"""
SanitizerStage: ITP Sanitizer and Charge Neutrality Preflight.

Runs after staging (htpolynet) and before any GROMACS grompp.

Implements:
- Atomtype extraction and de-duplication from all ITPs
- Conflict detection and namespace collision defense
- [ defaults ] / comb-rule validation and consistency
- Combined atomtypes file generation (versioned + current)
- Sanitized ITP generation (atomtypes sections removed)
- Charge neutrality preflight check
- Optional auto-correction for small charge imbalances
- system.top update with correct include order (defaults first)

Hardening (v2):
- Origin-based source_ff mapping (not global ctx.ff)
- Deterministic file ordering for reproducibility
- Proper file classification for conflict resolution
- Recursive include dependency discovery

Hardening (v3):
- PBC-aware dipole shift calculation (minimum-image unwrapping)
- Configurable grompp preprocessing timeout
- None-safe charge limit comparisons
- Duplicate [molecules] entry handling
- HTPolyNet/staged ITP sources in missing ITP validation
- Flexible forcefield directory resolution
- Configurable LJ outlier thresholds
- Preprocessor-aware moleculetype scanning
"""

from dataclasses import dataclass, asdict
from collections import deque
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple
import math
import hashlib
import os
import re
import shutil
import subprocess

from .base import BaseStage
from ..itp_sanitizer import IncludeResolution

if TYPE_CHECKING:
    from ..context import PipelineContext


class SanitizerError(Exception):
    """Error during sanitization stage."""
    pass


# Default timeout for grompp preprocessing (seconds)
DEFAULT_GROMPP_PP_TIMEOUT_S = 300

# Hardening v7: Safety cap for absurdly large GRO files (non-strict parse)
GRO_ABSURD_ATOMS = 5_000_000

# Hardening v6: Grompp timeout complexity estimation parameters
GROMPP_TIMEOUT_FLOOR_S = 60  # Minimum timeout for any system
GROMPP_TIMEOUT_PER_1K_LINES = 10  # +10s per 1000 lines of topology
GROMPP_TIMEOUT_PER_10_INCLUDES = 30  # +30s per 10 includes resolved
GROMPP_TIMEOUT_RETRY_MULTIPLIER = 2.0  # Retry with 2x timeout on first failure

# LJ outlier bound profiles (comb-rule 2/3; sigma/epsilon interpretation).
LJ_BOUNDS_PROFILES = {
    "aa": {
        "sigma_min_nm": 0.03,
        "sigma_max_nm": 1.0,
        "epsilon_max_kj_mol": 50.0,
    },
    "coarse_grained": {
        "sigma_min_nm": 0.0,
        "sigma_max_nm": 5.0,
        "epsilon_max_kj_mol": 5000.0,
    },
    "polarizable": {
        "sigma_min_nm": 0.0,
        "sigma_max_nm": 2.5,
        "epsilon_max_kj_mol": 500.0,
    },
}

# Always-fatal "units are almost certainly broken" bounds.
LJ_UNIT_SCREAMING_THRESHOLDS = {
    "sigma_max_nm": 10.0,
    "epsilon_max_kj_mol": 1e6,
}

# Hardening v6: Protected molecule patterns for charge correction
# These molecules should NEVER have their charges modified unless explicitly allowed
PROTECTED_SOLVENT_PATTERNS = frozenset({
    # Carbonate solvents (GPE)
    "PC", "EC", "DMC", "DEC", "EMC", "FEC", "VC",
    # Ether solvents
    "DOL", "DME", "TEGDME", "DIGLYME", "TRIGLYME",
    # Other common solvents
    "ACN", "AN", "DMSO", "NMP", "THF", "THP",
    # Water models
    "SOL", "HOH", "WAT", "TIP3P", "TIP4P", "SPC", "SPCE",
})

PROTECTED_ION_PATTERNS = frozenset({
    # Alkali cations
    "LI", "LI+", "NA", "NA+", "K", "K+",
    # Anions
    "TFSI", "TFSI-", "FSI", "FSI-", "PF6", "PF6-", "BF4", "BF4-",
    "CLO4", "CL", "CL-", "BR", "BR-", "I", "I-",
    # Salt molecule patterns
    "LITFSI", "LIFSI", "LIPF6", "LIBF4", "NACL", "KCL",
})

# Hardening v6: Include resolution priority defaults
# Canonical priorities:
# - forcefield_first (legacy "last"): search forcefield/system before user include dirs
# - sanitized_first (legacy "first"): allow user include dirs to override earlier
DEFAULT_INCLUDE_PRIORITY = "forcefield_first"
INCLUDE_PRIORITY_ALIASES = {
    "forcefield_first": "forcefield_first",
    "ff_first": "forcefield_first",
    "last": "forcefield_first",            # legacy alias
    "sanitized_first": "sanitized_first",
    "sanitizer_first": "sanitized_first",
    "first": "sanitized_first",            # legacy alias
}
CANONICAL_INCLUDE_PRIORITIES = frozenset({"sanitized_first", "forcefield_first"})
DEFAULT_ALLOW_INCLUDE_SHADOWING = False  # In strict mode, fail on shadowing
DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE = False

# Supported ptype tokens in [ atomtypes ].
PTYPE_TOKENS = frozenset({"A", "D", "V"})

# Sentinel markers for sanitizer-managed include block in system.top
# Hardening v6: Consistent marker wording
SANITIZER_BLOCK_BEGIN = "; >>> sanitizer managed block begin"
SANITIZER_BLOCK_END = "; <<< sanitizer managed block end"

# Deterministic priority for selecting conflicting moleculetype definitions.
# Lower value = higher priority.
MOLECULETYPE_SOURCE_PRIORITY = {
    "htpolynet": 0,
    "staged": 1,
    "molecule": 2,
    "molecules-lib": 2,  # Alias for readability in logs/docs
    "forcefield": 3,
}


# =============================================================================
# Parsing Utilities (v4 hardening)
# =============================================================================

def strip_gmx_comment(line: str) -> str:
    """
    Remove trailing '; comment' from a GROMACS line.
    
    Preserves preprocessor directives:
    - For #include: strips trailing '; comment' to allow parsing the target
    - For #ifdef/#define/etc: returns unchanged (comment may be part of macro)
    
    Args:
        line: Raw line from topology file
        
    Returns:
        Line with trailing comment removed (or unchanged for non-include preprocessor)
    """
    stripped = line.lstrip()
    if stripped.startswith("#"):
        # For #include: allow stripping trailing comment
        if stripped.startswith("#include"):
            idx = line.find(";")
            return line[:idx].rstrip() if idx >= 0 else line
        # Other preprocessor lines: preserve completely
        return line
    idx = line.find(";")
    return line[:idx] if idx >= 0 else line


def parse_section_name(line: str) -> Optional[str]:
    """
    Parse a GROMACS section header, supporting trailing comments.
    
    Examples:
        '[ atoms ]' -> 'atoms'
        '[ moleculetype ] ; comment' -> 'moleculetype'
        '[bonds]' -> 'bonds'
        'not a section' -> None
    
    Args:
        line: Raw line from topology file
        
    Returns:
        Section name in lowercase, or None if not a section header
    """
    stripped = strip_gmx_comment(line).strip()
    if stripped.startswith("[") and stripped.endswith("]"):
        return stripped[1:-1].strip().lower()
    return None


def parse_include_target(line: str) -> Optional[str]:
    """
    Parse #include directive, robust to trailing comments and both quote styles.
    
    Examples:
        '#include "topol.itp"' -> 'topol.itp'
        '#include <forcefield.itp>' -> 'forcefield.itp'
        '#include "path.itp" ; Include topology' -> 'path.itp'
    
    Args:
        line: Raw line from topology file
        
    Returns:
        Include target path, or None if not an include directive
    """
    stripped = strip_gmx_comment(line).strip()
    if not stripped.startswith("#include"):
        return None
    match = re.match(r'#include\s+[<"]([^">]+)[">]', stripped)
    return match.group(1) if match else None


def normalize_include_priority(raw_value: Optional[str]) -> str:
    """Normalize include-priority value to canonical semantics."""
    value = str(raw_value).strip().casefold() if raw_value is not None else ""
    if not value:
        return DEFAULT_INCLUDE_PRIORITY
    mapped = INCLUDE_PRIORITY_ALIASES.get(value)
    if mapped is not None:
        return mapped
    valid = ", ".join(["sanitized_first", "forcefield_first", "first", "last"])
    raise SanitizerError(
        f"Invalid --itp-include-dirs-priority value '{raw_value}'. "
        f"Valid values: {valid}"
    )


def _matches_protected_pattern(name: str, patterns: frozenset) -> bool:
    """
    Check if name matches any protected pattern using boundary-aware matching.
    
    Hardening v7: Short patterns (<=3 chars) require exact match after normalization
    to avoid false positives like "PC" matching "POLYMER_CHAIN".
    
    Longer patterns can match as substring but require word boundaries (non-alphanumeric
    or string start/end on both sides).
    
    Args:
        name: Molecule name to check
        patterns: Set of protected patterns (e.g., PROTECTED_SOLVENT_PATTERNS)
        
    Returns:
        True if name matches any pattern
    """
    # Local import keeps this helper self-contained for test extraction.
    import re as _re
    # Normalize: uppercase and remove NON-alphanumerics so "Li+" and "LI" compare equal.
    # This defends against false negatives when charges/punctuation are present.
    raw_upper = name.upper()
    normalized = _re.sub(r"[^A-Z0-9]", "", raw_upper)
    
    for pattern in patterns:
        # Normalize patterns the same way as names
        pat_upper = pattern.upper()
        pat_upper = _re.sub(r"[^A-Z0-9]", "", pat_upper)
        if not pat_upper:
            continue
        
        if len(pat_upper) <= 3:
            # Short patterns: exact match only after normalization
            # This prevents "PC" from matching "POLYMER_CHAIN" or "LIPID_GPC"
            if normalized == pat_upper:
                return True
        else:
            # Longer patterns: allow substring but check boundaries on ORIGINAL string.
            # Using the raw string preserves separators so boundary detection is reliable.
            boundary_re = rf"(?<![A-Z0-9]){_re.escape(pat_upper)}(?![A-Z0-9])"
            if _re.search(boundary_re, raw_upper):
                return True
    return False


@dataclass
class GROAtom:
    """Atom data from a GRO file line."""
    resnr: int
    resname: str
    atomname: str
    x: float
    y: float
    z: float


@dataclass
class GROParseResult:
    """Result of parsing a GRO file with box vectors."""
    atoms: List[GROAtom]
    box_x: float
    box_y: float
    box_z: float
    is_triclinic: bool
    title: str = ""



class SanitizerStage(BaseStage):
    """
    Stage 2.5: ITP Sanitizer (between htpolynet and gromacs)
    
    Reads:
    - IN/systems/<SYSTEM_ID>/gromacs/itp/*.itp (staged ITPs)
    - IN/systems/<SYSTEM_ID>/htpolynet/itp/*.itp (htpolynet ITPs)
    - IN/molecules/<NAME>/itp/*.itp (referenced molecule ITPs)
    - IN/forcefield/<FF>/gromacs/*.itp (forcefield atomtypes)
    
    Writes:
    - IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_<RUN_ID>.itp
    - IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_current.itp (link)
    - IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_<RUN_ID>/
    - IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_current/ (link)
    - IN/systems/<SYSTEM_ID>/gromacs/top/system.top (updated)
    """
    
    @property
    def name(self) -> str:
        return "sanitizer"
    
    @property
    def output_subdir(self) -> str:
        return "02_5_sanitizer"
    
    def run(self, ctx: "PipelineContext") -> bool:
        """
        Execute ITP sanitization and charge neutrality preflight.
        
        Returns:
            True if successful, False otherwise
        """
        from ..itp_sanitizer import (
            ItpSanitizer,
            ItpSanitizerError,
            ItpParser,
            generate_system_top,
        )
        from ..charge_neutrality import (
            ChargeNeutralityChecker,
            ChargeNeutralityError,
        )
        
        print("  Running ITP Sanitizer stage...")
        self.diagnostics: List[IncludeResolution] = []
        # Reset per-run cached paths for safe instance reuse.
        self._forcefield_dir = None
        self._sanitized_current_dir = None
        
        # Get paths
        system_gromacs_dir = ctx.get_input_path(
            "systems", ctx.system_id, "gromacs"
        )
        system_htpolynet_dir = ctx.get_input_path(
            "systems", ctx.system_id, "htpolynet"
        )
        molecules_dir = ctx.get_input_path("molecules")
        
        # Task F: Flexible forcefield directory resolution
        # Try both dashed and non-dashed directory names
        ff_base = ctx.ff.lower()
        ff_candidates = [
            ctx.get_input_path("forcefield", ff_base.replace("-", "")),
            ctx.get_input_path("forcefield", ff_base),
        ]
        forcefield_dir = None
        for ff_candidate in ff_candidates:
            if ff_candidate.exists():
                forcefield_dir = ff_candidate
                break
        if forcefield_dir is None:
            raise SanitizerError(
                f"Forcefield directory not found. Tried:\n"
                + "\n".join(f"  - {c}" for c in ff_candidates)
            )
        self._forcefield_dir = forcefield_dir
        
        print(f"  - System GROMACS dir: {system_gromacs_dir}")
        print(f"  - Forcefield: {ctx.ff} (dir: {forcefield_dir.name})")
        
        # Determine molecule counts early (ordered if possible)
        top_path = system_gromacs_dir / "top" / "system.top"
        ordered_molecules = self._get_ordered_molecules_from_top(top_path)
        if ordered_molecules:
            molecule_counts = {name: count for name, count in ordered_molecules}
            htpolynet_molecules_present = True
            print("  Using [molecules] from HTPolyNet-staged system.top")
        else:
            molecule_counts = self._get_molecule_counts(ctx)
            ordered_molecules = list(molecule_counts.items())
            htpolynet_molecules_present = False
            print("  Using [molecules] from manifest (pre-reaction counts)")
        
        needed_molecule_names = {
            name for name, count in ordered_molecules if count > 0
        }
        
        # Collect all ITP files with proper classification
        itp_paths: List[Path] = []
        ff_map: Dict[str, str] = {}  # file_path -> forcefield identifier
        class_map: Dict[str, str] = {}  # file_path -> file class
        
        # 1. Forcefield atomtypes (if exists) - highest priority
        ff_gromacs_dir = forcefield_dir / "gromacs"
        if ff_gromacs_dir.exists():
            ff_itps = self._sorted_paths_casefold(ff_gromacs_dir.glob("*.itp"))
            for itp in ff_itps:
                itp_paths.append(itp)
                ff_map[str(itp)] = ctx.ff  # Only FF files get ctx.ff
                class_map[str(itp)] = "forcefield"
                print(f"    + FF: {itp.name}")
            
            # Resolve includes from forcefield files
            for itp in ff_itps:
                include_dirs = self._build_include_search_paths(ctx, itp)
                includes = ItpParser.resolve_includes(
                    itp,
                    include_dirs=include_dirs,
                    strict=ctx.strict_include_resolution,
                    diagnostics=self.diagnostics,
                    check_shadowing=True,
                    allow_shadowing=ctx.allow_include_shadowing,
                )
                allow_escape = bool(
                    getattr(ctx, "allow_unsafe_include_escape", DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE)
                )
                allowed_roots = self._include_allowed_roots(
                    including_file=itp,
                    include_dirs=include_dirs,
                    ctx=ctx,
                )
                for inc in includes:
                    if not allow_escape and not self._include_path_allowed(inc.resolve(), allowed_roots):
                        roots_text = ", ".join(str(r) for r in allowed_roots) or "(none)"
                        msg = (
                            f"Blocked include escape '{inc}' while scanning {itp.name}.\n"
                            f"  Allowed roots: {roots_text}\n"
                            "Use --allow-unsafe-include-escape to override."
                        )
                        if ctx.strict_include_resolution:
                            raise SanitizerError(msg)
                        print(f"  [WARN] {msg}")
                        continue
                    if str(inc) not in ff_map:
                        itp_paths.append(inc)
                        ff_map[str(inc)] = ctx.ff
                        class_map[str(inc)] = "forcefield"
                        print(f"    + FF (included): {inc.name}")
        
        # 2. Staged system ITPs from gromacs/itp/
        staged_itp_dir = system_gromacs_dir / "itp"
        if staged_itp_dir.exists():
            staged_itps = self._sorted_paths_casefold(staged_itp_dir.glob("*.itp"))
            for itp in staged_itps:
                if str(itp) not in ff_map:  # Avoid duplicates
                    itp_paths.append(itp)
                    ff_map[str(itp)] = "UNKNOWN"  # Not from forcefield
                    class_map[str(itp)] = "staged"
                    print(f"    + Staged: {itp.name}")
        
        # 3. HTPolyNet produced ITPs
        htpolynet_itp_dir = system_htpolynet_dir / "itp"
        if htpolynet_itp_dir.exists():
            htpolynet_itps = self._sorted_paths_casefold(htpolynet_itp_dir.glob("*.itp"))
            for itp in htpolynet_itps:
                if str(itp) not in ff_map:  # Avoid duplicates
                    itp_paths.append(itp)
                    ff_map[str(itp)] = "UNKNOWN"  # From HTPolyNet, not FF library
                    class_map[str(itp)] = "htpolynet"
                    print(f"    + HTPolyNet: {itp.name}")
        
        # 4. Molecule ITPs (molecules library)
        needed_name_map = {n.casefold(): n for n in needed_molecule_names}
        if molecules_dir.exists():
            mol_dirs = sorted(molecules_dir.iterdir(), key=lambda p: p.name.casefold())
            for mol_dir in mol_dirs:
                if not mol_dir.is_dir():
                    continue
                mol_itp_dir = mol_dir / "itp"
                if mol_itp_dir.exists():
                    mol_itps = self._sorted_paths_casefold(mol_itp_dir.glob("*.itp"))
                    for itp in mol_itps:
                        if str(itp) in ff_map:
                            continue
                        mol_names = self._get_moleculetype_names(itp, ctx=ctx)
                        matched_names = sorted(
                            {
                                needed_name_map[m.casefold()]
                                for m in mol_names
                                if m.casefold() in needed_name_map
                            },
                            key=lambda n: n.casefold(),
                        )
                        if not matched_names:
                            continue
                        if str(itp) not in ff_map:
                            itp_paths.append(itp)
                            ff_map[str(itp)] = "UNKNOWN"
                            class_map[str(itp)] = "molecule"
                            print(
                                f"    + Molecules-lib {','.join(matched_names)}: {itp.name}"
                            )
        
        # Hardening v7: Expand include-closure for ALL ITPs (not just forcefield)
        # This ensures moleculetypes in included files are discovered
        print("  Expanding include-closure for all ITPs...")
        itp_paths, class_map, ff_map = self._expand_itp_include_closure(
            itp_paths, ctx, class_map, ff_map
        )
        
        # Build exact-name moleculetype candidates from all sources
        all_candidates: Dict[str, List[Path]] = {}
        casefold_index: Dict[str, Set[str]] = {}  # casefold -> set of exact names
        
        for itp in itp_paths:
            for mol_name in self._get_moleculetype_names(itp, ctx=ctx):
                all_candidates.setdefault(mol_name, []).append(itp)
                casefold_index.setdefault(mol_name.casefold(), set()).add(mol_name)
        
        # Fail-fast on case collisions (different exact names that casefold the same)
        for cf_key, exact_names in casefold_index.items():
            if len(exact_names) > 1:
                raise SanitizerError(
                    f"Case-insensitive moleculetype collision will break grompp:\n"
                    f"  Names: {sorted(exact_names)}\n"
                    f"  Rename one of these moleculetypes to avoid OS-dependent behavior."
                )

        # Select one deterministic winner per exact moleculetype name.
        all_moleculetypes: Dict[str, Path] = {}
        for mol_name in sorted(all_candidates, key=lambda n: (n.casefold(), n)):
            winner = self._select_moleculetype_winner(
                mol_name=mol_name,
                candidate_paths=all_candidates[mol_name],
                class_map=class_map,
                ctx=ctx,
            )
            all_moleculetypes[mol_name] = winner
        
        # Map needed molecules from the single winner map above (no contradictory paths).
        molecule_itp_paths: Dict[str, Path] = {}
        if needed_molecule_names:
            still_missing = []
            for n in sorted(needed_molecule_names, key=lambda x: x.casefold()):
                exact_names = sorted(
                    casefold_index.get(n.casefold(), set()),
                    key=lambda x: (x.casefold(), x),
                )
                if not exact_names:
                    still_missing.append(n)
                    continue
                winner = all_moleculetypes.get(exact_names[0])
                if winner is None:
                    still_missing.append(n)
                    continue
                molecule_itp_paths[n] = winner
            
            if still_missing:
                missing_sorted = sorted(still_missing, key=lambda x: x.casefold())
                # Count ITPs per source category
                ff_count = sum(1 for p, c in class_map.items() if c == "forcefield")
                staged_count = sum(1 for p, c in class_map.items() if c == "staged")
                htpoly_count = sum(1 for p, c in class_map.items() if c == "htpolynet")
                mol_count = sum(1 for p, c in class_map.items() if c == "molecule")
                
                raise SanitizerError(
                    f"Missing moleculetype(s) for [molecules]: {missing_sorted}\n"
                    f"Searched in:\n"
                    f"  - Forcefield: {forcefield_dir}/gromacs/ ({ff_count} files)\n"
                    f"  - Staged: {system_gromacs_dir}/itp/ ({staged_count} files)\n"
                    f"  - HTPolyNet: {system_htpolynet_dir}/itp/ ({htpoly_count} files)\n"
                    f"  - Molecules: {molecules_dir}/ ({mol_count} files)\n"
                    f"Suggestions:\n"
                    f"  - Verify HTPolyNet stage outputs\n"
                    f"  - Check #include paths in topology\n"
                    f"  - Run `grompp -pp` to see full preprocessed topology"
                )
        
        # Enforce case-insensitive uniqueness of ITP basenames (cross-platform reproducibility)
        self._assert_no_case_insensitive_duplicates(itp_paths, "ITP inputs")
        
        if not itp_paths:
            print("  [WARN] No ITP files found to sanitize")
            print("  Creating empty combined atomtypes file...")
            
            # Create minimal outputs
            self._create_minimal_outputs(ctx, system_gromacs_dir, allow_default_defaults=ctx.allow_default_defaults)
            return True
        
        print(f"  Found {len(itp_paths)} ITP files to process")
        
        # Get options from context
        allow_override = getattr(ctx, 'allow_override', False)
        atomtype_prefix = getattr(ctx, 'atomtype_prefix', None)
        strict_charge = getattr(ctx, 'strict_charge', False)
        
        # Initialize sanitizer
        sanitizer = ItpSanitizer(
            source_ff=ctx.ff,
            allow_override=allow_override,
            atomtype_prefix=atomtype_prefix,
            strict_charge=strict_charge,
        )
        
        try:
            # Run sanitization
            include_roots = [top_path] if top_path.exists() else []
            result = sanitizer.run(
                itp_paths=itp_paths,
                output_dir=system_gromacs_dir,
                run_id=ctx.run_id,
                ff_map=ff_map,
                class_map=class_map,
                ctx=ctx,
                include_roots=include_roots,
                allow_default_defaults=ctx.allow_default_defaults,
            )
            self._sanitized_current_dir = result.sanitized_current_dir
            
            print(f"  [OK] Sanitization complete:")
            print(f"       Atomtypes: {result.atomtype_count}")
            print(f"       Files sanitized: {len(result.sanitized_files)}")
            print(f"       Molecule ITPs: {len(result.molecule_itp_files)}")
            
            if result.conflicts_overridden:
                print(f"  [WARN] {len(result.conflicts_overridden)} conflicts overridden")
                for conflict in result.conflicts_overridden:
                    winner_file = conflict.winner.source_files[0] if conflict.winner and conflict.winner.source_files else "?"
                    print(f"         - {conflict.atomtype_name}: winner = {Path(winner_file).name}")
            
            if result.charge_warnings:
                print(f"  [WARN] {len(result.charge_warnings)} charge warnings (non-fatal)")
            
            if result.defaults_used:
                print(f"  [ defaults ] comb-rule: {result.defaults_used.comb_rule}")
            if getattr(result, "mixed_defaults_detected", False):
                defaults_report = getattr(result, "mixed_defaults_report", {})
                print("  [WARN] Mixed [defaults] tuples detected in sources.")
                print(
                    "  [WARN] Classification: "
                    f"{defaults_report.get('classification', 'unknown')} "
                    f"(fields: {defaults_report.get('differing_fields', [])})"
                )
                print(
                    "  [WARN] Primary policy: "
                    f"{defaults_report.get('chosen_policy', 'unknown')}, "
                    f"primary={defaults_report.get('primary_defaults_signature', 'n/a')}"
                )
                if getattr(ctx, "allow_mixed_defaults", False):
                    print("  [WARN] Continuing due to explicit unsafe override (--allow-mixed-defaults).")
                    print(
                        "  [WARN] Override reason: "
                        f"{getattr(ctx, 'allow_mixed_defaults_reason', '')}"
                    )
            if getattr(result, "nonbond_params_secondary_current_path", None):
                print(
                    "  [WARN] Added sanitizer-generated within-group preserve include: "
                    f"{result.nonbond_params_secondary_current_path}"
                )
            if getattr(result, "nonbond_params_cross_group_current_path", None):
                print(
                    "  [WARN] Added sanitizer-generated cross-group include: "
                    f"{result.nonbond_params_cross_group_current_path}"
                )
            prefix_summary = getattr(result, "prefix_policy_summary", {})
            if prefix_summary:
                print(
                    "  Prefix policy: "
                    f"{prefix_summary.get('policy')} "
                    f"-> {prefix_summary.get('applied_strategy')}"
                )
            if getattr(result, "prefix_injected_types_current_path", None):
                print(
                    "  [WARN] Added sanitizer-generated prefixed [*types] include: "
                    f"{result.prefix_injected_types_current_path}"
                )
            
        except ItpSanitizerError as e:
            print(f"  [FAIL] Sanitization failed: {e}")
            raise SanitizerError(str(e))

        if result.defaults_used is None:
            raise SanitizerError(
                "No [ defaults ] block found in inputs. Refusing to invent defaults.\n"
                "Provide explicit defaults in input ITPs or enable --allow-default-defaults."
            )
        if getattr(result, "defaults_guessed", False):
            print("  [WARN] Using fallback [ defaults ] per ctx.allow_default_defaults")
        
        cross_ff_conflicts = [
            c for c in result.conflicts_detected
            if c.conflict_type == "cross_ff_collision"
        ]
        if cross_ff_conflicts:
            msg = (
                f"Detected {len(cross_ff_conflicts)} cross-forcefield atomtype conflicts. "
                "This suggests mixed forcefield profiles."
            )
            if ctx.strict_forcefield_consistency:
                raise SanitizerError(msg)
            print(f"  [WARN] {msg}")
        
        # v4: Pass comb_rule for comb-rule-aware LJ outlier detection
        comb_rule = result.defaults_used.comb_rule if result.defaults_used else None
        outlier_warnings = self._warn_atomtype_outliers(result.combined_atomtypes_current_path, ctx, comb_rule)
        
        # v4 Issue 3: Invoke gen-pairs/[pairs] consistency check
        # Scan all source ITPs (not just sanitized_files) to catch FF/staged/htpolynet-derived pairs
        if result.defaults_used:
            self._warn_pairs_consistency(
                result.defaults_used,
                itp_paths,
                source_defaults_map=getattr(result, "source_defaults_map", None),
                ctx=ctx,
            )
        
        # === Charge Neutrality Check ===
        
        # Fail fast if no molecule counts available
        if not molecule_counts:
            print("  [ERROR] No molecule counts available from HTPolyNet or manifest")
            print("         Run HTPolyNet stage first or ensure manifest has composition data")
            raise SanitizerError("No molecule counts available for charge neutrality check")
        
        sanitized_molecule_itp_paths = self._map_sanitized_itps_to_molecules(
            result.sanitized_current_dir,
            needed_molecule_names,
        )
        
        if sanitized_molecule_itp_paths:
            print("  Running charge neutrality preflight...")
            
            # Parse target allowlist from comma-separated string
            target_allowlist = self._parse_csv_casefold_set(ctx.charge_fix_target_allowlist)
            
            # Part 2: Parse new comma-separated settings
            allowed_atomnames = None
            if ctx.charge_fix_target_atomnames:
                allowed_atomnames = set(
                    s.strip() for s in ctx.charge_fix_target_atomnames.split(",") if s.strip()
                )
            
            disallow_atomtypes = None
            if ctx.charge_fix_disallowed_atomtypes:
                disallow_atomtypes = set(
                    s.strip() for s in ctx.charge_fix_disallowed_atomtypes.split(",") if s.strip()
                )
            
            # Hardening v5 - Issue G: Pre-classify ionic moleculetypes
            ion_classifications = self._classify_ionic_moleculetypes(
                sanitized_molecule_itp_paths, ctx
            )

            class_by_casefold = {
                name.casefold(): cls for name, cls in ion_classifications.items()
            }
            protected_classes = {
                "protected_ion",
                "ionic",
                "monatomic_ion",
                "protected_solvent",
                "protected_polymer",
            }
            protected_detected = any(
                cls in protected_classes for cls in ion_classifications.values()
            )

            # Optional explicit correction target(s) from CLI. Checker supports one target;
            # we resolve deterministically and fail-fast on ambiguous strict requests.
            requested_target: Optional[str] = None
            requested_targets_raw = self._parse_csv_casefold_set(
                getattr(ctx, "charge_fix_target_molecules", None)
            )
            if requested_targets_raw:
                available_by_casefold = {
                    name.casefold(): name for name in sanitized_molecule_itp_paths.keys()
                }
                resolved_targets: List[str] = []
                unknown_targets: List[str] = []
                for token_cf in sorted(requested_targets_raw):
                    resolved = available_by_casefold.get(token_cf)
                    if resolved is None:
                        unknown_targets.append(token_cf)
                        continue
                    resolved_targets.append(resolved)
                if unknown_targets:
                    msg = (
                        "Unknown --charge-fix-target-molecules entries: "
                        f"{unknown_targets}. Available molecules: "
                        f"{sorted(sanitized_molecule_itp_paths.keys(), key=str.casefold)}"
                    )
                    if ctx.strict_charge_neutrality:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                if len(resolved_targets) > 1:
                    msg = (
                        "Multiple explicit charge-fix targets provided; only one target is supported "
                        f"per correction run: {resolved_targets}"
                    )
                    if ctx.strict_charge_neutrality:
                        raise SanitizerError(msg)
                    resolved_targets = sorted(set(resolved_targets), key=str.casefold)
                    requested_target = resolved_targets[0]
                    print(
                        f"  [WARN] {msg}. Using deterministic first target: '{requested_target}'."
                    )
                elif len(resolved_targets) == 1:
                    requested_target = resolved_targets[0]

            # Strict safety hardening: when protected species exist, require explicit target allowlist.
            if ctx.strict_charge_neutrality and protected_detected and target_allowlist is None:
                raise SanitizerError(
                    "Strict charge-neutrality mode requires an explicit "
                    "--charge-fix-target-allowlist when protected ions/solvents/polymers are present."
                )
            
            # Task C: Coalesce None to safe defaults at point of use
            max_delta_limit = (
                ctx.charge_fix_max_delta_per_atom
                if ctx.charge_fix_max_delta_per_atom is not None
                else 1e-5
            )
            max_total_limit = (
                ctx.charge_fix_max_total
                if ctx.charge_fix_max_total is not None
                else 1e-4
            )
            dipole_limit = (
                ctx.charge_fix_max_dipole_shift_debye
                if ctx.charge_fix_max_dipole_shift_debye is not None
                else 0.01
            )
            
            try:
                # Part 2: Use new context fields for threshold and correction limits
                checker = ChargeNeutralityChecker(
                    threshold=ctx.charge_neutrality_correct,  # Max |Q| for correction
                    max_delta_per_atom=max_delta_limit,
                    max_total_delta=max_total_limit,
                    max_dipole_shift_debye=dipole_limit,
                    strict_charge_physics=ctx.strict_charge_physics,
                    # Safe-subset config
                    correction_method=getattr(ctx, "charge_fix_method", "safe_subset"),
                    target_allowlist=target_allowlist,
                    allowed_atomnames=allowed_atomnames,
                    disallow_atomtypes=disallow_atomtypes,
                    # Part 2: New strict mode and threshold options
                    strict_mode=ctx.strict_charge_neutrality,
                    allow_ions=ctx.charge_fix_allow_ions,
                    warn_threshold=ctx.charge_neutrality_warn,
                    neutrality_tol=ctx.charge_neutrality_tol,
                    protected_polymers={
                        name for name, cls in ion_classifications.items() if cls == "protected_polymer"
                    },
                    polymer_net_charge_tol=getattr(ctx, "polymer_net_charge_tol", 1e-3),
                    polymer_correction_method=getattr(ctx, "charge_fix_polymer_method", "skip_if_small"),
                    polymer_exclusion_bonds=getattr(ctx, "charge_fix_polymer_exclusion_bonds", 2),
                    rounding_moleculetype_tol=getattr(ctx, "charge_fix_moleculetype_rounding_tol", 1e-4),
                    allow_non_rounding_correction=getattr(ctx, "charge_fix_allow_non_rounding", False),
                    allow_unparsed_atoms_lines=getattr(ctx, "allow_unparsed_atoms_lines", False),
                )
                correction_result = checker.run(
                    molecule_itp_paths=sanitized_molecule_itp_paths,
                    molecule_counts=molecule_counts,
                    output_dir=result.sanitized_current_dir,
                    run_id=ctx.run_id,
                    target_molecule=requested_target,
                )
                
                charge_summary = {
                    "computed_q": correction_result.original_q,
                    "corrected": correction_result.corrected,
                    "method": correction_result.correction_method,
                    "effective_method": correction_result.effective_method,
                    "method_desc": correction_result.method_desc,
                    "patched_itp": str(correction_result.patched_itp_path)
                        if correction_result.patched_itp_path else None,
                    "max_delta_applied": correction_result.max_delta_applied,
                    "total_delta": correction_result.total_delta,
                    "atoms_adjusted": correction_result.atoms_adjusted,
                    "target_molecule": correction_result.target_molecule,
                    "target_itp_moleculetype_name": correction_result.target_itp_moleculetype_name,
                    "patched_atoms": correction_result.patched_atoms,
                    # New fields for enhanced manifest logging (Issues A/B)
                    "safe_subset_used": correction_result.safe_subset_used,
                    "patched_atom_details": correction_result.patched_atom_details,
                    "patched_atom_details_truncated": correction_result.patched_atom_details_truncated,
                    "correction_refused": correction_result.correction_refused,
                    "refuse_reason": correction_result.refuse_reason,
                    "unparsed_atoms_lines": correction_result.unparsed_atoms_lines,
                    "unparsed_atoms_lines_truncated": correction_result.unparsed_atoms_lines_truncated,
                    "polymer_policy": correction_result.polymer_policy,
                    "rounding_only_passed": correction_result.rounding_only_passed,
                    "rounding_moleculetype_tol": correction_result.rounding_moleculetype_tol,
                    "per_moleculetype_report": correction_result.per_moleculetype_report,
                    "dominant_charge_contributors": correction_result.dominant_charge_contributors,
                    "hetero_atoms_modified": correction_result.hetero_atoms_modified,
                    "hetero_atoms_modified_count": correction_result.hetero_atoms_modified_count,
                    "hetero_atoms_modified_details": correction_result.hetero_atoms_modified_details,
                    # Record actual limits used (Task C)
                    "limits_used": {
                        "max_delta_per_atom": max_delta_limit,
                        "max_total": max_total_limit,
                        "max_dipole_shift_debye": dipole_limit,
                    },
                    "strict_mode_enforced": bool(ctx.strict_charge_neutrality),
                }
                electrolyte_warning = self._warn_electrolyte_fixed_charge_limitations(
                    correction_result=correction_result,
                    classifications=ion_classifications,
                    ctx=ctx,
                )
                if electrolyte_warning:
                    charge_summary["electrolyte_limitations"] = electrolyte_warning
                
                if correction_result.corrected:
                    # Task C: Use coalesced limits for comparison
                    if correction_result.max_delta_applied is not None and (
                        correction_result.max_delta_applied > max_delta_limit
                    ):
                        raise SanitizerError(
                            f"Charge correction exceeds per-atom limit: |Δq|={correction_result.max_delta_applied:.3e} "
                            f"> {max_delta_limit:.3e} (limit used)"
                        )
                    # v4 Issue 5: Use abs(total_delta) and guard for None
                    if (correction_result.total_delta is not None 
                        and max_total_limit is not None
                        and abs(correction_result.total_delta) > max_total_limit
                    ):
                        raise SanitizerError(
                            f"Charge correction exceeds total limit: |ΔQ|={abs(correction_result.total_delta):.3e} "
                            f"> {max_total_limit:.3e} (limit used)"
                        )
                    
                    # v4 Issue 5: Warn for uniform_all correction method
                    if (
                        getattr(correction_result, "effective_method", None) == "uniform_all"
                        or str(correction_result.correction_method or "").startswith("uniform_all:")
                    ):
                        print(
                            f"  [WARN] uniform_all correction applied (total_delta={correction_result.total_delta:.3e}). "
                            "Verify protected patterns/atomtypes were not affected."
                        )
                    
                    dipole_info = self._compute_dipole_shift(
                        system_gromacs_dir=system_gromacs_dir,
                        ordered_molecules=ordered_molecules,
                        molecule_itp_paths=sanitized_molecule_itp_paths,
                        correction_result=correction_result,
                        ctx=ctx,  # Hardening v5: explicit ctx passing
                    )
                    if dipole_info is not None:
                        charge_summary["dipole"] = {
                            "computed": dipole_info.computed,
                            "dipole_before_debye": dipole_info.dipole_before_debye,
                            "dipole_after_debye": dipole_info.dipole_after_debye,
                            "delta_dipole_debye": dipole_info.delta_dipole_debye,
                            "avg_delta_dipole_debye": dipole_info.delta_dipole_debye,
                            "max_delta_dipole_debye": dipole_info.max_delta_dipole_debye,
                            "coordinates_file": dipole_info.coordinates_file,
                            "skip_reasons": list(getattr(dipole_info, "skip_reasons", [])),
                        }
                        if dipole_info.computed and dipole_info.delta_dipole_debye is not None:
                            max_delta_val = (
                                dipole_info.max_delta_dipole_debye
                                if dipole_info.max_delta_dipole_debye is not None
                                else dipole_info.delta_dipole_debye
                            )
                            print(
                                "  Dipole shift summary: "
                                f"avg |Δμ|={dipole_info.delta_dipole_debye:.4f} D, "
                                f"max |Δμ|={max_delta_val:.4f} D"
                            )
                            # v4: Dipole check is heuristic by default
                            # Only hard-fail if strict_dipole_check=True AND no gating reasons
                            if dipole_info.delta_dipole_debye > dipole_limit:
                                msg = (
                                    f"Dipole shift |Δμ|={dipole_info.delta_dipole_debye:.4f} D exceeds "
                                    f"limit {dipole_limit:.4f} D"
                                )
                                # Check gating: skip enforcement for giant/percolated molecules
                                skip_reasons = getattr(dipole_info, 'skip_reasons', [])
                                can_enforce = not skip_reasons
                                if ctx.strict_dipole_check and can_enforce:
                                    raise SanitizerError(msg)
                                else:
                                    reason_str = f"; skip reasons: {skip_reasons}" if skip_reasons else " (heuristic mode)"
                                    print(f"  [WARN] {msg}{reason_str}")
                    
                    self._sync_corrected_itp(
                        correction_result=correction_result,
                        molecule_itp_paths=sanitized_molecule_itp_paths,
                        sanitized_dir=result.sanitized_dir,
                        sanitized_current_dir=result.sanitized_current_dir,
                    )
                else:
                    # v4: Use ctx.charge_neutrality_tol instead of hard-coded 1e-10
                    allowed_skip_statuses = {"polymer_acceptable_drift", "polymer_skipped_non_strict"}
                    if (
                        abs(correction_result.original_q) > ctx.charge_neutrality_tol
                        and correction_result.neutrality_status not in allowed_skip_statuses
                    ):
                        raise SanitizerError(
                            f"Charge neutrality correction failed: {correction_result.correction_method}"
                        )
                
                # Hardening v5 - Issue G: Post-correction ion protection check
                protection_audit = self._verify_ion_protection(
                    correction_result, ion_classifications, ctx
                )
                violations = protection_audit.get("violations", [])
                charge_summary["target_policy"] = protection_audit
                if violations:
                    for v in violations:
                        print(f"  [ERROR] {v}")
                    raise SanitizerError(
                        f"Charge protection policy violations: {violations}. "
                        "Use --charge-fix-allow-ions/--charge-fix-allow-solvents, and for polymer "
                        "targets provide --charge-fix-target-allowlist with "
                        "--charge-fix-polymer-method spread_safe."
                    )
                
                # Record in manifest
                if ctx.manifest:
                    ctx.manifest.set_sanitizer_output("charge_neutrality", charge_summary)
                    
            except ChargeNeutralityError as e:
                print(f"  [FAIL] {e}")
                raise SanitizerError(str(e))
        else:
            print("  [SKIP] Charge neutrality check (no molecule ITPs)")
        
        # === Validate GRO/TOP consistency ===
        # Fail fast if coordinate file doesn't match topology molecule counts
        gro_dir = system_gromacs_dir / "gro"
        gro_path = gro_dir / "system.gro"
        if gro_path.exists() and molecule_counts and sanitized_molecule_itp_paths:
            if not self._validate_gro_top_consistency(
                gro_path,
                system_top_path=top_path if top_path.exists() else None,
                molecule_counts=molecule_counts,
                molecule_itp_paths=sanitized_molecule_itp_paths,
                ctx=ctx,
                output_dir=self.get_output_dir(ctx),
            ):
                raise SanitizerError(
                    f"GRO/TOP atom count mismatch. The coordinate file does not match "
                    f"the topology. This can happen if HTPolyNet crosslinking changed "
                    f"the molecule composition but the GRO was not updated."
                )
        
        # === Update system.top ===
        
        top_dir = system_gromacs_dir / "top"
        top_dir.mkdir(parents=True, exist_ok=True)
        
        system_top_path = top_dir / "system.top"
        molecule_itp_include_paths = self._ordered_molecule_itp_paths(
            ordered_molecules,
            sanitized_molecule_itp_paths,
            result.sanitized_current_dir,
        )
        extra_includes: List[str] = []
        prefix_injected_current = getattr(result, "prefix_injected_types_current_path", None)
        if prefix_injected_current and Path(prefix_injected_current).exists():
            extra_includes.append(
                Path(os.path.relpath(prefix_injected_current, top_dir)).as_posix()
            )
        cross_group_nonbond_current = getattr(result, "nonbond_params_cross_group_current_path", None)
        if cross_group_nonbond_current and Path(cross_group_nonbond_current).exists():
            extra_includes.append(
                Path(os.path.relpath(cross_group_nonbond_current, top_dir)).as_posix()
            )
        secondary_nonbond_current = getattr(result, "nonbond_params_secondary_current_path", None)
        if secondary_nonbond_current and Path(secondary_nonbond_current).exists():
            extra_includes.append(
                Path(os.path.relpath(secondary_nonbond_current, top_dir)).as_posix()
            )
        
        # Check if HTPolyNet-staged system.top already exists
        # If so, preserve its [molecules] section and only update includes
        self._system_top_update_status = {
            "mode": "generated_new",
            "degraded": False,
        }
        if system_top_path.exists() and htpolynet_molecules_present:
            print(f"  Updating existing system.top (preserving [molecules] section)...")
            top_content = self._update_system_top_includes(
                existing_top=system_top_path,
                combined_atomtypes_path=result.combined_atomtypes_current_path,
                molecule_itp_files=molecule_itp_include_paths,
                defaults=result.defaults_used,
                top_dir=top_dir,
                allow_default_defaults=ctx.allow_default_defaults,
                ctx=ctx,
                extra_includes=extra_includes,
            )
        else:
            # Generate new system.top with proper include order
            top_content = generate_system_top(
                combined_atomtypes_path=result.combined_atomtypes_current_path,
                sanitized_itp_paths=result.sanitized_files,
                system_name=ctx.system_id,
                molecule_counts=molecule_counts,
                defaults=result.defaults_used,
                molecule_itp_files=molecule_itp_include_paths,
                top_dir=top_dir,
                allow_guess_defaults=ctx.allow_default_defaults,
                extra_includes=extra_includes,
            )

        defaults_sections_in_top = sum(
            1 for line in top_content.splitlines() if parse_section_name(line) == "defaults"
        )
        if defaults_sections_in_top > 1:
            raise SanitizerError(
                "Generated system.top contains multiple [ defaults ] sections. "
                "Refusing to proceed because GROMACS accepts only one global defaults tuple."
            )
        
        previous_top = ""
        if system_top_path.exists():
            previous_top = self.safe_read_text(system_top_path, strict=True)
        if previous_top != top_content:
            self._atomic_write_text(system_top_path, top_content, encoding="utf-8")
            print(f"  Wrote: {system_top_path}")
        else:
            print(f"  system.top unchanged: {system_top_path}")
        self._validate_system_top_includes(system_top_path, ctx)
        
        # === Record in manifest ===
        
        if ctx.manifest:
            include_priority_raw = getattr(
                ctx,
                "itp_include_dirs_priority",
                DEFAULT_INCLUDE_PRIORITY,
            )
            include_priority_norm = normalize_include_priority(include_priority_raw)
            ctx.manifest.set_sanitizer_output(
                "include_resolution_policy",
                {
                    "itp_include_dirs_priority": include_priority_raw,
                    "normalized_priority": include_priority_norm,
                    "allow_include_shadowing": bool(getattr(ctx, "allow_include_shadowing", False)),
                    "allow_unsafe_include_escape": bool(
                        getattr(ctx, "allow_unsafe_include_escape", DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE)
                    ),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "system_top_update",
                {
                    "path": str(system_top_path),
                    "status": getattr(self, "_system_top_update_status", {}),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "combined_atomtypes",
                {
                    "versioned": str(result.combined_atomtypes_path),
                    "current": str(result.combined_atomtypes_current_path),
                    "atomtype_count": result.atomtype_count,
                }
            )
            ctx.manifest.set_sanitizer_output(
                "sanitized_itps",
                {
                    "versioned_dir": str(result.sanitized_dir),
                    "current_dir": str(result.sanitized_current_dir),
                    "file_count": len(result.sanitized_files),
                    "molecule_itp_count": len(result.molecule_itp_files),
                    "files": [str(f) for f in result.sanitized_files],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "conflicts",
                {
                    "detected": len(result.conflicts_detected),
                    "overridden": len(result.conflicts_overridden),
                    "allow_override_flag": allow_override,
                    "atomtype_prefix": atomtype_prefix,
                    "overridden_details": [
                        {
                            "atomtype": c.atomtype_name,
                            "winner": Path(c.winner.source_files[0]).name if c.winner and c.winner.source_files else None,
                        }
                        for c in result.conflicts_overridden
                    ],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "defaults",
                {
                    "used": result.defaults_used is not None,
                    "nbfunc": result.defaults_used.nbfunc if result.defaults_used else None,
                    "comb_rule": result.defaults_used.comb_rule if result.defaults_used else None,
                    "source_file": Path(result.defaults_used.source_file).name if result.defaults_used else None,
                    "guessed": getattr(result, "defaults_guessed", False),
                    "allow_mixed_defaults": bool(getattr(ctx, "allow_mixed_defaults", False)),
                    "allow_mixed_defaults_reason": getattr(ctx, "allow_mixed_defaults_reason", None),
                    "mixed_defaults_cross_group_policy": getattr(
                        ctx, "mixed_defaults_cross_group_policy", None
                    ),
                    "mixed_defaults_cross_group_rule": getattr(
                        ctx, "mixed_defaults_cross_group_rule", None
                    ),
                    "mixed_defaults_cross_group_max_pairs": getattr(
                        ctx, "mixed_defaults_cross_group_max_pairs", None
                    ),
                    "mixed_defaults_cross_group_reason": getattr(
                        ctx, "mixed_defaults_cross_group_reason", None
                    ),
                    "mixed_defaults_preserve_within_group": bool(
                        getattr(ctx, "mixed_defaults_preserve_within_group", False)
                    ),
                    "unsafe_mixed_defaults_override_used": bool(
                        getattr(ctx, "allow_mixed_defaults", False)
                        and getattr(result, "mixed_defaults_detected", False)
                    ),
                    "primary_defaults_signature": getattr(result, "mixed_defaults_report", {}).get(
                        "primary_defaults_signature"
                    ),
                    "secondary_signatures": getattr(result, "mixed_defaults_report", {}).get(
                        "secondary_signatures", []
                    ),
                    "chosen_policy": getattr(result, "mixed_defaults_report", {}).get("chosen_policy"),
                    "classification": getattr(result, "mixed_defaults_report", {}).get("classification"),
                    "differing_fields": getattr(result, "mixed_defaults_report", {}).get(
                        "differing_fields", []
                    ),
                    "comb_rules_present": getattr(result, "mixed_defaults_report", {}).get(
                        "comb_rules_present", []
                    ),
                    "grouped_report": getattr(result, "mixed_defaults_report", {}).get("groups", []),
                    "source_defaults": [
                        {
                            "source_file": str(path),
                            "nbfunc": entry.nbfunc,
                            "comb_rule": entry.comb_rule,
                            "gen_pairs": entry.gen_pairs,
                            "fudge_lj": entry.fudge_lj,
                            "fudge_qq": entry.fudge_qq,
                        }
                        for path, entry in sorted(
                            getattr(result, "source_defaults_map", {}).items(),
                            key=lambda item: item[0].as_posix().casefold(),
                        )
                    ],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "mixed_defaults_within_group_preservation",
                {
                    "versioned": str(result.nonbond_params_secondary_path)
                    if getattr(result, "nonbond_params_secondary_path", None)
                    else None,
                    "current": str(result.nonbond_params_secondary_current_path)
                    if getattr(result, "nonbond_params_secondary_current_path", None)
                    else None,
                    "summary": getattr(result, "nonbond_params_secondary_summary", {}),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "mixed_defaults_cross_group_override",
                {
                    "versioned": str(result.nonbond_params_cross_group_path)
                    if getattr(result, "nonbond_params_cross_group_path", None)
                    else None,
                    "current": str(result.nonbond_params_cross_group_current_path)
                    if getattr(result, "nonbond_params_cross_group_current_path", None)
                    else None,
                    "summary": getattr(result, "nonbond_params_cross_group_summary", {}),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "prefix_implicit_topology",
                {
                    "policy_summary": getattr(result, "prefix_policy_summary", {}),
                    "injected_types_versioned": str(result.prefix_injected_types_path)
                    if getattr(result, "prefix_injected_types_path", None)
                    else None,
                    "injected_types_current": str(result.prefix_injected_types_current_path)
                    if getattr(result, "prefix_injected_types_current_path", None)
                    else None,
                },
            )
            if bool(getattr(ctx, "allow_mixed_defaults", False)) and bool(
                getattr(result, "mixed_defaults_detected", False)
            ):
                unsafe_overrides = ctx.manifest.get("unsafe_overrides", [])
                if not isinstance(unsafe_overrides, list):
                    unsafe_overrides = []
                unsafe_overrides.append(
                    {
                        "kind": "mixed_defaults_override",
                        "unsafe_override": True,
                        "reason": getattr(ctx, "allow_mixed_defaults_reason", None),
                        "primary_defaults_signature": getattr(
                            result, "mixed_defaults_report", {}
                        ).get("primary_defaults_signature"),
                        "secondary_signatures": getattr(
                            result, "mixed_defaults_report", {}
                        ).get("secondary_signatures", []),
                        "chosen_policy": getattr(result, "mixed_defaults_report", {}).get(
                            "chosen_policy"
                        ),
                    }
                )
                cross_group_summary = getattr(result, "nonbond_params_cross_group_summary", {})
                if (
                    isinstance(cross_group_summary, dict)
                    and cross_group_summary.get("policy") == "generate"
                ):
                    unsafe_overrides.append(
                        {
                            "kind": "mixed_defaults_cross_group_generate",
                            "unsafe_override": True,
                            "reason": cross_group_summary.get("reason"),
                            "rule": cross_group_summary.get("rule"),
                            "pair_count": cross_group_summary.get("pair_count"),
                            "file": cross_group_summary.get("current_path"),
                        }
                    )
                ctx.manifest.set("unsafe_overrides", unsafe_overrides)
            ctx.manifest.set_sanitizer_output(
                "charge_warnings",
                {
                    "count": len(result.charge_warnings),
                    "strict_mode": strict_charge,
                    "warnings": [str(w) for w in result.charge_warnings],
                }
            )
            
            # v5 Hardening: Include resolution report (Issue 1)
            ctx.manifest.set_sanitizer_output(
                "include_resolution",
                [asdict(d) for d in self.diagnostics]
            )
            ctx.manifest.save()
        
        # Log output for stage directory
        output_dir = self.get_output_dir(ctx)
        log_path = output_dir / "sanitizer.log"
        
        log_lines = [
            "ITP Sanitizer completed",
            f"Atomtypes: {result.atomtype_count}",
            f"Files: {len(result.sanitized_files)}",
            f"Molecule ITPs: {len(result.molecule_itp_files)}",
            f"Conflicts: {len(result.conflicts_detected)}",
            f"Overridden: {len(result.conflicts_overridden)}",
            f"Charge warnings: {len(result.charge_warnings)}",
        ]
        
        if result.defaults_used:
            log_lines.append(f"Defaults: nbfunc={result.defaults_used.nbfunc}, comb-rule={result.defaults_used.comb_rule}")
        defaults_report = getattr(result, "mixed_defaults_report", {})
        if defaults_report:
            log_lines.append(
                "Defaults policy: "
                f"classification={defaults_report.get('classification')}, "
                f"chosen_policy={defaults_report.get('chosen_policy')}, "
                f"primary={defaults_report.get('primary_defaults_signature')}"
            )
        nonbond_summary = getattr(result, "nonbond_params_secondary_summary", {})
        if nonbond_summary:
            log_lines.append(
                "Secondary nonbond_params preserve: "
                f"requested={nonbond_summary.get('requested')}, "
                f"applied={nonbond_summary.get('applied')}, "
                f"pairs={nonbond_summary.get('pair_count')}, "
                f"reason={nonbond_summary.get('reason')}"
            )
        cross_group_summary = getattr(result, "nonbond_params_cross_group_summary", {})
        if cross_group_summary:
            log_lines.append(
                "Cross-group nonbond_params: "
                f"policy={cross_group_summary.get('policy')}, "
                f"rule={cross_group_summary.get('rule')}, "
                f"applied={cross_group_summary.get('applied')}, "
                f"pairs={cross_group_summary.get('pair_count')}, "
                f"status={cross_group_summary.get('status')}"
            )
        prefix_summary = getattr(result, "prefix_policy_summary", {})
        if prefix_summary:
            log_lines.append(
                "Prefix implicit-topology policy: "
                f"policy={prefix_summary.get('policy')}, "
                f"strategy={prefix_summary.get('applied_strategy')}, "
                f"effective_prefix={prefix_summary.get('effective_prefix')}"
            )
        
        if result.conflicts_overridden:
            log_lines.append("\nOverridden conflicts (winner selected by priority):")
            for c in result.conflicts_overridden:
                winner_file = Path(c.winner.source_files[0]).name if c.winner and c.winner.source_files else "?"
                log_lines.append(f"  - {c.atomtype_name}: {winner_file}")
        
        self._atomic_write_text(log_path, "\n".join(log_lines) + "\n", encoding="utf-8")
        
        print("  [OK] ITP Sanitizer stage complete")
        return True
    
    def _get_molecule_counts(self, ctx: "PipelineContext") -> Dict[str, int]:
        """Get molecule counts from manifest or composition (fallback)."""
        counts = {}
        
        if ctx.manifest:
            manifest_data = ctx.manifest.data
            composition = manifest_data.get("composition", {})
            # Prefer 'counts' (canonical), fallback to 'molecule_counts'
            counts = composition.get("counts") or composition.get("molecule_counts", {})
        
        return counts

    def _sorted_paths_casefold(self, paths) -> List[Path]:
        """Deterministic path sort with casefold to avoid OS-dependent ordering."""
        return sorted(
            list(paths),
            key=lambda p: (p.as_posix().casefold(), p.as_posix()),
        )

    @staticmethod
    def _parse_csv_casefold_set(raw: Optional[str]) -> Optional[Set[str]]:
        """Parse comma-separated values into a casefold-normalized set."""
        if not raw:
            return None
        parsed = {token.strip().casefold() for token in raw.split(",") if token.strip()}
        return parsed or None

    def _assert_no_case_insensitive_duplicates(self, paths: List[Path], label: str) -> None:
        """Fail fast on case-insensitive basename duplicates (non-reproducible)."""
        seen: Dict[str, Path] = {}
        for path in paths:
            key = path.name.casefold()
            if key in seen and seen[key].resolve() != path.resolve():
                raise SanitizerError(
                    f"Case-insensitive filename collision in {label}:\n"
                    f"  - {seen[key]}\n  - {path}\n"
                    "Rename one of the files to avoid OS-dependent behavior."
                )
            seen[key] = path

    def safe_read_text(self, path: Path, strict: bool) -> str:
        """
        Read text with consistent strict/non-strict behavior.

        Strict mode raises SanitizerError with path + root cause.
        Non-strict mode logs a warning and returns empty content.
        """
        try:
            return path.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            msg = f"Failed to read text file: {path} ({e})"
            if strict:
                raise SanitizerError(msg)
            print(f"  [WARN] {msg}")
            return ""

    def _moleculetype_source_class(
        self,
        path: Path,
        class_map: Dict[str, str],
    ) -> str:
        """Get stable source class for an ITP path from class_map."""
        return (
            class_map.get(str(path))
            or class_map.get(str(path.resolve()))
            or "unknown"
        )

    def _moleculetype_source_priority(
        self,
        source_class: str,
    ) -> int:
        """Map source class to deterministic priority rank."""
        return MOLECULETYPE_SOURCE_PRIORITY.get(source_class, 99)

    def _select_moleculetype_winner(
        self,
        mol_name: str,
        candidate_paths: List[Path],
        class_map: Dict[str, str],
        ctx: "PipelineContext",
    ) -> Path:
        """
        Select one deterministic winner for a moleculetype across all sources.

        Priority: htpolynet > staged > molecules-lib > forcefield.
        On differing signatures: strict mode fails; otherwise warn and pick winner.
        """
        deduped: Dict[Path, Path] = {}
        for p in candidate_paths:
            deduped[p.resolve()] = p.resolve()
        ordered = sorted(
            deduped.values(),
            key=lambda p: (
                self._moleculetype_source_priority(
                    self._moleculetype_source_class(p, class_map)
                ),
                p.as_posix().casefold(),
                p.as_posix(),
            ),
        )
        winner = ordered[0]
        winner_source = self._moleculetype_source_class(winner, class_map)
        conflicts: List[Path] = []

        for candidate in ordered[1:]:
            if self._moleculetype_signatures_differ(
                winner,
                candidate,
                mol_name,
                ctx=ctx,
            ):
                conflicts.append(candidate)

        if conflicts:
            conflict_lines = [
                f"  - {winner} ({winner_source}, selected by priority)"
            ]
            for conflict in conflicts:
                conflict_lines.append(
                    f"  - {conflict} ({self._moleculetype_source_class(conflict, class_map)})"
                )
            strict_conflict = bool(
                getattr(ctx, "strict_forcefield_consistency", False)
                or getattr(ctx, "strict_include_resolution", False)
                or getattr(ctx, "strict_charge_neutrality", False)
            )
            msg = (
                f"Conflicting moleculetype '{mol_name}' definitions detected:\n"
                + "\n".join(conflict_lines)
                + "\nResolution options:\n"
                + "  1. Rename one moleculetype so RESP/CM5 variants are explicit.\n"
                + "  2. Keep only one charge model in include inputs.\n"
                + "  3. Use non-strict mode to accept deterministic priority selection."
            )
            if strict_conflict:
                raise SanitizerError(msg)
            print(f"  [WARN] {msg}")

        return winner

    def _get_moleculetype_names(
        self,
        itp_path: Path,
        include_dirs: Optional[List[Path]] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> List[str]:
        """
        Extract moleculetype names from an ITP, following #include chains.
        
        Hardening v7: Now supports proper include resolution using GROMACS -I semantics.
        
        Args:
            itp_path: Path to ITP file
            include_dirs: Optional list of include directories to search
            ctx: PipelineContext (used to build include_dirs if not provided)
        """
        names, _ = self._get_moleculetype_names_recursive(
            itp_path, set(), include_dirs=include_dirs, ctx=ctx
        )
        return names

    def _get_moleculetype_names_recursive(
        self,
        itp_path: Path,
        visited: Set[Path],
        include_dirs: Optional[List[Path]] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> Tuple[List[str], bool]:
        """
        Extract moleculetype names from an ITP, recursively following #include.
        
        Hardening v7: Now uses full GROMACS -I style include resolution instead
        of only resolving relative to current file's directory.
        
        Args:
            itp_path: Path to ITP file
            visited: Set of already-visited paths (prevents infinite loops)
            include_dirs: Optional list of include directories to search
            ctx: PipelineContext (used to build include_dirs via _build_include_search_paths)
            
        Returns:
            Tuple of (list of moleculetype names, uncertain flag)
            The uncertain flag is True if some #include could not be resolved.
        """
        if not itp_path.exists():
            return [], False
        resolved = itp_path.resolve()
        if resolved in visited:
            return [], False  # Already processed
        visited.add(resolved)
        
        # Build search paths if not provided
        if include_dirs is None and ctx is not None:
            include_dirs = self._build_include_search_paths(ctx, itp_path)
        elif include_dirs is None:
            # Fallback: only search relative to current file
            include_dirs = [itp_path.parent]
        
        names: List[str] = []
        uncertain = False
        current_section = ""
        strict_includes = bool(ctx and ctx.strict_include_resolution)
        allow_shadowing = bool(ctx and getattr(ctx, "allow_include_shadowing", False))
        content = self.safe_read_text(itp_path, strict=strict_includes)
        if not content and not strict_includes:
            uncertain = True

        for line in content.splitlines():
            stripped = line.strip()

            # Handle #include directives using parse_include_target
            if stripped.startswith("#include"):
                include_target = parse_include_target(line)
                if include_target:
                    # Resolve using search paths (GROMACS -I semantics)
                    include_path = self._resolve_include_in_dirs(
                        include_target,
                        itp_path,
                        include_dirs,
                        strict=strict_includes,
                        check_shadowing=True,
                        allow_shadowing=allow_shadowing,
                        ctx=ctx,
                    )
                    if include_path is not None:
                        inc_names, inc_uncertain = self._get_moleculetype_names_recursive(
                            include_path, visited,
                            include_dirs=include_dirs, ctx=ctx
                        )
                        names.extend(inc_names)
                        uncertain = uncertain or inc_uncertain
                    else:
                        # Could not resolve include
                        uncertain = True
                continue

            if not stripped or stripped.startswith(";"):
                continue
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                continue
            if current_section == "moleculetype":
                data = stripped.split(";")[0].strip()
                if data:
                    parts = data.split()
                    if parts:
                        names.append(parts[0])
                        current_section = ""
        
        return names, uncertain

    @staticmethod
    def _path_within_root(path: Path, root: Path) -> bool:
        """Return True when `path` is inside `root` (or equal)."""
        try:
            path.relative_to(root)
            return True
        except ValueError:
            return False

    def _include_allowed_roots(
        self,
        including_file: Path,
        include_dirs: List[Path],
        ctx: Optional["PipelineContext"] = None,
    ) -> List[Path]:
        """
        Build include-resolution boundary roots.

        Allowed roots are:
        - including file parent
        - configured include search roots
        - optional project root (if configured)
        """
        roots: List[Path] = []
        candidates: List[Path] = [including_file.parent] + list(include_dirs)

        project_root = getattr(ctx, "project_root", None) if ctx is not None else None
        if project_root:
            candidates.append(Path(project_root))

        seen: Set[Path] = set()
        for candidate in candidates:
            try:
                resolved = candidate.resolve()
            except OSError:
                continue
            if resolved in seen:
                continue
            seen.add(resolved)
            roots.append(resolved)
        return roots

    def _include_path_allowed(
        self,
        resolved_path: Path,
        allowed_roots: List[Path],
    ) -> bool:
        """Return True if include path is under at least one allowed root."""
        return any(self._path_within_root(resolved_path, root) for root in allowed_roots)
    
    def _resolve_include_in_dirs(
        self,
        include_target: str,
        including_file: Path,
        include_dirs: List[Path],
        *,
        strict: bool = False,
        check_shadowing: bool = False,
        allow_shadowing: bool = False,
        diagnostics: Optional[List[Any]] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> Optional[Path]:
        """
        Resolve an include target using GROMACS -I style search order.
        
        Shared helper for include resolution with optional shadowing checks.
        
        Search order:
        1. Including file's parent directory
        2. Each directory in include_dirs (in order)
        
        Args:
            include_target: The target path from #include directive
            including_file: Path of the file containing the #include
            include_dirs: List of directories to search
            
        Returns:
            Resolved absolute path, or None if not found
        """
        # Build search order: including file's parent first, then include_dirs
        search_order = [including_file.parent]
        for d in include_dirs:
            if d not in search_order and d.exists():
                search_order.append(d)

        allow_escape = bool(
            getattr(ctx, "allow_unsafe_include_escape", DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE)
            if ctx is not None
            else DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE
        )
        allowed_roots = self._include_allowed_roots(
            including_file=including_file,
            include_dirs=search_order,
            ctx=ctx,
        )

        attempted: List[str] = []
        found_paths: List[Path] = []
        chosen_from: Optional[Path] = None
        target_path = Path(include_target)

        if target_path.is_absolute():
            attempted.append(str(target_path))
            if target_path.exists():
                resolved = target_path.resolve()
                if not allow_escape and not self._include_path_allowed(resolved, allowed_roots):
                    roots_text = ", ".join(str(r) for r in allowed_roots) or "(none)"
                    msg = (
                        f"Blocked include target '{include_target}' in {including_file.name}: "
                        f"resolved path escapes allowed roots.\n"
                        f"  Resolved: {resolved}\n"
                        f"  Allowed roots: {roots_text}\n"
                        "Use --allow-unsafe-include-escape to override."
                    )
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                else:
                    found_paths.append(resolved)
                    chosen_from = resolved.parent
        else:
            for search_dir in search_order:
                candidate = (search_dir / include_target).resolve()
                attempted.append(str(candidate))
                if candidate.exists():
                    if not allow_escape and not self._include_path_allowed(candidate, allowed_roots):
                        roots_text = ", ".join(str(r) for r in allowed_roots) or "(none)"
                        msg = (
                            f"Blocked include escape '{include_target}' from {including_file.name}:\n"
                            f"  Candidate: {candidate}\n"
                            f"  Allowed roots: {roots_text}\n"
                            "Use --allow-unsafe-include-escape to override."
                        )
                        if strict:
                            raise SanitizerError(msg)
                        print(f"  [WARN] {msg}")
                        continue
                    found_paths.append(candidate)
                    if chosen_from is None:
                        chosen_from = search_dir.resolve()
                    if not check_shadowing:
                        break

        # Deduplicate while preserving order
        unique_found: List[Path] = []
        seen: Set[str] = set()
        for fp in found_paths:
            key = fp.as_posix()
            if key in seen:
                continue
            seen.add(key)
            unique_found.append(fp)

        resolved_path = unique_found[0] if unique_found else None
        shadowed = unique_found[1:] if len(unique_found) > 1 else []
        if check_shadowing and shadowed:
            shadowed_str = ", ".join(str(p) for p in shadowed)
            msg = (
                f"Include shadowing detected for '{include_target}' in {including_file.name}:\n"
                f"  Selected: {resolved_path}\n"
                f"  Shadowed: {shadowed_str}"
            )
            if strict and not allow_shadowing:
                raise SanitizerError(
                    f"Include shadowing violation (strict mode): {msg}\n"
                    "Use --allow-include-shadowing to bypass."
                )
            print(f"  [WARN] {msg}")

        if diagnostics is not None:
            diagnostics.append(
                IncludeResolution(
                    including_file=str(including_file),
                    include_target=include_target,
                    resolved_path=str(resolved_path) if resolved_path is not None else None,
                    resolved_from_dir=str(chosen_from) if chosen_from is not None else None,
                    attempted_paths=attempted,
                    shadowed_candidates=[str(p) for p in shadowed],
                )
            )

        if resolved_path is None and strict:
            searched = ", ".join(str(p) for p in search_order) or "(no include dirs)"
            raise SanitizerError(
                f"Unresolved #include '{include_target}' in {including_file}\n"
                f"  Searched: {searched}"
            )
        if resolved_path is None and not strict:
            searched = ", ".join(str(p) for p in search_order) or "(no include dirs)"
            print(
                f"  [WARN] Unresolved #include '{include_target}' in {including_file}\n"
                f"         Searched: {searched}"
            )

        return resolved_path
    
    def _expand_itp_include_closure(
        self,
        itp_paths: List[Path],
        ctx: "PipelineContext",
        class_map: Dict[str, str],
        ff_map: Dict[str, str],
    ) -> Tuple[List[Path], Dict[str, str], Dict[str, str]]:
        """
        BFS-expand includes for all ITPs to discover reachable files.
        
        Hardening v7: Ensures all included files are discovered, not just
        those in the immediate ITP list. Propagates class/FF maps to
        included files from their parent.
        
        Args:
            itp_paths: Initial list of ITP paths
            ctx: PipelineContext for include path building
            class_map: Dict of path_str -> classification (forcefield/staged/htpolynet/molecule)
            ff_map: Dict of path_str -> forcefield name
            
        Returns:
            (expanded_paths, updated_class_map, updated_ff_map)
        """
        visited: Set[Path] = set()
        expanded: List[Path] = []
        queue = deque(itp_paths)
        strict_includes = getattr(ctx, "strict_include_resolution", False)
        allow_shadowing = getattr(ctx, "allow_include_shadowing", False)
        
        while queue:
            itp = queue.popleft()
            resolved = itp.resolve()
            if resolved in visited:
                continue
            visited.add(resolved)
            expanded.append(resolved)

            # Normalize maps to resolved paths for deterministic downstream lookups.
            parent_class = (
                class_map.get(str(itp))
                or class_map.get(str(resolved))
                or "unknown"
            )
            parent_ff = (
                ff_map.get(str(itp))
                or ff_map.get(str(resolved))
                or "UNKNOWN"
            )
            class_map[str(resolved)] = parent_class
            ff_map[str(resolved)] = parent_ff
            
            if not itp.exists():
                continue
            
            # Build search paths for this file
            include_dirs = self._build_include_search_paths(ctx, itp)
            
            # Scan for includes
            content = self.safe_read_text(itp, strict=strict_includes)
            if not content:
                continue
            
            for line in content.splitlines():
                if not line.strip().startswith("#include"):
                    continue
                include_target = parse_include_target(line)
                if not include_target:
                    continue
                
                inc_path = self._resolve_include_in_dirs(
                    include_target,
                    itp,
                    include_dirs,
                    strict=strict_includes,
                    check_shadowing=True,
                    allow_shadowing=allow_shadowing,
                    diagnostics=self.diagnostics,
                    ctx=ctx,
                )
                if inc_path is None:
                    continue
                
                inc_resolved = inc_path.resolve()
                if inc_resolved in visited:
                    continue
                
                # Only add if not already tracked
                if str(inc_resolved) not in ff_map:
                    # Inherit class/ff from parent
                    class_map[str(inc_resolved)] = parent_class
                    ff_map[str(inc_resolved)] = parent_ff
                    print(f"    + Include ({parent_class}): {inc_resolved.name}")
                
                queue.append(inc_resolved)
        
        return expanded, class_map, ff_map
    
    def _moleculetype_signatures_differ(
        self,
        itp1: Path,
        itp2: Path,
        mol_name: str,
        ctx: Optional["PipelineContext"] = None,
    ) -> bool:
        """
        Check if moleculetype definitions differ between two ITPs.
        
        Uses full-table streaming signature from _compute_moleculetype_signature():
        - SHA256 over canonicalized [ atoms ] entries (type|name|charge@1e-4)
        - Section presence flags and cheap charge distribution stats
        
        Args:
            itp1: First ITP path
            itp2: Second ITP path  
            mol_name: Moleculetype name to compare
            
        Returns:
            True if signatures differ, False if they match
        """
        sig1 = self._compute_moleculetype_signature(itp1, mol_name, ctx=ctx)
        sig2 = self._compute_moleculetype_signature(itp2, mol_name, ctx=ctx)
        
        if sig1 is None or sig2 is None:
            # Could not compute signature - assume different to be safe
            return True
        
        return sig1 != sig2
    
    def _classify_ionic_moleculetypes(
        self,
        molecule_itp_paths: Dict[str, Path],
        ctx: Optional["PipelineContext"] = None,
    ) -> Dict[str, str]:
        """
        Classify moleculetypes for charge correction protection.
        
        Hardening v6: Extended to protect solvents AND ions.
        
        Classification rules (priority order):
        1. Pattern match against PROTECTED_ION_PATTERNS -> "protected_ion"
        2. Pattern match against PROTECTED_SOLVENT_PATTERNS -> "protected_solvent"
        3. n_atoms == 1 and |q| > 0 -> "monatomic_ion"
        4. abs(round(q_total)) >= 1 and close to integer -> "ionic"
        5. Name contains polymer/chain patterns -> "protected_polymer"
        6. Otherwise -> "neutral"
        
        Args:
            molecule_itp_paths: Dict of molecule_name -> itp_path
            ctx: PipelineContext for tolerance
            
        Returns:
            Dict of molecule_name -> classification
        """
        from ..charge_neutrality import ChargeParser
        
        tol = ctx.charge_neutrality_tol if ctx else 1e-6
        classifications: Dict[str, str] = {}
        
        for mol_name, itp_path in molecule_itp_paths.items():
            # Hardening v7: Use boundary-aware pattern matching
            # Priority 1: Check ion patterns
            if _matches_protected_pattern(mol_name, PROTECTED_ION_PATTERNS):
                classifications[mol_name] = "protected_ion"
                continue
            
            # Priority 2: Check solvent patterns
            if _matches_protected_pattern(mol_name, PROTECTED_SOLVENT_PATTERNS):
                classifications[mol_name] = "protected_solvent"
                continue
            
            # Priority 3-5: Charge-based classification
            try:
                per_mol = ChargeParser.parse_atoms_charges_per_molecule(itp_path)
                resolved_mol_name = mol_name if mol_name in per_mol else None
                if resolved_mol_name is None:
                    matches = [k for k in per_mol if k.casefold() == mol_name.casefold()]
                    if len(matches) == 1:
                        resolved_mol_name = matches[0]
                    elif len(matches) == 0 and len(per_mol) == 1:
                        resolved_mol_name = next(iter(per_mol.keys()))
                if resolved_mol_name is None:
                    # Check if this looks like a polymer/fragment
                    name_lower = mol_name.lower()
                    if any(p in name_lower for p in ('polymer', 'chain', 'frag', 'oligo', 'peo', 'peg')):
                        classifications[mol_name] = "protected_polymer"
                        continue
                    classifications[mol_name] = "unknown"
                    continue
                
                charges = per_mol[resolved_mol_name][0]
                n_atoms = len(charges)
                q_total = sum(charges)
                
                # Monatomic ion check
                if n_atoms == 1 and abs(q_total) > tol:
                    classifications[mol_name] = "monatomic_ion"
                    continue
                
                # Integer charge check
                rounded_q = round(q_total)
                if abs(rounded_q) >= 1 and abs(q_total - rounded_q) < tol:
                    classifications[mol_name] = "ionic"
                    continue
                
                # Check for polymer patterns in name
                name_lower = mol_name.lower()
                if any(p in name_lower for p in ('polymer', 'chain', 'frag', 'oligo', 'peo', 'peg')):
                    classifications[mol_name] = "protected_polymer"
                    continue
                
                # Default: neutral
                classifications[mol_name] = "neutral"
                    
            except Exception:
                classifications[mol_name] = "unknown"
        
        # Log classification summary
        prot_ions = [k for k, v in classifications.items() if v in ("protected_ion", "ionic", "monatomic_ion")]
        prot_solvents = [k for k, v in classifications.items() if v == "protected_solvent"]
        prot_polymers = [k for k, v in classifications.items() if v == "protected_polymer"]
        neutral = [k for k, v in classifications.items() if v == "neutral"]
        
        if prot_ions:
            print(f"  [CHARGE-FIX] Protected ions: {prot_ions}")
        if prot_solvents:
            print(f"  [CHARGE-FIX] Protected solvents: {prot_solvents}")
        if prot_polymers:
            print(f"  [CHARGE-FIX] Protected polymers: {prot_polymers}")
        if neutral:
            print(f"  [CHARGE-FIX] Correction candidates (neutral): {neutral}")
        
        return classifications

    def _warn_electrolyte_fixed_charge_limitations(
        self,
        correction_result,
        classifications: Dict[str, str],
        ctx: "PipelineContext",
    ) -> Optional[Dict[str, Any]]:
        """
        Emit targeted warnings for concentrated-electrolyte charge-fix edge cases.
        """
        target = correction_result.target_molecule
        target_cf = target.casefold() if target else None
        class_by_casefold = {name.casefold(): cls for name, cls in classifications.items()}
        target_class = class_by_casefold.get(target_cf or "", "unknown")

        protected_classes = {"protected_ion", "ionic", "monatomic_ion", "protected_solvent"}
        protected_detected = any(cls in protected_classes for cls in classifications.values())
        correction_refused = bool(getattr(correction_result, "correction_refused", False))
        forced_polymer_target = bool(
            correction_result.corrected and target_class == "protected_polymer"
        )
        if not ((protected_detected and correction_refused) or forced_polymer_target):
            return None

        advisory = (
            "Concentrated-electrolyte caution: fixed-charge non-polarizable force fields "
            "omit electronic polarizability and can over-bind ions."
        )
        if correction_refused and protected_detected:
            print(f"  [WARN] {advisory}")
            print(
                "  [WARN] Charge correction was refused due to protected ion/solvent policy. "
                "This often indicates model limitations, not just automation constraints."
            )
        if forced_polymer_target:
            print(f"  [WARN] {advisory}")
            print(
                "  [WARN] Correction pressure shifted onto polymer/backbone atoms while ions/solvents "
                "remained protected; review physical validity before production."
            )
        print(
            "  [WARN] Consider mitigation: ECC-style charge scaling or a polarizable force field model."
        )
        print(
            "  [WARN] Recommended safer pipeline settings for GPE automation: "
            "--charge-fix-polymer-exclusion-bonds 3, --polymer-net-charge-tol 0.05, "
            "and keep --strict-gro-top-check enabled."
        )

        return {
            "warning_emitted": True,
            "target_molecule": target,
            "target_classification": target_class,
            "correction_refused": correction_refused,
            "forced_polymer_target": forced_polymer_target,
            "fixed_charge_polarizability_caveat": True,
            "recommended_mitigations": [
                "ECC-style charge scaling",
                "polarizable force field models",
            ],
            "recommended_flags": {
                "charge_fix_polymer_exclusion_bonds": 3,
                "polymer_net_charge_tol": 0.05,
                "strict_gro_top_check": True,
            },
        }
    
    def _verify_ion_protection(
        self,
        correction_result,
        classifications: Dict[str, str],
        ctx: "PipelineContext",
    ) -> Dict[str, Any]:
        """
        Verify that charge correction obeys protected-species policy.

        Returns structured audit details for logs/manifest, including violations.
        """
        violations: List[str] = []

        # Get config
        allow_ions = ctx.charge_fix_allow_ions
        allow_solvents = getattr(ctx, "charge_fix_allow_solvents", False)
        max_delta = ctx.charge_fix_max_delta_per_atom if ctx.charge_fix_max_delta_per_atom is not None else 1e-5
        target_allowlist = self._parse_csv_casefold_set(ctx.charge_fix_target_allowlist)
        strict_mode = bool(getattr(ctx, "strict_charge_neutrality", False))
        polymer_mode_cfg = str(getattr(ctx, "charge_fix_polymer_method", "skip_if_small"))
        
        target_mol = correction_result.target_molecule
        target_mol_cf = target_mol.casefold() if target_mol else None
        class_by_casefold = {name.casefold(): cls for name, cls in classifications.items()}
        target_classification = class_by_casefold.get(target_mol_cf or "", "unknown")
        target_from_allowlist = bool(
            target_allowlist and target_mol_cf and target_mol_cf in target_allowlist
        )

        protected_molecules = {
            name: cls
            for name, cls in classifications.items()
            if cls in {
                "protected_ion",
                "ionic",
                "monatomic_ion",
                "protected_solvent",
                "protected_polymer",
            }
        }

        # Protected classification categories
        PROTECTED_CATEGORIES = {
            "protected_ion": (allow_ions, "--charge-fix-allow-ions"),
            "ionic": (allow_ions, "--charge-fix-allow-ions"),
            "monatomic_ion": (allow_ions, "--charge-fix-allow-ions"),
            "protected_solvent": (allow_solvents, "--charge-fix-allow-solvents"),
            "protected_polymer": (
                False,
                "--charge-fix-target-allowlist + --charge-fix-polymer-method spread_safe",
            ),  # Never auto-allow
        }

        if strict_mode and protected_molecules and target_allowlist is None:
            violations.append(
                "Strict mode with protected molecules requires explicit "
                "--charge-fix-target-allowlist."
            )

        # Check if target is protected
        if target_mol:
            classification = target_classification
            neutrality_status = str(getattr(correction_result, "neutrality_status", "") or "")
            polymer_skip_status = neutrality_status in {
                "polymer_acceptable_drift",
                "polymer_skipped_non_strict",
            }

            # If target is in explicit allowlist, it's OK
            if target_from_allowlist:
                print(f"  [CHARGE-FIX] Target '{target_mol}' ({classification}) is in explicit allowlist")
            elif classification in PROTECTED_CATEGORIES:
                is_allowed, flag_hint = PROTECTED_CATEGORIES[classification]
                if (
                    classification == "protected_polymer"
                    and not correction_result.corrected
                    and polymer_skip_status
                ):
                    is_allowed = True
                if not is_allowed:
                    violations.append(
                        f"Charge correction targeted protected species '{target_mol}' "
                        f"(classification={classification}). "
                        f"Use {flag_hint} to allow, or add to --charge-fix-target-allowlist."
                    )

        polymer_outcome = "not_applicable"
        effective_method = str(getattr(correction_result, "effective_method", "") or "")
        if target_mol and target_classification == "protected_polymer":
            if correction_result.corrected:
                polymer_outcome = "allowed"
                if not target_from_allowlist:
                    polymer_outcome = "refused"
                    violations.append(
                        f"Protected polymer target '{target_mol}' must be explicitly allowlisted."
                    )
                if polymer_mode_cfg != "spread_safe":
                    polymer_outcome = "refused"
                    violations.append(
                        f"Protected polymer target '{target_mol}' requires "
                        "--charge-fix-polymer-method spread_safe."
                    )
                if effective_method and effective_method != "spread_safe":
                    polymer_outcome = "refused"
                    violations.append(
                        f"Protected polymer target '{target_mol}' was corrected with "
                        f"unsafe method '{effective_method}' (expected spread_safe)."
                    )
            else:
                polymer_outcome = "refused" if correction_result.correction_refused else "skipped"

        # Check max delta per atom
        if correction_result.max_delta_applied is not None:
            if correction_result.max_delta_applied > max_delta:
                violations.append(
                    f"Max per-atom charge delta ({correction_result.max_delta_applied:.3e}) "
                    f"exceeds limit ({max_delta:.3e})"
                )
        
        # Log summary
        if correction_result.corrected:
            print(f"  [CHARGE-FIX SUMMARY]")
            print(f"    Target: {target_mol or 'unknown'}")
            print(f"    Classification: {class_by_casefold.get(target_mol_cf or '', 'unknown')}")
            print(f"    Method: {correction_result.correction_method}")
            print(f"    Atoms adjusted: {correction_result.atoms_adjusted}")
            if correction_result.max_delta_applied is not None:
                print(f"    Max delta: {correction_result.max_delta_applied:.3e}")
            if correction_result.total_delta is not None:
                print(f"    Total delta: {correction_result.total_delta:.3e}")
            classification = target_classification
            if classification in PROTECTED_CATEGORIES:
                print(f"    [WARN] Target '{target_mol}' classified as {classification}")

        audit: Dict[str, Any] = {
            "target_molecule": target_mol,
            "target_classification": target_classification,
            "all_classifications": classifications,
            "protected_molecules": protected_molecules,
            "strict_mode_enforced": strict_mode,
            "explicit_allowlist_provided": target_allowlist is not None,
            "target_from_allowlist": target_from_allowlist,
            "corrected": bool(correction_result.corrected),
            "correction_refused": bool(getattr(correction_result, "correction_refused", False)),
            "neutrality_status": getattr(correction_result, "neutrality_status", None),
            "polymer_correction_outcome": polymer_outcome,
            "polymer_safe_mode_configured": polymer_mode_cfg == "spread_safe",
            "polymer_safe_mode_effective": effective_method == "spread_safe",
            "effective_method": effective_method or None,
            "violations": violations,
            "config": {
                "allow_ions": allow_ions,
                "allow_solvents": allow_solvents,
                "max_delta_per_atom": max_delta,
                "target_allowlist": list(target_allowlist) if target_allowlist else None,
                "polymer_net_charge_tol": getattr(ctx, "polymer_net_charge_tol", 1e-3),
                "charge_fix_polymer_method": polymer_mode_cfg,
                "charge_fix_polymer_exclusion_bonds": getattr(ctx, "charge_fix_polymer_exclusion_bonds", 2),
            },
            "polymer_policy": getattr(correction_result, "polymer_policy", None),
        }

        # Record audit to manifest
        if ctx.manifest:
            ctx.manifest.set_sanitizer_output("charge_protection_audit", audit)

        return audit
    
    def _compute_moleculetype_signature(
        self,
        itp_path: Path,
        mol_name: str,
        ctx: Optional["PipelineContext"] = None,
    ) -> Optional[str]:
        """
        Compute a collision-resistant signature for one moleculetype.

        Uses streaming SHA-256 over canonicalized content for topology-relevant
        sections (atoms/bonds/angles/dihedrals/...); this catches materially
        different same-name definitions even when atoms-only signatures match.
        """
        if not itp_path.exists():
            return None

        strict_reads = bool(ctx and ctx.strict_include_resolution)
        content = self.safe_read_text(itp_path, strict=strict_reads)
        if not content:
            return None
        
        current_section = ""
        current_mol = None
        in_target_mol = False
        
        atom_count = 0
        sections_present: Set[str] = set()
        section_line_counts: Dict[str, int] = {}
        qmin: Optional[float] = None
        qmax: Optional[float] = None
        sum_abs_q = 0.0
        sum_q2 = 0.0
        hasher = hashlib.sha256()
        quant_step = 1e-4
        topology_sections = {
            "atoms",
            "bonds",
            "pairs",
            "angles",
            "dihedrals",
            "impropers",
            "constraints",
            "settles",
            "virtual_sites2",
            "virtual_sites3",
            "virtual_sites4",
            "virtual_sitesn",
            "exclusions",
            "position_restraints",
            "distance_restraints",
            "dihedral_restraints",
            "angle_restraints",
            "orientation_restraints",
        }
        
        for line in content.splitlines():
            stripped = line.strip()
            
            # Track section
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if current_section == "moleculetype":
                    current_mol = None
                    in_target_mol = False
                continue
            
            # Skip comments/blanks
            data = strip_gmx_comment(stripped).strip()
            if not data:
                continue
            
            # Track moleculetype name
            if current_section == "moleculetype":
                parts = data.split()
                if parts:
                    current_mol = parts[0]
                    in_target_mol = (current_mol == mol_name)
                continue

            if not in_target_mol:
                continue

            if data.startswith("#"):
                preproc_key = f"{current_section}:preproc"
                sections_present.add(preproc_key)
                section_line_counts[preproc_key] = section_line_counts.get(preproc_key, 0) + 1
                hasher.update(f"{preproc_key}|{data}\n".encode("utf-8"))
                continue

            if current_section in topology_sections:
                canonical_tokens = " ".join(data.split())
                sections_present.add(current_section)
                section_line_counts[current_section] = section_line_counts.get(current_section, 0) + 1
                hasher.update(f"{current_section}|{canonical_tokens}\n".encode("utf-8"))
            
            # Hash all atom entries in target molecule
            if current_section == "atoms":
                parts = data.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    atom_count += 1
                    charge_idx = 6 if len(parts) >= 8 else (5 if len(parts) >= 7 else None)
                    if charge_idx is None:
                        return None
                    try:
                        raw_q = float(parts[charge_idx])
                    except (ValueError, TypeError):
                        return None

                    qmin = raw_q if qmin is None else min(qmin, raw_q)
                    qmax = raw_q if qmax is None else max(qmax, raw_q)
                    sum_abs_q += abs(raw_q)
                    sum_q2 += raw_q * raw_q
                    quant_q = round(raw_q / quant_step) * quant_step
                    canonical = f"{parts[1]}|{parts[4]}|{quant_q:+.4f}\n"
                    hasher.update(canonical.encode("utf-8"))
        
        if atom_count == 0:
            return None

        qmin_val = qmin if qmin is not None else 0.0
        qmax_val = qmax if qmax is not None else 0.0
        digest = hasher.hexdigest()

        return (
            f"atoms={atom_count};"
            f"digest={digest};"
            f"qmin={qmin_val:+.6f};"
            f"qmax={qmax_val:+.6f};"
            f"sum_abs_q={sum_abs_q:.6f};"
            f"sum_q2={sum_q2:.6f};"
            f"sections={','.join(sorted(sections_present))};"
            f"section_counts={','.join(f'{k}:{section_line_counts[k]}' for k in sorted(section_line_counts))}"
        )

    def _get_ordered_molecules_from_top(self, top_path: Path) -> List[Tuple[str, int]]:
        """
        Parse [molecules] section preserving order, collapsing duplicates.
        
        Task D: GROMACS allows the same molecule name to appear on multiple
        lines in [molecules]. We collapse duplicates by summing counts into
        the first occurrence to avoid breaking include ordering.
        """
        if not top_path.exists():
            return []
        # Use list for order + dict for fast lookup of first occurrence index
        ordered: List[Tuple[str, int]] = []
        seen_indices: Dict[str, int] = {}  # name -> index in ordered list
        in_molecules = False
        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return []
        for line in content.splitlines():
            stripped = line.strip()
            section = parse_section_name(line)
            if section is not None:
                in_molecules = (section == "molecules")
                continue
            if not in_molecules:
                continue
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            data = stripped.split(";")[0].strip()
            if not data:
                continue
            parts = data.split()
            if len(parts) >= 2:
                try:
                    name = parts[0]
                    count = int(parts[1])
                    if name in seen_indices:
                        # Collapse: add count to existing entry
                        idx = seen_indices[name]
                        old_name, old_count = ordered[idx]
                        ordered[idx] = (old_name, old_count + count)
                    else:
                        # First occurrence: add new entry
                        seen_indices[name] = len(ordered)
                        ordered.append((name, count))
                except ValueError:
                    continue
        return ordered
    
    def _get_molecules_from_htpolynet_top(self, ctx: "PipelineContext") -> Optional[Dict[str, int]]:
        """
        Parse [molecules] section from HTPolyNet-staged system.top.
        
        Bug #2: After HTPolyNet crosslinks, the authoritative [molecules]
        section comes from the HTPolyNet-generated topology, not from
        the pre-reaction manifest counts.
        
        Returns:
            Dict of molecule_name -> count, or None if no HTPolyNet top exists
        """
        top_path = ctx.get_input_path(
            "systems", ctx.system_id, "gromacs", "top", "system.top"
        )
        
        if not top_path.exists():
            return None

        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return None
        
        molecules: Dict[str, int] = {}
        in_molecules = False
        
        for line in content.splitlines():
            stripped = line.strip()
            
            # Detect section headers
            section = parse_section_name(line)
            if section is not None:
                in_molecules = (section == "molecules")
                continue
            
            if in_molecules:
                # Skip comments and empty lines
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                
                parts = data.split()
                if len(parts) >= 2:
                    try:
                        mol_name = parts[0]
                        mol_count = int(parts[1])
                        # Accumulate in case same molecule appears multiple times
                        molecules[mol_name] = molecules.get(mol_name, 0) + mol_count
                    except ValueError:
                        continue
        
        return molecules if molecules else None

    def _map_sanitized_itps_to_molecules(
        self,
        sanitized_current_dir: Path,
        needed_molecule_names: set,
    ) -> Dict[str, Path]:
        """Map molecule names to sanitized ITPs (by moleculetype name)."""
        if not sanitized_current_dir.exists() or not needed_molecule_names:
            return {}
        name_map = {n.casefold(): n for n in needed_molecule_names}
        mapping: Dict[str, Path] = {}
        for itp_path in self._sorted_paths_casefold(sanitized_current_dir.glob("*.itp")):
            for mol_name in self._get_moleculetype_names(itp_path):
                key = mol_name.casefold()
                if key not in name_map:
                    continue
                canonical = name_map[key]
                existing = mapping.get(canonical)
                if existing and existing.resolve() != itp_path.resolve():
                    raise SanitizerError(
                        f"Multiple sanitized ITPs define moleculetype '{canonical}':\n"
                        f"  - {existing}\n  - {itp_path}"
                    )
                mapping[canonical] = itp_path
        missing = [n for n in needed_molecule_names if n not in mapping]
        if missing:
            raise SanitizerError(
                f"Missing sanitized ITPs for molecules: {sorted(missing, key=str.casefold)}"
            )
        return mapping

    def _ordered_molecule_itp_paths(
        self,
        ordered_molecules: List[Tuple[str, int]],
        molecule_itp_paths: Dict[str, Path],
        sanitized_current_dir: Path,
    ) -> List[Path]:
        """
        Return molecule ITP paths ordered by [molecules] section.
        
        Task D: Ensure each ITP is included at most once, preserving first-seen
        order. This handles cases where the same molecule appears multiple times
        in [molecules] (which is valid GROMACS syntax).
        """
        ordered_paths: List[Path] = []
        seen_paths: Set[Path] = set()  # Track already-included ITPs
        for name, count in ordered_molecules:
            if count <= 0:
                continue
            itp_path = molecule_itp_paths.get(name)
            if itp_path is None:
                raise SanitizerError(
                    f"No sanitized ITP found for molecule '{name}' needed by [molecules]."
                )
            resolved = itp_path.resolve()
            if resolved not in seen_paths:
                ordered_paths.append(itp_path)
                seen_paths.add(resolved)
        # Relaxed assertion: only check for case-insensitive basename collisions
        # among distinct paths (not duplicates of the same path)
        self._assert_no_case_insensitive_duplicates(ordered_paths, "sanitized molecule ITPs")
        return ordered_paths

    def _validate_system_top_includes(
        self, system_top_path: Path, ctx: Optional["PipelineContext"] = None
    ) -> None:
        """Ensure all #include targets resolve on disk."""
        if not system_top_path.exists():
            return
        strict_includes = getattr(ctx, "strict_include_resolution", False) if ctx else False
        include_dirs = [system_top_path.parent]
        if ctx is not None:
            include_dirs = self._build_include_search_paths(ctx, system_top_path)
        content = self.safe_read_text(system_top_path, strict=strict_includes)
        if not content:
            return
        for line in content.splitlines():
            include_target = parse_include_target(line)
            if not include_target:
                continue
            include_path = self._resolve_include_in_dirs(
                include_target,
                system_top_path,
                include_dirs,
                strict=strict_includes,
                check_shadowing=True,
                allow_shadowing=getattr(ctx, "allow_include_shadowing", False) if ctx else False,
                ctx=ctx,
            )
            if include_path is None or not include_path.exists():
                if strict_includes:
                    searched = ", ".join(str(p) for p in include_dirs) or "(no include dirs)"
                    raise SanitizerError(
                        f"system.top includes missing file: {include_target}\n"
                        f"  Searched: {searched}"
                    )

    def _parse_lj_atomtypes_row(
        self,
        data: str,
        source_name: str,
        line_no: int,
    ) -> Tuple[str, str, float, float, Optional[str]]:
        """
        Parse one [ atomtypes ] data row for LJ validation.

        Supports 6/7/8-column styles and shifted optional fields by:
        - taking LJ params from the final two numeric tokens
        - resolving ptype from the right-most {A,D,V} token before LJ params
        """
        parts = data.split()
        if len(parts) < 5:
            raise ValueError(f"Too few fields for [ atomtypes ] row: '{data}'")

        atomtype_name = parts[0]
        try:
            p1 = float(parts[-2])
            p2 = float(parts[-1])
        except ValueError as exc:
            raise ValueError(
                f"Non-numeric LJ params in atomtype '{atomtype_name}': "
                f"'{parts[-2]}', '{parts[-1]}'"
            ) from exc

        ptype_positions = [
            idx for idx, token in enumerate(parts[:-2]) if token.upper() in PTYPE_TOKENS
        ]
        warning_msg: Optional[str] = None
        if not ptype_positions:
            ptype = "A"
            warning_msg = (
                f"Ambiguous ptype in [ atomtypes ] row at {source_name}:{line_no}; "
                f"defaulting to ptype='A'. Raw: {data}"
            )
        else:
            chosen_idx = ptype_positions[-1]  # closest token to LJ params wins
            ptype = parts[chosen_idx].upper()
            if len(ptype_positions) > 1:
                warning_msg = (
                    f"Multiple ptype tokens in [ atomtypes ] row at {source_name}:{line_no}; "
                    f"using right-most token '{ptype}'. Raw: {data}"
                )
            elif chosen_idx != len(parts) - 3:
                warning_msg = (
                    f"Shifted ptype column detected at {source_name}:{line_no}; "
                    f"interpreting token '{ptype}' as ptype. Raw: {data}"
                )

        return atomtype_name, ptype, p1, p2, warning_msg

    def _warn_atomtype_outliers(
        self, combined_atomtypes_path: Path, ctx: Optional["PipelineContext"] = None,
        comb_rule: Optional[int] = None
    ) -> List[str]:
        """
        Validate LJ outliers with profile+policy controls.

        - comb-rule 2/3: sigma/epsilon profile bounds
        - comb-rule 1: robust statistical outlier check
        - unit-screaming thresholds always raise SanitizerError
        """
        warnings: List[str] = []
        if not combined_atomtypes_path.exists():
            return warnings

        policy = getattr(ctx, "lj_outlier_policy", "warn") if ctx else "warn"
        if getattr(ctx, "strict_lj_validation", False) and policy == "warn":
            # Backward-compatibility with legacy strict flag.
            policy = "error"

        profile = getattr(ctx, "lj_bounds_profile", "aa") if ctx else "aa"
        if profile == "custom":
            sigma_max = float(getattr(ctx, "lj_sigma_max_nm", 1.0))
            epsilon_max = float(getattr(ctx, "lj_epsilon_max_kj_mol", 50.0))
            sigma_min = 0.0
        else:
            profile_thresholds = LJ_BOUNDS_PROFILES.get(profile, LJ_BOUNDS_PROFILES["aa"])
            sigma_min = float(profile_thresholds.get("sigma_min_nm", 0.0))
            sigma_max = float(profile_thresholds.get("sigma_max_nm", 1.0))
            epsilon_max = float(profile_thresholds.get("epsilon_max_kj_mol", 50.0))

        # Backward-compatible threshold override hook.
        ctx_thresholds = getattr(ctx, "lj_outlier_thresholds", None) if ctx else None
        if isinstance(ctx_thresholds, dict):
            if "sigma_min_nm" in ctx_thresholds:
                sigma_min = float(ctx_thresholds["sigma_min_nm"])
            if "sigma_max_nm" in ctx_thresholds:
                sigma_max = float(ctx_thresholds["sigma_max_nm"])
            if "epsilon_max_kj_mol" in ctx_thresholds:
                epsilon_max = float(ctx_thresholds["epsilon_max_kj_mol"])

        atomtype_data: List[Dict[str, Any]] = []
        in_atomtypes = False

        content = self.safe_read_text(combined_atomtypes_path, strict=False)
        if not content:
            return warnings

        for line_no, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            section = parse_section_name(line)
            if section is not None:
                in_atomtypes = (section == "atomtypes")
                continue
            if not in_atomtypes:
                continue
            data = strip_gmx_comment(stripped).strip()
            if not data:
                continue
            parts = data.split()
            if len(parts) < 5:
                continue
            try:
                atomtype_name, ptype, p1, p2, parse_warning = self._parse_lj_atomtypes_row(
                    data,
                    combined_atomtypes_path.name,
                    line_no,
                )
            except ValueError:
                msg = (
                    f"Failed to parse [ atomtypes ] LJ row in {combined_atomtypes_path.name}:{line_no}: "
                    f"{data}"
                )
                print(f"  [WARN] {msg}")
                warnings.append(msg)
                continue
            if parse_warning:
                print(f"  [WARN] {parse_warning}")
                warnings.append(parse_warning)

            comment = stripped.split(";", 1)[1].strip() if ";" in stripped else ""
            provenance = None
            if comment.lower().startswith("from:"):
                provenance = comment[5:].strip()

            atomtype_data.append(
                {
                    "name": atomtype_name,
                    "ptype": ptype,
                    "p1": p1,
                    "p2": p2,
                    "provenance": provenance,
                    "raw": data,
                }
            )

        if not atomtype_data:
            return warnings

        report: Dict[str, Any] = {
            "policy": policy,
            "profile": profile,
            "comb_rule": comb_rule,
            "thresholds": {
                "sigma_min_nm": sigma_min,
                "sigma_max_nm": sigma_max,
                "epsilon_max_kj_mol": epsilon_max,
            },
            "unit_screaming_thresholds": dict(LJ_UNIT_SCREAMING_THRESHOLDS),
            "top_offenders": [],
            "warning_count": 0,
            "unit_error_count": 0,
        }

        # comb-rule=1: retain robust statistical fallback.
        if comb_rule == 1:
            print(
                "  [INFO] comb-rule=1: treating last two fields as raw LJ params "
                "(p1/p2), not sigma/epsilon"
            )
            robust_data = [
                (row["name"], float(row["p1"]), float(row["p2"]), str(row["raw"]))
                for row in atomtype_data
            ]
            warnings.extend(self._check_lj_outliers_robust(
                robust_data, combined_atomtypes_path.name
            ))
            report["warning_count"] = len(warnings)
            report["top_offenders"] = [
                {"message": msg} for msg in warnings[:10]
            ]
            if ctx and hasattr(ctx, "manifest") and ctx.manifest:
                ctx.manifest.set_sanitizer_output("lj_validation", report)
            if warnings and policy == "error":
                raise SanitizerError(
                    "LJ outlier policy=error triggered for comb-rule=1 robust outlier checks."
                )
            return warnings

        offenders: List[Dict[str, Any]] = []
        unit_errors: List[Dict[str, Any]] = []
        unit_sigma_max = float(LJ_UNIT_SCREAMING_THRESHOLDS["sigma_max_nm"])
        unit_eps_max = float(LJ_UNIT_SCREAMING_THRESHOLDS["epsilon_max_kj_mol"])

        for row in atomtype_data:
            name = str(row["name"])
            ptype = str(row.get("ptype", "A")).upper()
            sigma = float(row["p1"])
            epsilon = float(row["p2"])
            provenance = row.get("provenance")

            if ptype in {"D", "V"}:
                # Skip dummy/virtual sites to avoid false positives in modern topologies.
                continue

            if any((math.isnan(sigma), math.isnan(epsilon), math.isinf(sigma), math.isinf(epsilon))):
                msg = (
                    f"NaN/Inf LJ params in atomtype '{name}': sigma={sigma}, epsilon={epsilon} "
                    f"({provenance or combined_atomtypes_path.name})"
                )
                unit_errors.append(
                    {
                        "name": name,
                        "ptype": ptype,
                        "sigma_nm": sigma,
                        "epsilon_kj_mol": epsilon,
                        "provenance": provenance,
                        "kind": "nan_or_inf",
                        "message": msg,
                    }
                )
                continue

            if sigma > unit_sigma_max or epsilon > unit_eps_max:
                msg = (
                    f"Unit-screaming LJ value in atomtype '{name}': "
                    f"sigma={sigma:.6g} nm, epsilon={epsilon:.6g} kJ/mol "
                    f"(>{unit_sigma_max} nm or >{unit_eps_max} kJ/mol)."
                )
                unit_errors.append(
                    {
                        "name": name,
                        "ptype": ptype,
                        "sigma_nm": sigma,
                        "epsilon_kj_mol": epsilon,
                        "provenance": provenance,
                        "kind": "unit_screaming",
                        "message": msg,
                    }
                )
                continue

            if sigma < 0:
                msg = (
                    f"Negative sigma in atomtype '{name}' (sigma={sigma:.6g} nm, epsilon={epsilon:.6g} kJ/mol). "
                    "Some GROMACS forcefields use negative sigma as special encoding; review manually."
                )
                offenders.append(
                    {
                        "name": name,
                        "ptype": ptype,
                        "sigma_nm": sigma,
                        "epsilon_kj_mol": epsilon,
                        "provenance": provenance,
                        "kind": "negative_sigma_special",
                        "severity": "warn",
                        "score": 0.0,
                        "message": msg,
                    }
                )
                warnings.append(msg)
                continue

            if epsilon < 0:
                msg = (
                    f"Negative epsilon in atomtype '{name}' (epsilon={epsilon:.6g} kJ/mol) "
                    f"at {provenance or combined_atomtypes_path.name}."
                )
                offenders.append(
                    {
                        "name": name,
                        "ptype": ptype,
                        "sigma_nm": sigma,
                        "epsilon_kj_mol": epsilon,
                        "provenance": provenance,
                        "kind": "negative_epsilon",
                        "severity": "warn",
                        "score": 0.0,
                        "message": msg,
                    }
                )
                warnings.append(msg)
                continue

            if sigma == 0 and epsilon == 0:
                continue

            local_msgs: List[str] = []
            score_terms: List[float] = []
            if sigma_max > 0 and sigma > sigma_max:
                local_msgs.append(
                    f"sigma={sigma:.6g} nm > {sigma_max:.6g} nm"
                )
                score_terms.append(sigma / sigma_max)
            if epsilon_max > 0 and epsilon > epsilon_max:
                local_msgs.append(
                    f"epsilon={epsilon:.6g} kJ/mol > {epsilon_max:.6g} kJ/mol"
                )
                score_terms.append(epsilon / epsilon_max)
            if sigma_min > 0 and sigma > 0 and sigma < sigma_min:
                local_msgs.append(
                    f"sigma={sigma:.6g} nm < {sigma_min:.6g} nm"
                )
                score_terms.append((sigma_min / sigma) if sigma > 0 else 0.0)

            if local_msgs:
                msg = (
                    f"LJ profile outlier '{name}' ({provenance or combined_atomtypes_path.name}): "
                    + "; ".join(local_msgs)
                )
                offenders.append(
                    {
                        "name": name,
                        "ptype": ptype,
                        "sigma_nm": sigma,
                        "epsilon_kj_mol": epsilon,
                        "provenance": provenance,
                        "kind": "profile_outlier",
                        "severity": "warn",
                        "score": max(score_terms) if score_terms else 0.0,
                        "message": msg,
                    }
                )
                warnings.append(msg)

        offenders_sorted = sorted(
            offenders,
            key=lambda item: (
                -float(item.get("score", 0.0)),
                str(item.get("name", "")).casefold(),
                str(item.get("name", "")),
            ),
        )
        top_offenders = offenders_sorted[:10]
        report["top_offenders"] = top_offenders
        report["warning_count"] = len(warnings)
        report["unit_error_count"] = len(unit_errors)

        if top_offenders:
            print(
                "  [WARN] LJ outlier report: "
                f"profile={profile}, policy={policy}, offenders={len(offenders)}"
            )
            for item in top_offenders[:5]:
                print(f"  [WARN]   - {item.get('message')}")

        if unit_errors:
            report["top_unit_errors"] = unit_errors[:10]
            if ctx and hasattr(ctx, "manifest") and ctx.manifest:
                ctx.manifest.set_sanitizer_output("lj_validation", report)
            first = unit_errors[0]
            raise SanitizerError(
                "LJ unit-screaming threshold violated: "
                f"{first.get('message')}"
            )

        if ctx and hasattr(ctx, "manifest") and ctx.manifest:
            ctx.manifest.set_sanitizer_output("lj_validation", report)

        if policy == "off":
            return []

        if policy == "error" and offenders:
            raise SanitizerError(
                f"LJ outlier policy=error triggered ({len(offenders)} offender(s)); "
                "inspect sanitizer lj_validation report in manifest."
            )

        for msg in warnings[:20]:
            print(f"  [WARN] {msg}")
        if len(warnings) > 20:
            print(
                f"  [WARN] ... {len(warnings) - 20} additional LJ outlier warning(s) omitted."
            )

        return warnings
    
    def _check_lj_outliers_robust(
        self, atomtype_data: List[Tuple[str, float, float, str]], source_name: str
    ) -> List[str]:
        """
        Check for LJ parameter outliers using robust statistics (MAD z-score on log-scale).
        
        For comb-rule=1 where we cannot assume sigma/epsilon semantics.
        Uses median and median absolute deviation (MAD) to detect outliers.
        
        Args:
            atomtype_data: List of (name, p1, p2, raw_line) tuples
            source_name: Source file name for warnings
            
        Returns:
            List of warning messages
        """
        warnings: List[str] = []
        
        # Filter to positive values only for log-scale analysis
        p1_positive = [(name, abs(p1)) for name, p1, _, _ in atomtype_data if p1 != 0]
        p2_positive = [(name, abs(p2)) for name, _, p2, _ in atomtype_data if p2 != 0]
        
        # Check for negative values (warn-only, not "impossible")
        for name, p1, p2, _ in atomtype_data:
            if p1 < 0 or p2 < 0:
                msg = (
                    f"Negative LJ params in atomtype '{name}': p1={p1:.4e}, p2={p2:.4e} "
                    f"at {source_name}; comb-rule=1 is not sigma/epsilon, please review"
                )
                print(f"  [WARN] {msg}")
                warnings.append(msg)
        
        # Robust outlier detection using MAD z-score on log scale
        z_threshold = 5.0  # Conservative threshold for outliers
        
        for label, values in [("p1", p1_positive), ("p2", p2_positive)]:
            if len(values) < 3:
                continue  # Not enough data for statistics
            
            # Compute log-scale values
            log_vals = [math.log10(v) for _, v in values]
            
            # Median
            sorted_logs = sorted(log_vals)
            n = len(sorted_logs)
            median = sorted_logs[n // 2] if n % 2 == 1 else (sorted_logs[n//2 - 1] + sorted_logs[n//2]) / 2
            
            # MAD (Median Absolute Deviation)
            abs_devs = sorted([abs(lv - median) for lv in log_vals])
            mad = abs_devs[n // 2] if n % 2 == 1 else (abs_devs[n//2 - 1] + abs_devs[n//2]) / 2
            
            # Scale factor for MAD -> sigma-equivalent (1.4826 for normal distribution)
            mad_scaled = mad * 1.4826 if mad > 0 else 1.0
            
            # Check for outliers
            for (name, val), log_val in zip(values, log_vals):
                if mad_scaled > 0:
                    z_score = abs(log_val - median) / mad_scaled
                    if z_score > z_threshold:
                        msg = (
                            f"Statistical outlier in atomtype '{name}': "
                            f"{label}={val:.4e} (z={z_score:.1f}) at {source_name}"
                        )
                        print(f"  [WARN] {msg}")
                        warnings.append(msg)
        
        return warnings

    def _sync_corrected_itp(
        self,
        correction_result,
        molecule_itp_paths: Dict[str, Path],
        sanitized_dir: Path,
        sanitized_current_dir: Path,
    ) -> None:
        """Ensure corrected ITP matches the include path used by system.top."""
        if not correction_result.corrected or not correction_result.patched_itp_path:
            return
        target = correction_result.target_molecule
        if not target and correction_result.patched_itp_path:
            name = correction_result.patched_itp_path.stem
            if name.endswith("_corrected"):
                target = name[: -len("_corrected")]
        target_casefold_map = {name.casefold(): name for name in molecule_itp_paths}
        resolved_target = target_casefold_map.get(target.casefold()) if target else None
        if not resolved_target:
            raise SanitizerError(
                "Charge correction produced a patched ITP but the target molecule "
                "cannot be resolved to an included ITP."
            )
        original_path = molecule_itp_paths[resolved_target]
        if correction_result.patched_itp_path.resolve() != original_path.resolve():
            self._atomic_replace_file(correction_result.patched_itp_path, original_path)
        versioned_path = sanitized_dir / original_path.name
        if versioned_path.resolve() != original_path.resolve():
            self._atomic_replace_file(original_path, versioned_path)
        # Remove stray corrected filename to avoid confusion
        if correction_result.patched_itp_path.exists() and correction_result.patched_itp_path.resolve() != original_path.resolve():
            correction_result.patched_itp_path.unlink()

    def _atomic_replace_file(self, src: Path, dest: Path) -> None:
        """Atomic replace of dest with src content."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        temp_path = dest.parent / f".{dest.name}.tmp"
        try:
            shutil.copy2(src, temp_path)
            temp_path.replace(dest)
        finally:
            if temp_path.exists():
                temp_path.unlink()

    def _atomic_write_text(self, path: Path, content: str, encoding: str = "utf-8") -> None:
        """
        Write text content atomically via temp file + replace.
        
        Prevents partial files after crashes.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = path.parent / f".{path.name}.tmp"
        try:
            temp_path.write_text(content, encoding=encoding)
            temp_path.replace(path)
        finally:
            if temp_path.exists():
                temp_path.unlink()

    def _build_include_search_paths(self, ctx: "PipelineContext", base_path: Path) -> List[Path]:
        """
        Build restricted search paths for #include resolution.
        
        Hardening v5 - Issue D: Priority-aware include resolution.
        Priority values are normalized with deterministic aliases:
        - sanitized_first (alias: first)
        - forcefield_first (alias: last)
        
        Order when priority == "forcefield_first" (legacy alias: "last"):
        1. Current file's parent directory
        2. system_gromacs_dir/top and /itp
        3. system_htpolynet_dir/itp
        4. Forcefield directory
        5. Sanitized output directory
        6. Configured GROMACS share/data dirs from context
        7. User-provided ctx.itp_include_dirs (LAST)

        Order when priority == "sanitized_first" (legacy alias: "first"):
        1. Current file's parent directory
        2. User-provided ctx.itp_include_dirs (FIRST - override position)
        3. system_gromacs_dir/top and /itp
        4. system_htpolynet_dir/itp
        5. Forcefield directory
        6. Sanitized output directory
        7. Configured GROMACS share/data dirs from context

        Relative user paths are resolved against ctx.project_root.
        
        Returns:
            List of existing directories to search (in priority order)
        """
        paths: List[Path] = [base_path.parent]
        project_root = getattr(ctx, "project_root", None)
        project_root_path = Path(project_root).resolve() if project_root else None

        def _resolve_path(raw_path: Path) -> Path:
            if raw_path.is_absolute() or project_root_path is None:
                return raw_path
            return project_root_path / raw_path

        # Collect user paths (resolve relative paths against project_root for reproducibility)
        user_paths: List[Path] = []
        if ctx.itp_include_dirs:
            for d in ctx.itp_include_dirs.split(","):
                d = d.strip()
                if not d:
                    continue
                resolved = _resolve_path(Path(d))
                if resolved.exists():
                    user_paths.append(resolved)

        # Configured GROMACS share/data dirs from context
        gmx_share_paths: List[Path] = []
        for attr in (
            "gromacs_share_dirs",
            "gromacs_share_dir",
            "gmx_share_dirs",
            "gmx_share_dir",
            "gromacs_datadir",
            "gmx_datadir",
            "gmxdata",
        ):
            value = getattr(ctx, attr, None)
            if not value:
                continue
            if isinstance(value, (list, tuple)):
                raw_items = [str(item) for item in value if str(item).strip()]
            else:
                raw_items = [item.strip() for item in str(value).split(",") if item.strip()]
            for item in raw_items:
                resolved = _resolve_path(Path(item))
                if resolved.exists():
                    gmx_share_paths.append(resolved)

        # Normalize priority setting and reject invalid values with clear guidance.
        raw_priority = getattr(ctx, "itp_include_dirs_priority", DEFAULT_INCLUDE_PRIORITY)
        priority = normalize_include_priority(raw_priority)

        # If priority == "sanitized_first", add user paths early (position 1, after base_path)
        if priority == "sanitized_first" and user_paths:
            paths.extend(user_paths)

        # System paths
        system_gromacs_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs")
        system_htpolynet_dir = ctx.get_input_path("systems", ctx.system_id, "htpolynet")

        paths.extend([
            system_gromacs_dir / "top",
            system_gromacs_dir / "itp",
            system_htpolynet_dir / "itp",
        ])

        # Forcefield directory (cached if available)
        if hasattr(self, "_forcefield_dir") and self._forcefield_dir:
            paths.append(self._forcefield_dir / "gromacs")

        # Sanitized output directory
        if hasattr(self, "_sanitized_current_dir") and self._sanitized_current_dir:
            paths.append(self._sanitized_current_dir)

        # Configured GROMACS share/data paths
        paths.extend(gmx_share_paths)

        # If priority == "forcefield_first" (default), add user paths at end
        if priority == "forcefield_first" and user_paths:
            paths.extend(user_paths)

        deduped: List[Path] = []
        seen: Set[Path] = set()
        for path in paths:
            if not path.exists():
                continue
            try:
                key = path.resolve()
            except OSError:
                key = path
            if key in seen:
                continue
            seen.add(key)
            deduped.append(path)
        return deduped




    def _warn_pairs_consistency(
        self, defaults, itp_paths: List[Path],
        source_defaults_map: Optional[Dict[Path, Any]] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> None:
        """
        Check gen-pairs/[pairs] mismatch and defaults consistency.
        
        Hardening v5 - Issue E: Also checks comb-rule, fudgeLJ, fudgeQQ
        consistency across merged sources.
        
        Args:
            defaults: DefaultsEntry with gen_pairs setting (primary)
            itp_paths: List of molecule ITP files to scan for [pairs]
            source_defaults_map: Dict of source_path -> DefaultsEntry for consistency check
            ctx: PipelineContext for strict mode flag
        """
        if defaults is None:
            return
        
        gen_pairs_val = getattr(defaults, "gen_pairs", "no").lower()
        gen_pairs_yes = gen_pairs_val == "yes"
        
        # Scan for [pairs] sections
        pairs_count = 0
        strict_reads = bool(ctx and ctx.strict_include_resolution)
        for itp in itp_paths:
            content = self.safe_read_text(itp, strict=strict_reads)
            if not content:
                continue
            for line in content.splitlines():
                section = parse_section_name(line)
                if section == "pairs":
                    pairs_count += 1
                    break
        
        if gen_pairs_yes and pairs_count > 0:
            print(
                f"  [WARN] gen-pairs=yes but {pairs_count} ITP(s) have explicit [pairs] sections "
                "(may duplicate 1-4 interactions; #ifdef may change actual behavior)"
            )
        elif not gen_pairs_yes and pairs_count == 0:
            print(
                f"  [WARN] gen-pairs=no but no [pairs] sections found "
                "(1-4 interactions may be missing; #ifdef may change actual behavior)"
            )
        
        # Mixed [defaults] consistency is enforced in ItpSanitizer.validate_defaults().
        # Keep this check informational only to avoid duplicate/contradictory failures.
        if source_defaults_map and len(source_defaults_map) > 1:
            signatures = {
                (
                    getattr(src_defaults, "nbfunc", "?"),
                    getattr(src_defaults, "comb_rule", "?"),
                    getattr(src_defaults, "gen_pairs", "?"),
                    getattr(src_defaults, "fudge_lj", "?"),
                    getattr(src_defaults, "fudge_qq", "?"),
                )
                for src_defaults in source_defaults_map.values()
            }
            if len(signatures) > 1:
                print(
                    "  [WARN] Mixed [defaults] signatures detected in source map "
                    "(already validated during sanitizer defaults pass)."
                )



    def _find_coords_gro(self, system_gromacs_dir: Path) -> Optional[Path]:
        """Locate a coordinate file for dipole estimation."""
        coords_path = system_gromacs_dir / "coords" / "system.gro"
        if coords_path.exists():
            return coords_path
        gro_path = system_gromacs_dir / "gro" / "system.gro"
        if gro_path.exists():
            return gro_path
        return None

    def _read_gro_coords(self, gro_path: Path) -> Optional[List[Tuple[float, float, float]]]:
        """Read coordinates from a GRO file (nm). Legacy wrapper."""
        result = self._parse_gro_file(gro_path)
        if result is None:
            return None
        return [(a.x, a.y, a.z) for a in result.atoms]

    def _parse_gro_file(
        self,
        gro_path: Path,
        strict: bool = False,
        max_atoms: Optional[int] = None,
    ) -> Optional[GROParseResult]:
        """
        Parse a GRO file returning atoms with residue info and box vectors.
        
        Hardening v6: Robust parsing with truncation/corruption detection.
        Hardening v7: Streamed parsing (avoid loading huge files into memory).
        
        GRO format (fixed columns):
        - Line 1: title
        - Line 2: natoms (integer count)
        - Lines 3 to 2+natoms: atom lines
        - Last line: box vectors (3 or 9 floats)
        
        Atom line format (fixed columns):
        - Columns 1-5: residue number
        - Columns 6-10: residue name  
        - Columns 11-15: atom name
        - Columns 16-20: atom number
        - Columns 21-28: x coordinate (nm)
        - Columns 29-36: y coordinate (nm)
        - Columns 37-44: z coordinate (nm)
        
        Args:
            gro_path: Path to GRO file
            strict: If True, fail on any parsing issue; else return None
            max_atoms: Optional safety cap for declared atoms (non-strict early exit)
            
        Returns:
            GROParseResult with atoms, box dimensions, and triclinic flag
            
        Raises:
            SanitizerError: In strict mode on truncation or format errors
        """
        try:
            with open(gro_path, "r") as handle:
                title_line = handle.readline()
                count_line = handle.readline()
                if not title_line or not count_line:
                    if strict:
                        raise SanitizerError(
                            f"GRO file {gro_path} is too short (missing title/count line)"
                        )
                    return None
                title = title_line.strip()
                try:
                    natoms_declared = int(count_line.strip())
                except ValueError:
                    if strict:
                        raise SanitizerError(
                            f"GRO file {gro_path}: line 2 is not a valid atom count: {count_line!r}"
                        )
                    return None

                if natoms_declared < 0:
                    if strict:
                        raise SanitizerError(
                            f"GRO file {gro_path}: negative atom count declared: {natoms_declared}"
                        )
                    return None

                # Safety cap for absurdly large files (non-strict only)
                cap = max_atoms if max_atoms is not None else GRO_ABSURD_ATOMS
                if cap is not None and natoms_declared > cap:
                    msg = (
                        f"GRO file {gro_path}: declared atom count {natoms_declared} exceeds safety cap {cap}"
                    )
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                    return None

                # Parse atom lines (line numbers start at 3 for first atom)
                atoms: List[GROAtom] = []
                parse_errors: List[str] = []
                for i in range(natoms_declared):
                    line = handle.readline()
                    if not line:
                        parse_errors.append(f"Line {i + 3}: EOF before atom data")
                        break
                    if len(line) < 20:
                        parse_errors.append(
                            f"Line {i + 3}: too short for atom ({len(line)} chars)"
                        )
                        continue
                    try:
                        resnr = int(line[0:5].strip())
                        resname = line[5:10].strip()
                        atomname = line[10:15].strip()
                        x = float(line[20:28])
                        y = float(line[28:36])
                        z = float(line[36:44])
                        atoms.append(GROAtom(resnr, resname, atomname, x, y, z))
                    except (ValueError, IndexError) as e:
                        parse_errors.append(f"Line {i + 3}: {e}")
                        continue

                # Hardening v6: Verify natoms matches parsed atoms
                atoms_parsed = len(atoms)
                if atoms_parsed != natoms_declared:
                    msg = (
                        f"GRO file {gro_path}: atom count mismatch!\n"
                        f"  Declared: {natoms_declared}, Parsed: {atoms_parsed}\n"
                        f"  File may be truncated or corrupted."
                    )
                    if parse_errors:
                        msg += f"\n  Parse errors: {parse_errors[:5]}"
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                    if atoms_parsed == 0:
                        return None

                # Find the box line robustly by scanning remaining lines
                box_floats: List[float] = []
                last_nonempty: Optional[str] = None
                for line in handle:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    last_nonempty = stripped
                    parts = stripped.split()
                    try:
                        floats = [float(x) for x in parts]
                        if len(floats) == 3 or len(floats) == 9:
                            box_floats = floats
                            break
                    except ValueError:
                        continue

                if not box_floats and last_nonempty:
                    parts = last_nonempty.split()
                    try:
                        floats = [float(x) for x in parts]
                        if len(floats) == 3 or len(floats) == 9:
                            box_floats = floats
                    except ValueError:
                        pass

                if not box_floats:
                    msg = f"GRO file {gro_path}: no valid box vector line found"
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                    box_floats = [0.0, 0.0, 0.0]
        except Exception as e:
            if strict:
                raise SanitizerError(f"Cannot read GRO file {gro_path}: {e}")
            return None
        
        is_triclinic = len(box_floats) >= 9
        box_x = box_floats[0] if len(box_floats) >= 1 else 0.0
        box_y = box_floats[1] if len(box_floats) >= 2 else 0.0
        box_z = box_floats[2] if len(box_floats) >= 3 else 0.0
        
        return GROParseResult(
            atoms=atoms,
            box_x=box_x,
            box_y=box_y,
            box_z=box_z,
            is_triclinic=is_triclinic,
            title=title,
        )

    def _get_molecule_charges(self, itp_path: Path, mol_name: str) -> Optional[List[float]]:
        """Extract charges list for a specific moleculetype."""
        from ..charge_neutrality import ChargeParser
        per_mol = ChargeParser.parse_atoms_charges_per_molecule(itp_path)
        if mol_name in per_mol:
            return per_mol[mol_name][0]
        key = mol_name.casefold()
        matches = [name for name in per_mol if name.casefold() == key]
        if len(matches) == 1:
            return per_mol[matches[0]][0]
        return None

    def _compute_dipole_shift(
        self,
        system_gromacs_dir: Path,
        ordered_molecules: List[Tuple[str, int]],
        molecule_itp_paths: Dict[str, Path],
        correction_result,
        ctx: Optional["PipelineContext"] = None,
    ):
        """
        Compute approximate dipole shift for corrected molecule(s).
        
        PBC-aware implementation using residue grouping and topology-aware
        bond-graph unwrapping when bonds are available.
        
        v4: Adds gating logic for heuristic mode:
        - Giant molecule detection (>dipole_unwrap_max_atoms)
        - Net charge guard (skip enforcement if molecule is charged)
        - Percolation check via bond-loop image inconsistency (preferred)
          with span-based fallback only when bond graph is unavailable
        - Triclinic fallback to diagonal approximation
        """
        from ..charge_neutrality import ChargeParser, DipoleResult
        
        skip_reasons: List[str] = []
        
        if not correction_result.patched_itp_path:
            return None
        
        target = correction_result.target_molecule
        if not target and correction_result.patched_itp_path:
            name = correction_result.patched_itp_path.stem
            if name.endswith("_corrected"):
                target = name[: -len("_corrected")]
        if not target:
            return DipoleResult(computed=False)
        target_casefold_map = {name.casefold(): name for name in molecule_itp_paths}
        target_resolved = target_casefold_map.get(target.casefold())
        if not target_resolved:
            return DipoleResult(computed=False)
        target = target_resolved
        target_itp_name = (
            getattr(correction_result, "target_itp_moleculetype_name", None) or target
        )

        charges_before = self._get_molecule_charges(molecule_itp_paths[target], target_itp_name)
        charges_after = self._get_molecule_charges(correction_result.patched_itp_path, target_itp_name)
        if not charges_before or not charges_after or len(charges_before) != len(charges_after):
            return DipoleResult(computed=False)
        bond_source = correction_result.patched_itp_path or molecule_itp_paths[target]
        bond_graph = ChargeParser.parse_molecule_bond_graph(bond_source, target_itp_name)
        has_bond_graph = bool(bond_graph)
        
        # v4: Giant molecule gating (skip enforcement for large molecules)
        expected_natoms = len(charges_before)
        # Hardening v5: Use explicit ctx parameter (no more getattr(self, '_ctx'))
        dipole_max_atoms = ctx.dipole_unwrap_max_atoms if ctx else 512
        if expected_natoms > dipole_max_atoms:
            reason = f"giant_molecule ({expected_natoms} > {dipole_max_atoms} atoms)"
            skip_reasons.append(reason)
            print(f"  [INFO] Dipole enforcement gated: {reason}")
            gro_path = self._find_coords_gro(system_gromacs_dir)
            return DipoleResult(
                computed=False,
                skip_reasons=skip_reasons,
                coordinates_file=str(gro_path) if gro_path else None,
            )

        gro_path = self._find_coords_gro(system_gromacs_dir)
        if not gro_path:
            return DipoleResult(computed=False)

        # Use richer parser
        gro_data = self._parse_gro_file(gro_path)
        if gro_data is None or not gro_data.atoms:
            return DipoleResult(computed=False)

        # Skip if box dimensions are zero/invalid
        if gro_data.box_x <= 0 or gro_data.box_y <= 0 or gro_data.box_z <= 0:
            print("  [INFO] Dipole shift check skipped (reason: invalid box dimensions)")
            return DipoleResult(computed=False, coordinates_file=str(gro_path))
        
        # v4: Net charge guard (dipole is origin-dependent for charged molecules)
        charge_tol = ctx.charge_neutrality_tol if ctx else 1e-6
        net_before = sum(charges_before)
        net_after = sum(charges_after)
        if abs(net_before) > charge_tol or abs(net_after) > charge_tol:
            reason = f"charged_molecule (net_q={net_before:.4e}/{net_after:.4e})"
            skip_reasons.append(reason)
            print(f"  [INFO] Dipole enforcement gated: {reason}")
        
        # Hardening v6 - Issue C: Triclinic box policy
        # CRITICAL: Never apply orthorhombic MIC to triclinic boxes - it gives wrong results
        # Proper triclinic MIC requires box matrix inversion which is NOT implemented
        triclinic_policy_used: Optional[str] = None
        if gro_data.is_triclinic:
            # Determine effective policy
            policy = getattr(ctx, "dipole_triclinic_policy", "skip") if ctx else "skip"
            triclinic_policy_used = policy
            
            # If strict mode and no explicit policy, treat as "error"
            if ctx and ctx.strict_dipole_check and policy == "skip":
                # Check if user explicitly set the policy (vs default)
                # Default is "skip", so if strict we escalate to error unless user overrode
                if not hasattr(ctx, "_dipole_triclinic_policy_explicit"):
                    policy = "error"
                    triclinic_policy_used = "error (escalated from skip in strict mode)"
            
            if policy == "error":
                raise SanitizerError(
                    "Cannot compute accurate dipole shift for triclinic box.\n"
                    "Orthorhombic minimum-image convention (MIC) is mathematically invalid "
                    "for tilted boxes - it will give incorrect results.\n\n"
                    "Options:\n"
                    "  1. Re-box to orthorhombic: gmx editconf -bt cubic\n"
                    "  2. Skip check: --dipole-triclinic-policy skip\n"
                    "  3. Disable strict mode: --no-strict-dipole-check"
                )
            elif policy == "mic":
                # Hardening v6: REMOVED misleading "correct MIC" claim
                # Proper triclinic MIC requires matrix inversion - NOT IMPLEMENTED
                # Fail clearly rather than silently giving wrong results
                raise SanitizerError(
                    "Triclinic MIC policy='mic' is not supported.\n"
                    "Proper triclinic minimum-image convention requires box matrix inversion "
                    "which is not implemented in this version.\n\n"
                    "Options:\n"
                    "  1. Re-box to orthorhombic: gmx editconf -bt cubic\n"
                    "  2. Skip check: --dipole-triclinic-policy skip\n"
                    "  3. Use GROMACS 'gmx dipoles' for accurate triclinic dipole analysis"
                )
            else:  # policy == "skip"
                skip_reasons.append("triclinic_box (orthorhombic MIC invalid for tilted boxes)")
                print("  [INFO] Triclinic box detected; skipping dipole check (policy=skip)")
                return DipoleResult(
                    computed=False,
                    skip_reasons=skip_reasons,
                    coordinates_file=str(gro_path),
                )
        
        # Task A: Group atoms by (resnr, resname), match target by resname
        # Handle GRO 5-char truncation: try target and target[:5]
        target_resnames = {target.upper(), target[:5].upper()}
        
        # Build residue groups: (resnr, resname) -> list of atoms
        residue_groups: Dict[Tuple[int, str], List[GROAtom]] = {}
        for atom in gro_data.atoms:
            key = (atom.resnr, atom.resname.upper())
            if key not in residue_groups:
                residue_groups[key] = []
            residue_groups[key].append(atom)
        
        # Find groups matching target
        target_groups: List[List[GROAtom]] = []
        for (resnr, resname), atoms in residue_groups.items():
            if resname in target_resnames:
                target_groups.append(atoms)
        
        if not target_groups:
            # No matching residues found
            return DipoleResult(computed=False)
        
        group_max_atoms = (
            getattr(ctx, "dipole_group_max_atoms", dipole_max_atoms) if ctx else dipole_max_atoms
        )
        valid_groups: List[List[GROAtom]] = []
        for group in target_groups:
            group_size = len(group)
            if group_size > group_max_atoms:
                reason = f"group_too_large ({group_size} > {group_max_atoms} atoms)"
                if reason not in skip_reasons:
                    skip_reasons.append(reason)
                    print(f"  [INFO] Dipole enforcement gated: {reason}")
                continue
            if group_size == expected_natoms:
                valid_groups.append(group)
        
        if not valid_groups:
            # Atom count mismatch in all groups or all groups gated
            return DipoleResult(computed=False, skip_reasons=skip_reasons, coordinates_file=str(gro_path))
        
        # Compute dipole for each valid group with PBC unwrapping.
        # Prefer topology-aware bond graph unwrap; fallback to span-gated MIC unwrap.
        # v4: Issue 7 - Per-group vector delta: ||mu_after_vec - mu_before_vec||
        box = (gro_data.box_x, gro_data.box_y, gro_data.box_z)
        delta_magnitudes: List[float] = []  # Per-group vector delta norms
        mu_before_mags: List[float] = []  # For reporting
        mu_after_mags: List[float] = []
        percolation_threshold = (
            ctx.dipole_percolation_threshold if ctx else 0.45
        )
        
        for group in valid_groups:
            # Fast wrapped-span pre-gate: conservative early skip for obvious spanning groups.
            wrapped_span = self._wrapped_span(group)
            wrapped_span_ratio = max(
                wrapped_span[0] / box[0] if box[0] > 0 else 0,
                wrapped_span[1] / box[1] if box[1] > 0 else 0,
                wrapped_span[2] / box[2] if box[2] > 0 else 0,
            )
            if wrapped_span_ratio > percolation_threshold:
                print(
                    "  [INFO] Dipole wrapped-span suspicion: "
                    f"ratio={wrapped_span_ratio:.2f} > {percolation_threshold:.2f}. "
                    "Continuing with unwrap/percolation checks."
                )

            if has_bond_graph:
                coords_list, topo_percolated, topo_msg = self._unwrap_molecule_bond_graph(
                    group,
                    box,
                    bond_graph,
                )
                if topo_percolated:
                    reason = "percolated_topology"
                    if reason not in skip_reasons:
                        skip_reasons.append(reason)
                        detail = f" ({topo_msg})" if topo_msg else ""
                        print(f"  [INFO] Dipole enforcement gated: {reason}{detail}")
                    continue
            else:
                coords_list, span, percolated = self._unwrap_molecule_pbc_with_span(
                    group,
                    box,
                    early_exit_ratio=percolation_threshold,
                )
                if percolated:
                    max_span_ratio = max(
                        span[0] / box[0] if box[0] > 0 else 0,
                        span[1] / box[1] if box[1] > 0 else 0,
                        span[2] / box[2] if box[2] > 0 else 0,
                    )
                    reason = (
                        f"percolated_unwrapped_span (ratio={max_span_ratio:.2f} > "
                        f"{percolation_threshold})"
                    )
                    if reason not in skip_reasons:
                        skip_reasons.append(reason)
                        print(f"  [INFO] Dipole enforcement gated: {reason}")
                    continue
            if not coords_list:
                continue
            
            # Compute vector dipoles (returns (mx, my, mz) in Debye)
            mu_before_vec = self._dipole_vector_unwrapped(coords_list, charges_before)
            mu_after_vec = self._dipole_vector_unwrapped(coords_list, charges_after)
            
            # v4 Issue 7: Per-group vector delta norm
            delta_vec = (
                mu_after_vec[0] - mu_before_vec[0],
                mu_after_vec[1] - mu_before_vec[1],
                mu_after_vec[2] - mu_before_vec[2],
            )
            delta_mag = math.sqrt(delta_vec[0]**2 + delta_vec[1]**2 + delta_vec[2]**2)
            delta_magnitudes.append(delta_mag)
            
            # Also store magnitudes for reporting
            mu_before_mags.append(math.sqrt(sum(x**2 for x in mu_before_vec)))
            mu_after_mags.append(math.sqrt(sum(x**2 for x in mu_after_vec)))
        
        if not delta_magnitudes:
            return DipoleResult(
                computed=False,
                skip_reasons=skip_reasons,
                coordinates_file=str(gro_path),
            )
        
        # Report average and max delta
        avg_before = sum(mu_before_mags) / len(mu_before_mags)
        avg_after = sum(mu_after_mags) / len(mu_after_mags)
        avg_delta = sum(delta_magnitudes) / len(delta_magnitudes)
        max_delta = max(delta_magnitudes)
        
        # Use average delta as the primary metric
        return DipoleResult(
            dipole_before_debye=avg_before,
            dipole_after_debye=avg_after,
            delta_dipole_debye=avg_delta,
            max_delta_dipole_debye=max_delta,
            coordinates_file=str(gro_path),
            computed=True,
            skip_reasons=skip_reasons,
        )

    def _unwrap_molecule_pbc(
        self, atoms: List[GROAtom], box: Tuple[float, float, float]
    ) -> List[GROAtom]:
        """
        Apply minimum-image unwrapping to a molecule's atoms (orthorhombic).
        
        Task A: Uses first atom as reference, unwraps others relative to it.
        This ensures a molecule straddling a PBC boundary is contiguous.
        """
        if not atoms:
            return atoms
        
        Lx, Ly, Lz = box
        ref = atoms[0]
        unwrapped = [ref]
        
        for atom in atoms[1:]:
            dx = atom.x - ref.x
            dy = atom.y - ref.y
            dz = atom.z - ref.z
            
            # Minimum image convention
            dx -= round(dx / Lx) * Lx
            dy -= round(dy / Ly) * Ly
            dz -= round(dz / Lz) * Lz
            
            unwrapped.append(GROAtom(
                resnr=atom.resnr,
                resname=atom.resname,
                atomname=atom.atomname,
                x=ref.x + dx,
                y=ref.y + dy,
                z=ref.z + dz,
            ))
        
        return unwrapped
    
    def _unwrap_molecule_pbc_with_span(
        self,
        atoms: List[GROAtom],
        box: Tuple[float, float, float],
        early_exit_ratio: Optional[float] = None,
    ) -> Tuple[Optional[List[Tuple[float, float, float]]], Tuple[float, float, float], bool]:
        """
        Apply minimum-image unwrapping with span tracking.

        Returns compact coordinates (x, y, z) tuples instead of GROAtom objects.
        Supports optional early exit when span ratio exceeds a threshold.
        """
        if not atoms:
            return [], (0.0, 0.0, 0.0), False
        
        Lx, Ly, Lz = box
        ref = atoms[0]
        
        # Initialize min/max tracking with reference atom
        min_x = max_x = ref.x
        min_y = max_y = ref.y
        min_z = max_z = ref.z
        
        coords: List[Tuple[float, float, float]] = [(ref.x, ref.y, ref.z)]
        
        def _span_ratio() -> float:
            span_x = max_x - min_x
            span_y = max_y - min_y
            span_z = max_z - min_z
            return max(
                span_x / Lx if Lx > 0 else 0.0,
                span_y / Ly if Ly > 0 else 0.0,
                span_z / Lz if Lz > 0 else 0.0,
            )
        
        for atom in atoms[1:]:
            dx = atom.x - ref.x
            dy = atom.y - ref.y
            dz = atom.z - ref.z
            
            # Minimum image convention
            dx -= round(dx / Lx) * Lx
            dy -= round(dy / Ly) * Ly
            dz -= round(dz / Lz) * Lz
            
            new_x = ref.x + dx
            new_y = ref.y + dy
            new_z = ref.z + dz
            
            # Stream min/max update
            min_x = min(min_x, new_x)
            max_x = max(max_x, new_x)
            min_y = min(min_y, new_y)
            max_y = max(max_y, new_y)
            min_z = min(min_z, new_z)
            max_z = max(max_z, new_z)
            
            coords.append((new_x, new_y, new_z))

            if early_exit_ratio is not None and _span_ratio() > early_exit_ratio:
                span = (max_x - min_x, max_y - min_y, max_z - min_z)
                return None, span, True
        
        span = (max_x - min_x, max_y - min_y, max_z - min_z)
        return coords, span, False

    def _unwrap_molecule_bond_graph(
        self,
        atoms: List[GROAtom],
        box: Tuple[float, float, float],
        bond_graph: Dict[int, Set[int]],
    ) -> Tuple[Optional[List[Tuple[float, float, float]]], bool, Optional[str]]:
        """
        Unwrap coordinates via bond-graph traversal and detect topology percolation.

        For revisits through alternate graph paths, inconsistent implied image vectors
        indicate percolation/topological wrapping and the dipole check is skipped.
        """
        if not atoms:
            return [], False, None

        natoms = len(atoms)
        Lx, Ly, Lz = box
        adjacency: Dict[int, Set[int]] = {i: set() for i in range(1, natoms + 1)}
        for i, neighbors in bond_graph.items():
            if i < 1 or i > natoms:
                continue
            for j in neighbors:
                if j < 1 or j > natoms or i == j:
                    continue
                adjacency[i].add(j)
                adjacency[j].add(i)

        if not any(adjacency.values()):
            coords, _, _ = self._unwrap_molecule_pbc_with_span(atoms, box, early_exit_ratio=None)
            return coords, False, None

        unwrapped: Dict[int, Tuple[float, float, float]] = {}
        images: Dict[int, Tuple[int, int, int]] = {}
        tol = 1e-6

        for root in range(1, natoms + 1):
            if root in unwrapped:
                continue
            root_atom = atoms[root - 1]
            unwrapped[root] = (root_atom.x, root_atom.y, root_atom.z)
            images[root] = (0, 0, 0)
            queue = deque([root])

            while queue:
                i = queue.popleft()
                xi, yi, zi = unwrapped[i]
                img_i = images[i]
                atom_i = atoms[i - 1]

                for j in sorted(adjacency.get(i, set())):
                    atom_j = atoms[j - 1]
                    dx_w = atom_j.x - atom_i.x
                    dy_w = atom_j.y - atom_i.y
                    dz_w = atom_j.z - atom_i.z

                    sx = int(round(dx_w / Lx)) if Lx > 0 else 0
                    sy = int(round(dy_w / Ly)) if Ly > 0 else 0
                    sz = int(round(dz_w / Lz)) if Lz > 0 else 0

                    dx = dx_w - sx * Lx
                    dy = dy_w - sy * Ly
                    dz = dz_w - sz * Lz

                    implied_pos = (xi + dx, yi + dy, zi + dz)
                    implied_img = (img_i[0] - sx, img_i[1] - sy, img_i[2] - sz)

                    if j not in unwrapped:
                        unwrapped[j] = implied_pos
                        images[j] = implied_img
                        queue.append(j)
                        continue

                    existing_pos = unwrapped[j]
                    existing_img = images[j]
                    if implied_img != existing_img:
                        diff_img = (
                            implied_img[0] - existing_img[0],
                            implied_img[1] - existing_img[1],
                            implied_img[2] - existing_img[2],
                        )
                        if diff_img != (0, 0, 0):
                            return None, True, f"image_loop_shift={diff_img}"

                    dx_pos = implied_pos[0] - existing_pos[0]
                    dy_pos = implied_pos[1] - existing_pos[1]
                    dz_pos = implied_pos[2] - existing_pos[2]
                    if abs(dx_pos) <= tol and abs(dy_pos) <= tol and abs(dz_pos) <= tol:
                        continue

                    kx = int(round(dx_pos / Lx)) if Lx > 0 else 0
                    ky = int(round(dy_pos / Ly)) if Ly > 0 else 0
                    kz = int(round(dz_pos / Lz)) if Lz > 0 else 0
                    residual_ok = (
                        abs(dx_pos - kx * Lx) <= tol
                        and abs(dy_pos - ky * Ly) <= tol
                        and abs(dz_pos - kz * Lz) <= tol
                    )
                    if residual_ok and (kx, ky, kz) != (0, 0, 0):
                        return None, True, f"loop_lattice_shift={(kx, ky, kz)}"
                    return None, True, "inconsistent_bond_loop"

        coords: List[Tuple[float, float, float]] = []
        for idx in range(1, natoms + 1):
            coords.append(unwrapped.get(idx, (atoms[idx - 1].x, atoms[idx - 1].y, atoms[idx - 1].z)))
        return coords, False, None

    def _wrapped_span(self, atoms: List[GROAtom]) -> Tuple[float, float, float]:
        """Compute axis-aligned span from wrapped coordinates without allocations."""
        if not atoms:
            return (0.0, 0.0, 0.0)
        min_x = max_x = atoms[0].x
        min_y = max_y = atoms[0].y
        min_z = max_z = atoms[0].z
        for atom in atoms[1:]:
            x = atom.x
            y = atom.y
            z = atom.z
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)
            min_z = min(min_z, z)
            max_z = max(max_z, z)
        return (max_x - min_x, max_y - min_y, max_z - min_z)
    
    def _dipole_vector_unwrapped(
        self, coords: List[Tuple[float, float, float]], charges: List[float]
    ) -> Tuple[float, float, float]:
        """
        Compute dipole vector (Debye) using geometric center.
        
        v4 Issue 7: Returns vector components for proper delta calculation.
        
        Args:
            coords: List of (x, y, z) tuples in nm (should be PBC-unwrapped)
            charges: List of partial charges (e)
            
        Returns:
            Dipole vector (mx, my, mz) in Debye
        """
        if not coords or len(coords) != len(charges):
            return (0.0, 0.0, 0.0)
        
        # Geometric center for origin
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        
        # Compute dipole components relative to geometric center
        mux = muy = muz = 0.0
        for (x, y, z), q in zip(coords, charges):
            mux += q * (x - cx)
            muy += q * (y - cy)
            muz += q * (z - cz)
        
        # Convert from e*nm to Debye: 1 e*nm = 48.0321 Debye
        return (mux * 48.0321, muy * 48.0321, muz * 48.0321)

    def _dipole_magnitude(self, coords: List[Tuple[float, float, float]], charges: List[float]) -> float:
        """Compute dipole magnitude (Debye) using molecule-centered coordinates. Legacy."""
        return self._dipole_magnitude_unwrapped(coords, charges)

    def _dipole_magnitude_unwrapped(
        self, coords: List[Tuple[float, float, float]], charges: List[float]
    ) -> float:
        """
        Compute dipole magnitude (Debye) using geometric center.
        
        Physics Note (Issue E hardening):
        ---------------------------------
        The dipole moment μ = Σ q_i * r_i is ORIGIN-DEPENDENT for charged molecules.
        For a molecule with net charge Q ≠ 0, shifting the origin by vector R changes
        the dipole by ΔQ * R.
        
        For NEUTRAL molecules (which we're correcting toward), the dipole is
        origin-independent, so any choice of center gives the same result.
        
        We use GEOMETRIC CENTER because:
        1. It minimizes numerical error from large |r| values
        2. It provides consistent comparison before/after charge correction
        3. For the PURPOSE of checking dipole SHIFT (Δμ = μ_after - μ_before),
           the origin choice CANCELS OUT as long as it's consistent
        
        Alternative: Charge-weighted center could be used if comparing absolute
        dipole magnitudes to experimental values, but for SHIFT calculations
        (our use case), geometric center is equivalent and simpler.
        
        Units:
        - Input coords: nm (GROMACS native)
        - Input charges: elementary charges (e)
        - Output: Debye (1 e·nm = 48.0321 Debye)
        
        Args:
            coords: List of (x, y, z) tuples in nm (should be PBC-unwrapped)
            charges: List of partial charges (e)
            
        Returns:
            Dipole magnitude in Debye
        """
        if not coords or len(coords) != len(charges):
            return 0.0
        
        # Geometric center for origin
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        
        # Compute dipole components relative to geometric center
        mux = muy = muz = 0.0
        for (x, y, z), q in zip(coords, charges):
            mux += q * (x - cx)
            muy += q * (y - cy)
            muz += q * (z - cz)
        
        mu_nm = math.sqrt(mux * mux + muy * muy + muz * muz)
        # Convert from e*nm to Debye: 1 e*nm = 48.0321 Debye
        return mu_nm * 48.0321

    def _expected_atoms_from_grompp(
        self,
        gmx_exe: str,
        system_top_path: Path,
        gro_path: Path,
        output_dir: Path,
        ctx: "PipelineContext",
    ) -> Optional[int]:
        """
        Run grompp -pp to preprocess topology and parse expected atom count.
        
        Hardening v6: Complexity-based timeout estimation with retry logic.
        """
        workdir = output_dir / "grompp_preprocess"
        workdir.mkdir(parents=True, exist_ok=True)
        mdp_path = self._select_mdp_for_preprocess(ctx, workdir)
        pp_path = workdir / "system.preprocessed.top"
        tpr_path = workdir / "preprocess.tpr"
        cmd = [
            gmx_exe,
            "grompp",
            "-f",
            str(mdp_path),
            "-c",
            str(gro_path),
            "-p",
            str(system_top_path),
            "-pp",
            str(pp_path),
            "-o",
            str(tpr_path),
            "-maxwarn",
            str(ctx.grompp_maxwarn),
        ]
        
        # Hardening v6: Complexity-based timeout estimation
        complexity_metrics = self._estimate_topology_complexity(system_top_path, ctx)
        
        # Base timeout from ctx > env var > default
        base_timeout_s = getattr(ctx, "grompp_preprocess_timeout_s", None)
        if base_timeout_s is None:
            env_timeout = os.environ.get("GMX_GROMPP_PP_TIMEOUT_S", "")
            if env_timeout:
                try:
                    base_timeout_s = int(env_timeout)
                except ValueError:
                    base_timeout_s = DEFAULT_GROMPP_PP_TIMEOUT_S
            else:
                base_timeout_s = DEFAULT_GROMPP_PP_TIMEOUT_S
        
        # Compute complexity-adjusted timeout
        timeout_s = GROMPP_TIMEOUT_FLOOR_S
        timeout_s += complexity_metrics.get("lines", 0) // 1000 * GROMPP_TIMEOUT_PER_1K_LINES
        timeout_s += complexity_metrics.get("includes", 0) // 10 * GROMPP_TIMEOUT_PER_10_INCLUDES
        
        # Add GRO atom-based scaling
        gro_atoms = self._count_atoms_in_gro(gro_path)
        if gro_atoms is not None and gro_atoms > 10000:
            extra_time = ((gro_atoms - 10000) // 10000) * 60
            if gro_atoms > 100000:
                extra_time = int(extra_time * 1.5) + 600
            timeout_s += extra_time
        
        # Never reduce below user/env-set timeout or floor
        timeout_s = max(timeout_s, base_timeout_s, GROMPP_TIMEOUT_FLOOR_S)
        
        # Hardening v6: Run with retry logic
        retry_count = 0
        max_retries = 1  # Try once, then retry with longer timeout
        current_timeout = timeout_s
        
        while retry_count <= max_retries:
            try:
                result = subprocess.run(
                    cmd,
                    cwd=str(workdir),
                    capture_output=True,
                    text=True,
                    timeout=current_timeout,
                )
                break  # Success
            except subprocess.TimeoutExpired:
                retry_count += 1
                if retry_count <= max_retries:
                    # Retry with longer timeout
                    current_timeout = int(current_timeout * GROMPP_TIMEOUT_RETRY_MULTIPLIER)
                    print(
                        f"  [WARN] grompp preprocessing timed out after {timeout_s}s, "
                        f"retrying with {current_timeout}s timeout..."
                    )
                    continue
                else:
                    # Final failure
                    msg = (
                        f"grompp preprocessing timed out after {current_timeout}s (including retry).\n"
                        f"  Command: {' '.join(cmd)}\n"
                        f"  Complexity metrics: {complexity_metrics}\n"
                        f"  GRO atoms: {gro_atoms}\n"
                        f"  To increase timeout:\n"
                        f"    - Set environment variable: GMX_GROMPP_PP_TIMEOUT_S={int(current_timeout * 2)}\n"
                        f"    - Or add grompp_preprocess_timeout_s to context\n"
                        f"  Possible causes:\n"
                        f"    - Very large system (>100k atoms)\n"
                        f"    - Complex topology with many includes\n"
                        f"    - Slow filesystem or network-mounted storage"
                    )
                    print(f"  [ERROR] {msg}")
                    if ctx.strict_gro_top_check:
                        raise SanitizerError(
                            f"grompp preprocessing timed out (strict mode). {msg}"
                        )
                    return None
            except Exception as e:
                print(f"  [WARN] grompp preprocess failed: {e}")
                if ctx.strict_gro_top_check:
                    raise SanitizerError(f"grompp preprocess failed for GRO/TOP validation: {e}")
                return None
        
        if result.returncode != 0 or not pp_path.exists():
            print(f"  [WARN] grompp preprocess returned code {result.returncode}")
            if ctx.strict_gro_top_check:
                raise SanitizerError("grompp preprocess failed for GRO/TOP validation.")
            return None
        
        expected_atoms = self._parse_preprocessed_top(pp_path)
        if expected_atoms is None:
            print("  [WARN] Could not parse expected atoms from preprocessed topology.")
        
        # Hardening v6: Record complexity metrics in manifest
        if ctx.manifest:
            ctx.manifest.set_sanitizer_output("grompp_preprocess", {
                "timeout_used_s": current_timeout,
                "retry_count": retry_count,
                "complexity_metrics": complexity_metrics,
                "gro_atoms": gro_atoms,
                "expected_atoms": expected_atoms,
            })
        
        return expected_atoms
    
    def _estimate_topology_complexity(
        self, top_path: Path, ctx: "PipelineContext"
    ) -> Dict[str, int]:
        """
        Estimate topology complexity for timeout calculation.
        
        Hardening v6: Counts lines, includes, and other complexity indicators.
        
        Returns:
            Dict with 'lines', 'includes', 'moleculetypes' counts
        """
        metrics = {"lines": 0, "includes": 0, "moleculetypes": 0}
        
        if not top_path.exists():
            return metrics

        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return metrics
        lines = content.splitlines()
        metrics["lines"] = len(lines)

        for line in lines:
            if line.strip().startswith("#include"):
                metrics["includes"] += 1
            section = parse_section_name(line)
            if section == "moleculetype":
                metrics["moleculetypes"] += 1
        
        return metrics

    def _select_mdp_for_preprocess(self, ctx: "PipelineContext", workdir: Path) -> Path:
        """Select an existing MDP for preprocessing or create a minimal one."""
        mdp_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "mdp")
        if mdp_dir.exists():
            mdp_files = self._sorted_paths_casefold(mdp_dir.glob("*.mdp"))
            for preferred in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]:
                for mdp in mdp_files:
                    if mdp.name.lower() == preferred:
                        return mdp
            if mdp_files:
                return mdp_files[0]
        minimal_mdp = workdir / "minimal.mdp"
        self._atomic_write_text(
            minimal_mdp,
            "integrator = md\n"
            "nsteps = 0\n"
            "cutoff-scheme = Verlet\n"
            "coulombtype = PME\n"
            "rcoulomb = 1.0\n"
            "vdwtype = Cut-off\n"
            "rvdw = 1.0\n"
            "pbc = xyz\n"
        )
        return minimal_mdp

    def _parse_preprocessed_top(self, top_path: Path) -> Optional[int]:
        """Parse a preprocessed topology to compute expected atom count."""
        if not top_path.exists():
            return None
        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return None
        current_section = ""
        current_mol = None
        mol_atom_counts: Dict[str, int] = {}
        molecules_counts: Dict[str, int] = {}
        for line in content.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if current_section == "moleculetype":
                    current_mol = None
                continue
            if current_section == "moleculetype":
                data = stripped.split(";")[0].strip()
                if data:
                    parts = data.split()
                    if parts:
                        current_mol = parts[0]
                        mol_atom_counts.setdefault(current_mol, 0)
                continue
            if current_section == "atoms" and current_mol:
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                parts = data.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    mol_atom_counts[current_mol] += 1
                continue
            if current_section == "molecules":
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                parts = data.split()
                if len(parts) >= 2:
                    try:
                        molecules_counts[parts[0]] = molecules_counts.get(parts[0], 0) + int(parts[1])
                    except ValueError:
                        continue
        if not molecules_counts or not mol_atom_counts:
            return None
        total = 0
        for mol_name, count in molecules_counts.items():
            atom_count = mol_atom_counts.get(mol_name)
            if atom_count is None:
                return None
            total += atom_count * count
        return total
    
    def _update_system_top_includes(
        self,
        existing_top: Path,
        combined_atomtypes_path: Path,
        molecule_itp_files: List[Path],
        defaults: Optional["DefaultsEntry"] = None,
        top_dir: Optional[Path] = None,
        allow_default_defaults: bool = False,
        ctx: Optional["PipelineContext"] = None,
        extra_includes: Optional[List[str]] = None,
    ) -> str:
        """
        Update existing system.top using sentinel-block based replacement.
        
        v4 hardening: Non-destructive update that preserves all content outside
        the sentinel block. Maintains #define, #ifdef, extra includes, and custom
        sections from the original file.
        
        Sentinel markers:
        - ; >>> sanitizer managed block begin
        - ; <<< sanitizer managed block end
        
        If markers exist: replace ONLY content between them.
        If markers absent: insert block after forcefield include or after banner.
        
        Args:
            existing_top: Path to existing system.top
            combined_atomtypes_path: Path to combined atomtypes ITP
            molecule_itp_files: List of molecule ITP files
            defaults: DefaultsEntry for [ defaults ] section
            top_dir: Base directory for relative path calculation
            allow_default_defaults: Allow fallback defaults if none provided
            
        Returns:
            Updated system.top content
        """
        content = self.safe_read_text(existing_top, strict=True)
        lines = content.splitlines()
        self._system_top_update_status = {
            "mode": "existing_top_update",
            "degraded": False,
        }

        strict_top_update = (
            getattr(ctx, "strict_top_update", None)
            if ctx is not None
            else None
        )
        if strict_top_update is None:
            strict_top_update = getattr(ctx, "strict_include_resolution", False) if ctx else False

        # Check sentinel markers and validate pairing before any rewrite.
        begin_marker = SANITIZER_BLOCK_BEGIN
        end_marker = SANITIZER_BLOCK_END

        stack: List[int] = []
        pairs: List[Tuple[int, int]] = []
        unmatched_end: List[int] = []
        nested_begin: List[int] = []
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped == begin_marker:
                if stack:
                    nested_begin.append(i)
                stack.append(i)
            elif stripped == end_marker:
                if stack:
                    start = stack.pop()
                    pairs.append((start, i))
                else:
                    unmatched_end.append(i)

        unmatched_begin = stack[:]
        if unmatched_begin or unmatched_end or nested_begin:
            msg = (
                "Malformed sanitizer managed markers in system.top.\n"
                f"  Unmatched begin markers: {len(unmatched_begin)}\n"
                f"  Unmatched end markers: {len(unmatched_end)}\n"
                f"  Nested begin markers: {len(nested_begin)}"
            )
            if strict_top_update:
                raise SanitizerError(msg)
            print(f"  [WARN] {msg}")
            print("  [WARN] Rebuilding sanitizer managed block to avoid stale system.top.")

            cleaned_lines: List[str] = []
            in_managed = False
            for line in lines:
                stripped = line.strip()
                if stripped == begin_marker:
                    in_managed = True
                    continue
                if stripped == end_marker:
                    in_managed = False
                    continue
                if in_managed:
                    continue
                cleaned_lines.append(line)

            has_defaults, has_forcefield_include = self._top_has_defaults_or_forcefield(
                lines=cleaned_lines,
                ignore_ranges=None,
            )
            emit_defaults = not (has_defaults or has_forcefield_include)
            managed_lines = self._build_sanitizer_managed_block(
                combined_atomtypes_path,
                molecule_itp_files,
                defaults,
                top_dir,
                allow_default_defaults,
                emit_defaults=emit_defaults,
                extra_includes=extra_includes,
            )
            self._system_top_update_status = {
                "mode": "degraded_rebuild_from_malformed_markers",
                "degraded": True,
                "strict_top_update": bool(strict_top_update),
                "marker_issues": {
                    "unmatched_begin": len(unmatched_begin),
                    "unmatched_end": len(unmatched_end),
                    "nested_begin": len(nested_begin),
                },
            }
            rebuilt = self._insert_managed_block_at_safe_location(
                cleaned_lines,
                managed_lines,
                begin_marker,
                end_marker,
            )
            return '\n'.join(rebuilt) + '\n'

        # Detect pre-existing defaults/forcefield include outside managed block(s).
        # If found, do not inject a duplicate [ defaults ] section.
        has_defaults, has_forcefield_include = self._top_has_defaults_or_forcefield(
            lines=lines,
            ignore_ranges=pairs if pairs else None,
        )
        emit_defaults = not (has_defaults or has_forcefield_include)

        # Build the sanitizer-managed include block
        managed_lines = self._build_sanitizer_managed_block(
            combined_atomtypes_path,
            molecule_itp_files,
            defaults,
            top_dir,
            allow_default_defaults,
            emit_defaults=emit_defaults,
            extra_includes=extra_includes,
        )

        if len(pairs) == 1:
            begin_idx, end_idx = pairs[0]
            if begin_idx < end_idx:
                # Replace content between markers (preserve markers)
                result_lines = (
                    lines[:begin_idx + 1] +
                    managed_lines +
                    lines[end_idx:]
                )
                self._system_top_update_status = {
                    "mode": "replaced_single_managed_block",
                    "degraded": False,
                }
            else:
                # Defensive fallback: bad ordering
                result_lines = self._insert_managed_block_at_safe_location(
                    lines, managed_lines, begin_marker, end_marker
                )
                self._system_top_update_status = {
                    "mode": "fallback_insert_bad_marker_order",
                    "degraded": True,
                }
        elif len(pairs) > 1:
            msg = (
                f"Multiple sanitizer managed blocks detected ({len(pairs)}). "
                "Replacing all managed blocks with one clean block."
            )
            if strict_top_update:
                raise SanitizerError(msg)
            print(f"  [WARN] {msg}")

            # Remove ALL managed blocks deterministically, then insert one at the first block location
            pairs_sorted = sorted(pairs, key=lambda p: p[0])
            insertion_at = pairs_sorted[0][0]
            skip_indices: Set[int] = set()
            for start, end in pairs_sorted:
                skip_indices.update(range(start, end + 1))

            result_lines: List[str] = []
            inserted = False
            for idx, line in enumerate(lines):
                if not inserted and idx == insertion_at:
                    result_lines.append(begin_marker)
                    result_lines.extend(managed_lines)
                    result_lines.append(end_marker)
                    inserted = True
                if idx in skip_indices:
                    continue
                result_lines.append(line)
            self._system_top_update_status = {
                "mode": "collapsed_multiple_managed_blocks",
                "degraded": False,
                "blocks_found": len(pairs),
            }
        else:
            # Markers not found: insert at safe location
            result_lines = self._insert_managed_block_at_safe_location(
                lines, managed_lines, begin_marker, end_marker
            )
            self._system_top_update_status = {
                "mode": "inserted_new_managed_block",
                "degraded": False,
            }
        
        return '\n'.join(result_lines) + '\n'
    
    def _build_sanitizer_managed_block(
        self,
        combined_atomtypes_path: Path,
        molecule_itp_files: List[Path],
        defaults: Optional["DefaultsEntry"],
        top_dir: Optional[Path],
        allow_default_defaults: bool,
        emit_defaults: bool = True,
        extra_includes: Optional[List[str]] = None,
    ) -> List[str]:
        """Build the content that goes inside the sanitizer managed block."""
        managed_lines: List[str] = []

        # Defaults section
        if emit_defaults:
            if defaults:
                managed_lines.extend([
                    "[ defaults ]",
                    "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                    f"  {defaults.nbfunc}        {defaults.comb_rule}         {defaults.gen_pairs}       {defaults.fudge_lj}    {defaults.fudge_qq}",
                    "",
                ])
            elif allow_default_defaults:
                managed_lines.extend([
                    "[ defaults ]",
                    "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                    "  1        2         yes        0.5      0.8333  ; LJ, geometric, standard fudge",
                    "",
                ])
            else:
                raise SanitizerError(
                    "Missing [ defaults ] and allow_default_defaults is False. "
                    "Provide explicit defaults or enable --allow-default-defaults."
                )
        
        # Combined atomtypes include
        if top_dir is not None:
            combined_rel = Path(os.path.relpath(combined_atomtypes_path, top_dir))
        else:
            combined_rel = Path(combined_atomtypes_path.name)

        if emit_defaults:
            managed_lines.append("; Combined atomtypes (after defaults for correct LJ interpretation)")
        else:
            managed_lines.append("; Combined atomtypes")
        managed_lines.extend([
            f'#include "{combined_rel.as_posix()}"',
            "",
        ])

        if extra_includes:
            managed_lines.append("; Extra includes (sanitizer generated)")
            for inc in extra_includes:
                managed_lines.append(f'#include "{inc}"')
            managed_lines.append("")

        managed_lines.append("; Molecule definitions (sanitized)")
        
        # Molecule ITP includes
        for itp_path in molecule_itp_files:
            if top_dir is not None:
                rel_path = Path(os.path.relpath(itp_path, top_dir))
                managed_lines.append(f'#include "{rel_path.as_posix()}"')
            else:
                managed_lines.append(f'#include "itp_sanitized_current/{itp_path.name}"')
        
        managed_lines.append("")
        return managed_lines

    def _top_has_defaults_or_forcefield(
        self,
        lines: List[str],
        ignore_ranges: Optional[List[Tuple[int, int]]] = None,
    ) -> Tuple[bool, bool]:
        """
        Detect whether topology already defines defaults/forcefield include.

        Args:
            lines: Raw system.top lines
            ignore_ranges: Inclusive index ranges to ignore (e.g. old managed block)

        Returns:
            (has_defaults_section, has_forcefield_include)
        """
        skip_indices: Set[int] = set()
        if ignore_ranges:
            for start, end in ignore_ranges:
                skip_indices.update(range(start, end + 1))

        has_defaults = False
        has_forcefield = False
        for idx, line in enumerate(lines):
            if idx in skip_indices:
                continue
            section = parse_section_name(line)
            if section == "defaults":
                has_defaults = True
            include_target = parse_include_target(line)
            if include_target:
                target = include_target.lower()
                if "forcefield.itp" in target or ".ff/" in target:
                    has_forcefield = True
            if has_defaults and has_forcefield:
                break

        return has_defaults, has_forcefield
    
    def _insert_managed_block_at_safe_location(
        self,
        lines: List[str],
        managed_lines: List[str],
        begin_marker: str,
        end_marker: str,
    ) -> List[str]:
        """
        Insert the managed block at a safe location in system.top.
        
        Insertion priority:
        1. After forcefield include (#include "forcefield.itp" or *.ff/forcefield.itp)
        2. After leading comment banner (first non-comment line boundary)
        3. Never before #define directives
        """
        # Find insertion point
        insert_idx = 0
        forcefield_include_idx = None
        first_define_idx = None
        first_content_idx = None
        
        for i, line in enumerate(lines):
            stripped = line.strip()
            
            # Track forcefield include
            if stripped.startswith("#include"):
                include_target = parse_include_target(stripped)
                if include_target:
                    target_lower = include_target.lower()
                    if "forcefield.itp" in target_lower or ".ff/" in target_lower:
                        forcefield_include_idx = i
            
            # Track first #define
            if first_define_idx is None and stripped.startswith("#define"):
                first_define_idx = i
            
            # Track first non-comment content (section or include)
            if first_content_idx is None and stripped and not stripped.startswith(";"):
                first_content_idx = i
        
        # Determine insertion index
        if forcefield_include_idx is not None:
            # Insert after forcefield include
            insert_idx = forcefield_include_idx + 1
        elif first_define_idx is not None:
            # Insert after all leading defines
            # Find last consecutive define from first_define_idx
            last_define = first_define_idx
            for i in range(first_define_idx, len(lines)):
                stripped = lines[i].strip()
                if stripped.startswith("#define") or stripped.startswith("#") or not stripped:
                    last_define = i
                else:
                    break
            insert_idx = last_define + 1
        elif first_content_idx is not None:
            # Insert before first content but after comments
            insert_idx = first_content_idx
        else:
            # Empty or comment-only file: insert at end
            insert_idx = len(lines)
        
        # Build result with markers
        result_lines = lines[:insert_idx]
        result_lines.append("")  # Blank line before block
        result_lines.append(begin_marker)
        result_lines.extend(managed_lines)
        result_lines.append(end_marker)
        result_lines.append("")  # Blank line after block
        result_lines.extend(lines[insert_idx:])
        
        return result_lines
    
    def _create_minimal_outputs(
        self,
        ctx: "PipelineContext",
        system_gromacs_dir: Path,
        allow_default_defaults: bool = False,
    ) -> None:
        """Create minimal outputs when no ITPs are found."""
        from ..itp_sanitizer import generate_system_top
        
        if not allow_default_defaults:
            raise SanitizerError(
                "No ITPs found and allow_default_defaults is False. "
                "Refusing to generate a topology with guessed defaults."
            )
        
        # Empty combined atomtypes
        combined_path = system_gromacs_dir / f"combined_atomtypes_{ctx.run_id}.itp"
        self._atomic_write_text(
            combined_path,
            "; Combined atomtypes (empty - no source ITPs found)\n"
            "[ atomtypes ]\n"
        )
        
        combined_current = system_gromacs_dir / "combined_atomtypes_current.itp"
        if combined_current.exists() or combined_current.is_symlink():
            combined_current.unlink()
        self._atomic_replace_file(combined_path, combined_current)
        
        # Empty sanitized directory
        sanitized_dir = system_gromacs_dir / f"itp_sanitized_{ctx.run_id}"
        sanitized_dir.mkdir(parents=True, exist_ok=True)
        
        sanitized_current = system_gromacs_dir / "itp_sanitized_current"
        if sanitized_current.exists() or sanitized_current.is_symlink():
            if sanitized_current.is_symlink():
                sanitized_current.unlink()
            else:
                shutil.rmtree(sanitized_current)
        shutil.copytree(sanitized_dir, sanitized_current)
        
        # Minimal system.top
        top_dir = system_gromacs_dir / "top"
        top_dir.mkdir(parents=True, exist_ok=True)
        
        system_top = top_dir / "system.top"
        self._atomic_write_text(
            system_top,
            generate_system_top(
                combined_atomtypes_path=combined_current,
                sanitized_itp_paths=[],
                system_name=ctx.system_id,
                top_dir=top_dir,
                allow_guess_defaults=True,  # Minimal outputs case permits fallback defaults
            ),
            encoding="utf-8",
        )
        
        # Record in manifest
        if ctx.manifest:
            ctx.manifest.set_sanitizer_output("status", "minimal_outputs_created")
            ctx.manifest.save()
    
    def _count_atoms_in_gro(self, gro_path: Path) -> Optional[int]:
        """
        Count atoms in a GRO file from the second line.
        
        Returns:
            Number of atoms, or None if parsing fails
        """
        if not gro_path.exists():
            return None
        
        try:
            with open(gro_path, 'r') as f:
                _title = f.readline()  # First line: title
                count_line = f.readline().strip()  # Second line: atom count
                return int(count_line)
        except (ValueError, IOError):
            return None
    
    def _validate_gro_top_consistency(
        self,
        gro_path: Path,
        system_top_path: Optional[Path],
        molecule_counts: Dict[str, int],
        molecule_itp_paths: Dict[str, Path],
        ctx: "PipelineContext",
        output_dir: Path,
    ) -> bool:
        """
        Validate that GRO file atom count matches expected from topology [molecules].
        
        This prevents downstream grompp failures due to coordinate/topology mismatch.
        
        Args:
            gro_path: Path to GRO file
            molecule_counts: Dict of molecule_name -> count from [molecules]
            molecule_itp_paths: Dict of molecule_name -> itp_path
            ctx: Pipeline context for logging
            
        Returns:
            True if consistent, False if mismatch detected (caller should fail)
        """
        gro_atoms = self._count_atoms_in_gro(gro_path)
        if gro_atoms is None:
            print(f"  [WARN] Could not read atom count from {gro_path}")
            return True  # Can't validate, allow to proceed
        
        # Option A: use grompp -pp if available to preprocess topology
        expected_atoms = None
        gmx_exe = shutil.which("gmx")
        if gmx_exe and system_top_path:
            expected_atoms = self._expected_atoms_from_grompp(
                gmx_exe=gmx_exe,
                system_top_path=system_top_path,
                gro_path=gro_path,
                output_dir=output_dir,
                ctx=ctx,
            )
        
        # Option B: fallback parser
        if expected_atoms is None:
            expected_atoms, uncertain = self._expected_atoms_from_itps(
                molecule_counts, molecule_itp_paths
            )
            if uncertain:
                msg = "Preprocessor directives detected in [ atoms ]; GRO/TOP check is uncertain."
                if ctx.strict_gro_top_check:
                    raise SanitizerError(msg)
                print(f"  [WARN] {msg} Skipping GRO/TOP check.")
                return True
            if expected_atoms is None:
                print("  [WARN] Could not determine expected atom count from ITPs.")
                return True
        
        if gro_atoms != expected_atoms:
            print(f"  [ERROR] GRO/TOP atom count mismatch!")
            print(f"         GRO file: {gro_atoms} atoms")
            print(f"         Expected from [molecules]: {expected_atoms} atoms")
            print(f"         Molecule counts: {molecule_counts}")
            print(f"         This usually means the coordinate file doesn't match the topology.")
            return False
        
        print(f"  [OK] GRO/TOP consistency check: {gro_atoms} atoms match")
        return True
    
    def _expected_atoms_from_itps(
        self,
        molecule_counts: Dict[str, int],
        molecule_itp_paths: Dict[str, Path],
    ) -> Tuple[Optional[int], bool]:
        """Compute expected atom count from ITPs (fallback, preprocessor-blind)."""
        expected_atoms = 0
        uncertain = False
        unknown = []
        for mol_name, mol_count in molecule_counts.items():
            if mol_count <= 0:
                continue
            itp_path = molecule_itp_paths.get(mol_name)
            if not itp_path:
                unknown.append(mol_name)
                continue
            count, has_preproc = self._count_atoms_in_itp(itp_path, mol_name)
            if has_preproc:
                uncertain = True
            if count is None:
                unknown.append(mol_name)
                continue
            expected_atoms += count * mol_count
        if unknown:
            print(f"  [WARN] Cannot validate atom counts for molecules: {unknown}")
            return None, uncertain
        return expected_atoms, uncertain

    def _count_atoms_in_itp(self, itp_path: Path, mol_name: str) -> Tuple[Optional[int], bool]:
        """
        Count atoms in a molecule ITP file from the [atoms] section.
        
        Returns:
            Tuple of (atom_count, has_preprocessor_directives)
        """
        if not itp_path.exists():
            return None, False
        content = self.safe_read_text(itp_path, strict=False)
        if not content:
            return None, False
        current_section = ""
        current_mol = None
        atom_count = 0
        has_preproc = False
        for line in content.splitlines():
            stripped = line.strip()
            if stripped.startswith("#"):
                if current_section == "atoms":
                    has_preproc = True
                continue
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if current_section == "moleculetype":
                    current_mol = None
                continue
            if current_section == "moleculetype":
                data = stripped.split(";")[0].strip()
                if data and not data.startswith(";"):
                    parts = data.split()
                    if parts:
                        current_mol = parts[0]
                continue
            if current_section == "atoms":
                if current_mol != mol_name:
                    continue
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                parts = data.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    atom_count += 1
        return (atom_count if atom_count > 0 else None), has_preproc

__all__ = ["SanitizerStage", "SanitizerError", "GROAtom", "GROParseResult"]
