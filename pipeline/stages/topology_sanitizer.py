"""
Topology sanitizer helpers extracted from sanitizer.py.

Contains ITP/include/topology/forcefield sanitization and final output helpers.
"""

from dataclasses import dataclass, field
from collections import deque
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple
import math
import hashlib
import os
import re
import shutil
import uuid

from ..manifest import _best_effort_fsync, _best_effort_fsync_dir
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
DEFAULT_GRO_COORD_RESOLUTION_NM = 1e-3

# Hardening v6: Grompp timeout complexity estimation parameters
GROMPP_TIMEOUT_FLOOR_S = 60  # Minimum timeout for any system
GROMPP_TIMEOUT_PER_1K_LINES = 10  # +10s per 1000 lines of topology
GROMPP_TIMEOUT_PER_10_INCLUDES = 30  # +30s per 10 includes resolved
GROMPP_TIMEOUT_RETRY_MULTIPLIER = 2.0  # Retry with 2x timeout on first failure

# Moleculetype signature quantization (env-overridable)
DEFAULT_MOLECULE_SIG_QUANT_STEP = 1e-6
MOLECULE_SIG_QUANT_ENV = "GPE_MOLECULE_SIG_QUANT_STEP"
GPE_UNWRAP_LOOP_TOL_NM_ENV = "GPE_UNWRAP_LOOP_TOL_NM"
UNWRAP_LOOP_TOL_FACTOR = 10.0

FLOAT_FINDER_RE = re.compile(
    r"[+-]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[eE][+-]?\d+)?"
)

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


def _extract_floats_from_text(text: str) -> List[float]:
    """Extract floats from text, including glued tokens like '1.234-5.678'."""
    values: List[float] = []
    for match in FLOAT_FINDER_RE.finditer(text):
        token = match.group(0)
        try:
            values.append(float(token))
        except ValueError:
            continue
    return values


def _decimal_places_from_numeric_token(token: str) -> Optional[int]:
    """Infer decimal places from a numeric token; returns None when unavailable."""
    cleaned = token.strip()
    match = re.fullmatch(r"[+-]?\d+(?:\.(\d+))?(?:[eE][+-]?\d+)?", cleaned)
    if not match:
        return None
    frac = match.group(1)
    return len(frac) if frac is not None else 0


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
            if normalized == pat_upper:
                return True
            # Longer patterns: allow substring but check boundaries on ORIGINAL string.
            # Using the raw string preserves separators so boundary detection is reliable.
            boundary_re = rf"(?<![A-Z0-9]){_re.escape(pat_upper)}(?![A-Z0-9])"
            if _re.search(boundary_re, raw_upper):
                return True
    return False


@dataclass
class TopMoleculesParseResult:
    """Raw [ molecules ] parse result with uncertainty metadata."""
    ordered: List[Tuple[str, int]]
    uncertain: bool = False
    uncertainty_reasons: List[str] = field(default_factory=list)
    entries: List["TopMoleculesEntry"] = field(default_factory=list)


ATOMS_ROW_STATUS_PARSED = "parsed"
ATOMS_ROW_STATUS_UNSUPPORTED = "unsupported"
ATOMS_ROW_STATUS_INVALID = "invalid"

TOP_MOLECULES_ENTRY_STATUS_PARSED = "parsed"
TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED = "unresolved"
TOP_MOLECULES_ENTRY_STATUS_INVALID = "invalid"


@dataclass(frozen=True)
class AtomRecord:
    """Minimal IR for one supported [ atoms ] row."""
    atom_id: int
    atom_type: str
    resnr: int
    residue: str
    atom_name: str
    cgnr: Optional[int]
    charge: float
    mass: float


@dataclass
class MoleculeTypeIR:
    """Minimal IR for one moleculetype atoms table."""
    name: str
    atoms: List[AtomRecord] = field(default_factory=list)


@dataclass(frozen=True)
class TopMoleculesEntry:
    """Minimal IR for one [ molecules ] entry."""
    name: str
    count_token: str
    count_value: Optional[int]


@dataclass(frozen=True)
class AtomRowParseResult:
    """Explicit parse result for one [ atoms ] row."""
    status: str
    record: Optional[AtomRecord] = None
    reason: Optional[str] = None


@dataclass(frozen=True)
class TopMoleculesEntryParseResult:
    """Explicit parse result for one [ molecules ] entry."""
    status: str
    entry: Optional[TopMoleculesEntry] = None
    reason: Optional[str] = None


def _tokenize_atoms_row(data: str) -> List[str]:
    """Split an [ atoms ] row and explode glued numeric tokens when exact."""
    tokens: List[str] = []
    for token in data.split():
        numeric_chunks = FLOAT_FINDER_RE.findall(token)
        if len(numeric_chunks) >= 2 and "".join(numeric_chunks) == token:
            tokens.extend(numeric_chunks)
        else:
            tokens.append(token)
    return tokens


def _parse_int_token(token: str) -> Optional[int]:
    """Parse an integer token, returning None on failure."""
    try:
        return int(token)
    except ValueError:
        return None


def _parse_float_token(token: str) -> Optional[float]:
    """Parse a float token, returning None on failure."""
    try:
        return float(token)
    except ValueError:
        return None


def parse_atoms_row(line: str) -> AtomRowParseResult:
    """Parse a supported [ atoms ] row into AtomRecord."""
    data = strip_gmx_comment(line).strip()
    if not data:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_INVALID,
            reason="empty [ atoms ] row",
        )

    raw_tokens = data.split()
    if len(raw_tokens) < 7:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_INVALID,
            reason=f"too few raw tokens for [ atoms ] row ({len(raw_tokens)})",
        )

    tokens = _tokenize_atoms_row(data)
    if len(tokens) > 8:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_UNSUPPORTED,
            reason=f"unsupported [ atoms ] row width ({len(tokens)})",
        )
    if len(raw_tokens) == 7 and len(tokens) == 7 and _parse_int_token(raw_tokens[5]) is not None:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_UNSUPPORTED,
            reason="ambiguous raw 7-column [ atoms ] row with integer-like sixth token",
        )

    atom_id = _parse_int_token(tokens[0])
    resnr = _parse_int_token(tokens[2])
    if atom_id is None or resnr is None:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_INVALID,
            reason="atom id/resnr must be integers",
        )

    charge = _parse_float_token(tokens[-2])
    mass = _parse_float_token(tokens[-1])
    if charge is None or mass is None:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_INVALID,
            reason="charge/mass must be numeric",
        )
    if mass < 0:
        return AtomRowParseResult(
            status=ATOMS_ROW_STATUS_INVALID,
            reason="mass must be non-negative",
        )

    cgnr: Optional[int] = None
    if len(tokens) == 8:
        cgnr = _parse_int_token(tokens[5])
        if cgnr is None:
            return AtomRowParseResult(
                status=ATOMS_ROW_STATUS_INVALID,
                reason="cgnr must be an integer in 8-column [ atoms ] rows",
            )

    return AtomRowParseResult(
        status=ATOMS_ROW_STATUS_PARSED,
        record=AtomRecord(
            atom_id=atom_id,
            atom_type=tokens[1],
            resnr=resnr,
            residue=tokens[3],
            atom_name=tokens[4],
            cgnr=cgnr,
            charge=charge,
            mass=mass,
        ),
    )


def parse_top_molecules_entry(line: str) -> TopMoleculesEntryParseResult:
    """Parse one [ molecules ] entry without applying caller policy."""
    data = strip_gmx_comment(line).strip()
    if not data:
        return TopMoleculesEntryParseResult(
            status=TOP_MOLECULES_ENTRY_STATUS_INVALID,
            reason="empty [ molecules ] row",
        )
    if data.startswith("#"):
        return TopMoleculesEntryParseResult(
            status=TOP_MOLECULES_ENTRY_STATUS_INVALID,
            reason="preprocessor directives are handled by the caller",
        )

    parts = data.split()
    if len(parts) != 2:
        return TopMoleculesEntryParseResult(
            status=TOP_MOLECULES_ENTRY_STATUS_INVALID,
            reason=f"expected 2 fields in [ molecules ] row, got {len(parts)}",
        )

    name, count_token = parts
    try:
        count_value = int(count_token)
    except ValueError:
        return TopMoleculesEntryParseResult(
            status=TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED,
            entry=TopMoleculesEntry(
                name=name,
                count_token=count_token,
                count_value=None,
            ),
            reason=f"unresolved [ molecules ] count token '{count_token}'",
        )

    return TopMoleculesEntryParseResult(
        status=TOP_MOLECULES_ENTRY_STATUS_PARSED,
        entry=TopMoleculesEntry(
            name=name,
            count_token=count_token,
            count_value=count_value,
        ),
    )



class TopologySanitizerMixin:
    """Topology/ITP/include sanitization mixin."""
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

    @staticmethod
    def _strict_charge_fix_mode(ctx: Optional["PipelineContext"]) -> bool:
        """Return True when charge-fix related strict mode is enabled."""
        return bool(ctx and getattr(ctx, "strict_charge_neutrality", False))

    @staticmethod
    def _strict_moleculetype_conflict_mode(ctx: Optional["PipelineContext"]) -> bool:
        """Return True when same-name moleculetype conflicts must fail closed."""
        return bool(
            ctx
            and (
                getattr(ctx, "strict_forcefield_consistency", False)
                or getattr(ctx, "strict_include_resolution", False)
                or getattr(ctx, "strict_charge_neutrality", False)
            )
        )

    @staticmethod
    def _is_polymer_like_molecule_name(name: str) -> bool:
        """Conservative name-based polymer marker for charge-fix protection."""
        name_lower = name.lower()
        return any(
            token in name_lower
            for token in ("polymer", "chain", "frag", "oligo", "peo", "peg")
        )

    def _charge_fix_extra_protected_patterns(
        self,
        ctx: Optional["PipelineContext"],
    ) -> frozenset[str]:
        """Return user-configured extra protected residue/molecule patterns."""
        raw = getattr(ctx, "charge_fix_protect_resnames", None) if ctx is not None else None
        if not raw:
            return frozenset()
        return frozenset(token.strip() for token in raw.split(",") if token.strip())

    def _derive_include_search_paths(
        self,
        current_file: Path,
        include_dirs: Optional[List[Path]] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> List[Path]:
        """
        Rebuild include search paths for the currently scanned file.

        Recursive include parsing must derive search order from the current include
        context, not from the original parent file that started the walk.
        """
        if ctx is not None:
            paths = list(self._build_include_search_paths(ctx, current_file))
            for candidate in include_dirs or []:
                if candidate not in paths:
                    paths.append(candidate)
            return paths

        paths: List[Path] = [current_file.parent]
        for candidate in include_dirs or []:
            if candidate not in paths:
                paths.append(candidate)
        return paths

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
            strict_conflict = self._strict_moleculetype_conflict_mode(ctx)
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

        if include_dirs is None and ctx is None:
            # Fallback is intentionally narrow: call-sites should pass ctx/include_dirs
            # when complete include search context exists.
            print(
                "  [WARN] Include resolution context unavailable while scanning "
                f"{itp_path}. Falling back to local directory only."
            )
        include_dirs = self._derive_include_search_paths(
            itp_path,
            include_dirs=include_dirs,
            ctx=ctx,
        )
        
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
        strict_compare = self._strict_moleculetype_conflict_mode(ctx)
        sig1 = self._compute_moleculetype_signature(itp1, mol_name, ctx=ctx)
        sig2 = self._compute_moleculetype_signature(itp2, mol_name, ctx=ctx)
        
        if sig1 is None or sig2 is None:
            msg = (
                f"Unable to compute deterministic signature for same-name moleculetype "
                f"'{mol_name}' while comparing:\n"
                f"  - {itp1}\n"
                f"  - {itp2}"
            )
            if strict_compare:
                raise SanitizerError(msg)
            print(f"  [WARN][DEGRADED] {msg}. Treating definitions as conflicting.")
            return True
        if sig1 != sig2:
            return True

        # Secondary guard: signatures can quantize equal; compare per-atom charges.
        quant_step = self._molecule_sig_quant_step()
        charge_tol = max(quant_step / 2.0, 0.0)
        charges1, uncertain1 = self._extract_moleculetype_atom_charges(
            itp1, mol_name, ctx=ctx
        )
        charges2, uncertain2 = self._extract_moleculetype_atom_charges(
            itp2, mol_name, ctx=ctx
        )
        if uncertain1 or uncertain2:
            msg = (
                f"Uncertain same-name moleculetype comparison for '{mol_name}': "
                "charge extraction encountered preprocessing or ambiguous formatting.\n"
                f"  - {itp1}\n"
                f"  - {itp2}"
            )
            if strict_compare:
                raise SanitizerError(msg)
            print(f"  [WARN][DEGRADED] {msg}. Treating definitions as conflicting.")
            return True
        if charges1 is None or charges2 is None:
            msg = (
                f"Incomplete charge extraction for same-name moleculetype '{mol_name}' "
                f"while comparing:\n"
                f"  - {itp1}\n"
                f"  - {itp2}"
            )
            if strict_compare:
                raise SanitizerError(msg)
            print(f"  [WARN][DEGRADED] {msg}. Treating definitions as conflicting.")
            return True
        if len(charges1) != len(charges2):
            return True
        for q1, q2 in zip(charges1, charges2):
            if abs(q1 - q2) > charge_tol:
                return True
        return False

    @staticmethod
    def _quantize_moleculetype_charge(raw_q: float, quant_step: float) -> float:
        """
        Quantize charge deterministically for signature hashing.

        Decimal arithmetic avoids platform-dependent binary midpoint/tie behavior.
        """
        if not math.isfinite(raw_q) or not math.isfinite(quant_step) or quant_step <= 0:
            return raw_q
        try:
            step_dec = Decimal(str(quant_step))
            raw_dec = Decimal(str(raw_q))
            scaled = raw_dec / step_dec
            quantized = scaled.to_integral_value(rounding=ROUND_HALF_EVEN) * step_dec
            value = float(quantized)
        except (InvalidOperation, OverflowError, ValueError, ZeroDivisionError):
            return raw_q
        return 0.0 if value == 0.0 else value

    def _molecule_sig_quant_step(self) -> float:
        """Resolve moleculetype charge quantization step from env or default."""
        raw = os.environ.get(MOLECULE_SIG_QUANT_ENV, "").strip()
        if not raw:
            return DEFAULT_MOLECULE_SIG_QUANT_STEP
        try:
            value = float(raw)
        except ValueError:
            print(
                f"  [WARN] Invalid {MOLECULE_SIG_QUANT_ENV}={raw!r}; "
                f"using default {DEFAULT_MOLECULE_SIG_QUANT_STEP:g}"
            )
            return DEFAULT_MOLECULE_SIG_QUANT_STEP
        if not math.isfinite(value) or value <= 0:
            print(
                f"  [WARN] Non-positive/non-finite {MOLECULE_SIG_QUANT_ENV}={raw!r}; "
                f"using default {DEFAULT_MOLECULE_SIG_QUANT_STEP:g}"
            )
            return DEFAULT_MOLECULE_SIG_QUANT_STEP
        return value

    @staticmethod
    def _is_int_like(value: float, tol: float = 1e-9) -> bool:
        """Return True when value is numerically close to an integer."""
        return abs(value - round(value)) <= tol

    def _parse_moleculetype_ir(
        self,
        itp_path: Path,
        mol_name: str,
        ctx: Optional["PipelineContext"] = None,
    ) -> Tuple[Optional[MoleculeTypeIR], bool]:
        """Parse one target moleculetype into the minimal IR used in Stage 1."""
        if not itp_path.exists():
            return None, False
        strict_reads = bool(ctx and ctx.strict_include_resolution)
        content = self.safe_read_text(itp_path, strict=strict_reads)
        if not content:
            return None, False

        current_section = ""
        current_mol = None
        in_target_mol = False
        parse_uncertain = False
        atoms: List[AtomRecord] = []

        for line in content.splitlines():
            stripped = line.strip()
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if current_section == "moleculetype":
                    current_mol = None
                    in_target_mol = False
                continue

            data = strip_gmx_comment(stripped).strip()
            if not data:
                continue

            if current_section == "moleculetype":
                parts = data.split()
                if parts:
                    current_mol = parts[0]
                    in_target_mol = (current_mol == mol_name)
                continue

            if not in_target_mol:
                continue

            if data.startswith("#"):
                parse_uncertain = True
                continue
            if current_section != "atoms":
                continue

            row_result = parse_atoms_row(data)
            if row_result.status != ATOMS_ROW_STATUS_PARSED or row_result.record is None:
                parse_uncertain = True
                continue
            atoms.append(row_result.record)

        if not atoms:
            return None, parse_uncertain
        return MoleculeTypeIR(name=mol_name, atoms=atoms), parse_uncertain

    def _extract_atoms_row_charge_value(self, data: str) -> Tuple[Optional[float], bool]:
        """
        Extract charge from an [ atoms ] row.

        Compatibility wrapper around the Stage 1 row parser.
        """
        row_result = parse_atoms_row(data)
        if row_result.status == ATOMS_ROW_STATUS_PARSED and row_result.record is not None:
            return row_result.record.charge, False
        return None, row_result.status == ATOMS_ROW_STATUS_UNSUPPORTED

    def _extract_moleculetype_atom_charges(
        self,
        itp_path: Path,
        mol_name: str,
        ctx: Optional["PipelineContext"] = None,
    ) -> Tuple[Optional[List[float]], bool]:
        """Extract per-atom charges for one moleculetype with uncertainty signaling."""
        mol_ir, charges_uncertain = self._parse_moleculetype_ir(
            itp_path,
            mol_name,
            ctx=ctx,
        )
        if mol_ir is None or not mol_ir.atoms:
            return None, charges_uncertain
        return [atom.charge for atom in mol_ir.atoms], charges_uncertain
    
    def _classify_ionic_moleculetypes(
        self,
        molecule_itp_paths: Dict[str, Path],
        ctx: Optional["PipelineContext"] = None,
    ) -> Tuple[Dict[str, str], Dict[str, Any]]:
        """
        Classify moleculetypes for charge correction protection.
        
        Hardening v8: Explicitly reports degraded/unknown classifications instead
        of silently treating uncertain parsing as neutral.
        
        Classification rules (priority order):
        1. Pattern match against PROTECTED_ION_PATTERNS -> "protected_ion"
        2. Pattern match against PROTECTED_SOLVENT_PATTERNS -> "protected_solvent"
        3. Pattern match against ctx.charge_fix_protect_resnames -> "protected_configured"
        3. n_atoms == 1 and |q| > 0 -> "monatomic_ion"
        4. abs(round(q_total)) >= 1 and close to integer -> "ionic"
        5. Name contains polymer/chain patterns -> "protected_polymer"
        6. Otherwise -> "neutral"
        
        Args:
            molecule_itp_paths: Dict of molecule_name -> itp_path
            ctx: PipelineContext for tolerance
            
        Returns:
            Tuple of:
            - Dict of molecule_name -> classification
            - Classification audit with degraded/uncertain details
        """
        from ..charge_neutrality import ChargeNeutralityError, ChargeParser

        tol = ctx.charge_neutrality_tol if ctx else 1e-6
        classifications: Dict[str, str] = {}
        strict_mode = self._strict_charge_fix_mode(ctx)
        configured_patterns = self._charge_fix_extra_protected_patterns(ctx)
        classification_audit: Dict[str, Any] = {
            "configured_protected_patterns": sorted(configured_patterns, key=str.casefold),
            "details": {},
            "degraded_molecules": {},
        }

        def _record(
            mol_name: str,
            itp_path: Path,
            classification: str,
            *,
            status: str = "exact",
            source: str,
            reason: Optional[str] = None,
            atom_count: Optional[int] = None,
            total_charge: Optional[float] = None,
            has_preprocessor: bool = False,
            unparsed_atoms_line_count: int = 0,
        ) -> None:
            detail = {
                "classification": classification,
                "status": status,
                "source": source,
                "reason": reason,
                "itp_path": str(itp_path),
                "atom_count": atom_count,
                "total_charge": total_charge,
                "has_preprocessor": has_preprocessor,
                "unparsed_atoms_line_count": unparsed_atoms_line_count,
            }
            classifications[mol_name] = classification
            classification_audit["details"][mol_name] = detail
            if status != "exact":
                classification_audit["degraded_molecules"][mol_name] = detail

        for mol_name, itp_path in sorted(
            molecule_itp_paths.items(),
            key=lambda item: item[0].casefold(),
        ):
            looks_polymer = self._is_polymer_like_molecule_name(mol_name)

            if _matches_protected_pattern(mol_name, PROTECTED_ION_PATTERNS):
                _record(
                    mol_name,
                    itp_path,
                    "protected_ion",
                    source="protected_ion_pattern",
                )
                continue

            if _matches_protected_pattern(mol_name, PROTECTED_SOLVENT_PATTERNS):
                _record(
                    mol_name,
                    itp_path,
                    "protected_solvent",
                    source="protected_solvent_pattern",
                )
                continue

            if configured_patterns and _matches_protected_pattern(mol_name, configured_patterns):
                _record(
                    mol_name,
                    itp_path,
                    "protected_configured",
                    source="charge_fix_protect_resnames",
                    reason="Matched user-configured protected residue/molecule pattern.",
                )
                continue

            try:
                mol_charge = ChargeParser.get_molecule_charge(
                    itp_path,
                    name=mol_name,
                    strict=False,
                    allow_unparsed_atoms_lines=True,
                )

                issue_bits: List[str] = []
                if mol_charge.atom_count <= 0 or not mol_charge.charges:
                    issue_bits.append("missing or empty [ atoms ] charge data")
                if mol_charge.has_preprocessor:
                    issue_bits.append("preprocessor directives in [ atoms ]")
                if mol_charge.unparsed_atoms_lines:
                    issue_bits.append(
                        f"{len(mol_charge.unparsed_atoms_lines)} unparseable [ atoms ] data line(s)"
                    )

                if issue_bits:
                    classification = "protected_polymer" if looks_polymer else "unknown"
                    msg = (
                        f"Charge-fix classification for active molecule '{mol_name}' is degraded: "
                        + "; ".join(issue_bits)
                    )
                    if strict_mode:
                        raise SanitizerError(msg)
                    print(f"  [WARN][DEGRADED] {msg}")
                    _record(
                        mol_name,
                        itp_path,
                        classification,
                        status="degraded",
                        source="charge_parser_degraded",
                        reason="; ".join(issue_bits),
                        atom_count=mol_charge.atom_count,
                        total_charge=mol_charge.total_charge,
                        has_preprocessor=mol_charge.has_preprocessor,
                        unparsed_atoms_line_count=len(mol_charge.unparsed_atoms_lines or []),
                    )
                    continue

                n_atoms = mol_charge.atom_count
                q_total = mol_charge.total_charge

                if n_atoms == 1 and abs(q_total) > tol:
                    _record(
                        mol_name,
                        itp_path,
                        "monatomic_ion",
                        source="charge_parser",
                        atom_count=n_atoms,
                        total_charge=q_total,
                    )
                    continue

                rounded_q = round(q_total)
                if abs(rounded_q) >= 1 and abs(q_total - rounded_q) < tol:
                    _record(
                        mol_name,
                        itp_path,
                        "ionic",
                        source="charge_parser",
                        atom_count=n_atoms,
                        total_charge=q_total,
                    )
                    continue

                if looks_polymer:
                    _record(
                        mol_name,
                        itp_path,
                        "protected_polymer",
                        source="polymer_name_pattern",
                        atom_count=n_atoms,
                        total_charge=q_total,
                    )
                    continue

                _record(
                    mol_name,
                    itp_path,
                    "neutral",
                    source="charge_parser",
                    atom_count=n_atoms,
                    total_charge=q_total,
                )

            except ChargeNeutralityError as exc:
                classification = "protected_polymer" if looks_polymer else "unknown"
                msg = (
                    f"Charge-fix classification for active molecule '{mol_name}' is degraded: "
                    f"{str(exc).splitlines()[0]}"
                )
                if strict_mode:
                    raise SanitizerError(msg)
                print(f"  [WARN][DEGRADED] {msg}")
                _record(
                    mol_name,
                    itp_path,
                    classification,
                    status="degraded",
                    source="charge_parser_error",
                    reason=str(exc).splitlines()[0],
                )
                continue

            except Exception as exc:
                classification = "protected_polymer" if looks_polymer else "unknown"
                msg = (
                    f"Charge-fix classification for active molecule '{mol_name}' failed unexpectedly: "
                    f"{exc}"
                )
                if strict_mode:
                    raise SanitizerError(msg)
                print(f"  [WARN][DEGRADED] {msg}")
                _record(
                    mol_name,
                    itp_path,
                    classification,
                    status="degraded",
                    source="unexpected_error",
                    reason=str(exc),
                )
        
        # Log classification summary
        prot_ions = [k for k, v in classifications.items() if v in ("protected_ion", "ionic", "monatomic_ion")]
        prot_solvents = [k for k, v in classifications.items() if v == "protected_solvent"]
        prot_configured = [k for k, v in classifications.items() if v == "protected_configured"]
        prot_polymers = [k for k, v in classifications.items() if v == "protected_polymer"]
        neutral = [k for k, v in classifications.items() if v == "neutral"]
        degraded = sorted(classification_audit["degraded_molecules"], key=str.casefold)
        
        if prot_ions:
            print(f"  [CHARGE-FIX] Protected ions: {prot_ions}")
        if prot_solvents:
            print(f"  [CHARGE-FIX] Protected solvents: {prot_solvents}")
        if prot_configured:
            print(f"  [CHARGE-FIX] User-protected molecules: {prot_configured}")
        if prot_polymers:
            print(f"  [CHARGE-FIX] Protected polymers: {prot_polymers}")
        if neutral:
            print(f"  [CHARGE-FIX] Correction candidates (neutral): {neutral}")
        if degraded:
            print(f"  [WARN][DEGRADED] Charge-fix classification degraded for: {degraded}")

        return classifications, classification_audit

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

        protected_classes = {
            "protected_ion",
            "ionic",
            "monatomic_ion",
            "protected_solvent",
            "protected_configured",
        }
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
        classification_audit: Optional[Dict[str, Any]] = None,
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
        classification_details = (
            classification_audit.get("details", {})
            if isinstance(classification_audit, dict)
            else {}
        )
        detail_by_casefold = {
            name.casefold(): detail
            for name, detail in classification_details.items()
        }
        degraded_molecules = {
            name: detail
            for name, detail in classification_details.items()
            if str(detail.get("status", "exact")) != "exact"
        }
        
        target_mol = correction_result.target_molecule
        target_mol_cf = target_mol.casefold() if target_mol else None
        class_by_casefold = {name.casefold(): cls for name, cls in classifications.items()}
        target_classification = class_by_casefold.get(target_mol_cf or "", "unknown")
        target_detail = detail_by_casefold.get(target_mol_cf or "", {})
        target_classification_status = str(
            target_detail.get(
                "status",
                "degraded" if target_classification == "unknown" else "exact",
            )
        )
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
                "protected_configured",
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
            "protected_configured": (
                False,
                "--charge-fix-target-allowlist",
            ),
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
            if (
                not target_from_allowlist
                and (
                    classification == "unknown"
                    or target_classification_status != "exact"
                )
            ):
                violations.append(
                    f"Charge correction targeted '{target_mol}' with {target_classification_status} "
                    f"classification (classification={classification}). "
                    "Use --charge-fix-target-allowlist only when you explicitly accept this degraded target."
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
            "target_classification_status": target_classification_status,
            "all_classifications": classifications,
            "classification_audit": classification_audit,
            "degraded_molecules": degraded_molecules,
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
        quant_step = self._molecule_sig_quant_step()
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
                    row_result = parse_atoms_row(data)
                    if row_result.status != ATOMS_ROW_STATUS_PARSED or row_result.record is None:
                        unknown_key = "atoms:charge_unknown"
                        if row_result.status == ATOMS_ROW_STATUS_UNSUPPORTED:
                            unknown_key = "atoms:charge_uncertain"
                        sections_present.add(unknown_key)
                        section_line_counts[unknown_key] = section_line_counts.get(unknown_key, 0) + 1
                        hasher.update(f"{unknown_key}|{' '.join(parts[:5])}\n".encode("utf-8"))
                        continue
                    raw_q = row_result.record.charge

                    qmin = raw_q if qmin is None else min(qmin, raw_q)
                    qmax = raw_q if qmax is None else max(qmax, raw_q)
                    sum_abs_q += abs(raw_q)
                    sum_q2 += raw_q * raw_q
                    quant_q = self._quantize_moleculetype_charge(raw_q, quant_step)
                    canonical = (
                        f"{row_result.record.atom_type}|"
                        f"{row_result.record.atom_name}|"
                        f"{quant_q:+.10f}\n"
                    )
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
            f"charge_quant_step={quant_step:.12g};"
            f"sections={','.join(sorted(sections_present))};"
            f"section_counts={','.join(f'{k}:{section_line_counts[k]}' for k in sorted(section_line_counts))}"
        )

    def _get_ordered_molecules_from_top(self, top_path: Path) -> TopMoleculesParseResult:
        """
        Parse [molecules] section preserving order, collapsing duplicates.
        
        Task D: GROMACS allows the same molecule name to appear on multiple
        lines in [molecules]. We collapse duplicates by summing counts into
        the first occurrence to avoid breaking include ordering.
        """
        if not top_path.exists():
            return TopMoleculesParseResult(ordered=[])
        # Use list for order + dict for fast lookup of first occurrence index
        ordered: List[Tuple[str, int]] = []
        entries: List[TopMoleculesEntry] = []
        seen_indices: Dict[str, int] = {}  # name -> index in ordered list
        in_molecules = False
        uncertain_reasons: List[str] = []
        conditional_directives = frozenset(
            {"#if", "#ifdef", "#ifndef", "#elif", "#else", "#endif"}
        )
        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return TopMoleculesParseResult(ordered=[])
        for line_no, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            section = parse_section_name(line)
            if section is not None:
                in_molecules = (section == "molecules")
                continue
            if not in_molecules:
                continue
            if stripped.startswith("#"):
                directive = stripped.split()[0].lower()
                if directive in conditional_directives:
                    uncertain_reasons.append(
                        f"{directive} in [molecules] at {top_path.name}:{line_no}"
                    )
                continue
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            data = stripped.split(";")[0].strip()
            if not data:
                continue
            entry_result = parse_top_molecules_entry(data)
            if entry_result.entry is not None:
                entries.append(entry_result.entry)
            if entry_result.status == TOP_MOLECULES_ENTRY_STATUS_PARSED and entry_result.entry is not None:
                name = entry_result.entry.name
                count = entry_result.entry.count_value or 0
                if name in seen_indices:
                    # Collapse: add count to existing entry
                    idx = seen_indices[name]
                    old_name, old_count = ordered[idx]
                    ordered[idx] = (old_name, old_count + count)
                else:
                    # First occurrence: add new entry
                    seen_indices[name] = len(ordered)
                    ordered.append((name, count))
                continue
            if entry_result.status == TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED and entry_result.entry is not None:
                uncertain_reasons.append(
                    f"unresolved [molecules] count token '{entry_result.entry.count_token}' "
                    f"for '{entry_result.entry.name}' at {top_path.name}:{line_no}"
                )
                continue
            uncertain_reasons.append(
                f"invalid [molecules] row at {top_path.name}:{line_no}: "
                f"{entry_result.reason or data}"
            )
        return TopMoleculesParseResult(
            ordered=ordered,
            entries=entries,
            uncertain=bool(uncertain_reasons),
            uncertainty_reasons=uncertain_reasons,
        )
    
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
        ctx: Optional["PipelineContext"] = None,
    ) -> Dict[str, Path]:
        """Map molecule names to sanitized ITPs (by moleculetype name)."""
        if not sanitized_current_dir.exists() or not needed_molecule_names:
            return {}
        name_map = {n.casefold(): n for n in needed_molecule_names}
        mapping: Dict[str, Path] = {}
        for itp_path in self._sorted_paths_casefold(sanitized_current_dir.glob("*.itp")):
            for mol_name in self._get_moleculetype_names(itp_path, ctx=ctx):
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

        Conservative parser:
        - Requires ptype in the canonical position (third token from tail)
        - Requires tail shape: ... mass charge ptype sigma epsilon
        - Rejects ambiguous rows rather than guessing shifted columns
        """
        parts = data.split()
        if len(parts) < 6:
            raise ValueError(f"Too few fields for [ atomtypes ] row: '{data}'")

        atomtype_name = parts[0]
        lj_tail_floats = _extract_floats_from_text(" ".join(parts[-2:]))
        if len(lj_tail_floats) != 2:
            raise ValueError(
                f"Ambiguous LJ params in atomtype '{atomtype_name}' at "
                f"{source_name}:{line_no}. Expected exactly two float values in LJ tail, "
                f"got {len(lj_tail_floats)}. Raw: {data}"
            )
        p1, p2 = lj_tail_floats

        ptype = parts[-3].upper()
        if ptype not in PTYPE_TOKENS:
            raise ValueError(
                f"Ambiguous/missing ptype token at {source_name}:{line_no}. "
                f"Expected one of {sorted(PTYPE_TOKENS)} in column -3; got '{parts[-3]}'. "
                f"Raw: {data}"
            )

        extra_ptypes = [
            token for token in parts[:-3] if token.upper() in PTYPE_TOKENS
        ]
        if extra_ptypes:
            raise ValueError(
                f"Ambiguous row with multiple ptype-like tokens at {source_name}:{line_no}: "
                f"{extra_ptypes + [ptype]}. Raw: {data}"
            )

        charge_fields = _extract_floats_from_text(parts[-4]) if len(parts) >= 4 else []
        mass_fields = _extract_floats_from_text(parts[-5]) if len(parts) >= 5 else []
        if len(charge_fields) != 1 or len(mass_fields) != 1:
            raise ValueError(
                f"Malformed mass/charge columns before ptype in {source_name}:{line_no}. Raw: {data}"
            )

        return atomtype_name, ptype, p1, p2, None

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
            except ValueError as exc:
                msg = (
                    f"Ambiguous/malformed [ atomtypes ] LJ row in "
                    f"{combined_atomtypes_path.name}:{line_no}: {data}\n"
                    f"  Reason: {exc}"
                )
                if policy == "error":
                    raise SanitizerError(msg)
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
    ) -> Optional[Path]:
        """Ensure corrected ITP matches the include path used by system.top."""
        if not correction_result.corrected or not correction_result.patched_itp_path:
            return None
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
        run_token = getattr(correction_result, "run_id", None) or "run"
        backup_path = sanitized_dir / f"{original_path.stem}_pre_patch_{run_token}{original_path.suffix}"
        if correction_result.patched_itp_path.resolve() != original_path.resolve():
            # Required order for crash-safe semantics:
            # 1) backup original -> versioned backup path
            # 2) replace original with patched content
            self._atomic_replace_file(original_path, backup_path)
            self._atomic_replace_file(correction_result.patched_itp_path, original_path)
        # Remove stray corrected filename to avoid confusion
        if correction_result.patched_itp_path.exists() and correction_result.patched_itp_path.resolve() != original_path.resolve():
            correction_result.patched_itp_path.unlink()
        return backup_path

    def _atomic_replace_file(self, src: Path, dest: Path) -> None:
        """Atomic replace of dest with src content."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        temp_path = dest.parent / f".{dest.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
        try:
            with open(src, "rb") as src_handle, open(temp_path, "wb") as temp_handle:
                shutil.copyfileobj(src_handle, temp_handle)
                temp_handle.flush()
                _best_effort_fsync(temp_handle.fileno())
            os.replace(temp_path, dest)
            _best_effort_fsync_dir(dest.parent)
        finally:
            if temp_path.exists():
                temp_path.unlink()

    def _atomic_write_text(self, path: Path, content: str, encoding: str = "utf-8") -> None:
        """
        Write text content atomically via temp file + replace.
        
        Prevents partial files after crashes.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = path.parent / f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
        try:
            with open(temp_path, "w", encoding=encoding) as handle:
                handle.write(content)
                handle.flush()
                _best_effort_fsync(handle.fileno())
            os.replace(temp_path, path)
            _best_effort_fsync_dir(path.parent)
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
        
        Deterministic policy:
        1. If forcefield include block exists, insert immediately after its last line.
        2. Otherwise insert after leading comment/blank preamble.
        """
        insert_idx = 0
        forcefield_include_idx = None

        for i, line in enumerate(lines):
            stripped = line.strip()
            include_target = parse_include_target(stripped)
            if include_target:
                target_lower = include_target.lower()
                if "forcefield.itp" in target_lower or ".ff/" in target_lower:
                    forcefield_include_idx = i

        if forcefield_include_idx is not None:
            insert_idx = forcefield_include_idx + 1
        else:
            idx = 0
            while idx < len(lines):
                stripped = lines[idx].strip()
                if not stripped or stripped.startswith(";"):
                    idx += 1
                    continue
                break
            insert_idx = idx
        
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
    

