"""
ITP Sanitizer: Atomtype extraction, de-duplication, and conflict detection.

This module implements:
- Parsing ITP files to extract [ atomtypes ], [ atoms ], [ defaults ] and other sections
- Token-based atomtype parsing (robust across GROMACS ITP variants)
- Atomtype signature computation for conflict detection (excludes placeholder charge)
- Merging atomtypes from multiple sources with namespace collision defense
- [ defaults ] / comb-rule validation
- Generating sanitized ITP files with atomtypes blocks removed
- Optional namespace prefixing for atomtype names with safety checks
- Recursive #include resolution for dependency discovery

Per AGENT_SKILLS/06_ITP_SANITIZER.md:
- Run after staging and before any grompp
- Fail-fast on atomtype conflicts (different params for same name)
- Track source_ff for cross-forcefield contamination detection
- Allow override only with explicit --allow-override flag

Hardening (v2):
- Signature excludes charge column (often placeholder 0.0)
- Token-based parser replaces fragile regex
- Defaults validation ensures comb-rule consistency
- Prefix safety check detects implicit topology usage
- Deterministic file ordering for reproducibility
"""

import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
import shutil
from enum import Enum


class LJParamKind(Enum):
    """
    LJ parameter representation kind based on GROMACS comb-rule.
    
    - SIGMA_EPSILON: comb-rule 2 or 3, values are sigma (nm) and epsilon (kJ/mol)
    - C6_C12: comb-rule 1, values are C6 and C12 coefficients (very small magnitudes)
    """
    SIGMA_EPSILON = "sigma_epsilon"
    C6_C12 = "c6_c12"


# =============================================================================
# Constants and Tolerances (configurable via env vars)
# =============================================================================

# Mass comparison tolerances (Task B: robust abs+rel tolerance)
# Use max(abs_tol, rel_tol * max(|a|, |b|)) for mass comparisons
MASS_ABS_TOL = float(os.environ.get("ITP_MASS_ABS_TOL", "1e-3"))  # ~0.001 Da
MASS_REL_TOL = float(os.environ.get("ITP_MASS_REL_TOL", "1e-4"))  # 0.01%

# Legacy constant (for backward compat, use new tolerances in params_match)
MASS_TOLERANCE = 1e-4
LJ_TOLERANCE = 1e-6
CHARGE_TOLERANCE = 1e-4

# Relative tolerances for comb-rule aware comparison
LJ_REL_TOL_SIGMA_EPSILON = 1e-5  # 0.001% relative tolerance
LJ_REL_TOL_C6_C12 = 1e-4  # 0.01% for small C12 values
LJ_ABS_FLOOR_SIGMA = 1e-9  # absolute floor for sigma comparison
LJ_ABS_FLOOR_EPSILON = 1e-12  # absolute floor for epsilon
LJ_ABS_FLOOR_C6 = 1e-18  # absolute floor for C6
LJ_ABS_FLOOR_C12 = 1e-24  # absolute floor for C12 (can be very small)

# Sanity bounds for LJ parameters (warn on outliers) - Task L: configurable
SIGMA_MIN_NM = float(os.environ.get("ITP_SIGMA_MIN_NM", "0.01"))
SIGMA_MAX_NM = float(os.environ.get("ITP_SIGMA_MAX_NM", "2.0"))
SIGMA_BOUNDS_NM = (SIGMA_MIN_NM, SIGMA_MAX_NM)  # typical: 0.1-0.5 nm

EPSILON_MIN_KJ_MOL = float(os.environ.get("ITP_EPSILON_MIN_KJ_MOL", "0.0"))
EPSILON_MAX_KJ_MOL = float(os.environ.get("ITP_EPSILON_MAX_KJ_MOL", "1000.0"))
EPSILON_BOUNDS_KJ_MOL = (EPSILON_MIN_KJ_MOL, EPSILON_MAX_KJ_MOL)

C6_MIN = float(os.environ.get("ITP_C6_MIN", "0.0"))
C6_MAX = float(os.environ.get("ITP_C6_MAX", "1e-2"))
C6_BOUNDS = (C6_MIN, C6_MAX)  # typical ~1e-3 for C-C

C12_MIN = float(os.environ.get("ITP_C12_MIN", "0.0"))
C12_MAX = float(os.environ.get("ITP_C12_MAX", "1e-4"))
C12_BOUNDS = (C12_MIN, C12_MAX)  # typical ~1e-6 for C-C

# Valid ptype values (D=Dummy, V=Virtual are treated specially for LJ warnings)
VALID_PTYPES = {"A", "D", "S", "V"}  # Atom, Dummy, Shell, Virtual
DUMMY_PTYPES = {"D", "V"}  # ptype values where sigma=0, epsilon=0 is expected

# Type-based sections that reference atomtype names (Task H: prefix safety)
TYPE_BASED_SECTIONS = {
    "pairtypes", "nonbond_params", "bondtypes", "angletypes", "dihedraltypes",
    "constrainttypes", "cmaptypes", "virtual_sites2", "virtual_sites3",
    "virtual_sites4", "virtual_sitesn",
}

# Number of leading type-name tokens in type-based sections.
# Used when rewriting atomtype prefixes outside [ atoms ]/[ atomtypes ].
TYPE_SECTION_ATOMTYPE_TOKEN_COUNT = {
    "pairtypes": 2,
    "nonbond_params": 2,
    "bondtypes": 2,
    "constrainttypes": 2,
    "angletypes": 3,
    "dihedraltypes": 4,
    "impropertypes": 4,
    "cmaptypes": 5,
}

# Wildcard markers commonly used in type tables.
WILDCARD_TYPE_TOKENS = {"X", "x", "*"}

# File classification priorities (higher = preferred in conflicts)
FILE_PRIORITY = {
    "forcefield": 100,
    "staged": 50,
    "htpolynet": 40,
    "molecule": 30,
    "unknown": 0,
}


# =============================================================================
# Exceptions
# =============================================================================

class ItpSanitizerError(Exception):
    """Error during ITP sanitization."""
    pass


class AtomtypeConflictError(ItpSanitizerError):
    """Raised when atomtype definitions conflict."""
    pass


class NamespaceCollisionError(ItpSanitizerError):
    """Raised when same atomtype appears from different forcefields with different params."""
    pass


class DefaultsConflictError(ItpSanitizerError):
    """Raised when [ defaults ] blocks conflict."""
    pass


class PrefixUnsafeError(ItpSanitizerError):
    """Raised when atomtype prefixing would break implicit topology."""
    pass


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class DefaultsEntry:
    """
    Represents [ defaults ] section content.
    
    GROMACS format:
    ; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ
      1       2          yes        0.5      0.8333
    """
    nbfunc: int  # 1 = LJ (Lennard-Jones), 2 = Buckingham
    comb_rule: int  # 1, 2, or 3 (LJ combination rule)
    gen_pairs: str  # "yes" or "no"
    fudge_lj: float  # LJ 1-4 scaling
    fudge_qq: float  # Coulomb 1-4 scaling
    source_file: str = ""
    
    def signature(self) -> str:
        """Signature for conflict detection."""
        return f"{self.nbfunc}|{self.comb_rule}|{self.gen_pairs}|{self.fudge_lj:.4f}|{self.fudge_qq:.4f}"
    
    def matches(self, other: "DefaultsEntry") -> bool:
        """
        Check if defaults match (Task F: tolerant float comparison).
        
        Uses math.isclose for fudgeLJ/fudgeQQ to handle formatting differences
        like 0.8333 vs 0.833333 across toolchains.
        """
        if self.nbfunc != other.nbfunc:
            return False
        if self.comb_rule != other.comb_rule:
            return False
        if self.gen_pairs.lower() != other.gen_pairs.lower():
            return False
        # Use math.isclose for float comparison (rel_tol=1e-4 handles 0.8333 vs 0.833333)
        if not math.isclose(self.fudge_lj, other.fudge_lj, rel_tol=1e-4, abs_tol=1e-6):
            return False
        if not math.isclose(self.fudge_qq, other.fudge_qq, rel_tol=1e-4, abs_tol=1e-6):
            return False
        return True


def _fallback_defaults_entry() -> DefaultsEntry:
    """Return fallback defaults (used only when explicitly allowed)."""
    return DefaultsEntry(
        nbfunc=1,
        comb_rule=2,
        gen_pairs="yes",
        fudge_lj=0.5,
        fudge_qq=0.8333,
        source_file="FALLBACK_DEFAULTS",
    )


def _is_forcefield_defaults_source(
    source_file: str,
    class_map: Optional[Dict[str, str]] = None,
) -> bool:
    """
    Determine whether a defaults source is forcefield-owned.

    Priority markers:
    - explicit `forcefield.itp`
    - any path component ending with `.ff`
    - class_map classification set to `forcefield`
    """
    source_path = Path(source_file)
    try:
        resolved = source_path.resolve()
    except OSError:
        resolved = source_path

    source_text = source_path.as_posix().lower()
    if source_path.name.lower() == "forcefield.itp":
        return True
    if any(part.lower().endswith(".ff") for part in source_path.parts):
        return True
    if ".ff/" in source_text:
        return True

    if class_map:
        class_hint = (
            class_map.get(str(source_path))
            or class_map.get(str(resolved))
            or class_map.get(source_text)
            or class_map.get(resolved.as_posix())
        )
        if class_hint == "forcefield":
            return True

    return False


@dataclass
class AtomtypeEntry:
    """
    Represents an atomtype definition.
    
    GROMACS atomtype format (from [ atomtypes ] section):
    ; name  at.num  mass     charge   ptype  sigma      epsilon
    ca         6    12.01    0.0000   A      0.339967   0.359824
    
    Or older format:
    ; name  bond_type  at.num  mass  charge  ptype  sigma  epsilon
    
    Or simplified format (GAFF2/OPLS):
    ; name  mass  charge  ptype  sigma  epsilon
    
    NOTE: For comb-rule 1, last two columns are C6/C12 instead of sigma/epsilon.
    """
    name: str
    mass: float
    charge: float
    ptype: str  # Usually 'A' for atom
    sigma: float  # or C6 if comb_rule==1
    epsilon: float  # or C12 if comb_rule==1
    at_num: Optional[int] = None
    bond_type: Optional[str] = None
    source_ff: str = ""
    source_files: List[str] = field(default_factory=list)
    file_class: str = "unknown"  # forcefield, staged, htpolynet, molecule, unknown
    raw_line: str = ""
    # Comb-rule awareness fields
    lj_kind: LJParamKind = field(default=LJParamKind.SIGMA_EPSILON)
    comb_rule: Optional[int] = None  # 1, 2, or 3
    
    def signature(self, include_charge: bool = False) -> str:
        """
        Compute a signature for conflict detection.
        
        IMPORTANT: By default, signature EXCLUDES charge because in many forcefields
        the [ atomtypes ] charge is a placeholder (often 0.0) and real
        charges live in [ atoms ] sections.
        
        Args:
            include_charge: If True, include charge in signature
        
        Includes: name, ptype, mass (4 decimals), sigma (6 decimals), epsilon (6 decimals)
        """
        base = f"{self.name}|{self.ptype}|{self.mass:.4f}|{self.sigma:.6f}|{self.epsilon:.6f}"
        if include_charge:
            base += f"|{self.charge:.4f}"
        if self.comb_rule is not None:
            base += f"|comb{self.comb_rule}"
        return base
    
    def params_match(self, other: "AtomtypeEntry") -> bool:
        """
        Check if physical parameters match within tolerances.
        
        Uses comb-rule aware relative tolerances to avoid:
        - False positives from absolute tolerances on very small C12 values
        - False negatives from large absolute differences on sigma
        
        Returns False (mismatch) if comb-rules differ and cannot be compared.
        """
        # Task B: Use robust mass tolerance with both abs and rel components
        # This handles rounding differences like 12.0110 vs 12.0107 (carbon)
        mass_diff = abs(self.mass - other.mass)
        mass_max = max(abs(self.mass), abs(other.mass))
        mass_threshold = max(MASS_ABS_TOL, MASS_REL_TOL * mass_max)
        if mass_diff > mass_threshold:
            return False
        
        # Handle ptype=None (missing ptype) - count as match if both None or same value
        if self.ptype != other.ptype:
            return False
        
        # Comb-rule must match for comparison to be valid when BOTH are known
        if self.comb_rule is not None and other.comb_rule is not None:
            if self.comb_rule != other.comb_rule:
                # Can't reliably compare sigma/epsilon to C6/C12
                # Caller should fail-fast rather than guess
                # Issue A: Log diagnostic message for debugging 
                print(
                    f"  [DEBUG] Comb-rule mismatch for atomtype '{self.name}': "
                    f"comb-rule={self.comb_rule} vs {other.comb_rule} "
                    f"(files: {self.source_files[0] if self.source_files else '?'} vs "
                    f"{other.source_files[0] if other.source_files else '?'})"
                )
                return False
        elif self.comb_rule is None or other.comb_rule is None:
            # Unknown comb-rule: compare numerics but warn about uncertainty
            print(
                f"  [WARN] Comb-rule unknown for atomtype '{self.name}' during comparison; "
                "numeric parameters will be compared without comb-rule certainty."
            )
        
        # Use relative tolerance with absolute floor depending on representation
        # If comb-rule is unknown on either side, prefer the known side for tolerances
        if self.comb_rule is None and other.comb_rule is not None:
            lj_kind = other.lj_kind
        elif other.comb_rule is None and self.comb_rule is not None:
            lj_kind = self.lj_kind
        else:
            lj_kind = self.lj_kind

        if lj_kind == LJParamKind.C6_C12:
            rel_tol = LJ_REL_TOL_C6_C12
            floor1 = LJ_ABS_FLOOR_C6
            floor2 = LJ_ABS_FLOOR_C12
        else:
            rel_tol = LJ_REL_TOL_SIGMA_EPSILON
            floor1 = LJ_ABS_FLOOR_SIGMA
            floor2 = LJ_ABS_FLOOR_EPSILON
        
        if not self._vals_match(self.sigma, other.sigma, rel_tol, floor1):
            return False
        if not self._vals_match(self.epsilon, other.epsilon, rel_tol, floor2):
            return False
        
        return True
    
    @staticmethod
    def _vals_match(a: float, b: float, rel_tol: float, abs_floor: float) -> bool:
        """Compare values with relative tolerance and absolute floor."""
        diff = abs(a - b)
        max_val = max(abs(a), abs(b))
        threshold = max(abs_floor, max_val * rel_tol)
        return diff < threshold
    
    def charge_matches(self, other: "AtomtypeEntry", tolerance: float = CHARGE_TOLERANCE) -> bool:
        """Check if charges match within tolerance."""
        return abs(self.charge - other.charge) < tolerance
    
    def priority(self) -> int:
        """Get priority for override selection."""
        return FILE_PRIORITY.get(self.file_class, 0)
    
    def validate_lj_bounds(self) -> List[str]:
        """
        Validate LJ parameters against sanity bounds.
        
        Returns list of warning messages for outliers.
        
        Task C: Skip warnings for dummy/virtual atoms (ptype D/V) with LJ=(0,0),
        as this is normal for TIP4P-like models.
        """
        warnings = []
        
        # Task C: Virtual sites / dummies with LJ=0 are normal
        if self.ptype in DUMMY_PTYPES and self.sigma == 0 and self.epsilon == 0:
            return []  # No warnings for expected zero params
        
        if self.lj_kind == LJParamKind.SIGMA_EPSILON:
            if not (SIGMA_BOUNDS_NM[0] <= self.sigma <= SIGMA_BOUNDS_NM[1]):
                warnings.append(
                    f"Atomtype '{self.name}': sigma={self.sigma:.4f} nm outside bounds "
                    f"{SIGMA_BOUNDS_NM} (typical: 0.1-0.5 nm)"
                )
            if not (EPSILON_BOUNDS_KJ_MOL[0] <= self.epsilon <= EPSILON_BOUNDS_KJ_MOL[1]):
                warnings.append(
                    f"Atomtype '{self.name}': epsilon={self.epsilon:.4f} kJ/mol outside bounds "
                    f"{EPSILON_BOUNDS_KJ_MOL}"
                )
            # Only warn about sigma<=0 if NOT a dummy/virtual type
            if self.sigma <= 0 and self.ptype not in DUMMY_PTYPES:
                warnings.append(f"Atomtype '{self.name}': sigma={self.sigma} must be > 0")
            if self.epsilon < 0:
                warnings.append(f"Atomtype '{self.name}': epsilon={self.epsilon} must be >= 0")
        else:  # C6/C12
            c6, c12 = self.sigma, self.epsilon
            if not (C6_BOUNDS[0] <= c6 <= C6_BOUNDS[1]):
                warnings.append(
                    f"Atomtype '{self.name}': C6={c6:.4e} outside bounds {C6_BOUNDS}"
                )
            if not (C12_BOUNDS[0] <= c12 <= C12_BOUNDS[1]):
                warnings.append(
                    f"Atomtype '{self.name}': C12={c12:.4e} outside bounds {C12_BOUNDS}"
                )
            if c6 < 0:
                warnings.append(f"Atomtype '{self.name}': C6={c6} must be >= 0")
            if c12 < 0:
                warnings.append(f"Atomtype '{self.name}': C12={c12} must be >= 0")
        return warnings


@dataclass
class ChargeWarning:
    """Records a charge mismatch warning."""
    atomtype_name: str
    charges: List[float]
    source_files: List[str]
    
    def __str__(self) -> str:
        charges_str = ", ".join(f"{c:.4f}" for c in self.charges)
        files_str = ", ".join(Path(f).name for f in self.source_files)
        return f"Charge mismatch for '{self.atomtype_name}': [{charges_str}] from [{files_str}]"


@dataclass
class Conflict:
    """Represents a detected atomtype conflict."""
    atomtype_name: str
    entries: List[AtomtypeEntry]
    conflict_type: str  # "parameter_mismatch" or "cross_ff_collision"
    winner: Optional[AtomtypeEntry] = None  # Set when overridden
    
    def __str__(self) -> str:
        lines = [f"Conflict [{self.conflict_type}] for atomtype '{self.atomtype_name}':"]
        for e in self.entries:
            files = [Path(f).name for f in e.source_files[:3]]
            lines.append(
                f"  - {e.source_ff}: mass={e.mass:.4f}, σ={e.sigma:.6e}, ε={e.epsilon:.6e} (files: {files})"
            )
        return "\n".join(lines)


@dataclass
class IncludeResolution:
    """Records how a #include target was resolved (or not)."""
    including_file: str
    include_target: str
    resolved_path: Optional[str]
    resolved_from_dir: Optional[str]
    attempted_paths: List[str]
    shadowed_candidates: List[str] = field(default_factory=list)


@dataclass
class BondedTypeEntry:
    """A parsed bonded type-table row (bondtypes/angletypes/dihedraltypes)."""
    section: str
    atomtypes: Tuple[str, ...]
    funct: str
    params: Tuple[str, ...]
    source_file: str
    wildcard: bool = False


@dataclass 
class SanitizerResult:
    """Result of a sanitization run."""
    combined_atomtypes_path: Path
    combined_atomtypes_current_path: Path
    sanitized_dir: Path
    sanitized_current_dir: Path
    sanitized_files: List[Path]
    molecule_itp_files: List[Path]  # Files with [ moleculetype ] section
    conflicts_detected: List[Conflict]
    conflicts_overridden: List[Conflict]
    charge_warnings: List[ChargeWarning]
    atomtype_count: int
    source_files_processed: int
    defaults_used: Optional[DefaultsEntry] = None
    source_defaults_map: Dict[Path, DefaultsEntry] = field(default_factory=dict)
    mixed_defaults_detected: bool = False
    mixed_defaults_report: Dict[str, Any] = field(default_factory=dict)
    include_resolution: List[IncludeResolution] = field(default_factory=list)
    defaults_guessed: bool = False
    nonbond_params_secondary_path: Optional[Path] = None
    nonbond_params_secondary_current_path: Optional[Path] = None
    nonbond_params_secondary_summary: Dict[str, Any] = field(default_factory=dict)
    nonbond_params_cross_group_path: Optional[Path] = None
    nonbond_params_cross_group_current_path: Optional[Path] = None
    nonbond_params_cross_group_summary: Dict[str, Any] = field(default_factory=dict)
    prefix_policy_summary: Dict[str, Any] = field(default_factory=dict)
    prefix_injected_types_path: Optional[Path] = None
    prefix_injected_types_current_path: Optional[Path] = None


# =============================================================================
# ITP Parser
# =============================================================================

class ItpParser:
    """
    Parser for GROMACS ITP files.
    
    Handles extraction of sections demarcated by [ section_name ].
    Uses token-based parsing for robustness across ITP variants.
    """
    
    # Task A: Regex for section headers - tolerant to trailing comments
    # Examples that should match:
    #   [ atomtypes ]
    #   [ atomtypes ] ; Define my atoms here
    #   [  defaults  ]   ; trailing whitespace and comment
    SECTION_RE = re.compile(r'^\s*\[\s*(\w+)\s*\](?:\s*;.*)?\s*$')
    
    # Regex for #include directives
    INCLUDE_RE = re.compile(r'^\s*#include\s+["\']([^"\']+)["\']')
    
    @classmethod
    def parse_sections(cls, content: str) -> Dict[str, List[str]]:
        """
        Parse ITP content into sections.
        
        Returns dict mapping section name to list of content lines.
        Lines not in any section go under key "".
        """
        sections: Dict[str, List[str]] = {"": []}
        current_section = ""
        
        for line in content.splitlines():
            match = cls.SECTION_RE.match(line)
            if match:
                current_section = match.group(1).lower().strip()
                if current_section not in sections:
                    sections[current_section] = []
                sections[current_section].append(line)
            else:
                if current_section not in sections:
                    sections[current_section] = []
                sections[current_section].append(line)
        
        return sections
    
    @classmethod
    def parse_defaults(cls, lines: List[str], source_file: str = "") -> Optional[DefaultsEntry]:
        """
        Parse [ defaults ] section.
        
        Args:
            lines: Lines from defaults section (including header)
            source_file: Source file path for tracking
            
        Returns:
            DefaultsEntry or None if parsing fails
        """
        for line in lines:
            # Strip inline comments
            stripped = line.split(";")[0].strip()
            if not stripped or stripped.startswith("["):
                continue
            
            tokens = stripped.split()
            if len(tokens) >= 5:
                try:
                    return DefaultsEntry(
                        nbfunc=int(tokens[0]),
                        comb_rule=int(tokens[1]),
                        gen_pairs=tokens[2].lower(),
                        fudge_lj=float(tokens[3]),
                        fudge_qq=float(tokens[4]),
                        source_file=source_file,
                    )
                except (ValueError, IndexError):
                    continue
        return None
    
    @classmethod
    def parse_atomtypes(
        cls,
        lines: List[str],
        source_ff: str = "",
        source_file: str = "",
        file_class: str = "unknown",
        comb_rule: Optional[int] = None,
    ) -> List[AtomtypeEntry]:
        """
        Token-based parser for [ atomtypes ] section.
        
        Uses TAIL-STRUCTURE parsing for robustness:
        - Last 2 tokens: always LJ params (sigma/epsilon or C6/C12)
        - Token before that: ptype (validated against VALID_PTYPES)
        - Preceding tokens: charge(float), mass(float)
        - Remaining middle tokens: optional at_num, bond_type
        
        Args:
            lines: Lines from the atomtypes section (including header)
            source_ff: Forcefield identifier
            source_file: Source file path for tracking
            file_class: File classification for priority
            comb_rule: Comb-rule from [ defaults ] (1=C6/C12, 2/3=sigma/epsilon)
            
        Returns:
            List of AtomtypeEntry objects
            
        Raises:
            ItpSanitizerError: If parsing fails with actionable error
        """
        entries = []
        
        for line_num, raw_line in enumerate(lines, 1):
            # Strip inline comments
            line = raw_line.split(";")[0].strip()
            
            # Skip empty lines and section header
            if not line or line.startswith("["):
                continue
            
            tokens = line.split()
            n = len(tokens)
            
            # Need at least 6 tokens for simplest format
            if n < 6:
                continue
            
            try:
                entry = cls._parse_atomtype_tokens(
                    tokens, source_ff, source_file, file_class, raw_line, comb_rule
                )
                if entry:
                    entries.append(entry)
            except ValueError as e:
                raise ItpSanitizerError(
                    f"Failed to parse atomtype in {source_file}:\n"
                    f"  Line {line_num}: {raw_line.strip()}\n"
                    f"  Error: {e}"
                )
        
        return entries
    
    @classmethod
    def _parse_atomtype_tokens(
        cls,
        tokens: List[str],
        source_ff: str,
        source_file: str,
        file_class: str,
        raw_line: str,
        comb_rule: Optional[int] = None,
    ) -> Optional[AtomtypeEntry]:
        """
        Parse tokens into AtomtypeEntry using robust TAIL-STRUCTURE approach.
        
        Task D: Tolerant parser that handles missing ptype.
        
        Strategy:
        1. Last 2 tokens are always LJ params (sigma/epsilon or C6/C12)
        2. Search remaining tail tokens for single-letter ptype (A/D/S/V)
        3. If ptype not found, allow ptype=None and continue
        4. Parse charge, mass from remaining tokens
        
        Tail structure variants:
        - tokens[-1]: epsilon (or C12)
        - tokens[-2]: sigma (or C6)
        - tokens[-3]: ptype (if present) OR charge (if ptype missing)
        - tokens[-4]: charge (if ptype present) OR mass
        - tokens[-5]: mass (if ptype present)
        - tokens[0]: name
        - middle tokens: at_num and/or bond_type
        """
        n = len(tokens)
        if n < 5:  # Minimum: name mass charge sigma epsilon (no ptype)
            return None
        
        # Parse LJ params from tail (always last 2)
        try:
            lj_param2 = float(tokens[-1])  # epsilon or C12
            lj_param1 = float(tokens[-2])  # sigma or C6
        except (ValueError, IndexError) as e:
            print(
                f"  [WARN] Cannot parse LJ parameters from tail tokens: {raw_line.strip()}"
            )
            raise ValueError(f"Cannot parse LJ parameters from tokens[-1], tokens[-2]: {e}")
        
        # Task D: Try to find ptype - check if tokens[-3] is a valid ptype
        # NOTE: In rare edge-cases, an atomtype name could be "A/S/D/V".
        # Under standard GROMACS atomtypes formats, this is still safe because
        # the name is token[0] and ptype is only accepted from token[-3].
        ptype = None
        ptype_idx = None
        
        # Check if tokens[-3] is a valid ptype (single letter A/D/S/V)
        candidate = tokens[-3].upper() if n >= 3 else ""
        if candidate in VALID_PTYPES and n >= 6:
            ptype = candidate
            ptype_idx = n - 3
        else:
            # ptype might be missing - this is allowed per Task D
            # In this case: name [optional...] mass charge sigma epsilon (5 tokens min)
            ptype = None
            ptype_idx = None
        
        # Parse charge and mass based on whether ptype was found
        try:
            if ptype is not None:
                # Standard format: ... mass charge ptype sigma epsilon
                charge = float(tokens[-4])
                mass = float(tokens[-5])
                tail_consumed = 5  # mass, charge, ptype, sigma, epsilon
            else:
                # Missing ptype: ... mass charge sigma epsilon
                charge = float(tokens[-3])
                mass = float(tokens[-4])
                tail_consumed = 4  # mass, charge, sigma, epsilon
        except (ValueError, IndexError):
            # Guard: if ptype guess leads to parse failure, fall back to missing ptype
            if ptype is not None:
                try:
                    charge = float(tokens[-3])
                    mass = float(tokens[-4])
                    tail_consumed = 4
                    ptype = None
                    print(
                        f"  [WARN] Ptype inference failed; treating as missing: {raw_line.strip()}"
                    )
                except (ValueError, IndexError) as e:
                    raise ValueError(f"Cannot parse mass/charge: {e}")
            else:
                raise ValueError("Cannot parse mass/charge from tokens")
        
        # Name is always first
        name = tokens[0]
        
        # Determine LJ parameter kind based on comb-rule
        if comb_rule == 1:
            lj_kind = LJParamKind.C6_C12
        else:
            lj_kind = LJParamKind.SIGMA_EPSILON
        
        # Middle tokens (between name and tail-consumed tokens)
        # These can be at_num, bond_type, or both
        middle_tokens = tokens[1:n-tail_consumed]
        
        at_num = None
        bond_type = None
        
        if len(middle_tokens) == 0:
            pass
        elif len(middle_tokens) == 1:
            # Could be at_num (int) or bond_type (string)
            try:
                at_num = int(middle_tokens[0])
            except ValueError:
                bond_type = middle_tokens[0]
        elif len(middle_tokens) >= 2:
            # First is bond_type, second is at_num
            bond_type = middle_tokens[0]
            try:
                at_num = int(middle_tokens[1])
            except ValueError:
                # Not critical - could be unusual format
                pass
        
        return AtomtypeEntry(
            name=name,
            mass=mass,
            charge=charge,
            ptype=ptype,  # May be None if missing
            sigma=lj_param1,
            epsilon=lj_param2,
            at_num=at_num,
            bond_type=bond_type,
            source_ff=source_ff,
            source_files=[source_file] if source_file else [],
            file_class=file_class,
            raw_line=raw_line,
            lj_kind=lj_kind,
            comb_rule=comb_rule,
        )
    
    @classmethod
    def extract_atomtypes_from_file(
        cls,
        itp_path: Path,
        source_ff: str = "",
        file_class: str = "unknown",
        comb_rule: Optional[int] = None,
    ) -> List[AtomtypeEntry]:
        """
        Extract atomtype entries from an ITP file.
        
        Args:
            itp_path: Path to ITP file
            source_ff: Forcefield identifier
            file_class: File classification
            comb_rule: Comb-rule from [ defaults ] (1=C6/C12, 2/3=sigma/epsilon)
            
        Returns:
            List of AtomtypeEntry objects found in the file
        """
        if not itp_path.exists():
            return []
        
        content = itp_path.read_text()
        sections = cls.parse_sections(content)
        
        if 'atomtypes' not in sections:
            return []
        
        return cls.parse_atomtypes(
            sections['atomtypes'],
            source_ff=source_ff,
            source_file=str(itp_path),
            file_class=file_class,
            comb_rule=comb_rule,
        )
    
    @classmethod
    def extract_defaults_from_file(cls, itp_path: Path) -> Optional[DefaultsEntry]:
        """Extract [ defaults ] section from an ITP file."""
        if not itp_path.exists():
            return None
        
        content = itp_path.read_text()
        sections = cls.parse_sections(content)
        
        if 'defaults' not in sections:
            return None
        
        return cls.parse_defaults(sections['defaults'], str(itp_path))

    @classmethod
    def has_defaults_or_atomtypes(cls, itp_path: Path) -> bool:
        """Check if a file contains [ defaults ] or [ atomtypes ] sections."""
        if not itp_path.exists():
            return False
        content = itp_path.read_text()
        sections = cls.parse_sections(content)
        return ('defaults' in sections) or ('atomtypes' in sections)
    
    @classmethod
    def has_moleculetype_section(cls, itp_path: Path) -> bool:
        """Check if file contains [ moleculetype ] section."""
        if not itp_path.exists():
            return False
        
        content = itp_path.read_text()
        sections = cls.parse_sections(content)
        return 'moleculetype' in sections

    @classmethod
    def has_moleculetype_including_includes(
        cls,
        itp_path: Path,
        include_dirs: Optional[List[Path]] = None,
        strict: bool = False,
        check_shadowing: bool = False,
        allow_shadowing: bool = False,
    ) -> bool:
        """
        Check whether a file or any file in its include-closure contains [ moleculetype ].
        
        This is include-aware and follows GROMACS -I semantics using resolve_includes().
        """
        if not itp_path.exists():
            return False
        
        if cls.has_moleculetype_section(itp_path):
            return True
        
        includes = cls.resolve_includes(
            itp_path,
            include_dirs=include_dirs,
            strict=strict,
            diagnostics=None,
            check_shadowing=check_shadowing,
            allow_shadowing=allow_shadowing,
        )
        for inc in includes:
            if cls.has_moleculetype_section(inc):
                return True
        
        return False
    
    @classmethod
    def resolve_includes(
        cls,
        itp_path: Path,
        visited: Optional[Set[Path]] = None,
        include_dirs: Optional[List[Path]] = None,
        strict: bool = False,
        diagnostics: Optional[List["IncludeResolution"]] = None,
        check_shadowing: bool = False,
        allow_shadowing: bool = False,
    ) -> List[Path]:
        """
        Recursively resolve #include directives.
        
        Task I: Ensures deduplication and shadowing detection.
        
        Args:
            itp_path: Starting ITP file
            visited: Set of already-visited paths
            include_dirs: Additional directories to search
            strict: If True, fail-fast on missing includes
            diagnostics: Optional list to collect resolution details
            check_shadowing: If True, check for shadowed files (Issue 1)
            allow_shadowing: If True, allow shadowing without error even in strict mode
            
        Returns:
            List of all included files (depth-first order, deduplicated)
            
        Raises:
            ItpSanitizerError: If strict and include not found OR (strict and shadowing detected)
        """
        if visited is None:
            visited = set()
        if include_dirs is None:
            include_dirs = []
        
        itp_path = itp_path.resolve()
        if itp_path in visited:
            return []
        visited.add(itp_path)
        
        if not itp_path.exists():
            return []
        
        includes = []
        try:
            content = itp_path.read_text()
        except Exception:
            return []
        
        for line in content.splitlines():
            match = cls.INCLUDE_RE.match(line)
            if match:
                include_rel = match.group(1)

                # Search order: relative to including file, then include_dirs
                search_paths = [itp_path.parent] + list(include_dirs)
                search_paths = [p for p in search_paths if p.exists()]
                
                include_path = None
                resolved_from = None
                attempted_paths: List[str] = []
                found_paths: List[Path] = []
                
                # Check absolute path first (bypass search)
                include_candidate = Path(include_rel)
                if include_candidate.is_absolute():
                    candidate = include_candidate.resolve()
                    attempted_paths.append(str(candidate))
                    if candidate.exists():
                        found_paths.append(candidate)
                        include_path = candidate
                        resolved_from = str(candidate.parent)
                else:
                    # Scan all search paths
                    for search_dir in search_paths:
                        candidate = (search_dir / include_rel).resolve()
                        attempted_paths.append(str(candidate))
                        
                        if candidate.exists():
                            if not found_paths:
                                # First match wins
                                include_path = candidate
                                resolved_from = str(search_dir)
                            if check_shadowing:
                                found_paths.append(candidate)
                            else:
                                # Optimization: stop at first match if not checking shadowing
                                found_paths.append(candidate)
                                break

                # Handle Shadowing Detection (Issue 1)
                shadowed_by = []
                if check_shadowing and len(found_paths) > 1:
                    # Filter out same file reachable via different paths (symlinks or .. normalization)
                    unique_found = []
                    seen_resolved = set()
                    for p in found_paths:
                        r = p.resolve()
                        if r not in seen_resolved:
                            seen_resolved.add(r)
                            unique_found.append(p)
                    
                    if len(unique_found) > 1:
                        top_winner = unique_found[0]
                        shadowed = unique_found[1:]
                        shadowed_str = ", ".join(str(p) for p in shadowed)
                        msg = (
                            f"Include shadowing detected for '{include_rel}' in {itp_path.name}:\n"
                            f"  Used: {top_winner}\n"
                            f"  Shadowed: {shadowed_str}"
                        )
                        print(f"  [WARN] {msg}")
                        shadowed_by = [str(p) for p in shadowed]
                        
                        if strict and not allow_shadowing:
                            raise ItpSanitizerError(
                                f"Include shadowing violation (strict mode): {msg}\n"
                                "Use --allow-include-shadowing to bypass."
                            )

                if include_path is not None:
                    if diagnostics is not None:
                        diagnostics.append(IncludeResolution(
                            including_file=str(itp_path),
                            include_target=include_rel,
                            resolved_path=str(include_path),
                            resolved_from_dir=resolved_from,
                            attempted_paths=attempted_paths,
                            shadowed_candidates=shadowed_by,
                        ))
                    
                    if include_path not in visited:
                        includes.append(include_path)
                        # Recurse with same params
                        includes.extend(cls.resolve_includes(
                            include_path, 
                            visited, 
                            include_dirs, 
                            strict, 
                            diagnostics,
                            check_shadowing,
                            allow_shadowing
                        ))
                else:
                    if diagnostics is not None:
                        diagnostics.append(IncludeResolution(
                            including_file=str(itp_path),
                            include_target=include_rel,
                            resolved_path=None,
                            resolved_from_dir=None,
                            attempted_paths=attempted_paths,
                        ))
                    if strict:
                        raise ItpSanitizerError(
                            f"Missing include: '{include_rel}' from {itp_path.name}\n"
                            f"Attempted paths:\n"
                            + "\n".join(f"  - {p}" for p in attempted_paths) + "\n"
                            "Use --itp-include-dirs to add search paths, or disable strict mode."
                        )
                    print(
                        f"  [WARN] Include not found: {include_rel} (from {itp_path.name})\n"
                        f"         Attempted paths: {attempted_paths}"
                    )
    
        return includes
    
    @classmethod
    def detect_implicit_topology(cls, content: str) -> List[str]:
        """
        Detect sections that use implicit/library parameter lookup.
        
        These sections reference bondtypes/angletypes/dihedraltypes from
        forcefield library files. If we prefix atomtypes, these lookups
        will break unless we also update the library files.
        
        Returns:
            List of section names with implicit parameter usage
        """
        implicit_sections = []
        sections = cls.parse_sections(content)
        
        # Check bonded sections for implicit parameters
        # Format: ai aj [ak [al]] funct [params...]
        # If params are missing, it's an implicit lookup
        
        checks = [
            ("bonds", 2, 3),      # ai aj funct [params] -> need >=4 for explicit
            ("angles", 3, 4),    # ai aj ak funct [params] -> need >=5 for explicit
            ("dihedrals", 4, 5), # ai aj ak al funct [params] -> need >=6 for explicit
            ("impropers", 4, 5), # same as dihedrals
        ]
        
        for section_name, atom_count, min_explicit in checks:
            if section_name not in sections:
                continue
            
            for line in sections[section_name]:
                stripped = line.split(";")[0].strip()
                if not stripped or stripped.startswith("["):
                    continue
                
                tokens = stripped.split()
                # Has atom indices + funct but no explicit parameters
                if len(tokens) == atom_count + 1:
                    implicit_sections.append(section_name)
                    break
        
        return list(set(implicit_sections))
    
    @classmethod
    def strip_atomtypes_section(cls, content: str) -> str:
        """
        Remove the [ atomtypes ] section from ITP content.
        
        Preserves all other sections and blank lines.
        """
        lines = content.splitlines()
        result = []
        in_atomtypes = False
        
        for line in lines:
            match = cls.SECTION_RE.match(line)
            if match:
                section_name = match.group(1).lower()
                if section_name == 'atomtypes':
                    in_atomtypes = True
                    continue
                else:
                    in_atomtypes = False
            
            if not in_atomtypes:
                result.append(line)
        
        return '\n'.join(result)
    
    @classmethod
    def strip_defaults_section(cls, content: str) -> str:
        """Remove the [ defaults ] section from ITP content."""
        lines = content.splitlines()
        result = []
        in_defaults = False
        
        for line in lines:
            match = cls.SECTION_RE.match(line)
            if match:
                section_name = match.group(1).lower()
                if section_name == 'defaults':
                    in_defaults = True
                    continue
                else:
                    in_defaults = False
            
            if not in_defaults:
                result.append(line)
        
        return '\n'.join(result)
    
    @classmethod
    def apply_atomtype_prefix(
        cls,
        content: str,
        prefix: str,
        atomtype_names: Set[str],
    ) -> Tuple[str, int]:
        """
        Apply a prefix to atomtype names in ITP content.
        
        Updates both [ atomtypes ] and [ atoms ] sections.
        
        Args:
            content: ITP file content
            prefix: Prefix to add (e.g., "G2_")
            atomtype_names: Set of atomtype names to prefix
            
        Returns:
            Tuple of (modified content, count of replacements)
        """
        lines = content.splitlines()
        result = []
        replacement_count = 0
        current_section = ""

        def _split_comment(raw_line: str) -> Tuple[str, str]:
            idx = raw_line.find(";")
            if idx < 0:
                return raw_line, ""
            return raw_line[:idx], raw_line[idx:]

        def _rewrite_type_tokens(
            raw_line: str,
            token_count: int,
        ) -> Tuple[str, int]:
            pre_comment, comment = _split_comment(raw_line)
            leading_ws_match = re.match(r"^(\s*)", pre_comment)
            leading_ws = leading_ws_match.group(1) if leading_ws_match else ""
            payload = pre_comment.strip()
            if not payload:
                return raw_line, 0
            tokens = payload.split()
            if len(tokens) < token_count:
                return raw_line, 0
            local_replacements = 0
            for idx in range(token_count):
                token = tokens[idx]
                if token in atomtype_names:
                    tokens[idx] = f"{prefix}{token}"
                    local_replacements += 1
            if local_replacements == 0:
                return raw_line, 0
            rewritten = leading_ws + " ".join(tokens)
            if comment:
                if rewritten and not rewritten.endswith(" "):
                    rewritten += " "
                rewritten += comment.lstrip()
            return rewritten, local_replacements
        
        for line in lines:
            match = cls.SECTION_RE.match(line)
            if match:
                current_section = match.group(1).lower()
                result.append(line)
                continue
            
            # Skip comments and empty lines
            stripped = line.strip()
            if not stripped or stripped.startswith(';'):
                result.append(line)
                continue
            
            # Process atomtypes section
            if current_section == 'atomtypes':
                for name in atomtype_names:
                    if re.match(rf'^\s*{re.escape(name)}\s+', line):
                        line = re.sub(
                            rf'^(\s*){re.escape(name)}(\s+)',
                            rf'\g<1>{prefix}{name}\g<2>',
                            line,
                        )
                        replacement_count += 1
                        break
            
            # Process atoms section (type column is column 2, 1-indexed)
            elif current_section == 'atoms':
                parts = stripped.split()
                if len(parts) >= 2:
                    atomtype = parts[1]
                    if atomtype in atomtype_names:
                        # Replace the type field with prefixed version
                        pattern = rf'^(\s*\S+\s+){re.escape(atomtype)}(\s+)'
                        replacement = rf'\g<1>{prefix}{atomtype}\g<2>'
                        new_line = re.sub(pattern, replacement, line)
                        if new_line != line:
                            line = new_line
                            replacement_count += 1
            elif current_section in TYPE_SECTION_ATOMTYPE_TOKEN_COUNT:
                rewritten, local_count = _rewrite_type_tokens(
                    line, TYPE_SECTION_ATOMTYPE_TOKEN_COUNT[current_section]
                )
                if local_count > 0:
                    line = rewritten
                    replacement_count += local_count
            
            result.append(line)
        
        return '\n'.join(result), replacement_count


# =============================================================================
# ITP Sanitizer - Main Class
# =============================================================================

class ItpSanitizer:
    """
    Main sanitizer class for processing ITP files.
    
    Implements the sanitization workflow:
    1. Scan all input ITP files (in deterministic order)
    2. Extract and validate [ defaults ] sections
    3. Extract and register atomtypes
    4. Detect conflicts (same name, different params)
    5. Detect namespace collisions (cross-FF contamination)
    6. Handle charge mismatches as warnings
    7. Generate combined atomtypes file (with deterministic charge)
    8. Generate sanitized ITP files (atomtypes removed)
    9. Generate system.top with correct include order
    """
    
    def __init__(
        self,
        source_ff: str = "",
        allow_override: bool = False,
        atomtype_prefix: Optional[str] = None,
        strict_charge: bool = False,
    ):
        """
        Initialize sanitizer.
        
        Args:
            source_ff: Default forcefield identifier
            allow_override: Allow conflicts to be resolved by priority-based override
            atomtype_prefix: Optional prefix for atomtype names
            strict_charge: Fail on charge mismatches (default: warn only)
        """
        self.source_ff = source_ff
        self.allow_override = allow_override
        self.atomtype_prefix = atomtype_prefix
        self.strict_charge = strict_charge
        self.effective_atomtype_prefix: Optional[str] = atomtype_prefix
        
        # Registry: atomtype_name -> AtomtypeEntry
        self.registry: Dict[str, AtomtypeEntry] = {}
        
        # Track all source files
        self.source_files: List[Path] = []
        
        # Defaults from all files
        self.defaults_entries: List[DefaultsEntry] = []
        self.canonical_defaults: Optional[DefaultsEntry] = None
        self.defaults_report: Dict[str, Any] = {}
        
        # Conflicts and warnings
        self.conflicts: List[Conflict] = []
        self.overridden_conflicts: List[Conflict] = []
        self.charge_warnings: List[ChargeWarning] = []
        
        # Charge tracking for deterministic output
        self._charge_by_atomtype: Dict[str, List[Tuple[float, str]]] = {}  # name -> [(charge, file)]
        # Optional source-content overrides (used by prefix bake policy).
        self._content_overrides: Dict[Path, str] = {}
        self.prefix_policy_summary: Dict[str, Any] = {}
        self.prefix_injected_types_path: Optional[Path] = None
        self.prefix_injected_types_current_path: Optional[Path] = None
    
    def scan_itp_files(
        self,
        itp_paths: List[Path],
        ff_map: Optional[Dict[str, str]] = None,
        class_map: Optional[Dict[str, str]] = None,
        comb_rule_override: Optional[int] = None,
        collect_defaults: bool = True,
    ) -> Dict[str, AtomtypeEntry]:
        """
        Scan multiple ITP files and build atomtype registry.
        
        IMPORTANT: [ defaults ] are collected separately when needed.
        If comb_rule_override is provided, atomtypes are parsed using it.
        
        Files are processed in deterministic order (sorted by path with casefold).
        
        Args:
            itp_paths: List of ITP file paths
            ff_map: Mapping of file path string to forcefield ID
            class_map: Mapping of file path string to file class
            comb_rule_override: Optional comb-rule to use for atomtype parsing
            collect_defaults: If True, collect [ defaults ] entries while scanning
            
        Returns:
            Registry of atomtype name to entry
        """
        ff_map = ff_map or {}
        class_map = class_map or {}
        
        # Sort for reproducibility (casefold for cross-platform stability)
        sorted_paths = sorted(
            itp_paths,
            key=lambda p: (p.as_posix().casefold(), p.as_posix())
        )
        
        # === PHASE 1: Collect ALL defaults (optional) ===
        if collect_defaults:
            for itp_path in sorted_paths:
                if not itp_path.exists():
                    continue
                defaults = ItpParser.extract_defaults_from_file(itp_path)
                if defaults:
                    self.defaults_entries.append(defaults)
        
        # Get comb-rule for atomtype parsing
        comb_rule = comb_rule_override
        if comb_rule is None:
            comb_rule = self.canonical_defaults.comb_rule if self.canonical_defaults else None
        
        # === PHASE 2: Parse atomtypes with known comb-rule ===
        for itp_path in sorted_paths:
            if not itp_path.exists():
                print(f"  [WARN] ITP file not found: {itp_path}")
                continue
            
            self.source_files.append(itp_path)
            
            # Determine source FF and class
            file_ff = ff_map.get(str(itp_path), "UNKNOWN")
            file_class = class_map.get(str(itp_path), "unknown")
            
            # Extract atomtypes with comb-rule awareness
            entries = ItpParser.extract_atomtypes_from_file(
                itp_path, 
                source_ff=file_ff,
                file_class=file_class,
                comb_rule=comb_rule,
            )
            
            # Validate LJ bounds for each entry
            for entry in entries:
                for warning in entry.validate_lj_bounds():
                    print(f"  [WARN] {warning}")
            
            for entry in entries:
                self._register_atomtype(entry)
        
        return self.registry
    
    def _register_atomtype(self, entry: AtomtypeEntry) -> None:
        """
        Register an atomtype entry, detecting conflicts.
        
        Conflict detection is based on physical parameters (excluding charge).
        Charge mismatches are tracked separately as warnings.
        """
        name = entry.name
        
        # Track charge for later deterministic output
        if name not in self._charge_by_atomtype:
            self._charge_by_atomtype[name] = []
        self._charge_by_atomtype[name].append((
            entry.charge, 
            entry.source_files[0] if entry.source_files else ""
        ))
        
        if name not in self.registry:
            self.registry[name] = entry
            return
        
        existing = self.registry[name]
        
        # Check for physical parameter mismatch (signature excludes charge)
        if not existing.params_match(entry):
            # Determine conflict type
            if (existing.source_ff and entry.source_ff and 
                existing.source_ff != entry.source_ff and
                existing.source_ff != "UNKNOWN" and entry.source_ff != "UNKNOWN"):
                conflict_type = "cross_ff_collision"
            else:
                conflict_type = "parameter_mismatch"
            
            conflict = Conflict(
                atomtype_name=name,
                entries=[existing, entry],
                conflict_type=conflict_type,
            )
            
            self.conflicts.append(conflict)
            
            if self.allow_override:
                # Task G: Prefer entries matching active forcefield over priority alone
                existing_matches_ff = (
                    self.source_ff and existing.source_ff == self.source_ff
                )
                entry_matches_ff = (
                    self.source_ff and entry.source_ff == self.source_ff
                )
                
                if entry_matches_ff and not existing_matches_ff:
                    # New entry matches active FF, existing doesn't -> prefer new
                    winner = entry
                elif existing_matches_ff and not entry_matches_ff:
                    # Existing matches active FF, new doesn't -> prefer existing
                    winner = existing
                elif entry.priority() > existing.priority():
                    # Neither or both match FF -> fall back to priority
                    winner = entry
                else:
                    winner = existing  # Existing wins on tie
                
                conflict.winner = winner
                self.overridden_conflicts.append(conflict)
                self.registry[name] = winner
                
                winner_file = winner.source_files[0] if winner.source_files else "unknown"
                loser = entry if winner == existing else existing
                loser_file = loser.source_files[0] if loser.source_files else "unknown"
                print(
                    f"  [WARN] Atomtype conflict overridden: {name}\n"
                    f"    Winner: {winner.source_ff} from {Path(winner_file).name}\n"
                    f"    Loser:  {loser.source_ff} from {Path(loser_file).name}"
                )
        else:
            # Physical params match - merge source_files list
            existing.source_files.extend(entry.source_files)
            
            # Check charge mismatch (warning, not error)
            if not existing.charge_matches(entry):
                # Record warning
                warning = ChargeWarning(
                    atomtype_name=name,
                    charges=[existing.charge, entry.charge],
                    source_files=existing.source_files[:1] + entry.source_files[:1],
                )
                # Avoid duplicate warnings
                if not any(w.atomtype_name == name for w in self.charge_warnings):
                    self.charge_warnings.append(warning)
                    print(f"  [WARN] {warning}")
    
    def validate_defaults(
        self,
        allow_mixed_defaults: bool = False,
        allow_mixed_defaults_reason: Optional[str] = None,
        class_map: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """
        Validate [ defaults ] consistency across all files.

        Returns a structured report used for diagnostics/manifest output.
        """
        empty_report: Dict[str, Any] = {
            "mixed_defaults_detected": False,
            "classification": "none",
            "differing_fields": [],
            "comb_rules_present": [],
            "groups": [],
            "primary_defaults_signature": None,
            "secondary_signatures": [],
            "chosen_policy": "none",
            "allow_mixed_defaults": bool(allow_mixed_defaults),
            "allow_mixed_defaults_reason": None,
        }
        if not self.defaults_entries:
            self.defaults_report = dict(empty_report)
            return dict(empty_report)

        # Keep include-encounter order from collection phase for deterministic policy.
        entries_in_order = list(self.defaults_entries)
        order_index = {id(entry): idx for idx, entry in enumerate(entries_in_order)}

        groups: List[Dict[str, Any]] = []
        for entry in entries_in_order:
            placed = False
            for group in groups:
                rep: DefaultsEntry = group["rep"]
                if rep.matches(entry):
                    group["entries"].append(entry)
                    placed = True
                    break
            if not placed:
                groups.append({"rep": entry, "entries": [entry]})

        reps = [group["rep"] for group in groups]
        mixed_defaults = len(groups) > 1

        differing_fields: List[str] = []
        if len({rep.nbfunc for rep in reps}) > 1:
            differing_fields.append("nbfunc")
        if len({rep.comb_rule for rep in reps}) > 1:
            differing_fields.append("comb-rule")
        if len({rep.gen_pairs.lower() for rep in reps}) > 1:
            differing_fields.append("gen-pairs")

        ref_lj = reps[0].fudge_lj if reps else 0.0
        ref_qq = reps[0].fudge_qq if reps else 0.0
        if reps and not all(math.isclose(rep.fudge_lj, ref_lj, rel_tol=1e-4, abs_tol=1e-6) for rep in reps):
            differing_fields.append("fudgeLJ")
        if reps and not all(math.isclose(rep.fudge_qq, ref_qq, rel_tol=1e-4, abs_tol=1e-6) for rep in reps):
            differing_fields.append("fudgeQQ")

        comb_rules_present = {rep.comb_rule for rep in reps}
        has_representation_mismatch = bool(
            (1 in comb_rules_present) and bool(comb_rules_present.intersection({2, 3}))
        )
        if not mixed_defaults:
            classification = "consistent"
        elif has_representation_mismatch:
            classification = "representation-mismatch"
        elif differing_fields == ["comb-rule"]:
            classification = "comb-rule-only"
        elif "comb-rule" in differing_fields:
            classification = "comb-rule-and-other-fields"
        else:
            classification = "non-comb-rule-mismatch"

        # Primary defaults selection policy:
        # 1) forcefield.itp / *.ff defaults, 2) first encountered in include order, 3) previous selection.
        primary_defaults: Optional[DefaultsEntry] = None
        chosen_policy = "reuse_previous_defaults"
        forcefield_entries = [
            entry
            for entry in entries_in_order
            if _is_forcefield_defaults_source(entry.source_file, class_map=class_map)
        ]
        if forcefield_entries:
            primary_defaults = min(
                forcefield_entries,
                key=lambda entry: (
                    order_index.get(id(entry), 10**9),
                    Path(entry.source_file).as_posix().casefold(),
                    Path(entry.source_file).as_posix(),
                ),
            )
            chosen_policy = "prefer_forcefield_defaults"
        elif entries_in_order:
            primary_defaults = entries_in_order[0]
            chosen_policy = "first_encountered_include_order"
        elif self.canonical_defaults is not None:
            primary_defaults = self.canonical_defaults

        if primary_defaults is not None:
            self.canonical_defaults = primary_defaults

        primary_signature = primary_defaults.signature() if primary_defaults else None
        secondary_signatures = [
            group["rep"].signature()
            for group in groups
            if not (primary_defaults and group["rep"].matches(primary_defaults))
        ]

        group_reports: List[Dict[str, Any]] = []
        report_lines: List[str] = []
        for idx, group in enumerate(groups, 1):
            rep: DefaultsEntry = group["rep"]
            files_in_group = [
                entry.source_file
                for entry in sorted(
                    group["entries"],
                    key=lambda entry: (
                        order_index.get(id(entry), 10**9),
                        Path(entry.source_file).as_posix().casefold(),
                        Path(entry.source_file).as_posix(),
                    ),
                )
            ]
            group_reports.append(
                {
                    "signature": rep.signature(),
                    "nbfunc": rep.nbfunc,
                    "comb_rule": rep.comb_rule,
                    "gen_pairs": rep.gen_pairs,
                    "fudge_lj": rep.fudge_lj,
                    "fudge_qq": rep.fudge_qq,
                    "files": files_in_group,
                }
            )
            report_lines.append(
                f"  [{idx}] nbfunc={rep.nbfunc}, comb-rule={rep.comb_rule}, "
                f"gen-pairs={rep.gen_pairs}, fudgeLJ={rep.fudge_lj}, fudgeQQ={rep.fudge_qq}"
            )
            for source_file in files_in_group:
                report_lines.append(f"      - {source_file}")

        reason_clean = (allow_mixed_defaults_reason or "").strip() or None
        report: Dict[str, Any] = {
            "mixed_defaults_detected": mixed_defaults,
            "classification": classification,
            "differing_fields": list(differing_fields),
            "comb_rules_present": sorted(comb_rules_present),
            "groups": group_reports,
            "primary_defaults_signature": primary_signature,
            "secondary_signatures": secondary_signatures,
            "chosen_policy": chosen_policy,
            "allow_mixed_defaults": bool(allow_mixed_defaults),
            "allow_mixed_defaults_reason": reason_clean,
        }

        if mixed_defaults:
            if classification == "comb-rule-only":
                class_line = "Classification: mismatch is comb-rule only."
            elif classification == "representation-mismatch":
                class_line = (
                    "Classification: representation mismatch (comb-rule 1 mixed with 2/3). "
                    "This cannot be safely overridden."
                )
            else:
                class_line = (
                    "Classification: mismatch includes non-comb-rule fields "
                    f"({', '.join(differing_fields) if differing_fields else 'unknown'})."
                )
            msg = (
                "Conflicting [ defaults ] tuples detected across sources:\n"
                + "\n".join(report_lines)
                + f"\n{class_line}\n"
                + f"Primary selection policy: {chosen_policy}\n"
                + "[defaults] is global in GROMACS; mixed tuples can distort LJ mixing."
            )
            if classification == "representation-mismatch":
                remediation = (
                    "\nRemediation options:\n"
                    "  1) Use one representation family only (all comb-rule 1, or all comb-rule 2/3).\n"
                    "  2) Regenerate/convert parameters so atomtypes share a compatible LJ representation.\n"
                    "Unsafe override is intentionally disabled for this mismatch."
                )
                raise DefaultsConflictError(msg + remediation)
            if not allow_mixed_defaults:
                remediation = (
                    "\nRemediation options:\n"
                    "  1) Use a single consistent forcefield/[defaults] tuple.\n"
                    "  2) Use --allow-mixed-defaults with --allow-mixed-defaults-reason and accept risk.\n"
                    "  3) For comb-rule 2/3 mixes, use --mixed-defaults-cross-group-policy generate "
                    "to emit explicit cross-group [ nonbond_params ]."
                )
                raise DefaultsConflictError(msg + remediation)
            if reason_clean is None:
                raise DefaultsConflictError(
                    msg
                    + "\nUnsafe override requires --allow-mixed-defaults-reason for auditability."
                )
            print("  [WARN] ===========================================================")
            print("  [WARN] UNSAFE OVERRIDE: --allow-mixed-defaults is enabled.")
            print(f"  [WARN] Reason: {reason_clean}")
            print(f"  [WARN] {msg}")
            if classification == "comb-rule-only":
                print(
                    "  [WARN] Recommendation: use --mixed-defaults-cross-group-policy generate "
                    "to avoid silent cross-group LJ distortion."
                )
            print("  [WARN] ===========================================================")

        first = self.canonical_defaults
        if first is not None:
            print(
                f"  [ defaults ] primary: {Path(first.source_file).name} "
                f"(nbfunc={first.nbfunc}, comb-rule={first.comb_rule}, policy={chosen_policy})"
            )

        self.defaults_report = dict(report)
        return report
    
    def validate_conflicts(self) -> None:
        """
        Validate registry and raise on unresolved conflicts.
        
        Raises:
            AtomtypeConflictError: If conflicts exist and override not allowed
            NamespaceCollisionError: If cross-FF collision detected
        """
        if self.allow_override:
            return  # All conflicts were resolved
        
        for conflict in self.conflicts:
            if conflict.conflict_type == "cross_ff_collision":
                raise NamespaceCollisionError(
                    f"Cross-forcefield collision detected: {conflict}\n"
                    f"Use --allow-override to force resolution."
                )
            else:
                raise AtomtypeConflictError(
                    f"Atomtype conflict detected: {conflict}\n"
                    f"Use --allow-override to force resolution."
                )
    
    def validate_charge_strict(self) -> None:
        """Fail on charge warnings if strict mode enabled."""
        if self.strict_charge and self.charge_warnings:
            msgs = "\n".join(f"  - {w}" for w in self.charge_warnings)
            raise ItpSanitizerError(
                f"Charge mismatches detected (--strict-charge enabled):\n{msgs}"
            )
    
    def validate_prefix_safety(self, itp_paths: List[Path]) -> None:
        """
        Check that prefixing is safe for all ITP files.
        
        Task H: Fails if any file:
        1. Uses implicit/library-based parameter lookup (bonds without explicit params)
        2. Contains type-based sections (pairtypes, bondtypes, etc.) that reference atomtype names
        
        Raises:
            PrefixUnsafeError: If implicit topology or type-based sections detected
        """
        if not self.atomtype_prefix:
            return
        
        implicit_files = []
        type_based_files = []
        
        for itp_path in itp_paths:
            if not itp_path.exists():
                continue
            content = itp_path.read_text()
            sections = ItpParser.parse_sections(content)
            
            # Check for implicit topology (existing check)
            implicit = ItpParser.detect_implicit_topology(content)
            if implicit:
                implicit_files.append((itp_path, implicit))
            
            # Task H: Check for type-based sections
            found_type_based = [s for s in TYPE_BASED_SECTIONS if s in sections]
            if found_type_based:
                type_based_files.append((itp_path, found_type_based))
        
        # Report errors
        if implicit_files or type_based_files:
            msg_parts = []
            
            if implicit_files:
                msg_parts.append("Files with implicit topology (library-based lookup):")
                for path, sections in implicit_files[:3]:
                    msg_parts.append(f"  - {path.name}: implicit {sections}")
                if len(implicit_files) > 3:
                    msg_parts.append(f"  - ... and {len(implicit_files) - 3} more")
            
            if type_based_files:
                msg_parts.append("Files with type-based sections (reference atomtype names):")
                for path, sections in type_based_files[:3]:
                    msg_parts.append(f"  - {path.name}: {sections}")
                if len(type_based_files) > 3:
                    msg_parts.append(f"  - ... and {len(type_based_files) - 3} more")
            
            raise PrefixUnsafeError(
                f"Cannot apply atomtype prefix safely with implicit/library lookup present.\n"
                + "\n".join(msg_parts) +
                "\nPrefixing changes atomtype names and can break [ bondtypes ]/[ angletypes ]/[ dihedraltypes ] "
                "matching unless parameters are baked or prefixed type-tables are injected.\n"
                "Use --prefix-implicit-topology-policy=bake (preferred) or inject, or set policy=warn to continue "
                "without prefixing."
            )

    @staticmethod
    def _path_identity(path: Path) -> Path:
        try:
            return path.resolve()
        except OSError:
            return path

    @staticmethod
    def _pair_key(ai: str, aj: str) -> Tuple[str, str]:
        return tuple(sorted((ai, aj), key=lambda v: (v.casefold(), v)))  # type: ignore[return-value]

    @staticmethod
    def _split_inline_comment(raw_line: str) -> Tuple[str, str]:
        idx = raw_line.find(";")
        if idx < 0:
            return raw_line, ""
        return raw_line[:idx], raw_line[idx:]

    def _collect_used_atomtypes(self, scan_paths: Iterable[Path]) -> Set[str]:
        used: Set[str] = set()
        for path in sorted(
            list(scan_paths),
            key=lambda p: (p.as_posix().casefold(), p.as_posix()),
        ):
            if not path.exists():
                continue
            sections = ItpParser.parse_sections(path.read_text())
            atoms_section = sections.get("atoms", [])
            for line in atoms_section:
                data = line.split(";")[0].strip()
                if not data or data.startswith("["):
                    continue
                tokens = data.split()
                if len(tokens) < 2:
                    continue
                if not tokens[0].isdigit():
                    continue
                used.add(tokens[1])
        return used

    def _collect_explicit_nonbond_pair_keys(self, scan_paths: Iterable[Path]) -> Set[Tuple[str, str]]:
        pair_keys: Set[Tuple[str, str]] = set()
        for path in sorted(
            list(scan_paths),
            key=lambda p: (p.as_posix().casefold(), p.as_posix()),
        ):
            if not path.exists():
                continue
            sections = ItpParser.parse_sections(path.read_text())
            for section_name in ("nonbond_params", "pairtypes"):
                for line in sections.get(section_name, []):
                    data = line.split(";")[0].strip()
                    if not data or data.startswith("["):
                        continue
                    tokens = data.split()
                    if len(tokens) < 3:
                        continue
                    ai, aj = tokens[0], tokens[1]
                    if ai and aj:
                        pair_keys.add(self._pair_key(ai, aj))
        return pair_keys

    def _collect_bonded_type_tables(self, scan_paths: Iterable[Path]) -> Dict[str, List[BondedTypeEntry]]:
        section_token_counts = {
            "bondtypes": 2,
            "angletypes": 3,
            "dihedraltypes": 4,
        }
        tables: Dict[str, List[BondedTypeEntry]] = {
            "bondtypes": [],
            "angletypes": [],
            "dihedraltypes": [],
        }
        for path in sorted(
            list(scan_paths),
            key=lambda p: (p.as_posix().casefold(), p.as_posix()),
        ):
            if not path.exists():
                continue
            sections = ItpParser.parse_sections(path.read_text())
            for section_name, atom_count in section_token_counts.items():
                for line in sections.get(section_name, []):
                    data = line.split(";")[0].strip()
                    if not data or data.startswith("["):
                        continue
                    tokens = data.split()
                    if len(tokens) < atom_count + 2:
                        continue
                    atomtypes = tuple(tokens[:atom_count])
                    funct = tokens[atom_count]
                    params = tuple(tokens[atom_count + 1 :])
                    if not params:
                        continue
                    tables[section_name].append(
                        BondedTypeEntry(
                            section=section_name,
                            atomtypes=atomtypes,
                            funct=funct,
                            params=params,
                            source_file=str(path),
                            wildcard=any(token in WILDCARD_TYPE_TOKENS for token in atomtypes),
                        )
                    )
        return tables

    def _lookup_bonded_params(
        self,
        table_section: str,
        atomtypes: Tuple[str, ...],
        funct: str,
        tables: Dict[str, List[BondedTypeEntry]],
        allow_wildcards: bool,
        exact_only: bool,
    ) -> Tuple[Optional[Tuple[str, ...]], str]:
        entries = tables.get(table_section, [])
        matches: List[BondedTypeEntry] = []
        for entry in entries:
            if entry.funct != funct:
                continue
            if exact_only and entry.wildcard:
                continue
            if len(entry.atomtypes) != len(atomtypes):
                continue
            candidates = [entry.atomtypes]
            # Bonded lookups are typically symmetric with reversed ordering.
            candidates.append(tuple(reversed(entry.atomtypes)))
            for candidate in candidates:
                ok = True
                for a, c in zip(atomtypes, candidate):
                    if allow_wildcards and c in WILDCARD_TYPE_TOKENS:
                        continue
                    if a != c:
                        ok = False
                        break
                if ok:
                    matches.append(entry)
                    break
        if not matches:
            return None, "missing"
        unique_params = sorted(
            {m.params for m in matches},
            key=lambda p: tuple(str(x) for x in p),
        )
        if len(unique_params) != 1:
            return None, "ambiguous"
        return unique_params[0], "ok"

    def _build_atoms_type_map(self, content: str) -> Dict[int, str]:
        atom_map: Dict[int, str] = {}
        sections = ItpParser.parse_sections(content)
        for line in sections.get("atoms", []):
            data = line.split(";")[0].strip()
            if not data or data.startswith("["):
                continue
            tokens = data.split()
            if len(tokens) < 2:
                continue
            if not tokens[0].isdigit():
                continue
            atom_map[int(tokens[0])] = tokens[1]
        return atom_map

    def _bake_implicit_sections_in_content(
        self,
        content: str,
        source_path: Path,
        tables: Dict[str, List[BondedTypeEntry]],
    ) -> Tuple[str, Dict[str, Any]]:
        section_specs = {
            "bonds": (2, "bondtypes"),
            "angles": (3, "angletypes"),
            "dihedrals": (4, "dihedraltypes"),
            "impropers": (4, "dihedraltypes"),
        }
        stats: Dict[str, Any] = {
            "baked_bonds": 0,
            "baked_angles": 0,
            "baked_dihedrals": 0,
            "unresolved_bonds": 0,
            "unresolved_angles": 0,
            "unresolved_dihedrals": 0,
            "unresolved_items": [],
        }

        atom_map = self._build_atoms_type_map(content)
        lines = content.splitlines()
        result_lines: List[str] = []
        current_section = ""

        for line_no, raw_line in enumerate(lines, 1):
            match = ItpParser.SECTION_RE.match(raw_line)
            if match:
                current_section = match.group(1).lower().strip()
                result_lines.append(raw_line)
                continue

            spec = section_specs.get(current_section)
            if spec is None:
                result_lines.append(raw_line)
                continue

            atom_count, table_section = spec
            pre_comment, comment = self._split_inline_comment(raw_line)
            data = pre_comment.strip()
            if not data:
                result_lines.append(raw_line)
                continue
            tokens = data.split()
            # Implicit lookup format: atom indices + funct (no explicit params)
            if len(tokens) != atom_count + 1:
                result_lines.append(raw_line)
                continue

            try:
                atom_indices = [int(tok) for tok in tokens[:atom_count]]
            except ValueError:
                result_lines.append(raw_line)
                continue
            funct = tokens[atom_count]
            atomtypes: List[str] = []
            missing_atom = False
            for idx in atom_indices:
                atomtype = atom_map.get(idx)
                if atomtype is None:
                    missing_atom = True
                    break
                atomtypes.append(atomtype)
            if missing_atom:
                if current_section == "bonds":
                    stats["unresolved_bonds"] += 1
                elif current_section == "angles":
                    stats["unresolved_angles"] += 1
                else:
                    stats["unresolved_dihedrals"] += 1
                stats["unresolved_items"].append(
                    {
                        "file": str(source_path),
                        "line": line_no,
                        "section": current_section,
                        "reason": "missing_atom_index_mapping",
                    }
                )
                result_lines.append(raw_line)
                continue

            exact_only = current_section in {"dihedrals", "impropers"}
            params, lookup_reason = self._lookup_bonded_params(
                table_section=table_section,
                atomtypes=tuple(atomtypes),
                funct=funct,
                tables=tables,
                allow_wildcards=not exact_only,
                exact_only=exact_only,
            )
            if params is None:
                if current_section == "bonds":
                    stats["unresolved_bonds"] += 1
                elif current_section == "angles":
                    stats["unresolved_angles"] += 1
                else:
                    stats["unresolved_dihedrals"] += 1
                stats["unresolved_items"].append(
                    {
                        "file": str(source_path),
                        "line": line_no,
                        "section": current_section,
                        "reason": lookup_reason,
                        "atomtypes": list(atomtypes),
                        "funct": funct,
                    }
                )
                result_lines.append(raw_line)
                continue

            leading = re.match(r"^(\s*)", pre_comment).group(1) if pre_comment else ""
            rewritten = leading + " ".join(tokens + list(params))
            if comment:
                rewritten = rewritten.rstrip() + " " + comment.lstrip()
            result_lines.append(rewritten)
            if current_section == "bonds":
                stats["baked_bonds"] += 1
            elif current_section == "angles":
                stats["baked_angles"] += 1
            else:
                stats["baked_dihedrals"] += 1

        return "\n".join(result_lines) + "\n", stats

    def _write_prefix_injected_types(
        self,
        output_dir: Path,
        run_id: str,
        tables: Dict[str, List[BondedTypeEntry]],
        used_atomtypes: Set[str],
    ) -> Tuple[Optional[Path], Optional[Path], Dict[str, Any]]:
        summary: Dict[str, Any] = {
            "applied": False,
            "injected_counts": {"bondtypes": 0, "angletypes": 0, "dihedraltypes": 0},
            "wildcard_dihedral_entries_skipped": 0,
            "used_atomtype_count": len(used_atomtypes),
        }
        if not self.effective_atomtype_prefix:
            summary["reason"] = "prefix_disabled"
            return None, None, summary

        output_sections: Dict[str, List[BondedTypeEntry]] = {
            "bondtypes": [],
            "angletypes": [],
            "dihedraltypes": [],
        }
        for section_name in ("bondtypes", "angletypes", "dihedraltypes"):
            for entry in tables.get(section_name, []):
                if entry.wildcard:
                    if section_name == "dihedraltypes":
                        summary["wildcard_dihedral_entries_skipped"] += 1
                    continue
                if not all(atomtype in used_atomtypes for atomtype in entry.atomtypes):
                    continue
                output_sections[section_name].append(entry)

        total_entries = sum(len(v) for v in output_sections.values())
        if total_entries == 0:
            summary["reason"] = "no_exact_entries_for_used_types"
            return None, None, summary

        versioned = output_dir / f"prefix_injected_types_{run_id}.itp"
        current = output_dir / "prefix_injected_types_current.itp"
        lines: List[str] = [
            "; Generated by ITP Sanitizer: prefix implicit-topology inject policy",
            f"; Prefix: {self.effective_atomtype_prefix}",
            f"; Used atomtypes: {len(used_atomtypes)}",
            "",
        ]
        for section_name in ("bondtypes", "angletypes", "dihedraltypes"):
            entries = sorted(
                output_sections[section_name],
                key=lambda e: (
                    tuple((t.casefold(), t) for t in e.atomtypes),
                    e.funct,
                    tuple((p.casefold(), p) for p in e.params),
                    Path(e.source_file).as_posix().casefold(),
                    Path(e.source_file).as_posix(),
                ),
            )
            if not entries:
                continue
            lines.append(f"[ {section_name} ]")
            for entry in entries:
                prefixed_types = [f"{self.effective_atomtype_prefix}{t}" for t in entry.atomtypes]
                lines.append(
                    " ".join(prefixed_types + [entry.funct] + list(entry.params))
                )
            lines.append("")
            summary["injected_counts"][section_name] = len(entries)

        versioned.write_text("\n".join(lines).rstrip() + "\n")
        if current.exists() or current.is_symlink():
            current.unlink()
        try:
            current.symlink_to(versioned.name)
        except OSError:
            shutil.copy2(versioned, current)

        summary["applied"] = True
        summary["reason"] = "generated"
        summary["versioned_path"] = str(versioned)
        summary["current_path"] = str(current)
        return versioned, current, summary

    def _prepare_prefix_policy(
        self,
        scan_paths: List[Path],
        sanitize_paths: List[Path],
        output_dir: Path,
        run_id: str,
        policy: str,
        reason: Optional[str],
    ) -> Dict[str, Any]:
        summary: Dict[str, Any] = {
            "requested_prefix": self.atomtype_prefix,
            "effective_prefix": self.atomtype_prefix,
            "policy": policy,
            "reason": reason,
            "applied_strategy": "none",
            "implicit_files": [],
            "type_based_files": [],
            "baked_bonds": 0,
            "baked_angles": 0,
            "baked_dihedrals": 0,
            "unresolved_bonds": 0,
            "unresolved_angles": 0,
            "unresolved_dihedrals": 0,
            "unresolved_items": [],
            "inject_summary": {},
        }
        self._content_overrides = {}
        self.prefix_injected_types_path = None
        self.prefix_injected_types_current_path = None
        self.effective_atomtype_prefix = self.atomtype_prefix

        if not self.atomtype_prefix:
            summary["applied_strategy"] = "prefix_not_requested"
            summary["effective_prefix"] = None
            return summary

        implicit_files: List[Tuple[Path, List[str]]] = []
        type_based_files: List[Tuple[Path, List[str]]] = []
        for itp_path in scan_paths:
            if not itp_path.exists():
                continue
            content = itp_path.read_text()
            implicit = ItpParser.detect_implicit_topology(content)
            if implicit:
                implicit_files.append((itp_path, implicit))
            sections = ItpParser.parse_sections(content)
            found_type_based = [s for s in TYPE_BASED_SECTIONS if s in sections]
            if found_type_based:
                type_based_files.append((itp_path, found_type_based))

        summary["implicit_files"] = [
            {"path": str(path), "sections": sections}
            for path, sections in implicit_files
        ]
        summary["type_based_files"] = [
            {"path": str(path), "sections": sections}
            for path, sections in type_based_files
        ]

        if not implicit_files and not type_based_files:
            summary["applied_strategy"] = "safe_no_rewrite_needed"
            return summary

        if policy == "error":
            self.validate_prefix_safety(scan_paths)
            summary["applied_strategy"] = "error"
            return summary

        if policy == "warn":
            print("  [WARN] ===========================================================")
            print("  [WARN] Prefix implicit-topology policy=warn: disabling atomtype prefix.")
            print("  [WARN] Name-collision risk remains for mixed namespace topologies.")
            print("  [WARN] Suggested fixes:")
            print("  [WARN]   - Use --prefix-implicit-topology-policy=bake (preferred)")
            print("  [WARN]   - Or use --prefix-implicit-topology-policy=inject")
            print("  [WARN] ===========================================================")
            self.effective_atomtype_prefix = None
            summary["effective_prefix"] = None
            summary["applied_strategy"] = "warn_disable_prefix"
            return summary

        tables = self._collect_bonded_type_tables(scan_paths)
        if policy == "bake":
            unresolved: List[Dict[str, Any]] = []
            for source_path in sorted(
                sanitize_paths,
                key=lambda p: (p.as_posix().casefold(), p.as_posix()),
            ):
                if not source_path.exists():
                    continue
                original = source_path.read_text()
                rewritten, stats = self._bake_implicit_sections_in_content(
                    original,
                    source_path,
                    tables,
                )
                if rewritten != original:
                    self._content_overrides[self._path_identity(source_path)] = rewritten
                summary["baked_bonds"] += int(stats.get("baked_bonds", 0))
                summary["baked_angles"] += int(stats.get("baked_angles", 0))
                summary["baked_dihedrals"] += int(stats.get("baked_dihedrals", 0))
                summary["unresolved_bonds"] += int(stats.get("unresolved_bonds", 0))
                summary["unresolved_angles"] += int(stats.get("unresolved_angles", 0))
                summary["unresolved_dihedrals"] += int(stats.get("unresolved_dihedrals", 0))
                unresolved.extend(list(stats.get("unresolved_items", [])))
            unresolved_sorted = sorted(
                unresolved,
                key=lambda item: (
                    str(item.get("file", "")).casefold(),
                    int(item.get("line", 0)),
                    str(item.get("section", "")),
                ),
            )
            summary["unresolved_items"] = unresolved_sorted

            if summary["unresolved_bonds"] or summary["unresolved_angles"]:
                preview = unresolved_sorted[:6]
                details = "\n".join(
                    f"  - {item.get('file')}:{item.get('line')} [{item.get('section')}] "
                    f"reason={item.get('reason')}"
                    for item in preview
                )
                raise PrefixUnsafeError(
                    "prefix-implicit-topology-policy=bake could not resolve all implicit "
                    "bond/angle parameters uniquely.\n"
                    + details
                    + ("\n  - ... (truncated)" if len(unresolved_sorted) > len(preview) else "")
                )
            if summary["unresolved_dihedrals"]:
                print("  [WARN] ===========================================================")
                print(
                    "  [WARN] Prefix bake left unresolved implicit dihedrals "
                    f"({summary['unresolved_dihedrals']})."
                )
                print(
                    "  [WARN] These remain parameter-by-lookup and may need inject policy or "
                    "manual parameterization."
                )
                print("  [WARN] ===========================================================")
            summary["applied_strategy"] = "bake"
            return summary

        if policy == "inject":
            used_atomtypes = self._collect_used_atomtypes(scan_paths)
            injected_versioned, injected_current, inject_summary = self._write_prefix_injected_types(
                output_dir=output_dir,
                run_id=run_id,
                tables=tables,
                used_atomtypes=used_atomtypes,
            )
            self.prefix_injected_types_path = injected_versioned
            self.prefix_injected_types_current_path = injected_current
            summary["inject_summary"] = inject_summary
            summary["applied_strategy"] = (
                "inject" if inject_summary.get("applied") else "inject_no_entries"
            )
            if inject_summary.get("wildcard_dihedral_entries_skipped", 0):
                print(
                    "  [WARN] Prefix inject skipped wildcard dihedraltypes entries; "
                    "wildcard-only coverage may still require manual review."
                )
            return summary

        raise PrefixUnsafeError(
            f"Unsupported prefix_implicit_topology_policy: {policy!r}"
        )

    def _collect_defaults(self, itp_paths: List[Path]) -> None:
        """Collect [ defaults ] entries from paths in the provided deterministic order."""
        seen: Set[Path] = set()
        for itp_path in itp_paths:
            if not itp_path.exists():
                continue
            try:
                key = itp_path.resolve()
            except OSError:
                key = itp_path
            if key in seen:
                continue
            seen.add(key)
            defaults = ItpParser.extract_defaults_from_file(itp_path)
            if defaults:
                self.defaults_entries.append(defaults)

    def _build_include_dirs(
        self,
        itp_paths: List[Path],
        output_dir: Path,
        ctx: Optional["PipelineContext"] = None,
        include_dirs: Optional[List[Path]] = None,
        include_roots: Optional[List[Path]] = None,
    ) -> List[Path]:
        """
        Build include search paths for #include resolution (GROMACS -I style).
        
        Order (best-effort):
        1) itp_path.parent for each scanned file
        2) common system dirs (top, itp, htpolynet itp)
        3) forcefield dir (if known)
        4) sanitized output dir (if exists)
        5) user-provided include dirs
        6) gromacs share dirs from ctx (if configured)
        """
        paths: List[Path] = []
        seen: Set[Path] = set()
        project_root_path: Optional[Path] = None
        if ctx is not None and getattr(ctx, "project_root", None):
            try:
                project_root_path = Path(ctx.project_root).resolve()
            except OSError:
                project_root_path = Path(ctx.project_root)

        def _resolve_ctx_path(raw_path: Path) -> Path:
            if raw_path.is_absolute() or project_root_path is None:
                return raw_path
            return project_root_path / raw_path

        def _add_path(p: Path) -> None:
            try:
                resolved = p.resolve()
            except OSError:
                return
            if resolved in seen:
                return
            if resolved.exists() and resolved.is_dir():
                seen.add(resolved)
                paths.append(resolved)

        # 1) Parents of scanned files
        for itp_path in itp_paths:
            _add_path(itp_path.parent)
        if include_roots:
            for root in include_roots:
                _add_path(root.parent)

        # 2) Common system dirs
        if ctx is not None:
            try:
                system_gromacs_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs")
                system_htpolynet_dir = ctx.get_input_path("systems", ctx.system_id, "htpolynet")
                _add_path(system_gromacs_dir / "top")
                _add_path(system_gromacs_dir / "itp")
                _add_path(system_htpolynet_dir / "itp")
            except Exception:
                pass

            # 3) Forcefield dir (if known)
            ff_dir = getattr(ctx, "ff_dir", None)
            if ff_dir:
                _add_path(Path(ff_dir) / "gromacs")
            else:
                # Best-effort default from IN/forcefield
                ff_base = ctx.ff.lower()
                candidates = [
                    ctx.get_input_path("forcefield", ff_base.replace("-", "")),
                    ctx.get_input_path("forcefield", ff_base),
                ]
                for candidate in candidates:
                    if candidate.exists():
                        _add_path(candidate / "gromacs")
                        break

        # 4) Sanitized output dir (if exists)
        _add_path(output_dir / "itp_sanitized_current")

        # 5) User-provided include dirs (CLI)
        if ctx is not None and getattr(ctx, "itp_include_dirs", None):
            for d in ctx.itp_include_dirs.split(","):
                d = d.strip()
                if d:
                    _add_path(_resolve_ctx_path(Path(d)))
        if include_dirs:
            for d in include_dirs:
                _add_path(_resolve_ctx_path(d))

        # 6) GROMACS share dirs (only if available from ctx)
        if ctx is not None:
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
                    for item in value:
                        _add_path(_resolve_ctx_path(Path(item)))
                elif isinstance(value, str):
                    for item in value.split(","):
                        item = item.strip()
                        if item:
                            _add_path(_resolve_ctx_path(Path(item)))

        return paths
    
    def _determine_output_charge(self, name: str) -> float:
        """
        Determine the charge to output for an atomtype.
        
        Policy:
        - If all charges agree within tolerance: use that value
        - Otherwise: keep the first deterministic source value (warning already issued)
        """
        if name not in self._charge_by_atomtype:
            return 0.0
        
        charges = [c for c, _ in self._charge_by_atomtype[name]]
        if not charges:
            return 0.0
        
        # Check if all agree
        first = charges[0]
        if all(abs(c - first) < CHARGE_TOLERANCE for c in charges):
            return first
        
        # Disagreement: keep deterministic primary source value.
        return first
    
    def generate_combined_atomtypes(self) -> str:
        """
        Generate combined atomtypes ITP content.
        
        Charge output is deterministic:
        - If all sources agree on charge: use that value
        - If disagreement: keep first deterministic source value (and warning is logged)
        
        Task E: Labels LJ columns correctly based on comb-rule.
        """
        # Determine comb-rule for column labeling
        comb_rule = self.canonical_defaults.comb_rule if self.canonical_defaults else 2
        
        if comb_rule == 1:
            lj_header = "; name    mass      charge   ptype  C6           C12"
            lj_comment = f"; Comb-rule 1: LJ columns are C6 and C12 coefficients"
        else:
            lj_header = "; name    mass      charge   ptype  sigma        epsilon"
            lj_comment = f"; Comb-rule {comb_rule}: LJ columns are sigma (nm) and epsilon (kJ/mol)"
        
        lines = [
            "; Combined atomtypes generated by ITP Sanitizer",
            f"; Source files: {len(self.source_files)}",
            f"; Atomtype count: {len(self.registry)}",
            lj_comment,
            "",
            "[ atomtypes ]",
            lj_header,
        ]
        
        # Sort atomtypes for reproducibility
        for name in sorted(self.registry.keys()):
            entry = self.registry[name]
            
            # Apply prefix if configured
            output_name = name
            if self.effective_atomtype_prefix:
                output_name = f"{self.effective_atomtype_prefix}{name}"
            
            # Determine charge (deterministic)
            output_charge = self._determine_output_charge(name)
            
            # Handle ptype=None (Task D)
            output_ptype = entry.ptype if entry.ptype else "A"
            
            line = (
                f"{output_name:8s} {entry.mass:10.4f} {output_charge:8.4f}  "
                f"{output_ptype}  {entry.sigma:12.6e}  {entry.epsilon:12.6e}"
            )
            
            # Add source comment
            if entry.source_files:
                sources = ", ".join(Path(f).name for f in entry.source_files[:3])
                line += f"  ; from: {sources}"
            
            lines.append(line)
        
        return '\n'.join(lines) + '\n'

    def _generate_cross_group_nonbond_params(
        self,
        output_dir: Path,
        run_id: str,
        source_defaults_map: Dict[Path, DefaultsEntry],
        defaults_report: Dict[str, Any],
        scan_paths: List[Path],
        ctx: Optional["PipelineContext"] = None,
    ) -> Tuple[Optional[Path], Optional[Path], Dict[str, Any]]:
        allow_mixed = bool(getattr(ctx, "allow_mixed_defaults", False)) if ctx else False
        policy = (
            getattr(ctx, "mixed_defaults_cross_group_policy", None) if ctx else None
        )
        if policy is None:
            policy = "warn" if allow_mixed else "off"
        rule = getattr(ctx, "mixed_defaults_cross_group_rule", None) if ctx else None
        max_pairs = int(
            getattr(ctx, "mixed_defaults_cross_group_max_pairs", 20000) if ctx else 20000
        )
        reason = (
            (getattr(ctx, "mixed_defaults_cross_group_reason", None) if ctx else None) or None
        )

        summary: Dict[str, Any] = {
            "policy": policy,
            "rule": rule,
            "max_pairs": max_pairs,
            "reason": reason,
            "applied": False,
            "pair_count": 0,
            "candidate_pair_count": 0,
            "skipped_existing_count": 0,
            "skipped_invalid_count": 0,
            "preexisting_pair_count": 0,
        }
        if policy == "off":
            summary["status"] = "disabled"
            return None, None, summary

        mixed_detected = bool(defaults_report.get("mixed_defaults_detected", False))
        if not mixed_detected:
            summary["status"] = "no_mixed_defaults"
            return None, None, summary

        classification = str(defaults_report.get("classification", "unknown"))
        summary["classification"] = classification
        if classification != "comb-rule-only":
            if policy == "generate":
                raise ItpSanitizerError(
                    "mixed-defaults-cross-group-policy=generate is only supported for "
                    f"classification='comb-rule-only' (got '{classification}')."
                )
            summary["status"] = f"incompatible_classification:{classification}"
            if policy == "warn":
                print(
                    "  [WARN] mixed-defaults cross-group policy=warn requested, but classification "
                    f"is '{classification}'. No cross-group pair generation attempted."
                )
            return None, None, summary

        primary_defaults = self.canonical_defaults
        if primary_defaults is None:
            summary["status"] = "missing_primary_defaults"
            return None, None, summary
        if primary_defaults.nbfunc != 1:
            raise ItpSanitizerError(
                "mixed-defaults cross-group generation currently supports nbfunc=1 (LJ 12-6) only."
            )
        if primary_defaults.comb_rule not in {2, 3}:
            raise ItpSanitizerError(
                "mixed-defaults cross-group generation requires canonical comb-rule 2 or 3 "
                f"(got {primary_defaults.comb_rule})."
            )

        if rule is None:
            rule = "lorentz-berthelot" if primary_defaults.comb_rule == 2 else "geometric"
        summary["rule"] = rule

        primary_signature = defaults_report.get("primary_defaults_signature")
        secondary_signatures: List[str] = list(defaults_report.get("secondary_signatures", []))
        if not primary_signature or not secondary_signatures:
            summary["status"] = "no_secondary_signatures"
            return None, None, summary

        source_signature_map: Dict[Path, str] = {}
        for src_path, defaults in source_defaults_map.items():
            source_signature_map[self._path_identity(src_path)] = defaults.signature()

        used_atomtypes = self._collect_used_atomtypes(scan_paths)
        if not used_atomtypes:
            used_atomtypes = set(self.registry.keys())

        primary_types: Set[str] = set()
        secondary_types: Set[str] = set()
        for atomtype_name in sorted(self.registry.keys(), key=lambda s: (s.casefold(), s)):
            if atomtype_name not in used_atomtypes:
                continue
            entry = self.registry[atomtype_name]
            origin_signatures: Set[str] = set()
            for source_file in entry.source_files:
                source_key = self._path_identity(Path(source_file))
                signature = source_signature_map.get(source_key)
                if signature:
                    origin_signatures.add(signature)
            if not origin_signatures:
                continue
            if primary_signature in origin_signatures:
                primary_types.add(atomtype_name)
                continue
            if any(sig in secondary_signatures for sig in origin_signatures):
                secondary_types.add(atomtype_name)

        primary_list = sorted(primary_types, key=lambda s: (s.casefold(), s))
        secondary_list = sorted(secondary_types, key=lambda s: (s.casefold(), s))
        summary["primary_type_count"] = len(primary_list)
        summary["secondary_type_count"] = len(secondary_list)
        if not primary_list or not secondary_list:
            summary["status"] = "insufficient_cross_group_types"
            return None, None, summary

        candidate_pairs = len(primary_list) * len(secondary_list)
        summary["candidate_pair_count"] = candidate_pairs
        if candidate_pairs > max_pairs:
            msg = (
                "Refusing mixed-defaults cross-group [ nonbond_params ] generation: "
                f"{candidate_pairs} candidate pairs exceeds max {max_pairs}. "
                "Increase --mixed-defaults-cross-group-max-pairs if this is intentional."
            )
            if policy == "generate":
                raise ItpSanitizerError(msg)
            print(f"  [WARN] {msg}")
            summary["status"] = "over_limit_warn_only"
            return None, None, summary

        if policy == "warn":
            print("  [WARN] ===========================================================")
            print(
                "  [WARN] Mixed comb-rule topology detected (2 vs 3). "
                "Cross-group LJ currently falls back to one global [ defaults ] rule."
            )
            print(
                "  [WARN] This can distort coordination, RDFs, and transport properties "
                "for electrolyte/polymer systems."
            )
            print(
                "  [WARN] Recommendation: set --mixed-defaults-cross-group-policy generate "
                "with --mixed-defaults-cross-group-reason."
            )
            print(
                f"  [WARN] Candidate cross-group pairs in this system: {candidate_pairs} "
                f"(cap={max_pairs})."
            )
            print("  [WARN] ===========================================================")
            summary["status"] = "warn_only"
            return None, None, summary

        if policy == "generate":
            reason_clean = (reason or "").strip()
            if not reason_clean:
                raise ItpSanitizerError(
                    "mixed-defaults cross-group generate policy requires a non-empty reason."
                )
            summary["reason"] = reason_clean

        existing_pairs = self._collect_explicit_nonbond_pair_keys(scan_paths)
        summary["preexisting_pair_count"] = len(existing_pairs)

        param_lines: List[str] = []
        skipped_existing = 0
        skipped_invalid = 0
        for ai in primary_list:
            ei = self.registry.get(ai)
            if ei is None:
                continue
            sigma_i = float(ei.sigma)
            eps_i = float(ei.epsilon)
            for aj in secondary_list:
                ej = self.registry.get(aj)
                if ej is None:
                    continue
                pair_key = self._pair_key(ai, aj)
                if pair_key in existing_pairs:
                    skipped_existing += 1
                    continue

                sigma_j = float(ej.sigma)
                eps_j = float(ej.epsilon)
                # For generated sigma/epsilon pairs we require non-negative values.
                if sigma_i < 0 or sigma_j < 0 or eps_i < 0 or eps_j < 0:
                    skipped_invalid += 1
                    continue
                if rule == "geometric":
                    sigma_ij = math.sqrt(sigma_i * sigma_j)
                else:
                    sigma_ij = 0.5 * (sigma_i + sigma_j)
                eps_ij = math.sqrt(eps_i * eps_j)

                out_ai = (
                    f"{self.effective_atomtype_prefix}{ai}"
                    if self.effective_atomtype_prefix
                    else ai
                )
                out_aj = (
                    f"{self.effective_atomtype_prefix}{aj}"
                    if self.effective_atomtype_prefix
                    else aj
                )
                param_lines.append(
                    f"{out_ai:16s} {out_aj:16s} 1 {sigma_ij:.12e} {eps_ij:.12e}"
                )

        if not param_lines:
            summary["status"] = "no_pairs_generated"
            summary["skipped_existing_count"] = skipped_existing
            summary["skipped_invalid_count"] = skipped_invalid
            return None, None, summary

        versioned_path = output_dir / f"nonbond_params_cross_group_{run_id}.itp"
        current_path = output_dir / "nonbond_params_cross_group_current.itp"
        lines = [
            "; Generated by ITP Sanitizer: mixed-defaults cross-group LJ override",
            f"; policy={policy}",
            f"; rule={rule}",
            f"; reason={summary.get('reason')}",
            f"; canonical_comb_rule={primary_defaults.comb_rule}",
            "; For comb-rule 2/3, V/W are sigma (nm) and epsilon (kJ/mol).",
            "",
            "[ nonbond_params ]",
            "; i               j               func   V(sigma_nm)      W(epsilon_kJ_mol)",
        ]
        lines.extend(param_lines)
        versioned_path.write_text("\n".join(lines) + "\n")

        if current_path.exists() or current_path.is_symlink():
            current_path.unlink()
        try:
            current_path.symlink_to(versioned_path.name)
        except OSError:
            shutil.copy2(versioned_path, current_path)

        summary.update(
            {
                "applied": True,
                "status": "generated",
                "pair_count": len(param_lines),
                "skipped_existing_count": skipped_existing,
                "skipped_invalid_count": skipped_invalid,
                "versioned_path": str(versioned_path),
                "current_path": str(current_path),
            }
        )
        return versioned_path, current_path, summary

    def _generate_secondary_nonbond_params(
        self,
        output_dir: Path,
        run_id: str,
        source_defaults_map: Dict[Path, DefaultsEntry],
        defaults_report: Dict[str, Any],
        preserve_within_group: bool = False,
    ) -> Tuple[Optional[Path], Optional[Path], Dict[str, Any]]:
        """
        Optionally generate explicit [ nonbond_params ] for secondary comb-rule groups.
        """
        summary: Dict[str, Any] = {
            "requested": bool(preserve_within_group),
            "applied": False,
            "reason": "flag_disabled",
            "pair_count": 0,
            "atomtype_count": 0,
            "secondary_groups": [],
        }
        if not preserve_within_group:
            return None, None, summary

        mixed_detected = bool(defaults_report.get("mixed_defaults_detected", False))
        if not mixed_detected:
            summary["reason"] = "no_mixed_defaults"
            return None, None, summary

        classification = str(defaults_report.get("classification", "unknown"))
        if classification != "comb-rule-only":
            summary["reason"] = f"incompatible_classification:{classification}"
            print(
                "  [WARN] --mixed-defaults-preserve-within-group requested, but mixed defaults "
                f"classification is '{classification}' (needs comb-rule-only). Skipping."
            )
            return None, None, summary

        primary = self.canonical_defaults
        if primary is None:
            summary["reason"] = "missing_primary_defaults"
            return None, None, summary
        if primary.nbfunc != 1:
            summary["reason"] = f"unsupported_nbfunc:{primary.nbfunc}"
            print(
                "  [WARN] --mixed-defaults-preserve-within-group currently supports nbfunc=1 "
                f"(LJ 12-6) only. Primary nbfunc={primary.nbfunc}; skipping."
            )
            return None, None, summary

        primary_signature = defaults_report.get("primary_defaults_signature")
        secondary_signatures: List[str] = list(defaults_report.get("secondary_signatures", []))
        if not primary_signature or not secondary_signatures:
            summary["reason"] = "no_secondary_signatures"
            return None, None, summary

        source_signature_map: Dict[Path, str] = {}
        for src_path, defaults in source_defaults_map.items():
            try:
                key = src_path.resolve()
            except OSError:
                key = src_path
            source_signature_map[key] = defaults.signature()

        group_comb_rule: Dict[str, int] = {}
        for group in defaults_report.get("groups", []):
            signature = str(group.get("signature", ""))
            if not signature:
                continue
            comb_rule = group.get("comb_rule")
            if isinstance(comb_rule, int):
                group_comb_rule[signature] = comb_rule

        secondary_group_atomtypes: Dict[str, Set[str]] = {}
        for atomtype_name in sorted(self.registry.keys(), key=str.casefold):
            entry = self.registry[atomtype_name]
            origin_signatures: Set[str] = set()
            for source_file in entry.source_files:
                source_path = Path(source_file)
                try:
                    source_key = source_path.resolve()
                except OSError:
                    source_key = source_path
                signature = source_signature_map.get(source_key)
                if signature is not None:
                    origin_signatures.add(signature)

            if not origin_signatures:
                continue
            if primary_signature in origin_signatures:
                continue

            secondary_hits = sorted(
                [sig for sig in origin_signatures if sig in secondary_signatures],
                key=str.casefold,
            )
            if not secondary_hits:
                continue
            chosen_signature = secondary_hits[0]
            secondary_group_atomtypes.setdefault(chosen_signature, set()).add(atomtype_name)

        if not secondary_group_atomtypes:
            summary["reason"] = "no_secondary_atomtypes"
            print(
                "  [WARN] --mixed-defaults-preserve-within-group requested, but no atomtypes "
                "were uniquely attributable to a secondary defaults group."
            )
            return None, None, summary

        param_lines: List[str] = []
        pair_count = 0
        all_atomtypes: Set[str] = set()
        for signature in sorted(secondary_group_atomtypes.keys(), key=str.casefold):
            atomtypes = sorted(secondary_group_atomtypes[signature], key=str.casefold)
            all_atomtypes.update(atomtypes)
            comb_rule = group_comb_rule.get(signature)
            if comb_rule not in {2, 3}:
                print(
                    f"  [WARN] Secondary signature {signature} has unsupported comb-rule={comb_rule}; "
                    "expected 2 or 3 for sigma/epsilon. Skipping this group."
                )
                continue

            summary["secondary_groups"].append(
                {
                    "signature": signature,
                    "comb_rule": comb_rule,
                    "atomtypes": atomtypes,
                }
            )
            param_lines.append(
                f"; secondary signature: {signature} (comb-rule={comb_rule}, atomtypes={len(atomtypes)})"
            )

            for i, ai in enumerate(atomtypes):
                entry_i = self.registry.get(ai)
                if entry_i is None:
                    continue
                sigma_i = float(entry_i.sigma)
                epsilon_i = float(entry_i.epsilon)
                for aj in atomtypes[i:]:
                    entry_j = self.registry.get(aj)
                    if entry_j is None:
                        continue
                    sigma_j = float(entry_j.sigma)
                    epsilon_j = float(entry_j.epsilon)

                    if sigma_i < 0 or sigma_j < 0 or epsilon_i < 0 or epsilon_j < 0:
                        print(
                            f"  [WARN] Skipping nonbond_params pair ({ai}, {aj}) due to negative "
                            "sigma/epsilon values in secondary group."
                        )
                        continue

                    if comb_rule == 3:
                        sigma_ij = math.sqrt(sigma_i * sigma_j)
                    else:
                        sigma_ij = 0.5 * (sigma_i + sigma_j)
                    epsilon_ij = math.sqrt(epsilon_i * epsilon_j)

                    c6 = 4.0 * epsilon_ij * (sigma_ij ** 6)
                    c12 = 4.0 * epsilon_ij * (sigma_ij ** 12)
                    out_ai = (
                        f"{self.effective_atomtype_prefix}{ai}"
                        if self.effective_atomtype_prefix
                        else ai
                    )
                    out_aj = (
                        f"{self.effective_atomtype_prefix}{aj}"
                        if self.effective_atomtype_prefix
                        else aj
                    )
                    param_lines.append(
                        f"{out_ai:16s} {out_aj:16s} 1 {c6:.12e} {c12:.12e}"
                    )
                    pair_count += 1

        if pair_count == 0:
            summary["reason"] = "no_pairs_generated"
            return None, None, summary

        versioned_path = output_dir / f"nonbond_params_secondary_{run_id}.itp"
        current_path = output_dir / "nonbond_params_secondary_current.itp"
        content_lines = [
            "; Generated by ITP Sanitizer: secondary comb-rule within-group preservation",
            "; WARNING: Cross-group LJ interactions still use the global primary comb-rule.",
            ";          Manual cross-group nonbonded reparameterization may be required for accuracy.",
            "",
            "[ nonbond_params ]",
            "; i               j               func   C6             C12",
        ]
        content_lines.extend(param_lines)
        versioned_path.write_text("\n".join(content_lines) + "\n")

        if current_path.exists() or current_path.is_symlink():
            current_path.unlink()
        try:
            current_path.symlink_to(versioned_path.name)
        except OSError:
            shutil.copy2(versioned_path, current_path)

        print("  [WARN] ===========================================================")
        print(
            "  [WARN] mixed-defaults preserve-within-group generated explicit "
            f"[ nonbond_params ] for {pair_count} pair(s)."
        )
        print(
            "  [WARN] Cross-group LJ (secondary vs primary) is still governed by the "
            "global primary comb-rule and may need manual re-parameterization."
        )
        print("  [WARN] ===========================================================")

        summary.update(
            {
                "applied": True,
                "reason": "generated",
                "pair_count": pair_count,
                "atomtype_count": len(all_atomtypes),
                "versioned_path": str(versioned_path),
                "current_path": str(current_path),
            }
        )
        return versioned_path, current_path, summary
    
    def sanitize_itp(self, itp_path: Path) -> str:
        """
        Sanitize an ITP file by removing atomtypes and defaults sections.
        
        Args:
            itp_path: Path to ITP file
            
        Returns:
            Sanitized ITP content
        """
        content = self._content_overrides.get(
            self._path_identity(itp_path),
            itp_path.read_text(),
        )
        
        # Strip atomtypes
        content = ItpParser.strip_atomtypes_section(content)
        
        # Strip defaults (will be in system.top)
        content = ItpParser.strip_defaults_section(content)
        
        # Apply prefix if needed
        if self.effective_atomtype_prefix:
            content, _ = ItpParser.apply_atomtype_prefix(
                content,
                self.effective_atomtype_prefix,
                set(self.registry.keys()),
            )
        
        return content
    
    def run(
        self,
        itp_paths: List[Path],
        output_dir: Path,
        run_id: str,
        ff_map: Optional[Dict[str, str]] = None,
        class_map: Optional[Dict[str, str]] = None,
        ctx: Optional["PipelineContext"] = None,
        include_dirs: Optional[List[Path]] = None,
        include_roots: Optional[List[Path]] = None,
        allow_default_defaults: bool = False,
        sanitize_include_closure: bool = True,
    ) -> SanitizerResult:
        """
        Run full sanitization workflow.
        
        Args:
            itp_paths: List of ITP files to process
            output_dir: Directory for outputs
            run_id: Run identifier for versioning
            ff_map: Optional mapping of file path to forcefield
            class_map: Optional mapping of file path to class
            ctx: Optional pipeline context (for include resolution and policy)
            include_dirs: Optional include search dirs for #include resolution
            include_roots: Optional root files (e.g., system.top) to seed include closure
            allow_default_defaults: Allow fallback defaults if none found
            sanitize_include_closure: If True, sanitize include-closure files that contain
                                      [ defaults ] or [ atomtypes ] to prevent reintroduction
            
        Returns:
            SanitizerResult with output paths and statistics
        """
        print(f"  Scanning {len(itp_paths)} ITP files...")
        
        # Task J: Reset all instance state at start of run()
        # This prevents state leakage if instance is reused
        self.registry = {}
        self.source_files = []
        self.defaults_entries = []
        self.canonical_defaults = None
        self.defaults_report = {}
        self.conflicts = []
        self.overridden_conflicts = []
        self.charge_warnings = []
        self._charge_by_atomtype = {}
        self.effective_atomtype_prefix = self.atomtype_prefix
        self._content_overrides = {}
        self.prefix_policy_summary = {}
        self.prefix_injected_types_path = None
        self.prefix_injected_types_current_path = None
        
        # Sort paths for reproducibility (casefold for cross-platform stability)
        itp_paths = sorted(
            itp_paths,
            key=lambda p: (p.as_posix().casefold(), p.as_posix()),
        )

        # Detect basename collisions early (case-insensitive) to avoid silent overwrites
        basename_seen: Dict[str, Path] = {}
        for itp_path in itp_paths:
            key = itp_path.name.casefold()
            existing = basename_seen.get(key)
            if existing is not None:
                try:
                    same = existing.resolve() == itp_path.resolve()
                except OSError:
                    same = existing == itp_path
                if not same:
                    raise ItpSanitizerError(
                        "Case-insensitive ITP basename collision detected:\n"
                        f"  - {existing}\n"
                        f"  - {itp_path}\n"
                        "Rename one of the files to avoid output overwrites."
                    )
            else:
                basename_seen[key] = itp_path

        # Build include search paths (GROMACS -I style)
        include_dirs = self._build_include_dirs(
            itp_paths=itp_paths,
            output_dir=output_dir,
            ctx=ctx,
            include_dirs=include_dirs,
            include_roots=include_roots,
        )

        # Resolve include closure (from system.top or entry ITPs)
        include_diagnostics: List[IncludeResolution] = []
        include_roots = include_roots or []
        include_seed: List[Path] = []
        include_seed_seen: Set[Path] = set()

        def _path_identity(path: Path) -> Path:
            try:
                return path.resolve()
            except OSError:
                return path

        def _append_seed(path: Path) -> None:
            if not path.exists():
                return
            key = _path_identity(path)
            if key in include_seed_seen:
                return
            include_seed_seen.add(key)
            include_seed.append(path)

        # Prefer explicit topology roots first, then discovered ITP entries.
        for root in include_roots:
            _append_seed(root)
        for itp in itp_paths:
            _append_seed(itp)

        strict_includes = getattr(ctx, "strict_include_resolution", False) if ctx else False
        allow_shadowing = getattr(ctx, "allow_include_shadowing", False) if ctx else False
        include_paths: Set[Path] = set()
        include_order: List[Path] = []
        include_order_seen: Set[Path] = set()

        def _record_include_order(path: Path) -> None:
            key = _path_identity(path)
            if key in include_order_seen:
                return
            include_order_seen.add(key)
            include_order.append(path)

        for root in include_seed:
            _record_include_order(root)
            resolved_includes = ItpParser.resolve_includes(
                root,
                include_dirs=include_dirs,
                strict=strict_includes,
                diagnostics=include_diagnostics,
                check_shadowing=True,
                allow_shadowing=allow_shadowing,
            )
            for inc in resolved_includes:
                include_paths.add(inc)
                _record_include_order(inc)

        # Build full scan set (entry ITPs + include closure + include roots)
        scan_paths = set(itp_paths) | set(include_paths) | set(include_roots)
        scan_paths_list = sorted(
            scan_paths,
            key=lambda p: (p.as_posix().casefold(), p.as_posix()),
        )

        # Ensure ff_map/class_map cover scanned files (default to UNKNOWN)
        ff_map = ff_map or {}
        class_map = class_map or {}
        for p in scan_paths_list:
            if str(p) not in ff_map:
                ff_map[str(p)] = "UNKNOWN"
            if str(p) not in class_map:
                class_map[str(p)] = "unknown"

        # Determine which files will be sanitized (entry ITPs + include-closure defaults/atomtypes)
        sanitize_sources: Dict[Path, Path] = {}
        for itp_path in itp_paths:
            if itp_path.exists() and itp_path.suffix.lower() == ".itp":
                key = _path_identity(itp_path)
                if key not in sanitize_sources:
                    sanitize_sources[key] = itp_path
        if sanitize_include_closure:
            for scan_path in scan_paths_list:
                if not scan_path.exists():
                    continue
                if scan_path.suffix.lower() != ".itp":
                    continue
                if not ItpParser.has_defaults_or_atomtypes(scan_path):
                    continue
                key = _path_identity(scan_path)
                if key not in sanitize_sources:
                    sanitize_sources[key] = scan_path

        # Collect defaults across include-closure (prefer include encounter order)
        scan_lookup: Dict[Path, Path] = {}
        for scan_path in scan_paths_list:
            scan_lookup[_path_identity(scan_path)] = scan_path
        defaults_scan_order: List[Path] = []
        defaults_seen: Set[Path] = set()
        for ordered_path in include_order:
            key = _path_identity(ordered_path)
            candidate = scan_lookup.get(key)
            if candidate is None or key in defaults_seen:
                continue
            defaults_seen.add(key)
            defaults_scan_order.append(candidate)
        for scan_path in scan_paths_list:
            key = _path_identity(scan_path)
            if key in defaults_seen:
                continue
            defaults_seen.add(key)
            defaults_scan_order.append(scan_path)

        self._collect_defaults(defaults_scan_order)
        source_defaults_map: Dict[Path, DefaultsEntry] = {}
        for entry in self.defaults_entries:
            entry_path = Path(entry.source_file)
            try:
                entry_path = entry_path.resolve()
            except OSError:
                pass
            source_defaults_map[entry_path] = entry
        mixed_defaults_detected = False
        defaults_report: Dict[str, Any] = {}
        if self.defaults_entries:
            defaults_report = self.validate_defaults(
                allow_mixed_defaults=getattr(ctx, "allow_mixed_defaults", False) if ctx else False,
                allow_mixed_defaults_reason=(
                    getattr(ctx, "allow_mixed_defaults_reason", None) if ctx else None
                ),
                class_map=class_map,
            )
            mixed_defaults_detected = bool(defaults_report.get("mixed_defaults_detected", False))

        # Scan atomtypes across include-closure with known comb-rule (if any)
        comb_rule = self.canonical_defaults.comb_rule if self.canonical_defaults else None
        self.scan_itp_files(
            scan_paths_list,
            ff_map,
            class_map,
            comb_rule_override=comb_rule,
            collect_defaults=False,
        )

        # Task E: Comb-rule safety - fail-fast if no defaults but atomtypes found
        defaults_guessed = False
        if self.canonical_defaults is None and len(self.registry) > 0:
            if allow_default_defaults:
                print("  [WARN] No [ defaults ] found; using fallback defaults (nbfunc=1, comb-rule=2)")
                self.canonical_defaults = _fallback_defaults_entry()
                defaults_guessed = True
                source_defaults_map[Path(self.canonical_defaults.source_file)] = self.canonical_defaults
            else:
                scanned_list = [p.as_posix() for p in scan_paths_list]
                raise ItpSanitizerError(
                    "Missing [ defaults ] in reachable topology.\n"
                    f"Scanned {len(scanned_list)} files:\n  - "
                    + "\n  - ".join(scanned_list[:10])
                    + ("\n  - ... (truncated)" if len(scanned_list) > 10 else "")
                    + "\nProvide explicit [ defaults ] or enable --allow-default-defaults."
                )
        if not defaults_report and self.canonical_defaults is not None:
            defaults_report = {
                "mixed_defaults_detected": False,
                "classification": "fallback",
                "differing_fields": [],
                "comb_rules_present": [self.canonical_defaults.comb_rule],
                "groups": [
                    {
                        "signature": self.canonical_defaults.signature(),
                        "nbfunc": self.canonical_defaults.nbfunc,
                        "comb_rule": self.canonical_defaults.comb_rule,
                        "gen_pairs": self.canonical_defaults.gen_pairs,
                        "fudge_lj": self.canonical_defaults.fudge_lj,
                        "fudge_qq": self.canonical_defaults.fudge_qq,
                        "files": [self.canonical_defaults.source_file],
                    }
                ],
                "primary_defaults_signature": self.canonical_defaults.signature(),
                "secondary_signatures": [],
                "chosen_policy": (
                    "fallback_default_defaults" if defaults_guessed else "single_defaults_source"
                ),
                "allow_mixed_defaults": bool(getattr(ctx, "allow_mixed_defaults", False) if ctx else False),
                "allow_mixed_defaults_reason": (
                    getattr(ctx, "allow_mixed_defaults_reason", None) if ctx else None
                ),
            }
        
        print(f"  Found {len(self.registry)} unique atomtypes")
        print(f"  Detected {len(self.conflicts)} conflicts")
        if self.charge_warnings:
            print(f"  Charge warnings: {len(self.charge_warnings)}")
            print("  [WARN] Note: some forcefields store non-zero atomtype charges; "
                  "mismatches keep deterministic source values by policy.")
        
        # Validate conflicts (will raise on unresolved)
        self.validate_conflicts()
        
        # Validate charge strict mode
        self.validate_charge_strict()

        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        prefix_policy = (
            getattr(ctx, "prefix_implicit_topology_policy", "error") if ctx else "error"
        )
        prefix_reason = (
            getattr(ctx, "prefix_implicit_topology_reason", None) if ctx else None
        )
        self.prefix_policy_summary = self._prepare_prefix_policy(
            scan_paths=scan_paths_list,
            sanitize_paths=list(sanitize_sources.values()),
            output_dir=output_dir,
            run_id=run_id,
            policy=prefix_policy,
            reason=prefix_reason,
        )

        # Generate combined atomtypes
        combined_content = self.generate_combined_atomtypes()
        combined_path = output_dir / f"combined_atomtypes_{run_id}.itp"
        combined_path.write_text(combined_content)
        print(f"  Wrote: {combined_path.name}")
        
        # Create current symlink/copy
        combined_current = output_dir / "combined_atomtypes_current.itp"
        if combined_current.exists() or combined_current.is_symlink():
            combined_current.unlink()
        try:
            combined_current.symlink_to(combined_path.name)
        except OSError:
            # Symlinks may fail on Windows, use copy instead
            shutil.copy2(combined_path, combined_current)
        print(f"  Linked: {combined_current.name} -> {combined_path.name}")
        
        # Create sanitized directory
        sanitized_dir = output_dir / f"itp_sanitized_{run_id}"
        sanitized_dir.mkdir(parents=True, exist_ok=True)
        
        # Build sanitized output mapping (collision-safe)
        sanitized_files: List[Path] = []
        molecule_itp_files: List[Path] = []
        output_entries: List[Tuple[Path, Path]] = []  # (source_path, output_rel)
        output_map: Dict[str, Path] = {}  # casefold(rel) -> source_path
        sources_with_output: Set[Path] = set()

        def _safe_output_rel(include_target: str) -> Optional[Path]:
            rel = Path(include_target)
            if rel.is_absolute():
                return None
            if ".." in rel.parts:
                return None
            return rel

        def _register_output(source_path: Path, output_rel: Path) -> None:
            key = output_rel.as_posix().casefold()
            existing = output_map.get(key)
            source_id = _path_identity(source_path)
            if existing is not None:
                if _path_identity(existing) != source_id:
                    raise ItpSanitizerError(
                        "Sanitized output path collision (case-insensitive):\n"
                        f"  Output: {output_rel.as_posix()}\n"
                        f"  - {existing}\n"
                        f"  - {source_path}\n"
                        "Rename one of the inputs or disambiguate include paths."
                    )
            else:
                output_map[key] = source_path
                output_entries.append((source_path, output_rel))
            sources_with_output.add(source_id)

        # Register entry ITPs by basename (system.top includes these)
        entry_output_paths: Dict[Path, Path] = {}
        for itp_path in itp_paths:
            if not itp_path.exists():
                continue
            source_id = _path_identity(itp_path)
            output_rel = Path(itp_path.name)
            _register_output(itp_path, output_rel)
            entry_output_paths[source_id] = sanitized_dir / output_rel

        # Register include-target outputs for include-closure sanitization
        if sanitize_include_closure:
            for diag in include_diagnostics:
                if not diag.resolved_path:
                    continue
                resolved = Path(diag.resolved_path)
                source_id = _path_identity(resolved)
                if source_id not in sanitize_sources:
                    continue
                output_rel = _safe_output_rel(diag.include_target)
                if output_rel is None:
                    raise ItpSanitizerError(
                        "Include target cannot be mapped into sanitized output dir:\n"
                        f"  Include: {diag.include_target}\n"
                        f"  Resolved: {resolved}\n"
                        "Only relative include paths without '..' are supported for sanitization."
                    )
                _register_output(resolved, output_rel)

        # Ensure all sanitize sources have at least one output mapping
        missing_outputs = [
            str(source) for key, source in sanitize_sources.items()
            if key not in sources_with_output
        ]
        if missing_outputs:
            raise ItpSanitizerError(
                "Sanitization targets missing output mapping:\n  - "
                + "\n  - ".join(missing_outputs[:10])
                + ("\n  - ... (truncated)" if len(missing_outputs) > 10 else "")
            )

        # Sanitize and write files deterministically
        sanitized_cache: Dict[Path, str] = {}
        for source_path, output_rel in sorted(
            output_entries,
            key=lambda item: (item[1].as_posix().casefold(), item[1].as_posix()),
        ):
            if not source_path.exists():
                continue
            source_id = _path_identity(source_path)
            if source_id not in sanitized_cache:
                sanitized_cache[source_id] = self.sanitize_itp(source_path)
            out_path = sanitized_dir / output_rel
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(sanitized_cache[source_id])
            sanitized_files.append(out_path)

        # Moleculetype detection is include-aware
        for itp_path in itp_paths:
            if not itp_path.exists():
                continue
            source_id = _path_identity(itp_path)
            out_path = entry_output_paths.get(source_id)
            if out_path is None:
                continue
            if ItpParser.has_moleculetype_including_includes(
                itp_path,
                include_dirs=include_dirs,
                strict=strict_includes,
                check_shadowing=True,
                allow_shadowing=allow_shadowing,
            ):
                molecule_itp_files.append(out_path)
        
        print(f"  Wrote {len(sanitized_files)} sanitized ITPs to: {sanitized_dir.name}/")
        print(f"  Molecule ITPs (with [ moleculetype ]): {len(molecule_itp_files)}")
        
        # Create current directory symlink/copy
        sanitized_current = output_dir / "itp_sanitized_current"
        if sanitized_current.exists() or sanitized_current.is_symlink():
            if sanitized_current.is_symlink():
                sanitized_current.unlink()
            else:
                shutil.rmtree(sanitized_current)
        try:
            sanitized_current.symlink_to(sanitized_dir.name)
        except OSError:
            # Use copy for Windows
            shutil.copytree(sanitized_dir, sanitized_current)
        print(f"  Linked: {sanitized_current.name}/ -> {sanitized_dir.name}/")

        nonbond_params_cross_group_path = None
        nonbond_params_cross_group_current_path = None
        nonbond_params_cross_group_summary: Dict[str, Any] = {}
        if ctx is not None:
            (
                nonbond_params_cross_group_path,
                nonbond_params_cross_group_current_path,
                nonbond_params_cross_group_summary,
            ) = self._generate_cross_group_nonbond_params(
                output_dir=output_dir,
                run_id=run_id,
                source_defaults_map=source_defaults_map,
                defaults_report=defaults_report,
                scan_paths=scan_paths_list,
                ctx=ctx,
            )

        nonbond_params_secondary_path = None
        nonbond_params_secondary_current_path = None
        nonbond_params_secondary_summary: Dict[str, Any] = {}
        if ctx is not None:
            (
                nonbond_params_secondary_path,
                nonbond_params_secondary_current_path,
                nonbond_params_secondary_summary,
            ) = self._generate_secondary_nonbond_params(
                output_dir=output_dir,
                run_id=run_id,
                source_defaults_map=source_defaults_map,
                defaults_report=defaults_report,
                preserve_within_group=bool(
                    getattr(ctx, "mixed_defaults_preserve_within_group", False)
                ),
            )
        self.defaults_report = dict(defaults_report)
        
        return SanitizerResult(
            combined_atomtypes_path=combined_path,
            combined_atomtypes_current_path=combined_current,
            sanitized_dir=sanitized_dir,
            sanitized_current_dir=sanitized_current,
            sanitized_files=sanitized_files,
            molecule_itp_files=molecule_itp_files,
            conflicts_detected=self.conflicts,
            conflicts_overridden=self.overridden_conflicts,
            charge_warnings=self.charge_warnings,
            atomtype_count=len(self.registry),
            source_files_processed=len(self.source_files),
            defaults_used=self.canonical_defaults,
            source_defaults_map=source_defaults_map,
            mixed_defaults_detected=mixed_defaults_detected,
            mixed_defaults_report=defaults_report,
            include_resolution=include_diagnostics,
            defaults_guessed=defaults_guessed,
            nonbond_params_secondary_path=nonbond_params_secondary_path,
            nonbond_params_secondary_current_path=nonbond_params_secondary_current_path,
            nonbond_params_secondary_summary=nonbond_params_secondary_summary,
            nonbond_params_cross_group_path=nonbond_params_cross_group_path,
            nonbond_params_cross_group_current_path=nonbond_params_cross_group_current_path,
            nonbond_params_cross_group_summary=nonbond_params_cross_group_summary,
            prefix_policy_summary=self.prefix_policy_summary,
            prefix_injected_types_path=self.prefix_injected_types_path,
            prefix_injected_types_current_path=self.prefix_injected_types_current_path,
        )


# =============================================================================
# system.top Generation
# =============================================================================

def generate_system_top(
    combined_atomtypes_path: Path,
    sanitized_itp_paths: List[Path],
    system_name: str = "System",
    molecule_counts: Optional[Dict[str, int]] = None,
    defaults: Optional[DefaultsEntry] = None,
    molecule_itp_files: Optional[List[Path]] = None,
    top_dir: Optional[Path] = None,
    allow_guess_defaults: bool = False,
    validate_includes: bool = True,
    extra_includes: Optional[List[str]] = None,
) -> str:
    """
    Generate a system.top file with correct include order.
    
    Order:
    1. [ defaults ] section (mandatory for LJ interpretation)
    2. #include combined atomtypes (relative path to combined_atomtypes_path)
    3. Extra includes (POSRES, water models, etc.) - Task K
    4. #include for each molecule ITP (files with [ moleculetype ])
    5. [ system ] and [ molecules ]
    
    Args:
        combined_atomtypes_path: Path to combined atomtypes ITP (included via relative path)
        sanitized_itp_paths: List of all sanitized ITP paths
        system_name: System name for [ system ] section
        molecule_counts: Dict of molecule name to count for [ molecules ]
        defaults: DefaultsEntry to include (if None, fails unless allow_guess_defaults)
        molecule_itp_files: Subset of sanitized_itp_paths that contain [ moleculetype ]
        top_dir: Directory where system.top will be written (REQUIRED for correct relative paths)
        allow_guess_defaults: If True, use standard defaults when none provided (silent physics risk!)
        validate_includes: If True, verify all include targets exist on disk
        extra_includes: Task K - Optional list of extra include paths (e.g., "posre_polymer.itp", "tip3p.itp")
                       These are inserted after atomtypes but before molecule definitions.
    
    Returns:
        system.top content
        
    Raises:
        ItpSanitizerError: If defaults missing and allow_guess_defaults=False
        ItpSanitizerError: If top_dir is None
        ItpSanitizerError: If validate_includes=True and an include target is missing
    """
    # Fail-fast: require top_dir
    if top_dir is None:
        raise ItpSanitizerError(
            "generate_system_top requires top_dir to compute correct relative include paths. "
            "Pass the directory where system.top will be written."
        )
    
    molecule_counts = molecule_counts or {}
    
    lines = [
        "; System topology generated by ITP Sanitizer",
        "; IMPORTANT: [ defaults ] must appear before any LJ parameter interpretation",
        "",
    ]
    
    # Include defaults first - fail-fast if missing
    if defaults:
        lines.extend([
            "[ defaults ]",
            "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
            f"  {defaults.nbfunc}        {defaults.comb_rule}         {defaults.gen_pairs}       {defaults.fudge_lj}    {defaults.fudge_qq}",
            "",
        ])
    else:
        if not allow_guess_defaults:
            raise ItpSanitizerError(
                "No [ defaults ] section found in inputs. Refusing to invent physics.\n"
                "Either provide explicit defaults in input ITPs or enable --allow-default-defaults."
            )
        # Use standard defaults with warning logged
        print("  [WARN] No [ defaults ] found; using fallback (nbfunc=1, comb-rule=2)")
        defaults = _fallback_defaults_entry()
        lines.extend([
            "[ defaults ]",
            "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
            f"  {defaults.nbfunc}        {defaults.comb_rule}         {defaults.gen_pairs}       "
            f"{defaults.fudge_lj}    {defaults.fudge_qq}  ; FALLBACK: LJ, geometric, standard fudge",
            "",
        ])
    
    combined_rel = Path(os.path.relpath(combined_atomtypes_path, top_dir))

    # Track paths for validation
    include_paths = [(combined_atomtypes_path, combined_rel.as_posix())]
    
    lines.extend([
        "; Combined atomtypes (after defaults for correct LJ interpretation)",
        f'#include "{combined_rel.as_posix()}"',
        "",
    ])
    
    # Task K: Insert extra includes (POSRES, water models, etc.)
    if extra_includes:
        lines.append("; Extra includes (POSRES, water models, etc.)")
        for inc in extra_includes:
            lines.append(f'#include "{inc}"')
        lines.append("")
    
    lines.append("; Molecule definitions (only files with [ moleculetype ] section)")
    
    # Only include molecule ITPs if list is provided
    if molecule_itp_files:
        for itp_path in molecule_itp_files:
            rel_path = Path(os.path.relpath(itp_path, top_dir))
            lines.append(f'#include "{rel_path.as_posix()}"')
            include_paths.append((itp_path, rel_path.as_posix()))
    else:
        # Fallback: include all sanitized ITPs (legacy behavior)
        for itp_path in sanitized_itp_paths:
            rel_path = Path(os.path.relpath(itp_path, top_dir))
            lines.append(f'#include "{rel_path.as_posix()}"')
            include_paths.append((itp_path, rel_path.as_posix()))
    
    lines.extend([
        "",
        "[ system ]",
        f"{system_name}",
        "",
        "[ molecules ]",
        "; name    count",
    ])
    
    for mol_name, count in molecule_counts.items():
        lines.append(f"{mol_name:10s} {count}")
    
    content = '\n'.join(lines) + '\n'
    
    # Validate include paths exist on disk
    if validate_includes:
        for abs_path, rel_str in include_paths:
            if not abs_path.exists():
                raise ItpSanitizerError(
                    f"Include target does not exist: {rel_str}\n"
                    f"  Expected at: {abs_path}\n"
                    "Ensure all ITP files are created before generating system.top."
                )
    
    return content
