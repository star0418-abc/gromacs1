"""
Topology sanitizer helpers extracted from sanitizer.py.

Contains ITP/include/topology/forcefield sanitization and final output helpers.
"""

from dataclasses import dataclass, field
from collections import deque
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
from pathlib import Path, PurePosixPath
import sys
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Set, Tuple
import math
import hashlib
import os
import re
import shutil
import uuid

from ..manifest import _best_effort_fsync, _best_effort_fsync_dir
from ..itp_sanitizer import IncludeResolution, ItpParser

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
DEFAULT_MOLECULE_SIG_MASS_QUANT_STEP = 1e-6
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

LJ_ROBUST_LOG_MEDIAN_RELATIVE_FLOOR = 1e-12
LJ_ROBUST_LOG_ABSOLUTE_FLOOR = sys.float_info.min * 1024.0
LJ_ROBUST_Z_THRESHOLD = 5.0

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

PREPROCESSOR_CONDITIONAL_START_DIRECTIVES = frozenset({"#if", "#ifdef", "#ifndef"})
PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES = frozenset({"#elif", "#else"})
PREPROCESSOR_CONDITIONAL_END_DIRECTIVE = "#endif"
PREPROCESSOR_MACRO_DIRECTIVES = frozenset({"#define", "#undef"})
PREPROCESSOR_UNRESOLVED_SEMANTIC_DIRECTIVES = (
    PREPROCESSOR_CONDITIONAL_START_DIRECTIVES
    | PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES
    | PREPROCESSOR_MACRO_DIRECTIVES
    | frozenset({PREPROCESSOR_CONDITIONAL_END_DIRECTIVE})
)

# Supported ptype tokens in [ atomtypes ].
PTYPE_TOKENS = frozenset({"A", "D", "S", "V"})

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


def preprocessor_directive_token(line: str) -> Optional[str]:
    """Return normalized preprocessor directive token for a raw line."""
    stripped = line.lstrip()
    if not stripped.startswith("#"):
        return None
    parts = stripped.split(None, 1)
    if not parts:
        return None
    return parts[0].lower()


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


def _median(values: List[float]) -> float:
    """Deterministic median helper for small numeric lists."""
    if not values:
        raise ValueError("median requires at least one value")
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def _matches_protected_pattern(name: str, patterns: frozenset) -> bool:
    """
    Check if name matches any protected pattern using canonicalized boundaries.
    """
    import re as _re

    def _canonical_parts(value: str) -> Tuple[str, List[str], str]:
        raw_upper = value.upper().strip()
        tokens = [token for token in _re.split(r"[^A-Z0-9]+", raw_upper) if token]
        normalized = "".join(tokens)
        return raw_upper, tokens, normalized

    raw_upper, _, normalized = _canonical_parts(name)
    for pattern in patterns:
        _, pattern_tokens, pattern_normalized = _canonical_parts(pattern)
        if not pattern_normalized:
            continue

        if len(pattern_normalized) <= 3:
            if normalized == pattern_normalized:
                return True
            continue

        if normalized == pattern_normalized:
            return True

        token_pattern = (
            r"[^A-Z0-9]*".join(_re.escape(token) for token in pattern_tokens)
            if pattern_tokens
            else _re.escape(pattern_normalized)
        )
        boundary_re = rf"(?<![A-Z0-9]){token_pattern}(?![A-Z0-9])"
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


@dataclass(frozen=True)
class TopologyTruthSource:
    """Selected topology truth source for semantic reads."""
    path: Path
    source: str
    is_preprocessed: bool = False
    fallback_used: bool = False
    reason: Optional[str] = None
    candidate_paths: Tuple[str, ...] = ()


ATOMS_ROW_STATUS_PARSED = "parsed"
ATOMS_ROW_STATUS_UNSUPPORTED = "unsupported"
ATOMS_ROW_STATUS_INVALID = "invalid"

TOP_MOLECULES_ENTRY_STATUS_PARSED = "parsed"
TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED = "unresolved"
TOP_MOLECULES_ENTRY_STATUS_INVALID = "invalid"

MOLECULETYPE_IR_STATUS_EXACT = "exact"
MOLECULETYPE_IR_STATUS_UNCERTAIN = "uncertain"
MOLECULETYPE_IR_STATUS_MISSING = "missing"

MOLECULETYPE_COMPARE_EQUAL = "equal"
MOLECULETYPE_COMPARE_DIFFERENT = "different"
MOLECULETYPE_COMPARE_UNCERTAIN = "uncertain"


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
class LocalAtomTypeRecord:
    """Minimal local LJ record for atomtypes used by a target moleculetype."""
    name: str
    ptype: Optional[str]
    sigma: float
    epsilon: float


@dataclass(frozen=True)
class TopDefaultsOwnershipScan:
    """Existing system.top defaults/forcefield ownership scan result."""
    has_defaults: bool
    has_forcefield_include: bool
    owner_insert_after_idx: Optional[int] = None
    ambiguity_reasons: Tuple[str, ...] = ()


@dataclass
class MoleculeTypeIRParseResult:
    """Structured target-moleculetype parse result for same-name checks."""
    status: str
    molecule: Optional[MoleculeTypeIR] = None
    actual_name: Optional[str] = None
    signature: Optional[str] = None
    total_charge: Optional[float] = None
    uncertainty_reasons: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class MoleculeTypeComparisonResult:
    """Explicit same-name comparison result."""
    status: str
    reason: Optional[str] = None


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

    @staticmethod
    def _project_root_path(ctx: Optional["PipelineContext"]) -> Optional[Path]:
        """Return resolved project root if configured."""
        project_root = getattr(ctx, "project_root", None) if ctx is not None else None
        if not project_root:
            return None
        try:
            return Path(project_root).resolve()
        except OSError:
            return Path(project_root)

    def _resolve_include_root_path(
        self,
        raw_path: Path,
        ctx: Optional["PipelineContext"],
    ) -> Path:
        """Resolve include-root candidates against project_root when needed."""
        project_root_path = self._project_root_path(ctx)
        if raw_path.is_absolute() or project_root_path is None:
            return raw_path
        return project_root_path / raw_path

    def _collect_ctx_include_roots(self, ctx: Optional["PipelineContext"]) -> Tuple[List[Path], List[Path]]:
        """Return configured include roots split into user paths and GROMACS share/data paths."""
        if ctx is None:
            return [], []

        user_paths: List[Path] = []
        raw_user_paths = getattr(ctx, "itp_include_dirs", "")
        if raw_user_paths:
            for item in str(raw_user_paths).split(","):
                item = item.strip()
                if not item:
                    continue
                user_paths.append(self._resolve_include_root_path(Path(item), ctx))

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
                raw_items = [str(entry) for entry in value if str(entry).strip()]
            else:
                raw_items = [entry.strip() for entry in str(value).split(",") if entry.strip()]
            for item in raw_items:
                gmx_share_paths.append(self._resolve_include_root_path(Path(item), ctx))
        return user_paths, gmx_share_paths

    def _dedupe_existing_paths(self, paths: List[Path]) -> List[Path]:
        """Keep existing paths only, preserving deterministic order."""
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

    @staticmethod
    def _summarize_semantic_scan_issues(issues: List[str], max_items: int = 3) -> str:
        """Compact human-readable summary for degraded semantic scans."""
        if not issues:
            return ""
        summary = "; ".join(issues[:max_items])
        remaining = len(issues) - min(len(issues), max_items)
        if remaining > 0:
            summary += f"; ... {remaining} more"
        return summary

    def _record_semantic_scan_issue(
        self,
        *,
        issues: List[str],
        file_path: Path,
        line_no: int,
        directive: str,
        purpose: str,
        strict: bool,
    ) -> None:
        """Record or raise on unresolved preprocessor ambiguity during Python parsing."""
        reason = (
            f"{directive} at {file_path.name}:{line_no} makes {purpose} ambiguous "
            "without preprocessing"
        )
        if reason not in issues:
            issues.append(reason)
        if strict:
            raise SanitizerError(
                f"Cannot safely {purpose} in {file_path}: unresolved preprocessor directive "
                f"{directive} at line {line_no}. Use a preprocessed topology truth source or "
                "disable strict include resolution."
            )

    def _emit_semantic_scan_warning(
        self,
        *,
        file_path: Path,
        purpose: str,
        issues: List[str],
    ) -> None:
        """Emit one explicit degraded warning for conservative Python fallback parsing."""
        if not issues:
            return
        summary = self._summarize_semantic_scan_issues(issues)
        print(
            "  [WARN][DEGRADED] "
            f"Conservative Python fallback used for {purpose} in {file_path}: {summary}"
        )

    def _conservative_topology_lines(
        self,
        *,
        lines: List[str],
        file_path: Path,
        purpose: str,
        strict: bool,
        ignore_ranges: Optional[List[Tuple[int, int]]] = None,
    ) -> Tuple[List[Tuple[int, str]], List[str]]:
        """
        Return only topology lines that are unambiguous without preprocessing.

        Preprocessor directives make Python parsing ambiguous. Strict mode fails
        closed via `_record_semantic_scan_issue()`. Non-strict mode records one
        degraded warning and skips conditionally ambiguous content.
        """
        skip_indices: Set[int] = set()
        if ignore_ranges:
            for start, end in ignore_ranges:
                skip_indices.update(range(start, end + 1))

        active_lines: List[Tuple[int, str]] = []
        semantic_issues: List[str] = []
        conditional_depth = 0

        for idx, line in enumerate(lines):
            if idx in skip_indices:
                continue

            line_no = idx + 1
            directive = preprocessor_directive_token(line)

            if directive == "#include":
                if conditional_depth > 0:
                    self._record_semantic_scan_issue(
                        issues=semantic_issues,
                        file_path=file_path,
                        line_no=line_no,
                        directive=directive,
                        purpose=purpose,
                        strict=strict,
                    )
                    continue
                include_target = parse_include_target(line)
                if include_target is None:
                    self._record_semantic_scan_issue(
                        issues=semantic_issues,
                        file_path=file_path,
                        line_no=line_no,
                        directive="#include(macro-like)",
                        purpose=purpose,
                        strict=strict,
                    )
                    continue
                active_lines.append((line_no, line))
                continue

            if directive is not None:
                self._record_semantic_scan_issue(
                    issues=semantic_issues,
                    file_path=file_path,
                    line_no=line_no,
                    directive=directive,
                    purpose=purpose,
                    strict=strict,
                )
                if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                    conditional_depth += 1
                elif directive in PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES:
                    if conditional_depth == 0:
                        conditional_depth = 1
                elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE:
                    if conditional_depth > 0:
                        conditional_depth -= 1
                continue

            if conditional_depth > 0:
                continue

            active_lines.append((line_no, line))

        if conditional_depth > 0:
            self._record_semantic_scan_issue(
                issues=semantic_issues,
                file_path=file_path,
                line_no=len(lines) or 1,
                directive=f"unclosed conditional depth={conditional_depth}",
                purpose=purpose,
                strict=strict,
            )
        if semantic_issues and not strict:
            self._emit_semantic_scan_warning(
                file_path=file_path,
                purpose=purpose,
                issues=semantic_issues,
            )
        return active_lines, semantic_issues

    def _candidate_preprocessed_topology_paths(
        self,
        *,
        output_dir: Optional[Path] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> List[Path]:
        """Return existing candidate preprocessed topology paths in priority order."""
        candidates: List[Path] = []
        if ctx is not None:
            for attr in (
                "preprocessed_topology_path",
                "preprocessed_topology",
                "preprocessed_top_path",
                "topology_truth_source_path",
            ):
                value = getattr(ctx, attr, None)
                if not value:
                    continue
                candidate = self._resolve_include_root_path(Path(str(value)), ctx)
                candidates.append(candidate)
        if output_dir is not None and hasattr(self, "_grompp_preprocessed_top_path"):
            try:
                candidates.append(self._grompp_preprocessed_top_path(output_dir))
            except Exception:
                pass
        return self._dedupe_existing_paths(candidates)

    def _select_topology_truth_source(
        self,
        top_path: Path,
        *,
        output_dir: Optional[Path] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> TopologyTruthSource:
        """Prefer preprocessed topology when available; otherwise fall back explicitly."""
        candidate_paths = self._candidate_preprocessed_topology_paths(
            output_dir=output_dir,
            ctx=ctx,
        )
        if candidate_paths:
            selected = candidate_paths[0]
            return TopologyTruthSource(
                path=selected,
                source="preprocessed_topology",
                is_preprocessed=True,
                fallback_used=False,
                candidate_paths=tuple(str(path) for path in candidate_paths),
            )
        return TopologyTruthSource(
            path=top_path,
            source="python_fallback_raw_topology",
            is_preprocessed=False,
            fallback_used=True,
            reason="no preprocessed topology truth source available",
            candidate_paths=(),
        )

    def _derive_include_search_paths(
        self,
        current_file: Path,
        include_dirs: Optional[List[Path]] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> List[Path]:
        """
        Compatibility wrapper around the single include-priority implementation.
        """
        return list(
            self._build_include_search_paths(
                ctx,
                current_file,
                extra_include_dirs=include_dirs,
            )
        )

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

        names: List[str] = []
        uncertain = False
        current_section = ""
        strict_includes = bool(ctx and ctx.strict_include_resolution)
        allow_shadowing = bool(ctx and getattr(ctx, "allow_include_shadowing", False))
        content = self.safe_read_text(itp_path, strict=strict_includes)
        if not content and not strict_includes:
            uncertain = True
        semantic_issues: List[str] = []
        conditional_depth = 0

        for line_no, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            directive = preprocessor_directive_token(line)

            if directive == "#include":
                if conditional_depth > 0:
                    self._record_semantic_scan_issue(
                        issues=semantic_issues,
                        file_path=itp_path,
                        line_no=line_no,
                        directive=directive,
                        purpose="moleculetype name extraction",
                        strict=strict_includes,
                    )
                    uncertain = True
                    continue
                include_target = parse_include_target(line)
                if include_target is None:
                    self._record_semantic_scan_issue(
                        issues=semantic_issues,
                        file_path=itp_path,
                        line_no=line_no,
                        directive="#include(macro-like)",
                        purpose="moleculetype name extraction",
                        strict=strict_includes,
                    )
                    uncertain = True
                    continue
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
                        include_path,
                        visited,
                        include_dirs=include_dirs,
                        ctx=ctx,
                    )
                    names.extend(inc_names)
                    uncertain = uncertain or inc_uncertain
                else:
                    uncertain = True
                continue

            if directive is not None:
                self._record_semantic_scan_issue(
                    issues=semantic_issues,
                    file_path=itp_path,
                    line_no=line_no,
                    directive=directive,
                    purpose="moleculetype name extraction",
                    strict=strict_includes,
                )
                uncertain = True
                if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                    conditional_depth += 1
                elif directive in PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES:
                    if conditional_depth == 0:
                        conditional_depth = 1
                elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE:
                    if conditional_depth > 0:
                        conditional_depth -= 1
                continue
            if conditional_depth > 0:
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

        if conditional_depth > 0:
            self._record_semantic_scan_issue(
                issues=semantic_issues,
                file_path=itp_path,
                line_no=len(content.splitlines()) or 1,
                directive=f"unclosed conditional depth={conditional_depth}",
                purpose="moleculetype name extraction",
                strict=strict_includes,
            )
            uncertain = True
        if semantic_issues and not strict_includes:
            self._emit_semantic_scan_warning(
                file_path=itp_path,
                purpose="moleculetype name extraction",
                issues=semantic_issues,
            )
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
        Resolve an include target using the single per-file search-order builder.
        """
        search_order = self._derive_include_search_paths(
            including_file,
            include_dirs=include_dirs,
            ctx=ctx,
        )

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

            # Scan for includes
            content = self.safe_read_text(itp, strict=strict_includes)
            if not content:
                continue

            semantic_issues: List[str] = []
            conditional_depth = 0
            for line_no, line in enumerate(content.splitlines(), start=1):
                directive = preprocessor_directive_token(line)
                if directive == "#include":
                    if conditional_depth > 0:
                        self._record_semantic_scan_issue(
                            issues=semantic_issues,
                            file_path=itp,
                            line_no=line_no,
                            directive=directive,
                            purpose="include closure discovery",
                            strict=strict_includes,
                        )
                        continue
                    include_target = parse_include_target(line)
                    if not include_target:
                        self._record_semantic_scan_issue(
                            issues=semantic_issues,
                            file_path=itp,
                            line_no=line_no,
                            directive="#include(macro-like)",
                            purpose="include closure discovery",
                            strict=strict_includes,
                        )
                        continue

                    inc_path = self._resolve_include_in_dirs(
                        include_target,
                        itp,
                        [],
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

                    if str(inc_resolved) not in ff_map:
                        class_map[str(inc_resolved)] = parent_class
                        ff_map[str(inc_resolved)] = parent_ff
                        print(f"    + Include ({parent_class}): {inc_resolved.name}")

                    queue.append(inc_resolved)
                    continue

                if directive is not None:
                    self._record_semantic_scan_issue(
                        issues=semantic_issues,
                        file_path=itp,
                        line_no=line_no,
                        directive=directive,
                        purpose="include closure discovery",
                        strict=strict_includes,
                    )
                    if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                        conditional_depth += 1
                    elif directive in PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES:
                        if conditional_depth == 0:
                            conditional_depth = 1
                    elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE:
                        if conditional_depth > 0:
                            conditional_depth -= 1
                    continue
                if conditional_depth > 0:
                    continue

            if conditional_depth > 0:
                self._record_semantic_scan_issue(
                    issues=semantic_issues,
                    file_path=itp,
                    line_no=len(content.splitlines()) or 1,
                    directive=f"unclosed conditional depth={conditional_depth}",
                    purpose="include closure discovery",
                    strict=strict_includes,
                )
            if semantic_issues and not strict_includes:
                self._emit_semantic_scan_warning(
                    file_path=itp,
                    purpose="include closure discovery",
                    issues=semantic_issues,
                )
        
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
        """
        strict_compare = self._strict_moleculetype_conflict_mode(ctx)
        comparison = self._compare_moleculetype_ir(
            itp1,
            itp2,
            mol_name,
            ctx=ctx,
        )
        if comparison.status == MOLECULETYPE_COMPARE_UNCERTAIN:
            msg = comparison.reason or (
                f"Uncertain same-name moleculetype comparison for '{mol_name}'."
            )
            if strict_compare:
                raise SanitizerError(msg)
            print(f"  [WARN][DEGRADED] {msg}. Treating definitions as conflicting.")
            return True
        return comparison.status == MOLECULETYPE_COMPARE_DIFFERENT

    @staticmethod
    def _quantize_moleculetype_value(raw_value: float, quant_step: float) -> float:
        """
        Quantize a scalar deterministically for signature hashing.

        Decimal arithmetic avoids platform-dependent binary midpoint/tie behavior.
        """
        if not math.isfinite(raw_value) or not math.isfinite(quant_step) or quant_step <= 0:
            return raw_value
        try:
            step_dec = Decimal(str(quant_step))
            raw_dec = Decimal(str(raw_value))
            scaled = raw_dec / step_dec
            quantized = scaled.to_integral_value(rounding=ROUND_HALF_EVEN) * step_dec
            value = float(quantized)
        except (InvalidOperation, OverflowError, ValueError, ZeroDivisionError):
            return raw_value
        return 0.0 if value == 0.0 else value

    @staticmethod
    def _quantize_moleculetype_charge(raw_q: float, quant_step: float) -> float:
        """Compatibility wrapper for charge quantization."""
        return TopologySanitizerMixin._quantize_moleculetype_value(raw_q, quant_step)

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

    def _resolve_requested_moleculetype_name_in_content(
        self,
        content: str,
        requested_name: str,
        itp_path: Path,
    ) -> Tuple[Optional[str], List[str], Optional[str]]:
        """Resolve request to one real [ moleculetype ] name in content."""
        available: List[str] = []
        current_section = ""
        awaiting_moleculetype = False

        for raw in content.splitlines():
            section = parse_section_name(raw)
            if section is not None:
                current_section = section
                awaiting_moleculetype = section == "moleculetype"
                continue

            if current_section != "moleculetype" or not awaiting_moleculetype:
                continue

            data = strip_gmx_comment(raw.strip()).strip()
            if not data:
                continue
            parts = data.split()
            if not parts:
                continue
            available.append(parts[0])
            awaiting_moleculetype = False

        if not available:
            fallback = requested_name or itp_path.stem
            return fallback, [], None
        if requested_name in available:
            return requested_name, available, None

        requested_key = requested_name.casefold()
        casefold_matches = [
            name for name in available if name.casefold() == requested_key
        ]
        if len(casefold_matches) == 1:
            return casefold_matches[0], available, None
        if len(casefold_matches) > 1:
            return (
                None,
                available,
                (
                    f"ambiguous case-insensitive moleculetype match for '{requested_name}' "
                    f"in {itp_path}: {casefold_matches}"
                ),
            )
        if len(available) == 1:
            return available[0], available, None
        return (
            None,
            available,
            (
                f"moleculetype '{requested_name}' not found in {itp_path}; "
                f"available moleculetypes: {available}"
            ),
        )

    def _build_moleculetype_ir_signature(
        self,
        molecule: MoleculeTypeIR,
        topology_entries: List[Tuple[str, str]],
        local_atomtypes: Dict[str, LocalAtomTypeRecord],
    ) -> str:
        """Build deterministic same-name signature from typed IR plus local topology."""
        hasher = hashlib.sha256()
        quant_step = self._molecule_sig_quant_step()
        mass_quant_step = DEFAULT_MOLECULE_SIG_MASS_QUANT_STEP
        section_counts: Dict[str, int] = {}
        sections_present: Set[str] = set()
        qmin: Optional[float] = None
        qmax: Optional[float] = None
        sum_abs_q = 0.0
        sum_q2 = 0.0

        for atom in molecule.atoms:
            quant_q = self._quantize_moleculetype_charge(atom.charge, quant_step)
            qmin = quant_q if qmin is None else min(qmin, quant_q)
            qmax = quant_q if qmax is None else max(qmax, quant_q)
            sum_abs_q += abs(quant_q)
            sum_q2 += quant_q * quant_q
            quant_mass = self._quantize_moleculetype_value(atom.mass, mass_quant_step)
            cgnr_token = str(atom.cgnr) if atom.cgnr is not None else "-"
            canonical = (
                "atoms|"
                f"{atom.atom_id}|{atom.atom_type}|{atom.resnr}|{atom.residue}|"
                f"{atom.atom_name}|{cgnr_token}|{quant_q:+.10f}|{quant_mass:.6f}\n"
            )
            hasher.update(canonical.encode("utf-8"))

        section_counts["atoms"] = len(molecule.atoms)
        sections_present.add("atoms")

        for section, canonical_tokens in topology_entries:
            sections_present.add(section)
            section_counts[section] = section_counts.get(section, 0) + 1
            hasher.update(f"{section}|{canonical_tokens}\n".encode("utf-8"))

        if local_atomtypes:
            sections_present.add("atomtypes_local")
            section_counts["atomtypes_local"] = len(local_atomtypes)
            for atomtype_name in sorted(local_atomtypes, key=str.casefold):
                record = local_atomtypes[atomtype_name]
                sigma = self._quantize_moleculetype_value(record.sigma, quant_step)
                epsilon = self._quantize_moleculetype_value(record.epsilon, quant_step)
                ptype_token = record.ptype or "-"
                hasher.update(
                    (
                        "atomtypes_local|"
                        f"{record.name}|{ptype_token}|{sigma:.10f}|{epsilon:.10f}\n"
                    ).encode("utf-8")
                )

        qmin_val = qmin if qmin is not None else 0.0
        qmax_val = qmax if qmax is not None else 0.0
        return (
            f"atoms={len(molecule.atoms)};"
            f"digest={hasher.hexdigest()};"
            f"qmin={qmin_val:+.6f};"
            f"qmax={qmax_val:+.6f};"
            f"sum_abs_q={sum_abs_q:.6f};"
            f"sum_q2={sum_q2:.6f};"
            f"charge_quant_step={quant_step:.12g};"
            f"mass_quant_step={mass_quant_step:.12g};"
            f"sections={','.join(sorted(sections_present))};"
            f"section_counts={','.join(f'{k}:{section_counts[k]}' for k in sorted(section_counts))}"
        )

    def _describe_moleculetype_parse_result(
        self,
        result: MoleculeTypeIRParseResult,
        itp_path: Path,
        requested_name: str,
    ) -> str:
        """Format a concise parse status diagnostic for same-name comparisons."""
        if result.status == MOLECULETYPE_IR_STATUS_MISSING:
            reason = "; ".join(result.uncertainty_reasons) or "target moleculetype not found"
            return f"{itp_path}: missing '{requested_name}' ({reason})"
        if result.status == MOLECULETYPE_IR_STATUS_UNCERTAIN:
            reason = "; ".join(result.uncertainty_reasons) or "unsupported or ambiguous target parse"
            actual = result.actual_name or requested_name
            return f"{itp_path}: uncertain parse for '{actual}' ({reason})"
        actual = result.actual_name or requested_name
        atom_count = len(result.molecule.atoms) if result.molecule is not None else 0
        return f"{itp_path}: exact parse for '{actual}' ({atom_count} atoms)"

    def _parse_moleculetype_ir(
        self,
        itp_path: Path,
        mol_name: str,
        ctx: Optional["PipelineContext"] = None,
    ) -> MoleculeTypeIRParseResult:
        """Parse one target moleculetype into the Stage 1 IR with explicit status."""
        if not itp_path.exists():
            return MoleculeTypeIRParseResult(
                status=MOLECULETYPE_IR_STATUS_MISSING,
                uncertainty_reasons=[f"{itp_path} does not exist"],
            )
        strict_reads = bool(ctx and ctx.strict_include_resolution)
        content = self.safe_read_text(itp_path, strict=strict_reads)
        if not content:
            return MoleculeTypeIRParseResult(
                status=MOLECULETYPE_IR_STATUS_MISSING,
                uncertainty_reasons=[f"{itp_path} is empty or unreadable"],
            )

        resolved_name, available_names, resolution_error = (
            self._resolve_requested_moleculetype_name_in_content(
                content,
                mol_name,
                itp_path,
            )
        )
        if resolved_name is None:
            return MoleculeTypeIRParseResult(
                status=MOLECULETYPE_IR_STATUS_MISSING,
                uncertainty_reasons=[resolution_error or f"unable to resolve '{mol_name}'"],
            )

        current_section = ""
        current_mol: Optional[str] = None
        in_target_mol = False
        target_found = False
        awaiting_moleculetype = False
        atoms: List[AtomRecord] = []
        uncertainty_reasons: List[str] = []
        topology_entries: List[Tuple[str, str]] = []
        local_atomtypes: Dict[str, LocalAtomTypeRecord] = {}
        local_atomtype_issues: Dict[str, str] = {}
        used_atomtypes: Set[str] = set()
        topology_sections = {
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

        for line_no, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                awaiting_moleculetype = section == "moleculetype"
                if section == "moleculetype":
                    in_target_mol = False
                continue

            data = strip_gmx_comment(stripped).strip()
            if not data:
                continue

            if current_section == "atomtypes" and not data.startswith("#"):
                atomtype_name = data.split()[0]
                try:
                    parsed_name, ptype, sigma, epsilon, _ = self._parse_lj_atomtypes_row(
                        data,
                        itp_path.name,
                        line_no,
                    )
                except ValueError as exc:
                    local_atomtype_issues[atomtype_name] = (
                        f"malformed local [ atomtypes ] row for '{atomtype_name}' at line {line_no}: {exc}"
                    )
                    continue
                record = LocalAtomTypeRecord(
                    name=parsed_name,
                    ptype=ptype,
                    sigma=sigma,
                    epsilon=epsilon,
                )
                previous = local_atomtypes.get(parsed_name)
                if previous is not None and previous != record:
                    local_atomtype_issues[parsed_name] = (
                        f"conflicting local [ atomtypes ] definitions for '{parsed_name}'"
                    )
                    continue
                local_atomtypes[parsed_name] = record
                continue

            if current_section == "moleculetype":
                if data.startswith("#"):
                    continue
                if not awaiting_moleculetype:
                    continue
                parts = data.split()
                if not parts:
                    continue
                current_mol = parts[0]
                awaiting_moleculetype = False
                if current_mol == resolved_name:
                    if target_found:
                        uncertainty_reasons.append(
                            f"multiple [ moleculetype ] blocks for '{resolved_name}' in {itp_path}"
                        )
                        in_target_mol = False
                        continue
                    in_target_mol = True
                    target_found = True
                else:
                    in_target_mol = False
                continue

            if not in_target_mol:
                continue

            if data.startswith("#"):
                uncertainty_reasons.append(
                    f"preprocessor directive in target moleculetype '{resolved_name}' at line {line_no}"
                )
                continue

            if current_section == "atoms":
                row_result = parse_atoms_row(data)
                if row_result.status != ATOMS_ROW_STATUS_PARSED or row_result.record is None:
                    reason = row_result.reason or "unsupported [ atoms ] row"
                    uncertainty_reasons.append(
                        f"[ atoms ] line {line_no}: {reason}"
                    )
                    continue
                atoms.append(row_result.record)
                used_atomtypes.add(row_result.record.atom_type)
                continue

            if current_section in topology_sections:
                topology_entries.append((current_section, " ".join(data.split())))
                continue

        if not target_found:
            return MoleculeTypeIRParseResult(
                status=MOLECULETYPE_IR_STATUS_MISSING,
                uncertainty_reasons=[
                    (
                        f"target moleculetype '{mol_name}' not found in {itp_path}; "
                        f"available moleculetypes: {available_names or ['<none>']}"
                    )
                ],
            )

        for atomtype_name in sorted(used_atomtypes, key=str.casefold):
            if atomtype_name in local_atomtype_issues:
                uncertainty_reasons.append(local_atomtype_issues[atomtype_name])

        molecule = MoleculeTypeIR(name=resolved_name, atoms=atoms)
        total_charge = math.fsum(atom.charge for atom in atoms) if atoms else None
        if not atoms:
            uncertainty_reasons.append(
                f"target moleculetype '{resolved_name}' has no parseable [ atoms ] rows"
            )

        status = (
            MOLECULETYPE_IR_STATUS_EXACT
            if atoms and not uncertainty_reasons
            else MOLECULETYPE_IR_STATUS_UNCERTAIN
        )
        signature = None
        if status == MOLECULETYPE_IR_STATUS_EXACT:
            signature = self._build_moleculetype_ir_signature(
                molecule,
                topology_entries,
                {
                    name: record
                    for name, record in local_atomtypes.items()
                    if name in used_atomtypes
                },
            )

        return MoleculeTypeIRParseResult(
            status=status,
            molecule=molecule if atoms else None,
            actual_name=resolved_name,
            signature=signature,
            total_charge=total_charge,
            uncertainty_reasons=uncertainty_reasons,
        )

    def _compare_moleculetype_ir(
        self,
        itp1: Path,
        itp2: Path,
        mol_name: str,
        ctx: Optional["PipelineContext"] = None,
    ) -> MoleculeTypeComparisonResult:
        """Compare same-name moleculetype definitions using Stage 1 IR results."""
        parsed1 = self._parse_moleculetype_ir(itp1, mol_name, ctx=ctx)
        parsed2 = self._parse_moleculetype_ir(itp2, mol_name, ctx=ctx)

        if (
            parsed1.status != MOLECULETYPE_IR_STATUS_EXACT
            or parsed2.status != MOLECULETYPE_IR_STATUS_EXACT
        ):
            reason = (
                f"Uncertain same-name moleculetype comparison for '{mol_name}': "
                f"{self._describe_moleculetype_parse_result(parsed1, itp1, mol_name)}; "
                f"{self._describe_moleculetype_parse_result(parsed2, itp2, mol_name)}"
            )
            return MoleculeTypeComparisonResult(
                status=MOLECULETYPE_COMPARE_UNCERTAIN,
                reason=reason,
            )

        if parsed1.signature != parsed2.signature:
            return MoleculeTypeComparisonResult(
                status=MOLECULETYPE_COMPARE_DIFFERENT,
                reason=(
                    f"same-name moleculetype '{mol_name}' differs: deterministic IR signatures mismatch "
                    f"between {itp1} and {itp2}"
                ),
            )

        return MoleculeTypeComparisonResult(status=MOLECULETYPE_COMPARE_EQUAL)

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
        parsed = self._parse_moleculetype_ir(
            itp_path,
            mol_name,
            ctx=ctx,
        )
        if parsed.molecule is None or not parsed.molecule.atoms:
            return None, parsed.status == MOLECULETYPE_IR_STATUS_UNCERTAIN
        return (
            [atom.charge for atom in parsed.molecule.atoms],
            parsed.status == MOLECULETYPE_IR_STATUS_UNCERTAIN,
        )
    
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
            itp_moleculetype_name: Optional[str] = None,
            parse_status: Optional[str] = None,
            parse_reasons: Optional[List[str]] = None,
        ) -> None:
            detail = {
                "classification": classification,
                "status": status,
                "source": source,
                "reason": reason,
                "itp_path": str(itp_path),
                "itp_moleculetype_name": itp_moleculetype_name or mol_name,
                "atom_count": atom_count,
                "total_charge": total_charge,
                "has_preprocessor": has_preprocessor,
                "unparsed_atoms_line_count": unparsed_atoms_line_count,
                "parse_status": parse_status or status,
                "parse_reasons": list(parse_reasons or []),
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
                    itp_moleculetype_name=mol_name,
                )
                continue

            if _matches_protected_pattern(mol_name, PROTECTED_SOLVENT_PATTERNS):
                _record(
                    mol_name,
                    itp_path,
                    "protected_solvent",
                    source="protected_solvent_pattern",
                    itp_moleculetype_name=mol_name,
                )
                continue

            if configured_patterns and _matches_protected_pattern(mol_name, configured_patterns):
                _record(
                    mol_name,
                    itp_path,
                    "protected_configured",
                    source="charge_fix_protect_resnames",
                    reason="Matched user-configured protected residue/molecule pattern.",
                    itp_moleculetype_name=mol_name,
                )
                continue

            try:
                parsed = self._parse_moleculetype_ir(itp_path, mol_name, ctx=ctx)
                issue_bits = list(parsed.uncertainty_reasons)
                if parsed.status != MOLECULETYPE_IR_STATUS_EXACT or parsed.molecule is None:
                    classification = "protected_polymer" if looks_polymer else "unknown"
                    msg = (
                        f"Charge-fix classification for active molecule '{mol_name}' is degraded: "
                        + "; ".join(issue_bits or ["unsupported or ambiguous [ atoms ] parse"])
                    )
                    if strict_mode:
                        raise SanitizerError(msg)
                    print(f"  [WARN][DEGRADED] {msg}")
                    _record(
                        mol_name,
                        itp_path,
                        classification,
                        status="degraded",
                        source="moleculetype_ir_degraded",
                        reason="; ".join(issue_bits or ["unsupported or ambiguous [ atoms ] parse"]),
                        atom_count=(
                            len(parsed.molecule.atoms)
                            if parsed.molecule is not None
                            else 0
                        ),
                        total_charge=parsed.total_charge,
                        has_preprocessor=any(
                            "preprocessor directive" in item for item in issue_bits
                        ),
                        unparsed_atoms_line_count=len(
                            [item for item in issue_bits if item.startswith("[ atoms ] line ")]
                        ),
                        itp_moleculetype_name=parsed.actual_name or mol_name,
                        parse_status=parsed.status,
                        parse_reasons=issue_bits,
                    )
                    continue

                n_atoms = len(parsed.molecule.atoms)
                q_total = parsed.total_charge if parsed.total_charge is not None else 0.0

                if n_atoms == 1 and abs(q_total) > tol:
                    _record(
                        mol_name,
                        itp_path,
                        "monatomic_ion",
                        source="moleculetype_ir",
                        atom_count=n_atoms,
                        total_charge=q_total,
                        itp_moleculetype_name=parsed.actual_name or mol_name,
                        parse_status=parsed.status,
                    )
                    continue

                rounded_q = round(q_total)
                if abs(rounded_q) >= 1 and abs(q_total - rounded_q) < tol:
                    _record(
                        mol_name,
                        itp_path,
                        "ionic",
                        source="moleculetype_ir",
                        atom_count=n_atoms,
                        total_charge=q_total,
                        itp_moleculetype_name=parsed.actual_name or mol_name,
                        parse_status=parsed.status,
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
                        itp_moleculetype_name=parsed.actual_name or mol_name,
                        parse_status=parsed.status,
                    )
                    continue

                _record(
                    mol_name,
                    itp_path,
                    "neutral",
                    source="moleculetype_ir",
                    atom_count=n_atoms,
                    total_charge=q_total,
                    itp_moleculetype_name=parsed.actual_name or mol_name,
                    parse_status=parsed.status,
                )

            except SanitizerError:
                raise
            except Exception as exc:
                classification = "protected_polymer" if looks_polymer else "unknown"
                msg = (
                    f"Charge-fix classification for active molecule '{mol_name}' is degraded: "
                    f"{str(exc).splitlines()[0] if str(exc) else 'unexpected classification error'}"
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
                    reason=str(exc).splitlines()[0] if str(exc) else "unexpected classification error",
                    itp_moleculetype_name=mol_name,
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

    def _derive_charge_fix_checker_allowlist(
        self,
        classifications: Dict[str, str],
        classification_audit: Optional[Dict[str, Any]],
        explicit_allowlist: Optional[Set[str]],
    ) -> Set[str]:
        """
        Derive the checker allowlist without re-opening protected/degraded auto-targets.

        Explicit user allowlists pass through unchanged. Without an explicit allowlist,
        keep the legacy default target vocabulary but only for exact neutral molecules.
        """
        if explicit_allowlist is not None:
            return set(explicit_allowlist)

        from ..charge_neutrality import DEFAULT_CORRECTION_TARGET_ALLOWLIST

        default_allowlist = {
            token.casefold() for token in DEFAULT_CORRECTION_TARGET_ALLOWLIST
        }
        detail_by_casefold = {
            name.casefold(): detail
            for name, detail in (
                classification_audit.get("details", {}).items()
                if isinstance(classification_audit, dict)
                else []
            )
        }
        effective: Set[str] = set()
        for alias_name, classification in classifications.items():
            detail = detail_by_casefold.get(alias_name.casefold(), {})
            status = str(detail.get("status", "degraded"))
            if classification != "neutral" or status != "exact":
                continue
            actual_name = str(detail.get("itp_moleculetype_name") or alias_name)
            alias_key = alias_name.casefold()
            actual_key = actual_name.casefold()
            if alias_key in default_allowlist or actual_key in default_allowlist:
                effective.add(alias_name)
                effective.add(actual_name)
        return effective

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
        Compute a deterministic signature for one moleculetype using the IR parser.
        """
        parsed = self._parse_moleculetype_ir(itp_path, mol_name, ctx=ctx)
        if parsed.status != MOLECULETYPE_IR_STATUS_EXACT:
            return None
        return parsed.signature

    def _get_ordered_molecules_from_truth_source(
        self,
        top_path: Path,
        *,
        output_dir: Optional[Path] = None,
        ctx: Optional["PipelineContext"] = None,
    ) -> Tuple[TopMoleculesParseResult, TopologyTruthSource]:
        """Prefer preprocessed topology when available, otherwise use conservative raw parsing."""
        truth_source = self._select_topology_truth_source(
            top_path,
            output_dir=output_dir,
            ctx=ctx,
        )
        if truth_source.is_preprocessed and hasattr(self, "_parse_preprocessed_molecules"):
            ordered = self._parse_preprocessed_molecules(truth_source.path)
            if ordered:
                return TopMoleculesParseResult(ordered=ordered), truth_source
            print(
                "  [WARN][DEGRADED] Preferred preprocessed topology truth source was present "
                f"but could not be parsed for [molecules]: {truth_source.path}. "
                "Falling back to conservative Python parsing of system.top."
            )
            truth_source = TopologyTruthSource(
                path=top_path,
                source="python_fallback_raw_topology",
                is_preprocessed=False,
                fallback_used=True,
                reason=f"preferred preprocessed topology could not be parsed: {truth_source.path}",
                candidate_paths=truth_source.candidate_paths,
            )
        return self._get_ordered_molecules_from_top(top_path), truth_source

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
        conditional_depth = 0
        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return TopMoleculesParseResult(ordered=[])
        for line_no, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            directive = preprocessor_directive_token(line)
            section = parse_section_name(line)
            if section is not None:
                in_molecules = (section == "molecules")
                continue
            if not in_molecules:
                continue
            if directive in PREPROCESSOR_UNRESOLVED_SEMANTIC_DIRECTIVES:
                if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                    conditional_depth += 1
                elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE and conditional_depth > 0:
                    conditional_depth -= 1
                uncertain_reasons.append(
                    f"{directive} in [molecules] at {top_path.name}:{line_no}"
                )
                continue
            if directive is not None or conditional_depth > 0:
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
        if conditional_depth > 0:
            uncertain_reasons.append(
                f"unclosed conditional block in [molecules] at {top_path.name}"
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
        scan_lines, _ = self._conservative_topology_lines(
            lines=content.splitlines(),
            file_path=system_top_path,
            purpose="system.top include validation",
            strict=strict_includes,
        )
        for _line_no, line in scan_lines:
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

    def _format_lj_row_location(
        self,
        source_name: str,
        line_no: int,
        provenance: Optional[str] = None,
    ) -> str:
        """Format a stable source location for LJ diagnostics."""
        location = f"{source_name}:{line_no}"
        if provenance:
            return f"{location} (from {provenance})"
        return location

    def _scan_direct_include_for_defaults_ownership(
        self,
        include_path: Path,
        strict: bool,
    ) -> Tuple[bool, List[str]]:
        """
        Inspect one directly included file for obvious defaults ownership.

        Returns:
            (owns_defaults, ambiguity_reasons)
        """
        try:
            content = self.safe_read_text(include_path, strict=True)
        except SanitizerError as exc:
            return False, [f"{include_path.name}: {exc}"]

        if not content:
            return False, []

        lines = content.splitlines()
        _, issues = self._conservative_topology_lines(
            lines=lines,
            file_path=include_path,
            purpose="defaults ownership scan",
            strict=strict,
        )
        if issues:
            return False, list(issues)

        current_section: Optional[str] = None
        for line_no, line in enumerate(lines, start=1):
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if section == "defaults":
                    return True, []
                if section != "defaults":
                    break
                continue

            if current_section is not None:
                continue

            include_target = parse_include_target(line)
            if include_target:
                return (
                    False,
                    [
                        f"{include_path.name}:{line_no} contains nested pre-section include "
                        f"'{include_target}'"
                    ],
                )

        if ItpParser.extract_defaults_from_file(include_path) is not None:
            return True, []
        return False, []

    def _scan_existing_top_defaults_ownership(
        self,
        *,
        system_top_path: Path,
        lines: List[str],
        ignore_ranges: Optional[List[Tuple[int, int]]] = None,
        ctx: Optional["PipelineContext"] = None,
        strict: bool = False,
    ) -> TopDefaultsOwnershipScan:
        """Inspect existing system.top for defaults ownership and insertion anchoring."""
        skip_indices: Set[int] = set()
        if ignore_ranges:
            for start, end in ignore_ranges:
                skip_indices.update(range(start, end + 1))

        include_dirs = [system_top_path.parent]
        if ctx is not None:
            include_dirs = self._build_include_search_paths(ctx, system_top_path)

        allow_shadowing = getattr(ctx, "allow_include_shadowing", False) if ctx else False
        has_defaults = False
        has_forcefield_include = False
        owner_insert_after_idx: Optional[int] = None
        ambiguity_reasons: List[str] = []
        conditional_depth = 0
        current_section: Optional[str] = None

        def _note_owner(idx: int) -> None:
            nonlocal owner_insert_after_idx
            candidate = idx + 1
            if owner_insert_after_idx is None or candidate > owner_insert_after_idx:
                owner_insert_after_idx = candidate

        for idx, line in enumerate(lines):
            if idx in skip_indices:
                continue

            directive = preprocessor_directive_token(line)
            if directive == "#include":
                if conditional_depth > 0:
                    continue
                include_target = parse_include_target(line)
                if include_target is None:
                    continue
                include_target_lower = include_target.lower()
                if "forcefield.itp" in include_target_lower or ".ff/" in include_target_lower:
                    has_forcefield_include = True
                    _note_owner(idx)
                    continue
                if current_section is not None:
                    continue
                resolved = self._resolve_include_in_dirs(
                    include_target,
                    system_top_path,
                    include_dirs,
                    strict=False,
                    check_shadowing=True,
                    allow_shadowing=allow_shadowing,
                    ctx=ctx,
                )
                if resolved is None:
                    ambiguity_reasons.append(
                        f"{system_top_path.name}:{idx + 1} direct pre-section include "
                        f"'{include_target}' could not be resolved"
                    )
                    continue
                owns_defaults, nested_issues = self._scan_direct_include_for_defaults_ownership(
                    resolved,
                    strict=strict,
                )
                if owns_defaults:
                    has_defaults = True
                    _note_owner(idx)
                    continue
                for issue in nested_issues:
                    ambiguity_reasons.append(
                        f"{system_top_path.name}:{idx + 1} include '{include_target}' "
                        f"keeps [ defaults ] ownership ambiguous: {issue}"
                    )
                continue

            if directive is not None:
                if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                    conditional_depth += 1
                elif directive in PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES:
                    if conditional_depth == 0:
                        conditional_depth = 1
                elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE:
                    if conditional_depth > 0:
                        conditional_depth -= 1
                continue

            if conditional_depth > 0:
                continue

            section = parse_section_name(line)
            if section is not None:
                if section == "defaults":
                    has_defaults = True
                    current_section = section
                    _note_owner(idx)
                    continue
                break

            if current_section == "defaults":
                _note_owner(idx)

        return TopDefaultsOwnershipScan(
            has_defaults=has_defaults,
            has_forcefield_include=has_forcefield_include,
            owner_insert_after_idx=owner_insert_after_idx,
            ambiguity_reasons=tuple(ambiguity_reasons),
        )

    def _parse_lj_atomtypes_row(
        self,
        data: str,
        source_name: str,
        line_no: int,
    ) -> Tuple[str, Optional[str], float, float, Optional[str]]:
        """
        Parse one [ atomtypes ] data row for LJ validation.

        Conservative parser:
        - Supports tail shapes:
          * ... mass charge ptype sigma epsilon
          * ... mass charge sigma epsilon  (missing ptype)
        - Requires ptype, when present, in the canonical third token from tail
        - Rejects ambiguous rows rather than guessing shifted columns
        """
        parts = data.split()
        if len(parts) < 5:
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

        def _is_integer_like_token(token: str) -> bool:
            return re.fullmatch(r"[+-]?\d+", token) is not None

        def _is_explicit_charge_literal(token: str) -> bool:
            stripped = token.strip()
            if not stripped:
                return False
            if stripped[0] in "+-":
                return True
            try:
                numeric = float(stripped)
            except ValueError:
                return False
            return numeric == 0.0

        ptype: Optional[str]
        middle_tokens: List[str]
        if len(parts) >= 6 and parts[-3].upper() in PTYPE_TOKENS:
            ptype = parts[-3].upper()
            if len(parts) == 6:
                leading_token = parts[1]
                charge_token = parts[2]
                if _is_integer_like_token(leading_token) and not _is_explicit_charge_literal(charge_token):
                    raise ValueError(
                        f"Ambiguous 6-field explicit-ptype row at {source_name}:{line_no}: "
                        "the numeric fields before ptype could be interpreted as either "
                        "'mass charge' or 'at_num mass'. Use a signed/zero-like charge token "
                        "or include an additional non-numeric field before mass/charge. "
                        f"Raw: {data}"
                    )
            middle_tokens = parts[1:-5]
            charge_fields = _extract_floats_from_text(parts[-4])
            mass_fields = _extract_floats_from_text(parts[-5])
        else:
            ptype = None
            if len(parts) == 5:
                leading_token = parts[1]
                charge_token = parts[2]
                if _is_integer_like_token(leading_token) and not _is_explicit_charge_literal(charge_token):
                    raise ValueError(
                        f"Ambiguous 5-field missing-ptype row at {source_name}:{line_no}: "
                        "the numeric fields before the LJ tail could be interpreted as either "
                        "'mass charge' or 'at_num mass'. Use a signed/zero-like charge token "
                        "or include an additional non-numeric field before mass/charge. "
                        f"Raw: {data}"
                    )
            middle_tokens = parts[1:-4]
            charge_fields = _extract_floats_from_text(parts[-3]) if len(parts) >= 4 else []
            mass_fields = _extract_floats_from_text(parts[-4]) if len(parts) >= 5 else []

        extra_ptypes = [token for token in middle_tokens if token.upper() in PTYPE_TOKENS]
        if extra_ptypes:
            if ptype is None:
                raise ValueError(
                    f"Ambiguous row without explicit ptype at {source_name}:{line_no}: "
                    f"middle token(s) look like ptype markers {extra_ptypes}. Raw: {data}"
                )
            raise ValueError(
                f"Ambiguous row with multiple ptype-like tokens at {source_name}:{line_no}: "
                f"{extra_ptypes + [ptype]}. Raw: {data}"
            )

        if len(charge_fields) != 1 or len(mass_fields) != 1:
            raise ValueError(
                f"Malformed mass/charge columns in {source_name}:{line_no}. Raw: {data}"
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
                warnings.append(msg)
                continue
            if parse_warning:
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
                    "line_no": line_no,
                    "source_name": combined_atomtypes_path.name,
                }
            )

        if not atomtype_data and not warnings:
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
            "top_parse_issues": [],
            "warning_count": 0,
            "parse_issue_count": 0,
            "unit_error_count": 0,
        }

        report["parse_issue_count"] = len(warnings)
        report["top_parse_issues"] = warnings[:10]
        if not atomtype_data:
            report["warning_count"] = len(warnings)
            report["top_offenders"] = [{"message": msg} for msg in warnings[:10]]
            if ctx and hasattr(ctx, "manifest") and ctx.manifest:
                ctx.manifest.set_sanitizer_output("lj_validation", report)
            if policy == "off":
                return []
            for msg in warnings[:20]:
                print(f"  [WARN] {msg}")
            if len(warnings) > 20:
                print(
                    f"  [WARN] ... {len(warnings) - 20} additional LJ outlier warning(s) omitted."
                )
            return warnings

        offenders: List[Dict[str, Any]] = []
        unit_errors: List[Dict[str, Any]] = []
        unit_sigma_max = float(LJ_UNIT_SCREAMING_THRESHOLDS["sigma_max_nm"])
        unit_eps_max = float(LJ_UNIT_SCREAMING_THRESHOLDS["epsilon_max_kj_mol"])

        # comb-rule=1: retain robust statistical fallback with finite-value guardrails.
        if comb_rule == 1:
            if policy != "off":
                print(
                    "  [INFO] comb-rule=1: treating last two fields as raw LJ params "
                    "(p1/p2), not sigma/epsilon"
                )
            for row in atomtype_data:
                name = str(row["name"])
                p1 = float(row["p1"])
                p2 = float(row["p2"])
                location = self._format_lj_row_location(
                    str(row.get("source_name", combined_atomtypes_path.name)),
                    int(row.get("line_no", 0)),
                    row.get("provenance"),
                )
                if any((math.isnan(p1), math.isnan(p2), math.isinf(p1), math.isinf(p2))):
                    msg = (
                        f"Non-finite comb-rule=1 LJ params in atomtype '{name}' at {location}: "
                        f"p1={p1}, p2={p2}."
                    )
                    unit_errors.append(
                        {
                            "name": name,
                            "ptype": row.get("ptype"),
                            "sigma_nm": p1,
                            "epsilon_kj_mol": p2,
                            "provenance": row.get("provenance"),
                            "kind": "nan_or_inf",
                            "message": msg,
                        }
                    )
            warnings.extend(self._check_lj_outliers_robust(atomtype_data))
            report["warning_count"] = len(warnings)
            report["warning_examples"] = warnings[:10]
            report["top_offenders"] = [{"message": msg} for msg in warnings[:10]]
            report["unit_error_count"] = len(unit_errors)
            if unit_errors:
                report["top_unit_errors"] = unit_errors[:10]
            if ctx and hasattr(ctx, "manifest") and ctx.manifest:
                ctx.manifest.set_sanitizer_output("lj_validation", report)
            if unit_errors:
                first = unit_errors[0]
                raise SanitizerError(
                    "LJ unit-screaming threshold violated: "
                    f"{first.get('message')}"
                )
            if policy == "off":
                return []
            if warnings and policy == "error":
                raise SanitizerError(
                    "LJ outlier policy=error triggered for comb-rule=1 robust outlier checks. "
                    f"First issue: {warnings[0]}"
                )
            for msg in warnings[:20]:
                print(f"  [WARN] {msg}")
            if len(warnings) > 20:
                print(
                    f"  [WARN] ... {len(warnings) - 20} additional LJ outlier warning(s) omitted."
                )
            return warnings

        for row in atomtype_data:
            name = str(row["name"])
            ptype = str(row.get("ptype") or "A").upper()
            sigma = float(row["p1"])
            epsilon = float(row["p2"])
            provenance = row.get("provenance")
            location = self._format_lj_row_location(
                str(row.get("source_name", combined_atomtypes_path.name)),
                int(row.get("line_no", 0)),
                provenance,
            )

            if ptype in {"D", "V"}:
                # Skip dummy/virtual sites to avoid false positives in modern topologies.
                continue

            if any((math.isnan(sigma), math.isnan(epsilon), math.isinf(sigma), math.isinf(epsilon))):
                msg = (
                    f"NaN/Inf LJ params in atomtype '{name}' at {location}: "
                    f"sigma={sigma}, epsilon={epsilon}."
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
                    f"Unit-screaming LJ value in atomtype '{name}' at {location}: "
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
                if epsilon <= 0:
                    msg = (
                        f"Negative sigma in atomtype '{name}' at {location}: "
                        f"sigma={sigma:.6g} nm, epsilon={epsilon:.6g} kJ/mol. "
                        "Legacy/special encodings normally still use epsilon>0, so this row "
                        "is treated as suspicious."
                    )
                    kind = "negative_sigma_nonpositive_epsilon"
                else:
                    msg = (
                        f"Negative sigma in atomtype '{name}' at {location}: "
                        f"sigma={sigma:.6g} nm, epsilon={epsilon:.6g} kJ/mol. "
                        "Some GROMACS forcefields use negative sigma as a special encoding; "
                        "review manually."
                    )
                    kind = "negative_sigma_special"
                offenders.append(
                    {
                        "name": name,
                        "ptype": ptype,
                        "sigma_nm": sigma,
                        "epsilon_kj_mol": epsilon,
                        "provenance": provenance,
                        "kind": kind,
                        "severity": "warn",
                        "score": 0.0,
                        "message": msg,
                    }
                )
                warnings.append(msg)
                continue

            if epsilon < 0:
                msg = (
                    f"Negative epsilon in atomtype '{name}' at {location}: "
                    f"sigma={sigma:.6g} nm, epsilon={epsilon:.6g} kJ/mol."
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
                    f"LJ profile outlier '{name}' at {location}: "
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

        if policy != "off" and top_offenders:
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
                f"first offender: {top_offenders[0].get('message') if top_offenders else 'n/a'}"
            )

        for msg in warnings[:20]:
            print(f"  [WARN] {msg}")
        if len(warnings) > 20:
            print(
                f"  [WARN] ... {len(warnings) - 20} additional LJ outlier warning(s) omitted."
            )

        return warnings
    
    def _check_lj_outliers_robust(
        self,
        atomtype_data: List[Dict[str, Any]],
    ) -> List[str]:
        """
        Check for LJ parameter outliers using robust statistics (MAD z-score on log-scale).
        
        For comb-rule=1 where we cannot assume sigma/epsilon semantics.
        Uses median and median absolute deviation (MAD) to detect outliers.
        
        Args:
            atomtype_data: Parsed row dictionaries with name/p1/p2/location metadata
            
        Returns:
            List of warning messages
        """
        warnings: List[str] = []
        
        # Filter to positive values only for log-scale analysis
        def _row_location(row: Dict[str, Any]) -> str:
            return self._format_lj_row_location(
                str(row.get("source_name", "?")),
                int(row.get("line_no", 0)),
                row.get("provenance"),
            )

        p1_positive = [
            (row, abs(float(row["p1"])))
            for row in atomtype_data
            if math.isfinite(float(row["p1"])) and float(row["p1"]) != 0.0
        ]
        p2_positive = [
            (row, abs(float(row["p2"])))
            for row in atomtype_data
            if math.isfinite(float(row["p2"])) and float(row["p2"]) != 0.0
        ]
        
        # Check for negative values (warn-only, not "impossible")
        for row in atomtype_data:
            name = str(row["name"])
            p1 = float(row["p1"])
            p2 = float(row["p2"])
            if p1 < 0 or p2 < 0:
                msg = (
                    f"Negative comb-rule=1 LJ params in atomtype '{name}' at {_row_location(row)}: "
                    f"p1={p1:.4e}, p2={p2:.4e}. comb-rule=1 is not sigma/epsilon; review manually."
                )
                warnings.append(msg)
        
        # Robust outlier detection using MAD z-score on log scale
        for label, values in [("p1", p1_positive), ("p2", p2_positive)]:
            if len(values) < 3:
                continue  # Not enough data for statistics

            positive_values = [value for _, value in values]
            median_positive = _median(positive_values)
            near_zero_floor = max(
                LJ_ROBUST_LOG_ABSOLUTE_FLOOR,
                median_positive * LJ_ROBUST_LOG_MEDIAN_RELATIVE_FLOOR,
            )
            filtered_values = [
                (row, value)
                for row, value in values
                if value >= near_zero_floor
            ]
            excluded_values = [
                (row, value)
                for row, value in values
                if value < near_zero_floor
            ]
            if excluded_values:
                examples = ", ".join(
                    f"{row['name']}@{_row_location(row)}={value:.4e}"
                    for row, value in excluded_values[:3]
                )
                msg = (
                    f"Near-zero comb-rule=1 {label} value(s) excluded from log-domain MAD: "
                    f"floor={near_zero_floor:.4e} derived from median |{label}|={median_positive:.4e}; "
                    f"examples: {examples}"
                )
                warnings.append(msg)
            if len(filtered_values) < 3:
                continue

            # Compute log-scale values
            log_vals = [math.log10(value) for _, value in filtered_values]
            median = _median(log_vals)
            abs_devs = [abs(log_val - median) for log_val in log_vals]
            mad = _median(abs_devs)

            # Scale factor for MAD -> sigma-equivalent (1.4826 for normal distribution)
            mad_scaled = mad * 1.4826 if mad > 0 else 1.0

            # Check for outliers
            for (row, val), log_val in zip(filtered_values, log_vals):
                if mad_scaled > 0:
                    z_score = abs(log_val - median) / mad_scaled
                    if z_score > LJ_ROBUST_Z_THRESHOLD:
                        msg = (
                            f"Statistical LJ outlier '{row['name']}' for {label} at {_row_location(row)}: "
                            f"value={val:.4e} (log-MAD z={z_score:.1f})"
                        )
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

    def _build_include_search_paths(
        self,
        ctx: Optional["PipelineContext"],
        base_path: Path,
        extra_include_dirs: Optional[List[Path]] = None,
    ) -> List[Path]:
        """
        Build per-file include search paths.

        This is the single source of include-priority semantics:
        - `sanitized_first` places configured user include dirs immediately after the
          currently included file's parent directory.
        - `forcefield_first` keeps configured user include dirs after system/forcefield
          roots.

        `extra_include_dirs` are compatibility-only append roots. They never define a
        second priority mechanism.
        """
        paths: List[Path] = [base_path.parent]
        if ctx is not None:
            user_paths, gmx_share_paths = self._collect_ctx_include_roots(ctx)
            raw_priority = getattr(ctx, "itp_include_dirs_priority", DEFAULT_INCLUDE_PRIORITY)
            priority = normalize_include_priority(raw_priority)

            system_gromacs_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs")
            system_htpolynet_dir = ctx.get_input_path("systems", ctx.system_id, "htpolynet")
            configured_paths: List[Path] = [
                system_gromacs_dir / "top",
                system_gromacs_dir / "itp",
                system_htpolynet_dir / "itp",
            ]

            if hasattr(self, "_forcefield_dir") and self._forcefield_dir:
                configured_paths.append(self._forcefield_dir / "gromacs")
            if hasattr(self, "_sanitized_current_dir") and self._sanitized_current_dir:
                configured_paths.append(self._sanitized_current_dir)
            configured_paths.extend(gmx_share_paths)

            if priority == "sanitized_first":
                paths.extend(user_paths)
                paths.extend(configured_paths)
            else:
                paths.extend(configured_paths)
                paths.extend(user_paths)

        for candidate in extra_include_dirs or []:
            paths.append(self._resolve_include_root_path(candidate, ctx))

        return self._dedupe_existing_paths(paths)




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

    @staticmethod
    def _normalize_top_include_target(include_target: str) -> str:
        """Normalize a system.top include target for deterministic emission."""
        normalized = str(include_target or "").strip().replace("\\", "/")
        if not normalized:
            return ""
        return PurePosixPath(normalized).as_posix()

    def _relative_top_include_target(
        self,
        include_path: Path,
        *,
        top_dir: Optional[Path],
        fallback_dir: Optional[str] = None,
    ) -> str:
        """Return a normalized include target relative to system.top."""
        if top_dir is not None:
            rel_path = Path(os.path.relpath(include_path, top_dir)).as_posix()
            return self._normalize_top_include_target(rel_path)
        if fallback_dir:
            return self._normalize_top_include_target(
                f"{fallback_dir.rstrip('/')}/{include_path.name}"
            )
        return self._normalize_top_include_target(include_path.name)

    def _dedupe_include_targets(self, include_targets: Iterable[str]) -> List[str]:
        """Deduplicate include targets while preserving first-seen order."""
        ordered: List[str] = []
        seen: Set[str] = set()
        for include_target in include_targets:
            normalized = self._normalize_top_include_target(include_target)
            if not normalized or normalized in seen:
                continue
            seen.add(normalized)
            ordered.append(normalized)
        return ordered

    def _ordered_output_extra_includes(self, top_dir: Path, result: Any) -> List[str]:
        """
        Build deterministic sanitizer-generated sidecar includes for system.top.

        Order is explicit:
        1. Prefix-injected bonded/type tables
        2. Mixed-defaults cross-group nonbond_params override
        3. Mixed-defaults within-group preservation nonbond_params
        """
        ordered_targets: List[str] = []
        for attr_name in (
            "prefix_injected_types_current_path",
            "nonbond_params_cross_group_current_path",
            "nonbond_params_secondary_current_path",
        ):
            raw_path = getattr(result, attr_name, None)
            if not raw_path:
                continue
            include_path = Path(raw_path)
            if not include_path.exists():
                continue
            ordered_targets.append(
                self._relative_top_include_target(
                    include_path,
                    top_dir=top_dir,
                )
            )
        return self._dedupe_include_targets(ordered_targets)

    def _insert_managed_block_at_index(
        self,
        lines: List[str],
        managed_lines: List[str],
        begin_marker: str,
        end_marker: str,
        insert_idx: int,
    ) -> List[str]:
        """Insert a managed block at an exact line index with stable spacing."""
        result_lines = list(lines[:insert_idx])
        if result_lines and result_lines[-1].strip():
            result_lines.append("")
        result_lines.append(begin_marker)
        result_lines.extend(managed_lines)
        result_lines.append(end_marker)
        if insert_idx < len(lines) and lines[insert_idx].strip():
            result_lines.append("")
        result_lines.extend(lines[insert_idx:])
        return result_lines

    def _managed_block_payload_is_safe(self, payload_lines: List[str]) -> bool:
        """Return True when removed malformed-block payload is managed-block shaped."""
        in_defaults = False
        for line in payload_lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            section = parse_section_name(line)
            if section is not None:
                in_defaults = section == "defaults"
                if in_defaults:
                    continue
                return False
            if parse_include_target(line):
                in_defaults = False
                continue
            if preprocessor_directive_token(line) is not None:
                return False
            if in_defaults:
                continue
            return False
        return True

    def _cleanup_malformed_managed_block(
        self,
        lines: List[str],
        begin_marker: str,
        end_marker: str,
    ) -> Optional[List[str]]:
        """
        Remove malformed managed-block fragments only when the removed payload is safe.

        This prevents non-strict rebuild mode from silently deleting user topology
        content such as `[ system ]` / `[ molecules ]` sections.
        """
        cleaned_lines: List[str] = []
        removed_payload: List[str] = []
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
                removed_payload.append(line)
                continue
            cleaned_lines.append(line)
        if not self._managed_block_payload_is_safe(removed_payload):
            return None
        return cleaned_lines




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

        _, top_preprocessor_issues = self._conservative_topology_lines(
            lines=lines,
            file_path=existing_top,
            purpose="system.top managed-block update",
            strict=bool(strict_top_update),
        )

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
            cleaned_lines = self._cleanup_malformed_managed_block(
                lines,
                begin_marker,
                end_marker,
            )
            if cleaned_lines is None:
                raise SanitizerError(
                    msg
                    + "\nCannot safely rebuild non-strict system.top because malformed managed-block "
                    "fragments overlap non-managed topology content."
                )
            print(f"  [WARN] {msg}")
            print("  [WARN] Rebuilding sanitizer managed block from safe managed-content fragments only.")

            ownership_scan = self._scan_existing_top_defaults_ownership(
                system_top_path=existing_top,
                lines=cleaned_lines,
                ignore_ranges=None,
                ctx=ctx,
                strict=bool(strict_top_update),
            )
            if ownership_scan.ambiguity_reasons:
                ambiguity_msg = (
                    "system.top defaults ownership is ambiguous outside sanitizer-managed "
                    "content.\n  - "
                    + "\n  - ".join(ownership_scan.ambiguity_reasons)
                )
                if strict_top_update:
                    raise SanitizerError(ambiguity_msg)
                print(f"  [WARN][DEGRADED] {ambiguity_msg}")
            emit_defaults = not (
                ownership_scan.has_defaults or ownership_scan.has_forcefield_include
            )
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
                "preprocessor_ambiguity": list(top_preprocessor_issues),
                "defaults_ownership_ambiguity": list(ownership_scan.ambiguity_reasons),
            }
            rebuilt = self._insert_managed_block_at_safe_location(
                cleaned_lines,
                managed_lines,
                begin_marker,
                end_marker,
                owner_insert_after_idx=ownership_scan.owner_insert_after_idx,
            )
            return '\n'.join(rebuilt) + '\n'

        if len(pairs) > 1 and strict_top_update:
            raise SanitizerError(
                f"Multiple sanitizer managed blocks detected ({len(pairs)}). "
                "Replacing all managed blocks with one clean block."
            )

        # Detect pre-existing defaults/forcefield include outside managed block(s).
        # If found, do not inject a duplicate [ defaults ] section.
        ownership_scan = self._scan_existing_top_defaults_ownership(
            system_top_path=existing_top,
            lines=lines,
            ignore_ranges=pairs if pairs else None,
            ctx=ctx,
            strict=bool(strict_top_update),
        )
        if ownership_scan.ambiguity_reasons:
            ambiguity_msg = (
                "system.top defaults ownership is ambiguous outside sanitizer-managed "
                "content.\n  - "
                + "\n  - ".join(ownership_scan.ambiguity_reasons)
            )
            if strict_top_update:
                raise SanitizerError(ambiguity_msg)
            print(f"  [WARN][DEGRADED] {ambiguity_msg}")
        emit_defaults = not (
            ownership_scan.has_defaults or ownership_scan.has_forcefield_include
        )

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
                    "degraded": bool(
                        top_preprocessor_issues or ownership_scan.ambiguity_reasons
                    ),
                    "preprocessor_ambiguity": list(top_preprocessor_issues),
                    "defaults_ownership_ambiguity": list(ownership_scan.ambiguity_reasons),
                }
            else:
                # Defensive fallback: bad ordering
                result_lines = self._insert_managed_block_at_safe_location(
                    lines,
                    managed_lines,
                    begin_marker,
                    end_marker,
                    owner_insert_after_idx=ownership_scan.owner_insert_after_idx,
                )
                self._system_top_update_status = {
                    "mode": "fallback_insert_bad_marker_order",
                    "degraded": True,
                    "preprocessor_ambiguity": list(top_preprocessor_issues),
                    "defaults_ownership_ambiguity": list(ownership_scan.ambiguity_reasons),
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

            cleaned_lines = [
                line for idx, line in enumerate(lines) if idx not in skip_indices
            ]
            removed_before_insertion = sum(1 for idx in skip_indices if idx < insertion_at)
            cleaned_insert_idx = insertion_at - removed_before_insertion
            result_lines = self._insert_managed_block_at_index(
                cleaned_lines,
                managed_lines,
                begin_marker,
                end_marker,
                cleaned_insert_idx,
            )
            self._system_top_update_status = {
                "mode": "collapsed_multiple_managed_blocks",
                "degraded": True,
                "blocks_found": len(pairs),
                "preprocessor_ambiguity": list(top_preprocessor_issues),
                "defaults_ownership_ambiguity": list(ownership_scan.ambiguity_reasons),
            }
        else:
            # Markers not found: insert at safe location
            result_lines = self._insert_managed_block_at_safe_location(
                lines,
                managed_lines,
                begin_marker,
                end_marker,
                owner_insert_after_idx=ownership_scan.owner_insert_after_idx,
            )
            self._system_top_update_status = {
                "mode": "inserted_new_managed_block",
                "degraded": bool(
                    top_preprocessor_issues or ownership_scan.ambiguity_reasons
                ),
                "preprocessor_ambiguity": list(top_preprocessor_issues),
                "defaults_ownership_ambiguity": list(ownership_scan.ambiguity_reasons),
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
        emitted_includes: Set[str] = set()

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
        combined_rel = self._relative_top_include_target(
            combined_atomtypes_path,
            top_dir=top_dir,
        )

        if emit_defaults:
            managed_lines.append("; Combined atomtypes (after defaults for correct LJ interpretation)")
        else:
            managed_lines.append("; Combined atomtypes")
        managed_lines.extend([
            f'#include "{combined_rel}"',
            "",
        ])
        emitted_includes.add(combined_rel)

        normalized_extra_includes = self._dedupe_include_targets(extra_includes or [])
        extra_to_emit = [
            include_target
            for include_target in normalized_extra_includes
            if include_target not in emitted_includes
        ]
        if extra_to_emit:
            managed_lines.append("; Extra includes (sanitizer generated)")
            for include_target in extra_to_emit:
                managed_lines.append(f'#include "{include_target}"')
                emitted_includes.add(include_target)
            managed_lines.append("")

        managed_lines.append("; Molecule definitions (sanitized, [molecules]-ordered)")

        molecule_targets: List[str] = []
        seen_molecule_paths: Set[Path] = set()
        for itp_path in molecule_itp_files:
            resolved = itp_path.resolve()
            if resolved in seen_molecule_paths:
                continue
            seen_molecule_paths.add(resolved)
            molecule_targets.append(
                self._relative_top_include_target(
                    itp_path,
                    top_dir=top_dir,
                    fallback_dir="itp_sanitized_current",
                )
            )
        for include_target in molecule_targets:
            if include_target in emitted_includes:
                continue
            managed_lines.append(f'#include "{include_target}"')
            emitted_includes.add(include_target)

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
        conditional_depth = 0
        for idx, line in enumerate(lines):
            if idx in skip_indices:
                continue
            directive = preprocessor_directive_token(line)
            if directive == "#include":
                if conditional_depth > 0:
                    continue
                include_target = parse_include_target(line)
                if include_target:
                    target = include_target.lower()
                    if "forcefield.itp" in target or ".ff/" in target:
                        has_forcefield = True
                if has_defaults and has_forcefield:
                    break
                continue
            if directive is not None:
                if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                    conditional_depth += 1
                elif directive in PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES:
                    if conditional_depth == 0:
                        conditional_depth = 1
                elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE:
                    if conditional_depth > 0:
                        conditional_depth -= 1
                continue
            if conditional_depth > 0:
                continue
            section = parse_section_name(line)
            if section == "defaults":
                has_defaults = True
            if has_defaults and has_forcefield:
                break

        return has_defaults, has_forcefield
    
    def _insert_managed_block_at_safe_location(
        self,
        lines: List[str],
        managed_lines: List[str],
        begin_marker: str,
        end_marker: str,
        owner_insert_after_idx: Optional[int] = None,
    ) -> List[str]:
        """
        Insert the managed block at a safe location in system.top.
        
        Deterministic policy:
        1. If existing topology content already owns defaults (direct [ defaults ],
           direct pre-section include with direct defaults, or unconditional
           forcefield include), insert immediately after that ownership prefix.
        2. Otherwise insert before the first non-preamble line, where the preamble
           is limited to leading comments, blanks, and top-level macro definitions.
        3. Never insert after the first topology section header.
        """
        if owner_insert_after_idx is not None:
            return self._insert_managed_block_at_index(
                lines,
                managed_lines,
                begin_marker,
                end_marker,
                owner_insert_after_idx,
            )

        first_section_idx = len(lines)
        for idx, line in enumerate(lines):
            if parse_section_name(line) is not None:
                first_section_idx = idx
                break

        forcefield_include_idx = None
        conditional_depth = 0

        for idx, line in enumerate(lines[:first_section_idx]):
            directive = preprocessor_directive_token(line)
            if directive == "#include":
                if conditional_depth > 0:
                    continue
                include_target = parse_include_target(line)
                if include_target is None:
                    continue
                target_lower = include_target.lower()
                if "forcefield.itp" in target_lower or ".ff/" in target_lower:
                    forcefield_include_idx = idx
                continue
            if directive is not None:
                if directive in PREPROCESSOR_CONDITIONAL_START_DIRECTIVES:
                    conditional_depth += 1
                elif directive in PREPROCESSOR_CONDITIONAL_BRANCH_DIRECTIVES:
                    if conditional_depth == 0:
                        conditional_depth = 1
                elif directive == PREPROCESSOR_CONDITIONAL_END_DIRECTIVE:
                    if conditional_depth > 0:
                        conditional_depth -= 1
                continue

        if forcefield_include_idx is not None:
            insert_idx = forcefield_include_idx + 1
        else:
            insert_idx = 0
            while insert_idx < first_section_idx:
                stripped = lines[insert_idx].strip()
                directive = preprocessor_directive_token(lines[insert_idx])
                if (
                    not stripped
                    or stripped.startswith(";")
                    or directive in PREPROCESSOR_MACRO_DIRECTIVES
                ):
                    insert_idx += 1
                    continue
                break

        return self._insert_managed_block_at_index(
            lines,
            managed_lines,
            begin_marker,
            end_marker,
            insert_idx,
        )

    def _atomic_publish_current_directory(
        self,
        source_dir: Path,
        current_dir: Path,
    ) -> None:
        """Crash-safe publish of a current directory via temp copy + replace."""
        current_dir.parent.mkdir(parents=True, exist_ok=True)
        temp_dir = current_dir.parent / f".{current_dir.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
        backup_dir = current_dir.parent / f".{current_dir.name}.bak.{os.getpid()}.{uuid.uuid4().hex}"
        current_exists = current_dir.exists() or current_dir.is_symlink()
        try:
            shutil.copytree(source_dir, temp_dir, symlinks=True)
            _best_effort_fsync_dir(temp_dir)
            if current_exists:
                os.replace(current_dir, backup_dir)
                _best_effort_fsync_dir(current_dir.parent)
            os.replace(temp_dir, current_dir)
            _best_effort_fsync_dir(current_dir.parent)
        except Exception:
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
            if (backup_dir.exists() or backup_dir.is_symlink()) and not (
                current_dir.exists() or current_dir.is_symlink()
            ):
                os.replace(backup_dir, current_dir)
                _best_effort_fsync_dir(current_dir.parent)
            raise
        finally:
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
            if backup_dir.is_symlink():
                backup_dir.unlink()
            elif backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)
    
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
        self._atomic_publish_current_directory(sanitized_dir, sanitized_current)
        
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
    

