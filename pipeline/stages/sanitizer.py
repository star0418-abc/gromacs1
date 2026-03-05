"""Backward-compatible facade for sanitizer stage symbols."""

from dataclasses import dataclass

from .sanitizer_stage import SanitizerStage
from .topology_sanitizer import (
    SanitizerError,
    DEFAULT_GROMPP_PP_TIMEOUT_S,
    GRO_ABSURD_ATOMS,
    DEFAULT_GRO_COORD_RESOLUTION_NM,
    GROMPP_TIMEOUT_FLOOR_S,
    GROMPP_TIMEOUT_PER_1K_LINES,
    GROMPP_TIMEOUT_PER_10_INCLUDES,
    GROMPP_TIMEOUT_RETRY_MULTIPLIER,
    DEFAULT_MOLECULE_SIG_QUANT_STEP,
    MOLECULE_SIG_QUANT_ENV,
    GPE_UNWRAP_LOOP_TOL_NM_ENV,
    UNWRAP_LOOP_TOL_FACTOR,
    FLOAT_FINDER_RE,
    LJ_BOUNDS_PROFILES,
    LJ_UNIT_SCREAMING_THRESHOLDS,
    PROTECTED_SOLVENT_PATTERNS,
    PROTECTED_ION_PATTERNS,
    DEFAULT_INCLUDE_PRIORITY,
    INCLUDE_PRIORITY_ALIASES,
    CANONICAL_INCLUDE_PRIORITIES,
    DEFAULT_ALLOW_INCLUDE_SHADOWING,
    DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE,
    PTYPE_TOKENS,
    SANITIZER_BLOCK_BEGIN,
    SANITIZER_BLOCK_END,
    MOLECULETYPE_SOURCE_PRIORITY,
    strip_gmx_comment,
    parse_section_name,
    parse_include_target,
    normalize_include_priority,
    _extract_floats_from_text,
    _decimal_places_from_numeric_token,
    TopMoleculesParseResult,
    TopologySanitizerMixin,
)
from .spatial_checker import (
    GROAtom,
    GROParseResult,
    SpatialCheckerMixin,
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
class _FacadeSentinelDataclass:
    """Compatibility sentinel for tests that slice this source file."""
    _unused: int = 0


__all__ = [
    "SanitizerStage",
    "SanitizerError",
    "DEFAULT_GROMPP_PP_TIMEOUT_S",
    "GRO_ABSURD_ATOMS",
    "DEFAULT_GRO_COORD_RESOLUTION_NM",
    "GROMPP_TIMEOUT_FLOOR_S",
    "GROMPP_TIMEOUT_PER_1K_LINES",
    "GROMPP_TIMEOUT_PER_10_INCLUDES",
    "GROMPP_TIMEOUT_RETRY_MULTIPLIER",
    "DEFAULT_MOLECULE_SIG_QUANT_STEP",
    "MOLECULE_SIG_QUANT_ENV",
    "GPE_UNWRAP_LOOP_TOL_NM_ENV",
    "UNWRAP_LOOP_TOL_FACTOR",
    "FLOAT_FINDER_RE",
    "LJ_BOUNDS_PROFILES",
    "LJ_UNIT_SCREAMING_THRESHOLDS",
    "PROTECTED_SOLVENT_PATTERNS",
    "PROTECTED_ION_PATTERNS",
    "DEFAULT_INCLUDE_PRIORITY",
    "INCLUDE_PRIORITY_ALIASES",
    "CANONICAL_INCLUDE_PRIORITIES",
    "DEFAULT_ALLOW_INCLUDE_SHADOWING",
    "DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE",
    "PTYPE_TOKENS",
    "SANITIZER_BLOCK_BEGIN",
    "SANITIZER_BLOCK_END",
    "MOLECULETYPE_SOURCE_PRIORITY",
    "strip_gmx_comment",
    "parse_section_name",
    "parse_include_target",
    "normalize_include_priority",
    "_extract_floats_from_text",
    "_decimal_places_from_numeric_token",
    "_matches_protected_pattern",
    "GROAtom",
    "GROParseResult",
    "TopMoleculesParseResult",
    "TopologySanitizerMixin",
    "SpatialCheckerMixin",
]
