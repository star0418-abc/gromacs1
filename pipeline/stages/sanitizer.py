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
    Check if name matches any protected pattern using canonicalized boundaries.
    """
    import re as _re
    def _canonical_parts(value: str):
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
