"""
AnalysisStage: Computes RDF and coordination numbers.

Reads pair definitions from analysis_pairs.yaml and computes:
- RDF via gmx rdf
- Coordination numbers via gmx rdf -cn (GROMACS native)

Key features:
- Explicit group vs selection expression handling (with index.ndx support)
- Box-aware rmax capping using minimum box altitude for triclinic cells (strict invariant)
- Robust first_minimum detection with smoothing
- Full manifest provenance
"""

import difflib
import hashlib
import math
import re
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

import yaml

from .base import BaseStage
from ..gromacs_cmd import resolve_gmx_command
from ..utils import parse_bool as _parse_bool  # Shared canonical boolean parser

if TYPE_CHECKING:
    from ..context import PipelineContext


# Keywords/operators that indicate a gmx selection expression (not a simple group name)
_GMX_SELECTION_KEYWORDS = {
    "and", "or", "not", "resname", "resid", "name", "type", "atomnr",
    "group", "within", "same", "plus", "merge", "permute", "mol", "all",
    "byres", "byatom", "around", "insphere", "inbox", "exclusive", "inclusive",
}

_DEFAULT_RDF_MAX_R = 1.5
_DEFAULT_RDF_BIN = 0.002
_DEFAULT_CN_CUTOFF = 0.35
_DEFAULT_STRICT_ANALYSIS = True
_DEFAULT_GMX_TIMEOUT_S = 600
_MAX_FILENAME_LEN = 64
_HASH_SUFFIX_LEN = 8
_DEFAULT_PEAK_HEIGHT_THRESHOLD = 0.08
_DEFAULT_MIN_DEPTH = 0.05
_DEFAULT_BASELINE_WINDOW_NM = 0.10
_DEFAULT_TRUNCATION_TAIL_POINTS = 3
_DEFAULT_RECOVERY_WINDOW_NM = 0.08
_DEFAULT_RECOVERY_HEIGHT = 0.03
_DEFAULT_RECOVERY_POSITIVE_POINTS = 4
_DEFAULT_LOCAL_WINDOW_NM = 0.0
_DEFAULT_MAX_SEP_NM = 0.40


def _safe_filename(name: str, max_len: int = _MAX_FILENAME_LEN) -> str:
    """
    Sanitize a string for use as a filename component.
    
    - Keeps only [A-Za-z0-9._-]
    - Replaces other chars with '_'
    - Collapses multiple consecutive dots (prevents '..')
    - Strips leading dots and underscores
    - Strips trailing underscores
    - Prevents empty result (fallback to 'unnamed')
    - Truncates to max_len
    """
    if not name:
        return "unnamed"
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    # Collapse multiple consecutive dots (prevent path traversal like ..)
    sanitized = re.sub(r"\.{2,}", ".", sanitized)
    # Strip leading dots and underscores, trailing underscores
    sanitized = sanitized.lstrip("._").rstrip("_")
    if not sanitized:
        return "unnamed"
    if len(sanitized) > max_len:
        sanitized = sanitized[:max_len].rstrip("_.")
    return sanitized or "unnamed"


def _build_pair_output_basename(
    pair_name: str,
    group1_spec: Dict[str, str],
    group2_spec: Dict[str, str],
) -> Tuple[str, str]:
    """
    Build deterministic output basename for a pair.

    Format: <safe_prefix>__<stable_hash8>
    Hash includes name + normalized group specs to prevent collisions even if names repeat.
    """
    prefix_max_len = max(8, _MAX_FILENAME_LEN - (_HASH_SUFFIX_LEN + 2))
    safe_prefix = _safe_filename(pair_name, max_len=prefix_max_len)
    hash_input = (
        f"name={pair_name}|"
        f"g1={group1_spec.get('mode','')}:{group1_spec.get('value','')}|"
        f"g2={group2_spec.get('mode','')}:{group2_spec.get('value','')}"
    )
    suffix = hashlib.sha1(hash_input.encode("utf-8")).hexdigest()[:_HASH_SUFFIX_LEN]
    return safe_prefix, f"{safe_prefix}__{suffix}"


def _is_nonempty_file(path: Path) -> bool:
    """Return True if file exists and contains at least one byte."""
    return path.exists() and path.is_file() and path.stat().st_size > 0


# _parse_bool is now imported from ..utils (prevents bool("false") == True bug)


def parse_analysis_pairs(yaml_path: Path) -> Dict[str, Any]:
    """
    Parse analysis_pairs.yaml file.
    
    Args:
        yaml_path: Path to analysis_pairs.yaml
        
    Returns:
        Parsed configuration dict
    """
    if not yaml_path.exists():
        raise FileNotFoundError(f"Analysis config not found: {yaml_path}")
    
    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f)

    if config is None:
        return {}
    if not isinstance(config, dict):
        raise ValueError(f"Invalid analysis config (expected mapping): {yaml_path}")
    return config


def _is_selection_expression(s: str) -> bool:
    """
    Check if string looks like a gmx selection expression vs a simple group name.
    This is only used for backwards-compatible inference (string-only schema).
    """
    s_lower = s.lower().strip()
    if not s_lower:
        return False

    # Obvious selection syntax
    if any(ch in s_lower for ch in ("(", ")", "<", ">", "=", "\"", "'")):
        return True

    # "group ..." or other keyword operators imply selection expression
    tokens = re.split(r"\s+", s_lower)
    for tok in tokens:
        if tok in _GMX_SELECTION_KEYWORDS:
            return True

    return False


def _to_group_selection(group_name: str) -> str:
    """Convert a group name to a gmx selection string with safe quoting."""
    if '"' in group_name:
        raise ValueError(f'Group name contains a double quote: {group_name!r}')
    return f'group "{group_name}"'


def _extract_floats(text: str) -> List[float]:
    """Extract floats from a string, supporting scientific notation."""
    matches = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    values: List[float] = []
    for m in matches:
        try:
            values.append(float(m))
        except ValueError:
            continue
    return values


def _split_gmx_command(gmx_cmd: str) -> List[str]:
    """
    Split a possibly wrapped gmx command string into argv tokens.

    Examples:
      "gmx" -> ["gmx"]
      "mpirun -np 4 gmx_mpi" -> ["mpirun", "-np", "4", "gmx_mpi"]
    """
    base = shlex.split(str(gmx_cmd))
    if not base:
        raise ValueError(f"Invalid empty gmx command: {gmx_cmd!r}")
    return base


def _estimate_bin_width_from_r(r_values: List[float]) -> Optional[float]:
    """Estimate RDF bin width from monotonic r-values."""
    if len(r_values) < 2:
        return None
    diffs: List[float] = []
    for i in range(1, len(r_values)):
        dr = r_values[i] - r_values[i - 1]
        if dr > 0:
            diffs.append(dr)
    if not diffs:
        return None
    # Median-like robust choice from sorted positive deltas.
    diffs.sort()
    mid = len(diffs) // 2
    if len(diffs) % 2 == 1:
        return diffs[mid]
    return 0.5 * (diffs[mid - 1] + diffs[mid])


def _nm_to_point_count(length_nm: float, bin_width_nm: Optional[float]) -> int:
    """Convert a physical window length to an integer point count."""
    if bin_width_nm is None or bin_width_nm <= 0:
        return 1
    if length_nm <= 0:
        return 1
    return max(1, int(round(length_nm / bin_width_nm)))


def _moving_average(values: List[float], window: int) -> List[float]:
    """
    Apply fixed-width symmetric moving average smoothing.

    Uses reflect-style index folding so every output point uses exactly `window`
    samples, including near the array boundaries.
    """
    if window <= 1 or len(values) < window:
        return values

    n = len(values)
    half_w = window // 2
    result: List[float] = []

    def _reflect_index(idx: int) -> int:
        # Fold indices into [0, n-1] by reflection, preserving symmetry.
        while idx < 0 or idx >= n:
            if idx < 0:
                idx = -idx - 1
            else:
                idx = (2 * n) - idx - 1
        return idx

    for i in range(n):
        acc = 0.0
        for off in range(-half_w, half_w + 1):
            acc += values[_reflect_index(i + off)]
        result.append(acc / float(window))
    return result


def _dot3(u: List[float], v: List[float]) -> float:
    """3D dot product."""
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2])


def _cross3(u: List[float], v: List[float]) -> List[float]:
    """3D cross product."""
    return [
        (u[1] * v[2]) - (u[2] * v[1]),
        (u[2] * v[0]) - (u[0] * v[2]),
        (u[0] * v[1]) - (u[1] * v[0]),
    ]


def _norm3(v: List[float]) -> float:
    """3D Euclidean norm."""
    return math.sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]))


def _triclinic_min_altitude(a: List[float], b: List[float], c: List[float]) -> Optional[float]:
    """
    Compute the minimum triclinic box altitude (distance between opposite faces).

    h_a = V / ||b x c||, h_b = V / ||c x a||, h_c = V / ||a x b||
    where V = |a . (b x c)|.
    """
    bc = _cross3(b, c)
    ca = _cross3(c, a)
    ab = _cross3(a, b)
    n_bc = _norm3(bc)
    n_ca = _norm3(ca)
    n_ab = _norm3(ab)
    if n_bc <= 0.0 or n_ca <= 0.0 or n_ab <= 0.0:
        return None

    volume = abs(_dot3(a, bc))
    if volume <= 0.0 or not math.isfinite(volume):
        return None

    h_a = volume / n_bc
    h_b = volume / n_ca
    h_c = volume / n_ab
    h_min = min(h_a, h_b, h_c)
    if h_min <= 0.0 or not math.isfinite(h_min):
        return None
    return h_min


def parse_ndx_groups(ndx_path: Path) -> Tuple[List[str], set]:
    """
    Parse .ndx file and return (ordered list, set) of group names.
    """
    groups: List[str] = []
    with open(ndx_path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith("[") and line.endswith("]"):
                name = line[1:-1].strip()
                if name:
                    groups.append(name)
    return groups, set(groups)


def _resolve_index_path(ctx: "PipelineContext") -> Tuple[Optional[Path], List[str]]:
    """
    Determine index.ndx path following deterministic priority.
    Returns (path or None, warnings).
    """
    warnings: List[str] = []
    # (a) manifest-provided path (if exists)
    if ctx.manifest is not None:
        for key in ("index_ndx_path", "gromacs_index_ndx", "index_path", "ndx_path"):
            val = ctx.manifest.get(key)
            if val:
                candidate = Path(val)
                if candidate.exists():
                    try:
                        if ctx.path_manager is not None:
                            candidate = ctx.path_manager.validate_input_path(candidate)
                        return candidate, warnings
                    except Exception as e:
                        warnings.append(f"Manifest index path not under IN/: {candidate} ({e})")
                else:
                    warnings.append(f"Manifest index path not found: {candidate}")

    # (b) default index.ndx under IN
    ndx_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "ndx")
    default_ndx = ndx_dir / "index.ndx"
    if default_ndx.exists():
        return default_ndx, warnings

    # (c) single *.ndx under ndx_dir
    if ndx_dir.exists():
        ndx_candidates = sorted(ndx_dir.glob("*.ndx"))
        if len(ndx_candidates) == 1:
            return ndx_candidates[0], warnings
        if len(ndx_candidates) > 1:
            msg = "Multiple .ndx files found under IN; please keep exactly one or rename to index.ndx:\n"
            msg += "\n".join(f"  - {p}" for p in ndx_candidates)
            raise ValueError(msg)

    return None, warnings


def _looks_like_expression(s: str) -> bool:
    return _is_selection_expression(s)


def _normalize_group_spec(
    raw_value: Any,
    field_name: str,
    pair_name: str,
    warnings: List[str],
) -> Dict[str, str]:
    """
    Normalize group spec into dict with keys: mode, value.
    Accepts dict (canonical) or string (legacy).
    """
    if isinstance(raw_value, dict):
        mode = str(raw_value.get("mode", "")).strip().lower()
        value = str(raw_value.get("value", "")).strip()
        if mode not in {"group", "expr"} or not value:
            raise ValueError(
                f"Invalid {field_name} spec for pair '{pair_name}': expected mode=group|expr and value"
            )
        return {"mode": mode, "value": value, "source": "object"}

    if isinstance(raw_value, str):
        value = raw_value.strip()
        if not value:
            raise ValueError(f"Empty {field_name} value for pair '{pair_name}'")
        warnings.append(
            f"[DEPRECATION] Pair '{pair_name}' uses string {field_name}. "
            f"Prefer {field_name}: {{mode: group|expr, value: \"...\"}}"
        )
        return {"mode": "auto", "value": value, "source": "string"}

    raise ValueError(f"Invalid {field_name} for pair '{pair_name}': expected mapping or string")


def _resolve_selection(
    spec: Dict[str, str],
    index_groups: Optional[set],
    index_path: Optional[Path],
    strict: bool,
    field_name: str,
    pair_name: str,
) -> Tuple[str, str, Optional[str]]:
    """
    Resolve a group/expr spec into a gmx selection string.
    Returns (selection_string, resolved_mode, group_name_if_group).
    """
    mode = spec["mode"]
    value = spec["value"]

    if mode == "group":
        if index_path is None or index_groups is None:
            raise ValueError(
                f"index.ndx is required for group selection {field_name}='{value}' "
                f"(pair '{pair_name}'). Place index.ndx under IN/systems/<SYSTEM_ID>/gromacs/ndx/."
            )
        if value not in index_groups:
            suggestions = difflib.get_close_matches(value, sorted(index_groups), n=5)
            hint = f" Closest matches: {', '.join(suggestions)}" if suggestions else ""
            raise ValueError(
                f"Group '{value}' not found in index {index_path}.{hint}"
            )
        return _to_group_selection(value), "group", value

    if mode == "expr":
        return value, "expr", None

    # auto mode (legacy strings)
    if index_groups is not None and value in index_groups:
        return _to_group_selection(value), "group", value

    if _looks_like_expression(value):
        # Expression provided explicitly (no wrapping)
        return value, "expr", None

    # Ambiguous -> fail fast
    if strict:
        raise ValueError(
            f"Ambiguous selection '{value}' for {field_name} in pair '{pair_name}': "
            "not an index group and not a clear expression. "
            "Use explicit {mode: group|expr} or add group to index.ndx."
        )

    raise ValueError(
        f"Ambiguous selection '{value}' for {field_name} in pair '{pair_name}' "
        "(strict=false still requires explicit mode)."
    )


def find_first_minimum_details(
    xvg_path: Path,
    search_start_nm: float = 0.2,
    smooth_window: int = 5,
    peak_height_threshold: float = _DEFAULT_PEAK_HEIGHT_THRESHOLD,
    min_sep_nm: float = 0.05,
    min_depth: float = _DEFAULT_MIN_DEPTH,
    baseline_window_nm: float = _DEFAULT_BASELINE_WINDOW_NM,
    truncation_tail_points: int = _DEFAULT_TRUNCATION_TAIL_POINTS,
    recovery_window_nm: float = _DEFAULT_RECOVERY_WINDOW_NM,
    recovery_height: float = _DEFAULT_RECOVERY_HEIGHT,
    recovery_positive_points: int = _DEFAULT_RECOVERY_POSITIVE_POINTS,
    local_window_nm: float = _DEFAULT_LOCAL_WINDOW_NM,
    max_sep_nm: float = _DEFAULT_MAX_SEP_NM,
    refine_points: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Find first-shell cutoff as the first minimum after the first significant peak.

    Returns diagnostics:
    - cutoff_r (None if not found)
    - peak_r / peak_g (if found)
    - baseline, baseline_unreliable
    - warnings list and machine-readable reason
    """
    details: Dict[str, Any] = {
        "cutoff_method": "first_minimum_smoothed",
        "cutoff_r": None,
        "peak_r": None,
        "peak_g": None,
        "min_idx_smooth": None,
        "min_idx_raw": None,
        "g_rmin_smooth": None,
        "g_rmin_raw": None,
        "recovery_passed": False,
        "fallback_used": None,
        "baseline": None,
        "baseline_unreliable": False,
        "reliability_flags": [],
        "search_start_nm": float(search_start_nm),
        "smooth_window_requested": int(smooth_window),
        "smooth_window_used": None,
        "peak_height_threshold_requested": float(peak_height_threshold),
        "peak_height_threshold_used": None,
        "min_sep_nm": float(min_sep_nm),
        "min_depth": float(min_depth),
        "baseline_window_nm": float(baseline_window_nm),
        "truncation_tail_points": int(truncation_tail_points),
        "recovery_window_nm": float(recovery_window_nm),
        "recovery_height": float(recovery_height),
        "recovery_positive_points": int(recovery_positive_points),
        "local_window_nm": float(local_window_nm),
        "max_sep_nm": float(max_sep_nm),
        "refine_points_requested": None if refine_points is None else int(refine_points),
        "refine_points_used": None,
        "refine_mode": "none",
        "warnings": [],
        "reason": None,
    }

    if not xvg_path.exists():
        details["reason"] = "rdf_missing"
        details["warnings"].append("rdf_file_missing")
        details["reliability_flags"].append("rdf_missing")
        return details

    r_values: List[float] = []
    g_values: List[float] = []

    with open(xvg_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(("#", "@")) or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    r = float(parts[0])
                    g = float(parts[1])
                    r_values.append(r)
                    g_values.append(g)
                except ValueError:
                    continue

    if len(r_values) < 3:
        details["reason"] = "insufficient_points"
        details["warnings"].append("rdf_data_too_short")
        details["reliability_flags"].append("rdf_data_too_short")
        return details

    # Validate smoothing window: odd and <= number of points
    window = int(smooth_window)
    if window < 1:
        details["warnings"].append(f"smooth_window_invalid={window}; using 1")
        window = 1
    if window % 2 == 0:
        details["warnings"].append(f"smooth_window_even={window}; using {window + 1}")
        window += 1
    if window > len(g_values):
        clamped = len(g_values) if len(g_values) % 2 == 1 else max(1, len(g_values) - 1)
        details["warnings"].append(f"smooth_window_clamped={window}->{clamped}")
        window = clamped
    details["smooth_window_used"] = window

    # Apply smoothing
    g_smooth = _moving_average(g_values, window)
    bin_width_est = _estimate_bin_width_from_r(r_values)

    # Refinement is explicit opt-in. Default is no refinement.
    refine_points_value = 0 if refine_points is None else int(refine_points)
    if refine_points_value < 0:
        details["warnings"].append(f"refine_points_invalid={refine_points_value}; using 0")
        refine_points_value = 0
    details["refine_points_used"] = refine_points_value
    details["refine_mode"] = "smooth_only" if refine_points_value > 0 else "none"

    recovery_window_nm_value = max(0.0, float(recovery_window_nm))
    if recovery_window_nm_value != float(recovery_window_nm):
        details["warnings"].append(
            f"recovery_window_nm_invalid={recovery_window_nm}; using {recovery_window_nm_value}"
        )
    details["recovery_window_nm"] = recovery_window_nm_value
    recovery_window_points = _nm_to_point_count(recovery_window_nm_value, bin_width_est)

    recovery_height_value = max(0.0, float(recovery_height))
    if recovery_height_value != float(recovery_height):
        details["warnings"].append(f"recovery_height_invalid={recovery_height}; using {recovery_height_value}")
    details["recovery_height"] = recovery_height_value

    recovery_positive_points_value = int(recovery_positive_points)
    if recovery_positive_points_value < 1:
        details["warnings"].append(
            f"recovery_positive_points_invalid={recovery_positive_points_value}; using 1"
        )
        recovery_positive_points_value = 1
    details["recovery_positive_points"] = recovery_positive_points_value

    local_window_nm_value = max(0.0, float(local_window_nm))
    if local_window_nm_value != float(local_window_nm):
        details["warnings"].append(f"local_window_nm_invalid={local_window_nm}; using {local_window_nm_value}")
    details["local_window_nm"] = local_window_nm_value
    local_window_points = _nm_to_point_count(local_window_nm_value, bin_width_est)

    max_sep_nm_value = float(max_sep_nm)
    if max_sep_nm_value <= 0.0:
        details["warnings"].append(f"max_sep_nm_invalid={max_sep_nm_value}; disabling max-separation cap")
        max_sep_nm_value = 0.0
    details["max_sep_nm"] = max_sep_nm_value

    # Baseline around search_start (prefer centered window, then right-side fallback)
    baseline_half = max(0.0, float(baseline_window_nm)) / 2.0
    baseline_vals = [
        g for r, g in zip(r_values, g_smooth)
        if (search_start_nm - baseline_half) <= r <= (search_start_nm + baseline_half)
    ]
    if not baseline_vals:
        baseline_vals = [
            g for r, g in zip(r_values, g_smooth)
            if search_start_nm <= r <= (search_start_nm + max(0.0, float(baseline_window_nm)))
        ]

    if baseline_vals:
        baseline = sum(baseline_vals) / len(baseline_vals)
    else:
        baseline = 1.0
        details["baseline_unreliable"] = True
        details["warnings"].append("baseline_window_empty_using_1.0")
        details["reliability_flags"].append("baseline_unreliable")
    details["baseline"] = baseline

    threshold = float(peak_height_threshold)
    if threshold <= 0:
        details["warnings"].append(f"peak_height_threshold_invalid={threshold}; using 0.0")
        threshold = 0.0
    elif threshold > 1.0:
        # Backward-compatible interpretation for legacy factor-based settings (e.g. 1.2x baseline)
        converted = max(0.0, baseline * (threshold - 1.0))
        details["warnings"].append(
            f"peak_height_threshold_legacy_factor={threshold}; converted_to_additive={converted:.6f}"
        )
        threshold = converted
    details["peak_height_threshold_used"] = threshold

    min_depth_value = max(0.0, float(min_depth))
    if min_depth_value != float(min_depth):
        details["warnings"].append(f"min_depth_invalid={min_depth}; using {min_depth_value}")
    details["min_depth"] = min_depth_value

    # Find first significant peak (local maximum above baseline + threshold)
    peak_idx = None
    for i in range(1, len(r_values) - 1):
        if r_values[i] < search_start_nm:
            continue
        if g_smooth[i] >= g_smooth[i - 1] and g_smooth[i] >= g_smooth[i + 1]:
            if (g_smooth[i] - baseline) >= threshold:
                peak_idx = i
                details["peak_r"] = r_values[i]
                details["peak_g"] = g_smooth[i]
                break

    if peak_idx is None:
        details["reason"] = "no_first_peak"
        details["warnings"].append("no_first_peak")
        details["reliability_flags"].append("no_first_peak")
        return details

    peak_r = float(details["peak_r"])
    peak_g = float(details["peak_g"])
    min_search_r = peak_r + float(min_sep_nm)
    max_search_r = r_values[-1]
    if max_sep_nm_value > 0.0:
        max_search_r = min(max_search_r, peak_r + max_sep_nm_value)
    if max_search_r <= min_search_r:
        details["reason"] = "search_window_empty"
        details["reliability_flags"].append("search_window_empty")
        details["warnings"].append(
            f"search_window_empty_after_peak: min_r={min_search_r:.6f}, max_r={max_search_r:.6f}"
        )
        return details

    selected_min_idx_smooth: Optional[int] = None
    selected_recovery_passed = False
    local_min_found_in_range = False

    # Find first physically meaningful minimum after the first peak.
    for i in range(peak_idx + 1, len(r_values) - 1):
        r_i = r_values[i]
        if r_i < min_search_r:
            continue
        if r_i > max_search_r:
            break

        is_local_min = (
            g_smooth[i] <= g_smooth[i - 1]
            and g_smooth[i] <= g_smooth[i + 1]
            and (g_smooth[i] < g_smooth[i - 1] or g_smooth[i] < g_smooth[i + 1])
        )
        if not is_local_min:
            continue
        local_min_found_in_range = True

        depth = peak_g - g_smooth[i]
        if depth < min_depth_value:
            continue

        # Recovery test (A): RDF should rise after the candidate minimum.
        recovery_a = False
        rise_end = min(len(g_smooth) - 1, i + recovery_window_points)
        if rise_end > i:
            local_rise = max(g_smooth[j] - g_smooth[i] for j in range(i + 1, rise_end + 1))
            recovery_a = local_rise >= recovery_height_value

        # Recovery test (B): positive derivative for K consecutive bins.
        recovery_b = False
        if i + recovery_positive_points_value < len(g_smooth):
            recovery_b = True
            for j in range(i, i + recovery_positive_points_value):
                if (g_smooth[j + 1] - g_smooth[j]) <= 0.0:
                    recovery_b = False
                    break

        recovery_passed = recovery_a or recovery_b

        # Optional basin robustness: candidate must be minimum in local neighborhood.
        basin_ok = True
        if local_window_nm_value > 0.0:
            left = max(0, i - local_window_points)
            right = min(len(g_smooth) - 1, i + local_window_points)
            basin_min = min(g_smooth[left:right + 1])
            basin_ok = g_smooth[i] <= (basin_min + 1e-12)
        if not basin_ok:
            continue

        selected_min_idx_smooth = i
        selected_recovery_passed = recovery_passed
        break

    if selected_min_idx_smooth is None:
        tail_count = max(3, int(truncation_tail_points))
        tail_start = max(peak_idx + 1, len(g_smooth) - tail_count)
        tail = g_smooth[tail_start:]
        # Treat only strictly descending tails as truncation; flat shoulders are handled separately.
        if len(tail) >= 3 and all(tail[j] > tail[j + 1] for j in range(len(tail) - 1)):
            details["reason"] = "truncated_before_first_minimum"
            details["reliability_flags"].append("truncated_before_first_minimum")
            details["warnings"].append(
                "truncated_before_first_minimum; increase rdf.max_r/search_rmax to capture the first shell minimum"
            )
            return details

        details["reason"] = "no_first_minimum"
        details["reliability_flags"].append("no_first_minimum")
        if not local_min_found_in_range:
            details["reliability_flags"].append("shoulder_no_minimum")
            details["warnings"].append("shoulder_no_minimum")
        details["warnings"].append("no_first_minimum_after_peak")
        return details

    # Use the smoothed minimum index for stability. Refinement is smooth-only (no raw argmin).
    min_idx_raw = selected_min_idx_smooth

    details["min_idx_smooth"] = selected_min_idx_smooth
    details["min_idx_raw"] = min_idx_raw
    details["g_rmin_smooth"] = g_smooth[selected_min_idx_smooth]
    details["g_rmin_raw"] = g_values[min_idx_raw]
    details["recovery_passed"] = selected_recovery_passed
    if not selected_recovery_passed:
        details["reliability_flags"].append("no_recovery_after_minimum")
        details["warnings"].append("no_recovery_after_minimum")
    details["cutoff_r"] = r_values[min_idx_raw]
    details["reason"] = "ok"
    return details


def find_first_minimum(
    xvg_path: Path,
    search_start_nm: float = 0.2,
    smooth_window: int = 5,
    peak_threshold: float = _DEFAULT_PEAK_HEIGHT_THRESHOLD,
    min_sep_nm: float = 0.05,
    min_depth: float = _DEFAULT_MIN_DEPTH,
    baseline_window_nm: float = _DEFAULT_BASELINE_WINDOW_NM,
    truncation_tail_points: int = _DEFAULT_TRUNCATION_TAIL_POINTS,
    recovery_window_nm: float = _DEFAULT_RECOVERY_WINDOW_NM,
    recovery_height: float = _DEFAULT_RECOVERY_HEIGHT,
    recovery_positive_points: int = _DEFAULT_RECOVERY_POSITIVE_POINTS,
    local_window_nm: float = _DEFAULT_LOCAL_WINDOW_NM,
    max_sep_nm: float = _DEFAULT_MAX_SEP_NM,
    refine_points: Optional[int] = None,
) -> Optional[float]:
    """
    Backward-compatible wrapper: returns only cutoff value.
    """
    details = find_first_minimum_details(
        xvg_path,
        search_start_nm=search_start_nm,
        smooth_window=smooth_window,
        peak_height_threshold=peak_threshold,
        min_sep_nm=min_sep_nm,
        min_depth=min_depth,
        baseline_window_nm=baseline_window_nm,
        truncation_tail_points=truncation_tail_points,
        recovery_window_nm=recovery_window_nm,
        recovery_height=recovery_height,
        recovery_positive_points=recovery_positive_points,
        local_window_nm=local_window_nm,
        max_sep_nm=max_sep_nm,
        refine_points=refine_points,
    )
    return details.get("cutoff_r")


def parse_cn_from_xvg(xvg_path: Path, cutoff_nm: float) -> Optional[float]:
    """
    Parse CN(r) curve from gmx rdf -cn output and extract CN at cutoff.
    
    The cn.xvg contains cumulative CN as a function of r.
    We find the value at or just before the cutoff.
    
    Args:
        xvg_path: Path to CN .xvg file from gmx rdf -cn
        cutoff_nm: Integration cutoff (nm)
        
    Returns:
        Coordination number at cutoff, or None if parsing fails
    """
    if not xvg_path.exists():
        return None
    
    r_values: List[float] = []
    cn_values: List[float] = []
    
    with open(xvg_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(("#", "@")) or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    r = float(parts[0])
                    cn = float(parts[1])
                    r_values.append(r)
                    cn_values.append(cn)
                except ValueError:
                    continue
    
    if not r_values:
        return None
    
    # Find CN at cutoff (interpolate if needed)
    for i, r in enumerate(r_values):
        if r >= cutoff_nm:
            if i == 0:
                return cn_values[0]
            # Linear interpolation
            r_prev, r_curr = r_values[i-1], r
            cn_prev, cn_curr = cn_values[i-1], cn_values[i]
            if r_curr == r_prev:
                return cn_curr
            frac = (cutoff_nm - r_prev) / (r_curr - r_prev)
            return cn_prev + frac * (cn_curr - cn_prev)
    
    # Cutoff beyond data range
    return cn_values[-1] if cn_values else None


def get_box_dimensions(
    tpr_path: Path,
    gmx_cmd: str,
    timeout_s: int = _DEFAULT_GMX_TIMEOUT_S,
) -> Optional[Tuple[float, float, float]]:
    """
    Get box dimensions from TPR file using gmx check.
    
    Args:
        tpr_path: Path to .tpr file
        
    Returns:
        (Lx, Ly, Lz) in nm, or None if extraction fails
    """
    try:
        base = _split_gmx_command(gmx_cmd)
        cmd = base + ["check", "-s", str(tpr_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )

        if result.returncode != 0:
            return None

        lines = result.stderr.split("\n") + result.stdout.split("\n")
        for i, line in enumerate(lines):
            line_lower = line.lower().strip()
            if not re.match(r"^box\b", line_lower):
                continue

            # Handle triclinic (3x3) matrix
            if "3x3" in line_lower:
                floats: List[float] = []
                if ":" in line:
                    floats.extend(_extract_floats(line.split(":", 1)[1]))
                j = i + 1
                while len(floats) < 9 and j < len(lines):
                    floats.extend(_extract_floats(lines[j]))
                    j += 1
                if len(floats) >= 9:
                    lx, ly, lz = floats[0], floats[4], floats[8]
                else:
                    continue
            else:
                payload = line.split(":", 1)[1] if ":" in line else line
                floats = _extract_floats(payload)
                if len(floats) < 3:
                    continue
                lx, ly, lz = floats[0], floats[1], floats[2]

            if lx > 0.0 and ly > 0.0 and lz > 0.0:
                return (lx, ly, lz)

        return None

    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return None


def get_min_box_height_from_tpr(
    tpr_path: Path,
    gmx_cmd: str,
    timeout_s: int = _DEFAULT_GMX_TIMEOUT_S,
) -> Optional[float]:
    """
    Extract minimum safe box repeat distance from TPR via `gmx check -s`.

    - Triclinic: uses minimum cell altitude from full 3x3 box matrix
    - Orthorhombic: uses min(Lx, Ly, Lz)
    """
    try:
        base = _split_gmx_command(gmx_cmd)
        cmd = base + ["check", "-s", str(tpr_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )

        if result.returncode != 0:
            return None

        lines = result.stderr.split("\n") + result.stdout.split("\n")
        for i, line in enumerate(lines):
            line_lower = line.lower().strip()
            if not re.match(r"^box\b", line_lower):
                continue

            if "3x3" in line_lower:
                floats: List[float] = []
                payload = line.split(":", 1)[1] if ":" in line else ""
                floats.extend(_extract_floats(payload))
                j = i + 1
                while len(floats) < 9 and j < len(lines):
                    floats.extend(_extract_floats(lines[j]))
                    j += 1
                if len(floats) < 9:
                    continue
                a = floats[0:3]
                b = floats[3:6]
                c = floats[6:9]
                h_min = _triclinic_min_altitude(a, b, c)
                if h_min is not None:
                    return h_min
                continue

            payload = line.split(":", 1)[1] if ":" in line else line
            floats = _extract_floats(payload)
            if len(floats) < 3:
                continue
            lx, ly, lz = floats[0], floats[1], floats[2]
            if lx > 0.0 and ly > 0.0 and lz > 0.0:
                return min(lx, ly, lz)

        return None

    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return None


def get_min_box_length_from_traj(
    xtc_path: Path,
    tpr_path: Path,
    gmx_cmd: str,
    timeout_s: int = _DEFAULT_GMX_TIMEOUT_S,
    dt_ps: Optional[float] = None,
) -> Optional[float]:
    """
    Extract minimum safe box repeat distance from trajectory via `gmx traj -ob`.

    Supports `-ob` output formats:
    - time + 3 cols: Lx, Ly, Lz
    - time + 9 cols: triclinic vectors (ax ay az bx by bz cx cy cz)
    """
    tmp_path: Optional[Path] = None
    try:
        base = _split_gmx_command(gmx_cmd)
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".xvg",
            prefix="traj_box_",
            delete=False,
        ) as tf:
            tmp_path = Path(tf.name)

        cmd = base + [
            "traj",
            "-f",
            str(xtc_path),
            "-s",
            str(tpr_path),
            "-ob",
            str(tmp_path),
            "-xvg",
            "none",
        ]
        if dt_ps is not None and float(dt_ps) > 0.0:
            cmd.extend(["-dt", str(float(dt_ps))])

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            input="0\n",
            timeout=timeout_s,
        )
        if result.returncode != 0:
            print(f"  [WARN] Failed to extract trajectory box heights with gmx traj (exit={result.returncode})")
            return None
        if tmp_path is None or not tmp_path.exists():
            print("  [WARN] gmx traj did not produce -ob output file")
            return None

        min_box_height: Optional[float] = None
        with open(tmp_path, "r") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith(("#", "@")):
                    continue
                parts = line.split()
                try:
                    vals = [float(x) for x in parts]
                except ValueError:
                    continue

                frame_min: Optional[float] = None
                if len(vals) == 4:
                    frame_min = min(vals[1], vals[2], vals[3])
                elif len(vals) == 10:
                    a = vals[1:4]
                    b = vals[4:7]
                    c = vals[7:10]
                    frame_min = _triclinic_min_altitude(a, b, c)

                if frame_min is None:
                    continue
                if frame_min <= 0.0:
                    continue
                if min_box_height is None or frame_min < min_box_height:
                    min_box_height = frame_min

        if min_box_height is None:
            print("  [WARN] Could not parse valid trajectory box heights from gmx traj -ob output")
        return min_box_height

    except (subprocess.TimeoutExpired, FileNotFoundError, ValueError, Exception) as e:
        print(f"  [WARN] Trajectory box extraction failed: {e}")
        return None
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass


class AnalysisStage(BaseStage):
    """
    Stage 4: Analysis (RDF/CN)
    
    Reads:
    - IN/systems/<SYSTEM_ID>/gromacs/analysis_pairs.yaml - pair definitions
    - IN/systems/<SYSTEM_ID>/gromacs/ndx/index.ndx (if present)
    - OUT_GMX/<RUN_ID>/03_gromacs/md/ - production trajectory
    
    Writes:
    - OUT_GMX/<RUN_ID>/analysis/rdf/ - radial distribution functions
    - OUT_GMX/<RUN_ID>/analysis/cn/ - coordination numbers (via gmx rdf -cn)
    
    Key improvements:
    - CN via GROMACS native -cn option (not hardcoded water density)
    - Explicit group/expr handling with index support
    - rmax capped using 0.5*min box altitude for periodic artifact prevention
    - Robust first_minimum with smoothing
    """
    
    @property
    def name(self) -> str:
        return "analysis"
    
    @property
    def output_subdir(self) -> str:
        return "analysis"
    
    def run(self, ctx: "PipelineContext") -> bool:
        """Execute analysis stage."""
        print(f"  Running Analysis stage...")
        
        # Parse analysis config
        config_path = ctx.get_input_path(
            "systems", ctx.system_id, "gromacs", "analysis_pairs.yaml"
        )
        
        try:
            config = parse_analysis_pairs(config_path)
            print(f"  - Config: {config_path.name}")
        except FileNotFoundError as e:
            print(f"  [ERROR] {e}")
            return False
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False

        analysis_settings = config.get("analysis_settings", {}) if isinstance(config, dict) else {}
        if not isinstance(analysis_settings, dict):
            print(f"  [ERROR] analysis_settings must be a mapping")
            return False
        strict_analysis = _parse_bool(analysis_settings.get("strict"), _DEFAULT_STRICT_ANALYSIS)
        try:
            gmx_timeout_s = int(analysis_settings.get("gmx_timeout_s", _DEFAULT_GMX_TIMEOUT_S))
        except (TypeError, ValueError):
            print(f"  [ERROR] analysis_settings.gmx_timeout_s must be an integer")
            return False
        if gmx_timeout_s <= 0:
            print(f"  [ERROR] analysis_settings.gmx_timeout_s must be > 0")
            return False
        traj_box_dt_ps_raw = analysis_settings.get("traj_box_dt_ps")
        traj_box_dt_ps: Optional[float] = None
        if traj_box_dt_ps_raw is not None:
            try:
                traj_box_dt_ps = float(traj_box_dt_ps_raw)
            except (TypeError, ValueError):
                print(f"  [ERROR] analysis_settings.traj_box_dt_ps must be numeric when provided")
                return False
            if traj_box_dt_ps <= 0.0:
                print(f"  [ERROR] analysis_settings.traj_box_dt_ps must be > 0 when provided")
                return False
        
        # Verify trajectory exists (deterministic preference)
        md_dir = ctx.get_output_path("03_gromacs", "md")
        md_xtc = md_dir / "md.xtc"
        md_tpr = md_dir / "md.tpr"
        
        # Also check for .trr if .xtc not found
        if not md_xtc.exists():
            md_xtc = md_dir / "md.trr"
        
        if not md_xtc.exists():
            print(f"  [ERROR] No trajectory found in {md_dir}")
            print(f"          Expected {md_dir / 'md.xtc'} or {md_dir / 'md.trr'}")
            print(f"          Run gmx_prod first to generate a production trajectory.")
            return False
        
        if not md_tpr.exists():
            print(f"  [ERROR] TPR file not found: {md_tpr}")
            return False

        print(f"  - Trajectory: {md_xtc.name}")
        print(f"  - TPR: {md_tpr.name}")
        print(f"  - GROMACS timeout: {gmx_timeout_s}s")
        
        # Resolve index.ndx (if available) and parse groups
        try:
            index_path, index_warnings = _resolve_index_path(ctx)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False

        for w in index_warnings:
            print(f"  [WARN] {w}")

        index_groups_list: List[str] = []
        index_groups_set: Optional[set] = None
        if index_path is not None:
            try:
                index_groups_list, index_groups_set = parse_ndx_groups(index_path)
                print(f"  - Index: {index_path.name} ({len(index_groups_list)} groups)")
            except Exception as e:
                print(f"  [ERROR] Failed to parse index.ndx: {e}")
                return False
        else:
            print(f"  - Index: (none)")

        gmx_cmd = resolve_gmx_command(ctx)

        # Get box dimensions for conservative rmax capping.
        box_dims_tpr = get_box_dimensions(md_tpr, gmx_cmd, timeout_s=gmx_timeout_s)
        min_box_height_tpr = get_min_box_height_from_tpr(md_tpr, gmx_cmd, timeout_s=gmx_timeout_s)
        if box_dims_tpr:
            print(f"  - TPR box: {box_dims_tpr[0]:.3f} x {box_dims_tpr[1]:.3f} x {box_dims_tpr[2]:.3f} nm")
        else:
            print(f"  [WARN] Could not determine TPR box dimensions")
        if min_box_height_tpr is not None:
            print(f"  - TPR min box height (altitude): {min_box_height_tpr:.3f} nm")
        else:
            print(f"  [WARN] Could not determine TPR min box height")

        min_box_height_traj = get_min_box_length_from_traj(
            md_xtc,
            md_tpr,
            gmx_cmd,
            timeout_s=gmx_timeout_s,
            dt_ps=traj_box_dt_ps,
        )
        if min_box_height_traj is not None:
            dt_info = f", dt={traj_box_dt_ps:g} ps" if traj_box_dt_ps is not None else ""
            print(f"  - Trajectory min box height: {min_box_height_traj:.3f} nm{dt_info}")
        else:
            print(f"  [WARN] Could not determine trajectory min box height; using TPR-only cap if available")
        
        # Create output directories
        output_dir = self.get_output_dir(ctx)
        rdf_dir = output_dir / "rdf"
        cn_dir = output_dir / "cn"
        rdf_dir.mkdir(parents=True, exist_ok=True)
        cn_dir.mkdir(parents=True, exist_ok=True)
        
        # Process pairs
        pairs = config.get("pairs", [])
        if not isinstance(pairs, list):
            print(f"  [ERROR] analysis_pairs.yaml 'pairs' must be a list")
            return False

        rdf_settings = config.get("rdf_settings", {})
        cn_settings = config.get("cn_settings", {})
        if not isinstance(rdf_settings, dict) or not isinstance(cn_settings, dict):
            print(f"  [ERROR] rdf_settings/cn_settings must be mappings")
            return False

        warnings: List[str] = []
        
        analysis_results = []

        # Optional audit warning: readable prefix collisions (hash suffix still guarantees unique outputs)
        prefix_to_names: Dict[str, List[str]] = defaultdict(list)
        for pair_entry in pairs:
            if isinstance(pair_entry, dict):
                pair_name = str(pair_entry.get("name", "unnamed"))
                safe_prefix = _safe_filename(
                    pair_name,
                    max_len=max(8, _MAX_FILENAME_LEN - (_HASH_SUFFIX_LEN + 2)),
                )
                prefix_to_names[safe_prefix].append(pair_name)
        for prefix, names in sorted(prefix_to_names.items()):
            unique_names = sorted(set(names))
            if len(unique_names) > 1:
                msg = (
                    f"Output prefix collision for '{prefix}': {unique_names}; "
                    "hash suffix preserves unique filenames"
                )
                warnings.append(msg)
        
        for pair in pairs:
            if not isinstance(pair, dict):
                print(f"  [ERROR] Invalid pair entry (expected mapping): {pair!r}")
                return False

            name = str(pair.get("name", "unnamed"))
            description = pair.get("description", "")
            pair_warnings: List[str] = []

            # Backward-compat: support A/B
            if "group1" not in pair and "A" in pair:
                pair["group1"] = pair.get("A")
                warnings.append(f"[DEPRECATION] Pair '{name}' uses legacy key 'A'; use 'group1'")
            if "group2" not in pair and "B" in pair:
                pair["group2"] = pair.get("B")
                warnings.append(f"[DEPRECATION] Pair '{name}' uses legacy key 'B'; use 'group2'")

            if "group1" not in pair or "group2" not in pair:
                print(f"  [ERROR] Pair '{name}' missing group1/group2")
                return False

            try:
                group1_spec = _normalize_group_spec(pair.get("group1"), "group1", name, warnings)
                group2_spec = _normalize_group_spec(pair.get("group2"), "group2", name, warnings)
            except ValueError as e:
                print(f"  [ERROR] {e}")
                return False

            # Define labels for display/manifest (original values from spec)
            group1_label = group1_spec["value"]
            group2_label = group2_spec["value"]

            safe_prefix, output_basename = _build_pair_output_basename(name, group1_spec, group2_spec)

            print(f"\n  Analyzing pair: {name}")

            # Per-pair overrides
            pair_rdf_overrides = pair.get("rdf")
            if pair_rdf_overrides is not None and not isinstance(pair_rdf_overrides, dict):
                print(f"  [ERROR] Pair '{name}' rdf must be a mapping")
                return False
            pair_cn_overrides = pair.get("cn")
            if pair_cn_overrides is not None and not isinstance(pair_cn_overrides, dict):
                print(f"  [ERROR] Pair '{name}' cn must be a mapping")
                return False

            has_pair_max_r = isinstance(pair_rdf_overrides, dict) and "max_r" in pair_rdf_overrides
            has_global_max_r = "max_r" in rdf_settings
            if has_pair_max_r:
                rmax_source = "pair"
            elif has_global_max_r:
                rmax_source = "global"
            else:
                rmax_source = "default"

            pair_rdf_settings = dict(rdf_settings)
            pair_cn_settings = dict(cn_settings)
            if isinstance(pair_rdf_overrides, dict):
                pair_rdf_settings.update(pair_rdf_overrides)
            if isinstance(pair_cn_overrides, dict):
                pair_cn_settings.update(pair_cn_overrides)

            # Resolve selections
            try:
                ref_sel, ref_mode, _ = _resolve_selection(
                    group1_spec, index_groups_set, index_path, strict_analysis, "group1", name
                )
                sel_sel, sel_mode, _ = _resolve_selection(
                    group2_spec, index_groups_set, index_path, strict_analysis, "group2", name
                )
            except ValueError as e:
                print(f"  [ERROR] {e}")
                return False

            print(f"    {group1_spec['value']} - {group2_spec['value']}")
            print(f"    Resolved: ref={ref_sel!r} ({ref_mode}), sel={sel_sel!r} ({sel_mode})")
            print(f"    Output basename: {output_basename}")
            
            # Compute effective rmax (capped to half minimum box height)
            try:
                requested_max_r = float(pair_rdf_settings.get("max_r", _DEFAULT_RDF_MAX_R))
                bin_width = float(pair_rdf_settings.get("bin_width", _DEFAULT_RDF_BIN))
            except (TypeError, ValueError):
                print(f"  [ERROR] Invalid rdf settings for pair '{name}'")
                return False
            if bin_width <= 0:
                print(f"  [ERROR] Invalid bin_width for pair '{name}': {bin_width}")
                return False
            if requested_max_r <= 0:
                print(f"  [ERROR] Invalid max_r for pair '{name}': {requested_max_r}")
                return False
            
            effective_max_r = requested_max_r
            rmax_capped = False
            no_box_cap_applied = False
            rmax_margin_nm = 0.02
            safe_rmax_basis = "min_altitude"
            safe_rmax_tpr: Optional[float] = (
                (0.5 * min_box_height_tpr - rmax_margin_nm) if min_box_height_tpr is not None else None
            )
            safe_rmax_traj: Optional[float] = (
                (0.5 * min_box_height_traj - rmax_margin_nm) if min_box_height_traj is not None else None
            )
            safe_rmax_used: Optional[float] = None
            box_source_used: Optional[str] = None

            if safe_rmax_tpr is not None and safe_rmax_traj is not None:
                safe_rmax_used = min(safe_rmax_tpr, safe_rmax_traj)
                box_source_used = "traj" if safe_rmax_traj <= safe_rmax_tpr else "tpr"
            elif safe_rmax_traj is not None:
                safe_rmax_used = safe_rmax_traj
                box_source_used = "traj"
            elif safe_rmax_tpr is not None:
                safe_rmax_used = safe_rmax_tpr
                box_source_used = "tpr"

            if safe_rmax_used is not None:
                if safe_rmax_used <= 0:
                    msg = f"    [ERROR] Box too small for RDF (safe_rmax={safe_rmax_used:.4f} nm)"
                    print(msg)
                    if strict_analysis:
                        return False
                    analysis_results.append({
                        "name": name,
                        "status": "skipped",
                        "pair_skipped_reason": "box_too_small",
                    })
                    continue

                min_meaningful = max(3 * bin_width, 0.05)
                if safe_rmax_used < min_meaningful:
                    msg = (
                        f"    [ERROR] Box too small for meaningful RDF "
                        f"(safe_rmax={safe_rmax_used:.4f} nm < {min_meaningful:.4f} nm)"
                    )
                    print(msg)
                    if strict_analysis:
                        return False
                    analysis_results.append({
                        "name": name,
                        "status": "skipped",
                        "pair_skipped_reason": "box_too_small",
                    })
                    continue

                if requested_max_r > safe_rmax_used:
                    effective_max_r = safe_rmax_used
                    rmax_capped = True
                    print(
                        f"    [WARN] rmax capped: {requested_max_r:.3f} -> "
                        f"{effective_max_r:.3f} nm (box limit, source={box_source_used}, basis={safe_rmax_basis})"
                    )
            else:
                if rmax_source == "default":
                    msg = (
                        "    [ERROR] Box dimensions unknown and max_r not provided in config; "
                        "cannot choose safe rmax."
                    )
                    print(msg)
                    if strict_analysis:
                        return False
                    analysis_results.append({
                        "name": name,
                        "status": "skipped",
                        "pair_skipped_reason": "box_unknown_no_rmax",
                    })
                    continue
                no_box_cap_applied = True
                if strict_analysis:
                    print("    [WARN] Box unknown; using explicit max_r without capping (strict=true)")
                else:
                    print("    [WARN] Box unknown; using max_r without capping (strict=false)")
            print(
                f"    rmax: requested={requested_max_r:.3f} nm, "
                f"used={effective_max_r:.3f} nm (source={rmax_source}, cap_source={box_source_used})"
            )
            
            # Run RDF with CN output
            rdf_output = rdf_dir / f"rdf_{output_basename}.xvg"
            cn_xvg_output = cn_dir / f"cn_{output_basename}.xvg"
            
            rdf_result = self._compute_rdf_with_cn(
                ctx, md_xtc, md_tpr, ref_sel, sel_sel,
                rdf_output, cn_xvg_output,
                bin_width=bin_width,
                max_r=effective_max_r,
                index_path=index_path,
                gmx_cmd=gmx_cmd,
                gmx_timeout_s=gmx_timeout_s,
            )
            
            if not rdf_result["success"]:
                print(
                    f"    [WARN] RDF/CN computation failed for {name}: "
                    f"{rdf_result.get('error', 'unknown')}"
                )
                if rdf_result.get("stderr"):
                    print(f"    stderr: {rdf_result['stderr'][:300]}")
                analysis_results.append({
                    "name": name,
                    "output_basename": output_basename,
                    "group1": group1_label,
                    "group2": group2_label,
                    "resolved_ref": ref_sel,
                    "resolved_sel": sel_sel,
                    "status": "failed",
                    "error": rdf_result.get("error", "unknown"),
                    "warnings": pair_warnings,
                })
                continue
            
            cn_method = rdf_result.get("cn_method", "gmx_rdf_native")
            print(f"    RDF: {rdf_output.name}")
            
            # Determine CN cutoff
            default_cutoff_raw = pair_cn_settings.get("cutoff_nm", _DEFAULT_CN_CUTOFF)
            try:
                default_cutoff = float(default_cutoff_raw)
            except (TypeError, ValueError):
                print(f"  [ERROR] Invalid cn_settings.cutoff_nm for pair '{name}': {default_cutoff_raw!r}")
                return False
            cutoff_config = pair.get("cn_cutoff")
            if cutoff_config is None:
                cutoff_config = pair_cn_settings.get("cutoff", pair_cn_settings.get("cutoff_nm", default_cutoff))
            
            cutoff_mode = "fixed"
            cutoff_unreliable = False
            is_dryrun_pair = bool(rdf_result.get("dry_run"))
            cutoff_details: Dict[str, Any] = {
                "cutoff_method": "fixed_numeric",
                "peak_r": None,
                "peak_g": None,
                "cutoff_r": None,
                "min_idx_smooth": None,
                "min_idx_raw": None,
                "g_rmin_smooth": None,
                "g_rmin_raw": None,
                "recovery_passed": False,
                "fallback_used": None,
                "refine_mode": "none",
                "baseline": None,
                "baseline_unreliable": False,
                "reason": "fixed_numeric",
                "reliability_flags": [],
                "warnings": [],
            }
            
            cutoff_mode_token = str(cutoff_config).strip().lower() if isinstance(cutoff_config, str) else None
            if cutoff_mode_token in {"first_minimum", "first_minimum_smoothed"}:
                cutoff_mode = "first_minimum_smoothed"
                # Use cn_settings for first-shell cutoff parameters
                try:
                    smooth_window_nm_raw = pair_cn_settings.get("smooth_window_nm")
                    if smooth_window_nm_raw is not None:
                        smooth_window_nm = float(smooth_window_nm_raw)
                        smooth_window = _nm_to_point_count(smooth_window_nm, bin_width)
                        if smooth_window % 2 == 0:
                            smooth_window += 1
                        if "smooth_window" in pair_cn_settings:
                            pair_warnings.append("smooth_window_nm_overrides_smooth_window")
                    else:
                        smooth_window = int(pair_cn_settings.get("smooth_window", 5))
                    peak_height_threshold = float(
                        pair_cn_settings.get(
                            "peak_height_threshold",
                            pair_cn_settings.get("peak_threshold", _DEFAULT_PEAK_HEIGHT_THRESHOLD),
                        )
                    )
                    min_sep_nm = float(pair_cn_settings.get("min_sep_nm", 0.05))
                    min_depth = float(pair_cn_settings.get("min_depth", _DEFAULT_MIN_DEPTH))
                    baseline_window_nm = float(
                        pair_cn_settings.get("baseline_window_nm", _DEFAULT_BASELINE_WINDOW_NM)
                    )
                    search_start = float(pair_cn_settings.get("search_start_nm", 0.2))
                    truncation_tail_points = int(
                        pair_cn_settings.get(
                            "truncation_tail_points",
                            _DEFAULT_TRUNCATION_TAIL_POINTS,
                        )
                    )
                    recovery_window_nm = float(
                        pair_cn_settings.get("recovery_window_nm", _DEFAULT_RECOVERY_WINDOW_NM)
                    )
                    recovery_height = float(
                        pair_cn_settings.get("recovery_height", _DEFAULT_RECOVERY_HEIGHT)
                    )
                    recovery_positive_points = int(
                        pair_cn_settings.get(
                            "recovery_positive_points",
                            _DEFAULT_RECOVERY_POSITIVE_POINTS,
                        )
                    )
                    local_window_nm = float(
                        pair_cn_settings.get("local_window_nm", _DEFAULT_LOCAL_WINDOW_NM)
                    )
                    max_sep_nm = float(pair_cn_settings.get("max_sep_nm", _DEFAULT_MAX_SEP_NM))
                    refine_points_raw = pair_cn_settings.get("refine_points")
                    refine_points = int(refine_points_raw) if refine_points_raw is not None else None
                except (TypeError, ValueError):
                    print(f"  [ERROR] Invalid cn_settings for pair '{name}'")
                    return False

                if is_dryrun_pair:
                    cutoff = float(default_cutoff)
                    cutoff_unreliable = True
                    cutoff_details.update({
                        "cutoff_method": "first_minimum_smoothed",
                        "cutoff_r": cutoff,
                        "reason": "dry_run_no_data",
                        "fallback_used": "dry_run_default_cutoff",
                        "reliability_flags": ["dry_run_no_rdf", "fallback_default_cutoff"],
                        "warnings": ["dry_run_no_rdf_for_first_minimum; using default cutoff_nm"],
                    })
                    pair_warnings.extend(cutoff_details["warnings"])
                    print("    [WARN] dry-run: first_minimum detection skipped; using default cutoff")
                else:
                    cutoff_details = find_first_minimum_details(
                        rdf_output,
                        search_start_nm=search_start,
                        smooth_window=smooth_window,
                        peak_height_threshold=peak_height_threshold,
                        min_sep_nm=min_sep_nm,
                        min_depth=min_depth,
                        baseline_window_nm=baseline_window_nm,
                        truncation_tail_points=truncation_tail_points,
                        recovery_window_nm=recovery_window_nm,
                        recovery_height=recovery_height,
                        recovery_positive_points=recovery_positive_points,
                        local_window_nm=local_window_nm,
                        max_sep_nm=max_sep_nm,
                        refine_points=refine_points,
                    )

                    if (
                        cutoff_details.get("reason") == "truncated_before_first_minimum"
                        and rmax_capped
                    ):
                        guidance = (
                            "truncated_before_first_minimum; RDF rmax is box-limited by safe_rmax "
                            "(min-altitude basis). Use explicit cn_cutoff, analyze after longer/stabler NPT, "
                            "or use a larger box/system for RDF analysis."
                        )
                        existing_warnings = list(cutoff_details.get("warnings", []))
                        cutoff_details["warnings"] = [
                            w for w in existing_warnings
                            if "increase rdf.max_r/search_rmax" not in w
                        ]
                        cutoff_details["warnings"].append(guidance)

                    cutoff = cutoff_details.get("cutoff_r")
                    for cw in cutoff_details.get("warnings", []):
                        pair_warnings.append(cw)
                        print(f"    [WARN] {cw}")

                    if cutoff is None:
                        reason = cutoff_details.get("reason", "unknown")
                        msg = f"    [WARN] Could not find first minimum ({reason})"
                        if strict_analysis:
                            print(f"{msg}; strict=true -> failing analysis")
                            return False
                        print(f"{msg}; using default cutoff")
                        cutoff = float(default_cutoff)
                        cutoff_unreliable = True
                        fallback_warning = f"fallback_to_default_cutoff_nm={cutoff:.6f}"
                        cutoff_details["warnings"] = list(cutoff_details.get("warnings", [])) + [
                            fallback_warning
                        ]
                        cutoff_details["reliability_flags"] = list(
                            cutoff_details.get("reliability_flags", [])
                        ) + ["fallback_default_cutoff"]
                        pair_warnings.append(fallback_warning)
                    else:
                        g_rmin_raw = cutoff_details.get("g_rmin_raw")
                        if isinstance(g_rmin_raw, (int, float)):
                            print(f"    Detected first minimum at {cutoff:.3f} nm (g_raw={float(g_rmin_raw):.3f})")
                        else:
                            print(f"    Detected first minimum at {cutoff:.3f} nm")
            else:
                try:
                    cutoff = float(cutoff_config)
                except (TypeError, ValueError):
                    print(f"  [ERROR] Invalid cn_cutoff for pair '{name}': {cutoff_config!r}")
                    return False
                cutoff_details["cutoff_r"] = cutoff

            if cutoff_details.get("baseline_unreliable", False):
                cutoff_unreliable = True
            
            # Get CN from gmx rdf -cn output
            cn_value: Optional[float] = None
            cn_dat_output: Optional[Path] = None
            if is_dryrun_pair:
                pair_warnings.append("dry_run_no_rdf_cn_generated")
                print("    [DRY-RUN] Skipping CN parse/write (no RDF/CN files generated)")
            else:
                if _is_nonempty_file(cn_xvg_output):
                    cn_value = parse_cn_from_xvg(cn_xvg_output, cutoff)
                else:
                    msg = f"CN output missing/empty: {cn_xvg_output}"
                    pair_warnings.append(msg)
                    print(f"    [WARN] {msg}")
                    analysis_results.append({
                        "name": name,
                        "output_basename": output_basename,
                        "group1": group1_label,
                        "group2": group2_label,
                        "resolved_ref": ref_sel,
                        "resolved_sel": sel_sel,
                        "status": "failed",
                        "error": "cn_output_missing_or_empty",
                        "warnings": pair_warnings,
                    })
                    continue

                # Write CN output (.dat summary)
                cn_dat_output = cn_dir / f"cn_{output_basename}.dat"
                try:
                    with open(cn_dat_output, "w") as f:
                        f.write(f"# Coordination Number: {name}\n")
                        f.write(f"# {description}\n")
                        f.write(f"# Group 1: {group1_label} (resolved: {ref_sel})\n")
                        f.write(f"# Group 2: {group2_label} (resolved: {sel_sel})\n")
                        f.write(f"# Cutoff: {cutoff:.4f} nm (mode: {cutoff_mode})\n")
                        f.write(f"# Method: {cn_method}\n")
                        f.write(f"# \n")
                        if cn_value is not None:
                            f.write(f"CN = {cn_value:.3f}\n")
                        else:
                            f.write(f"CN = N/A (computation failed)\n")
                except Exception as e:
                    msg = f"Failed to write CN summary file: {e}"
                    pair_warnings.append(msg)
                    print(f"    [WARN] {msg}")
                    analysis_results.append({
                        "name": name,
                        "output_basename": output_basename,
                        "group1": group1_label,
                        "group2": group2_label,
                        "resolved_ref": ref_sel,
                        "resolved_sel": sel_sel,
                        "status": "failed",
                        "error": "cn_dat_write_failed",
                        "warnings": pair_warnings,
                    })
                    continue

                if not _is_nonempty_file(cn_dat_output):
                    msg = f"CN summary missing/empty: {cn_dat_output}"
                    pair_warnings.append(msg)
                    print(f"    [WARN] {msg}")
                    analysis_results.append({
                        "name": name,
                        "output_basename": output_basename,
                        "group1": group1_label,
                        "group2": group2_label,
                        "resolved_ref": ref_sel,
                        "resolved_sel": sel_sel,
                        "status": "failed",
                        "error": "cn_dat_missing_or_empty",
                        "warnings": pair_warnings,
                    })
                    continue

            cn_str = f"{cn_value:.3f}" if cn_value is not None else "N/A"

            cutoff_reliability = "unreliable" if cutoff_unreliable else "reliable"
            print(
                f"    Cutoff: {cutoff:.3f} nm ({cutoff_mode}, {cutoff_reliability}) | "
                f"CN: {cn_str}"
            )
            
            analysis_results.append({
                "name": name,
                "output_prefix": safe_prefix,
                "output_basename": output_basename,
                "group1": group1_spec.get("value"),
                "group2": group2_spec.get("value"),
                "resolved_ref": ref_sel,
                "resolved_sel": sel_sel,
                "description": description,
                "rdf_settings": {
                    "bin_width": bin_width,
                    "requested_max_r": requested_max_r,
                    "effective_max_r": effective_max_r,
                    "rmax_capped": rmax_capped,
                    "rmax_source": rmax_source,
                    "box_source_used": box_source_used,
                    "safe_rmax_basis": safe_rmax_basis,
                    "min_box_height_tpr_nm": min_box_height_tpr,
                    "min_box_height_traj_nm": min_box_height_traj,
                    # Backward-compat keys
                    "min_box_tpr": min_box_height_tpr,
                    "min_box_traj": min_box_height_traj,
                    "safe_rmax_tpr": safe_rmax_tpr,
                    "safe_rmax_traj": safe_rmax_traj,
                    "safe_rmax": safe_rmax_used,
                    "safe_rmax_used": safe_rmax_used,
                    "no_box_cap_applied": no_box_cap_applied,
                },
                "cn_method": cn_method,
                "cutoff_mode": cutoff_mode,
                "cutoff_method": cutoff_details.get("cutoff_method"),
                "cutoff_nm": cutoff,
                "r_min": cutoff,
                "peak_r": cutoff_details.get("peak_r"),
                "peak_g": cutoff_details.get("peak_g"),
                "cutoff_r": cutoff_details.get("cutoff_r"),
                "g_rmin_raw": cutoff_details.get("g_rmin_raw"),
                "g_rmin_smooth": cutoff_details.get("g_rmin_smooth"),
                "baseline_unreliable": cutoff_details.get("baseline_unreliable", False),
                "cutoff_unreliable": cutoff_unreliable,
                "cutoff_reason": cutoff_details.get("reason"),
                "reliability_flags": cutoff_details.get("reliability_flags", []),
                "quality_flags": cutoff_details.get("reliability_flags", []),
                "fallback_used": bool(cutoff_details.get("fallback_used")),
                "smoothing_params": {
                    "smooth_window_requested": cutoff_details.get("smooth_window_requested"),
                    "smooth_window_used": cutoff_details.get("smooth_window_used"),
                    "smooth_window_nm": pair_cn_settings.get("smooth_window_nm"),
                    "peak_height_threshold_used": cutoff_details.get("peak_height_threshold_used"),
                    "min_sep_nm": cutoff_details.get("min_sep_nm"),
                    "min_depth": cutoff_details.get("min_depth"),
                    "baseline_window_nm": cutoff_details.get("baseline_window_nm"),
                    "recovery_window_nm": cutoff_details.get("recovery_window_nm"),
                    "recovery_height": cutoff_details.get("recovery_height"),
                    "recovery_positive_points": cutoff_details.get("recovery_positive_points"),
                    "local_window_nm": cutoff_details.get("local_window_nm"),
                    "refine_points_used": cutoff_details.get("refine_points_used"),
                    "refine_mode": cutoff_details.get("refine_mode"),
                },
                "cutoff_details": cutoff_details,
                "cn_value": cn_value,
                "output_files": {
                    "rdf_xvg": str(rdf_output) if _is_nonempty_file(rdf_output) else None,
                    "cn_xvg": str(cn_xvg_output) if _is_nonempty_file(cn_xvg_output) else None,
                    "cn_dat": str(cn_dat_output) if (cn_dat_output is not None and _is_nonempty_file(cn_dat_output)) else None,
                },
                "warnings": pair_warnings,
                "status": "completed",
            })
        
        completed_pairs = sum(1 for r in analysis_results if r.get("status") == "completed")
        failed_pairs = sum(1 for r in analysis_results if r.get("status") == "failed")
        skipped_pairs = sum(1 for r in analysis_results if r.get("status") == "skipped")
        pair_error_count = failed_pairs + skipped_pairs

        if strict_analysis and pair_error_count > 0:
            overall_status = "failed"
        elif pair_error_count > 0:
            overall_status = "completed_with_errors"
        else:
            overall_status = "completed"

        # Record in manifest
        if ctx.manifest is not None:
            ctx.manifest.set("analysis", {
                "status": "failed" if overall_status == "failed" else "completed",
                "overall_status": overall_status,
                "config_path": str(config_path),
                "trajectory": str(md_xtc),
                "tpr": str(md_tpr),
                "box_dimensions_nm": list(box_dims_tpr) if box_dims_tpr else None,
                "safe_rmax_basis": "min_altitude",
                "box_min_height_tpr_nm": min_box_height_tpr,
                "box_min_height_traj_nm": min_box_height_traj,
                # Backward-compat keys
                "box_min_tpr_nm": min_box_height_tpr,
                "box_min_traj_nm": min_box_height_traj,
                "box_source": (
                    "traj+tpr"
                    if (min_box_height_tpr is not None and min_box_height_traj is not None)
                    else (
                        "traj"
                        if min_box_height_traj is not None
                        else ("tpr" if min_box_height_tpr is not None else None)
                    )
                ),
                "box_ok": bool(min_box_height_tpr is not None or min_box_height_traj is not None),
                "index_path": str(index_path) if index_path else None,
                "analysis_strict": strict_analysis,
                "gmx_timeout_s": gmx_timeout_s,
                "pair_counts": {
                    "total": len(analysis_results),
                    "completed": completed_pairs,
                    "failed": failed_pairs,
                    "skipped": skipped_pairs,
                },
                "warnings": warnings,
                "pairs": analysis_results,
            })
        
        for w in warnings:
            print(f"  [WARN] {w}")

        if strict_analysis and pair_error_count > 0:
            print(
                f"\n  [ERROR] Analysis strict mode failed: "
                f"{failed_pairs} pair(s) failed, {skipped_pairs} pair(s) skipped"
            )
            return False

        if pair_error_count > 0:
            print(
                f"\n  [WARN] Analysis completed with errors: "
                f"{failed_pairs} pair(s) failed, {skipped_pairs} pair(s) skipped"
            )
        print(f"\n  [OK] Analysis complete ({len(analysis_results)} pairs)")
        return True
    
    def _compute_rdf_with_cn(
        self,
        ctx: "PipelineContext",
        trajectory: Path,
        tpr: Path,
        ref_sel: str,
        sel_sel: str,
        rdf_output: Path,
        cn_output: Path,
        bin_width: float,
        max_r: float,
        index_path: Optional[Path] = None,
        gmx_cmd: Optional[str] = None,
        gmx_timeout_s: int = _DEFAULT_GMX_TIMEOUT_S,
    ) -> Dict[str, Any]:
        """
        Compute RDF and CN using gmx rdf with -cn option.
        
        Returns:
            Dict with 'success', 'cn_method', 'error', 'stderr' keys
        """
        # Build gmx rdf command with -cn for native CN output
        gmx_cmd = gmx_cmd or resolve_gmx_command(ctx)
        try:
            base = _split_gmx_command(gmx_cmd)
        except ValueError as e:
            print(f"    [ERROR] {e}")
            return {"success": False, "error": "invalid_gmx_cmd"}
        cmd = base + [
            "rdf",
            "-f", str(trajectory),
            "-s", str(tpr),
            "-o", str(rdf_output),
            "-cn", str(cn_output),
            "-bin", str(bin_width),
            "-rmax", str(max_r),
            "-ref", ref_sel,
            "-sel", sel_sel,
        ]

        if index_path is not None:
            cmd.extend(["-n", str(index_path)])
        
        # Use shlex.join for proper quoting in log
        cmd_str = shlex.join(cmd)
        print(f"    $ {cmd_str}")
        
        if ctx.dry_run:
            print(f"    [DRY-RUN] Planned command only; RDF/CN files not generated")
            ctx.log_command(cmd_str, self.name, exit_code=0)
            return {"success": True, "cn_method": "dryrun", "dry_run": True}
        
        try:
            result = subprocess.run(
                cmd,
                cwd=str(rdf_output.parent),
                capture_output=True,
                text=True,
                timeout=gmx_timeout_s,
            )
            
            ctx.log_command(
                cmd_str, self.name,
                exit_code=result.returncode,
                stderr=result.stderr[-1000:] if result.stderr else None,
            )
            
            if result.returncode != 0:
                # Provide actionable diagnostic
                stderr_tail = result.stderr[-500:] if result.stderr else ""
                print(f"    [ERROR] gmx rdf failed (exit={result.returncode})")
                print(f"    Command: {cmd_str}")
                return {
                    "success": False,
                    "error": f"exit_code={result.returncode}",
                    "stderr": stderr_tail,
                }
            
            if not rdf_output.exists():
                return {
                    "success": False,
                    "error": "rdf_output_not_created",
                }

            if not _is_nonempty_file(cn_output):
                return {
                    "success": False,
                    "error": "cn_output_missing_or_empty",
                }
            
            return {
                "success": True,
                "cn_method": "gmx_rdf_native",
            }
            
        except subprocess.TimeoutExpired:
            print(f"    [ERROR] gmx rdf timed out after {gmx_timeout_s}s")
            return {"success": False, "error": "timeout"}
        except FileNotFoundError:
            print(f"    [ERROR] GROMACS executable not found: {gmx_cmd}")
            return {"success": False, "error": "gmx_not_found"}
        except Exception as e:
            print(f"    [ERROR] gmx rdf exception: {e}")
            return {"success": False, "error": str(e)}
