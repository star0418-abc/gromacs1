"""
MDP Dynamic Patching Engine.

Copies MDP templates from IN/ to OUT/ and patches values based on CLI options:
- Temperature (ref_t, gen_temp)
- Pressure (ref_p)
- nsteps
- tc-grps mode (System vs Polymer/NonPolymer split)
- comm-grps and comm-mode for COM removal (Issue B)
- gen_vel handling for stage transitions (Issue C)

Safety features:
- Fail-fast on blind tau_t duplication (Issue A)
- Vector length validation (tc-grps vs tau_t/ref_t)
- Override application tracking to eliminate silent no-ops (Issue D)
- Comment formatting preserves original structure (Issue E)
"""

import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING, Union

if TYPE_CHECKING:
    from .context import PipelineContext


# MDP parameters that take temperature values
TEMPERATURE_PARAMS = ["ref_t", "ref-t", "gen_temp", "gen-temp"]

# MDP parameters that take pressure values
PRESSURE_PARAMS = ["ref_p", "ref-p"]

# MDP parameters for nsteps
NSTEPS_PARAMS = ["nsteps"]

# tc-grps related parameters (use canonical hyphenated form)
TC_GRPS_PARAMS = ["tc-grps", "tc_grps", "tau_t", "tau-t"]

# COM removal parameters
COMM_PARAMS = ["comm-mode", "comm_mode", "comm-grps", "comm_grps", "nstcomm"]

# Velocity generation parameters
GENVEL_PARAMS = ["gen_vel", "gen-vel", "gen_temp", "gen-temp", "gen_seed", "gen-seed"]

# Continuation parameters (gen_vel must be consistent with continuation)
CONTINUATION_PARAMS = ["continuation"]

# Unified temperature tolerance for gen_temp/ref_t consistency (K)
# Tighter than previous (was 10K for multi-group, inconsistent)
GEN_TEMP_REF_T_TOLERANCE_K = 0.5

# Heuristic keywords for percolated polymer/gel network COM-risk detection
PERCOLATED_RISK_KEYWORDS = ("polymer", "gel", "network", "matrix")

# Default tau_p values for barostats
# Parrinello-Rahman uses larger value for gel/semi-solid stability
TAU_P_DEFAULT_BERENDSEN = "1.0"
TAU_P_DEFAULT_PR_LIQUID = "2.0"
TAU_P_DEFAULT_PR_GEL = "5.0"  # Safer default for gel/semi-solid systems

# Compressibility defaults for pressure coupling (bar^-1)
# Liquids are often near water-like compressibility; crosslinked gels/solids should
# typically be much less compressible to avoid PR box ringing and LINCS blowups.
COMPRESSIBILITY_DEFAULT_LIQUID = "4.5e-5"
COMPRESSIBILITY_DEFAULT_GEL = "1.0e-6"
COMPRESSIBILITY_WATER_LIKE = 4.5e-5
COMPRESSIBILITY_WATER_LIKE_REL_TOL = 0.25

# Dimensionality map shared by ref_p / compressibility validation and injection.
PCOUPLTYPE_VALUE_COUNTS = {
    "isotropic": 1,
    "semiisotropic": 2,  # Requested spelling in this project
    "semisotropic": 2,   # Common GROMACS spelling
    "anisotropic": 6,
    "surface-tension": 2,
}

# Canonical GROMACS key mappings (prefer hyphenated form per GROMACS style)
CANONICAL_KEYS = {
    "tc_grps": "tc-grps",
    "tau_t": "tau-t",
    "ref_t": "ref-t",
    "ref_p": "ref-p",
    "tau_p": "tau-p",
    "gen_temp": "gen-temp",
    "gen_vel": "gen-vel",
    "gen_seed": "gen-seed",
    "comm_mode": "comm-mode",
    "comm_grps": "comm-grps",
}

# Semantic key groups: each semantic parameter may have multiple alias variants.
# When patching, we prefer the variant already in the template to avoid duplicates.
# Format: semantic_name -> (variant1, variant2, ...) with first as canonical default
SEMANTIC_KEY_GROUPS = {
    "gen_vel": ("gen-vel", "gen_vel"),
    "gen_seed": ("gen-seed", "gen_seed"),
    "gen_temp": ("gen-temp", "gen_temp"),
    "tc_grps": ("tc-grps", "tc_grps"),
    "tau_t": ("tau-t", "tau_t"),
    "ref_t": ("ref-t", "ref_t"),
    "ref_p": ("ref-p", "ref_p"),
    "tau_p": ("tau-p", "tau_p"),
    "comm_grps": ("comm-grps", "comm_grps"),
    "comm_mode": ("comm-mode", "comm_mode"),
    "compressibility": ("compressibility",),
    "continuation": ("continuation",),  # No alias, but include for consistency
    "pcoupl": ("pcoupl",),
    "pcoupltype": ("pcoupltype",),
    "tcoupl": ("tcoupl",),
}

def normalize_key(key: str) -> str:
    """Normalize keys for consistent comparisons."""
    return key.strip().lower()


def _canonicalize_key(key: str) -> str:
    """Convert key to canonical GROMACS form (hyphenated)."""
    key_norm = normalize_key(key)
    return CANONICAL_KEYS.get(key_norm, key_norm)

# Issue B: Reverse mapping from any variant to its semantic key
# Enables proper override normalization
VARIANT_TO_SEMANTIC = {}
for _semantic, _variants in SEMANTIC_KEY_GROUPS.items():
    for _v in _variants:
        VARIANT_TO_SEMANTIC[normalize_key(_v)] = _semantic

# Flattened set of all known semantic key variants (normalized)
SEMANTIC_VARIANTS = {
    normalize_key(v)
    for _variants in SEMANTIC_KEY_GROUPS.values()
    for v in _variants
}


class MDPParseError(Exception):
    """Raised when MDP parsing fails."""
    pass


class MDPPatchError(Exception):
    """Raised when MDP patching fails."""
    pass


class MDPValidationError(Exception):
    """Raised when MDP validation fails (Issue D)."""
    pass


class MDPValidationWarning(UserWarning):
    """Warning for non-fatal MDP validation issues."""
    pass


class OverrideStatus(str, Enum):
    """Outcome states for applying a single override."""
    APPLIED = "applied"
    SKIPPED_WITH_WARNING = "skipped_with_warning"
    FAILED = "failed"
    MISSING_IN_TEMPLATE = "missing_in_template"


@dataclass
class OverrideOutcome:
    """Result of attempting to apply one override."""
    status: OverrideStatus
    message: Optional[str] = None


@dataclass
class DiagnosticSink:
    """
    Collects warnings/errors and applies strict escalation in one place.
    """
    strict: bool = False
    stage: str = ""
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def _with_stage(self, msg: str) -> str:
        if not self.stage:
            return msg
        if f"Stage '{self.stage}'" in msg:
            return msg
        return f"Stage '{self.stage}': {msg}"

    def warn(self, msg: str, escalate_in_strict: bool = True) -> None:
        message = self._with_stage(msg)
        if self.strict and escalate_in_strict:
            self.errors.append(message)
        else:
            self.warnings.append(message)

    def info(self, msg: str) -> None:
        self.warn(f"[INFO] {msg}", escalate_in_strict=False)

    def error(self, msg: str) -> None:
        self.errors.append(self._with_stage(msg))


def find_existing_key(params: Dict[str, str], semantic_name: str) -> Optional[str]:
    """
    Find which variant of a semantic key exists in params.
    
    Args:
        params: Dictionary of MDP parameters
        semantic_name: Semantic key name (e.g., 'ref_t')
        
    Returns:
        The variant key that exists in params, or None if no variant found.
        If multiple variants exist, returns the first one found.
    """
    semantic_norm = normalize_key(semantic_name)
    variants = SEMANTIC_KEY_GROUPS.get(semantic_norm)
    if variants is None:
        # Not a known semantic key group, check directly
        return semantic_norm if semantic_norm in params else None
    
    for variant in variants:
        variant_norm = normalize_key(variant)
        if variant_norm in params:
            return variant_norm
    return None


def set_semantic_param(
    patched_params: Dict[str, str],
    original_params: Dict[str, str],
    semantic_name: str,
    value: str,
) -> str:
    """
    Set a semantic parameter, respecting existing key style (underscore vs hyphen).
    
    This ensures we modify the existing key in-place rather than appending
    a new key with different style, which would create duplicates.
    
    Args:
        patched_params: Parameters being built (will be modified)
        original_params: Original template parameters (read-only, for style detection)
        semantic_name: Semantic key name (e.g., 'ref_t')
        value: Value to set
        
    Returns:
        The chosen key name that was used.
    """
    # Find existing key in original template
    semantic_norm = normalize_key(semantic_name)
    existing_key = find_existing_key(original_params, semantic_norm)
    
    if existing_key is not None:
        # Use the same style as template
        chosen_key = normalize_key(existing_key)
    else:
        # No existing key, use canonical default (first variant)
        variants = SEMANTIC_KEY_GROUPS.get(semantic_norm)
        chosen_key = normalize_key(variants[0] if variants else semantic_norm)
    
    # Remove any other variants from patched_params to avoid duplicates
    variants = SEMANTIC_KEY_GROUPS.get(semantic_norm, (semantic_norm,))
    for variant in variants:
        variant_norm = normalize_key(variant)
        if variant_norm != chosen_key and variant_norm in patched_params:
            del patched_params[variant_norm]
    
    patched_params[chosen_key] = value
    return chosen_key


def get_semantic_value(params: Dict[str, str], semantic_name: str) -> Optional[str]:
    """
    Get the value of a semantic parameter from any of its variant keys.
    
    Args:
        params: Dictionary of MDP parameters
        semantic_name: Semantic key name (e.g., 'ref_t')
        
    Returns:
        The value if found, or None.
    """
    existing_key = find_existing_key(params, semantic_name)
    if existing_key is not None:
        return params.get(existing_key)
    return None


def normalize_group_list(groups: str) -> str:
    """
    Normalize a group list for deterministic comparison and output.
    
    Tokenizes by whitespace, drops empty tokens, joins with single spaces.
    Use this for clean MDP output formatting.
    
    Args:
        groups: Space-separated group names (may have extra whitespace)
        
    Returns:
        Cleaned group string with single spaces between tokens.
    """
    if not groups:
        return ""
    tokens = groups.split()
    return " ".join(t.strip() for t in tokens if t.strip())


def groups_equal(a: str, b: str) -> bool:
    """
    Case-insensitive, whitespace-insensitive group list comparison.
    
    Args:
        a, b: Group list strings to compare
        
    Returns:
        True if groups match after normalization and casefolding.
    """
    return normalize_group_list(a).casefold() == normalize_group_list(b).casefold()


def parse_ndx_groups(ndx_path: Path) -> List[str]:
    """
    Parse group names from a GROMACS index file.
    
    Args:
        ndx_path: Path to .ndx file
        
    Returns:
        List of group names in order of appearance.
        
    Raises:
        MDPParseError: If file cannot be parsed.
    """
    if not ndx_path.exists():
        raise MDPParseError(f"Index file not found: {ndx_path}")
    
    groups = []
    try:
        with open(ndx_path, "r") as f:
            for line in f:
                line = line.strip()
                # Group headers are [ GroupName ]
                if line.startswith("[") and line.endswith("]"):
                    group_name = line[1:-1].strip()
                    if group_name:
                        groups.append(group_name)
    except Exception as e:
        raise MDPParseError(f"Failed to parse index file {ndx_path}: {e}")
    
    return groups


def parse_ndx_group_sizes(ndx_path: Path) -> Dict[str, int]:
    """
    Parse group sizes (atom counts) from a GROMACS index file.

    Args:
        ndx_path: Path to .ndx file

    Returns:
        Mapping: group_name -> atom count

    Raises:
        MDPParseError: If file cannot be parsed.
    """
    if not ndx_path.exists():
        raise MDPParseError(f"Index file not found: {ndx_path}")

    sizes: Dict[str, int] = {}
    current_group: Optional[str] = None

    try:
        with open(ndx_path, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("[") and stripped.endswith("]"):
                    current_group = stripped[1:-1].strip()
                    if current_group and current_group not in sizes:
                        sizes[current_group] = 0
                    continue
                if current_group is None:
                    continue
                if stripped.startswith(";"):
                    continue
                sizes[current_group] = sizes.get(current_group, 0) + len(stripped.split())
    except Exception as e:
        raise MDPParseError(f"Failed to parse index group sizes from {ndx_path}: {e}")

    return sizes


def _split_group_tokens(groups: Optional[str]) -> List[str]:
    """Split a group list string into normalized tokens."""
    if not groups:
        return []
    return [g for g in normalize_group_list(groups).split() if g]


def _decide_auto_tc_grps_mode(ctx: Optional["PipelineContext"]) -> Dict[str, Any]:
    """
    Decide whether auto mode should use split or System tc-grps.

    Split is selected only when:
    - index groups `Polymer` and `NonPolymer` both exist
    - each group meets both absolute atom threshold and atom-fraction threshold
    """
    decision: Dict[str, Any] = {
        "requested_mode": "auto",
        "applied_mode": "system",
        "tc_grps": "System",
        "reason": "auto_mode_defaulted_to_system",
        "warnings": [],
        "group_sizes": {},
        "group_fractions": {},
        "thresholds": {
            "min_atoms": 50,
            "min_fraction": 0.01,
        },
        "ndx_path": None,
    }
    if ctx is None:
        decision["warnings"].append("tc-grps auto mode: missing context; using System.")
        return decision

    ndx_path_raw = getattr(ctx, "active_ndx_path", None) or getattr(ctx, "tc_grps_ndx_path", None)
    if not ndx_path_raw:
        decision["reason"] = "auto_mode_missing_ndx"
        decision["warnings"].append(
            "tc-grps auto mode: no index.ndx available to validate Polymer/NonPolymer groups; using System."
        )
        return decision

    decision["ndx_path"] = str(ndx_path_raw)
    try:
        group_sizes = parse_ndx_group_sizes(Path(str(ndx_path_raw)))
    except MDPParseError as e:
        decision["reason"] = "auto_mode_ndx_parse_failed"
        decision["warnings"].append(
            f"tc-grps auto mode: failed to parse ndx groups from {ndx_path_raw}: {e}. Using System."
        )
        return decision

    min_atoms = getattr(ctx, "tc_grps_min_atoms", 50)
    min_fraction = getattr(ctx, "tc_grps_min_fraction", 0.01)
    try:
        min_atoms = int(min_atoms)
    except (TypeError, ValueError):
        min_atoms = 50
    if min_atoms < 1:
        min_atoms = 1
    try:
        min_fraction = float(min_fraction)
    except (TypeError, ValueError):
        min_fraction = 0.01
    if min_fraction <= 0.0:
        min_fraction = 0.01
    if min_fraction >= 1.0:
        min_fraction = 0.99

    decision["thresholds"] = {
        "min_atoms": min_atoms,
        "min_fraction": min_fraction,
    }

    poly = group_sizes.get("Polymer")
    nonpoly = group_sizes.get("NonPolymer")
    if poly is None or nonpoly is None:
        decision["reason"] = "auto_mode_missing_polymer_groups"
        decision["warnings"].append(
            "tc-grps auto mode: required ndx groups 'Polymer' and 'NonPolymer' were not both found; using System."
        )
        decision["available_groups"] = sorted(group_sizes.keys())
        return decision

    total_atoms = group_sizes.get("System")
    if total_atoms is None or total_atoms <= 0:
        total_atoms = max(poly + nonpoly, sum(v for v in group_sizes.values() if v > 0))
    if total_atoms <= 0:
        total_atoms = poly + nonpoly

    poly_frac = float(poly) / float(total_atoms) if total_atoms > 0 else 0.0
    nonpoly_frac = float(nonpoly) / float(total_atoms) if total_atoms > 0 else 0.0

    decision["group_sizes"] = {
        "Polymer": int(poly),
        "NonPolymer": int(nonpoly),
        "System": int(total_atoms),
    }
    decision["group_fractions"] = {
        "Polymer": poly_frac,
        "NonPolymer": nonpoly_frac,
    }

    too_small = (
        poly < min_atoms
        or nonpoly < min_atoms
        or poly_frac < min_fraction
        or nonpoly_frac < min_fraction
    )
    if too_small:
        decision["reason"] = "auto_mode_groups_too_small"
        decision["warnings"].append(
            "tc-grps auto mode: Polymer/NonPolymer groups are too small for safe split "
            f"(Polymer={poly}, NonPolymer={nonpoly}, min_atoms={min_atoms}, "
            f"min_fraction={min_fraction:.4f}); using System."
        )
        return decision

    decision["applied_mode"] = "split"
    decision["tc_grps"] = "Polymer NonPolymer"
    decision["reason"] = "auto_mode_split_applied"
    return decision


def _has_percolated_risk(groups: Optional[str]) -> bool:
    """
    Heuristic: treat groups containing polymer/gel/network/matrix as risky
    for per-group COM removal in periodic percolated systems.
    """
    tokens = _split_group_tokens(groups)
    for token in tokens:
        lower = token.casefold()
        if any(key in lower for key in PERCOLATED_RISK_KEYWORDS):
            return True
    return False


def _is_pbc_effectively_on(params: Dict[str, str]) -> bool:
    """
    Return whether PBC is effectively enabled.

    If pbc is missing, assume standard xyz PBC is active.
    """
    pbc_raw = params.get("pbc")
    if pbc_raw is None:
        return True
    pbc = pbc_raw.strip().lower()
    return pbc not in {"no", "none", "off", "false", "0"}


def _resolve_pcoupltype_for_dimensionality(
    params: Dict[str, str],
    sink: Optional[DiagnosticSink] = None,
) -> str:
    """
    Resolve pcoupltype and normalize to a supported dimensionality key.

    Defaults to isotropic if missing or unsupported.
    """
    raw = get_semantic_value(params, "pcoupltype")
    if not raw:
        if sink is not None:
            sink.warn(
                "pcoupl is enabled but pcoupltype is not set; assuming isotropic dimensionality.",
                escalate_in_strict=False,
            )
        return "isotropic"

    key = raw.strip().lower()
    if key in PCOUPLTYPE_VALUE_COUNTS:
        return key

    if sink is not None:
        sink.warn(
            f"Unknown pcoupltype='{raw}'; defaulting dimensionality to isotropic.",
            escalate_in_strict=False,
        )
    return "isotropic"


def _looks_water_like_compressibility(values_raw: str) -> bool:
    """
    Heuristic detection for water-like compressibility vectors (~4.5e-5 bar^-1).
    """
    tokens = values_raw.split()
    if not tokens:
        return False
    parsed = []
    for token in tokens:
        try:
            parsed.append(float(token))
        except ValueError:
            return False
    tol = abs(COMPRESSIBILITY_WATER_LIKE) * COMPRESSIBILITY_WATER_LIKE_REL_TOL
    return all(abs(v - COMPRESSIBILITY_WATER_LIKE) <= tol for v in parsed)


def _build_compressibility_vector(default_value: str, count: int) -> str:
    """Build compressibility vector string with required dimensionality."""
    return " ".join([default_value] * count)


def _resolve_checkpoint_known_state(
    checkpoint_available: Optional[bool],
    stage: str,
    sink: Optional[DiagnosticSink] = None,
) -> Optional[bool]:
    """
    Return checkpoint availability as True/False/None without forcing None->False.

    Rule:
    - True/False: trusted
    - None: unknown, keep unknown and emit an info warning for traceability.
    """
    if checkpoint_available is None:
        if sink is not None:
            sink.info(
                f"checkpoint availability is unknown; "
                f"checkpoint-specific continuation safeguards cannot be fully verified."
            )
        return None
    return bool(checkpoint_available)


def _validate_gen_vel_continuation_consistency(
    params: Dict[str, str], stage: str, strict: bool = False,
    ctx: Optional["PipelineContext"] = None,
    checkpoint_available: Optional[bool] = None,
) -> Tuple[List[str], List[str]]:
    """
    Validate gen_vel/continuation consistency to prevent LINCS instability.
    
    In GROMACS:
    - continuation=no resets constraints/velocities at step 0
    - continuation=yes assumes velocities from checkpoint
    
    Inconsistent combinations can cause LINCS warnings or blowups,
    especially for stressed GPE systems.
    
    P0.2 / Safety A: Enhanced support for expert restart from .gro with velocities.
    P0.B / Safety B: Clear messaging about input requirements for continuation.
    
    Args:
        params: MDP parameters
        stage: Stage name for context
        strict: If True, escalate warnings to errors
        ctx: Pipeline context for allow_gro_velocity_restart flag
        checkpoint_available: True/False if known; None if unknown.
            Unknown is preserved (never silently treated as False).
        
    Returns:
        Tuple of (errors, warnings) lists.
    """
    errors = []
    warn_list = []
    
    gen_vel = get_semantic_value(params, "gen_vel")
    continuation = get_semantic_value(params, "continuation")
    
    # Normalize to lowercase for comparison
    gen_vel_lower = gen_vel.lower() if gen_vel else None
    continuation_lower = continuation.lower() if continuation else None
    
    # Extract flags from context
    allow_gro_velocity_restart = getattr(ctx, "allow_gro_velocity_restart", False) if ctx else False
    input_gro_has_velocities = getattr(ctx, "input_gro_has_velocities", None) if ctx else None
    
    # =========================================================================
    # P0.B: Continuation input requirement messaging
    # =========================================================================
    if continuation_lower == "yes":
        info_msg = (
            f"[INFO] Stage '{stage}' uses continuation=yes.\n"
            f"  Input source MUST provide velocities/state from either:\n"
            f"    - Checkpoint file (.cpt) from previous stage, OR\n"
            f"    - Input .gro with embedded velocities\n"
            f"  The MDP patcher cannot verify file contents. Ensure your workflow\n"
            f"  provides the appropriate input before running grompp/mdrun."
        )
        # Log to warnings list (INFO level, not error)
        warn_list.append(info_msg)
    
    # =========================================================================
    # Case 1: gen_vel=no AND continuation=no (P0.2 - expert restart guard)
    # =========================================================================
    if gen_vel_lower == "no" and continuation_lower == "no":
        # Improved diagnostic message with clear risk explanation
        base_msg = (
            f"Stage '{stage}' has gen_vel=no and continuation=no.\n\n"
            f"RISK: This mode requires the input .gro to contain velocities.\n"
            f"      The MDP patcher CANNOT verify whether the .gro has velocities.\n"
            f"      If velocities are missing, the run may start with:\n"
            f"        - Zero velocities (near-0K behavior)\n"
            f"        - Unintended dynamics / equilibration artifacts\n\n"
            f"SOLUTIONS:\n"
            f"  1. Use checkpoint continuation: continuation=yes with .cpt file\n"
            f"  2. Generate velocities: gen_vel=yes with gen_temp\n"
            f"  3. Expert restart: --allow-gro-velocity-restart (YOU must ensure .gro has velocities)"
        )
        
        # Checkpoint-context guard: if checkpoint available, always error
        # (prevents accidentally dropping a checkpoint when one exists)
        if checkpoint_available is True:
            checkpoint_error = (
                f"{base_msg}\n\n"
                f"ADDITIONAL ERROR: A checkpoint file (.cpt) is available for this stage.\n"
                f"  The (gen_vel=no, continuation=no) combo is disallowed when a checkpoint exists,\n"
                f"  even with --allow-gro-velocity-restart, to prevent accidentally dropping\n"
                f"  the checkpoint. Use continuation=yes to continue from the checkpoint."
            )
            errors.append(checkpoint_error)
        elif allow_gro_velocity_restart:
            # Expert mode: input .gro velocity availability must be checked.
            if input_gro_has_velocities is False:
                errors.append(
                    f"Stage '{stage}' requested expert restart "
                    f"(gen_vel=no, continuation=no, --allow-gro-velocity-restart), "
                    f"but the selected input .gro is known to lack velocities.\n"
                    f"This can silently start from zero velocities; the thermostat may inject "
                    f"large energy, causing hotspots/constraint failures.\n"
                    f"Solutions:\n"
                    f"  1) Provide a .gro that includes velocities, or\n"
                    f"  2) Set gen_vel=yes with gen_temp, or\n"
                    f"  3) Use continuation=yes with a checkpoint (.cpt)."
                )
            elif input_gro_has_velocities is None:
                warn_list.append(
                    f"[WARNING] Stage '{stage}' using expert restart mode "
                    f"(gen_vel=no, continuation=no) with unknown .gro velocity status.\n"
                    f"If the input .gro lacks velocities, thermostat startup can inject large "
                    f"energy and trigger hotspots/constraint failures.\n"
                    f"Verify velocities exist, or switch to gen_vel=yes with gen_temp, "
                    f"or continuation=yes with a checkpoint."
                )
            else:
                warn_list.append(
                    f"[INFO] Stage '{stage}' expert restart acknowledged: "
                    f"gen_vel=no, continuation=no with validated .gro velocities."
                )
        elif checkpoint_available is None:
            msg = (
                f"{base_msg}\n\n"
                f"Checkpoint availability is unknown. "
                f"The validator cannot enforce checkpoint-preference safeguards."
            )
            if strict:
                errors.append(msg)
            else:
                warn_list.append(f"[WARNING] {msg}")
        else:
            # Default behavior: WARNING (non-strict) or ERROR (strict)
            if strict:
                errors.append(base_msg)
            else:
                warn_list.append(f"[WARNING] {base_msg}")
    
    # =========================================================================
    # Case 2: gen_vel=yes AND continuation=yes (always contradictory)
    # =========================================================================
    if gen_vel_lower == "yes" and continuation_lower == "yes":
        # This is ALWAYS an error - contradictory settings
        msg = (
            f"Stage '{stage}' has gen_vel=yes and continuation=yes.\n"
            f"This is contradictory: gen_vel=yes generates new velocities,\n"
            f"but continuation=yes implies using checkpoint velocities.\n"
            f"Set gen_vel=no for true continuation, or continuation=no for fresh start.\n\n"
            f"This combination is ALWAYS an error regardless of strict mode."
        )
        # Always error, never just a warning
        errors.append(msg)
    
    return errors, warn_list


def _dedupe_semantic_keys(
    params: Dict[str, str],
    patched_values: Optional[Dict[str, str]] = None,
    sink: Optional[DiagnosticSink] = None,
) -> Dict[str, str]:
    """
    Remove duplicate alias keys for semantic parameters.
    
    SAFETY NET: The primary deduplication happens at patch time
    (copy_and_patch_mdp Task A fix). This function catches any
    remaining edge cases that slip through.
    
    For each semantic key group, ensure only one variant exists.
    Priority for which variant to keep:
      1. Variant that was patched in this run (present in patched_values)
      2. First variant found (in order of SEMANTIC_KEY_GROUPS tuple)
    
    Args:
        params: Dictionary of MDP parameters
        patched_values: Optional dict of values patched in this run.
            If provided, prefer keeping variants present here.
        
    Returns:
        New dictionary with duplicates removed.
    """
    result = dict(params)
    patched_values = {normalize_key(k): v for k, v in (patched_values or {}).items()}
    
    for semantic_name, variants in SEMANTIC_KEY_GROUPS.items():
        # Find which variants exist
        existing = [normalize_key(v) for v in variants if normalize_key(v) in result]
        
        if len(existing) > 1:
            # Multiple variants exist - determine which to keep
            # Priority: prefer variant that was patched in this run
            keep_key = None
            for v in existing:
                if v in patched_values:
                    keep_key = v
                    break
            
            # Fallback: keep first if no patched variant found
            if keep_key is None:
                keep_key = existing[0]
            
            for variant in existing:
                if variant == keep_key:
                    continue
                # Merge values if different (prefer kept key)
                if sink is not None:
                    sink.warn(
                        f"MDP contains duplicate semantic key variants: {existing}. "
                        f"Keeping '{keep_key}' (patched={keep_key in patched_values}) "
                        f"and removing '{variant}'."
                    )
                del result[variant]
    
    return result


def parse_mdp(mdp_path: Path) -> Dict[str, str]:
    """
    Parse an MDP file into a dictionary.
    
    Args:
        mdp_path: Path to MDP file
        
    Returns:
        Dictionary mapping parameter names to values (as strings)
    """
    params = {}
    
    if not mdp_path.exists():
        raise MDPParseError(f"MDP file not found: {mdp_path}")
    
    with open(mdp_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            # Strip comments
            if ";" in line:
                line = line.split(";", 1)[0]
            
            line = line.strip()
            if not line:
                continue
            
            # Parse key = value
            if "=" in line:
                key, value = line.split("=", 1)
                key = normalize_key(key)  # MDP params are case-insensitive
                value = value.strip()
                params[key] = value
    
    return params


def write_mdp(
    mdp_path: Path,
    params: Dict[str, str],
    original_content: Optional[str] = None,
    preserve_comments: bool = True,
    patched_values: Optional[Dict[str, str]] = None,
    strict: bool = False,
    sink: Optional[DiagnosticSink] = None,
) -> None:
    """
    Write an MDP file, preserving comments if possible.
    
    If original_content is provided, patches in-place.
    Otherwise, writes a fresh file.
    
    Args:
        mdp_path: Output path
        params: Parameters to write
        original_content: Original file content (for preserving comments)
        preserve_comments: Whether to preserve original comments/structure
        patched_values: Values patched in this run (for dedupe priority)
        strict: If True, fail fast on duplicate semantic key variants in template
    
    Issue E fix: Proper comment formatting to avoid "; ;comment"
    Task A fix: Dedupe semantic keys before writing to avoid duplicates
    Task D fix: Prefer last-wins on duplicate semantic keys (or fail fast in strict mode)
    """
    # Task A: Dedupe semantic keys before writing, preferring patched variants
    params = _dedupe_semantic_keys(params, patched_values=patched_values, sink=sink)
    
    if original_content and preserve_comments:
        # Issue 6: Line-level semantic dedupe - ensure no alias duplicates remain
        # First pass: build map of semantic_name -> list of line indices (in order)
        lines = original_content.split("\n")
        semantic_occurrences = {}  # semantic_name -> list[(idx, key)]
        
        for idx, line in enumerate(lines):
            stripped = line.split(";", 1)[0].strip()
            if "=" in stripped:
                key = normalize_key(stripped.split("=", 1)[0])
                semantic_name = VARIANT_TO_SEMANTIC.get(key)
                if semantic_name:
                    semantic_occurrences.setdefault(semantic_name, []).append((idx, key))
        
        duplicates = {
            s: occs for s, occs in semantic_occurrences.items() if len(occs) > 1
        }
        
        if duplicates:
            if strict:
                details = []
                for semantic_name, occs in duplicates.items():
                    line_info = ", ".join(f"{idx+1}:{key}" for idx, key in occs)
                    details.append(f"  - {semantic_name}: {line_info}")
                raise MDPValidationError(
                    "MDP template contains duplicate semantic key variants:\n"
                    + "\n".join(details)
                )
            else:
                for semantic_name, occs in duplicates.items():
                    line_info = ", ".join(f"{idx+1}:{key}" for idx, key in occs)
                    keep_line = occs[-1][0] + 1
                    if sink is not None:
                        sink.warn(
                            f"MDP contains duplicate semantic key variants for '{semantic_name}': "
                            f"{line_info}. Keeping last occurrence at line {keep_line} "
                            f"and removing earlier duplicates (last-wins)."
                        )
        
        lines_to_skip = set()  # line indices to skip (duplicate aliases)
        for occs in duplicates.values():
            for idx, _key in occs[:-1]:
                lines_to_skip.add(idx)
        
        # Second pass: rebuild lines, patching and skipping duplicates
        new_lines = []
        params_written = set()  # Track params by semantic name
        params_written_exact = set()  # Track exact keys written
        
        for idx, line in enumerate(lines):
            if idx in lines_to_skip:
                # Skip duplicate alias lines (or optionally comment them out)
                # new_lines.append(f"; [REMOVED DUPLICATE] {line}")
                continue
                
            # Check if this line contains a param assignment
            stripped = line.split(";", 1)[0].strip()
            if "=" in stripped:
                key = normalize_key(stripped.split("=", 1)[0])
                
                # Check semantic equivalence
                semantic_name = VARIANT_TO_SEMANTIC.get(key)
                
                # If we have params to patch and this key matches (exact or semantic)
                if key in params or (semantic_name and any(
                    normalize_key(v) in params for v in SEMANTIC_KEY_GROUPS.get(semantic_name, ())
                )):
                    # Find the param value (check exact key first, then semantic variants)
                    param_value = params.get(key)
                    if param_value is None and semantic_name:
                        for v in SEMANTIC_KEY_GROUPS.get(semantic_name, ()):
                            v_norm = normalize_key(v)
                            if v_norm in params:
                                param_value = params[v_norm]
                                break
                    
                    if param_value is not None:
                        # Replace value, preserve comment (Issue E: fix formatting)
                        comment_part = ""
                        if ";" in line:
                            raw_comment = line.split(";", 1)[1]
                            # Strip leading semicolons and spaces from comment to avoid "; ;"
                            clean_comment = raw_comment.lstrip("; ")
                            if clean_comment:
                                comment_part = f"  ; {clean_comment}"
                        new_lines.append(f"{key} = {param_value}{comment_part}")
                        params_written_exact.add(key)
                        if semantic_name:
                            params_written.add(semantic_name)
                        continue
            
            new_lines.append(line)
        
        # Append any new params not in original (checking semantic equivalence)
        for key, value in params.items():
            # Skip if already written (exact match)
            if key in params_written_exact:
                continue
            # Skip if semantic equivalent already written
            for sname, variants in SEMANTIC_KEY_GROUPS.items():
                if key in [normalize_key(v) for v in variants] and sname in params_written:
                    break
            else:
                new_lines.append(f"{key} = {value}")
        
        content = "\n".join(new_lines)
    else:
        # Write fresh file
        lines = ["; Generated MDP file", ""]
        for key, value in sorted(params.items()):
            lines.append(f"{key} = {value}")
        content = "\n".join(lines)
    
    mdp_path.parent.mkdir(parents=True, exist_ok=True)
    mdp_path.write_text(content)


def _validate_vector_lengths(
    params: Dict[str, str], stage: str, strict: bool = False
) -> Tuple[List[str], List[str]]:
    """
    Validate that vector parameters have consistent lengths.
    
    Issue A/D: Ensure len(tc-grps) == len(tau_t) == len(ref_t)
    Task E: Returns (errors, warnings) tuple. comm-grps mismatch is a warning by default.
    
    Args:
        params: MDP parameters
        stage: Stage name for context
        strict: If True, escalate warnings to errors
        
    Returns:
        Tuple of (errors, warnings) where each is a list of messages
    """
    errors = []
    warn_list = []  # renamed to avoid shadowing warnings module
    
    # Find tc-grps using semantic key helper
    tc_grps_val = get_semantic_value(params, "tc_grps")
    if tc_grps_val is None:
        return errors, warn_list  # No tc-grps, nothing to validate
    
    n_groups = len(tc_grps_val.split())
    
    # Check tau_t
    tau_t_val = get_semantic_value(params, "tau_t")
    if tau_t_val:
        n_tau = len(tau_t_val.split())
        if n_tau != n_groups:
            errors.append(
                f"Vector length mismatch: tc-grps has {n_groups} groups, "
                f"but tau_t has {n_tau} values"
            )
    
    # Check ref_t
    ref_t_val = get_semantic_value(params, "ref_t")
    if ref_t_val:
        n_ref = len(ref_t_val.split())
        if n_ref != n_groups:
            errors.append(
                f"Vector length mismatch: tc-grps has {n_groups} groups, "
                f"but ref_t has {n_ref} values"
            )
    
    # Check comm-grps only when COM removal is enabled.
    comm_mode_val = get_semantic_value(params, "comm_mode")
    comm_mode_none = bool(comm_mode_val and comm_mode_val.strip().lower() == "none")
    comm_grps_val = get_semantic_value(params, "comm_grps")
    if comm_grps_val and comm_grps_val.lower() != "system" and not comm_mode_none:
        n_comm = len(comm_grps_val.split())
        if n_comm != n_groups and n_groups > 1:
            msg = (
                f"comm-grps has {n_comm} groups but tc-grps has {n_groups} groups. "
                f"This may cause COM drift / 'flying ice cube' artifacts. "
                f"RECOMMENDED: Use comm-grps-policy auto/system, or match comm-grps to tc-grps. "
                f"Expert option: comm-grps-policy none disables COM removal and requires drift monitoring."
            )
            if strict:
                errors.append(msg)
            else:
                warn_list.append(msg)
    
    return errors, warn_list


def _validate_ref_p_dimensionality(
    params: Dict[str, str], stage: str, strict: bool = False
) -> Tuple[List[str], List[str]]:
    """
    Validate ref_p dimensionality matches pcoupltype (Issue C).
    
    Rules:
      - pcoupl=no: skip check (pressure coupling disabled)
      - pcoupltype=isotropic: ref_p must be 1 value
      - pcoupltype=semisotropic: ref_p must be 2 values
      - pcoupltype=anisotropic: ref_p must be 6 values (upper-triangular tensor)
    
    Args:
        params: MDP parameters
        stage: Stage name for context
        strict: If True, escalate warnings to errors
        
    Returns:
        Tuple of (errors, warnings) lists.
    """
    errors = []
    warn_list = []
    
    # Get pcoupl value
    pcoupl = get_semantic_value(params, "pcoupl")
    if not pcoupl or pcoupl.lower() == "no":
        # Pressure coupling disabled, skip check
        return errors, warn_list
    
    local_sink = DiagnosticSink(strict=strict, stage=stage)
    pcoupltype_lower = _resolve_pcoupltype_for_dimensionality(params, sink=local_sink)
    
    # Get ref_p
    ref_p = get_semantic_value(params, "ref_p")
    if not ref_p:
        # ref_p not specified - GROMACS will use defaults, skip validation
        return errors, warn_list
    
    ref_p_values = ref_p.split()
    n_ref_p = len(ref_p_values)
    
    expected = PCOUPLTYPE_VALUE_COUNTS.get(pcoupltype_lower)
    
    if expected is None:
        # Unknown pcoupltype - emit info and skip
        warn_list.append(
            f"Stage '{stage}': Unknown pcoupltype='{pcoupltype_lower}', skipping ref_p validation."
        )
        return errors, warn_list
    
    if n_ref_p != expected:
        msg = (
            f"Stage '{stage}': ref_p dimensionality mismatch.\n"
            f"  pcoupltype={pcoupltype_lower} requires {expected} value(s), "
            f"but ref_p has {n_ref_p}: '{ref_p}'.\n"
            f"  Solutions:\n"
            f"    - For isotropic: ref_p = 1.0\n"
            f"    - For semiisotropic: ref_p = 1.0 1.0 (xy, z)\n"
            f"    - For anisotropic: ref_p = xx yy zz xy xz yz (6 values)"
        )
        if strict:
            errors.append(msg)
        else:
            warn_list.append(f"[WARNING] {msg}")
    warn_list.extend(local_sink.warnings)
    if strict:
        errors.extend(local_sink.errors)
    else:
        warn_list.extend(local_sink.errors)

    return errors, warn_list


def _validate_compressibility_dimensionality(
    params: Dict[str, str], stage: str, strict: bool = False
) -> Tuple[List[str], List[str]]:
    """
    Validate compressibility dimensionality matches pcoupltype.
    """
    errors: List[str] = []
    warn_list: List[str] = []

    pcoupl = get_semantic_value(params, "pcoupl")
    if not pcoupl or pcoupl.lower() == "no":
        return errors, warn_list

    compressibility = get_semantic_value(params, "compressibility")
    if not compressibility:
        return errors, warn_list

    local_sink = DiagnosticSink(strict=strict, stage=stage)
    pcoupltype = _resolve_pcoupltype_for_dimensionality(params, sink=local_sink)
    expected = PCOUPLTYPE_VALUE_COUNTS.get(pcoupltype)
    if expected is None:
        warn_list.extend(local_sink.warnings)
        return errors, warn_list

    n_values = len(compressibility.split())
    if n_values != expected:
        msg = (
            f"Stage '{stage}': compressibility dimensionality mismatch.\n"
            f"  pcoupltype={pcoupltype} requires {expected} value(s), "
            f"but compressibility has {n_values}: '{compressibility}'.\n"
            f"  Use 1 (isotropic), 2 (semiisotropic), or 6 (anisotropic) values."
        )
        if strict:
            errors.append(msg)
        else:
            warn_list.append(f"[WARNING] {msg}")

    warn_list.extend(local_sink.warnings)
    if strict:
        errors.extend(local_sink.errors)
    else:
        warn_list.extend(local_sink.errors)
    return errors, warn_list


def _validate_comm_mode_pbc_safety(
    params: Dict[str, str],
    stage: str,
    strict: bool = False,
    ctx: Optional["PipelineContext"] = None,
) -> Tuple[List[str], List[str]]:
    """
    Validate that comm-mode=Angular is not used with periodic boundaries by default.
    """
    errors: List[str] = []
    warn_list: List[str] = []

    comm_mode = get_semantic_value(params, "comm_mode")
    if not comm_mode or comm_mode.strip().lower() != "angular":
        return errors, warn_list
    if not _is_pbc_effectively_on(params):
        return errors, warn_list

    allow = getattr(ctx, "allow_comm_mode_angular_with_pbc", False) if ctx else False
    pbc_raw = params.get("pbc", "xyz")
    msg = (
        f"Stage '{stage}': comm-mode=Angular with pbc='{pbc_raw}' is unsafe/ill-defined under PBC "
        f"and can yield incorrect momentum removal."
    )
    if allow:
        warn_list.append(
            f"[WARNING] {msg} Expert override keep is enabled "
            f"(allow_comm_mode_angular_with_pbc=True); proceed at your own risk."
        )
    elif strict:
        errors.append(
            f"{msg} Strict validation forbids this combination. "
            f"Use comm-mode=Linear or disable PBC."
        )
    else:
        warn_list.append(f"[WARNING] {msg} The patcher will force Linear where policy logic applies.")
    return errors, warn_list


def _validate_gen_vel_consistency(
    params: Dict[str, str],
    stage: str,
    strict: bool = False,
    ctx: Optional["PipelineContext"] = None,
) -> Tuple[List[str], List[str]]:
    """
    Validate gen_vel/gen_temp consistency.
    
    Issue C: If gen_vel=yes, ensure gen_temp matches ref_t.
    P0.3: If gen_vel=yes, ensure gen_temp is present and numeric.
    P1.7: Don't silently swallow numeric parse failures.
    Task E: Returns (errors, warnings) tuple. gen_seed missing is warning by default.
    Task F: When ref_t has multiple values, require gen_temp close to EACH value
            or all ref_t values consistent. Do NOT use average-based pass condition.
    Task B: Allow ctx.allow_gen_vel_with_ref_t_spread to downgrade strict error to warning.
    
    Args:
        params: MDP parameters
        stage: Stage name for context
        strict: If True, escalate warnings to errors
        
    Returns:
        Tuple of (errors, warnings) where each is a list of messages
    """
    errors = []
    warn_list = []
    
    gen_vel = get_semantic_value(params, "gen_vel")
    allow_ref_t_spread = getattr(ctx, "allow_gen_vel_with_ref_t_spread", False) if ctx else False
    if gen_vel and gen_vel.lower() == "yes":
        gen_temp = get_semantic_value(params, "gen_temp")
        ref_t = get_semantic_value(params, "ref_t")
        
        # =====================================================================
        # P0.3: Check gen_temp is present when gen_vel=yes
        # =====================================================================
        if not gen_temp:
            msg = (
                f"Stage '{stage}': gen_vel=yes but gen_temp is not set.\n"
                f"GROMACS will use gen_temp=0 by default, generating zero velocities.\n"
                f"This will cause near-0K behavior or undefined dynamics.\n"
                f"Set gen_temp to match ref_t (e.g., gen_temp = 300)."
            )
            if strict:
                errors.append(msg)
            else:
                warn_list.append(f"[WARNING] {msg}")
        
        # =====================================================================
        # P0.3 + P1.7: Check gen_temp is numeric (catch "300K" style)
        # =====================================================================
        if gen_temp:
            try:
                gen_t = float(gen_temp)
            except ValueError:
                # P1.7: Don't swallow - report the parse failure
                msg = (
                    f"Stage '{stage}': gen_temp='{gen_temp}' is not a valid number.\n"
                    f"Expected numeric value in Kelvin (e.g., 300, not '300K').\n"
                    f"GROMACS will fail or use incorrect temperature."
                )
                if strict:
                    errors.append(msg)
                else:
                    warn_list.append(f"[WARNING] {msg}")
                gen_t = None
        else:
            gen_t = None
        
        # =====================================================================
        # Validate gen_temp vs ref_t consistency (if both are valid)
        # =====================================================================
        if gen_t is not None and ref_t:
            ref_t_values = ref_t.split()
            # Issue 5: Use unified tolerance for both single and multi-group cases
            tolerance = GEN_TEMP_REF_T_TOLERANCE_K
            
            # P1.7: Parse ref_t with error handling
            try:
                ref_ts = [float(x) for x in ref_t_values]
            except ValueError:
                msg = (
                    f"Stage '{stage}': ref_t='{ref_t}' contains non-numeric values.\n"
                    f"Expected numeric values in Kelvin (e.g., '300 300', not '300K 300K')."
                )
                if strict:
                    errors.append(msg)
                else:
                    warn_list.append(f"[WARNING] {msg}")
                ref_ts = None
            
            if ref_ts is not None:
                if len(ref_ts) == 1:
                    # Single ref_t value - straightforward check
                    ref_t_single = ref_ts[0]
                    if abs(gen_t - ref_t_single) > tolerance:
                        errors.append(
                            f"gen_temp ({gen_temp}) differs from ref_t ({ref_t}) by "
                            f"{abs(gen_t - ref_t_single):.2f}K (tolerance={tolerance}K); "
                            f"velocities will be generated at wrong temperature"
                        )
                else:
                    # Task F: Multiple ref_t values - check for thermal gradient trap
                    # Check if ref_t values are consistent with each other
                    ref_t_spread = max(ref_ts) - min(ref_ts)
                    all_consistent = ref_t_spread <= tolerance
                    
                    # Check if gen_temp is close to ALL ref_t values
                    all_close = all(abs(gen_t - r) <= tolerance for r in ref_ts)
                    
                    if not all_consistent:
                        # Thermal gradient in ref_t - dangerous with gen_vel=yes
                        msg = (
                            f"Stage '{stage}': gen_vel=yes but ref_t has divergent values {ref_ts} "
                            f"(spread={ref_t_spread:.1f}K > {tolerance}K tolerance).\n"
                            f"gen_temp is a scalar; velocities are generated at ONE temperature "
                            f"({gen_t}K), then the thermostat drives each group to its own ref_t. "
                            f"This can cause thermal shock and instabilities.\n"
                            f"Fixes:\n"
                            f"  1) Use tc-grps-mode system (single thermostat group), or\n"
                            f"  2) Use a single ref_t during equilibration and split later, or\n"
                            f"  3) Set gen_vel=no and continue from a checkpoint, or\n"
                            f"  4) Advanced: explicitly set gen_temp and set "
                            f"ctx.allow_gen_vel_with_ref_t_spread=True to acknowledge risk."
                        )
                        if strict and not allow_ref_t_spread:
                            errors.append(msg)
                        else:
                            warn_list.append(f"[WARNING] {msg}")
                    elif not all_close:
                        # ref_t is consistent but gen_temp doesn't match
                        msg = (
                            f"gen_temp ({gen_temp}K) not consistent with ref_t values "
                            f"{ref_ts} (tolerance={tolerance}K)."
                        )
                        if strict:
                            errors.append(msg)
                        else:
                            warn_list.append(msg)
        
        # Task E: gen_seed missing is a warning by default (not error)
        gen_seed = get_semantic_value(params, "gen_seed")
        if not gen_seed or gen_seed == "-1":
            msg = (
                f"gen_vel=yes but gen_seed not set or is -1; "
                f"velocities will be non-reproducible. "
                f"Use --gen-seed <N> for reproducibility."
            )
            if strict:
                errors.append(msg)
            else:
                warn_list.append(msg)
    
    return errors, warn_list


def _normalize_tc_grps_vectors(
    patched_params: Dict[str, str],
    original_values: Dict[str, str],
    patched_values: Dict[str, str],
    ctx: Optional["PipelineContext"],
    original_params: Optional[Dict[str, str]] = None,
    sink: Optional[DiagnosticSink] = None,
) -> None:
    """
    Final normalization pass to align ref_t/tau_t with tc-grps length.
    
    - If tc-grps has N groups:
      - ref_t length 1 → replicate to N
      - ref_t length not in {1, N} → error
      - tau_t length 1 → replicate ONLY if allow_tau_t_autofill
      - tau_t length not in {1, N} → error
    """
    tc_grps_val = get_semantic_value(patched_params, "tc_grps")
    if not tc_grps_val:
        return
    
    n_groups = len(tc_grps_val.split())
    if n_groups <= 1:
        return
    
    style_ref = original_params if original_params else patched_params
    allow_autofill = getattr(ctx, "allow_tau_t_autofill", False) if ctx else False
    
    # Normalize ref_t
    ref_t_val = get_semantic_value(patched_params, "ref_t")
    if ref_t_val:
        ref_tokens = ref_t_val.split()
        if len(ref_tokens) == 1:
            new_ref_t = " ".join([ref_tokens[0]] * n_groups)
            existing_key = find_existing_key(patched_params, "ref_t")
            if existing_key:
                original_values.setdefault(existing_key, ref_t_val)
            chosen_key = set_semantic_param(patched_params, style_ref, "ref_t", new_ref_t)
            patched_values[normalize_key(chosen_key)] = new_ref_t
        elif len(ref_tokens) != n_groups:
            raise MDPPatchError(
                f"ref_t has {len(ref_tokens)} values but tc-grps has {n_groups} groups. "
                f"Provide 1 value (will replicate) or {n_groups} values."
            )
    
    # Normalize tau_t
    tau_t_val = get_semantic_value(patched_params, "tau_t")
    if tau_t_val:
        tau_tokens = tau_t_val.split()
        if len(tau_tokens) == 1:
            if allow_autofill:
                single_val = tau_tokens[0]
                new_tau_t = " ".join([single_val] * n_groups)
                if sink is not None:
                    sink.warn(
                        f"Final normalization duplicated tau_t={single_val} for {n_groups} groups "
                        f"(allow_tau_t_autofill=True). Consider explicit per-group tau_t values."
                    )
                existing_key = find_existing_key(patched_params, "tau_t")
                if existing_key:
                    original_values.setdefault(existing_key, tau_t_val)
                chosen_key = set_semantic_param(patched_params, style_ref, "tau_t", new_tau_t)
                patched_values[normalize_key(chosen_key)] = new_tau_t
            else:
                raise MDPPatchError(
                    f"tau_t has 1 value but tc-grps has {n_groups} groups. "
                    f"Provide {n_groups} values or enable allow_tau_t_autofill."
                )
        elif len(tau_tokens) != n_groups:
            raise MDPPatchError(
                f"tau_t has {len(tau_tokens)} values but tc-grps has {n_groups} groups. "
                f"Provide {n_groups} values."
            )


def copy_and_patch_mdp(
    template_path: Path,
    output_path: Path,
    overrides: Dict[str, Any],
    validate: bool = True,
    ctx: Optional["PipelineContext"] = None,
    stage: str = "",
    checkpoint_available: Optional[bool] = None,
    return_diagnostics: bool = False,
) -> Union[
    Tuple[Dict[str, str], Dict[str, str]],
    Tuple[Dict[str, str], Dict[str, str], List[str], List[str]],
]:
    """
    Copy MDP template and apply patches.
    
    Args:
        template_path: Path to template MDP in IN/
        output_path: Path to output MDP in OUT/
        overrides: Dict of parameter overrides {param_name: new_value}
        validate: If True, fail if required patches cannot be applied
        ctx: Pipeline context for accessing safety options
        stage: Stage name for context-aware validation
        checkpoint_available: If True/False known; None means unknown.
            If None, availability is treated as unknown (never silently assumed false).
        return_diagnostics: If True, also return (errors, warnings) diagnostic lists.
        
    Returns:
        Tuple of (original_values, patched_values), optionally plus diagnostics.
        
    Raises:
        MDPPatchError: If patching fails
        MDPValidationError: If validate=True and overrides didn't apply (Issue D)
    """
    if not template_path.exists():
        raise MDPPatchError(f"MDP template not found: {template_path}")
    
    # Read original content
    original_content = template_path.read_text()
    original_params = parse_mdp(template_path)
    
    # Record original values for patched params
    original_values = {}
    patched_values = {}

    # Build patched params
    patched_params = dict(original_params)

    strict = getattr(ctx, "strict_mdp_validation", False) if ctx else False
    sink = DiagnosticSink(strict=strict, stage=stage)

    # Track how each override resolved to avoid false "unapplied" errors.
    applied_overrides: Dict[str, OverrideOutcome] = {}

    # Phase 1: normalize overrides into a stable list (do NOT rely on dict order)
    normalized_overrides = []
    for param, value in overrides.items():
        param_norm = normalize_key(param)
        normalized_overrides.append((param, value, param_norm))
    
    # Phase 1b: categorize overrides by dependency
    mode_keys = {
        "tc_grps_mode", "tc-grps-mode", "tc_grps", "tc-grps",
        "comm_grps_policy", "comm-grps-policy",
    }
    scalar_temp_keys = {"temperature", "pressure", "gen_temp", "gen-temp"}
    vector_keys = {"ref_t", "ref-t", "tau_t", "tau-t", "ref_p", "ref-p", "tau_p", "tau-p"}
    
    mode_items = []
    scalar_items = []
    vector_items = []
    other_items = []
    
    for param, value, param_norm in normalized_overrides:
        if param_norm in mode_keys:
            mode_items.append((param, value, param_norm))
        elif param_norm in scalar_temp_keys:
            scalar_items.append((param, value, param_norm))
        elif param_norm in vector_keys:
            vector_items.append((param, value, param_norm))
        else:
            other_items.append((param, value, param_norm))
    
    def _order_items(items: List[Tuple[str, Any, str]]) -> List[Tuple[str, Any, str]]:
        return sorted(items, key=lambda x: (x[2], str(x[0])))
    
    def _order_mode_items(items: List[Tuple[str, Any, str]]) -> List[Tuple[str, Any, str]]:
        priority = [
            "tc_grps_mode", "tc-grps-mode", "tc_grps", "tc-grps",
            "comm_grps_policy", "comm-grps-policy",
        ]
        by_key: Dict[str, List[Tuple[str, Any, str]]] = {}
        for item in items:
            by_key.setdefault(item[2], []).append(item)
        for key in by_key:
            by_key[key].sort(key=lambda x: str(x[0]))
        ordered = []
        for key in priority:
            ordered.extend(by_key.pop(key, []))
        for key in sorted(by_key.keys()):
            ordered.extend(by_key[key])
        return ordered

    def _order_scalar_items(items: List[Tuple[str, Any, str]]) -> List[Tuple[str, Any, str]]:
        priority = ["temperature", "pressure", "gen_temp", "gen-temp"]
        by_key: Dict[str, List[Tuple[str, Any, str]]] = {}
        for item in items:
            by_key.setdefault(item[2], []).append(item)
        for key in by_key:
            by_key[key].sort(key=lambda x: str(x[0]))
        ordered = []
        for key in priority:
            ordered.extend(by_key.pop(key, []))
        for key in sorted(by_key.keys()):
            ordered.extend(by_key[key])
        return ordered
    
    def _apply_override(param: str, value: Any, param_norm: str) -> OverrideOutcome:
        # Handle temperature override (affects multiple params)
        # Task C: Also inject missing ref_t/gen_temp when user explicitly sets temperature
        if param_norm == "temperature":
            any_temp_applied = False
            for temp_param in TEMPERATURE_PARAMS:
                temp_key = normalize_key(temp_param)
                if temp_key in patched_params:
                    original_values.setdefault(temp_key, patched_params[temp_key])
                    patched_params[temp_key] = str(value)
                    patched_values[temp_key] = str(value)
                    any_temp_applied = True
            
            # Task C: Inject ref_t if no variant exists (user explicitly wants this temperature)
            if find_existing_key(patched_params, "ref_t") is None:
                chosen_key = set_semantic_param(patched_params, original_params, "ref_t", str(value))
                patched_values[normalize_key(chosen_key)] = str(value)
                any_temp_applied = True
            
            # Task C: Inject gen_temp if gen_vel=yes and gen_temp is missing
            gen_vel_val = get_semantic_value(patched_params, "gen_vel")
            if gen_vel_val and gen_vel_val.lower() == "yes":
                if find_existing_key(patched_params, "gen_temp") is None:
                    chosen_key = set_semantic_param(patched_params, original_params, "gen_temp", str(value))
                    patched_values[normalize_key(chosen_key)] = str(value)
                    any_temp_applied = True

            if any_temp_applied:
                return OverrideOutcome(OverrideStatus.APPLIED)
            return OverrideOutcome(
                OverrideStatus.MISSING_IN_TEMPLATE,
                "temperature override could not target any template temperature parameter",
            )
        
        # Handle pressure override
        if param_norm == "pressure":
            applied = False
            for pres_param in PRESSURE_PARAMS:
                pres_key = normalize_key(pres_param)
                if pres_key in patched_params:
                    original_values.setdefault(pres_key, patched_params[pres_key])
                    patched_params[pres_key] = str(value)
                    patched_values[pres_key] = str(value)
                    applied = True
            if applied:
                return OverrideOutcome(OverrideStatus.APPLIED)
            return OverrideOutcome(
                OverrideStatus.MISSING_IN_TEMPLATE,
                "pressure override could not be applied because ref_p is absent from template",
            )
        
        # Handle nsteps override
        if param_norm == "nsteps":
            # Issue E: Validate nsteps > 0 for production stages
            try:
                nsteps_val = int(value)
            except (TypeError, ValueError) as exc:
                stage_name = stage or "<unknown>"
                raise MDPPatchError(
                    f"Stage '{stage_name}': invalid value for parameter 'nsteps': {value!r}. "
                    f"Expected an integer number of MD steps."
                ) from exc
            if stage in ["nvt", "npt", "md"] and nsteps_val <= 0:
                raise MDPPatchError(
                    f"nsteps must be > 0 for production stage '{stage}', got {nsteps_val}"
                )
            if "nsteps" in patched_params:
                original_values.setdefault("nsteps", patched_params["nsteps"])
            patched_params["nsteps"] = str(nsteps_val)
            patched_values["nsteps"] = str(nsteps_val)
            return OverrideOutcome(OverrideStatus.APPLIED)
        
        # Handle tc-grps mode (Issue A: safe split)
        if param_norm in {"tc_grps_mode", "tc-grps-mode"}:
            if str(value).lower() == "split":
                return _apply_tc_grps_split(
                    patched_params, original_values, patched_values, ctx,
                    original_params=original_params,
                    sink=sink,
                )
            return OverrideOutcome(OverrideStatus.APPLIED)
        
        # Handle comm-grps policy (Issue B)
        if param_norm in {"comm_grps_policy", "comm-grps-policy"}:
            return _apply_comm_grps_policy(
                patched_params, original_values, patched_values, value, ctx,
                original_params=original_params,
                sink=sink,
            )
        
        # Handle gen_vel override (Issue C) - Task A: use semantic key helpers
        if param_norm == "gen_vel":
            existing_key = find_existing_key(original_params, "gen_vel")
            if existing_key:
                original_values.setdefault(existing_key, patched_params.get(existing_key, ""))
            chosen_key = set_semantic_param(patched_params, original_params, "gen_vel", str(value))
            patched_values[normalize_key(chosen_key)] = str(value)
            return OverrideOutcome(OverrideStatus.APPLIED)
        
        # Handle gen_seed override for reproducibility - Task A: use semantic key helpers
        if param_norm in {"gen_seed", "gen-seed"}:
            existing_key = find_existing_key(original_params, "gen_seed")
            if existing_key:
                original_values.setdefault(existing_key, patched_params.get(existing_key, ""))
            chosen_key = set_semantic_param(patched_params, original_params, "gen_seed", str(value))
            patched_values[normalize_key(chosen_key)] = str(value)
            return OverrideOutcome(OverrideStatus.APPLIED)
        
        # Handle npt_barostat override (Berendsen vs Parrinello-Rahman)
        # Task B: Remove C-rescale (it's a thermostat, not a barostat)
        if param_norm == "npt_barostat":
            barostat = str(value)
            if barostat not in ["Berendsen", "Parrinello-Rahman"]:
                raise MDPPatchError(
                    f"Invalid npt_barostat value: {barostat}. "
                    f"Must be 'Berendsen' or 'Parrinello-Rahman'. "
                    f"Note: C-rescale is a thermostat (temperature coupling), not a barostat."
                )
            # Update pcoupl
            if "pcoupl" in patched_params:
                original_values.setdefault("pcoupl", patched_params["pcoupl"])
            patched_params["pcoupl"] = barostat
            patched_values["pcoupl"] = barostat
            
            # Issue 3: tau_p defaults with system-type awareness for gel stability
            # Check for both underscore and hyphen variants
            tau_p_existing = find_existing_key(original_params, "tau_p")
            if tau_p_existing is None:
                # No tau_p in template - set default based on barostat and system type
                if barostat == "Berendsen":
                    default_tau_p = TAU_P_DEFAULT_BERENDSEN
                else:
                    # Parrinello-Rahman: use system type to determine tau_p
                    # Gel/solid-like systems need larger damping to prevent oscillations
                    system_type = getattr(ctx, "system_type", None) if ctx else None
                    if system_type in ("gel", "solid"):
                        default_tau_p = TAU_P_DEFAULT_PR_GEL
                    elif system_type == "liquid":
                        default_tau_p = TAU_P_DEFAULT_PR_LIQUID
                    else:
                        # Unknown system type - use safer gel default
                        default_tau_p = TAU_P_DEFAULT_PR_GEL
                        sink.warn(
                            f"Parrinello-Rahman barostat with unknown system type: "
                            f"using safer tau_p={default_tau_p} ps (gel default). "
                            f"Use --system-type liquid for tau_p={TAU_P_DEFAULT_PR_LIQUID} ps."
                        )
                chosen_key = set_semantic_param(patched_params, original_params, "tau_p", default_tau_p)
                patched_values[normalize_key(chosen_key)] = default_tau_p
            else:
                # tau_p exists in template - preserve it with warning
                current_tau_p = original_params.get(tau_p_existing, "")
                expected_tau_p = TAU_P_DEFAULT_BERENDSEN if barostat == "Berendsen" else TAU_P_DEFAULT_PR_GEL
                sink.warn(
                    f"Preserving existing tau_p={current_tau_p} from template. "
                    f"Barostat '{barostat}' typically uses tau_p={expected_tau_p} ps. "
                    f"If simulation becomes unstable, consider adjusting tau_p in the template."
                )

            # PR compressibility safety: inject if missing; preserve if existing.
            if barostat == "Parrinello-Rahman":
                pcoupltype = _resolve_pcoupltype_for_dimensionality(patched_params, sink=sink)
                n_comp_values = PCOUPLTYPE_VALUE_COUNTS.get(pcoupltype, 1)
                system_type = (getattr(ctx, "system_type", None) if ctx else None) or "unknown"
                system_type_norm = str(system_type).strip().lower()

                comp_existing = find_existing_key(patched_params, "compressibility")
                if comp_existing is None:
                    if system_type_norm == "liquid":
                        comp_default = COMPRESSIBILITY_DEFAULT_LIQUID
                    else:
                        comp_default = COMPRESSIBILITY_DEFAULT_GEL
                    comp_value = _build_compressibility_vector(comp_default, n_comp_values)
                    chosen_comp_key = set_semantic_param(
                        patched_params, original_params, "compressibility", comp_value
                    )
                    patched_values[normalize_key(chosen_comp_key)] = comp_value
                    sink.warn(
                        "Parrinello-Rahman selected and template has no compressibility. "
                        f"Injected compressibility={comp_value} for pcoupltype={pcoupltype}. "
                        "For gels/solids this uses a safer low value to reduce volume ringing/LINCS risk.",
                        escalate_in_strict=False,
                    )
                else:
                    comp_value = patched_params.get(comp_existing, "")
                    if (
                        system_type_norm in {"gel", "solid"}
                        and _looks_water_like_compressibility(comp_value)
                    ):
                        sink.warn(
                            "Parrinello-Rahman with gel/solid system is using water-like "
                            f"compressibility ({comp_value}). This can cause severe volume ringing "
                            "and LINCS failures. Consider a smaller value (e.g., 1.0e-6) or provide "
                            "an explicit compressibility override.",
                            escalate_in_strict=False,
                        )
            return OverrideOutcome(OverrideStatus.APPLIED)
        
        # Direct parameter override - with semantic key awareness (Task A fix)
        # If param belongs to a semantic group, update the existing variant
        # to avoid creating duplicates that would be lost during dedupe.
        semantic_name = VARIANT_TO_SEMANTIC.get(param_norm)
        
        if semantic_name:
            # Find which variant currently exists in patched_params
            existing_key = find_existing_key(patched_params, semantic_name)
            if existing_key:
                # Update the EXISTING variant (preserve template style)
                original_values.setdefault(existing_key, patched_params[existing_key])
                patched_params[existing_key] = str(value)
                patched_values[normalize_key(existing_key)] = str(value)
                return OverrideOutcome(OverrideStatus.APPLIED)
            else:
                # No existing variant - use preferred (first in group)
                variants = SEMANTIC_KEY_GROUPS.get(semantic_name, (param_norm,))
                preferred_key = normalize_key(variants[0])
                patched_params[preferred_key] = str(value)
                patched_values[preferred_key] = str(value)
                return OverrideOutcome(OverrideStatus.APPLIED)
        else:
            # Not a semantic key - simple override with strict unknown-key handling
            if param_norm in original_params:
                original_values.setdefault(param_norm, patched_params.get(param_norm, ""))
            else:
                if strict:
                    raise MDPValidationError(
                        f"Unknown override key '{param}'. This key is not present in the template "
                        f"and is not a known semantic alias. Likely typo.\n"
                        f"Available params: {', '.join(sorted(original_params.keys()))}"
                    )
                sink.warn(
                    f"Override key '{param}' not found in template; adding new line. "
                    f"Check for typos if this was not intentional."
                )
            patched_params[param_norm] = str(value)
            patched_values[param_norm] = str(value)
            return OverrideOutcome(OverrideStatus.APPLIED)
    
    # Phase 2: apply mode-changing overrides first (tc-grps split/unified)
    for param, value, param_norm in _order_mode_items(mode_items):
        applied_overrides[param_norm] = _apply_override(param, value, param_norm)
    
    # Phase 3: apply scalar temperature-like overrides
    for param, value, param_norm in _order_scalar_items(scalar_items):
        applied_overrides[param_norm] = _apply_override(param, value, param_norm)
    
    # Phase 4: apply explicit vector overrides (ref_t/tau_t/etc.)
    for param, value, param_norm in _order_items(vector_items):
        applied_overrides[param_norm] = _apply_override(param, value, param_norm)
    
    # Phase 5: apply everything else
    for param, value, param_norm in _order_items(other_items):
        applied_overrides[param_norm] = _apply_override(param, value, param_norm)
    
    # Phase 6: final normalization pass after all overrides
    _normalize_tc_grps_vectors(
        patched_params,
        original_values,
        patched_values,
        ctx,
        original_params=original_params,
        sink=sink,
    )

    # Global safety guard for templates/overrides that end up with Angular under PBC.
    _guard_angular_comm_mode_with_pbc(
        patched_params,
        original_values,
        patched_values,
        original_params,
        ctx,
        sink,
        strict,
    )

    # Validate override outcomes: only true failures/missing params error.
    if validate:
        failed = []
        missing = []
        for key, outcome in applied_overrides.items():
            if outcome.status == OverrideStatus.FAILED:
                failed.append(f"{key}: {outcome.message or 'failed'}")
            elif outcome.status == OverrideStatus.MISSING_IN_TEMPLATE:
                missing.append(f"{key}: {outcome.message or 'not found in template'}")
            elif outcome.status == OverrideStatus.SKIPPED_WITH_WARNING and outcome.message:
                sink.warn(outcome.message, escalate_in_strict=False)
        if failed or missing:
            details = []
            if missing:
                details.append("Missing template parameters:\n  - " + "\n  - ".join(missing))
            if failed:
                details.append("Failed overrides:\n  - " + "\n  - ".join(failed))
            raise MDPValidationError(
                "One or more overrides could not be applied.\n"
                + "\n".join(details)
                + f"\nAvailable params: {', '.join(sorted(original_params.keys()))}"
            )
    
    # Write patched MDP
    write_mdp(
        output_path,
        patched_params,
        original_content,
        preserve_comments=True,
        patched_values=patched_values,
        strict=strict,
        sink=sink,
    )

    # Issue D: Post-write validation
    if validate:
        validation_errors = []
        validation_warnings = []
        
        # Vector length validation
        vec_errors, vec_warnings = _validate_vector_lengths(patched_params, stage, strict=strict)
        validation_errors.extend(vec_errors)
        validation_warnings.extend(vec_warnings)
        
        # gen_vel consistency validation
        gv_errors, gv_warnings = _validate_gen_vel_consistency(
            patched_params, stage, strict=strict, ctx=ctx
        )
        validation_errors.extend(gv_errors)
        validation_warnings.extend(gv_warnings)
        
        # Issue C: ref_p dimensionality validation
        refp_errors, refp_warnings = _validate_ref_p_dimensionality(patched_params, stage, strict=strict)
        validation_errors.extend(refp_errors)
        validation_warnings.extend(refp_warnings)

        comp_errors, comp_warnings = _validate_compressibility_dimensionality(
            patched_params, stage, strict=strict
        )
        validation_errors.extend(comp_errors)
        validation_warnings.extend(comp_warnings)

        comm_errors, comm_warnings = _validate_comm_mode_pbc_safety(
            patched_params, stage, strict=strict, ctx=ctx
        )
        validation_errors.extend(comm_errors)
        validation_warnings.extend(comm_warnings)
        
        # Checkpoint availability is tri-state: True/False/Unknown (None).
        effective_checkpoint = _resolve_checkpoint_known_state(
            checkpoint_available, stage=stage, sink=sink
        )
        gvc_errors, gvc_warnings = _validate_gen_vel_continuation_consistency(
            patched_params, stage, strict=strict, ctx=ctx,
            checkpoint_available=effective_checkpoint,
        )
        validation_errors.extend(gvc_errors)
        validation_warnings.extend(gvc_warnings)
        
        # Verify written values match expected
        written_params = parse_mdp(output_path)
        
        # Post-write semantic validation (independent of patched_values)
        for semantic_name, variants in SEMANTIC_KEY_GROUPS.items():
            existing = [normalize_key(v) for v in variants if normalize_key(v) in written_params]
            if len(existing) > 1:
                variant_values = {v: written_params[v] for v in existing}
                msg = (
                    f"Semantic key '{semantic_name}' has multiple variants in output: "
                    f"{variant_values}. Only one variant should remain."
                )
                if strict:
                    validation_errors.append(msg)
                else:
                    validation_warnings.append(msg)
        
        for key, expected in patched_values.items():
            key_norm = normalize_key(key)
            # Check both the key and its semantic variants
            actual = written_params.get(key_norm)
            if actual is None:
                # Try semantic key lookup
                semantic_name = VARIANT_TO_SEMANTIC.get(key_norm)
                if semantic_name:
                    for v in SEMANTIC_KEY_GROUPS.get(semantic_name, ()):
                        v_norm = normalize_key(v)
                        if v_norm in written_params:
                            actual = written_params.get(v_norm)
                            break
            if actual != expected:
                validation_errors.append(
                    f"Post-write mismatch: {key_norm} expected '{expected}', got '{actual}'"
                )
        
        # Collect warnings into sink (do not rely on Python warning filters)
        for w in validation_warnings:
            sink.warn(w, escalate_in_strict=False)

        if validation_errors:
            raise MDPValidationError(
                f"MDP validation failed after patching {output_path.name}:\n" +
                "\n".join(f"  - {e}" for e in validation_errors)
            )

    if validate and sink.errors:
        raise MDPValidationError(
            f"MDP diagnostics failed for {output_path.name}:\n"
            + "\n".join(f"  - {e}" for e in sink.errors)
        )

    if return_diagnostics:
        return original_values, patched_values, list(sink.errors), list(sink.warnings)
    return original_values, patched_values


def _apply_tc_grps_split(
    patched_params: Dict[str, str],
    original_values: Dict[str, str],
    patched_values: Dict[str, str],
    ctx: Optional["PipelineContext"],
    original_params: Optional[Dict[str, str]] = None,
    sink: Optional[DiagnosticSink] = None,
) -> OverrideOutcome:
    """
    Apply tc-grps split mode safely.
    
    Issue A: Do NOT blindly duplicate tau_t. Require explicit per-group values
    or explicit opt-in to autofill.
    Issue 2: Do NOT hardcode Polymer/NonPolymer - use configurable groups.
    Task A: Use semantic key helpers to preserve key style.
    
    Args:
        patched_params: Parameters being built
        original_values: Record of original values before patching
        patched_values: Record of patched values
        ctx: Pipeline context
        original_params: Original template parameters (for style detection)
        
    Returns:
        OverrideOutcome describing whether split was applied or intentionally skipped.
    """
    # Use original_params for style detection, fall back to patched_params
    style_ref = original_params if original_params else patched_params
    strict = getattr(ctx, "strict_mdp_validation", False) if ctx else False
    auto_mode_requested = (getattr(ctx, "tc_grps_mode", "system") == "auto") if ctx else False
    if sink is None:
        sink = DiagnosticSink(strict=strict)
    
    # Issue 2: Determine tc-grps group names
    # Priority: 1) explicit ctx.tc_grps_groups, 2) auto-detect from ndx, 3) legacy fallback.
    tc_grps_value = None
    split_group_names: List[str] = []
    split_source = "unknown"
    tc_grps_ndx_path = getattr(ctx, "tc_grps_ndx_path", None) if ctx else None
    available_groups: Optional[List[str]] = None
    ndx_group_sizes: Optional[Dict[str, int]] = None

    if tc_grps_ndx_path:
        try:
            ndx_path = Path(tc_grps_ndx_path)
            available_groups = parse_ndx_groups(ndx_path)
            ndx_group_sizes = parse_ndx_group_sizes(ndx_path)
        except MDPParseError as e:
            sink.warn(f"Failed to parse ndx for tc-grps: {e}", escalate_in_strict=False)

    tc_grps_groups = getattr(ctx, "tc_grps_groups", None) if ctx else None
    if tc_grps_groups:
        tc_grps_value = normalize_group_list(tc_grps_groups)
        split_group_names = tc_grps_value.split()
        split_source = "explicit"

    if not tc_grps_value and available_groups is not None:
        non_system_groups = [g for g in available_groups if not groups_equal(g, "System")]
        if len(non_system_groups) >= 2:
            split_group_names = non_system_groups[:2]
            tc_grps_value = " ".join(split_group_names)
            split_source = "auto_ndx"
            sink.warn(
                f"tc-grps split: auto-detected groups from {Path(tc_grps_ndx_path).name}: "
                f"{split_group_names}. Use --tc-grps-groups to override if incorrect.",
                escalate_in_strict=False,
            )
        else:
            sink.warn(
                f"tc-grps split: ndx file {Path(tc_grps_ndx_path).name} has insufficient groups for split "
                f"(found: {available_groups}). Need at least 2 non-System groups.",
                escalate_in_strict=False,
            )

    if not tc_grps_value:
        if strict:
            raise MDPPatchError(
                f"tc-grps split requires explicit group configuration in strict mode.\n"
                f"Solutions:\n"
                f"  1. Use --tc-grps-groups 'GroupA GroupB' with your index group names\n"
                f"  2. Use --tc-grps-ndx-path /path/to/index.ndx for auto-detection\n"
                f"  3. Disable strict mode to use legacy 'Polymer NonPolymer' default\n"
                f"Note: Hardcoded group names may fail if your .ndx uses different labels."
            )
        tc_grps_value = "Polymer NonPolymer"
        split_group_names = ["Polymer", "NonPolymer"]
        split_source = "legacy_default"
        sink.warn(
            "tc-grps split: using legacy default 'Polymer NonPolymer'. "
            "Ensure these groups exist in your .ndx file. "
            "For HTPolyNet/GPE workflows with different labels, use --tc-grps-groups 'YourGroup1 YourGroup2'.",
            escalate_in_strict=False,
        )
    
    n_groups = len(split_group_names)

    def _fallback_to_system_due_to_small_groups(reason: str) -> OverrideOutcome:
        sink.warn(
            f"{reason} Falling back to tc-grps=System for thermostat safety.",
            escalate_in_strict=False,
        )
        existing_tc_grps = find_existing_key(patched_params, "tc_grps")
        if existing_tc_grps:
            original_values.setdefault(existing_tc_grps, patched_params[existing_tc_grps])
        chosen_tc = set_semantic_param(patched_params, style_ref, "tc_grps", "System")
        patched_values[normalize_key(chosen_tc)] = "System"

        for semantic_name in ("tau_t", "ref_t"):
            vec_val = get_semantic_value(patched_params, semantic_name)
            if not vec_val:
                continue
            vec_tokens = vec_val.split()
            if len(vec_tokens) <= 1:
                continue
            existing_vec_key = find_existing_key(patched_params, semantic_name)
            if existing_vec_key:
                original_values.setdefault(existing_vec_key, patched_params[existing_vec_key])
            collapsed = vec_tokens[0]
            chosen_vec_key = set_semantic_param(patched_params, style_ref, semantic_name, collapsed)
            patched_values[normalize_key(chosen_vec_key)] = collapsed
            sink.warn(
                f"Collapsed {semantic_name} to '{collapsed}' because tc-grps split was disabled.",
                escalate_in_strict=False,
            )
        return OverrideOutcome(OverrideStatus.APPLIED)

    min_atoms = getattr(ctx, "tc_grps_min_atoms", 50) if ctx else 50
    min_fraction = getattr(ctx, "tc_grps_min_fraction", 0.01) if ctx else 0.01
    allow_small_tc_grps = getattr(ctx, "allow_small_tc_grps", False) if ctx else False
    try:
        min_atoms = int(min_atoms)
    except (TypeError, ValueError):
        min_atoms = 50
    if min_atoms < 1:
        sink.warn(
            f"Invalid tc_grps_min_atoms={min_atoms}; using minimum threshold 1.",
            escalate_in_strict=False,
        )
        min_atoms = 1
    try:
        min_fraction = float(min_fraction)
    except (TypeError, ValueError):
        min_fraction = 0.01
    if min_fraction <= 0.0 or min_fraction >= 1.0:
        sink.warn(
            f"Invalid tc_grps_min_fraction={min_fraction}; using 0.01.",
            escalate_in_strict=False,
        )
        min_fraction = 0.01

    if ndx_group_sizes is not None and split_group_names:
        casefold_lookup = {name.casefold(): (name, size) for name, size in ndx_group_sizes.items()}
        missing_size = [g for g in split_group_names if g.casefold() not in casefold_lookup]
        tiny = []
        system_item = casefold_lookup.get("system")
        system_size = int(system_item[1]) if system_item and int(system_item[1]) > 0 else 0
        if system_size <= 0:
            system_size = sum(max(0, int(size)) for _, size in casefold_lookup.values())
        for group in split_group_names:
            item = casefold_lookup.get(group.casefold())
            if item is None:
                continue
            canonical, size = item
            fraction = (float(size) / float(system_size)) if system_size > 0 else 0.0
            if size < min_atoms or fraction < min_fraction:
                tiny.append((canonical, size, fraction))

        if tiny:
            tiny_desc = ", ".join(
                f"{name}={size} ({fraction:.4f})" for name, size, fraction in tiny
            )
            risk_msg = (
                "tc-grps split selects tiny thermostat group(s) below thresholds "
                f"(min_atoms={min_atoms}, min_fraction={min_fraction:.4f}): "
                f"{tiny_desc}. Tiny groups amplify thermostat noise (low DOF) and can trigger "
                f"temperature spikes/instability."
            )
            fix_msg = (
                "Fixes: merge ions into a larger thermostat group, use tc-grps=System, "
                "increase --tc-grps-min-atoms, or opt in with --allow-small-tc-grps."
            )
            if split_source == "explicit":
                if allow_small_tc_grps:
                    sink.warn(
                        f"{risk_msg} Proceeding only because allow_small_tc_grps=True. {fix_msg}",
                        escalate_in_strict=False,
                    )
                elif strict:
                    raise MDPPatchError(f"{risk_msg} {fix_msg}")
                else:
                    return _fallback_to_system_due_to_small_groups(f"{risk_msg} {fix_msg}")
            else:
                return _fallback_to_system_due_to_small_groups(f"{risk_msg} {fix_msg}")

        if missing_size:
            sink.warn(
                f"Could not verify atom counts for tc-grps groups {missing_size} from ndx; "
                f"small-group DOF safety check is partial.",
                escalate_in_strict=False,
            )
    elif split_source == "explicit" and ndx_group_sizes is None:
        sink.warn(
            "tc-grps split groups were explicit but no readable ndx was available for atom-count checks. "
            "Small-group DOF safety could not be verified.",
            escalate_in_strict=False,
        )

    # Handle tau_t (Issue A: no blind duplication)
    existing_tau_t = find_existing_key(patched_params, "tau_t")
    current_tau_t = get_semantic_value(patched_params, "tau_t") or ""
    current_tau_t = current_tau_t.strip()
    n_current_tau = len(current_tau_t.split()) if current_tau_t else 0
    pending_tau_t: Optional[str] = None
    
    # Custom per-group tau_t values (any number of groups)
    tau_t_values = getattr(ctx, "tc_grps_tau_t_values", None) if ctx else None
    if tau_t_values:
        if len(tau_t_values) != n_groups:
            raise MDPPatchError(
                f"tc-grps split expects {n_groups} tau_t values but got {len(tau_t_values)}. "
                f"Groups: {split_group_names}"
            )
        pending_tau_t = " ".join(str(v) for v in tau_t_values)
    else:
        # Check if explicit per-group tau_t provided (works for 2 groups)
        tau_t_polymer = getattr(ctx, "tau_t_polymer", None) if ctx else None
        tau_t_nonpolymer = getattr(ctx, "tau_t_nonpolymer", None) if ctx else None
        allow_autofill = getattr(ctx, "allow_tau_t_autofill", False) if ctx else False
        
        if tau_t_polymer is not None and tau_t_nonpolymer is not None and n_groups == 2:
            # Explicit per-group tau_t provided - use them
            pending_tau_t = f"{tau_t_polymer} {tau_t_nonpolymer}"
        elif n_current_tau == 1:
            # Single tau_t value - need to expand for n groups
            if allow_autofill:
                # Explicit opt-in: duplicate with warning
                single_val = current_tau_t.split()[0]
                pending_tau_t = " ".join([single_val] * n_groups)
                sink.warn(
                    f"tc-grps split with allow_tau_t_autofill=True: "
                    f"duplicating tau_t={single_val} for {n_groups} groups ({split_group_names}). "
                    f"This may cause non-ergodic energy partitioning in GPE systems. "
                    f"Consider providing explicit per-group tau_t values."
                )
            else:
                err_msg = (
                    f"tc-grps split requires tau_t for each of {n_groups} groups ({split_group_names}), "
                    f"but template has single value tau_t={current_tau_t}.\n"
                    f"Solutions:\n"
                    f"  1. Provide --tau-t-values <list> matching your tc-grps groups\n"
                    f"  2. Provide --tau-t-polymer <ps> --tau-t-nonpolymer <ps> (for 2 groups)\n"
                    f"  3. Use --allow-tau-t-autofill to duplicate (not recommended)\n"
                    f"  4. Edit the MDP template to have {n_groups} tau_t values\n"
                    f"Note: Different groups typically have different dominant timescales; "
                    f"identical tau_t can cause hidden temperature partitioning."
                )
                if auto_mode_requested:
                    return _fallback_to_system_due_to_small_groups(
                        err_msg + " Auto mode requires a safe split; reverting to System."
                    )
                raise MDPPatchError(err_msg)
        elif n_current_tau == n_groups:
            # Already has correct number of values - keep them
            pass
        elif n_current_tau > 0 and n_current_tau != n_groups:
            # Issue 2: Validate tau_t/ref_t count matches tc-grps count
            if auto_mode_requested:
                return _fallback_to_system_due_to_small_groups(
                    f"tc-grps split vector mismatch (tau_t has {n_current_tau}, groups={n_groups}). "
                    "Auto mode requires aligned vectors."
                )
            if strict:
                raise MDPPatchError(
                    f"tau_t has {n_current_tau} values but tc-grps split uses {n_groups} groups "
                    f"({split_group_names}). Please fix the MDP template."
                )
            else:
                return OverrideOutcome(
                    OverrideStatus.SKIPPED_WITH_WARNING,
                    f"tc-grps split skipped: tau_t has {n_current_tau} values but split needs "
                    f"{n_groups} groups ({split_group_names}). "
                    f"Fix template vector lengths, provide --tau-t-values, or enable strict mode to fail fast.",
                )
        # Issue E: n_current_tau == 0 - no tau_t in template
        # Check if tcoupl is enabled; if so, this will fail at grompp
        if n_current_tau == 0:
            tcoupl = get_semantic_value(patched_params, "tcoupl")
            if tcoupl and tcoupl.lower() != "no":
                msg = (
                    f"tc-grps split with tcoupl={tcoupl} but tau_t is missing.\n"
                    f"  Groups: {split_group_names}\n"
                    f"  This will cause grompp to fail.\n"
                    f"  Solution: Add tau_t values to the MDP template or use --tau-t-values."
                )
                if auto_mode_requested:
                    return _fallback_to_system_due_to_small_groups(
                        msg + " Auto mode requires complete thermostat vectors."
                    )
                if strict:
                    raise MDPPatchError(msg)
                else:
                    sink.warn(f"[WARNING] {msg}", escalate_in_strict=False)

    # Apply tc-grps only after all skip/fail checks to avoid partial state changes.
    existing_tc_grps = find_existing_key(patched_params, "tc_grps")
    if existing_tc_grps:
        original_values[existing_tc_grps] = patched_params[existing_tc_grps]
    chosen_tc_key = set_semantic_param(patched_params, style_ref, "tc_grps", tc_grps_value)
    patched_values[normalize_key(chosen_tc_key)] = tc_grps_value

    if pending_tau_t is not None:
        if existing_tau_t:
            original_values[existing_tau_t] = current_tau_t
        chosen_tau_key = set_semantic_param(patched_params, style_ref, "tau_t", pending_tau_t)
        patched_values[normalize_key(chosen_tau_key)] = pending_tau_t

    # NOTE: ref_t replication is intentionally deferred to _normalize_tc_grps_vectors().
    return OverrideOutcome(OverrideStatus.APPLIED)


def _resolve_available_ndx_groups(
    ctx: Optional["PipelineContext"],
    sink: Optional[DiagnosticSink] = None,
) -> Optional[List[str]]:
    """Resolve available index group names from context-provided ndx path."""
    ndx_path_raw = None
    if ctx is not None:
        ndx_path_raw = getattr(ctx, "active_ndx_path", None) or getattr(ctx, "tc_grps_ndx_path", None)
    if not ndx_path_raw:
        return None
    try:
        return parse_ndx_groups(Path(str(ndx_path_raw)))
    except MDPParseError as e:
        if sink is not None:
            sink.warn(f"Could not parse index groups from '{ndx_path_raw}': {e}", escalate_in_strict=False)
        return None


def _validate_comm_groups_exist(
    group_tokens: List[str],
    available_ndx_groups: Optional[List[str]],
    strict: bool,
    sink: Optional[DiagnosticSink],
) -> None:
    """Validate that comm-grps tokens exist in ndx when comm-grps is multi-group."""
    if len(group_tokens) <= 1:
        return
    if available_ndx_groups is None:
        if sink is not None:
            sink.warn(
                f"Unable to verify comm-grps groups {group_tokens} because no readable .ndx path was provided. "
                f"Provide --tc-grps-ndx-path or ensure stage injects active ndx path.",
                escalate_in_strict=False,
            )
        return
    available_lookup = {g.casefold() for g in available_ndx_groups}
    missing = [g for g in group_tokens if g.casefold() not in available_lookup]
    if not missing:
        return
    msg = (
        f"comm-grps groups missing in index file: {missing}. "
        f"Available groups: {available_ndx_groups}. "
        f"Remediation: update tc/comm groups to existing labels or supply an ndx containing these groups."
    )
    if strict:
        raise MDPPatchError(msg)
    if sink is not None:
        sink.warn(msg, escalate_in_strict=False)


def _guard_angular_comm_mode_with_pbc(
    patched_params: Dict[str, str],
    original_values: Dict[str, str],
    patched_values: Dict[str, str],
    style_ref: Dict[str, str],
    ctx: Optional["PipelineContext"],
    sink: Optional[DiagnosticSink],
    strict: bool,
) -> None:
    """
    Enforce comm-mode safety: Angular is unsafe under PBC by default.
    """
    comm_mode = get_semantic_value(patched_params, "comm_mode")
    if not comm_mode or comm_mode.strip().lower() != "angular":
        return
    if not _is_pbc_effectively_on(patched_params):
        return

    allow_angular = getattr(ctx, "allow_comm_mode_angular_with_pbc", False) if ctx else False
    pbc_value = patched_params.get("pbc", "xyz")
    if allow_angular:
        if sink is not None:
            sink.warn(
                "HIGH VISIBILITY WARNING: comm-mode=Angular with active PBC "
                f"(pbc='{pbc_value}') was explicitly allowed via "
                "allow_comm_mode_angular_with_pbc=True. This setting can produce "
                "incorrect results; validate carefully.",
                escalate_in_strict=False,
            )
        return

    msg = (
        f"comm-mode=Angular with active PBC (pbc='{pbc_value}') is unsafe/ill-defined and may "
        f"produce incorrect dynamics."
    )
    if strict:
        raise MDPPatchError(
            f"{msg} Strict mode forbids this combination. "
            f"Use comm-mode=Linear, disable PBC, or enable allow_comm_mode_angular_with_pbc for expert override."
        )

    existing_comm_mode = find_existing_key(patched_params, "comm_mode")
    if existing_comm_mode:
        original_values.setdefault(existing_comm_mode, patched_params[existing_comm_mode])
    chosen_mode_key = set_semantic_param(patched_params, style_ref, "comm_mode", "Linear")
    patched_values[normalize_key(chosen_mode_key)] = "Linear"
    if sink is not None:
        sink.warn(
            f"{msg} Overriding comm-mode to Linear by default for safety.",
            escalate_in_strict=False,
        )


def _apply_comm_grps_policy(
    patched_params: Dict[str, str],
    original_values: Dict[str, str],
    patched_values: Dict[str, str],
    policy: str,
    ctx: Optional["PipelineContext"],
    original_params: Optional[Dict[str, str]] = None,
    sink: Optional[DiagnosticSink] = None,
) -> OverrideOutcome:
    """
    Apply comm-grps policy for COM removal.
    
    Issue B: When tc-grps is split, global COM removal can cause artifacts.
    Task A: Use semantic key helpers to preserve key style.
    Task D: Add warning about ndx group availability.
    
    Args:
        patched_params: Parameters being built
        original_values: Record of original values
        patched_values: Record of patched values
        policy: Policy name (match_tc, system, explicit)
        ctx: Pipeline context
        original_params: Original template parameters (for style detection)
        
    Returns:
        OverrideOutcome describing whether policy was applied or failed.
    """
    strict = getattr(ctx, "strict_mdp_validation", False) if ctx else False
    comm_mode = getattr(ctx, "comm_mode", "Linear") if ctx else "Linear"
    style_ref = original_params if original_params else patched_params
    if sink is None:
        sink = DiagnosticSink(strict=strict)

    requested_policy = normalize_key(str(policy))
    tc_grps = get_semantic_value(patched_params, "tc_grps")
    percolated_risk = _has_percolated_risk(tc_grps)
    available_ndx_groups = _resolve_available_ndx_groups(ctx, sink=sink)

    # auto: choose safer system-wide COM removal when percolated risk is detected.
    # Rationale: system removes only overall drift and avoids per-group COM on
    # potentially ill-defined periodic network fragments.
    effective_policy = requested_policy
    if requested_policy == "auto":
        if percolated_risk:
            effective_policy = "system"
            sink.warn(
                "comm-grps-policy=auto detected polymer/gel/network/matrix thermostat groups; "
                "using system-wide COM removal (comm-grps=System).",
                escalate_in_strict=False,
            )
        else:
            effective_policy = "match_tc"
            sink.info("comm-grps-policy=auto selected match_tc (no percolated-risk keywords detected).")

    if effective_policy == "match_tc":
        # Match comm-grps to tc-grps
        if not tc_grps:
            return OverrideOutcome(
                OverrideStatus.MISSING_IN_TEMPLATE,
                "comm-grps-policy=match_tc requested but tc-grps is not set in template.",
            )
        existing_comm_grps = find_existing_key(patched_params, "comm_grps")
        if existing_comm_grps:
            original_values[existing_comm_grps] = patched_params[existing_comm_grps]
        existing_comm_mode = find_existing_key(patched_params, "comm_mode")
        if existing_comm_mode:
            original_values[existing_comm_mode] = patched_params[existing_comm_mode]
        chosen_grps_key = set_semantic_param(patched_params, style_ref, "comm_grps", tc_grps)
        patched_values[normalize_key(chosen_grps_key)] = tc_grps
        chosen_mode_key = set_semantic_param(patched_params, style_ref, "comm_mode", comm_mode)
        patched_values[normalize_key(chosen_mode_key)] = comm_mode
        _guard_angular_comm_mode_with_pbc(
            patched_params, original_values, patched_values, style_ref, ctx, sink, strict
        )

        group_names = _split_group_tokens(tc_grps)
        if percolated_risk:
            sink.warn(
                "STRONG WARNING: comm-grps-policy=match_tc applies per-group COM removal to "
                f"groups {group_names} containing polymer/gel/network/matrix-like names. "
                "For percolated periodic networks this COM can be ill-defined and may inject "
                "nonphysical forces. Prefer --comm-grps-policy system (or none for expert use).",
                escalate_in_strict=False,
            )
        _validate_comm_groups_exist(group_names, available_ndx_groups, strict, sink)
        return OverrideOutcome(OverrideStatus.APPLIED)

    if effective_policy == "system":
        existing_comm_grps = find_existing_key(patched_params, "comm_grps")
        if existing_comm_grps:
            original_values[existing_comm_grps] = patched_params[existing_comm_grps]
        existing_comm_mode = find_existing_key(patched_params, "comm_mode")
        if existing_comm_mode:
            original_values[existing_comm_mode] = patched_params[existing_comm_mode]
        chosen_grps_key = set_semantic_param(patched_params, style_ref, "comm_grps", "System")
        patched_values[normalize_key(chosen_grps_key)] = "System"
        chosen_mode_key = set_semantic_param(patched_params, style_ref, "comm_mode", comm_mode)
        patched_values[normalize_key(chosen_mode_key)] = comm_mode
        _guard_angular_comm_mode_with_pbc(
            patched_params, original_values, patched_values, style_ref, ctx, sink, strict
        )
        sink.warn(
            "comm-grps-policy=system removes overall center-of-mass drift only. "
            "Relative drift between sub-groups can still exist, but this is usually the safer default.",
            escalate_in_strict=False,
        )
        return OverrideOutcome(OverrideStatus.APPLIED)

    if effective_policy == "none":
        existing_comm_grps = find_existing_key(patched_params, "comm_grps")
        if existing_comm_grps:
            original_values[existing_comm_grps] = patched_params[existing_comm_grps]
        existing_comm_mode = find_existing_key(patched_params, "comm_mode")
        if existing_comm_mode:
            original_values[existing_comm_mode] = patched_params[existing_comm_mode]
        chosen_grps_key = set_semantic_param(patched_params, style_ref, "comm_grps", "System")
        patched_values[normalize_key(chosen_grps_key)] = "System"
        chosen_mode_key = set_semantic_param(patched_params, style_ref, "comm_mode", "None")
        patched_values[normalize_key(chosen_mode_key)] = "None"
        sink.warn(
            "comm-grps-policy=none disables COM removal (comm-mode=None; comm-grps set to System for clarity). "
            "Total/system drift may accumulate; monitor trajectory drift explicitly.",
            escalate_in_strict=False,
        )
        return OverrideOutcome(OverrideStatus.APPLIED)

    if effective_policy == "explicit":
        # Use explicitly provided comm-grps
        explicit_grps = getattr(ctx, "comm_grps_explicit", None) if ctx else None
        if explicit_grps:
            existing_comm_grps = find_existing_key(patched_params, "comm_grps")
            if existing_comm_grps:
                original_values[existing_comm_grps] = patched_params[existing_comm_grps]
            chosen_grps_key = set_semantic_param(patched_params, style_ref, "comm_grps", explicit_grps)
            patched_values[normalize_key(chosen_grps_key)] = explicit_grps
            chosen_mode_key = set_semantic_param(patched_params, style_ref, "comm_mode", comm_mode)
            patched_values[normalize_key(chosen_mode_key)] = comm_mode
            _guard_angular_comm_mode_with_pbc(
                patched_params, original_values, patched_values, style_ref, ctx, sink, strict
            )
            _validate_comm_groups_exist(_split_group_tokens(explicit_grps), available_ndx_groups, strict, sink)
            return OverrideOutcome(OverrideStatus.APPLIED)
        else:
            raise MDPPatchError(
                f"comm-grps-policy=explicit but no --comm-grps provided"
            )

    return OverrideOutcome(
        OverrideStatus.FAILED,
        f"Unsupported comm-grps-policy '{policy}'. Supported values: auto, match_tc, system, none, explicit.",
    )


def validate_mdp_consistency(
    mdp_path: Path,
    expected_values: Dict[str, str],
    validate_vectors: bool = True,
    stage: str = "",
    strict: bool = False,
    ctx: Optional["PipelineContext"] = None,
    checkpoint_available: Optional[bool] = None,
) -> Tuple[List[str], List[str]]:
    """
    Validate that MDP file contains expected values and is internally consistent.
    
    Task G: Now returns (errors, warnings) tuple for consistent handling.
    
    Args:
        mdp_path: Path to MDP file
        expected_values: Expected parameter values
        validate_vectors: If True, validate vector length consistency
        stage: Stage name for context-aware validation
        strict: If True, escalate warnings to errors
        ctx: Pipeline context (optional) for advanced validation toggles
        checkpoint_available: If True, a checkpoint file exists for this stage.
        
    Returns:
        Tuple of (errors, warnings) - lists of error/warning messages
    """
    errors = []
    warn_list = []
    actual_params = parse_mdp(mdp_path)
    
    for param, expected in expected_values.items():
        # Use semantic key lookup for param matching
        param_norm = normalize_key(param)
        actual = actual_params.get(param_norm)
        if actual is None:
            # Try semantic variant lookup
            semantic_name = VARIANT_TO_SEMANTIC.get(param_norm)
            if semantic_name:
                for v in SEMANTIC_KEY_GROUPS.get(semantic_name, ()):
                    v_norm = normalize_key(v)
                    if v_norm in actual_params:
                        actual = actual_params.get(v_norm)
                        break
        if actual is None:
            errors.append(f"Parameter '{param}' not found in {mdp_path.name}")
        elif str(actual).strip() != str(expected).strip():
            errors.append(
                f"Parameter '{param}' mismatch in {mdp_path.name}: "
                f"expected '{expected}', got '{actual}'"
            )
    
    if validate_vectors:
        vec_errors, vec_warnings = _validate_vector_lengths(actual_params, stage, strict=strict)
        errors.extend(vec_errors)
        warn_list.extend(vec_warnings)
    
    # Validate gen_vel consistency
    gv_errors, gv_warnings = _validate_gen_vel_consistency(
        actual_params, stage, strict=strict, ctx=ctx
    )
    errors.extend(gv_errors)
    warn_list.extend(gv_warnings)
    
    # Issue C: ref_p dimensionality validation
    refp_errors, refp_warnings = _validate_ref_p_dimensionality(actual_params, stage, strict=strict)
    errors.extend(refp_errors)
    warn_list.extend(refp_warnings)

    comp_errors, comp_warnings = _validate_compressibility_dimensionality(
        actual_params, stage, strict=strict
    )
    errors.extend(comp_errors)
    warn_list.extend(comp_warnings)

    comm_errors, comm_warnings = _validate_comm_mode_pbc_safety(
        actual_params, stage, strict=strict, ctx=ctx
    )
    errors.extend(comm_errors)
    warn_list.extend(comm_warnings)
    
    # P1.5: gen_vel/continuation consistency with checkpoint awareness
    local_sink = DiagnosticSink(strict=strict, stage=stage)
    effective_checkpoint = _resolve_checkpoint_known_state(
        checkpoint_available, stage=stage, sink=local_sink
    )
    gvc_errors, gvc_warnings = _validate_gen_vel_continuation_consistency(
        actual_params, stage, strict=strict, ctx=ctx, checkpoint_available=effective_checkpoint
    )
    errors.extend(gvc_errors)
    warn_list.extend(gvc_warnings)
    warn_list.extend(local_sink.warnings)
    if strict:
        errors.extend(local_sink.errors)
    else:
        warn_list.extend(local_sink.errors)
    
    return errors, warn_list


def get_mdp_overrides_for_stage(ctx: "PipelineContext", stage: str) -> Dict[str, Any]:
    """
    Get MDP overrides for a specific stage based on context.
    
    P0.1: Close the stage policy loop for BOTH continuation AND gen_vel.
    
    Override precedence (highest → lowest):
      1. Explicit user force flags (force_gen_vel)
      2. Expert allow flags (allow_continuation_override, allow_gro_velocity_restart)
      3. Stage defaults (NVT=fresh start, NPT/MD=continuation)
    
    Args:
        ctx: Pipeline context with CLI options
        stage: Stage name (em, nvt, npt, md)
        
    Returns:
        Dictionary of overrides to apply
    """
    overrides = {}
    
    # Temperature override (applies to NVT, NPT, MD)
    if ctx.temperature is not None and stage in ["nvt", "npt", "md"]:
        overrides["temperature"] = ctx.temperature
    
    # Pressure override (applies to NPT, MD)
    if ctx.pressure is not None and stage in ["npt", "md"]:
        overrides["pressure"] = ctx.pressure
    
    # nsteps override
    nsteps_map = {
        "em": ctx.nsteps_em,
        "nvt": ctx.nsteps_nvt,
        "npt": ctx.nsteps_npt,
        "md": ctx.nsteps_md,
    }
    if stage in nsteps_map and nsteps_map[stage] is not None:
        overrides["nsteps"] = nsteps_map[stage]
    
    # tc-grps mode (applies to NVT, NPT, MD)
    if stage in ["nvt", "npt", "md"]:
        requested_mode = getattr(ctx, "tc_grps_mode", "auto")
        if requested_mode == "split":
            overrides["tc_grps_mode"] = "split"
            overrides["comm_grps_policy"] = ctx.comm_grps_policy
            setattr(
                ctx,
                "tc_grps_auto_decision",
                {
                    "requested_mode": "split",
                    "applied_mode": "split",
                    "reason": "explicit_split_requested",
                    "warnings": [],
                },
            )
        elif requested_mode == "auto":
            auto_decision = _decide_auto_tc_grps_mode(ctx)
            setattr(ctx, "tc_grps_auto_decision", auto_decision)
            for warning in auto_decision.get("warnings", []):
                print(f"  [MDP WARNING] {warning}")
            if auto_decision.get("applied_mode") == "split":
                overrides["tc_grps_mode"] = "split"
                overrides["comm_grps_policy"] = ctx.comm_grps_policy
            else:
                # Enforce deterministic fallback regardless of template defaults.
                overrides["tc_grps_mode"] = "system"
        else:
            setattr(
                ctx,
                "tc_grps_auto_decision",
                {
                    "requested_mode": requested_mode,
                    "applied_mode": "system",
                    "reason": "explicit_system_requested",
                    "warnings": [],
                },
            )
    
    # =========================================================================
    # P0.1: Stage policy for continuation AND gen_vel
    # Override precedence:
    #   1. force_gen_vel (explicit force flag) → gen_vel=yes
    #   2. allow_continuation_override → skip auto-setting continuation
    #   3. Stage defaults (NVT=fresh, NPT/MD=continuation)
    # =========================================================================
    force_gen_vel = getattr(ctx, "force_gen_vel", False)
    allow_continuation_override = getattr(ctx, "allow_continuation_override", False)
    allow_gro_velocity_restart = getattr(ctx, "allow_gro_velocity_restart", False)
    
    if stage == "nvt":
        # NVT: Fresh start by default
        # continuation=no, gen_vel=yes
        if not allow_continuation_override:
            overrides["continuation"] = "no"
        if force_gen_vel:
            overrides["gen_vel"] = "yes"
        else:
            # Stage default for NVT: gen_vel=yes (generate velocities)
            overrides["gen_vel"] = "yes"
            
    elif stage in ["npt", "md"]:
        # NPT/MD: Continuation from previous stage by default
        # continuation=yes, gen_vel=no
        # Expert restart intent: when allow_gro_velocity_restart is enabled,
        # prefer continuation=no unless user explicitly bypasses auto-continuation.
        if force_gen_vel:
            # Issue D: force_gen_vel MUST set continuation=no to avoid contradiction
            # gen_vel=yes + continuation=yes is always invalid
            overrides["gen_vel"] = "yes"
            overrides["continuation"] = "no"
            # Ensure gen_temp is set for velocity generation
            temperature = getattr(ctx, "temperature", None)
            if temperature is not None:
                overrides["gen_temp"] = str(temperature)
        else:
            # Stage default: continuation=yes, gen_vel=no
            # Expert restart mode shifts default continuation to "no".
            if not allow_continuation_override:
                overrides["continuation"] = "no" if allow_gro_velocity_restart else "yes"
            overrides["gen_vel"] = "no"
    
    # gen_seed override for reproducibility (applies to NVT where gen_vel=yes)
    gen_seed = getattr(ctx, "gen_seed", None)
    if gen_seed is not None and stage == "nvt":
        overrides["gen_seed"] = gen_seed
    
    # npt_barostat override (applies to NPT only, not production MD)
    npt_barostat = getattr(ctx, "npt_barostat", None)
    if npt_barostat is not None and stage == "npt":
        overrides["npt_barostat"] = npt_barostat
    
    return overrides


def resolve_mdp_overrides(
    ctx: "PipelineContext",
    stage: str,
    extra_overrides: Optional[Dict[str, Any]] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    Resolve MDP overrides with explicit precedence.

    Precedence: stage_overrides < extra_overrides

    Returns:
        Tuple of (stage_overrides, extra_overrides, effective_overrides)
    """
    stage_overrides = get_mdp_overrides_for_stage(ctx, stage)
    extra = dict(extra_overrides) if extra_overrides else {}
    effective = dict(stage_overrides)
    effective.update(extra)
    return stage_overrides, extra, effective


def prepare_stage_mdp(
    ctx: "PipelineContext",
    stage: str,
    output_dir: Path,
    extra_overrides: Optional[Dict[str, Any]] = None,
    effective_overrides: Optional[Dict[str, Any]] = None,
    checkpoint_available: Optional[bool] = None,
) -> Tuple[Path, Dict[str, str], Dict[str, str]]:
    """
    Prepare MDP file for a GROMACS stage.
    
    Copies template from IN/, patches with CLI overrides, writes to stage dir.
    
    Args:
        ctx: Pipeline context
        stage: Stage name (em, nvt, npt, md)
        output_dir: Stage output directory
        extra_overrides: Optional overrides to merge after stage defaults
        checkpoint_available: If True, a checkpoint file exists for this stage.
            Passed to validation to enable checkpoint-aware consistency checks.
        
    Returns:
        Tuple of (patched_mdp_path, original_values, patched_values)
        
    Raises:
        MDPPatchError: If template not found or patching fails
        MDPValidationError: If validation fails (Issue D)
    """
    # Template path
    template_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "mdp")
    template_path = template_dir / f"{stage}.mdp"
    
    # Output path
    output_path = output_dir / f"{stage}.mdp"
    
    # Get overrides
    if effective_overrides is None:
        _, _, overrides = resolve_mdp_overrides(ctx, stage, extra_overrides=extra_overrides)
    else:
        overrides = dict(effective_overrides)
    
    if not template_path.exists():
        raise MDPPatchError(
            f"MDP template not found: {template_path}\n"
            f"Please create the template file before running the pipeline."
        )
    
    # Copy and patch (with validation) - Issue A: pass checkpoint_available
    original_values, patched_values, patch_errors, patch_warnings = copy_and_patch_mdp(
        template_path, output_path, overrides,
        validate=True, ctx=ctx, stage=stage,
        checkpoint_available=checkpoint_available,
        return_diagnostics=True,
    )
    
    # Task G: Final validation - ALWAYS run validate_mdp_consistency
    # even when no patches were applied (to catch template inconsistencies)
    strict = getattr(ctx, "strict_mdp_validation", False)
    validation_errors, validation_warnings = validate_mdp_consistency(
        output_path,
        patched_values,
        validate_vectors=True,
        stage=stage,
        strict=strict,
        ctx=ctx,
        checkpoint_available=checkpoint_available,
    )

    # Copy/patch and final validation both return explicit diagnostics lists.
    all_errors = list(patch_errors) + list(validation_errors)
    all_warnings = list(patch_warnings) + list(validation_warnings)

    for w in all_warnings:
        print(f"  [MDP WARNING] {w}")

    if all_errors:
        raise MDPValidationError(
            f"Final MDP validation failed for {stage}:\n" +
            "\n".join(f"  - {e}" for e in all_errors)
        )
    
    # Record in manifest
    if ctx.manifest:
        ctx.manifest.set_mdp_patch(
            stage,
            original_values=original_values,
            patched_values=patched_values,
        )
    
    return output_path, original_values, patched_values
