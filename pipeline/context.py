"""
PipelineContext: Holds runtime context for pipeline execution.
"""

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

from .path_manager import PathManager
from .manifest import ManifestWriter


SUPPORTED_FF_CHARGE_PAIRS = {
    ("GAFF2", "CM5"),
    ("OPLS-AA", "RESP"),
}
FF_CHARGE_COMPAT_DOC_SECTION = "AGENT_SKILLS/03_ASSET_RESOLUTION.md#compatibility-matrix"


MISSING = object()


@dataclass(frozen=True)
class ContextArgSpec:
    ctx_field: str
    arg_name: str
    caster: Callable[[str, Any], Any]
    default_value: Any = MISSING
    numeric_or_bool: bool = False


def _type_error(name: str, expected: str, value: Any) -> ValueError:
    return ValueError(
        f"Invalid value for '{name}': expected {expected}, got {value!r} "
        f"(type={type(value).__name__})."
    )


def _as_required_str(name: str, value: Any) -> str:
    if value is None:
        raise ValueError(f"Missing required argument '{name}'.")
    if not isinstance(value, str):
        raise _type_error(name, "string", value)
    stripped = value.strip()
    if not stripped:
        raise ValueError(f"Argument '{name}' cannot be empty.")
    return stripped


def _as_optional_str(name: str, value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, str):
        return value
    raise _type_error(name, "string or null", value)


def _as_bool(name: str, value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "on"}:
            return True
        if normalized in {"false", "0", "no", "off"}:
            return False
    raise _type_error(name, "bool or one of true/false/1/0/yes/no/on/off", value)


def _as_int(name: str, value: Any, *, min_value: Optional[int] = None) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, bool):
        raise _type_error(name, "int (bool is not allowed)", value)
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        text = value.strip()
        if not text:
            raise _type_error(name, "integer string", value)
        if text[0] in {"+", "-"}:
            digits = text[1:]
        else:
            digits = text
        if not digits or not digits.isdigit():
            raise _type_error(name, "integer string (e.g. '8')", value)
        parsed = int(text, 10)
    else:
        raise _type_error(name, "int or integer string", value)

    if min_value is not None and parsed < min_value:
        raise ValueError(
            f"Invalid value for '{name}': expected >= {min_value}, got {parsed}."
        )
    return parsed


def _as_float(name: str, value: Any, *, min_value: Optional[float] = None) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, bool):
        raise _type_error(name, "float (bool is not allowed)", value)
    if isinstance(value, (int, float)):
        parsed = float(value)
    elif isinstance(value, str):
        text = value.strip()
        if not text:
            raise _type_error(name, "numeric string", value)
        try:
            parsed = float(text)
        except ValueError as exc:
            raise _type_error(name, "float or numeric string (e.g. '0.2')", value) from exc
    else:
        raise _type_error(name, "float/int or numeric string", value)

    if not math.isfinite(parsed):
        raise ValueError(f"Invalid value for '{name}': expected finite float, got {value!r}.")
    if min_value is not None and parsed < min_value:
        raise ValueError(
            f"Invalid value for '{name}': expected >= {min_value}, got {parsed}."
        )
    return parsed


def _as_optional_path(name: str, value: Any) -> Optional[Path]:
    if value is None:
        return None
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        return Path(stripped)
    raise _type_error(name, "Path, string path, or null", value)


def _as_optional_float_list(name: str, value: Any) -> Optional[list]:
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        tokens = [token.strip() for token in stripped.split(",") if token.strip()]
    elif isinstance(value, (list, tuple)):
        tokens = list(value)
    else:
        raise _type_error(name, "comma-separated string or list of floats", value)
    if not tokens:
        return None
    parsed: list = []
    for idx, token in enumerate(tokens):
        num = _as_float(f"{name}[{idx}]", token)
        if num is None:
            raise _type_error(f"{name}[{idx}]", "float", token)
        if num <= 0:
            raise ValueError(
                f"Invalid value for '{name}[{idx}]': expected > 0, got {num}."
            )
        parsed.append(num)
    return parsed


def _as_optional_threshold_dict(name: str, value: Any) -> Optional[dict]:
    if value is None:
        return None

    raw_dict: Any = None
    if isinstance(value, dict):
        raw_dict = value
    elif isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        try:
            raw_dict = json.loads(stripped)
        except json.JSONDecodeError:
            pairs = [part.strip() for part in stripped.split(",") if part.strip()]
            parsed_pairs: dict[str, Any] = {}
            for pair in pairs:
                if "=" not in pair:
                    raise _type_error(
                        name,
                        "JSON object or comma-separated key=value pairs",
                        value,
                    )
                key, raw_val = pair.split("=", 1)
                parsed_pairs[key.strip()] = raw_val.strip()
            raw_dict = parsed_pairs
    else:
        raise _type_error(name, "dict, JSON object string, or key=value list", value)

    if not isinstance(raw_dict, dict):
        raise _type_error(name, "dictionary-like object", value)

    coerced: dict[str, float] = {}
    for key, raw_val in raw_dict.items():
        key_name = str(key)
        parsed = _as_float(f"{name}.{key_name}", raw_val)
        if parsed is None:
            raise _type_error(f"{name}.{key_name}", "float", raw_val)
        coerced[key_name] = parsed
    return coerced


def _as_choice(*choices: str) -> Callable[[str, Any], str]:
    allowed = set(choices)

    def _cast(name: str, value: Any) -> str:
        parsed = _as_required_str(name, value)
        if parsed not in allowed:
            raise ValueError(
                f"Invalid value for '{name}': expected one of {sorted(allowed)}, got {parsed!r}."
            )
        return parsed

    return _cast


def _as_optional_choice(*choices: str) -> Callable[[str, Any], Optional[str]]:
    allowed = set(choices)

    def _cast(name: str, value: Any) -> Optional[str]:
        parsed = _as_optional_str(name, value)
        if parsed is None:
            return None
        if parsed not in allowed:
            raise ValueError(
                f"Invalid value for '{name}': expected one of {sorted(allowed)}, got {parsed!r}."
            )
        return parsed

    return _cast


def _as_stage_name_list(name: str, value: Any) -> list[str]:
    """Parse a comma-separated stage list into deterministic order-preserving tokens."""
    if value is None:
        return []
    if isinstance(value, str):
        raw_tokens = [tok.strip() for tok in value.split(",")]
    elif isinstance(value, (list, tuple, set)):
        raw_tokens = [str(tok).strip() for tok in value]
    else:
        raise _type_error(name, "comma-separated string or list of stage names", value)

    parsed: list[str] = []
    for token in raw_tokens:
        if not token:
            continue
        parsed.append(token)
    return parsed


def _invert_bool(name: str, value: Any) -> bool:
    parsed = _as_bool(name, value)
    if parsed is None:
        raise _type_error(name, "bool", value)
    return not parsed


def _as_timeout_or_none(name: str, value: Any) -> Optional[int]:
    parsed = _as_int(name, value, min_value=0)
    if parsed == 0:
        return None
    return parsed


_CONTEXT_ARG_SPECS = (
    ContextArgSpec("resume", "resume", _as_bool, True, numeric_or_bool=True),
    ContextArgSpec("force", "force", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "qc_policy",
        "qc_policy",
        _as_choice("off", "warn", "error"),
        "off",
    ),
    ContextArgSpec(
        "qc_density_rel_tol",
        "qc_density_rel_tol",
        _as_float,
        0.05,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "qc_enable_after_stages",
        "qc_enable_after_stages",
        _as_stage_name_list,
        ("gmx_eq",),
    ),
    ContextArgSpec(
        "all_qc_stop_on_warn",
        "all_qc_stop_on_warn",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_charge_model_mismatch",
        "allow_charge_model_mismatch",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_model_mismatch_reason",
        "charge_model_mismatch_reason",
        _as_optional_str,
        None,
    ),
    ContextArgSpec("verbose", "verbose", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("dry_run", "dry_run", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("strict_context_args", "strict_context_args", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("allow_override", "allow_override", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("atomtype_prefix", "atomtype_prefix", _as_optional_str, None),
    ContextArgSpec("strict_charge", "strict_charge", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "allow_unsanitized_grompp",
        "allow_unsanitized_grompp",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("charge_fix_method", "charge_fix_method", _as_choice("safe_subset", "uniform_all"), "safe_subset"),
    ContextArgSpec(
        "charge_fix_target_allowlist",
        "charge_fix_target_allowlist",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "polymer_net_charge_tol",
        "polymer_net_charge_tol",
        _as_float,
        1e-3,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_polymer_method",
        "charge_fix_polymer_method",
        _as_choice("skip_if_small", "spread_safe"),
        "skip_if_small",
    ),
    ContextArgSpec(
        "charge_fix_polymer_exclusion_bonds",
        "charge_fix_polymer_exclusion_bonds",
        lambda n, v: _as_int(n, v, min_value=0),
        2,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_max_delta_per_atom",
        "charge_fix_max_delta_per_atom",
        lambda n, v: _as_float(n, v, min_value=0.0),
        1e-5,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_max_total",
        "charge_fix_max_total",
        lambda n, v: _as_float(n, v, min_value=0.0),
        1e-4,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_max_dipole_shift_debye",
        "charge_fix_max_dipole_shift",
        lambda n, v: _as_float(n, v, min_value=0.0),
        0.01,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_moleculetype_rounding_tol",
        "charge_fix_moleculetype_rounding_tol",
        lambda n, v: _as_float(n, v, min_value=0.0),
        1e-4,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_allow_non_rounding",
        "charge_fix_allow_non_rounding",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("strict_charge_physics", "strict_charge_physics", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("strict_gro_top_check", "strict_gro_top_check", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("allow_default_defaults", "allow_default_defaults", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("allow_mixed_defaults", "allow_mixed_defaults", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "allow_mixed_defaults_reason",
        "allow_mixed_defaults_reason",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "mixed_defaults_cross_group_policy",
        "mixed_defaults_cross_group_policy",
        _as_optional_choice("off", "warn", "generate"),
        None,
    ),
    ContextArgSpec(
        "mixed_defaults_cross_group_rule",
        "mixed_defaults_cross_group_rule",
        _as_optional_choice("lorentz-berthelot", "geometric"),
        None,
    ),
    ContextArgSpec(
        "mixed_defaults_cross_group_max_pairs",
        "mixed_defaults_cross_group_max_pairs",
        lambda n, v: _as_int(n, v, min_value=1),
        20000,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "mixed_defaults_cross_group_reason",
        "mixed_defaults_cross_group_reason",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "mixed_defaults_preserve_within_group",
        "mixed_defaults_preserve_within_group",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "prefix_implicit_topology_policy",
        "prefix_implicit_topology_policy",
        _as_choice("error", "warn", "bake", "inject"),
        "error",
    ),
    ContextArgSpec(
        "prefix_implicit_topology_reason",
        "prefix_implicit_topology_reason",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "strict_forcefield_consistency",
        "strict_forcefield_consistency",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("temperature", "temperature", _as_float, None, numeric_or_bool=True),
    ContextArgSpec("pressure", "pressure", _as_float, None, numeric_or_bool=True),
    ContextArgSpec(
        "tc_grps_mode",
        "tc_grps_mode",
        _as_choice("auto", "system", "split"),
        "auto",
    ),
    ContextArgSpec("nsteps_em", "nsteps_em", _as_int, None, numeric_or_bool=True),
    ContextArgSpec("nsteps_nvt", "nsteps_nvt", _as_int, None, numeric_or_bool=True),
    ContextArgSpec("nsteps_npt", "nsteps_npt", _as_int, None, numeric_or_bool=True),
    ContextArgSpec("nsteps_md", "nsteps_md", _as_int, None, numeric_or_bool=True),
    ContextArgSpec(
        "gmx_eq_metrics_window_ps",
        "gmx_eq_metrics_window_ps",
        _as_float,
        20.0,
        numeric_or_bool=True,
    ),
    ContextArgSpec("tau_t_polymer", "tau_t_polymer", _as_float, None, numeric_or_bool=True),
    ContextArgSpec("tau_t_nonpolymer", "tau_t_nonpolymer", _as_float, None, numeric_or_bool=True),
    ContextArgSpec("allow_tau_t_autofill", "allow_tau_t_autofill", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "comm_grps_policy",
        "comm_grps_policy",
        _as_choice("auto", "match_tc", "system", "none", "explicit"),
        "match_tc",
    ),
    ContextArgSpec(
        "comm_mode",
        "comm_mode",
        _as_choice("Linear", "Angular", "None"),
        "Linear",
    ),
    ContextArgSpec("comm_grps_explicit", "comm_grps_explicit", _as_optional_str, None),
    ContextArgSpec(
        "allow_comm_mode_angular_with_pbc",
        "allow_comm_mode_angular_with_pbc",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("force_gen_vel", "force_gen_vel", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("gen_seed", "gen_seed", _as_int, None, numeric_or_bool=True),
    ContextArgSpec(
        "npt_barostat",
        "npt_barostat",
        _as_choice("Berendsen", "Parrinello-Rahman"),
        "Berendsen",
    ),
    ContextArgSpec("gmx_nt", "gmx_nt", _as_int, None, numeric_or_bool=True),
    ContextArgSpec("gmx_gpu_id", "gmx_gpu_id", _as_optional_str, None),
    ContextArgSpec("omp_num_threads", "omp_num_threads", _as_int, None, numeric_or_bool=True),
    ContextArgSpec(
        "htpolynet_inject_add_missing",
        "htpolynet_inject_add_missing",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("cleanup_success", "cleanup_success", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("cleanup_delete_trr", "cleanup_delete_trr", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "analysis_requires_velocities",
        "analysis_requires_velocities",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("grompp_maxwarn", "grompp_maxwarn", _as_int, 0, numeric_or_bool=True),
    ContextArgSpec("allow_velocity_reset", "allow_velocity_reset", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("allow_resume_repair", "allow_resume_repair", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("strict_restart_policy", "strict_restart_policy", _as_bool, None, numeric_or_bool=True),
    ContextArgSpec("resume_ignore_runtime", "resume_ignore_runtime", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("allow_placeholder", "allow_placeholder", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "allow_placeholder_propagate",
        "allow_placeholder_propagate",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_placeholder_stage_to_gromacs",
        "allow_placeholder_stage_to_gromacs",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_placeholder_gromacs_compile",
        "allow_placeholder_gromacs_compile",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_triclinic_unsafe_diagonal_approx",
        "allow_triclinic_unsafe_diagonal_approx",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("skip_gelation_precheck", "skip_gelation_precheck", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "allow_gelation_precheck_fail",
        "allow_gelation_precheck_fail",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_missing_htpolynet_qc",
        "allow_missing_htpolynet_qc",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("min_conversion", "min_conversion", _as_float, None, numeric_or_bool=True),
    ContextArgSpec(
        "htpolynet_timeout_sec",
        "htpolynet_timeout_sec",
        _as_timeout_or_none,
        43200,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "unsafe_allow_out_fallback",
        "unsafe_allow_out_fallback",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("strict_publish", "allow_partial_publish", _invert_bool, True, numeric_or_bool=True),
    ContextArgSpec(
        "allow_density_reduction",
        "allow_density_reduction",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "density_floor_fraction",
        "density_floor_fraction",
        _as_float,
        0.85,
        numeric_or_bool=True,
    ),
    ContextArgSpec("allow_default_recipe", "allow_default_recipe", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec("min_polymer_chains", "min_polymer_chains", _as_int, 100, numeric_or_bool=True),
    ContextArgSpec(
        "min_total_molecules_polymer",
        "min_total_molecules_polymer",
        _as_int,
        5000,
        numeric_or_bool=True,
    ),
    ContextArgSpec("strict_polymer_check", "strict_polymer_check", _as_bool, True, numeric_or_bool=True),
    ContextArgSpec("packmol_pdb_scale", "packmol_pdb_scale", _as_float, None, numeric_or_bool=True),
    ContextArgSpec(
        "packmol_edge_margin_nm",
        "packmol_edge_margin",
        _as_float,
        0.2,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "strict_packmol_units",
        "strict_packmol_units",
        _as_bool,
        True,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_max_box_expand_fraction",
        "packmol_max_box_expand_fraction",
        _as_float,
        0.03,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_min_contact_nm",
        "packmol_min_contact_nm",
        _as_float,
        0.18,
        numeric_or_bool=True,
    ),
    ContextArgSpec("packmol_random_seed", "packmol_random_seed", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "packmol_allow_density_floor_violation",
        "packmol_allow_density_floor_violation",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_preassembly_mode",
        "packmol_preassembly_mode",
        _as_choice("none", "li_near_polymer", "li_solvent_shell"),
        "none",
    ),
    ContextArgSpec(
        "packmol_preassembly_li_fraction",
        "packmol_preassembly_li_fraction",
        _as_float,
        1.0,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_preassembly_retry_count",
        "packmol_preassembly_retry_count",
        lambda n, v: _as_int(n, v, min_value=0),
        2,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_preassembly_sample_cap",
        "packmol_preassembly_sample_cap",
        lambda n, v: _as_int(n, v, min_value=1),
        4000,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_li_inner_exclusion_nm",
        "packmol_li_inner_exclusion_nm",
        _as_float,
        0.19,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_li_polymer_cutoff_nm",
        "packmol_li_polymer_cutoff_nm",
        _as_float,
        0.35,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_li_solvent_cutoff_nm",
        "packmol_li_solvent_cutoff_nm",
        _as_float,
        0.35,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_li_tfsi_cutoff_nm",
        "packmol_li_tfsi_cutoff_nm",
        _as_float,
        0.40,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_target_li_polymer_fraction",
        "packmol_target_li_polymer_fraction",
        _as_float,
        0.40,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_target_li_solvent_fraction",
        "packmol_target_li_solvent_fraction",
        _as_float,
        0.50,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "packmol_max_li_tfsi_close_fraction",
        "packmol_max_li_tfsi_close_fraction",
        _as_float,
        0.30,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_composition_changing_retries",
        "allow_composition_changing_retries",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("strict_mdp_validation", "strict_mdp_validation", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "allow_gen_vel_with_ref_t_spread",
        "allow_gen_vel_with_ref_t_spread",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "strict_charge_neutrality",
        "strict_charge_neutrality",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_neutrality_tol",
        "charge_neutrality_tol",
        _as_float,
        1e-6,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_neutrality_warn",
        "charge_neutrality_warn",
        _as_float,
        1e-4,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_neutrality_correct",
        "charge_neutrality_correct",
        _as_float,
        1e-4,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_fix_target_molecules",
        "charge_fix_target_molecules",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "charge_fix_target_atomnames",
        "charge_fix_target_atomnames",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "charge_fix_disallowed_atomtypes",
        "charge_fix_disallowed_atomtypes",
        _as_optional_str,
        None,
    ),
    ContextArgSpec("charge_fix_allow_ions", "charge_fix_allow_ions", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "charge_fix_source_counts",
        "charge_fix_source_counts",
        _as_choice("auto", "manifest", "system_top"),
        "auto",
    ),
    ContextArgSpec(
        "grompp_preprocess_timeout_s",
        "grompp_preprocess_timeout_s",
        lambda n, v: _as_int(n, v, min_value=1),
        None,
        numeric_or_bool=True,
    ),
    ContextArgSpec("itp_include_dirs", "itp_include_dirs", _as_optional_str, None),
    ContextArgSpec("tc_grps_groups", "tc_grps_groups", _as_optional_str, None),
    ContextArgSpec("tc_grps_ndx_path", "tc_grps_ndx_path", _as_optional_path, None),
    ContextArgSpec(
        "tc_grps_tau_t_values",
        "tc_grps_tau_t_values",
        _as_optional_float_list,
        None,
    ),
    ContextArgSpec(
        "tc_grps_min_atoms",
        "tc_grps_min_atoms",
        lambda n, v: _as_int(n, v, min_value=1),
        50,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "tc_grps_min_fraction",
        "tc_grps_min_fraction",
        _as_float,
        0.01,
        numeric_or_bool=True,
    ),
    ContextArgSpec("allow_small_tc_grps", "allow_small_tc_grps", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "system_type",
        "system_type",
        _as_optional_choice("gel", "solid", "liquid"),
        None,
    ),
    ContextArgSpec(
        "allow_continuation_override",
        "allow_continuation_override",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_gro_velocity_restart",
        "allow_gro_velocity_restart",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("posres_reference_gro", "posres_reference_gro", _as_optional_path, None),
    ContextArgSpec(
        "strict_posres_reference",
        "strict_posres_reference",
        _as_bool,
        True,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "strict_dipole_check",
        "strict_dipole_check",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "dipole_unwrap_max_atoms",
        "dipole_unwrap_max_atoms",
        lambda n, v: _as_int(n, v, min_value=1),
        512,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "dipole_group_max_atoms",
        "dipole_group_max_atoms",
        lambda n, v: _as_int(n, v, min_value=1),
        2048,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "dipole_percolation_threshold",
        "dipole_percolation_threshold",
        _as_float,
        0.45,
        numeric_or_bool=True,
    ),
    ContextArgSpec("ff_dir", "ff_dir", _as_optional_path, None),
    ContextArgSpec(
        "dipole_triclinic_policy",
        "dipole_triclinic_policy",
        _as_choice("skip", "error", "mic"),
        "skip",
    ),
    ContextArgSpec(
        "itp_include_dirs_priority",
        "itp_include_dirs_priority",
        _as_choice("sanitized_first", "forcefield_first", "first", "last"),
        "forcefield_first",
    ),
    ContextArgSpec(
        "charge_fix_protect_resnames",
        "charge_fix_protect_resnames",
        _as_optional_str,
        None,
    ),
    ContextArgSpec(
        "charge_fix_max_abs_delta_per_atom",
        "charge_fix_max_abs_delta_per_atom",
        lambda n, v: _as_float(n, v, min_value=0.0),
        1e-5,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "charge_protect_moleculetype_patterns",
        "charge_protect_moleculetype_patterns",
        _as_optional_str,
        None,
    ),
    ContextArgSpec("charge_protect_atomtypes", "charge_protect_atomtypes", _as_optional_str, None),
    ContextArgSpec(
        "strict_include_resolution",
        "strict_include_resolution",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec("strict_top_update", "strict_top_update", _as_bool, False, numeric_or_bool=True),
    ContextArgSpec(
        "strict_lj_validation",
        "strict_lj_validation",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "lj_outlier_policy",
        "lj_outlier_policy",
        _as_choice("off", "warn", "error"),
        "warn",
    ),
    ContextArgSpec(
        "lj_bounds_profile",
        "lj_bounds_profile",
        _as_choice("aa", "coarse_grained", "polarizable", "custom"),
        "aa",
    ),
    ContextArgSpec(
        "lj_sigma_max_nm",
        "lj_sigma_max_nm",
        lambda n, v: _as_float(n, v, min_value=0.0),
        None,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "lj_epsilon_max_kj_mol",
        "lj_epsilon_max_kj_mol",
        lambda n, v: _as_float(n, v, min_value=0.0),
        None,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "lj_outlier_thresholds",
        "lj_outlier_thresholds",
        _as_optional_threshold_dict,
        None,
    ),
    ContextArgSpec(
        "charge_fix_allow_solvents",
        "charge_fix_allow_solvents",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_include_shadowing",
        "allow_include_shadowing",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
    ContextArgSpec(
        "allow_unsafe_include_escape",
        "allow_unsafe_include_escape",
        _as_bool,
        False,
        numeric_or_bool=True,
    ),
)

_CONTEXT_SPEC_ARG_NAMES = {spec.arg_name for spec in _CONTEXT_ARG_SPECS}
_NUMERIC_OR_BOOL_FIELDS = {
    spec.ctx_field for spec in _CONTEXT_ARG_SPECS if spec.numeric_or_bool
}

_BASE_CONSUMED_ARG_NAMES = {
    "ff",
    "charge",
    "system",
    "run_id",
    "stage",
    "ff_key",
    "project_root",
    "allow_auto_run_id",
    "tau_t_values",
    "ff_charge_pair_supported",
    "ff_charge_doc_section",
}

_PIPELINE_ARG_PREFIXES = (
    "all_",
    "allow_",
    "analysis_",
    "charge_",
    "charge_fix_",
    "cleanup_",
    "comm_",
    "density_",
    "dipole_",
    "ff_",
    "gmx_",
    "grompp_",
    "htpolynet_",
    "itp_",
    "lj_",
    "min_",
    "npt_",
    "nsteps_",
    "omp_",
    "packmol_",
    "polymer_",
    "posres_",
    "qc_",
    "resume_",
    "strict_",
    "tau_t_",
    "tc_grps_",
    "unsafe_",
)

_PIPELINE_ARG_NAMES = {
    "temperature",
    "pressure",
    "gen_seed",
    "npt_barostat",
    "system_type",
}


def _is_pipeline_related_key(key: str) -> bool:
    if key in _CONTEXT_SPEC_ARG_NAMES or key in _PIPELINE_ARG_NAMES:
        return True
    return key.startswith(_PIPELINE_ARG_PREFIXES)


@dataclass
class PipelineContext:
    """
    Runtime context passed to all pipeline stages.
    
    Contains configuration options and shared resources.
    """
    
    # Configuration options
    ff: str  # Forcefield: GAFF2 or OPLS-AA (display name)
    charge: str  # Charge method: RESP or CM5
    system_id: str  # System identifier
    run_id: str  # Run identifier
    stage: str  # Target stage or 'all'
    ff_key: str = ""  # Normalized forcefield key for path matching (e.g., OPLSAA)
    resume: bool = True  # Skip completed stages
    force: bool = False  # Re-run even if complete
    qc_policy: str = "off"  # Dispatcher quality-gate policy: off | warn | error
    qc_density_rel_tol: float = 0.05  # Relative density tolerance for QC vs expected
    qc_enable_after_stages: list[str] = field(default_factory=lambda: ["gmx_eq"])
    all_qc_stop_on_warn: bool = False  # Stop early in all-mode after QC warning
    allow_charge_model_mismatch: bool = False
    charge_model_mismatch_reason: Optional[str] = None
    ff_charge_pair_supported: bool = True
    
    # Shared resources (initialized post-creation)
    path_manager: Optional[PathManager] = field(default=None, repr=False)
    manifest: Optional[ManifestWriter] = field(default=None, repr=False)
    
    # Runtime state
    project_root: Optional[Path] = None
    verbose: bool = False
    dry_run: bool = False
    strict_context_args: bool = False  # Fail on unconsumed pipeline-related args in from_args()
    run_interrupted: bool = False  # Set by dispatcher on KeyboardInterrupt
    
    # ITP Sanitizer options
    allow_override: bool = False
    atomtype_prefix: Optional[str] = None
    strict_charge: bool = False  # Fail on atomtype charge mismatches
    allow_unsanitized_grompp: bool = False  # Allow grompp without sanitizer outputs (not recommended)
    
    # GROMACS MDP patching options
    temperature: Optional[float] = None  # Override ref_t and gen_temp
    pressure: Optional[float] = None  # Override ref_p
    tc_grps_mode: str = "auto"  # 'auto', 'system', or 'split'
    nsteps_em: Optional[int] = None
    nsteps_nvt: Optional[int] = None
    nsteps_npt: Optional[int] = None
    nsteps_md: Optional[int] = None
    gmx_eq_metrics_window_ps: float = 20.0  # Last-window averaging window for gmx_eq metrics
    
    # tc-grps split safety options (Issue A)
    tau_t_polymer: Optional[float] = None  # tau_t for Polymer group (ps)
    tau_t_nonpolymer: Optional[float] = None  # tau_t for NonPolymer group (ps)
    allow_tau_t_autofill: bool = False  # Allow duplicate tau_t when split (legacy)
    
    # COM removal policy for split tc-grps (Issue B)
    comm_grps_policy: str = "match_tc"  # "auto" | "match_tc" | "system" | "none" | "explicit"
    comm_mode: str = "Linear"  # "Linear" | "Angular" | "None"
    comm_grps_explicit: Optional[str] = None  # Only used if policy="explicit"
    allow_comm_mode_angular_with_pbc: bool = False  # Expert escape hatch (unsafe under PBC)
    
    # Stage transition velocity handling (Issue C)
    force_gen_vel: bool = False  # Force gen_vel=yes even in production
    
    # Reproducibility: velocity generation seed
    gen_seed: Optional[int] = None  # CLI-provided seed for gen_seed (overrides template default)
    
    # NPT barostat selection (equilibration vs production-like)
    npt_barostat: str = "Berendsen"  # "Berendsen" (equil stable) or "Parrinello-Rahman" (correct fluctuations)
    
    # Resource scheduling / pinning
    gmx_nt: Optional[int] = None  # mdrun -nt
    gmx_gpu_id: Optional[str] = None  # mdrun -gpu_id
    omp_num_threads: Optional[int] = None  # OMP_NUM_THREADS
    htpolynet_inject_add_missing: bool = False  # Allow creating gromacs.ntomp/gpu_id in HTPolyNet config injection
    htpolynet_injection_method: Optional[str] = None  # Runtime: config | env | none
    
    # Cleanup options
    cleanup_success: bool = False
    cleanup_delete_trr: bool = False  # Explicit opt-in to delete TRR (velocities needed for Green-Kubo)
    analysis_requires_velocities: bool = False  # Signal that analysis needs velocities (e.g., Green-Kubo)
    
    # Strict reproducibility options (Part 1)
    grompp_maxwarn: int = 0  # Default strict: no warnings allowed
    allow_velocity_reset: bool = False  # Require checkpoint for continuation
    allow_resume_repair: bool = False  # Allow rebuilding missing md.tpr from md.cpt
    strict_restart_policy: Optional[bool] = None  # If None, fallback to strict_mdp_validation
    resume_ignore_runtime: bool = False  # Ignore runtime-only fingerprint keys during resume compare
    
    # HTPolyNet options
    allow_placeholder: bool = False  # Allow placeholder outputs when htpolynet not found
    allow_placeholder_propagate: bool = False  # Allow placeholder outputs to be published to IN/htpolynet
    allow_placeholder_stage_to_gromacs: bool = False  # Explicit opt-in to stage propagated placeholders to IN/gromacs
    allow_placeholder_gromacs_compile: bool = False  # Disable placeholder poison pill (unsafe)
    allow_triclinic_unsafe_diagonal_approx: bool = False  # Unsafe triclinic -> diagonal collapse policy
    skip_gelation_precheck: bool = False  # Skip conservative gelation feasibility precheck
    allow_gelation_precheck_fail: bool = False  # Allow stage to continue after failed precheck
    allow_missing_htpolynet_qc: bool = False  # Allow missing conversion/gel-fraction parse
    min_conversion: Optional[float] = None  # Minimum crosslinking conversion threshold (0.0-1.0)
    htpolynet_timeout_sec: Optional[int] = 43200  # Subprocess timeout (default: 12h; None=no limit)
    unsafe_allow_out_fallback: bool = False  # Allow OUT/ fallback for initial structure (breaks provenance)
    
    # PACKMOL strict reproducibility options
    strict_publish: bool = True  # Fail stage on publish errors (strict mode)
    allow_density_reduction: bool = False  # Allow retries to reduce density below floor
    density_floor_fraction: float = 0.85  # Minimum density as fraction of target
    allow_composition_changing_retries: bool = False  # A5: Allow molecule count changes on retry (DANGEROUS)
    
    # PACKMOL polymer percolation guardrails
    allow_default_recipe: bool = False  # Allow demo recipe when composition.yaml missing
    min_polymer_chains: int = 100  # Minimum polymer chains for percolation
    min_total_molecules_polymer: int = 5000  # Minimum total for polymer systems
    strict_polymer_check: bool = True  # Fail on polymer threshold violations
    
    # PACKMOL unit handling
    packmol_pdb_scale: Optional[float] = None  # Override: 0.1 for Å, 1.0 for nm
    packmol_edge_margin_nm: float = 0.2  # Edge margin for molecule placement
    strict_packmol_units: bool = True  # Fail-fast on ambiguous unit inference
    packmol_max_box_expand_fraction: float = 0.03  # Hard cap for auto box expansion during PDB->GRO
    packmol_min_contact_nm: float = 0.18  # Li-heavy minimum safe contact distance
    packmol_random_seed: bool = False  # Explicit opt-in for non-deterministic PACKMOL seed
    packmol_allow_density_floor_violation: bool = False  # Unsafe override for achieved density floor

    # PACKMOL optional preassembly hints
    packmol_preassembly_mode: str = "none"  # none | li_near_polymer | li_solvent_shell
    packmol_preassembly_li_fraction: float = 1.0
    packmol_preassembly_retry_count: int = 2
    packmol_preassembly_sample_cap: int = 4000
    packmol_li_inner_exclusion_nm: float = 0.19
    packmol_li_polymer_cutoff_nm: float = 0.35
    packmol_li_solvent_cutoff_nm: float = 0.35
    packmol_li_tfsi_cutoff_nm: float = 0.40
    packmol_target_li_polymer_fraction: float = 0.40
    packmol_target_li_solvent_fraction: float = 0.50
    packmol_max_li_tfsi_close_fraction: float = 0.30
    
    # Sanitizer strict reproducibility / physics validation options
    charge_fix_method: str = "safe_subset"  # "safe_subset" or "uniform_all"
    charge_fix_target_allowlist: Optional[str] = None  # Comma-separated allowed targets
    polymer_net_charge_tol: float = 1e-3  # Acceptable per-molecule polymer net-charge drift (e)
    charge_fix_polymer_method: str = "skip_if_small"  # skip_if_small | spread_safe
    charge_fix_polymer_exclusion_bonds: int = 2  # Exclude atoms within N bonds of hetero atoms
    charge_fix_max_delta_per_atom: Optional[float] = 1e-5  # Max |Δq| per atom for charge correction
    charge_fix_max_total: Optional[float] = 1e-4  # Max total |ΔQ| (relaxed for FF rounding)
    charge_fix_max_dipole_shift_debye: Optional[float] = 0.01  # Max |Δμ| in Debye
    charge_fix_moleculetype_rounding_tol: float = 1e-4  # Max |Q_mol-round(Q_mol)| treated as rounding drift
    charge_fix_allow_non_rounding: bool = False  # Explicit unsafe override for non-rounding drift fixes
    strict_charge_physics: bool = False  # Fail on dipole shift violations
    strict_gro_top_check: bool = False  # Fail on uncertain atom counts (preprocessor ambiguity)
    allow_default_defaults: bool = False  # Allow fallback [ defaults ] if none found in ITPs
    allow_mixed_defaults: bool = False  # Explicit unsafe override for mixed [defaults] tuples
    allow_mixed_defaults_reason: Optional[str] = None  # Required audit reason for mixed defaults override
    mixed_defaults_cross_group_policy: Optional[str] = None  # None=>auto (off unless allow_mixed_defaults, then warn)
    mixed_defaults_cross_group_rule: Optional[str] = None  # None=>auto from canonical comb-rule
    mixed_defaults_cross_group_max_pairs: int = 20000  # Pair-generation cap for cross-group nonbond_params
    mixed_defaults_cross_group_reason: Optional[str] = None  # Required when policy=generate
    mixed_defaults_preserve_within_group: bool = False  # Guard rail: preserve secondary within-group LJ via nonbond_params
    prefix_implicit_topology_policy: str = "error"  # error | warn | bake | inject
    prefix_implicit_topology_reason: Optional[str] = None  # Required when policy in {bake,inject}
    strict_forcefield_consistency: bool = False  # Fail on forcefield mixing (e.g., GAFF2+OPLS)
    
    # Part 2: New charge neutrality hardening options
    strict_charge_neutrality: bool = False  # Fail-fast on preprocessor, ambiguous parsing
    charge_neutrality_tol: float = 1e-6  # Tolerance for |Q| < tol → neutral
    charge_neutrality_warn: float = 1e-4  # Threshold for warning
    charge_neutrality_correct: float = 1e-4  # Max |Q| for auto-correction
    charge_fix_target_molecules: Optional[str] = None  # Comma-separated target molecules
    charge_fix_target_atomnames: Optional[str] = None  # Comma-separated atom name patterns
    charge_fix_disallowed_atomtypes: Optional[str] = None  # Comma-separated disallowed types
    charge_fix_allow_ions: bool = False  # Allow correction on detected ions
    charge_fix_source_counts: str = "auto"  # "auto", "manifest", or "system_top"
    
    # Task B: Configurable grompp preprocessing timeout
    grompp_preprocess_timeout_s: Optional[int] = None  # Timeout for grompp -pp (env: GMX_GROMPP_PP_TIMEOUT_S)
    
    # Task C: Configurable ITP include directories for #include resolution
    itp_include_dirs: Optional[str] = None  # Comma-separated paths for #include resolution
    
    # MDP validation strictness (Task E: strict vs non-strict mode)
    strict_mdp_validation: bool = False  # Escalate MDP validation warnings to errors
    
    # MDP Patcher Hardening: allow gen_vel with ref_t spread (unsafe opt-in)
    allow_gen_vel_with_ref_t_spread: bool = False  # Expert override for thermal gradient risk
    
    # MDP Patcher Hardening: tc-grps split configurability (Issue 2)
    tc_grps_groups: Optional[str] = None  # Custom group names, e.g. "PEO LiTFSI"
    tc_grps_ndx_path: Optional[Path] = None  # Path to ndx file for group detection
    tc_grps_tau_t_values: Optional[list] = None  # Per-group tau_t values (split mode)
    tc_grps_min_atoms: int = 50  # Tiny-group DOF safety threshold for split tc-grps
    tc_grps_min_fraction: float = 0.01  # Tiny-group atom-fraction threshold for split tc-grps
    allow_small_tc_grps: bool = False  # Expert opt-in to allow tiny split groups
    
    # MDP Patcher Hardening: System type for barostat tau_p defaults (Issue 3)
    system_type: Optional[str] = None  # "gel", "solid", "liquid", or None
    
    # MDP Patcher Hardening: Continuation semantics override (Issue 1)
    allow_continuation_override: bool = False  # Allow gen_vel/continuation mismatch
    
    # MDP Patcher Hardening: Expert restart from .gro with velocities (P0.2)
    # Allows (gen_vel=no, continuation=no) combo for restart from .gro that contains velocities.
    # RISK: MDP patcher cannot verify if .gro contains velocities; user must ensure this.
    allow_gro_velocity_restart: bool = False
    input_gro_has_velocities: Optional[bool] = None  # Stage-resolved .gro velocity availability
    active_ndx_path: Optional[str] = None  # Stage-resolved ndx path for group validation
    
    # POSRES Reference: explicit -r for position restraints (new)
    posres_reference_gro: Optional[Path] = None  # Explicit -r reference structure for POSRES
    strict_posres_reference: bool = True  # Fail if POSRES active but no explicit reference
    
    # Sanitizer v4: Dipole check hardening (heuristic by default)
    strict_dipole_check: bool = False  # Only enforce dipole limits if True AND preconditions pass
    dipole_unwrap_max_atoms: int = 512  # Skip dipole enforcement for molecules with more atoms
    dipole_group_max_atoms: int = 2048  # Per-group atom cap before any unwrap work
    dipole_percolation_threshold: float = 0.45  # Span ratio threshold for percolation gating
    
    # Sanitizer v4: Forcefield override
    ff_dir: Optional[Path] = None  # Explicit forcefield directory path
    
    # Sanitizer v5: Triclinic box policy for dipole/PBC unwrap (Issue C)
    # "skip" = skip dipole check with reason; "error" = fail; "mic" = correct triclinic MIC
    dipole_triclinic_policy: str = "skip"
    
    # Sanitizer v5: Include path priority (Issue D)
    # Canonical values:
    # - forcefield_first (legacy alias: "last")
    # - sanitized_first (legacy alias: "first")
    itp_include_dirs_priority: str = "forcefield_first"
    
    # Sanitizer v5: Ion protection for charge-fix (Issue G)
    charge_fix_protect_resnames: Optional[str] = None  # Comma-separated resnames to protect
    charge_fix_max_abs_delta_per_atom: float = 1e-5  # Max delta per atom before fail
    
    # Sanitizer v4: Charge protection patterns
    charge_protect_moleculetype_patterns: Optional[str] = None  # Comma-separated regex patterns
    charge_protect_atomtypes: Optional[str] = None  # Comma-separated atomtype denylist
    
    # Sanitizer v4: Include resolution strictness
    strict_include_resolution: bool = False  # Fail on unresolved includes (default: warn)
    strict_top_update: bool = False  # Fail on malformed sanitizer managed markers in system.top
    
    # Hardening v6: LJ parameter validation strictness
    strict_lj_validation: bool = False  # Fail on LJ parameter outliers (default: warn)
    lj_outlier_policy: str = "warn"  # off | warn | error
    lj_bounds_profile: str = "aa"  # aa | coarse_grained | polarizable | custom
    lj_sigma_max_nm: Optional[float] = None  # Used for lj_bounds_profile=custom
    lj_epsilon_max_kj_mol: Optional[float] = None  # Used for lj_bounds_profile=custom
    lj_outlier_thresholds: Optional[dict] = None  # Custom thresholds {sigma_max_nm, epsilon_max_kj_mol, etc.}
    
    # Hardening v6: Solvent/polymer protection for charge correction
    charge_fix_allow_solvents: bool = False  # Allow correction on detected solvents
    
    # Hardening v6: Include shadowing detection
    allow_include_shadowing: bool = False  # Allow shadowing without warning
    allow_unsafe_include_escape: bool = False  # Allow includes to resolve outside configured roots
    
    def __post_init__(self):
        """Validate options after initialization."""
        valid_ffs = {"GAFF2", "OPLS-AA"}
        if self.ff not in valid_ffs:
            raise ValueError(f"Invalid forcefield: {self.ff}. Must be one of {valid_ffs}")
        
        valid_charges = {"RESP", "CM5"}
        if self.charge not in valid_charges:
            raise ValueError(f"Invalid charge method: {self.charge}. Must be one of {valid_charges}")

        pair = (self.ff, self.charge)
        pair_supported = pair in SUPPORTED_FF_CHARGE_PAIRS
        self.ff_charge_pair_supported = pair_supported
        if not pair_supported and not self.allow_charge_model_mismatch:
            supported = ", ".join(
                f"{ff}/{charge}" for ff, charge in sorted(SUPPORTED_FF_CHARGE_PAIRS)
            )
            raise ValueError(
                f"Unsupported forcefield/charge-model pair {self.ff}/{self.charge}. "
                f"Supported pair(s): {supported}. "
                f"See {FF_CHARGE_COMPAT_DOC_SECTION}."
            )
        if not pair_supported and self.allow_charge_model_mismatch:
            reason = (self.charge_model_mismatch_reason or "").strip()
            if not reason:
                raise ValueError(
                    "Unsafe forcefield/charge override requires charge_model_mismatch_reason."
                )
            self.charge_model_mismatch_reason = reason

        if self.allow_mixed_defaults:
            reason = (self.allow_mixed_defaults_reason or "").strip()
            if not reason:
                raise ValueError(
                    "Unsafe mixed-defaults override requires allow_mixed_defaults_reason."
                )
            self.allow_mixed_defaults_reason = reason
        else:
            self.allow_mixed_defaults_reason = None

        if self.mixed_defaults_cross_group_policy is None:
            self.mixed_defaults_cross_group_policy = (
                "warn" if self.allow_mixed_defaults else "off"
            )
        if self.mixed_defaults_cross_group_policy != "off" and not self.allow_mixed_defaults:
            raise ValueError(
                "mixed_defaults_cross_group_policy requires allow_mixed_defaults "
                "(or set policy='off')."
            )
        if self.mixed_defaults_cross_group_max_pairs <= 0:
            raise ValueError("mixed_defaults_cross_group_max_pairs must be > 0.")
        if self.mixed_defaults_cross_group_policy == "generate":
            reason = (self.mixed_defaults_cross_group_reason or "").strip()
            if not reason:
                raise ValueError(
                    "mixed_defaults_cross_group_policy='generate' requires "
                    "mixed_defaults_cross_group_reason."
                )
            self.mixed_defaults_cross_group_reason = reason
        else:
            self.mixed_defaults_cross_group_reason = None

        if self.mixed_defaults_preserve_within_group and not self.allow_mixed_defaults:
            raise ValueError(
                "mixed_defaults_preserve_within_group requires allow_mixed_defaults."
            )

        if self.prefix_implicit_topology_policy in {"bake", "inject"}:
            reason = (self.prefix_implicit_topology_reason or "").strip()
            if not reason:
                raise ValueError(
                    "prefix_implicit_topology_policy in {bake,inject} requires "
                    "prefix_implicit_topology_reason."
                )
            self.prefix_implicit_topology_reason = reason
        else:
            self.prefix_implicit_topology_reason = None

        # Backward compatibility: strict_lj_validation escalates to error unless policy is explicit.
        if self.strict_lj_validation and self.lj_outlier_policy == "warn":
            self.lj_outlier_policy = "error"

        if self.lj_bounds_profile == "custom":
            if self.lj_sigma_max_nm is None or self.lj_sigma_max_nm <= 0:
                raise ValueError(
                    "lj_bounds_profile='custom' requires lj_sigma_max_nm > 0."
                )
            if self.lj_epsilon_max_kj_mol is None or self.lj_epsilon_max_kj_mol <= 0:
                raise ValueError(
                    "lj_bounds_profile='custom' requires lj_epsilon_max_kj_mol > 0."
                )
        else:
            if self.lj_sigma_max_nm is not None and self.lj_sigma_max_nm <= 0:
                raise ValueError("lj_sigma_max_nm must be > 0 when provided.")
            if self.lj_epsilon_max_kj_mol is not None and self.lj_epsilon_max_kj_mol <= 0:
                raise ValueError("lj_epsilon_max_kj_mol must be > 0 when provided.")

        if self.qc_policy not in {"off", "warn", "error"}:
            raise ValueError(
                f"Invalid qc_policy: {self.qc_policy!r}. Must be one of off/warn/error."
            )
        if not (0.0 < float(self.qc_density_rel_tol) < 1.0):
            raise ValueError("qc_density_rel_tol must be in (0, 1).")
        if self.qc_policy == "error" and not (0.0 < float(self.qc_density_rel_tol) < 0.5):
            raise ValueError(
                "qc_policy='error' requires qc_density_rel_tol in (0, 0.5)."
            )
        if float(self.gmx_eq_metrics_window_ps) <= 0.0:
            raise ValueError("gmx_eq_metrics_window_ps must be > 0.")

    def _validate_no_string_in_numeric_fields(self) -> None:
        """
        Fail fast if numeric/bool context fields still contain raw strings.
        """
        offenders = []
        for field_name in sorted(_NUMERIC_OR_BOOL_FIELDS):
            value = getattr(self, field_name, None)
            if isinstance(value, str):
                offenders.append((field_name, value))
        if offenders:
            details = ", ".join(f"{name}={raw!r}" for name, raw in offenders)
            raise ValueError(
                "Context type validation failed: numeric/bool fields contain strings: "
                f"{details}. Check wrapper/orchestrator argument typing."
            )
    
    @classmethod
    def from_args(cls, args, project_root: Path) -> "PipelineContext":
        """
        Create context from parsed CLI arguments.
        
        Args:
            args: Parsed argparse namespace
            project_root: Project root directory
            
        Returns:
            Initialized PipelineContext
        """
        args_dict = dict(args) if isinstance(args, dict) else vars(args).copy()
        ff = _as_required_str("ff", args_dict.get("ff"))
        charge = _as_required_str("charge", args_dict.get("charge"))
        system_id = _as_required_str("system", args_dict.get("system"))
        run_id = _as_required_str("run_id", args_dict.get("run_id"))
        stage = _as_required_str("stage", args_dict.get("stage"))
        ff_key = _as_required_str(
            "ff_key",
            args_dict.get("ff_key", ff.replace("-", "")),
        )

        project_root_path = Path(project_root)
        path_manager = PathManager(project_root_path)

        # Ensure output directory exists
        run_dir = path_manager.ensure_output_dir(run_id)
        manifest_path = run_dir / "manifest.json"
        manifest = ManifestWriter(manifest_path)

        kwargs: dict[str, Any] = {
            "ff": ff,
            "ff_key": ff_key,
            "charge": charge,
            "system_id": system_id,
            "run_id": run_id,
            "stage": stage,
            "path_manager": path_manager,
            "manifest": manifest,
            "project_root": project_root_path,
        }

        for spec in _CONTEXT_ARG_SPECS:
            raw_value = args_dict.get(spec.arg_name, MISSING)
            if raw_value is MISSING:
                if spec.default_value is MISSING:
                    continue
                kwargs[spec.ctx_field] = spec.default_value
                continue
            if raw_value is None:
                kwargs[spec.ctx_field] = None
            else:
                kwargs[spec.ctx_field] = spec.caster(spec.arg_name, raw_value)

        ctx = cls(**kwargs)
        ctx._validate_no_string_in_numeric_fields()

        consumed_arg_names = set(_BASE_CONSUMED_ARG_NAMES)
        consumed_arg_names.update(_CONTEXT_SPEC_ARG_NAMES)
        unconsumed_pipeline_keys = [
            key
            for key in sorted(args_dict)
            if key not in consumed_arg_names and _is_pipeline_related_key(key)
        ]
        if unconsumed_pipeline_keys:
            formatted = ", ".join(
                f"{name}={args_dict.get(name)!r}" for name in unconsumed_pipeline_keys
            )
            message = (
                "Unconsumed pipeline-related context args detected: "
                f"{formatted}. Add mappings to _CONTEXT_ARG_SPECS."
            )
            if bool(ctx.strict_context_args):
                raise ValueError(message)
            if bool(ctx.verbose):
                print(f"[WARN] {message}", file=sys.stderr)

        # Record options in manifest from normalized context
        manifest.set_options(
            ff=ctx.ff,
            charge=ctx.charge,
            system_id=ctx.system_id,
            run_id=ctx.run_id,
            stage=ctx.stage,
            resume=ctx.resume,
            force=ctx.force,
        )

        supported_pairs = [
            {"ff": ff_name, "charge_model": charge_name}
            for ff_name, charge_name in sorted(SUPPORTED_FF_CHARGE_PAIRS)
        ]
        pair_supported = bool(ctx.ff_charge_pair_supported)
        mismatch_allowed = bool(ctx.allow_charge_model_mismatch)
        mismatch_reason = ctx.charge_model_mismatch_reason
        manifest.set(
            "forcefield_charge_compatibility",
            {
                "ff": ctx.ff,
                "charge_model": ctx.charge,
                "pair_supported": pair_supported,
                "supported_pairs": supported_pairs,
                "doc_section": FF_CHARGE_COMPAT_DOC_SECTION,
                "unsafe_override": bool((not pair_supported) and mismatch_allowed),
                "unsafe_override_reason": mismatch_reason if (not pair_supported) else None,
                "validated_parameter_set": pair_supported,
            },
        )

        if not pair_supported and mismatch_allowed:
            unsafe = manifest.get("unsafe_overrides", [])
            if not isinstance(unsafe, list):
                unsafe = []
            unsafe.append(
                {
                    "kind": "forcefield_charge_model_mismatch",
                    "unsafe_override": True,
                    "reason": mismatch_reason,
                    "ff": ctx.ff,
                    "charge_model": ctx.charge,
                }
            )
            manifest.set("unsafe_overrides", unsafe)

        return ctx
    
    def get_input_path(self, *parts: str) -> Path:
        """Convenience wrapper for path_manager.get_input_path()."""
        if self.path_manager is None:
            raise RuntimeError("PathManager not initialized")
        return self.path_manager.get_input_path(*parts)
    
    def get_output_path(self, *parts: str) -> Path:
        """Convenience wrapper for path_manager.get_output_path()."""
        if self.path_manager is None:
            raise RuntimeError("PathManager not initialized")
        return self.path_manager.get_output_path(self.run_id, *parts)
    
    def ensure_output_dir(self, *parts: str) -> Path:
        """Convenience wrapper for path_manager.ensure_output_dir()."""
        if self.path_manager is None:
            raise RuntimeError("PathManager not initialized")
        return self.path_manager.ensure_output_dir(self.run_id, *parts)
    
    def log_command(self, command, stage: str, **kwargs) -> None:
        """Convenience wrapper for manifest.log_command()."""
        if self.manifest is not None:
            self.manifest.log_command(command, stage, **kwargs)
