"""
CLI entrypoint for the MD simulation pipeline.

Usage:
    python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage all
"""

import argparse
import re
import sys
import uuid
from datetime import datetime
from pathlib import Path

from .context import PipelineContext
from .dispatcher import StageDispatcher
from .stages import STAGE_ORDER


# Mapping from CLI-friendly display names to internal keys (no dashes for path matching)
FORCEFIELD_DISPLAY_TO_KEY = {
    "GAFF2": "GAFF2",
    "OPLS-AA": "OPLSAA",
}

# Forcefield/charge-model compatibility matrix (validated default-supported pairs).
# Any pair outside this set requires explicit unsafe override.
SUPPORTED_FF_CHARGE_PAIRS = {
    ("GAFF2", "CM5"),
    ("OPLS-AA", "RESP"),
}
FF_CHARGE_COMPAT_DOC_SECTION = "AGENT_SKILLS/03_ASSET_RESOLUTION.md#compatibility-matrix"


def normalize_forcefield(ff_display: str) -> str:
    """
    Convert display forcefield name to internal key.
    
    Args:
        ff_display: CLI-friendly name (e.g., "OPLS-AA")
        
    Returns:
        Internal key for path matching (e.g., "OPLSAA")
    """
    return FORCEFIELD_DISPLAY_TO_KEY.get(ff_display, ff_display.replace("-", ""))


def generate_run_id() -> str:
    """Generate a unique run ID based on timestamp and UUID."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    short_uuid = uuid.uuid4().hex[:6]
    return f"run_{timestamp}_{short_uuid}"


def _format_supported_pairs(pairs) -> str:
    """Format supported (ff, charge) pairs for user-facing diagnostics."""
    return ", ".join(f"{ff}/{charge}" for ff, charge in sorted(pairs))


def validate_args(parser: argparse.ArgumentParser, args) -> None:
    """
    Validate co-dependent and mutually-exclusive CLI arguments.
    
    Uses parser.error() for hard failures (exits with code 2).
    This is called after parse_args() and before any execution.
    
    Args:
        parser: ArgumentParser instance for error reporting
        args: Parsed arguments namespace
    """
    # Default for parsed tau_t values
    if not hasattr(args, "tc_grps_tau_t_values"):
        args.tc_grps_tau_t_values = None
    # =========================================================================
    # A1) tc-grps split requires both tau values (or autofill flag)
    # =========================================================================
    if args.tc_grps_mode == "split":
        # Parse custom groups (whitespace-delimited tokens)
        group_list = []
        if args.tc_grps_groups:
            group_list = [g for g in args.tc_grps_groups.split() if g]
            if len(group_list) < 2:
                parser.error(
                    "--tc-grps-groups must contain at least 2 group names for split mode.\n"
                    "Example: --tc-grps-groups 'Li TFSI Solvent'"
                )
        
        # Parse tau_t_values if provided (comma-separated)
        tau_t_values = None
        if args.tau_t_values:
            raw_vals = [v.strip() for v in args.tau_t_values.split(",") if v.strip()]
            if not raw_vals:
                parser.error("--tau-t-values provided but no values were parsed.")
            try:
                tau_t_values = [float(v) for v in raw_vals]
            except ValueError:
                parser.error(
                    f"--tau-t-values contains non-numeric entries: '{args.tau_t_values}'.\n"
                    "Example: --tau-t-values '0.1,0.5,2.0'"
                )
            if any(v <= 0 for v in tau_t_values):
                parser.error("--tau-t-values must all be > 0.")
        
        has_poly = args.tau_t_polymer is not None
        has_nonpoly = args.tau_t_nonpolymer is not None
        has_legacy = has_poly or has_nonpoly
        
        if tau_t_values is not None and has_legacy:
            parser.error(
                "--tau-t-values cannot be combined with --tau-t-polymer/--tau-t-nonpolymer.\n"
                "Use one method to avoid ambiguous group mapping."
            )
        
        if group_list:
            # Custom split groups
            if tau_t_values is not None:
                if len(tau_t_values) != len(group_list):
                    parser.error(
                        "--tau-t-values count must match --tc-grps-groups count.\n"
                        f"Groups ({len(group_list)}): {group_list}\n"
                        f"tau_t values ({len(tau_t_values)}): {tau_t_values}"
                    )
                args.tc_grps_tau_t_values = tau_t_values
            else:
                # Optional exception: legacy tau_t for Polymer/NonPolymer
                if (len(group_list) == 2 and group_list == ["Polymer", "NonPolymer"]
                        and has_poly and has_nonpoly):
                    print(
                        "[INFO] Using legacy --tau-t-polymer/--tau-t-nonpolymer with "
                        "--tc-grps-groups 'Polymer NonPolymer'.",
                        file=sys.stderr,
                    )
                elif args.allow_tau_t_autofill:
                    print(
                        "[WARN] tc-grps split with custom groups but no --tau-t-values provided. "
                        "Autofill may duplicate a single tau_t from the MDP template, which can be "
                        "physically inappropriate for heterogeneous groups.",
                        file=sys.stderr,
                    )
                else:
                    parser.error(
                        "Custom tc-grps split requires --tau-t-values with one value per group.\n"
                        "Example: --tc-grps-groups 'Li TFSI Solvent' --tau-t-values '0.1,0.5,2.0'\n"
                        "If you intentionally want duplication, add --allow-tau-t-autofill (not recommended)."
                    )
        else:
            # Legacy Polymer/NonPolymer split
            if tau_t_values is not None:
                parser.error(
                    "--tau-t-values requires --tc-grps-groups for custom split groups.\n"
                    "For legacy Polymer/NonPolymer split, use --tau-t-polymer/--tau-t-nonpolymer."
                )
            if not args.allow_tau_t_autofill:
                if has_poly != has_nonpoly:  # XOR - only one provided
                    parser.error(
                        "Legacy Polymer/NonPolymer split requires BOTH --tau-t-polymer and "
                        "--tau-t-nonpolymer, or use --allow-tau-t-autofill (not recommended).\n"
                        "Example: --tc-grps-mode split --tau-t-polymer 0.5 --tau-t-nonpolymer 0.1\n"
                        "For custom splits, use --tc-grps-groups + --tau-t-values."
                    )
                if not has_poly and not has_nonpoly:
                    parser.error(
                        "Legacy Polymer/NonPolymer split requires --tau-t-polymer and "
                        "--tau-t-nonpolymer.\n"
                        "Example: --tc-grps-mode split --tau-t-polymer 0.5 --tau-t-nonpolymer 0.1\n"
                        "For custom splits, use --tc-grps-groups + --tau-t-values."
                    )
        
        # Warn about implicit match_tc expectation in legacy mode
        if args.tc_grps_mode == "split" and not args.tc_grps_groups:
            print(
                "[WARN] tc-grps split without --tc-grps-groups uses legacy "
                "'Polymer NonPolymer' group names. Ensure your .ndx contains these groups.",
                file=sys.stderr,
            )
            # Optional early ndx inspection
            if args.tc_grps_ndx_path:
                ndx_path = Path(args.tc_grps_ndx_path)
                if ndx_path.exists():
                    try:
                        group_names = []
                        for line in ndx_path.read_text().splitlines():
                            m = re.match(r'^\s*\[\s*(.+?)\s*\]\s*$', line)
                            if m:
                                group_names.append(m.group(1))
                        if "Polymer" not in group_names or "NonPolymer" not in group_names:
                            print(
                                f"[WARN] Index file {ndx_path} does not contain both "
                                "'Polymer' and 'NonPolymer' groups. Consider --tc-grps-groups "
                                "to match your actual group names.",
                                file=sys.stderr,
                            )
                    except OSError:
                        print(
                            f"[WARN] Failed to read ndx file for early group check: {ndx_path}",
                            file=sys.stderr,
                        )
    
    # =========================================================================
    # A2) Validate tau values are positive floats (>0)
    # =========================================================================
    if args.tau_t_polymer is not None and args.tau_t_polymer <= 0:
        parser.error(f"--tau-t-polymer must be > 0, got {args.tau_t_polymer}")
    if args.tau_t_nonpolymer is not None and args.tau_t_nonpolymer <= 0:
        parser.error(f"--tau-t-nonpolymer must be > 0, got {args.tau_t_nonpolymer}")
    if getattr(args, "tc_grps_min_atoms", 50) is not None and args.tc_grps_min_atoms <= 0:
        parser.error("--tc-grps-min-atoms must be > 0.")
    if getattr(args, "tc_grps_min_fraction", 0.01) is not None:
        frac = float(args.tc_grps_min_fraction)
        if frac <= 0.0 or frac >= 1.0:
            parser.error("--tc-grps-min-fraction must be in (0, 1).")
    
    # =========================================================================
    # A3) comm-grps-policy explicit requires --comm-grps
    # =========================================================================
    if args.comm_grps_policy == "explicit" and not args.comm_grps_explicit:
        parser.error(
            "--comm-grps-policy explicit requires --comm-grps <groups> to be set.\n"
            "Example: --comm-grps-policy explicit --comm-grps 'Polymer NonPolymer'"
        )
    
    # =========================================================================
    # A4) --comm-grps without policy=explicit is an error (strict contract)
    # =========================================================================
    if args.comm_grps_explicit and args.comm_grps_policy != "explicit":
        parser.error(
            f"--comm-grps was set but --comm-grps-policy is '{args.comm_grps_policy}'.\n"
            "Either use --comm-grps-policy explicit, or remove --comm-grps."
        )
    
    # =========================================================================
    # B) Resume/force validation: --force + --no-resume is redundant/ambiguous
    # =========================================================================
    if args.force and not args.resume:
        parser.error(
            "--force and --no-resume cannot both be set.\n"
            "--force means 'rerun selected stage(s) even if done.ok exists'.\n"
            "--no-resume means 'do not skip any stages'.\n"
            "These are redundant; choose one based on your intent."
        )

    # =========================================================================
    # B2) Restart recovery modes must be unambiguous
    # =========================================================================
    if args.allow_velocity_reset and args.allow_gro_velocity_restart:
        parser.error(
            "--allow-velocity-reset and --allow-gro-velocity-restart are mutually exclusive.\n"
            "Choose ONE restart mode:\n"
            "  --allow-velocity-reset: regenerate velocities from scratch (gen_vel=yes).\n"
            "  --allow-gro-velocity-restart: restart from .gro velocities "
            "(gen_vel=no; requires valid .gro velocities).\n"
            "Use only one flag based on your recovery intent."
        )

    # =========================================================================
    # C) Numeric guardrails for CLI overrides and thresholds
    # =========================================================================
    if getattr(args, "gmx_nt", None) is not None and args.gmx_nt <= 0:
        parser.error("--gmx-nt must be > 0.")
    if getattr(args, "omp_num_threads", None) is not None and args.omp_num_threads <= 0:
        parser.error("--omp-num-threads must be > 0.")
    if getattr(args, "grompp_maxwarn", 0) < 0:
        parser.error("--grompp-maxwarn must be >= 0.")

    for flag_name in ("nsteps_em", "nsteps_nvt", "nsteps_npt", "nsteps_md"):
        nsteps_value = getattr(args, flag_name, None)
        if nsteps_value is not None and nsteps_value <= 0:
            parser.error(f"--{flag_name.replace('_', '-')} must be > 0.")

    # =========================================================================
    # C0) Dispatcher QC policy/threshold validation
    # =========================================================================
    qc_policy = getattr(args, "qc_policy", None)
    if qc_policy is None:
        qc_policy = "warn" if args.stage == "all" else "off"
    args.qc_policy = qc_policy

    qc_density_rel_tol = getattr(args, "qc_density_rel_tol", 0.05)
    if qc_density_rel_tol is None:
        qc_density_rel_tol = 0.05
    try:
        qc_density_rel_tol = float(qc_density_rel_tol)
    except (TypeError, ValueError):
        parser.error("--qc-density-rel-tol must be a float in (0, 1).")
    if qc_density_rel_tol <= 0.0 or qc_density_rel_tol >= 1.0:
        parser.error("--qc-density-rel-tol must be in (0, 1).")
    args.qc_density_rel_tol = qc_density_rel_tol

    qc_enable_raw = getattr(args, "qc_enable_after_stages", "gmx_eq")
    if isinstance(qc_enable_raw, str):
        parsed_qc_stages = [tok.strip() for tok in qc_enable_raw.split(",") if tok.strip()]
    elif isinstance(qc_enable_raw, (list, tuple, set)):
        parsed_qc_stages = [str(tok).strip() for tok in qc_enable_raw if str(tok).strip()]
    else:
        parser.error(
            "--qc-enable-after-stages must be a comma-separated stage list "
            "(e.g., 'gmx_eq,htpolynet')."
        )

    valid_qc_stages = set(STAGE_ORDER)
    invalid_qc = [name for name in parsed_qc_stages if name not in valid_qc_stages]
    if invalid_qc:
        parser.error(
            "--qc-enable-after-stages contains invalid stage(s): "
            f"{', '.join(invalid_qc)}. Valid values: {', '.join(STAGE_ORDER)}"
        )
    args.qc_enable_after_stages = parsed_qc_stages

    if qc_policy == "error":
        if qc_density_rel_tol <= 0.0 or qc_density_rel_tol >= 0.5:
            parser.error(
                "--qc-policy=error requires --qc-density-rel-tol in (0, 0.5). "
                f"Got {qc_density_rel_tol}."
            )

    neutrality_thresholds = (
        ("--charge-neutrality-tol", getattr(args, "charge_neutrality_tol", None)),
        ("--charge-neutrality-warn", getattr(args, "charge_neutrality_warn", None)),
        ("--charge-neutrality-correct", getattr(args, "charge_neutrality_correct", None)),
    )
    for flag_name, threshold in neutrality_thresholds:
        if threshold is not None and threshold < 0:
            parser.error(f"{flag_name} must be >= 0.")

    charge_tol = getattr(args, "charge_neutrality_tol", None)
    charge_warn = getattr(args, "charge_neutrality_warn", None)
    charge_correct = getattr(args, "charge_neutrality_correct", None)
    if None not in (charge_tol, charge_warn, charge_correct):
        if not (charge_tol <= charge_warn <= charge_correct):
            parser.error(
                "Charge-neutrality thresholds must satisfy:\n"
                "  --charge-neutrality-tol <= --charge-neutrality-warn <= --charge-neutrality-correct\n"
                f"Got: tol={charge_tol}, warn={charge_warn}, correct={charge_correct}\n"
                "Adjust the thresholds so they increase from pass tolerance "
                "to warning threshold to correction limit."
            )
    
    # =========================================================================
    # D) --run-id explicitness: require --run-id or --allow-auto-run-id
    # =========================================================================
    if args.run_id is None and not getattr(args, 'allow_auto_run_id', False):
        parser.error(
            "--run-id is required for reproducibility and auditing.\n"
            "Provide an explicit run ID, e.g., --run-id run_001\n"
            "Or use --allow-auto-run-id to permit auto-generated IDs."
        )
    if args.run_id is not None and getattr(args, 'allow_auto_run_id', False):
        print(
            "[WARN] --allow-auto-run-id is ignored because --run-id was provided.",
            file=sys.stderr,
        )

    # =========================================================================
    # D2) Forcefield/charge compatibility matrix (fail-fast by default)
    # =========================================================================
    requested_pair = (args.ff, args.charge)
    pair_supported = requested_pair in SUPPORTED_FF_CHARGE_PAIRS
    args.ff_charge_pair_supported = pair_supported
    args.ff_charge_doc_section = FF_CHARGE_COMPAT_DOC_SECTION

    if not pair_supported:
        supported_for_ff = [p for p in SUPPORTED_FF_CHARGE_PAIRS if p[0] == args.ff]
        if supported_for_ff:
            suggestion = _format_supported_pairs(supported_for_ff)
        else:
            suggestion = _format_supported_pairs(SUPPORTED_FF_CHARGE_PAIRS)

        if not getattr(args, "allow_charge_model_mismatch", False):
            parser.error(
                "Unsupported forcefield/charge-model pair: "
                f"{args.ff}/{args.charge}.\n"
                f"Supported pair(s): {suggestion}\n"
                f"See {FF_CHARGE_COMPAT_DOC_SECTION}.\n"
                "If you intentionally want an unvalidated pair, set "
                "--allow-charge-model-mismatch and provide "
                "--charge-model-mismatch-reason."
            )

        reason = (getattr(args, "charge_model_mismatch_reason", "") or "").strip()
        if not reason:
            parser.error(
                "Unsafe override requires --charge-model-mismatch-reason.\n"
                "Provide a brief reason for auditability when using "
                "--allow-charge-model-mismatch."
            )
        args.charge_model_mismatch_reason = reason
        print(
            "[WARN] ===========================================================\n"
            "[WARN] UNSAFE OVERRIDE: unvalidated forcefield/charge-model pair\n"
            f"[WARN] Requested: {args.ff}/{args.charge}\n"
            f"[WARN] Supported pair(s): {suggestion}\n"
            f"[WARN] This is no longer a validated parameter set. "
            "Nonbonded terms and charges may be incompatible.\n"
            "[WARN] ===========================================================",
            file=sys.stderr,
        )
    else:
        if getattr(args, "allow_charge_model_mismatch", False):
            print(
                "[WARN] --allow-charge-model-mismatch ignored: requested pair is supported.",
                file=sys.stderr,
            )
        if getattr(args, "charge_model_mismatch_reason", None):
            print(
                "[WARN] --charge-model-mismatch-reason ignored: requested pair is supported.",
                file=sys.stderr,
            )

    # =========================================================================
    # D3) Mixed [defaults] override must include an explicit audit reason
    # =========================================================================
    mixed_defaults_reason = (getattr(args, "allow_mixed_defaults_reason", "") or "").strip()
    if getattr(args, "allow_mixed_defaults", False):
        if not mixed_defaults_reason:
            parser.error(
                "Unsafe override requires --allow-mixed-defaults-reason.\n"
                "Provide a brief reason for auditability when using --allow-mixed-defaults."
            )
        args.allow_mixed_defaults_reason = mixed_defaults_reason
    else:
        if mixed_defaults_reason:
            print(
                "[WARN] --allow-mixed-defaults-reason ignored: --allow-mixed-defaults is not set.",
                file=sys.stderr,
            )
        args.allow_mixed_defaults_reason = None

    if getattr(args, "mixed_defaults_preserve_within_group", False) and not getattr(
        args, "allow_mixed_defaults", False
    ):
        parser.error(
            "--mixed-defaults-preserve-within-group requires --allow-mixed-defaults.\n"
            "This guard rail only applies when mixed [defaults] override is explicitly enabled."
        )

    cross_group_policy = getattr(args, "mixed_defaults_cross_group_policy", None)
    if cross_group_policy is None:
        cross_group_policy = "warn" if getattr(args, "allow_mixed_defaults", False) else "off"
    args.mixed_defaults_cross_group_policy = cross_group_policy

    if cross_group_policy != "off" and not getattr(args, "allow_mixed_defaults", False):
        parser.error(
            "--mixed-defaults-cross-group-policy requires --allow-mixed-defaults "
            "(or set policy=off)."
        )
    if getattr(args, "mixed_defaults_cross_group_max_pairs", 0) <= 0:
        parser.error("--mixed-defaults-cross-group-max-pairs must be > 0.")

    cross_group_reason = (
        getattr(args, "mixed_defaults_cross_group_reason", "") or ""
    ).strip()
    if cross_group_policy == "generate":
        if not cross_group_reason:
            parser.error(
                "--mixed-defaults-cross-group-policy=generate requires "
                "--mixed-defaults-cross-group-reason."
            )
        args.mixed_defaults_cross_group_reason = cross_group_reason
    else:
        if cross_group_reason:
            print(
                "[WARN] --mixed-defaults-cross-group-reason ignored unless "
                "--mixed-defaults-cross-group-policy=generate.",
                file=sys.stderr,
            )
        args.mixed_defaults_cross_group_reason = None
        if getattr(args, "mixed_defaults_cross_group_rule", None):
            print(
                "[WARN] --mixed-defaults-cross-group-rule is ignored unless "
                "--mixed-defaults-cross-group-policy=generate.",
                file=sys.stderr,
            )

    prefix_policy = getattr(args, "prefix_implicit_topology_policy", "error")
    prefix_reason = (getattr(args, "prefix_implicit_topology_reason", "") or "").strip()
    if prefix_policy in {"bake", "inject"}:
        if not prefix_reason:
            parser.error(
                f"--prefix-implicit-topology-policy={prefix_policy} requires "
                "--prefix-implicit-topology-reason."
            )
        args.prefix_implicit_topology_reason = prefix_reason
    else:
        if prefix_reason:
            print(
                "[WARN] --prefix-implicit-topology-reason ignored unless policy is bake or inject.",
                file=sys.stderr,
            )
        args.prefix_implicit_topology_reason = None

    # =========================================================================
    # E) posres reference path must exist (if provided)
    # =========================================================================
    if getattr(args, "posres_reference_gro", None):
        posres_path = Path(args.posres_reference_gro)
        if not posres_path.exists():
            parser.error(
                f"--posres-reference path does not exist: {posres_path}"
            )

    # =========================================================================
    # E2) packmol-pdb-scale validation
    # =========================================================================
    if getattr(args, "packmol_pdb_scale", None) is not None:
        scale = args.packmol_pdb_scale
        if scale <= 0:
            parser.error("--packmol-pdb-scale must be > 0.")
        # Contract: only 0.1 or 1.0 are supported reliably
        if abs(scale - 0.1) > 1e-6 and abs(scale - 1.0) > 1e-6:
            parser.error(
                "--packmol-pdb-scale must be 0.1 (Å→nm) or 1.0 (nm→nm).\n"
                "Other values are not supported by downstream unit handling."
            )

    # =========================================================================
    # E2b) PACKMOL preassembly validation
    # =========================================================================
    if getattr(args, "packmol_preassembly_li_fraction", 1.0) < 0.0 or getattr(args, "packmol_preassembly_li_fraction", 1.0) > 1.0:
        parser.error("--packmol-preassembly-li-fraction must be within [0, 1].")
    if getattr(args, "packmol_preassembly_retry_count", 0) < 0:
        parser.error("--packmol-preassembly-retry-count must be >= 0.")
    if getattr(args, "packmol_preassembly_sample_cap", 0) <= 0:
        parser.error("--packmol-preassembly-sample-cap must be > 0.")
    for flag_name in (
        "packmol_li_polymer_cutoff_nm",
        "packmol_li_solvent_cutoff_nm",
        "packmol_li_tfsi_cutoff_nm",
    ):
        if getattr(args, flag_name, 0.0) <= 0:
            parser.error(f"--{flag_name.replace('_', '-')} must be > 0.")
    for flag_name in (
        "packmol_target_li_polymer_fraction",
        "packmol_target_li_solvent_fraction",
        "packmol_max_li_tfsi_close_fraction",
    ):
        val = getattr(args, flag_name, 0.0)
        if val < 0.0 or val > 1.0:
            parser.error(f"--{flag_name.replace('_', '-')} must be within [0, 1].")

    # =========================================================================
    # E3) dipole percolation threshold validation
    # =========================================================================
    if getattr(args, "dipole_percolation_threshold", None) is not None:
        threshold = args.dipole_percolation_threshold
        if threshold <= 0 or threshold >= 1:
            parser.error(
                "--dipole-percolation-threshold must be between 0 and 1 (exclusive)."
            )
    if getattr(args, "dipole_group_max_atoms", None) is not None:
        if args.dipole_group_max_atoms <= 0:
            parser.error("--dipole-group-max-atoms must be > 0.")
    if getattr(args, "dipole_unwrap_max_atoms", None) is not None:
        if args.dipole_unwrap_max_atoms <= 0:
            parser.error("--dipole-unwrap-max-atoms must be > 0.")
    if getattr(args, "polymer_net_charge_tol", None) is not None:
        if args.polymer_net_charge_tol < 0:
            parser.error("--polymer-net-charge-tol must be >= 0.")
    if getattr(args, "charge_fix_moleculetype_rounding_tol", None) is not None:
        if args.charge_fix_moleculetype_rounding_tol < 0:
            parser.error("--charge-fix-moleculetype-rounding-tol must be >= 0.")
    if getattr(args, "charge_fix_polymer_exclusion_bonds", None) is not None:
        if args.charge_fix_polymer_exclusion_bonds < 0:
            parser.error("--charge-fix-polymer-exclusion-bonds must be >= 0.")
    if getattr(args, "grompp_preprocess_timeout_s", None) is not None:
        if args.grompp_preprocess_timeout_s <= 0:
            parser.error("--grompp-preprocess-timeout-s must be > 0.")
    if getattr(args, "gmx_eq_metrics_window_ps", None) is not None:
        if args.gmx_eq_metrics_window_ps <= 0:
            parser.error("--gmx-eq-metrics-window-ps must be > 0.")
    if getattr(args, "lj_sigma_max_nm", None) is not None and args.lj_sigma_max_nm <= 0:
        parser.error("--lj-sigma-max-nm must be > 0.")
    if (
        getattr(args, "lj_epsilon_max_kj_mol", None) is not None
        and args.lj_epsilon_max_kj_mol <= 0
    ):
        parser.error("--lj-epsilon-max-kj-mol must be > 0.")
    if getattr(args, "lj_bounds_profile", "aa") == "custom":
        if getattr(args, "lj_sigma_max_nm", None) is None:
            parser.error("--lj-bounds-profile=custom requires --lj-sigma-max-nm.")
        if getattr(args, "lj_epsilon_max_kj_mol", None) is None:
            parser.error("--lj-bounds-profile=custom requires --lj-epsilon-max-kj-mol.")

    # =========================================================================
    # E4) HTPolyNet safety option validation
    # =========================================================================
    if getattr(args, "min_conversion", None) is not None:
        min_conversion = float(args.min_conversion)
        if min_conversion < 0.0 or min_conversion > 1.0:
            parser.error("--min-conversion must be in [0.0, 1.0].")

    timeout_sec = getattr(args, "htpolynet_timeout_sec", None)
    if timeout_sec is not None:
        if timeout_sec < 0:
            parser.error("--htpolynet-timeout must be >= 0 seconds (0 disables timeout).")
        if timeout_sec == 0:
            args.htpolynet_timeout_sec = None

    if getattr(args, "allow_placeholder_stage_to_gromacs", False) and not getattr(
        args, "allow_placeholder_propagate", False
    ):
        parser.error(
            "--allow-placeholder-stage-to-gromacs requires --allow-placeholder-propagate."
        )
    if getattr(args, "allow_placeholder_gromacs_compile", False) and not getattr(
        args, "allow_placeholder_propagate", False
    ):
        parser.error(
            "--allow-placeholder-gromacs-compile requires --allow-placeholder-propagate."
        )

    # =========================================================================
    # F) GPU id validation: sanitize and reject invalid patterns
    # =========================================================================
    if args.gmx_gpu_id is not None:
        gpu_id = args.gmx_gpu_id.strip()
        # Reject whitespace within the string
        if ' ' in gpu_id or '\t' in gpu_id:
            parser.error(
                f"--gmx-gpu-id contains whitespace: '{args.gmx_gpu_id}'.\n"
                "Use comma-separated values without spaces.\n"
                "Valid examples: '0', '0123', '0,1', '0,1,2'"
            )
        # Validate pattern: single int or comma-separated ints (no trailing comma)
        if gpu_id and not re.match(r'^[0-9]+(,[0-9]+)*$', gpu_id):
            parser.error(
                f"--gmx-gpu-id has invalid format: '{args.gmx_gpu_id}'.\n"
                "Expected single GPU ID or comma-separated IDs.\n"
                "Valid examples: '0', '0123', '0,1', '0,1,2'"
            )
        # Store sanitized version (stripped)
        args.gmx_gpu_id = gpu_id if gpu_id else None


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="PACKMOL → HTPOLYNET → GROMACS → RDF/CN Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all stages with GAFF2 forcefield (explicit run-id required)
  python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage all
  
  # Run only PACKMOL stage
  python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage packmol
  
  # Force re-run of a stage (archives prior outputs)
  python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage packmol --force
  
  # Allow auto-generated run ID (for quick tests)
  python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --allow-auto-run-id --stage all
  
  # tc-grps split mode with explicit tau_t values
  python run_pipeline.py --ff GAFF2 --charge CM5 --system SYS_001 --run-id run_001 --stage gmx_eq \\
      --tc-grps-mode split --tau-t-polymer 0.5 --tau-t-nonpolymer 0.1

Stage order: {stages}

Forcefield internal naming:
  GAFF2   → GAFF2 (unchanged)
  OPLS-AA → OPLSAA (dashes removed for path matching)
        """.format(stages=" → ".join(STAGE_ORDER))
    )
    
    # Required arguments
    parser.add_argument(
        "--ff",
        choices=["GAFF2", "OPLS-AA"],
        required=True,
        help="Forcefield to use (OPLS-AA internally normalized to OPLSAA)",
    )
    parser.add_argument(
        "--charge",
        choices=["RESP", "CM5"],
        required=True,
        help="Charge method",
    )
    parser.add_argument(
        "--allow-charge-model-mismatch",
        action="store_true",
        default=False,
        help=(
            "Unsafe override: allow forcefield/charge-model pairs outside the "
            "validated compatibility matrix."
        ),
    )
    parser.add_argument(
        "--charge-model-mismatch-reason",
        type=str,
        default=None,
        metavar="TEXT",
        help="Required audit reason when --allow-charge-model-mismatch is used.",
    )
    
    # Stage selection
    valid_stages = STAGE_ORDER + ["all"]
    parser.add_argument(
        "--stage",
        choices=valid_stages,
        required=True,
        help="Stage to run (or 'all' for complete pipeline)",
    )
    parser.add_argument(
        "--qc-policy",
        choices=["off", "warn", "error"],
        default=None,
        help=(
            "Dispatcher quality-gate policy. Default is 'warn' for --stage all, "
            "otherwise 'off'."
        ),
    )
    parser.add_argument(
        "--qc-density-rel-tol",
        type=float,
        default=0.05,
        metavar="FRAC",
        help="Relative density tolerance for QC against expected density (default: 0.05).",
    )
    parser.add_argument(
        "--qc-enable-after-stages",
        type=str,
        default="gmx_eq",
        metavar="CSV",
        help=(
            "Comma-separated stages where dispatcher QC is evaluated if stage_metrics.json "
            "is present (default: gmx_eq)."
        ),
    )
    parser.add_argument(
        "--all-qc-stop-on-warn",
        action="store_true",
        default=False,
        help="In --stage all mode, stop early after dispatcher QC warning.",
    )
    
    # System ID - now required
    parser.add_argument(
        "--system",
        required=True,
        help="System ID (e.g., SYS_001). Required for path resolution.",
    )
    
    # Run ID - required unless --allow-auto-run-id
    parser.add_argument(
        "--run-id",
        dest="run_id",
        default=None,
        help="Run ID for output directory. Required unless --allow-auto-run-id is set.",
    )
    parser.add_argument(
        "--allow-auto-run-id",
        action="store_true",
        default=False,
        help="Allow auto-generated run ID if --run-id not provided (for quick tests)",
    )
    
    # Resume/force behavior
    # NOTE: --resume flag removed (default is True). Only --no-resume toggles it off.
    parser.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help="Do not skip any stages (default: resume=True, skip completed stages)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force re-run even if done.ok exists (archives prior outputs)",
    )
    
    # Set default for resume (since we removed --resume flag)
    parser.set_defaults(resume=True)
    
    # Additional options
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )
    parser.add_argument(
        "--strict-context-args",
        action="store_true",
        default=False,
        help="Fail if pipeline-related CLI args are present but not consumed by PipelineContext.from_args().",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Project root directory (default: auto-detect)",
    )
    
    # === ITP Sanitizer options ===
    parser.add_argument(
        "--allow-override",
        action="store_true",
        default=False,
        help="Allow atomtype conflicts to be resolved by override (with warning)",
    )
    parser.add_argument(
        "--atomtype-prefix",
        type=str,
        default=None,
        metavar="PREFIX",
        help="Optional prefix for atomtype names (e.g., 'G2_' for GAFF2)",
    )
    parser.add_argument(
        "--prefix-implicit-topology-policy",
        choices=["error", "warn", "bake", "inject"],
        default="error",
        help=(
            "Policy when atomtype prefixing would break implicit/library parameter lookup: "
            "error (default), warn (disable prefixing and continue), bake (materialize explicit bonded params), "
            "or inject (emit prefixed [*types] parameter tables)."
        ),
    )
    parser.add_argument(
        "--prefix-implicit-topology-reason",
        type=str,
        default=None,
        metavar="TEXT",
        help="Required audit reason when --prefix-implicit-topology-policy is bake or inject.",
    )
    parser.add_argument(
        "--strict-charge",
        action="store_true",
        default=False,
        help="Fail on atomtype charge mismatches (default: warn only)",
    )
    parser.add_argument(
        "--allow-unsanitized-grompp",
        action="store_true",
        default=False,
        help="Allow grompp to proceed without sanitizer outputs (NOT recommended).",
    )
    parser.add_argument(
        "--charge-fix-method",
        choices=["safe_subset", "uniform_all"],
        default="safe_subset",
        help="Charge correction method: 'safe_subset' (default, only low-|q| atoms) or 'uniform_all'",
    )
    parser.add_argument(
        "--charge-fix-target-allowlist",
        type=str,
        default=None,
        metavar="MOLS",
        help="Comma-separated molecule names allowed as correction targets (e.g., 'PC,EC,DMC')",
    )
    parser.add_argument(
        "--polymer-net-charge-tol",
        type=float,
        default=1e-3,
        metavar="E",
        help=(
            "Acceptable per-molecule net charge drift for protected polymers. "
            "If |q_polymer| <= tol, correction is skipped with audit (default: 1e-3)."
        ),
    )
    parser.add_argument(
        "--charge-fix-polymer-method",
        choices=["skip_if_small", "spread_safe"],
        default="skip_if_small",
        help=(
            "Protected polymer correction policy: skip small drifts by default, "
            "or use spread_safe when explicit correction is allowlisted."
        ),
    )
    parser.add_argument(
        "--charge-fix-polymer-exclusion-bonds",
        type=int,
        default=2,
        metavar="N",
        help=(
            "For spread_safe polymer correction, exclude atoms within N bonds of hetero atoms "
            "(O/N/S/B/F/P). Default: 2."
        ),
    )
    parser.add_argument(
        "--charge-fix-max-delta-per-atom",
        type=float,
        default=1e-5,
        metavar="E",
        help="Max |Δq| per atom for charge correction (default: 1e-5)",
    )
    parser.add_argument(
        "--charge-fix-max-total",
        type=float,
        default=1e-4,
        metavar="E",
        help="Max |ΔQ| total correction for charge neutrality (default: 1e-4)",
    )
    parser.add_argument(
        "--charge-fix-max-dipole-shift",
        type=float,
        default=0.01,
        metavar="DEBYE",
        help="Max |Δμ| dipole shift for corrected molecule(s) in Debye (default: 0.01)",
    )
    parser.add_argument(
        "--charge-fix-moleculetype-rounding-tol",
        type=float,
        default=1e-4,
        metavar="E",
        help=(
            "Rounding-only gate tolerance for moleculetype charges: "
            "|Q_mol-round(Q_mol)| <= tol (default: 1e-4)."
        ),
    )
    parser.add_argument(
        "--charge-fix-allow-non-rounding",
        action="store_true",
        default=False,
        help=(
            "Unsafe override: allow auto-correction even when moleculetype charge drift "
            "exceeds rounding-only tolerance."
        ),
    )
    parser.add_argument(
        "--strict-charge-physics",
        action="store_true",
        default=False,
        help="Fail if charge correction violates physics guardrails (default: warn)",
    )
    parser.add_argument(
        "--strict-gro-top-check",
        action="store_true",
        default=False,
        help="Fail if GRO/TOP atom count check is uncertain (default: warn/skip)",
    )
    parser.add_argument(
        "--allow-default-defaults",
        action="store_true",
        default=False,
        help="Allow fallback [ defaults ] if none found in inputs",
    )
    parser.add_argument(
        "--allow-mixed-defaults",
        action="store_true",
        default=False,
        help=(
            "Unsafe override: allow mixed [defaults] tuples across sources "
            "(comb-rule/fudge mismatch)."
        ),
    )
    parser.add_argument(
        "--allow-mixed-defaults-reason",
        type=str,
        default=None,
        metavar="TEXT",
        help="Required audit reason when --allow-mixed-defaults is used.",
    )
    parser.add_argument(
        "--mixed-defaults-preserve-within-group",
        action="store_true",
        default=False,
        help=(
            "Guard rail for mixed comb-rule runs: emit explicit [ nonbond_params ] for "
            "within-secondary-group LJ pairs (cross-group pairs still use global defaults)."
        ),
    )
    parser.add_argument(
        "--mixed-defaults-cross-group-policy",
        choices=["off", "warn", "generate"],
        default=None,
        help=(
            "Policy for mixed comb-rule (2 vs 3) cross-group LJ handling: "
            "off, warn, or generate explicit cross-group [ nonbond_params ]. "
            "Default: off unless --allow-mixed-defaults is set, then warn."
        ),
    )
    parser.add_argument(
        "--mixed-defaults-cross-group-rule",
        choices=["lorentz-berthelot", "geometric"],
        default=None,
        help=(
            "Cross-group sigma mixing rule when policy=generate. "
            "Default is chosen from canonical comb-rule: 2->lorentz-berthelot, 3->geometric."
        ),
    )
    parser.add_argument(
        "--mixed-defaults-cross-group-max-pairs",
        type=int,
        default=20000,
        metavar="N",
        help="Maximum generated cross-group nonbond pair count (default: 20000).",
    )
    parser.add_argument(
        "--mixed-defaults-cross-group-reason",
        type=str,
        default=None,
        metavar="TEXT",
        help="Required audit reason when --mixed-defaults-cross-group-policy=generate.",
    )
    parser.add_argument(
        "--strict-forcefield-consistency",
        action="store_true",
        default=False,
        help="Fail on mixed forcefield atomtype conflicts",
    )
    parser.add_argument(
        "--strict-top-update",
        action="store_true",
        default=False,
        help="Fail if system.top has malformed sanitizer managed block markers (default: warn and leave unchanged)",
    )
    parser.add_argument(
        "--strict-include-resolution",
        action="store_true",
        default=False,
        help="Fail on unresolved sanitizer #include targets (default: warn).",
    )
    parser.add_argument(
        "--allow-include-shadowing",
        action="store_true",
        default=False,
        help="Allow duplicate include candidates (shadowing) without raising strict-mode errors.",
    )
    parser.add_argument(
        "--itp-include-dirs-priority",
        choices=["sanitized_first", "forcefield_first", "first", "last"],
        default="forcefield_first",
        help=(
            "Include-dir precedence: forcefield_first/sanitized_first "
            "(legacy aliases: last/first)."
        ),
    )
    parser.add_argument(
        "--allow-unsafe-include-escape",
        action="store_true",
        default=False,
        help=(
            "Allow #include targets that resolve outside configured include roots. "
            "Unsafe; disabled by default."
        ),
    )
    parser.add_argument(
        "--grompp-preprocess-timeout-s",
        dest="grompp_preprocess_timeout_s",
        type=int,
        default=None,
        metavar="SEC",
        help="Base timeout in seconds for sanitizer grompp -pp preprocessing (overrides env var).",
    )
    parser.add_argument(
        "--strict-lj-validation",
        dest="strict_lj_validation",
        action="store_true",
        default=False,
        help="Fail on sanitizer LJ outlier checks (default: warn).",
    )
    parser.add_argument(
        "--no-strict-lj-validation",
        dest="strict_lj_validation",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--lj-outlier-thresholds",
        type=str,
        default=None,
        metavar="JSON_OR_KV",
        help="Override LJ outlier thresholds via JSON or key=value pairs (e.g. 'sigma_max_nm=2.5,epsilon_max_kj_mol=500').",
    )
    parser.add_argument(
        "--lj-outlier-policy",
        choices=["off", "warn", "error"],
        default="warn",
        help="LJ outlier handling policy (default: warn). Unit-screaming thresholds always hard-fail.",
    )
    parser.add_argument(
        "--lj-bounds-profile",
        choices=["aa", "coarse_grained", "polarizable", "custom"],
        default="aa",
        help=(
            "LJ outlier bound profile: aa (default), coarse_grained, polarizable, or custom."
        ),
    )
    parser.add_argument(
        "--lj-sigma-max-nm",
        type=float,
        default=None,
        metavar="NM",
        help="Custom LJ sigma upper bound in nm (used only when --lj-bounds-profile=custom).",
    )
    parser.add_argument(
        "--lj-epsilon-max-kj-mol",
        type=float,
        default=None,
        metavar="KJ",
        help=(
            "Custom LJ epsilon upper bound in kJ/mol (used only when --lj-bounds-profile=custom)."
        ),
    )
    parser.add_argument(
        "--dipole-triclinic-policy",
        choices=["skip", "error", "mic"],
        default="skip",
        help="Triclinic dipole policy: skip (default), error, or mic.",
    )
    parser.add_argument(
        "--strict-dipole-check",
        action="store_true",
        default=False,
        help="Enable strict dipole-check enforcement when preconditions pass.",
    )
    parser.add_argument(
        "--no-strict-dipole-check",
        dest="strict_dipole_check",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--dipole-unwrap-max-atoms",
        type=int,
        default=512,
        metavar="N",
        help="Skip dipole enforcement for molecule groups above N atoms (default: 512).",
    )
    parser.add_argument(
        "--charge-fix-allow-solvents",
        action="store_true",
        default=False,
        help="Allow charge correction on protected solvent molecules (unsafe override).",
    )
    parser.add_argument(
        "--charge-fix-protect-resnames",
        type=str,
        default=None,
        metavar="LIST",
        help="Comma-separated residue names to protect from charge-fix edits.",
    )
    parser.add_argument(
        "--dipole-percolation-threshold",
        type=float,
        default=0.45,
        metavar="FRAC",
        help="Percolation gating threshold for dipole check span ratio (default: 0.45)",
    )
    parser.add_argument(
        "--dipole-group-max-atoms",
        type=int,
        default=2048,
        metavar="N",
        help="Skip dipole work for any residue group larger than N atoms (default: 2048)",
    )
    
    # === New Charge Neutrality Hardening Options (Part 2) ===
    parser.add_argument(
        "--strict-charge-neutrality",
        action="store_true",
        default=False,
        help="Fail-fast on preprocessor directives, ambiguous parsing, or count mismatches (default: warn and skip)",
    )
    parser.add_argument(
        "--charge-neutrality-tol",
        type=float,
        default=1e-6,
        metavar="TOL",
        help="Tolerance for neutrality pass: |Q| < TOL → neutral (default: 1e-6)",
    )
    parser.add_argument(
        "--charge-neutrality-warn",
        type=float,
        default=1e-4,
        metavar="WARN",
        help="Threshold for charge warning (default: 1e-4)",
    )
    parser.add_argument(
        "--charge-neutrality-correct",
        type=float,
        default=1e-4,
        metavar="CORRECT",
        help="Max |Q| for auto-correction (default: 1e-4)",
    )
    parser.add_argument(
        "--charge-fix-target-molecules",
        type=str,
        default=None,
        metavar="MOLS",
        help="Comma-separated molecule names to target for correction (overrides allowlist)",
    )
    parser.add_argument(
        "--charge-fix-target-atomnames",
        type=str,
        default=None,
        metavar="NAMES",
        help="Comma-separated atom name patterns allowed for correction (e.g., 'H,HC,HA')",
    )
    parser.add_argument(
        "--charge-fix-disallowed-atomtypes",
        type=str,
        default=None,
        metavar="TYPES",
        help="Comma-separated atomtype patterns to never adjust (e.g., 'O,N,S')",
    )
    parser.add_argument(
        "--charge-fix-allow-ions",
        action="store_true",
        default=False,
        help="Allow correction on molecules detected as ions (advisory override)",
    )
    parser.add_argument(
        "--charge-fix-source-counts",
        choices=["auto", "manifest", "system_top"],
        default="auto",
        help="Source for molecule counts: 'auto' (prefer system.top), 'manifest', or 'system_top' (default: auto)",
    )
    parser.add_argument(
        "--itp-include-dirs",
        type=str,
        default=None,
        metavar="DIRS",
        help="Comma-separated directories to search for #include targets in ITP files",
    )

    # === GROMACS MDP Patching options ===
    parser.add_argument(
        "-T", "--temperature",
        dest="temperature",
        type=float,
        default=None,
        metavar="KELVIN",
        help="Temperature override for ref_t and gen_temp in MDP files",
    )
    parser.add_argument(
        "-P", "--pressure",
        dest="pressure",
        type=float,
        default=None,
        metavar="BAR",
        help="Pressure override for ref_p in MDP files",
    )
    parser.add_argument(
        "--tc-grps-mode",
        choices=["auto", "system", "split"],
        default="auto",
        help=(
            "tc-grps mode: 'auto' (recommended; split Polymer/NonPolymer when safe), "
            "'system' (force System), or 'split' (force split)."
        ),
    )
    parser.add_argument(
        "--nsteps-em",
        type=int,
        default=None,
        help="Override nsteps for energy minimization",
    )
    parser.add_argument(
        "--nsteps-nvt",
        type=int,
        default=None,
        help="Override nsteps for NVT equilibration",
    )
    parser.add_argument(
        "--nsteps-npt",
        type=int,
        default=None,
        help="Override nsteps for NPT equilibration",
    )
    parser.add_argument(
        "--nsteps-md",
        type=int,
        default=None,
        help="Override nsteps for production MD",
    )
    parser.add_argument(
        "--gmx-eq-metrics-window-ps",
        type=float,
        default=20.0,
        metavar="PS",
        help="Window size in ps for gmx_eq thermodynamic metric averages (default: 20.0).",
    )
    
    # === tc-grps split safety options (Issue A) ===
    parser.add_argument(
        "--tau-t-polymer",
        type=float,
        default=None,
        metavar="PS",
        help="tau_t for Polymer group when tc-grps split (ps). Required with --tc-grps-mode split.",
    )
    parser.add_argument(
        "--tau-t-nonpolymer",
        type=float,
        default=None,
        metavar="PS",
        help="tau_t for NonPolymer group when tc-grps split (ps). Required with --tc-grps-mode split.",
    )
    parser.add_argument(
        "--tau-t-values",
        type=str,
        default=None,
        metavar="LIST",
        help="Comma-separated tau_t list for split groups (e.g., '0.1,0.5,2.0'). "
             "Count must match --tc-grps-groups.",
    )
    parser.add_argument(
        "--allow-tau-t-autofill",
        action="store_true",
        default=False,
        help="Allow duplicating scalar tau_t for split mode (legacy, not recommended for GPE)",
    )
    
    # === COM removal policy options (Issue B) ===
    parser.add_argument(
        "--comm-grps-policy",
        choices=["auto", "match_tc", "system", "none", "explicit"],
        default="match_tc",
        help=(
            "COM group policy: 'auto' (risk-aware), 'match_tc' (match tc-grps, default), "
            "'system' (global COM removal), 'none' (comm-mode=None), or 'explicit'"
        ),
    )
    parser.add_argument(
        "--comm-mode",
        choices=["Linear", "Angular", "None"],
        default="Linear",
        help="COM removal mode (default: Linear)",
    )
    parser.add_argument(
        "--allow-comm-mode-angular-with-pbc",
        action="store_true",
        default=False,
        help="Expert override: allow comm-mode=Angular while PBC is active (unsafe; default blocks/overrides).",
    )
    parser.add_argument(
        "--comm-grps",
        dest="comm_grps_explicit",
        type=str,
        default=None,
        metavar="GROUPS",
        help="Explicit comm-grps value. Only valid with --comm-grps-policy explicit.",
    )
    
    # === Stage velocity handling (Issue C) ===
    parser.add_argument(
        "--force-gen-vel",
        action="store_true",
        default=False,
        help="Force gen_vel=yes even in production stage (not recommended)",
    )
    
    # === Reproducibility: velocity generation seed ===
    parser.add_argument(
        "--gen-seed",
        dest="gen_seed",
        type=int,
        default=None,
        metavar="SEED",
        help="Deterministic seed for velocity generation (overrides template's gen_seed; only affects NVT)",
    )
    
    # === NPT barostat selection (C-rescale removed - it's a thermostat concept) ===
    parser.add_argument(
        "--npt-barostat",
        dest="npt_barostat",
        choices=["Berendsen", "Parrinello-Rahman"],
        default="Berendsen",
        help="Barostat (pcoupl) for NPT equilibration: 'Berendsen' (stable for stressed gels, default) or 'Parrinello-Rahman' (correct fluctuations)",
    )

    # === Resource scheduling / pinning options ===
    parser.add_argument(
        "--gmx-nt",
        type=int,
        default=None,
        metavar="THREADS",
        help="Number of threads for gmx mdrun (-nt)",
    )
    parser.add_argument(
        "--gmx-gpu-id",
        type=str,
        default=None,
        metavar="ID",
        help="GPU ID(s) for gmx mdrun. Single ID '0' or comma-separated '0,1,2'. No spaces.",
    )
    parser.add_argument(
        "--omp-num-threads",
        type=int,
        default=None,
        metavar="N",
        help="OpenMP threads per rank (sets mdrun -ntomp and OMP_NUM_THREADS)",
    )
    parser.add_argument(
        "--htpolynet-inject-add-missing",
        action="store_true",
        default=False,
        help="Allow HTPolyNet config injection to create missing gromacs.ntomp/gpu_id keys.",
    )

    # === Cleanup options ===
    parser.add_argument(
        "--cleanup-success",
        action="store_true",
        default=False,
        help="Enable post-success cleanup (TRR kept by default; use --cleanup-delete-trr to remove)",
    )
    parser.add_argument(
        "--cleanup-delete-trr",
        action="store_true",
        default=False,
        help="Explicitly allow TRR deletion during cleanup (velocities needed for Green-Kubo)",
    )
    parser.add_argument(
        "--analysis-requires-velocities",
        action="store_true",
        default=False,
        help="Signal that analysis requires velocities (e.g., Green-Kubo transport); prevents TRR deletion",
    )

    # === Strict reproducibility options ===
    parser.add_argument(
        "--grompp-maxwarn",
        type=int,
        default=0,
        metavar="N",
        help="Maximum warnings allowed by grompp (default: 0, strict)",
    )
    parser.add_argument(
        "--allow-velocity-reset",
        action="store_true",
        default=False,
        help="Allow proceeding without checkpoint (velocities regenerated from scratch)",
    )
    parser.add_argument(
        "--allow-resume-repair",
        action="store_true",
        default=False,
        help="Allow production resume repair when md.cpt exists but md.tpr is missing (rebuild md.tpr, then continue).",
    )
    parser.add_argument(
        "--strict-restart-policy",
        action="store_true",
        default=None,
        help="Enforce strict checkpoint/restart policy independently of --strict-mdp-validation.",
    )
    parser.add_argument(
        "--resume-ignore-runtime",
        action="store_true",
        default=False,
        help="Ignore runtime-only fingerprint keys (e.g., mdrun thread/GPU command base) during production resume comparison.",
    )
    parser.add_argument(
        "--strict-mdp-validation",
        action="store_true",
        default=False,
        help="Escalate MDP validation warnings to errors (comm-grps mismatch, gen_seed missing, thermal gradient)",
    )
    
    # === MDP Patcher Hardening options ===
    parser.add_argument(
        "--tc-grps-groups",
        type=str,
        default=None,
        metavar="GROUPS",
        help="Custom tc-grps group names for split mode (e.g., 'PEO LiTFSI'). Overrides default 'Polymer NonPolymer'.",
    )
    parser.add_argument(
        "--tc-grps-ndx-path",
        type=str,
        default=None,
        metavar="PATH",
        help="Path to .ndx file for auto-detecting tc-grps group names.",
    )
    parser.add_argument(
        "--tc-grps-min-atoms",
        type=int,
        default=50,
        metavar="N",
        help="Minimum atom count per thermostat group in split mode (tiny groups fallback to System; default: 50).",
    )
    parser.add_argument(
        "--tc-grps-min-fraction",
        type=float,
        default=0.01,
        metavar="FRAC",
        help=(
            "Minimum atom fraction per thermostat group in split mode/auto mode "
            "(tiny groups fallback to System; default: 0.01)."
        ),
    )
    parser.add_argument(
        "--allow-small-tc-grps",
        action="store_true",
        default=False,
        help="Expert override: permit split tc-grps even when a group is below --tc-grps-min-atoms.",
    )
    parser.add_argument(
        "--system-type",
        choices=["gel", "solid", "liquid"],
        default=None,
        help=(
            "System type hint for PR barostat safety defaults "
            "(tau_p and compressibility warnings/injection)."
        ),
    )
    parser.add_argument(
        "--allow-continuation-override",
        action="store_true",
        default=False,
        help="Allow gen_vel/continuation mismatch (not recommended, can cause LINCS instability).",
    )
    parser.add_argument(
        "--allow-gro-velocity-restart",
        action="store_true",
        default=False,
        help="Allow gen_vel=no + continuation=no for expert restart from .gro with velocities. "
             "RISK: if .gro lacks velocities, thermostat startup can trigger hotspots/constraint failures.",
    )
    
    # === POSRES Reference Options ===
    parser.add_argument(
        "--posres-reference",
        dest="posres_reference_gro",
        type=str,
        default=None,
        metavar="PATH",
        help="Explicit reference structure for POSRES (-r to grompp). Required when MDP uses -DPOSRES.",
    )
    parser.add_argument(
        "--no-strict-posres-reference",
        dest="strict_posres_reference",
        action="store_false",
        default=True,
        help="Allow sticky auto-pinned POSRES reference in OUT_GMX/<RUN_ID>/03_gromacs/_shared/ when explicit --posres-reference is not set.",
    )

    # === PACKMOL strict reproducibility options ===
    parser.add_argument(
        "--allow-partial-publish",
        action="store_true",
        default=False,
        help="Allow PACKMOL stage to succeed even if publish fails (non-strict mode)",
    )
    parser.add_argument(
        "--allow-density-reduction",
        action="store_true",
        default=False,
        help="Allow PACKMOL retries to reduce density below floor (default 85%%)",
    )
    parser.add_argument(
        "--allow-default-recipe",
        action="store_true",
        default=False,
        help="Allow using default demo recipe when composition.yaml is missing (not for production)",
    )
    parser.add_argument(
        "--min-polymer-chains",
        type=int,
        default=100,
        metavar="N",
        help="Minimum polymer chain count for percolation (default: 100)",
    )
    parser.add_argument(
        "--min-total-molecules-polymer",
        type=int,
        default=5000,
        metavar="N",
        help="Minimum total molecules for polymer-like systems (default: 5000)",
    )
    parser.add_argument(
        "--no-strict-polymer-check",
        dest="strict_polymer_check",
        action="store_false",
        default=True,
        help="Warn instead of fail when polymer percolation thresholds not met",
    )
    parser.add_argument(
        "--packmol-pdb-scale",
        type=float,
        default=None,
        metavar="SCALE",
        help="Override PDB unit detection: 0.1 for Å→nm, 1.0 for nm→nm",
    )
    parser.add_argument(
        "--no-strict-packmol-units",
        dest="strict_packmol_units",
        action="store_false",
        default=True,
        help="Allow ambiguous unit inference to fallback to scale=1.0 (nm) instead of failing.",
    )
    parser.add_argument(
        "--strict-gro-conversion",
        dest="strict_gro_conversion",
        action="store_true",
        default=False,
        help="Fail-fast if PDB->GRO conversion cannot run or editconf fails.",
    )
    parser.add_argument(
        "--no-strict-gro-conversion",
        dest="strict_gro_conversion",
        action="store_false",
        help="Allow non-strict behavior for PDB->GRO conversion failures.",
    )
    parser.add_argument(
        "--packmol-edge-margin",
        type=float,
        default=0.2,
        metavar="NM",
        help="Edge margin for PACKMOL molecule placement (nm, default: 0.2)",
    )
    parser.add_argument(
        "--packmol-preassembly-mode",
        choices=["none", "li_near_polymer", "li_solvent_shell"],
        default="none",
        help="Optional initial-placement bias mode for Li solvation quality (default: none).",
    )
    parser.add_argument(
        "--packmol-preassembly-li-fraction",
        type=float,
        default=1.0,
        metavar="FRAC",
        help="Fraction of Li molecules targeted by preassembly hints (default: 1.0).",
    )
    parser.add_argument(
        "--packmol-preassembly-retry-count",
        type=int,
        default=2,
        metavar="N",
        help="Max post-pack preassembly repair retries with adjusted seed/hints (default: 2).",
    )
    parser.add_argument(
        "--packmol-preassembly-sample-cap",
        type=int,
        default=4000,
        metavar="N",
        help="Max sampled atoms per role for preassembly distance metrics (default: 4000).",
    )
    parser.add_argument(
        "--packmol-li-polymer-cutoff-nm",
        type=float,
        default=0.35,
        metavar="NM",
        help="Cutoff for Li-to-polymer oxygen proximity metric (default: 0.35 nm).",
    )
    parser.add_argument(
        "--packmol-li-solvent-cutoff-nm",
        type=float,
        default=0.35,
        metavar="NM",
        help="Cutoff for Li-to-solvent oxygen proximity metric (default: 0.35 nm).",
    )
    parser.add_argument(
        "--packmol-li-tfsi-cutoff-nm",
        type=float,
        default=0.40,
        metavar="NM",
        help="Cutoff for Li-to-TFSI close-contact proxy metric (default: 0.40 nm).",
    )
    parser.add_argument(
        "--packmol-target-li-polymer-fraction",
        type=float,
        default=0.40,
        metavar="FRAC",
        help="Target minimum fraction of Li within Li-polymer cutoff in li_near_polymer mode.",
    )
    parser.add_argument(
        "--packmol-target-li-solvent-fraction",
        type=float,
        default=0.50,
        metavar="FRAC",
        help="Target minimum fraction of Li within Li-solvent cutoff in li_solvent_shell mode.",
    )
    parser.add_argument(
        "--packmol-max-li-tfsi-close-fraction",
        type=float,
        default=0.30,
        metavar="FRAC",
        help="Target maximum fraction of Li within Li-TFSI close-contact cutoff.",
    )
    parser.add_argument(
        "--allow-composition-changing-retries",
        action="store_true",
        default=False,
        help="Allow PACKMOL retries to reduce molecule counts (DANGEROUS: changes recipe meaning)",
    )

    # === HTPolyNet options ===
    parser.add_argument(
        "--allow-placeholder",
        action="store_true",
        default=False,
        help="Allow HTPolyNet placeholder outputs when tool not installed (demo only)",
    )
    parser.add_argument(
        "--allow-placeholder-propagate",
        action="store_true",
        default=False,
        help="Allow publishing placeholder outputs to IN/systems/<SYSTEM_ID>/htpolynet (demo only).",
    )
    parser.add_argument(
        "--allow-placeholder-stage-to-gromacs",
        action="store_true",
        default=False,
        help="Explicitly stage propagated placeholders to IN/.../gromacs (unsafe demo override).",
    )
    parser.add_argument(
        "--allow-placeholder-gromacs-compile",
        action="store_true",
        default=False,
        help="Disable placeholder poison pill that intentionally breaks grompp (unsafe).",
    )
    parser.add_argument(
        "--allow-triclinic-unsafe-diagonal-approx",
        action="store_true",
        default=False,
        help="Unsafe override: allow triclinic GRO boxes via diagonal approximation for HTPolyNet config.",
    )
    parser.add_argument(
        "--skip-gelation-precheck",
        action="store_true",
        default=False,
        help="Skip conservative gelation feasibility precheck before running HTPolyNet.",
    )
    parser.add_argument(
        "--allow-gelation-precheck-fail",
        action="store_true",
        default=False,
        help="Downgrade failed gelation precheck to warning and continue (unsafe).",
    )
    parser.add_argument(
        "--allow-missing-htpolynet-qc",
        action="store_true",
        default=False,
        help="Allow successful runs with missing conversion/gel-fraction parse (unsafe).",
    )
    parser.add_argument(
        "--min-conversion",
        type=float,
        default=None,
        metavar="FRAC",
        help="Minimum crosslinking conversion threshold (0.0-1.0)",
    )
    parser.add_argument(
        "--htpolynet-timeout",
        dest="htpolynet_timeout_sec",
        type=int,
        default=43200,
        metavar="SEC",
        help="HTPolyNet subprocess timeout in seconds (default: 43200; set 0 to disable).",
    )
    parser.add_argument(
        "--unsafe-allow-out-fallback",
        action="store_true",
        default=False,
        help="Allow reading initial structure from OUT_GMX when IN/published input is missing (unsafe).",
    )

    args = parser.parse_args(argv)
    
    # Validate co-dependent and mutually-exclusive arguments
    validate_args(parser, args)
    
    # Auto-generate run ID if allowed and not provided
    if args.run_id is None and args.allow_auto_run_id:
        args.run_id = generate_run_id()
        # Print to stderr to avoid polluting stdout for scripting
        print(f"Generated run ID: {args.run_id}", file=sys.stderr)
    
    # Normalize forcefield name for internal use
    args.ff_key = normalize_forcefield(args.ff)
    
    return args


def get_project_root(override: Path = None) -> Path:
    """
    Determine the project root directory.

    Priority:
      1) explicit override from --project-root
      2) first parent of this file that looks like repo root
         (contains run_pipeline.py, .git, or IN/)
      3) legacy fallback (this_file.parent.parent)
    """
    if override is not None:
        return override.resolve()

    this_file = Path(__file__).resolve()

    # Walk upward and pick the first parent that matches repo-root markers.
    for candidate in this_file.parents:
        if (
            (candidate / "run_pipeline.py").exists()
            or (candidate / ".git").exists()
            or (candidate / "IN").is_dir()
        ):
            return candidate

    # Backward-compatible fallback for older layouts.
    return this_file.parent.parent


def main(argv=None) -> int:
    """
    Main entry point.
    
    Returns:
        0 on success, 1 on error
    """
    args = None
    try:
        args = parse_args(argv)
        
        project_root = get_project_root(args.project_root)
        # Only print project root in verbose mode to reduce noise
        if args.verbose:
            print(f"Project root: {project_root}")
        
        # Create pipeline context
        ctx = PipelineContext.from_args(args, project_root)
        
        # Run the pipeline
        dispatcher = StageDispatcher(ctx)
        success = dispatcher.run()

        if success:
            return 0
        if getattr(ctx, "run_interrupted", False):
            return 130
        return 1
        
    except KeyboardInterrupt:
        print("\n\nInterrupted by user.")
        return 130
    except Exception as e:
        print(f"\n[FATAL ERROR] {e}")
        # Use args.verbose if available, otherwise check sys.argv as fallback
        show_traceback = (args is not None and args.verbose)
        if show_traceback:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
