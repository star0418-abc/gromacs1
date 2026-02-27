"""
GROMACS stages: Energy minimization, equilibration, and production.

Implements:
- MDP dynamic patching (T/P/nsteps/tc-grps)
- Resource pinning (threads, GPU)
- Resume/skip with checkpoint logic
- Crash context preservation
- Strict reproducibility with stage_state.json fingerprinting (Part 1)
"""

import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import socket
import subprocess
import signal
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

from .base import BaseStage

if TYPE_CHECKING:
    from ..context import PipelineContext


# Import utilities
from ..mdp_patcher import (
    prepare_stage_mdp,
    get_mdp_overrides_for_stage,
    resolve_mdp_overrides,
    parse_mdp,
    get_semantic_value,
    parse_ndx_group_sizes,
    GEN_TEMP_REF_T_TOLERANCE_K,
    MDPParseError,
    MDPPatchError,
)
from ..resource_scheduler import (
    build_mdrun_command, build_grompp_command,
    get_execution_env, get_resource_summary,
)
from ..crash_context import handle_stage_failure
from ..gromacs_cmd import resolve_gmx_command


def _gro_has_velocities(gro_path: Path) -> Optional[bool]:
    """
    Check if a GRO file contains velocity columns.
    
    GRO format uses fixed-width columns:
    - Residue number: columns 1-5 (5 chars)
    - Residue name: columns 6-10 (5 chars)  
    - Atom name: columns 11-15 (5 chars)
    - Atom number: columns 16-20 (5 chars)
    - Positions X,Y,Z: columns 21-44 (3×8 chars = 24 chars)
    - Velocities VX,VY,VZ: columns 45-68 (3×8 chars = 24 chars, optional)
    
    A line with velocities has at least 44 chars (positions) + 24 chars (velocities) = 68 chars.
    This implementation is O(1) memory: it reads title, natoms, and only the first atom line.
    
    Args:
        gro_path: Path to GRO file
        
    Returns:
        True if the first atom line has parseable velocity columns.
        False if file is parseable but velocities are absent/invalid.
        None if file cannot be read or is malformed/ambiguous.
    """
    if not gro_path.exists():
        return None

    try:
        with open(gro_path, "r", encoding="utf-8", errors="replace") as fh:
            # Title/comment
            title = fh.readline()
            if title == "":
                return None

            # Atom count
            natoms_line = fh.readline()
            if natoms_line == "":
                return None
            try:
                atom_count = int(natoms_line.strip())
            except ValueError:
                return None
            if atom_count <= 0:
                return False

            # GRO velocities are all-or-none in normal output. Check only the first atom line.
            first_atom_line = fh.readline()
            if first_atom_line == "":
                return None
            line = first_atom_line.rstrip("\r\n")
            if len(line) < 68:
                return False
            try:
                float(line[44:52])
                float(line[52:60])
                float(line[60:68])
            except ValueError:
                return False
            return True
    except OSError:
        return None
    except Exception:
        return None


def _define_tokens_enable_posres(define_value: str) -> bool:
    """Return True if an MDP define RHS enables any POSRES* macro."""
    raw = (define_value or "").strip()
    if not raw:
        return False
    try:
        tokens = shlex.split(raw, posix=True)
    except ValueError:
        tokens = raw.split()

    i = 0
    while i < len(tokens):
        tok = tokens[i].strip().strip("'\"")
        macro = None
        if tok == "-D":
            if i + 1 < len(tokens):
                macro = tokens[i + 1].strip().strip("'\"")
                if macro.startswith("-D"):
                    macro = macro[2:]
                i += 1
        elif tok.startswith("-D"):
            macro = tok[2:]
        if macro and "POSRES" in macro.upper():
            return True
        i += 1
    return False


def _detect_posres_in_mdp(mdp_path: Path) -> bool:
    """
    Check if MDP enables position restraints via any POSRES* define macro.
    
    Handles various valid GROMACS syntaxes:
    - define = -DPOSRES
    - define=-DPOSRES
    - define = -D POSRES
    - define = -DSOME -DPOSRES
    
    Args:
        mdp_path: Path to MDP file
        
    Returns:
        True if define contains POSRES macro tokens (e.g. -DPOSRES_HEAVY), False otherwise.
    """
    if not mdp_path.exists():
        return False

    try:
        with open(mdp_path, "r", encoding="utf-8", errors="replace") as fh:
            for raw_line in fh:
                # Strip comments (anything after ';')
                line = raw_line.split(";", 1)[0].strip()
                if not line or "=" not in line:
                    continue

                # Parse key = value
                key, _, val = line.partition("=")
                if key.strip().lower() != "define":
                    continue

                if _define_tokens_enable_posres(val):
                    return True
    except OSError:
        return False

    return False


def _get_posres_reference(
    ctx: "PipelineContext",
    stage: str,
    current_structure: Path,
    mdp_path: Path,
) -> Optional[Path]:
    """
    Determine the appropriate POSRES reference structure.
    
    Resolution order:
    1. ctx.posres_reference_gro if explicitly provided
    2. Sticky auto-pinned reference in OUT_GMX/<RUN_ID>/03_gromacs/_shared/posres_reference.gro
       (only when strict_posres_reference=False)
    
    Args:
        ctx: Pipeline context
        stage: Stage name for error messages
        current_structure: Current -c structure path
        mdp_path: Path to MDP file to check for POSRES
        
    Returns:
        Path to reference structure if POSRES active and reference determined,
        None if POSRES not active.
        
    Raises:
        ValueError: If POSRES active with strict mode and no explicit reference.
    """
    if not _detect_posres_in_mdp(mdp_path):
        return None  # No POSRES, no -r needed
    
    # Explicit reference from CLI
    posres_ref = getattr(ctx, 'posres_reference_gro', None) or getattr(ctx, "posres_reference", None)
    if posres_ref:
        ref_path = Path(posres_ref)
        if not ref_path.exists():
            raise FileNotFoundError(f"POSRES reference not found: {ref_path}")
        print(f"  - POSRES reference: {ref_path}")
        return ref_path
    
    # Strict mode: require explicit reference
    strict = getattr(ctx, 'strict_posres_reference', True)
    if strict:
        raise ValueError(
            f"Stage '{stage}' has POSRES enabled but no explicit reference structure. "
            f"Use --posres-reference PATH to specify, or --no-strict-posres-reference "
            f"to enable sticky auto-pinning under OUT_GMX/<RUN_ID>/03_gromacs/_shared/."
        )
    
    # Non-strict: sticky auto-pinned reference in run output.
    # This prevents stage-to-stage reference creep when POSRES is active.
    if ctx.dry_run:
        shared_dir = ctx.get_output_path("03_gromacs", "_shared")
    else:
        shared_dir = ctx.ensure_output_dir("03_gromacs", "_shared")
    pinned_ref = shared_dir / "posres_reference.gro"
    if pinned_ref.exists():
        print(f"  [WARN] POSRES active without explicit --posres-reference; reusing pinned reference: {pinned_ref}")
        return pinned_ref

    if ctx.dry_run:
        print(
            "  [WARN] POSRES active without explicit --posres-reference; dry-run cannot pin "
            "reference file, using current structure for preview only."
        )
        return current_structure

    shutil.copy2(current_structure, pinned_ref)
    print(
        f"  [WARN] POSRES active without explicit --posres-reference; pinned initial reference to {pinned_ref}"
    )
    return pinned_ref


def _sha256_file(path: Path) -> str:
    """Compute SHA256 hash for a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _normalize_path_for_record(ctx: "PipelineContext", path: Optional[Path]) -> Optional[str]:
    """Normalize a path for manifest/stage_state (prefer project-relative)."""
    if path is None:
        return None
    try:
        resolved = Path(path).resolve()
    except Exception:
        return str(path)
    project_root = getattr(ctx, "project_root", None)
    if project_root:
        try:
            return str(resolved.relative_to(Path(project_root).resolve()))
        except Exception:
            pass
    return str(resolved)


def _resolve_recorded_path(ctx: "PipelineContext", raw_path: Optional[str]) -> Optional[Path]:
    """Resolve a path recorded in stage_state (absolute or project-relative)."""
    if not raw_path:
        return None
    candidate = Path(str(raw_path))
    if candidate.is_absolute():
        return candidate
    project_root = getattr(ctx, "project_root", None)
    if project_root:
        return Path(project_root) / candidate
    return candidate


def _resolve_posres_info(
    ctx: "PipelineContext",
    stage: str,
    current_structure: Path,
    mdp_path: Path,
) -> Tuple[Optional[Path], Dict[str, Any]]:
    """Resolve POSRES reference and return info for logging/manifest."""
    active = _detect_posres_in_mdp(mdp_path)
    if not active:
        return None, {"active": False, "reference_path": None, "reason": "not_active"}
    ref_path = _get_posres_reference(ctx, stage, current_structure, mdp_path)
    posres_ref = getattr(ctx, "posres_reference_gro", None) or getattr(ctx, "posres_reference", None)
    if posres_ref:
        reason = "user_provided"
    else:
        reason = "auto_pinned"
    info = {
        "active": True,
        "reference_path": _normalize_path_for_record(ctx, ref_path),
        "reason": reason,
    }
    return ref_path, info


def _parse_float_list(raw: str) -> Optional[List[float]]:
    """Parse whitespace-separated floats; return None on failure."""
    if raw is None:
        return None
    values = []
    for token in str(raw).split():
        try:
            values.append(float(token))
        except ValueError:
            return None
    return values


def _mdp_requests_custom_groups(params: Dict[str, str]) -> bool:
    """
    Detect if MDP requests custom index groups beyond default System.

    Conservative: any non-empty group list not equal to "System" triggers.
    """
    group_keys = [
        "tc-grps", "tc_grps",
        "comm-grps", "comm_grps",
        "energygrps",
        "energygrp-excl", "energygrp_excl",
        "freezegrps",
        "acc-grps", "acc_grps",
    ]
    for key in group_keys:
        if key not in params:
            continue
        raw = params.get(key, "").strip()
        if not raw:
            continue
        tokens = raw.split()
        if not tokens:
            continue
        if len(tokens) > 1:
            return True
        if tokens[0].casefold() != "system":
            return True
    return False


def _assert_ndx_for_mdp(mdp_path: Path, ndx_path: Optional[Path]) -> None:
    """Fail fast if MDP requests custom groups but no .ndx is available."""
    try:
        params = parse_mdp(mdp_path)
    except MDPParseError as e:
        raise ValueError(f"Failed to parse MDP for group checks: {e}") from e
    if _mdp_requests_custom_groups(params) and ndx_path is None:
        raise ValueError(
            "MDP requests custom groups but no index.ndx was found. "
            "Provide an index.ndx or adjust MDP to use default groups."
        )


def _assert_velocity_reset_mdp(
    mdp_path: Path,
    stage: str,
    expected_gen_temp: Optional[str],
) -> None:
    """Ensure velocity reset settings are present and coherent in final MDP."""
    params = parse_mdp(mdp_path)
    continuation = get_semantic_value(params, "continuation")
    gen_vel = get_semantic_value(params, "gen_vel")
    gen_temp = get_semantic_value(params, "gen_temp")

    if continuation is None or continuation.strip().lower() != "no":
        raise MDPPatchError(
            f"Stage '{stage}': velocity reset requires continuation=no in MDP."
        )
    if gen_vel is None or gen_vel.strip().lower() != "yes":
        raise MDPPatchError(
            f"Stage '{stage}': velocity reset requires gen_vel=yes in MDP."
        )
    if gen_temp is None:
        raise MDPPatchError(
            f"Stage '{stage}': velocity reset requires gen_temp to be set in MDP."
        )
    if expected_gen_temp is not None:
        expected_vals = _parse_float_list(str(expected_gen_temp))
        actual_vals = _parse_float_list(str(gen_temp))
        if not expected_vals or not actual_vals:
            raise MDPPatchError(
                f"Stage '{stage}': gen_temp is non-numeric in MDP: '{gen_temp}'."
            )
        if abs(expected_vals[0] - actual_vals[0]) > GEN_TEMP_REF_T_TOLERANCE_K:
            raise MDPPatchError(
                f"Stage '{stage}': gen_temp mismatch after velocity reset "
                f"(expected {expected_vals[0]}K, got {actual_vals[0]}K)."
            )


def _resolve_velocity_reset_overrides(
    ctx: "PipelineContext",
    stage: str,
    template_path: Path,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Resolve overrides to ensure coherent velocity reset semantics.
    
    Enforces: continuation=no, gen_vel=yes, gen_temp present.
    """
    strict = getattr(ctx, "strict_mdp_validation", False)
    params = parse_mdp(template_path)
    overrides: Dict[str, Any] = {
        "continuation": "no",
        "gen_vel": "yes",
    }
    warnings: List[str] = []

    tmpl_gen_vel = get_semantic_value(params, "gen_vel")
    tmpl_cont = get_semantic_value(params, "continuation")
    if tmpl_gen_vel and tmpl_cont:
        if tmpl_gen_vel.strip().lower() == "yes" and tmpl_cont.strip().lower() == "yes":
            msg = (
                f"Template MDP has gen_vel=yes and continuation=yes (contradictory). "
                f"Velocity reset will override to continuation=no, gen_vel=yes."
            )
            if strict:
                raise MDPPatchError(msg)
            warnings.append(msg)

    gen_temp_source = None
    cli_temperature = getattr(ctx, "temperature", None)
    if cli_temperature is not None:
        gen_temp = str(cli_temperature)
        gen_temp_source = "cli_temperature"
    else:
        ref_t_raw = get_semantic_value(params, "ref_t")
        if ref_t_raw:
            ref_t_values = _parse_float_list(ref_t_raw)
            if not ref_t_values:
                raise MDPPatchError(
                    f"ref_t='{ref_t_raw}' is not numeric in {template_path.name}. "
                    f"Set a numeric ref_t or use --temperature."
                )
            if len(ref_t_values) > 1:
                max_dev = max(abs(v - ref_t_values[0]) for v in ref_t_values)
                if max_dev > GEN_TEMP_REF_T_TOLERANCE_K:
                    msg = (
                        f"ref_t has multiple values {ref_t_values} in {template_path.name}. "
                        f"Velocity reset needs a single gen_temp."
                    )
                    if strict:
                        raise MDPPatchError(
                            msg + " Provide --temperature or set a single ref_t."
                        )
                    warnings.append(
                        msg + f" Using the first value {ref_t_values[0]}K."
                    )
            gen_temp = str(ref_t_values[0])
            gen_temp_source = "ref_t"
        else:
            gen_temp_raw = get_semantic_value(params, "gen_temp")
            if not gen_temp_raw:
                raise MDPPatchError(
                    f"Velocity reset requires a target temperature, but {template_path.name} has "
                    f"neither ref_t nor gen_temp. Set ref_t/gen_temp or use --temperature."
                )
            gen_temp_values = _parse_float_list(gen_temp_raw)
            if not gen_temp_values:
                raise MDPPatchError(
                    f"gen_temp='{gen_temp_raw}' is not numeric in {template_path.name}. "
                    f"Set a numeric gen_temp or use --temperature."
                )
            warnings.append(
                f"{template_path.name}: ref_t missing; using template gen_temp={gen_temp_values[0]}K."
            )
            gen_temp = str(gen_temp_values[0])
            gen_temp_source = "template_gen_temp_fallback"

    overrides["gen_temp"] = gen_temp
    info = {
        "gen_temp_source": gen_temp_source,
        "overrides": overrides,
        "warnings": warnings,
    }
    return overrides, info


def decide_prev_state_policy(
    stage_name: str,
    prev_cpt_path: Path,
    strict_mode: bool,
    allow_velocity_reset: bool,
    allow_gro_velocity_restart: bool,
    template_mdp_path: Optional[Path] = None,
    ctx: Optional["PipelineContext"] = None,
    prev_gro_path: Optional[Path] = None,
) -> Tuple[bool, Optional[Dict[str, Any]], str, List[str], Dict[str, Any]]:
    """
    Centralized decision function for checkpoint/velocity handling.
    
    Used by both NPT (missing nvt.cpt) and Production MD (missing npt.cpt) stages.
    
    Decision matrix (matches README):
    | prev_cpt exists | allow_velocity_reset | allow_gro_velocity_restart | strict | Result |
    |-----------------|---------------------|---------------------------|--------|--------|
    | True            | —                   | —                         | —      | velocity_mode="checkpoint" |
    | False           | True                | False                     | —      | velocity_mode="velocity_reset" |
    | False           | False               | True                      | —      | velocity_mode="gro_velocity_restart" (validated) |
    | False           | True                | True                      | —      | ERROR (mutually exclusive) |
    | False           | False               | False                     | True   | ERROR |
    | False           | False               | False                     | False  | Auto-fallback only after GRO velocity validation |
    
    Args:
        stage_name: Stage name for error messages (e.g., "npt", "md")
        prev_cpt_path: Path to the previous stage checkpoint file
        strict_mode: If True, missing checkpoint without flags is an error
        allow_velocity_reset: If True, regenerate velocities with gen_vel=yes
        allow_gro_velocity_restart: If True, trust .gro velocities (gen_vel=no, continuation=no)
        template_mdp_path: Path to MDP template (needed for velocity reset overrides)
        ctx: Pipeline context (needed for velocity reset temperature resolution)
        prev_gro_path: Path to input .gro file (for velocity validation when using gro_velocity_restart)
        
    Returns:
        Tuple of:
        - checkpoint_available (bool): Whether checkpoint exists
        - extra_overrides (dict or None): MDP overrides to apply
        - velocity_mode (str): "checkpoint" | "velocity_reset" | "gro_velocity_restart" | "error"
        - warnings (list): Warning/info strings to print
        - policy_metadata (dict): observability metadata for velocity validation
        
    Raises:
        ValueError: If both allow_velocity_reset and allow_gro_velocity_restart are True
    """
    warnings_list: List[str] = []
    policy_metadata: Dict[str, Any] = {
        "gro_velocity_validation_performed": False,
        "gro_velocity_validation_outcome": "not_checked",
        "gro_velocity_validation_path": _normalize_path_for_record(ctx, prev_gro_path),
        "auto_fallback": False,
    }
    
    # Log decision context for observability
    cpt_exists = prev_cpt_path.exists()
    print(f"  [RESTART] Policy decision for '{stage_name}':")
    print(f"    prev_cpt exists: {cpt_exists} ({prev_cpt_path.name})")
    print(f"    allow_velocity_reset: {allow_velocity_reset}")
    print(f"    allow_gro_velocity_restart: {allow_gro_velocity_restart}")
    
    # Case 1: Checkpoint exists - normal continuation
    if cpt_exists:
        print(f"    → velocity_mode: checkpoint (using existing checkpoint)")
        return True, {"continuation": "yes", "gen_vel": "no"}, "checkpoint", [], policy_metadata
    
    # Case 2: Checkpoint missing - need a recovery strategy
    # First, check for mutually exclusive flags (fail-fast on ambiguous config)
    if allow_velocity_reset and allow_gro_velocity_restart:
        raise ValueError(
            "Configuration error: --allow-velocity-reset and --allow-gro-velocity-restart "
            "are mutually exclusive.\n\n"
            "  --allow-velocity-reset:\n"
            "    Regenerates velocities from scratch using gen_vel=yes.\n"
            "    Use when: You want fresh Maxwell-Boltzmann velocities at gen_temp.\n"
            "    Risk: May cause thermal shock and destroy long-relaxation history.\n\n"
            "  --allow-gro-velocity-restart:\n"
            "    Trusts velocities embedded in the .gro file (gen_vel=no).\n"
            "    Use when: You have a .gro with validated velocities from a prior run.\n"
            "    Risk: If .gro lacks velocities, simulation starts near 0K.\n\n"
            "Choose ONE based on your recovery intent."
        )
    
    if allow_velocity_reset:
        # Regenerate velocities: gen_vel=yes, continuation=no
        print(f"    → velocity_mode: velocity_reset (regenerating velocities)")
        warnings_list.append(
            f"Previous checkpoint missing ({prev_cpt_path.name}); velocities will be regenerated.\n"
            f"  [HIGH-RISK] Cross-linked GPE warning: velocity reset can destroy long-relaxation history\n"
            f"  and trigger thermal shock (hotspots/LINCS instability). Use only as explicit recovery."
        )
        
        # Resolve velocity reset overrides
        if template_mdp_path is None or not template_mdp_path.exists():
            return False, None, "error", [
                f"Velocity reset requires MDP template but template not found: {template_mdp_path}"
            ], policy_metadata
        
        try:
            extra_overrides, reset_info = _resolve_velocity_reset_overrides(
                ctx, stage_name, template_mdp_path
            )
        except (MDPPatchError, MDPParseError) as e:
            return False, None, "error", [str(e)], policy_metadata
        
        warnings_list.extend(reset_info.get("warnings", []))
        return False, extra_overrides, "velocity_reset", warnings_list, policy_metadata

    if allow_gro_velocity_restart:
        # Trust .gro velocities: gen_vel=no, continuation=no
        # But first, validate that the GRO actually contains velocities
        if prev_gro_path is not None:
            policy_metadata["gro_velocity_validation_performed"] = True
            gro_has_velocities = _gro_has_velocities(prev_gro_path)
            if gro_has_velocities is False:
                policy_metadata["gro_velocity_validation_outcome"] = "missing_velocities"
                print(f"    → velocity_mode: error (GRO lacks velocities)")
                return False, None, "error", [
                    f"Stage '{stage_name}': GRO velocity restart requested but "
                    f"GRO file lacks velocities: {prev_gro_path}\n\n"
                    f"The .gro file must contain parseable velocity columns on all atom lines.\n"
                    f"Without velocities, simulation would start near 0K (catastrophic).\n\n"
                    f"SOLUTIONS:\n"
                    f"  1. Use --allow-velocity-reset instead to regenerate velocities\n"
                    f"  2. Provide a .gro file that contains velocities\n"
                    f"  3. Re-run the previous stage to generate a proper checkpoint"
                ], policy_metadata
            if gro_has_velocities is None:
                policy_metadata["gro_velocity_validation_outcome"] = "unreadable_or_invalid"
                print(f"    → velocity_mode: error (unable to validate GRO velocities)")
                return False, None, "error", [
                    f"Stage '{stage_name}': GRO velocity restart requested but velocity presence "
                    f"could not be validated from {prev_gro_path}.\n"
                    f"Use --allow-velocity-reset, or provide a valid GRO with parseable velocities."
                ], policy_metadata
            policy_metadata["gro_velocity_validation_outcome"] = "has_velocities"
            print(f"    → velocity_mode: gro_velocity_restart (GRO has velocities: validated)")
            warnings_list.append(
                f"Previous checkpoint missing ({prev_cpt_path.name}); using expert restart from .gro.\n"
                f"  INFO: GRO velocity validation passed ({prev_gro_path.name})."
            )
        else:
            # No GRO path provided for validation - fall back to warning
            policy_metadata["gro_velocity_validation_outcome"] = "skipped_no_prev_gro_path"
            print(f"    → velocity_mode: gro_velocity_restart (GRO velocities: not validated)")
            warnings_list.append(
                f"Previous checkpoint missing ({prev_cpt_path.name}); using expert restart from .gro.\n"
                f"  WARN: GRO velocity validation skipped (no path provided).\n"
                f"  If velocities are missing, the run may start with zero velocities (near-0K behavior)."
            )
        
        extra_overrides = {
            "continuation": "no",
            "gen_vel": "no",
        }
        return False, extra_overrides, "gro_velocity_restart", warnings_list, policy_metadata

    # No flags set - behavior depends on strict mode
    if strict_mode:
        print(f"    → velocity_mode: error (strict mode, no recovery flag set)")
        error_msg = (
            f"Stage '{stage_name}': Previous checkpoint missing ({prev_cpt_path.name}).\n"
            f"  This is a hard error in strict mode because velocities are required for continuation.\n\n"
            f"  SOLUTIONS:\n"
            f"    1. --allow-velocity-reset    Regenerate velocities (gen_vel=yes)\n"
            f"       Risk: May cause thermal shock / destroy relaxation history.\n\n"
            f"    2. --allow-gro-velocity-restart    Trust .gro velocities (expert mode)\n"
            f"       Risk: YOU must ensure the .gro file contains velocities.\n"
            f"       If velocities are missing, simulation starts near 0K.\n\n"
            f"  NOTE: These flags are mutually exclusive. Choose ONE.\n"
            f"  If strict restart policy is enabled, rerun without --strict-restart-policy "
            f"  (or without --strict-mdp-validation when restart policy inherits from it)."
        )
        return False, None, "error", [error_msg], policy_metadata
    
    # Non-strict: automatic fallback to .gro restart is allowed only after velocity validation.
    if prev_gro_path is None:
        print(f"    → velocity_mode: error (non-strict fallback blocked: no prev_gro_path)")
        return False, None, "error", [
            f"Stage '{stage_name}': checkpoint is missing ({prev_cpt_path.name}) and no restart flag was provided.\n"
            f"Non-strict auto-fallback now requires GRO velocity validation, but no prev_gro_path was supplied.\n\n"
            f"SOLUTIONS:\n"
            f"  1. Provide a valid prev_gro_path so fallback can validate velocities\n"
            f"  2. Use --allow-gro-velocity-restart explicitly\n"
            f"  3. Use --allow-velocity-reset explicitly (high-risk thermal shock path)"
        ], policy_metadata

    policy_metadata["gro_velocity_validation_performed"] = True
    gro_has_velocities = _gro_has_velocities(prev_gro_path)
    if gro_has_velocities is None:
        policy_metadata["gro_velocity_validation_outcome"] = "unreadable_or_invalid"
        print(f"    → velocity_mode: error (non-strict fallback blocked: unreadable GRO)")
        return False, None, "error", [
            f"Stage '{stage_name}': checkpoint is missing ({prev_cpt_path.name}), and automatic non-strict fallback "
            f"requires a parseable GRO with velocities.\n"
            f"Could not validate velocities from {prev_gro_path}.\n\n"
            f"SOLUTIONS:\n"
            f"  1. Fix/provide a valid GRO with velocities\n"
            f"  2. Use --allow-velocity-reset explicitly\n"
            f"  3. Re-run the previous stage to regenerate checkpoint + GRO"
        ], policy_metadata
    if gro_has_velocities is False:
        policy_metadata["gro_velocity_validation_outcome"] = "missing_velocities"
        print(f"    → velocity_mode: error (non-strict fallback blocked: GRO lacks velocities)")
        return False, None, "error", [
            f"Stage '{stage_name}': checkpoint is missing ({prev_cpt_path.name}), and automatic non-strict fallback "
            f"cannot proceed because {prev_gro_path} lacks velocities.\n\n"
            f"SOLUTIONS:\n"
            f"  1. Use --allow-velocity-reset explicitly (high-risk thermal shock path)\n"
            f"  2. Provide a GRO file with valid velocities\n"
            f"  3. Re-run the previous stage to produce a checkpoint"
        ], policy_metadata

    policy_metadata["gro_velocity_validation_outcome"] = "has_velocities"
    policy_metadata["auto_fallback"] = True
    print(f"    → velocity_mode: gro_velocity_restart (non-strict auto-fallback after validation)")
    warnings_list.append(
        f"Previous checkpoint missing ({prev_cpt_path.name}); non-strict auto-fallback selected .gro velocity restart.\n"
        f"  WARN: Auto-fallback is using validated GRO velocities from {prev_gro_path.name} "
        f"(gen_vel=no, continuation=no)."
    )
    extra_overrides = {
        "continuation": "no",
        "gen_vel": "no",
    }
    return False, extra_overrides, "gro_velocity_restart", warnings_list, policy_metadata


def _get_git_info(project_root: Optional[Path]) -> Dict[str, Optional[str]]:
    """Return git commit and dirty flag if available."""
    info = {"commit": None, "dirty": None}
    if project_root is None:
        return info
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        info["commit"] = commit
        info["dirty"] = "true" if dirty else "false"
    except Exception:
        return info
    return info


def _split_command_tokens(cmd: Any) -> List[str]:
    """Return command tokens from list/tuple/string input."""
    if cmd is None:
        return []
    if isinstance(cmd, (list, tuple)):
        return [str(part) for part in cmd if part is not None and str(part).strip()]
    raw = str(cmd).strip()
    if not raw:
        return []
    try:
        return shlex.split(raw, posix=True)
    except ValueError:
        return raw.split()


def _resolve_gmx_cmd_used(
    *,
    ctx: Optional["PipelineContext"],
    grompp_cmd: Any = None,
    mdrun_cmd_base: Any = None,
) -> Optional[str]:
    """Resolve the exact GROMACS frontend command used by stage commands."""
    for cmd in (grompp_cmd, mdrun_cmd_base):
        tokens = _split_command_tokens(cmd)
        if tokens:
            return tokens[0]
    return resolve_gmx_command(ctx) if ctx is not None else "gmx"


def _get_gromacs_version(gmx_cmd_used: Optional[str]) -> Optional[str]:
    """Return GROMACS version string from the specific executable used."""
    if not gmx_cmd_used:
        return None
    try:
        result = subprocess.run(
            [gmx_cmd_used, "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in (result.stdout or "").splitlines():
            if "GROMACS" in line:
                return line.strip()
        for line in (result.stderr or "").splitlines():
            if "GROMACS" in line:
                return line.strip()
    except Exception:
        return None
    return None


def _get_gromacs_provenance(
    *,
    ctx: Optional["PipelineContext"],
    grompp_cmd: Any = None,
    mdrun_cmd_base: Any = None,
) -> Dict[str, Optional[str]]:
    """Collect reproducibility metadata for the exact gmx command in use."""
    gmx_cmd_used = _resolve_gmx_cmd_used(
        ctx=ctx,
        grompp_cmd=grompp_cmd,
        mdrun_cmd_base=mdrun_cmd_base,
    )
    gmx_bin_resolved = shutil.which(gmx_cmd_used) if gmx_cmd_used else None
    gmx_version = _get_gromacs_version(gmx_cmd_used)
    return {
        "gmx_cmd_used": gmx_cmd_used,
        "gmx_bin_resolved": gmx_bin_resolved,
        "gromacs_version": gmx_version,
    }


def _stage_state_payload(
    *,
    stage: str,
    gro_path_record: Optional[str],
    top_path_record: Optional[str],
    ndx_path_record: Optional[str],
    mdp_sha256: str,
    gro_sha256: str,
    top_sha256: str,
    ndx_sha256: Optional[str],
    md_tpr_sha256: Optional[str],
    md_cpt_sha256: Optional[str],
    stage_overrides: Dict[str, Any],
    extra_overrides: Dict[str, Any],
    effective_overrides: Dict[str, Any],
    cli_overrides: List[str],
    grompp_maxwarn: int,
    grompp_cmd: str,
    mdrun_cmd_base: str,
    gmx_cmd_used: Optional[str],
    gmx_bin_resolved: Optional[str],
    gromacs_version: Optional[str],
    git_info: Dict[str, Optional[str]],
    allow_velocity_reset: bool,
    allow_gro_velocity_restart: bool = False,
    velocity_mode: Optional[str] = None,
) -> Dict[str, Any]:
    """Build fingerprint payload for deterministic resume checks."""
    return {
        "stage": stage,
        "gro_path": gro_path_record,
        "top_path": top_path_record,
        "ndx_path": ndx_path_record,
        "mdp_sha256": mdp_sha256,
        "gro_sha256": gro_sha256,
        "top_sha256": top_sha256,
        "ndx_sha256": ndx_sha256,
        "md_tpr_sha256": md_tpr_sha256,
        "md_cpt_sha256": md_cpt_sha256,
        "stage_overrides": stage_overrides,
        "extra_overrides": extra_overrides,
        "effective_overrides": effective_overrides,
        "cli_overrides": cli_overrides,
        "grompp_maxwarn": grompp_maxwarn,
        "grompp_command": grompp_cmd,
        "mdrun_command_base": mdrun_cmd_base,
        "gmx_cmd_used": gmx_cmd_used,
        "gmx_bin_resolved": gmx_bin_resolved,
        "gromacs_version": gromacs_version,
        "git_commit": git_info.get("commit"),
        "git_dirty": git_info.get("dirty"),
        "allow_velocity_reset": allow_velocity_reset,
        "allow_gro_velocity_restart": allow_gro_velocity_restart,
        "velocity_mode": velocity_mode,
    }


def _build_stage_state(
    ctx: "PipelineContext",
    stage: str,
    mdp_path: Path,
    gro_path: Path,
    top_path: Path,
    ndx_path: Optional[Path],
    grompp_cmd: Any,
    mdrun_cmd_base: Any,
    grompp_maxwarn: int,
    posres_info: Optional[Dict[str, Any]] = None,
    stage_overrides: Optional[Dict[str, Any]] = None,
    extra_overrides: Optional[Dict[str, Any]] = None,
    effective_overrides: Optional[Dict[str, Any]] = None,
    velocity_mode: Optional[str] = None,
) -> Dict[str, Any]:
    """Build stage_state.json data."""
    mdp_sha = _sha256_file(mdp_path)
    gro_sha = _sha256_file(gro_path)
    top_sha = _sha256_file(top_path)
    ndx_sha = _sha256_file(ndx_path) if ndx_path else None
    md_tpr_sha = None
    md_cpt_sha = None
    if stage == "md":
        md_tpr = mdp_path.parent / "md.tpr"
        md_cpt = mdp_path.parent / "md.cpt"
        md_tpr_sha = _sha256_file(md_tpr) if md_tpr.exists() else None
        md_cpt_sha = _sha256_file(md_cpt) if md_cpt.exists() else None
    stage_overrides = stage_overrides or get_mdp_overrides_for_stage(ctx, stage)
    extra_overrides = extra_overrides or {}
    if effective_overrides is None:
        effective_overrides = dict(stage_overrides)
        effective_overrides.update(extra_overrides)
    overrides_list = [f"{k}={effective_overrides[k]}" for k in sorted(effective_overrides.keys())]
    if isinstance(grompp_cmd, str):
        grompp_cmd_str = grompp_cmd
    else:
        grompp_cmd_str = " ".join(str(c) for c in grompp_cmd)
    if isinstance(mdrun_cmd_base, str):
        mdrun_cmd_base_str = mdrun_cmd_base
    else:
        mdrun_cmd_base_str = " ".join(str(c) for c in mdrun_cmd_base)
    gmx_provenance = _get_gromacs_provenance(
        ctx=ctx,
        grompp_cmd=grompp_cmd,
        mdrun_cmd_base=mdrun_cmd_base,
    )
    git_info = _get_git_info(ctx.project_root)
    payload = _stage_state_payload(
        stage=stage,
        gro_path_record=_normalize_path_for_record(ctx, gro_path),
        top_path_record=_normalize_path_for_record(ctx, top_path),
        ndx_path_record=_normalize_path_for_record(ctx, ndx_path) if ndx_path else None,
        mdp_sha256=mdp_sha,
        gro_sha256=gro_sha,
        top_sha256=top_sha,
        ndx_sha256=ndx_sha,
        md_tpr_sha256=md_tpr_sha,
        md_cpt_sha256=md_cpt_sha,
        stage_overrides=stage_overrides,
        extra_overrides=extra_overrides,
        effective_overrides=effective_overrides,
        cli_overrides=overrides_list,
        grompp_maxwarn=grompp_maxwarn,
        grompp_cmd=grompp_cmd_str,
        mdrun_cmd_base=mdrun_cmd_base_str,
        gmx_cmd_used=gmx_provenance.get("gmx_cmd_used"),
        gmx_bin_resolved=gmx_provenance.get("gmx_bin_resolved"),
        gromacs_version=gmx_provenance.get("gromacs_version"),
        git_info=git_info,
        allow_velocity_reset=ctx.allow_velocity_reset,
        allow_gro_velocity_restart=getattr(ctx, "allow_gro_velocity_restart", False),
        velocity_mode=velocity_mode,
    )
    fingerprint = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()
    data = {
        "stage": stage,
        "run_id": ctx.run_id,
        "system_id": ctx.system_id,
        "created_at": datetime.now().isoformat(),
        "mdp_path": str(mdp_path),
        "gro_path": str(gro_path),
        "top_path": str(top_path),
        "ndx_path": str(ndx_path) if ndx_path else None,
        "mdp_sha256": mdp_sha,
        "gro_sha256": gro_sha,
        "top_sha256": top_sha,
        "ndx_sha256": ndx_sha,
        "md_tpr_sha256": md_tpr_sha,
        "md_cpt_sha256": md_cpt_sha,
        "stage_overrides": stage_overrides,
        "extra_overrides": extra_overrides,
        "effective_overrides": effective_overrides,
        "cli_overrides": overrides_list,
        "grompp_maxwarn": grompp_maxwarn,
        "grompp_command": grompp_cmd_str,
        "mdrun_command_base": mdrun_cmd_base_str,
        "gmx_cmd_used": gmx_provenance.get("gmx_cmd_used"),
        "gmx_bin_resolved": gmx_provenance.get("gmx_bin_resolved"),
        "gromacs_version": gmx_provenance.get("gromacs_version"),
        "git": git_info,
        "allow_velocity_reset": ctx.allow_velocity_reset,
        "allow_gro_velocity_restart": getattr(ctx, "allow_gro_velocity_restart", False),
        "velocity_mode": velocity_mode,
        "fingerprint": fingerprint,
        "fingerprint_payload": payload,
    }
    if posres_info:
        data["posres"] = posres_info
    return data


def _diff_payloads(old: Dict[str, Any], new: Dict[str, Any]) -> List[str]:
    """Return human-readable diffs between two payload dicts."""
    diffs = []
    keys = sorted(set(old.keys()) | set(new.keys()))
    for key in keys:
        if old.get(key) != new.get(key):
            diffs.append(f"{key}: {old.get(key)!r} -> {new.get(key)!r}")
    return diffs


def _is_restart_policy_strict(ctx: "PipelineContext") -> bool:
    """Return strictness for checkpoint/restart policy with backward-compatible fallback."""
    strict_restart_policy = getattr(ctx, "strict_restart_policy", None)
    if strict_restart_policy is None:
        return bool(getattr(ctx, "strict_mdp_validation", False))
    return bool(strict_restart_policy)


def _payload_for_resume_compare(
    payload: Dict[str, Any],
    ignore_runtime: bool,
) -> Tuple[Dict[str, Any], List[str]]:
    """
    Optionally drop runtime-only keys from fingerprint comparison.

    Keeps exact matching by default; opt-in ignore mode only excludes keys that can
    change due to hardware/thread/pinning choices.
    """
    runtime_only_keys = {"mdrun_command_base"}
    if not ignore_runtime:
        return dict(payload), []
    filtered: Dict[str, Any] = {}
    ignored: List[str] = []
    for key, value in payload.items():
        if key in runtime_only_keys:
            ignored.append(key)
            continue
        filtered[key] = value
    ignored.sort()
    return filtered, ignored


def _atomic_write_json(path: Path, data: Dict[str, Any], *, sort_keys: bool = False) -> None:
    """Write JSON atomically (temp file + fsync + replace)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(data, indent=2, sort_keys=sort_keys)
    temp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=str(path.parent),
            prefix=f".{path.name}.tmp.",
            suffix=".json",
            delete=False,
        ) as fh:
            fh.write(payload)
            fh.flush()
            os.fsync(fh.fileno())
            temp_path = Path(fh.name)
        os.replace(temp_path, path)
        temp_path = None

        # Best effort: flush directory metadata for durability.
        dir_flags = os.O_RDONLY
        if hasattr(os, "O_DIRECTORY"):
            dir_flags |= os.O_DIRECTORY
        try:
            dir_fd = os.open(str(path.parent), dir_flags)
            try:
                os.fsync(dir_fd)
            finally:
                os.close(dir_fd)
        except OSError:
            pass
    finally:
        if temp_path is not None and temp_path.exists():
            try:
                temp_path.unlink()
            except OSError:
                pass


def _write_stage_state(path: Path, data: Dict[str, Any], ctx: "PipelineContext") -> None:
    """Write stage_state.json unless dry-run."""
    if ctx.dry_run:
        return
    _atomic_write_json(path, data, sort_keys=False)


def _record_stage_manifest(ctx: "PipelineContext", stage_key: str, data: Dict[str, Any]) -> None:
    """Record stage-specific metadata into manifest."""
    if ctx.manifest is None:
        return
    stages = ctx.manifest.get("stages", {})
    stages[stage_key] = data
    ctx.manifest.set("stages", stages)


def _collect_thermostat_audit(
    ctx: "PipelineContext",
    stage_key: str,
    mdp_path: Path,
    ndx_path: Optional[Path],
) -> Dict[str, Any]:
    """Collect thermostat group/coupling metadata for manifest auditability."""
    params = parse_mdp(mdp_path)
    tcoupl = get_semantic_value(params, "tcoupl")
    tc_grps_raw = get_semantic_value(params, "tc_grps") or "System"
    tau_t_raw = get_semantic_value(params, "tau_t")
    ref_t_raw = get_semantic_value(params, "ref_t")
    tc_grps_tokens = tc_grps_raw.split()
    group_sizes: Dict[str, Optional[int]] = {}
    size_lookup_warnings: List[str] = []
    if ndx_path is not None:
        try:
            ndx_sizes = parse_ndx_group_sizes(ndx_path)
            for grp in tc_grps_tokens:
                group_sizes[grp] = ndx_sizes.get(grp)
        except MDPParseError as e:
            size_lookup_warnings.append(f"failed_to_parse_ndx_group_sizes:{e}")
    else:
        for grp in tc_grps_tokens:
            if grp.casefold() == "system":
                group_sizes[grp] = None
    audit = {
        "stage": stage_key,
        "thermostat_type": tcoupl,
        "tc_grps": tc_grps_raw,
        "tc_grps_tokens": tc_grps_tokens,
        "tc_grps_mode_effective": "split" if len(tc_grps_tokens) > 1 else "system",
        "tau_t": tau_t_raw,
        "tau_t_tokens": tau_t_raw.split() if tau_t_raw else [],
        "ref_t": ref_t_raw,
        "ref_t_tokens": ref_t_raw.split() if ref_t_raw else [],
        "group_sizes": group_sizes,
        "ndx_path": str(ndx_path) if ndx_path else None,
        "auto_mode_decision": getattr(ctx, "tc_grps_auto_decision", None),
        "warnings": size_lookup_warnings,
    }
    return audit


def _record_thermostat_manifest(
    ctx: "PipelineContext",
    stage_key: str,
    thermostat: Dict[str, Any],
) -> None:
    """Persist thermostat audit records in manifest."""
    if ctx.manifest is None:
        return
    existing = ctx.manifest.get("thermostat", {})
    if not isinstance(existing, dict):
        existing = {}
    by_stage = existing.get("by_stage", {})
    if not isinstance(by_stage, dict):
        by_stage = {}
    by_stage[stage_key] = thermostat
    existing["by_stage"] = by_stage
    ctx.manifest.set("thermostat", existing)


def _archive_stage_outputs(stage_dir: Path, deffnm: str, ctx: "PipelineContext") -> None:
    """Archive previous outputs before a forced rerun."""
    if ctx.dry_run:
        print("  [DRY-RUN] Skipping archive of prior outputs")
        return
    if not stage_dir.exists():
        return
    patterns = [
        f"{deffnm}.*",
        "*.tpr",
        "*.cpt",
        "*.log",
        "*.edr",
        "*.xtc",
        "*.trr",
        "mdrun*.log",
        "done.ok",
        "stage_state.json",
    ]
    files: List[Path] = []
    for pattern in patterns:
        files.extend(stage_dir.glob(pattern))
    unique_files = []
    seen = set()
    for path in files:
        if not path.exists() or "_archive" in path.parts:
            continue
        if path in seen:
            continue
        seen.add(path)
        unique_files.append(path)
    if not unique_files:
        return
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    archive_dir = stage_dir / "_archive" / timestamp
    archive_dir.mkdir(parents=True, exist_ok=True)
    for path in unique_files:
        shutil.move(str(path), str(archive_dir / path.name))
    print(f"  [FORCE] Archived {len(unique_files)} items to {archive_dir}")


def _relpath_under_run(run_root: Path, path: Path) -> str:
    """Return path relative to run root if possible."""
    return os.path.relpath(path, run_root)


def _file_entry(run_root: Path, path: Path, role: str) -> Dict[str, Any]:
    """Build a file entry with sha256 and size."""
    return {
        "path": _relpath_under_run(run_root, path),
        "role": role,
        "sha256": _sha256_file(path),
        "bytes": path.stat().st_size,
    }


def _infer_output_role(path: Path) -> Optional[str]:
    """Infer output role from file suffix/name."""
    suffix = path.suffix.lower()
    if suffix == ".tpr":
        return "tpr"
    if suffix == ".cpt":
        return "cpt"
    if suffix == ".gro":
        return "gro"
    if suffix == ".xtc":
        return "xtc"
    if suffix == ".edr":
        return "edr"
    if suffix == ".log":
        return "log"
    return None


def _build_repro_json(
    *,
    ctx: "PipelineContext",
    stage: str,
    stage_dir: Path,
    mdp_path: Path,
    mdp_orig: Dict[str, Any],
    mdp_patched: Dict[str, Any],
    gro_path: Path,
    top_path: Path,
    ndx_path: Optional[Path],
    cpt_path: Optional[Path],
    grompp_cmd: str,
    mdrun_cmd: str,
    resume_allowed: bool,
) -> Dict[str, Any]:
    """Build repro.json payload for a stage."""
    run_root = ctx.get_output_path()
    host = {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
    }
    git_info = _get_git_info(ctx.project_root)
    gmx_provenance = _get_gromacs_provenance(
        ctx=ctx,
        grompp_cmd=grompp_cmd,
        mdrun_cmd_base=mdrun_cmd,
    )
    inputs: List[Dict[str, Any]] = []
    for path, role in [
        (gro_path, "gro"),
        (top_path, "top"),
        (ndx_path, "ndx"),
        (mdp_path, "mdp"),
        (cpt_path, "cpt"),
    ]:
        if path is None:
            continue
        if not path.exists():
            continue
        inputs.append(_file_entry(run_root, path, role))
    inputs.sort(key=lambda x: x["path"])

    diff_preview_items = []
    for key in sorted(set(mdp_orig.keys()) | set(mdp_patched.keys())):
        if mdp_orig.get(key) != mdp_patched.get(key):
            diff_preview_items.append(
                f"{key}: {mdp_orig.get(key)} -> {mdp_patched.get(key)}"
            )
    diff_preview = "; ".join(diff_preview_items[:12]) if diff_preview_items else ""

    outputs: List[Dict[str, Any]] = []
    if stage_dir.exists():
        for file_path in sorted(stage_dir.iterdir()):
            if file_path.is_dir():
                continue
            role = _infer_output_role(file_path)
            if role is None:
                continue
            outputs.append(_file_entry(run_root, file_path, role))

    return {
        "stage": stage,
        "timestamp_utc": datetime.utcnow().isoformat() + "Z",
        "host": host,
        "git": {
            "commit": git_info.get("commit"),
            "is_dirty": True if git_info.get("dirty") == "true" else False,
        },
        "gromacs": {
            "version": gmx_provenance.get("gromacs_version"),
            "gmx_bin": gmx_provenance.get("gmx_bin_resolved"),
            "gmx_cmd_used": gmx_provenance.get("gmx_cmd_used"),
            "gmx_bin_resolved": gmx_provenance.get("gmx_bin_resolved"),
        },
        "inputs": inputs,
        "mdp": {
            "path": _relpath_under_run(run_root, mdp_path),
            "sha256": _sha256_file(mdp_path),
            "patch": {"orig": mdp_orig, "patched": mdp_patched},
            "diff_preview": diff_preview,
        },
        "commands": {"grompp": grompp_cmd, "mdrun": mdrun_cmd},
        "policies": {
            "resume_requested": ctx.resume,
            "resume_allowed": resume_allowed,
            "force": ctx.force,
            "allow_velocity_reset": ctx.allow_velocity_reset,
            "grompp_maxwarn": ctx.grompp_maxwarn,
        },
        "outputs": outputs,
    }


def _write_repro_json(ctx: "PipelineContext", stage_dir: Path, data: Dict[str, Any]) -> None:
    """Write repro.json unless dry-run."""
    if ctx.dry_run:
        return
    repro_path = stage_dir / "repro.json"
    _atomic_write_json(repro_path, data, sort_keys=True)


def _write_run_report(ctx: "PipelineContext") -> None:
    """Aggregate repro.json files and manifest into run_report.json."""
    if ctx.dry_run:
        return
    run_root = ctx.get_output_path()
    report_path = run_root / "run_report.json"
    stages_data = []
    errors: List[Dict[str, str]] = []
    for stage_dir in [
        run_root / "03_gromacs" / "em",
        run_root / "03_gromacs" / "nvt",
        run_root / "03_gromacs" / "npt",
        run_root / "03_gromacs" / "md",
    ]:
        repro_path = stage_dir / "repro.json"
        if repro_path.exists():
            try:
                stages_data.append(json.loads(repro_path.read_text()))
            except Exception as e:
                errors.append({
                    "path": str(repro_path),
                    "error": f"{type(e).__name__}: {e}",
                })
                continue
    manifest_data = ctx.manifest.data if ctx.manifest else {}
    summary = {
        "decisions": manifest_data.get("stages", {}),
    }
    report = {
        "run_id": ctx.run_id,
        "system_id": ctx.system_id,
        "generated_at_utc": datetime.utcnow().isoformat() + "Z",
        "manifest": manifest_data,
        "stages": stages_data,
        "errors": errors,
        "summary": summary,
    }
    _atomic_write_json(report_path, report, sort_keys=True)


def verify_staged_inputs(ctx: "PipelineContext") -> Tuple[Path, Path]:
    """
    Verify that staged inputs exist for GROMACS.
    
    Returns:
        Tuple of (gro_path, top_path)
        
    Raises:
        FileNotFoundError: If required inputs are missing
    """
    gromacs_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs")

    # Look for staged structure
    gro_candidates = [
        gromacs_dir / "gro" / "system_current.gro",
        gromacs_dir / "gro" / "system.gro",
    ]
    
    gro_path = None
    for candidate in gro_candidates:
        if candidate.exists():
            gro_path = candidate
            break
    
    if gro_path is None:
        raise FileNotFoundError(
            f"No staged GRO file found in {gromacs_dir / 'gro'}. "
            "Run htpolynet and sanitizer stages first."
        )
    
    # Look for topology
    top_candidates = [
        gromacs_dir / "top" / "system_current.top",
        gromacs_dir / "top" / "system.top",
    ]
    
    top_path = None
    for candidate in top_candidates:
        if candidate.exists():
            top_path = candidate
            break
    
    if top_path is None:
        raise FileNotFoundError(
            f"No staged TOP file found in {gromacs_dir / 'top'}. "
            "Run htpolynet and sanitizer stages first."
        )

    preflight: Dict[str, Any] = {
        "gro_path": _normalize_path_for_record(ctx, gro_path),
        "top_path": _normalize_path_for_record(ctx, top_path),
        "gro_atom_count": None,
        "top_includes_checked": [],
        "top_missing_includes": [],
        "molecules_section_found": False,
        "status": "pending",
    }

    try:
        with open(gro_path, "r", encoding="utf-8", errors="replace") as gro_fh:
            gro_fh.readline()  # title
            atom_count_raw = gro_fh.readline()
        if atom_count_raw == "":
            raise ValueError(f"{gro_path}: missing atom-count line in GRO file.")
        try:
            atom_count = int(atom_count_raw.strip())
        except ValueError as exc:
            raise ValueError(
                f"{gro_path}: GRO atom-count line is not an integer: '{atom_count_raw.strip()}'."
            ) from exc
        if atom_count <= 0:
            raise ValueError(
                f"{gro_path}: GRO atom-count must be > 0 (got {atom_count})."
            )
        preflight["gro_atom_count"] = atom_count
    except OSError as exc:
        raise ValueError(f"Failed to read GRO for preflight check: {gro_path}: {exc}") from exc

    include_re = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]')
    molecules_found = False
    top_missing_includes: List[str] = []
    top_checked_includes: List[Dict[str, Any]] = []
    try:
        with open(top_path, "r", encoding="utf-8", errors="replace") as top_fh:
            for lineno, raw_line in enumerate(top_fh, start=1):
                line_no_comment = raw_line.split(";", 1)[0].strip()
                if not line_no_comment:
                    continue

                include_match = include_re.match(line_no_comment)
                if include_match:
                    include_target = include_match.group(1).strip()
                    include_path = Path(include_target)
                    if not include_path.is_absolute():
                        include_path = (top_path.parent / include_path).resolve()
                    exists = include_path.exists()
                    top_checked_includes.append({
                        "line": lineno,
                        "include": include_target,
                        "resolved_path": _normalize_path_for_record(ctx, include_path),
                        "exists": exists,
                    })
                    if not exists:
                        top_missing_includes.append(
                            f"line {lineno}: {include_target} -> {include_path}"
                        )
                    continue

                if line_no_comment.startswith("[") and line_no_comment.endswith("]"):
                    section = line_no_comment.strip("[]").strip().lower()
                    if section == "molecules":
                        molecules_found = True
                    continue
    except OSError as exc:
        raise ValueError(f"Failed to read TOP for preflight check: {top_path}: {exc}") from exc

    preflight["top_includes_checked"] = top_checked_includes
    preflight["top_missing_includes"] = top_missing_includes
    preflight["molecules_section_found"] = molecules_found

    if top_missing_includes:
        preflight["status"] = "failed"
        if ctx.manifest:
            ctx.manifest.set("topology_preflight", preflight)
        missing_summary = "\n  - " + "\n  - ".join(top_missing_includes)
        raise FileNotFoundError(
            "Topology preflight failed: unresolved #include file(s) in system TOP.\n"
            f"TOP: {top_path}\n"
            f"Missing include(s):{missing_summary}\n"
            "Fix include paths under IN/systems/<SYSTEM_ID>/gromacs/top/ (or referenced subdirectories)."
        )

    if not molecules_found:
        preflight["status"] = "failed"
        if ctx.manifest:
            ctx.manifest.set("topology_preflight", preflight)
        raise ValueError(
            f"Topology preflight failed: {top_path} does not contain a [ molecules ] section."
        )

    preflight["status"] = "ok"
    if ctx.manifest:
        ctx.manifest.set("topology_preflight", preflight)

    # Task D: Verify sanitizer outputs exist (contract: gromacs/combined_atomtypes_current.itp)
    allow_unsanitized = getattr(ctx, "allow_unsanitized_grompp", False)
    combined_canonical = gromacs_dir / "combined_atomtypes_current.itp"
    combined_legacy = gromacs_dir / "itp" / "combined_atomtypes_current.itp"
    sanitized_canonical = gromacs_dir / "itp_sanitized_current"
    sanitized_legacy = gromacs_dir / "itp" / "itp_sanitized_current"

    combined_path = combined_canonical if combined_canonical.exists() else None
    combined_legacy_found = False
    if combined_path is None and combined_legacy.exists():
        combined_path = combined_legacy
        combined_legacy_found = True

    sanitized_path = sanitized_canonical if sanitized_canonical.exists() else None
    sanitized_legacy_found = False
    if sanitized_path is None and sanitized_legacy.exists():
        sanitized_path = sanitized_legacy
        sanitized_legacy_found = True

    missing = []
    if combined_path is None:
        missing.append(str(combined_canonical))
    if sanitized_path is None:
        missing.append(str(sanitized_canonical))

    if missing:
        msg = (
            "ITP Sanitizer outputs missing. The Sanitizer must run before grompp.\n"
            f"Expected:\n  - {combined_canonical}\n  - {sanitized_canonical}\n"
            "Run: python run_pipeline.py ... --stage sanitizer"
        )
        if allow_unsanitized:
            print(f"  [WARN] {msg}")
        else:
            raise FileNotFoundError(msg)

    if combined_path:
        if combined_legacy_found:
            print(
                "  [WARN] Using legacy sanitizer path for combined_atomtypes_current.itp. "
                f"Please migrate to {combined_canonical}."
            )
        print(f"  - Sanitizer output verified: {combined_path}")
    if sanitized_path:
        if sanitized_legacy_found:
            print(
                "  [WARN] Using legacy sanitizer path for itp_sanitized_current. "
                f"Please migrate to {sanitized_canonical}."
            )
        print(f"  - Sanitizer output verified: {sanitized_path}")

    # Record sanitizer verification in manifest
    if ctx.manifest:
        ctx.manifest.set("sanitizer_verified", combined_path is not None and sanitized_path is not None)
        ctx.manifest.set("sanitizer_check", {
            "combined_atomtypes_path": _normalize_path_for_record(ctx, combined_path),
            "combined_atomtypes_legacy": combined_legacy_found,
            "sanitized_itp_dir": _normalize_path_for_record(ctx, sanitized_path),
            "sanitized_itp_legacy": sanitized_legacy_found,
            "allow_unsanitized_grompp": allow_unsanitized,
        })

    print(f"  - Selected GRO: {gro_path}")
    print(f"  - Selected TOP: {top_path}")
    print(f"  - Preflight: GRO atom count={preflight['gro_atom_count']}, includes={len(top_checked_includes)}")
    return gro_path, top_path


def find_ndx_file(ctx: "PipelineContext") -> Optional[Path]:
    """Find index file if it exists."""
    ndx_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "ndx")
    
    if ndx_dir.exists():
        index_ndx = ndx_dir / "index.ndx"
        if index_ndx.exists():
            return index_ndx
        ndx_files = sorted(ndx_dir.glob("*.ndx"))
        if len(ndx_files) == 1:
            return ndx_files[0]
        if len(ndx_files) > 1:
            raise ValueError(
                f"Multiple .ndx files found in {ndx_dir}. "
                "Please name the intended file 'index.ndx' or remove extras."
            )
    
    return None


def _tail_text(text: str, max_lines: int = 60, max_chars: int = 2000) -> str:
    """Return a bounded tail excerpt for error surfacing."""
    if not text:
        return ""
    lines = text.splitlines()
    excerpt = "\n".join(lines[-max_lines:])
    if len(excerpt) > max_chars:
        excerpt = excerpt[-max_chars:]
    return excerpt.strip()


def _safe_float(value: Any) -> Optional[float]:
    """Best-effort float conversion."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _mean_std(values: List[float]) -> Tuple[Optional[float], Optional[float]]:
    """Return mean and population std for a list."""
    if not values:
        return None, None
    avg = sum(values) / float(len(values))
    var = sum((x - avg) ** 2 for x in values) / float(len(values))
    return avg, math.sqrt(max(0.0, var))


def _window_stats(series: List[Tuple[float, float]], window_ps: float) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Compute avg/std over the last window_ps in a time series."""
    if not series:
        return None, None, None
    end_t = series[-1][0]
    start_t = max(series[0][0], end_t - window_ps)
    vals = [value for t_ps, value in series if t_ps >= start_t]
    avg, std = _mean_std(vals)
    return avg, std, max(0.0, end_t - start_t)


def _parse_energy_menu(menu_text: str) -> List[Tuple[int, str]]:
    """Parse gmx energy menu lines into (index, label)."""
    entries: List[Tuple[int, str]] = []
    seen: set[int] = set()
    for raw_line in (menu_text or "").splitlines():
        line = raw_line.rstrip()
        match = re.match(r"^\s*(\d+)\s+(.+?)\s*$", line)
        if not match:
            continue
        idx = int(match.group(1))
        if idx in seen:
            continue
        label = match.group(2).strip()
        if not label:
            continue
        entries.append((idx, label))
        seen.add(idx)
    return entries


def _select_energy_index(entries: List[Tuple[int, str]], candidates: List[str]) -> Optional[int]:
    """Find a deterministic term index from menu entries."""
    lowered = [(idx, label.lower()) for idx, label in entries]
    for candidate in candidates:
        target = candidate.lower()
        for idx, label in lowered:
            if label == target:
                return idx
    for candidate in candidates:
        target = candidate.lower()
        for idx, label in lowered:
            if label.startswith(target):
                return idx
    for candidate in candidates:
        target = candidate.lower()
        for idx, label in lowered:
            if target in label:
                return idx
    return None


def _parse_xvg_series(path: Path) -> List[Tuple[float, float]]:
    """Parse 2-column XVG-like data (time,value)."""
    out: List[Tuple[float, float]] = []
    if not path.exists():
        return out
    try:
        for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw_line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            t_ps = _safe_float(parts[0])
            val = _safe_float(parts[1])
            if t_ps is None or val is None:
                continue
            out.append((t_ps, val))
    except OSError:
        return []
    return out


def _run_energy_command(
    ctx: "PipelineContext",
    cmd: List[str],
    stage_dir: Path,
    input_text: str,
    timeout_s: int = 120,
) -> subprocess.CompletedProcess:
    """Run gmx energy non-interactively with stdin selection."""
    cmd_str = " ".join(str(c) for c in cmd)
    if ctx.dry_run:
        return subprocess.CompletedProcess(cmd, 0, "", "")

    env = get_execution_env(ctx)
    result = subprocess.run(
        cmd,
        cwd=str(stage_dir),
        env=env,
        input=input_text,
        capture_output=True,
        text=True,
        timeout=timeout_s,
        check=False,
    )
    ctx.log_command(
        cmd_str,
        "gmx_eq",
        exit_code=result.returncode,
        stdout=result.stdout[-2000:] if result.stdout else None,
        stderr=result.stderr[-2000:] if result.stderr else None,
    )
    return result


def _resolve_expected_density_g_cm3(ctx: "PipelineContext") -> Optional[float]:
    """Best-effort expected density from prior manifest composition metadata."""
    if ctx.manifest is None:
        return None
    manifest_data = ctx.manifest.data
    if not isinstance(manifest_data, dict):
        return None
    composition = manifest_data.get("composition", {})
    if not isinstance(composition, dict):
        return None
    for key in ("target_density_g_cm3", "density_g_cm3", "achieved_density_g_cm3"):
        value = _safe_float(composition.get(key))
        if value is not None and value > 0:
            return value
    return None


def _build_gmx_eq_metrics(ctx: "PipelineContext", npt_dir: Path, window_ps: float) -> Dict[str, Any]:
    """
    Extract thermodynamic averages from NPT energies into stage_metrics schema.
    """
    metrics: Dict[str, Any] = {
        "stage": "gmx_eq",
        "temperature_K": {"avg": None, "std": None},
        "pressure_bar": {"avg": None, "std": None},
        "density_kg_m3": {"avg": None, "std": None},
        "density_g_cm3": {"avg": None, "std": None},
        "volume_nm3": {"avg": None, "std": None},
        "last_window_ps": float(window_ps),
        "notes": [],
    }
    expected_density = _resolve_expected_density_g_cm3(ctx)
    if expected_density is not None:
        metrics["expected_density_g_cm3"] = expected_density

    if ctx.dry_run:
        metrics["notes"].append("dry_run_no_metrics")
        return metrics

    edr_path = npt_dir / "npt.edr"
    if not edr_path.exists():
        metrics["notes"].extend(["failed_to_extract_metrics", "missing_npt.edr"])
        return metrics

    gmx_cmd = resolve_gmx_command(ctx)
    gmx_tokens = _split_command_tokens(gmx_cmd)
    if not gmx_tokens:
        gmx_tokens = ["gmx"]

    menu_cmd = gmx_tokens + ["energy", "-f", str(edr_path), "-xvg", "none"]
    menu_result = _run_energy_command(ctx, menu_cmd, npt_dir, "0\n")
    menu_text = f"{menu_result.stdout or ''}\n{menu_result.stderr or ''}"
    menu_entries = _parse_energy_menu(menu_text)
    if not menu_entries:
        metrics["notes"].extend(["failed_to_extract_metrics", "energy_menu_parse_failed"])
        return metrics

    selected = {
        "temperature_K": _select_energy_index(menu_entries, ["Temperature"]),
        "pressure_bar": _select_energy_index(menu_entries, ["Pressure"]),
        "density_kg_m3": _select_energy_index(menu_entries, ["Density"]),
        "volume_nm3": _select_energy_index(menu_entries, ["Volume"]),
    }
    missing_terms = [key for key, idx in selected.items() if idx is None]
    if missing_terms:
        metrics["notes"].append(f"missing_energy_terms:{','.join(sorted(missing_terms))}")

    actual_windows: List[float] = []
    extracted_any = False
    with tempfile.TemporaryDirectory(prefix=".eq_metrics.", dir=str(npt_dir)) as tmp_dir:
        tmp_root = Path(tmp_dir)
        for metric_key, index in selected.items():
            if index is None:
                continue
            xvg_path = tmp_root / f"{metric_key}.xvg"
            probe_cmd = gmx_tokens + [
                "energy",
                "-f",
                str(edr_path),
                "-o",
                str(xvg_path),
                "-xvg",
                "none",
            ]
            result = _run_energy_command(ctx, probe_cmd, npt_dir, f"{index}\n0\n")
            if result.returncode != 0 or not xvg_path.exists():
                metrics["notes"].append(f"failed_energy_extract:{metric_key}")
                continue
            series = _parse_xvg_series(xvg_path)
            avg, std, actual_window = _window_stats(series, window_ps)
            if avg is None or std is None or actual_window is None:
                metrics["notes"].append(f"empty_series:{metric_key}")
                continue
            metrics[metric_key] = {"avg": avg, "std": std}
            actual_windows.append(actual_window)
            extracted_any = True

    if metrics["density_kg_m3"]["avg"] is not None:
        density_avg = metrics["density_kg_m3"]["avg"] / 1000.0
        density_std = metrics["density_kg_m3"]["std"] / 1000.0
        metrics["density_g_cm3"] = {"avg": density_avg, "std": density_std}

    if actual_windows:
        metrics["last_window_ps"] = min(actual_windows)
    if not extracted_any:
        metrics["notes"].append("failed_to_extract_metrics")

    return metrics


def _write_gmx_eq_metrics(ctx: "PipelineContext", npt_dir: Path) -> Optional[Path]:
    """Write gmx_eq stage_metrics.json after successful NPT."""
    if ctx.dry_run:
        return None

    window_ps_raw = getattr(ctx, "gmx_eq_metrics_window_ps", 20.0)
    window_ps = _safe_float(window_ps_raw)
    if window_ps is None or window_ps <= 0:
        window_ps = 20.0

    metrics = _build_gmx_eq_metrics(ctx, npt_dir, window_ps=window_ps)
    primary_path = npt_dir / "stage_metrics.json"
    _atomic_write_json(primary_path, metrics, sort_keys=True)

    # Stable dispatcher lookup path (requested at stage granularity).
    eq_metrics_dir = ctx.ensure_output_dir("03_gromacs", "gmx_eq")
    alias_path = eq_metrics_dir / "stage_metrics.json"
    _atomic_write_json(alias_path, metrics, sort_keys=True)
    return alias_path


def _gromacs_failure_hints(stderr: str, stdout: str) -> List[str]:
    """Best-effort hints for common grompp/mdrun staging issues."""
    combined = f"{stderr}\n{stdout}".lower()
    hints: List[str] = []
    if "number of coordinates in coordinate file does not match topology" in combined:
        hints.append(
            "Likely cause: GRO/TOP composition mismatch after staging/crosslinking. "
            "Suggested fix: regenerate/stage matching system.gro + system.top and rerun sanitizer."
        )
    if "no such file or directory" in combined and ("#include" in combined or ".itp" in combined):
        hints.append(
            "Likely cause: missing topology include. Suggested fix: verify #include targets under "
            "IN/systems/<SYSTEM_ID>/gromacs/top and rerun preflight."
        )
    if "atom index" in combined or ("group" in combined and "not found" in combined):
        hints.append(
            "Likely cause: MDP/index group mismatch. Suggested fix: check tc-grps/comm-grps/energygrps "
            "against index.ndx group names."
        )
    return hints


def run_gromacs_command(
    ctx: "PipelineContext",
    cmd: List[str],
    stage_dir: Path,
    stage_name: str,
    check: bool = True,
    stream_logs: bool = False,
    log_prefix: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """
    Run a GROMACS command with proper environment and logging.
    
    Args:
        ctx: Pipeline context
        cmd: Command to run
        stage_dir: Stage directory (for crash context)
        stage_name: Stage name for logging
        check: If True, raise on non-zero exit
        
    Returns:
        CompletedProcess result
    """
    cmd_str = " ".join(str(c) for c in cmd)
    print(f"  $ {cmd_str}")
    
    if ctx.dry_run:
        print("  [DRY-RUN] Command not executed")
        ctx.log_command(cmd_str, stage_name, exit_code=0)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    
    env = get_execution_env(ctx)
    
    if not stream_logs:
        result = subprocess.run(
            cmd,
            cwd=str(stage_dir),
            env=env,
            capture_output=True,
            text=True,
        )

        # Log command
        ctx.log_command(
            cmd_str, stage_name,
            exit_code=result.returncode,
            stdout=result.stdout[-2000:] if result.stdout else None,
            stderr=result.stderr[-2000:] if result.stderr else None,
        )

        if check and result.returncode != 0:
            # Handle crash context
            handle_stage_failure(stage_dir, cmd, result.returncode)
            stderr_tail = _tail_text(result.stderr or "")
            stdout_tail = _tail_text(result.stdout or "")
            hints = _gromacs_failure_hints(result.stderr or "", result.stdout or "")
            msg_lines = [
                f"Command failed with exit code {result.returncode}.",
                f"Command: {cmd_str}",
            ]
            if stderr_tail:
                msg_lines.append("stderr tail:")
                msg_lines.append(stderr_tail)
            if stdout_tail:
                msg_lines.append("stdout tail:")
                msg_lines.append(stdout_tail)
            if hints:
                msg_lines.append("Likely cause / Suggested fix:")
                msg_lines.extend(f"- {hint}" for hint in hints)
            raise RuntimeError("\n".join(msg_lines))

        return result

    # Stream stdout/stderr to logs for long-running commands (mdrun)
    logs_dir = stage_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    prefix = log_prefix or "mdrun"
    stdout_path = logs_dir / f"{prefix}.stdout.log"
    stderr_path = logs_dir / f"{prefix}.stderr.log"

    with open(stdout_path, "w") as stdout_f, open(stderr_path, "w") as stderr_f:
        popen_kwargs = {
            "cwd": str(stage_dir),
            "env": env,
            "stdout": stdout_f,
            "stderr": stderr_f,
            "text": True,
        }
        if os.name == "posix":
            popen_kwargs["start_new_session"] = True
        try:
            proc = subprocess.Popen(cmd, **popen_kwargs)
        except TypeError:
            popen_kwargs.pop("start_new_session", None)
            proc = subprocess.Popen(cmd, **popen_kwargs)

        original_handlers = {}

        def _forward_signal(signum, _frame):
            if proc.poll() is not None:
                return
            try:
                if hasattr(os, "killpg"):
                    os.killpg(proc.pid, signum)
                else:
                    proc.send_signal(signum)
            except Exception:
                try:
                    proc.send_signal(signum)
                except Exception:
                    pass

        for sig in (signal.SIGTERM, signal.SIGINT):
            try:
                original_handlers[sig] = signal.getsignal(sig)
                signal.signal(sig, _forward_signal)
            except Exception:
                continue

        try:
            returncode = proc.wait()
        finally:
            for sig, handler in original_handlers.items():
                try:
                    signal.signal(sig, handler)
                except Exception:
                    continue

    # Log command with log paths
    ctx.log_command(
        cmd_str,
        stage_name,
        exit_code=returncode,
        stdout=str(stdout_path),
        stderr=str(stderr_path),
    )

    if check and returncode != 0:
        handle_stage_failure(
            stage_dir,
            cmd,
            returncode,
            log_paths=[stdout_path, stderr_path],
        )
        raise RuntimeError(f"Command failed with exit code {returncode}")

    return subprocess.CompletedProcess(cmd, returncode, "", "")


def cleanup_large_files(stage_dir: Path, ctx: "PipelineContext", stage_tag: str) -> List[str]:
    """
    Remove large files (.trr) if cleanup policy permits.
    
    TRR cleanup policy is stage-aware and gated on:
    - ctx.cleanup_success must be True
    - ctx.cleanup_delete_trr must be True (explicit opt-in)

    Policy:
    - em/nvt/npt: delete TRR when cleanup_delete_trr=True (equilibration TRR is disposable)
    - md: keep TRR if analysis_requires_velocities=True, else delete when cleanup_delete_trr=True
    
    Returns:
        List of deleted file names (empty if nothing deleted)
    """
    if ctx.dry_run:
        return []
    if not ctx.cleanup_success:
        return []

    # Gate TRR deletion: explicit opt-in required
    if not getattr(ctx, 'cleanup_delete_trr', False):
        print(f"  [CLEANUP] TRR files retained (use --cleanup-delete-trr to remove)")
        return []

    analysis_needs_vel = bool(getattr(ctx, "analysis_requires_velocities", False))
    if stage_tag == "md" and analysis_needs_vel:
        print("  [CLEANUP] Keeping production MD TRR (analysis requires velocities).")
        if ctx.manifest:
            existing = ctx.manifest.get("trr_cleanup", {"records": []})
            records = list(existing.get("records", [])) if isinstance(existing, dict) else []
            records.append({
                "stage_tag": stage_tag,
                "stage_dir": str(stage_dir),
                "deleted_files": [],
                "cleanup_delete_trr": True,
                "analysis_requires_velocities": analysis_needs_vel,
                "decision": "kept_for_analysis",
            })
            ctx.manifest.set("trr_cleanup", {"records": records})
        return []

    print(f"  [CLEANUP] TRR deletion enabled for stage '{stage_tag}' (explicit opt-in).")
    deleted = []
    failed_deletes: List[Dict[str, str]] = []
    for trr_file in stage_dir.glob("*.trr"):
        try:
            print(f"  [CLEANUP] Removing: {trr_file.name}")
            trr_file.unlink()
            deleted.append(trr_file.name)
        except OSError as e:
            err_msg = f"{type(e).__name__}: {e}"
            print(f"  [CLEANUP][WARN] Failed to remove {trr_file.name}: {err_msg}")
            failed_deletes.append({
                "file": trr_file.name,
                "error": err_msg,
            })

    if ctx.manifest:
        existing = ctx.manifest.get("trr_cleanup", {"records": []})
        records = list(existing.get("records", [])) if isinstance(existing, dict) else []
        decision = "failed_to_delete" if failed_deletes else ("deleted" if deleted else "nothing_to_delete")
        record: Dict[str, Any] = {
            "stage_tag": stage_tag,
            "stage_dir": str(stage_dir),
            "deleted_files": deleted,
            "cleanup_delete_trr": True,
            "analysis_requires_velocities": analysis_needs_vel,
            "decision": decision,
        }
        if failed_deletes:
            record["failed_deletes"] = failed_deletes
        records.append(record)
        ctx.manifest.set("trr_cleanup", {"records": records})
    
    return deleted


class GmxEmStage(BaseStage):
    """
    Stage 3a: GROMACS Energy Minimization
    
    Reads:
    - IN/systems/<SYSTEM_ID>/gromacs/mdp/em.mdp - EM parameters
    - IN/systems/<SYSTEM_ID>/gromacs/gro/system.gro - staged structure
    - IN/systems/<SYSTEM_ID>/gromacs/top/system.top - staged topology
    
    Writes:
    - OUT_GMX/<RUN_ID>/03_gromacs/em/ - energy minimized system
    """
    
    @property
    def name(self) -> str:
        return "gmx_em"
    
    @property
    def output_subdir(self) -> str:
        return "03_gromacs/em"
    
    def run(self, ctx: "PipelineContext") -> bool:
        """Execute GROMACS energy minimization."""
        print(f"  Running GROMACS Energy Minimization...")
        
        # Output directory
        output_dir = self.get_output_dir(ctx)
        print(f"  - Output dir: {output_dir}")
        done_sentinel = output_dir / "done.ok"

        # Skip if already complete
        if done_sentinel.exists() and ctx.resume and not ctx.force:
            print(f"  [SKIP] EM already complete")
            _record_stage_manifest(ctx, "gmx_em", {
                "decision": "skip",
                "reason": "done.ok exists and resume=true",
                "stage_dir": str(output_dir),
            })
            return True

        # Force rerun: archive prior outputs
        if ctx.force:
            _archive_stage_outputs(output_dir, "em", ctx)

        # Verify inputs
        try:
            gro_path, top_path = verify_staged_inputs(ctx)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [ERROR] {e}")
            return False

        # Resolve ndx before MDP patching so comm-grps policy checks can validate group presence.
        try:
            ndx_path = find_ndx_file(ctx)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        ctx.active_ndx_path = str(ndx_path) if ndx_path else None
        ctx.input_gro_has_velocities = _gro_has_velocities(gro_path)
        
        # Prepare MDP
        try:
            stage_overrides, extra_overrides, effective_overrides = resolve_mdp_overrides(ctx, "em")
            mdp_path, orig_vals, patched_vals = prepare_stage_mdp(
                ctx,
                "em",
                output_dir,
                effective_overrides=effective_overrides,
                checkpoint_available=None,
            )
            print(f"  - MDP: {mdp_path.name}")
            if patched_vals:
                print(f"    Patched: {patched_vals}")
        except MDPPatchError as e:
            print(f"  [ERROR] {e}")
            return False
        try:
            _assert_ndx_for_mdp(mdp_path, ndx_path)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        
        # POSRES reference handling
        try:
            posres_ref, posres_info = _resolve_posres_info(ctx, "em", gro_path, mdp_path)
        except (ValueError, FileNotFoundError) as e:
            print(f"  [ERROR] {e}")
            return False

        # Run grompp
        tpr_path = output_dir / "em.tpr"
        grompp_cmd = build_grompp_command(
            ctx=ctx,
            mdp_path=str(mdp_path),
            structure_path=str(gro_path),
            topology_path=str(top_path),
            output_path=str(tpr_path),
            index_path=str(ndx_path) if ndx_path else None,
            reference_structure_path=str(posres_ref) if posres_ref else None,
            maxwarn=ctx.grompp_maxwarn,
        )
        
        mdrun_cmd_base = build_mdrun_command(ctx, "em")
        stage_state = _build_stage_state(
            ctx,
            "em",
            mdp_path,
            gro_path,
            top_path,
            ndx_path,
            grompp_cmd,
            mdrun_cmd_base,
            ctx.grompp_maxwarn,
            posres_info=posres_info,
            stage_overrides=stage_overrides,
            extra_overrides=extra_overrides,
            effective_overrides=effective_overrides,
        )

        try:
            run_gromacs_command(ctx, grompp_cmd, output_dir, self.name)
        except RuntimeError:
            return False
        
        # Run mdrun
        mdrun_cmd = mdrun_cmd_base
        
        try:
            run_gromacs_command(ctx, mdrun_cmd, output_dir, self.name, stream_logs=True, log_prefix="mdrun")
        except RuntimeError:
            return False

        if not ctx.dry_run:
            done_sentinel.write_text("EM complete\n")
            _write_stage_state(output_dir / "stage_state.json", stage_state, ctx)

        # Record resources in manifest
        if ctx.manifest:
            ctx.manifest.set("resources", get_resource_summary(ctx))

        _record_stage_manifest(ctx, "gmx_em", {
            "decision": "run",
            "reason": "completed",
            "stage_dir": str(output_dir),
            "mdp_path": str(mdp_path),
            "mdp_sha256": stage_state["mdp_sha256"],
            "mdp_original_values": orig_vals,
            "mdp_patched_values": patched_vals,
            "inputs": {
                "gro": {"path": str(gro_path), "sha256": stage_state["gro_sha256"]},
                "top": {"path": str(top_path), "sha256": stage_state["top_sha256"]},
                "ndx": {
                    "path": str(ndx_path) if ndx_path else None,
                    "sha256": stage_state["ndx_sha256"],
                },
            },
            "grompp_maxwarn": ctx.grompp_maxwarn,
            "grompp_command": stage_state["grompp_command"],
            "mdrun_command": " ".join(mdrun_cmd),
            "stage_state_path": str(output_dir / "stage_state.json"),
            "posres": posres_info,
        })

        repro = _build_repro_json(
            ctx=ctx,
            stage="em",
            stage_dir=output_dir,
            mdp_path=mdp_path,
            mdp_orig=orig_vals,
            mdp_patched=patched_vals,
            gro_path=gro_path,
            top_path=top_path,
            ndx_path=ndx_path,
            cpt_path=None,
            grompp_cmd=stage_state["grompp_command"],
            mdrun_cmd=" ".join(mdrun_cmd),
            resume_allowed=False,
        )

        # Cleanup BEFORE writing repro.json so it reflects final on-disk state
        deleted_files = cleanup_large_files(output_dir, ctx, stage_tag="em")
        if deleted_files:
            repro["cleanup"] = {
                "trr_deleted": deleted_files,
                "enabled": True,
            }

        _write_repro_json(ctx, output_dir, repro)
        _write_run_report(ctx)
        
        print(f"  [OK] Energy minimization complete")
        return True


class GmxEqStage(BaseStage):
    """
    Stage 3b: GROMACS Equilibration (NVT + NPT)
    
    Two substages with separate sentinels:
    - NVT: OUT_GMX/<RUN_ID>/03_gromacs/nvt/done.ok
    - NPT: OUT_GMX/<RUN_ID>/03_gromacs/npt/done.ok
    """
    
    @property
    def name(self) -> str:
        return "gmx_eq"
    
    @property
    def output_subdir(self) -> str:
        return "03_gromacs/nvt"  # Primary for sentinel
    
    def run(self, ctx: "PipelineContext") -> bool:
        """Execute GROMACS equilibration (NVT + NPT)."""
        print(f"  Running GROMACS Equilibration...")
        
        # Run NVT
        nvt_success = self._run_nvt(ctx)
        if not nvt_success:
            return False
        
        # Run NPT
        npt_success = self._run_npt(ctx)
        if not npt_success:
            return False
        
        print(f"  [OK] Equilibration complete (NVT + NPT)")
        return True
    
    def _run_nvt(self, ctx: "PipelineContext") -> bool:
        """Run NVT equilibration substage."""
        nvt_dir = ctx.ensure_output_dir("03_gromacs", "nvt")
        nvt_sentinel = nvt_dir / "done.ok"
        
        # Check if NVT already complete
        if nvt_sentinel.exists() and ctx.resume and not ctx.force:
            print(f"  [SKIP] NVT already complete")
            _record_stage_manifest(ctx, "gmx_eq_nvt", {
                "decision": "skip",
                "reason": "done.ok exists and resume=true",
                "stage_dir": str(nvt_dir),
            })
            return True

        if ctx.force:
            _archive_stage_outputs(nvt_dir, "nvt", ctx)
        
        print(f"  Running NVT equilibration...")
        
        # Get input: EM output
        em_dir = ctx.get_output_path("03_gromacs", "em")
        em_gro = em_dir / "em.gro"
        
        if not em_gro.exists():
            print(f"  [ERROR] EM output not found: {em_gro}")
            return False
        
        # Get topology
        try:
            _, top_path = verify_staged_inputs(ctx)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [ERROR] {e}")
            return False

        try:
            ndx_path = find_ndx_file(ctx)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        ctx.active_ndx_path = str(ndx_path) if ndx_path else None
        ctx.input_gro_has_velocities = _gro_has_velocities(em_gro)
        
        # Prepare MDP
        try:
            stage_overrides, extra_overrides, effective_overrides = resolve_mdp_overrides(ctx, "nvt")
            mdp_path, orig_vals, patched_vals = prepare_stage_mdp(
                ctx,
                "nvt",
                nvt_dir,
                effective_overrides=effective_overrides,
                checkpoint_available=None,
            )
            print(f"    MDP: {mdp_path.name}")
            if patched_vals:
                print(f"    Patched: {patched_vals}")
        except MDPPatchError as e:
            print(f"  [ERROR] {e}")
            return False
        try:
            _assert_ndx_for_mdp(mdp_path, ndx_path)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        thermostat_audit = _collect_thermostat_audit(ctx, "gmx_eq_nvt", mdp_path, ndx_path)
        
        # POSRES reference handling
        try:
            posres_ref, posres_info = _resolve_posres_info(ctx, "nvt", em_gro, mdp_path)
        except (ValueError, FileNotFoundError) as e:
            print(f"  [ERROR] {e}")
            return False

        # Run grompp
        tpr_path = nvt_dir / "nvt.tpr"
        grompp_cmd = build_grompp_command(
            ctx=ctx,
            mdp_path=str(mdp_path),
            structure_path=str(em_gro),
            topology_path=str(top_path),
            output_path=str(tpr_path),
            index_path=str(ndx_path) if ndx_path else None,
            reference_structure_path=str(posres_ref) if posres_ref else None,
            maxwarn=ctx.grompp_maxwarn,
        )

        mdrun_cmd_base = build_mdrun_command(ctx, "nvt")
        stage_state = _build_stage_state(
            ctx,
            "nvt",
            mdp_path,
            em_gro,
            top_path,
            ndx_path,
            grompp_cmd,
            mdrun_cmd_base,
            ctx.grompp_maxwarn,
            posres_info=posres_info,
            stage_overrides=stage_overrides,
            extra_overrides=extra_overrides,
            effective_overrides=effective_overrides,
        )
        
        try:
            run_gromacs_command(ctx, grompp_cmd, nvt_dir, self.name)
        except RuntimeError:
            return False
        
        # Run mdrun
        mdrun_cmd = mdrun_cmd_base

        try:
            run_gromacs_command(ctx, mdrun_cmd, nvt_dir, self.name, stream_logs=True, log_prefix="mdrun")
        except RuntimeError:
            return False
        
        # Mark NVT complete
        if not ctx.dry_run:
            nvt_sentinel.write_text("NVT equilibration complete\n")
            _write_stage_state(nvt_dir / "stage_state.json", stage_state, ctx)

        _record_stage_manifest(ctx, "gmx_eq_nvt", {
            "decision": "run",
            "reason": "completed",
            "stage_dir": str(nvt_dir),
            "mdp_path": str(mdp_path),
            "mdp_sha256": stage_state["mdp_sha256"],
            "mdp_original_values": orig_vals,
            "mdp_patched_values": patched_vals,
            "inputs": {
                "gro": {"path": str(em_gro), "sha256": stage_state["gro_sha256"]},
                "top": {"path": str(top_path), "sha256": stage_state["top_sha256"]},
                "ndx": {
                    "path": str(ndx_path) if ndx_path else None,
                    "sha256": stage_state["ndx_sha256"],
                },
            },
            "grompp_maxwarn": ctx.grompp_maxwarn,
            "grompp_command": stage_state["grompp_command"],
            "mdrun_command": " ".join(mdrun_cmd),
            "stage_state_path": str(nvt_dir / "stage_state.json"),
            "posres": posres_info,
            "thermostat": thermostat_audit,
        })
        _record_thermostat_manifest(ctx, "gmx_eq_nvt", thermostat_audit)

        repro = _build_repro_json(
            ctx=ctx,
            stage="nvt",
            stage_dir=nvt_dir,
            mdp_path=mdp_path,
            mdp_orig=orig_vals,
            mdp_patched=patched_vals,
            gro_path=em_gro,
            top_path=top_path,
            ndx_path=ndx_path,
            cpt_path=None,
            grompp_cmd=stage_state["grompp_command"],
            mdrun_cmd=" ".join(mdrun_cmd),
            resume_allowed=False,
        )

        # Cleanup BEFORE writing repro.json so it reflects final on-disk state
        deleted_files = cleanup_large_files(nvt_dir, ctx, stage_tag="nvt")
        if deleted_files:
            repro["cleanup"] = {
                "trr_deleted": deleted_files,
                "enabled": True,
            }

        _write_repro_json(ctx, nvt_dir, repro)
        _write_run_report(ctx)
        print(f"    [OK] NVT complete")
        return True
    
    def _run_npt(self, ctx: "PipelineContext") -> bool:
        """Run NPT equilibration substage."""
        npt_dir = ctx.ensure_output_dir("03_gromacs", "npt")
        npt_sentinel = npt_dir / "done.ok"
        
        # Check if NPT already complete
        if npt_sentinel.exists() and ctx.resume and not ctx.force:
            print(f"  [SKIP] NPT already complete")
            _record_stage_manifest(ctx, "gmx_eq_npt", {
                "decision": "skip",
                "reason": "done.ok exists and resume=true",
                "stage_dir": str(npt_dir),
            })
            return True

        if ctx.force:
            _archive_stage_outputs(npt_dir, "npt", ctx)
        
        print(f"  Running NPT equilibration...")
        
        # Get input: NVT output
        nvt_dir = ctx.get_output_path("03_gromacs", "nvt")
        nvt_gro = nvt_dir / "nvt.gro"
        nvt_cpt = nvt_dir / "nvt.cpt"
        
        if not nvt_gro.exists():
            print(f"  [ERROR] NVT output not found: {nvt_gro}")
            return False
        
        # Use shared decision function for checkpoint/velocity handling
        template_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "mdp")
        template_path = template_dir / "npt.mdp"
        
        cpt_available, extra_overrides, velocity_mode, policy_warnings, policy_metadata = decide_prev_state_policy(
            stage_name="npt",
            prev_cpt_path=nvt_cpt,
            strict_mode=_is_restart_policy_strict(ctx),
            allow_velocity_reset=ctx.allow_velocity_reset,
            allow_gro_velocity_restart=getattr(ctx, "allow_gro_velocity_restart", False),
            template_mdp_path=template_path,
            ctx=ctx,
            prev_gro_path=nvt_gro,
        )
        
        # Handle error case
        if velocity_mode == "error":
            for msg in policy_warnings:
                print(f"  [ERROR] {msg}")
            return False
        
        # Print warnings
        for w in policy_warnings:
            print(f"  [WARN] {w}")
        
        # Build velocity reset info for auditing
        velocity_reset_info = {
            "velocity_mode": velocity_mode,
            "warnings": policy_warnings,
            "policy_metadata": policy_metadata,
        }
        if extra_overrides:
            velocity_reset_info["overrides"] = extra_overrides

        # Get topology
        try:
            _, top_path = verify_staged_inputs(ctx)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [ERROR] {e}")
            return False

        try:
            ndx_path = find_ndx_file(ctx)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        ctx.active_ndx_path = str(ndx_path) if ndx_path else None
        ctx.input_gro_has_velocities = _gro_has_velocities(nvt_gro)
        
        # Prepare MDP
        try:
            stage_overrides, extra_overrides_resolved, effective_overrides = resolve_mdp_overrides(
                ctx, "npt", extra_overrides=extra_overrides
            )
            mdp_path, orig_vals, patched_vals = prepare_stage_mdp(
                ctx,
                "npt",
                npt_dir,
                effective_overrides=effective_overrides,
                checkpoint_available=cpt_available,
            )
            print(f"    MDP: {mdp_path.name}")
            if patched_vals:
                print(f"    Patched: {patched_vals}")
            if velocity_mode != "checkpoint":
                print(f"    Velocity mode: {velocity_mode}")
                if extra_overrides:
                    print(f"    Velocity overrides: {extra_overrides}")
        except MDPPatchError as e:
            print(f"  [ERROR] {e}")
            return False
        # Validate velocity reset semantics against final patched NPT MDP.
        if velocity_mode == "velocity_reset":
            if not extra_overrides or "gen_temp" not in extra_overrides:
                print("  [ERROR] Stage 'npt': velocity reset requires resolved gen_temp override.")
                return False
            try:
                expected_gen_temp = extra_overrides["gen_temp"]
                _assert_velocity_reset_mdp(mdp_path, "npt", expected_gen_temp)
            except MDPPatchError as e:
                print(f"  [ERROR] {e}")
                return False

        try:
            _assert_ndx_for_mdp(mdp_path, ndx_path)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        thermostat_audit = _collect_thermostat_audit(ctx, "gmx_eq_npt", mdp_path, ndx_path)
        
        # POSRES reference handling
        try:
            posres_ref, posres_info = _resolve_posres_info(ctx, "npt", nvt_gro, mdp_path)
        except (ValueError, FileNotFoundError) as e:
            print(f"  [ERROR] {e}")
            return False

        # Run grompp with checkpoint for velocities
        tpr_path = npt_dir / "npt.tpr"
        grompp_cmd = build_grompp_command(
            ctx=ctx,
            mdp_path=str(mdp_path),
            structure_path=str(nvt_gro),
            topology_path=str(top_path),
            output_path=str(tpr_path),
            checkpoint_path=str(nvt_cpt) if nvt_cpt.exists() else None,
            index_path=str(ndx_path) if ndx_path else None,
            reference_structure_path=str(posres_ref) if posres_ref else None,
            maxwarn=ctx.grompp_maxwarn,
        )

        mdrun_cmd_base = build_mdrun_command(ctx, "npt")
        stage_state = _build_stage_state(
            ctx,
            "npt",
            mdp_path,
            nvt_gro,
            top_path,
            ndx_path,
            grompp_cmd,
            mdrun_cmd_base,
            ctx.grompp_maxwarn,
            posres_info=posres_info,
            stage_overrides=stage_overrides,
            extra_overrides=extra_overrides_resolved,
            effective_overrides=effective_overrides,
            velocity_mode=velocity_mode,
        )
        
        try:
            run_gromacs_command(ctx, grompp_cmd, npt_dir, self.name)
        except RuntimeError:
            return False
        
        # Run mdrun
        mdrun_cmd = mdrun_cmd_base
        
        try:
            run_gromacs_command(ctx, mdrun_cmd, npt_dir, self.name, stream_logs=True, log_prefix="mdrun")
        except RuntimeError:
            return False

        stage_metrics_path: Optional[Path] = None
        if not ctx.dry_run:
            try:
                stage_metrics_path = _write_gmx_eq_metrics(ctx, npt_dir)
                if stage_metrics_path is not None:
                    print(f"    Metrics: {stage_metrics_path}")
            except Exception as exc:
                print(f"    [WARN] Failed to write gmx_eq stage metrics: {type(exc).__name__}: {exc}")
        
        # Mark NPT complete
        if not ctx.dry_run:
            npt_sentinel.write_text("NPT equilibration complete\n")
            _write_stage_state(npt_dir / "stage_state.json", stage_state, ctx)
        
        # Record resources
        if ctx.manifest:
            ctx.manifest.set("resources", get_resource_summary(ctx))

        _record_stage_manifest(ctx, "gmx_eq_npt", {
            "decision": "run",
            "reason": "completed",
            "stage_dir": str(npt_dir),
            "mdp_path": str(mdp_path),
            "mdp_sha256": stage_state["mdp_sha256"],
            "mdp_original_values": orig_vals,
            "mdp_patched_values": patched_vals,
            "inputs": {
                "gro": {"path": str(nvt_gro), "sha256": stage_state["gro_sha256"]},
                "top": {"path": str(top_path), "sha256": stage_state["top_sha256"]},
                "ndx": {
                    "path": str(ndx_path) if ndx_path else None,
                    "sha256": stage_state["ndx_sha256"],
                },
                "checkpoint": {
                    "path": str(nvt_cpt) if nvt_cpt.exists() else None,
                    "required": True,
                    "allow_velocity_reset": ctx.allow_velocity_reset,
                },
            },
            "grompp_maxwarn": ctx.grompp_maxwarn,
            "grompp_command": stage_state["grompp_command"],
            "mdrun_command": " ".join(mdrun_cmd),
            "stage_state_path": str(npt_dir / "stage_state.json"),
            "stage_metrics_path": str(stage_metrics_path) if stage_metrics_path else None,
            "posres": posres_info,
            "velocity_reset": velocity_reset_info,
            "thermostat": thermostat_audit,
        })
        _record_thermostat_manifest(ctx, "gmx_eq_npt", thermostat_audit)

        repro = _build_repro_json(
            ctx=ctx,
            stage="npt",
            stage_dir=npt_dir,
            mdp_path=mdp_path,
            mdp_orig=orig_vals,
            mdp_patched=patched_vals,
            gro_path=nvt_gro,
            top_path=top_path,
            ndx_path=ndx_path,
            cpt_path=nvt_cpt if cpt_available else None,
            grompp_cmd=stage_state["grompp_command"],
            mdrun_cmd=" ".join(mdrun_cmd),
            resume_allowed=False,
        )
        
        # Cleanup BEFORE writing repro.json so it reflects final on-disk state
        deleted_files = cleanup_large_files(npt_dir, ctx, stage_tag="npt")
        if deleted_files:
            repro["cleanup"] = {
                "trr_deleted": deleted_files,
                "enabled": True,
            }
        
        _write_repro_json(ctx, npt_dir, repro)
        _write_run_report(ctx)

        print(f"    [OK] NPT complete")
        return True


class GmxProdStage(BaseStage):
    """
    Stage 3c: GROMACS Production MD
    
    Supports checkpoint resume:
    - If md/*.cpt exists: mdrun -cpi md.cpt -append
    - Else: grompp -c npt.gro -t npt.cpt
    """
    
    @property
    def name(self) -> str:
        return "gmx_prod"
    
    @property
    def output_subdir(self) -> str:
        return "03_gromacs/md"
    
    def run(self, ctx: "PipelineContext") -> bool:
        """Execute GROMACS production MD."""
        print(f"  Running GROMACS Production MD...")
        
        output_dir = self.get_output_dir(ctx)
        done_sentinel = output_dir / "done.ok"

        if done_sentinel.exists() and ctx.resume and not ctx.force:
            print(f"  [SKIP] MD already complete")
            _record_stage_manifest(ctx, "gmx_prod", {
                "decision": "skip",
                "reason": "done.ok exists and resume=true",
                "stage_dir": str(output_dir),
            })
            return True

        if ctx.force:
            _archive_stage_outputs(output_dir, "md", ctx)

        md_cpt = output_dir / "md.cpt"
        md_tpr = output_dir / "md.tpr"
        stage_state_path = output_dir / "stage_state.json"

        # Inputs
        try:
            _, top_path = verify_staged_inputs(ctx)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [ERROR] {e}")
            return False

        npt_dir = ctx.get_output_path("03_gromacs", "npt")
        npt_gro = npt_dir / "npt.gro"
        npt_cpt = npt_dir / "npt.cpt"

        if not npt_gro.exists():
            print(f"  [ERROR] NPT output not found: {npt_gro}")
            print("  Run gmx_eq stage first.")
            return False

        velocity_reset_info = None
        extra_overrides = None
        mdp_path = output_dir / "md.mdp"

        repair_performed = False
        repair_details: Optional[Dict[str, Any]] = None
        if md_cpt.exists() and not md_tpr.exists():
            print("  [WARN] Partial MD outputs detected: md.cpt exists but md.tpr is missing.")
            if not ctx.resume:
                print("  [ERROR] Resume repair requires resume mode.")
                print("  Use default resume behavior or remove partial outputs before rerunning.")
                return False
            if not getattr(ctx, "allow_resume_repair", False):
                print("  [ERROR] Resume repair is disabled.")
                print("  md.cpt exists but md.tpr is missing; do not use --force for this case.")
                print("  Use --allow-resume-repair to rebuild md.tpr and continue, or restore md.tpr manually.")
                return False
            if not mdp_path.exists():
                print("  [ERROR] Cannot repair resume: md.mdp is missing.")
                print("  Restore md.mdp from this run, or use --force for a clean fresh production start.")
                return False

            old_state_for_repair: Optional[Dict[str, Any]] = None
            if stage_state_path.exists():
                try:
                    loaded_state = json.loads(stage_state_path.read_text())
                    if isinstance(loaded_state, dict):
                        old_state_for_repair = loaded_state
                except Exception as e:
                    print(f"  [ERROR] Failed to read stage_state.json for repair: {e}")
                    return False

            repair_gro = npt_gro
            repair_top = top_path
            repair_ndx: Optional[Path] = None
            if old_state_for_repair:
                old_gro = _resolve_recorded_path(ctx, old_state_for_repair.get("gro_path"))
                if old_gro and old_gro.exists():
                    repair_gro = old_gro
                old_top = _resolve_recorded_path(ctx, old_state_for_repair.get("top_path"))
                if old_top and old_top.exists():
                    repair_top = old_top
                old_ndx = _resolve_recorded_path(ctx, old_state_for_repair.get("ndx_path"))
                if old_ndx and old_ndx.exists():
                    repair_ndx = old_ndx
            if repair_ndx is None:
                try:
                    repair_ndx = find_ndx_file(ctx)
                except ValueError as e:
                    print(f"  [ERROR] {e}")
                    return False

            repair_posres_ref: Optional[Path] = None
            if old_state_for_repair and isinstance(old_state_for_repair.get("posres"), dict):
                old_posres = old_state_for_repair["posres"]
                old_ref = _resolve_recorded_path(ctx, old_posres.get("reference_path"))
                if old_ref and old_ref.exists():
                    repair_posres_ref = old_ref
                elif old_posres.get("active"):
                    print(
                        "  [WARN] stage_state POSRES reference is missing on disk; "
                        "resolving reference from current CLI policy."
                    )
            if repair_posres_ref is None:
                try:
                    repair_posres_ref, _ = _resolve_posres_info(ctx, "md", repair_gro, mdp_path)
                except (ValueError, FileNotFoundError) as e:
                    print(f"  [ERROR] Cannot repair md.tpr: {e}")
                    return False

            repair_maxwarn = ctx.grompp_maxwarn
            if old_state_for_repair and isinstance(old_state_for_repair.get("grompp_maxwarn"), int):
                if int(old_state_for_repair["grompp_maxwarn"]) >= 0:
                    repair_maxwarn = int(old_state_for_repair["grompp_maxwarn"])
            repair_grompp_cmd = build_grompp_command(
                ctx=ctx,
                mdp_path=str(mdp_path),
                structure_path=str(repair_gro),
                topology_path=str(repair_top),
                output_path=str(md_tpr),
                checkpoint_path=str(md_cpt),
                index_path=str(repair_ndx) if repair_ndx else None,
                reference_structure_path=str(repair_posres_ref) if repair_posres_ref else None,
                maxwarn=repair_maxwarn,
            )
            print("  [REPAIR] Rebuilding md.tpr from existing md.cpt for continuation.")
            try:
                run_gromacs_command(ctx, repair_grompp_cmd, output_dir, self.name)
            except RuntimeError:
                return False
            if not ctx.dry_run and not md_tpr.exists():
                print("  [ERROR] Resume repair failed: md.tpr was not created.")
                return False
            repair_performed = True
            repair_details = {
                "enabled": True,
                "tpr_rebuilt": True,
                "mdp_path": str(mdp_path),
                "gro_path": str(repair_gro),
                "top_path": str(repair_top),
                "ndx_path": str(repair_ndx) if repair_ndx else None,
                "used_checkpoint": str(md_cpt),
                "grompp_maxwarn": repair_maxwarn,
            }
        elif md_tpr.exists() and not md_cpt.exists():
            print("  [ERROR] Partial MD outputs detected: md.tpr exists but md.cpt is missing.")
            print("  Continuation cannot be reconstructed automatically without md.cpt.")
            print("  Options:")
            print("    1. Restore md.cpt from backup and rerun with resume.")
            print("    2. Use --force to archive current md outputs and start fresh.")
            print("    3. If intentionally restarting without checkpoint, remove md.tpr and choose")
            print("       --allow-velocity-reset or --allow-gro-velocity-restart explicitly.")
            return False

        resume_available = md_cpt.exists() and md_tpr.exists()

        # Resume path: do NOT rewrite md.mdp or require npt.cpt
        if resume_available:
            if not ctx.resume:
                print("  [ERROR] MD outputs exist but --resume is false.")
                print("  Use --resume to continue or --force to restart.")
                return False
            if not stage_state_path.exists():
                print("  [ERROR] Missing stage_state.json for resume validation.")
                print("  Restore stage_state.json from this run, or use --force to restart cleanly.")
                return False
            if not mdp_path.exists():
                print("  [ERROR] Missing md.mdp for resume validation.")
                print("  Restore md.mdp from this run, or use --force to restart cleanly.")
                return False
            try:
                old_state = json.loads(stage_state_path.read_text())
            except Exception as e:
                print(f"  [ERROR] Failed to read stage_state.json: {e}")
                return False
            old_payload = old_state.get("fingerprint_payload")
            if not isinstance(old_payload, dict):
                print("  [ERROR] stage_state.json missing fingerprint_payload.")
                print("  Restore a valid stage_state.json, or use --force to restart cleanly.")
                return False

            resume_fallback_warnings: List[str] = []

            def _resolve_resume_input(
                payload_key: str,
                state_key: str,
                default_path: Optional[Path],
            ) -> Optional[Path]:
                raw_value = old_payload.get(payload_key)
                if raw_value is None and payload_key not in old_payload:
                    raw_value = old_state.get(state_key)
                candidate = _resolve_recorded_path(ctx, raw_value) if raw_value else None
                if candidate and candidate.exists():
                    return candidate
                if raw_value:
                    resume_fallback_warnings.append(
                        f"{payload_key} recorded in stage_state could not be resolved on disk ({raw_value}); "
                        f"falling back to current input."
                    )
                return default_path

            current_ndx: Optional[Path] = None

            def _current_ndx() -> Optional[Path]:
                nonlocal current_ndx
                if current_ndx is not None:
                    return current_ndx
                try:
                    current_ndx = find_ndx_file(ctx)
                except ValueError as e:
                    raise RuntimeError(str(e)) from e
                return current_ndx

            try:
                resume_gro = _resolve_resume_input("gro_path", "gro_path", npt_gro)
                resume_top = _resolve_resume_input("top_path", "top_path", top_path)
                if "ndx_path" in old_payload and old_payload.get("ndx_path") is None:
                    resume_ndx = None
                else:
                    resume_ndx = _resolve_resume_input("ndx_path", "ndx_path", _current_ndx())
            except RuntimeError as e:
                print(f"  [ERROR] {e}")
                return False

            for warn_msg in resume_fallback_warnings:
                print(f"  [WARN] {warn_msg}")

            stage_overrides, extra_overrides_resolved, effective_overrides = resolve_mdp_overrides(
                ctx, "md", extra_overrides=None
            )
            mdrun_cmd_base = build_mdrun_command(ctx, "md")
            grompp_cmd_str = old_state.get("grompp_command", "")
            stage_state = _build_stage_state(
                ctx,
                "md",
                mdp_path,
                resume_gro,
                resume_top,
                resume_ndx,
                grompp_cmd_str,
                mdrun_cmd_base,
                ctx.grompp_maxwarn,
                posres_info=old_state.get("posres"),
                stage_overrides=stage_overrides,
                extra_overrides=extra_overrides_resolved,
                effective_overrides=effective_overrides,
            )

            ignore_runtime = bool(getattr(ctx, "resume_ignore_runtime", False))
            old_for_compare, ignored_keys = _payload_for_resume_compare(old_payload, ignore_runtime)
            new_for_compare, _ = _payload_for_resume_compare(
                stage_state["fingerprint_payload"],
                ignore_runtime,
            )
            compare_keys = sorted(set(old_for_compare.keys()) & set(new_for_compare.keys()))
            schema_only_old = sorted(set(old_for_compare.keys()) - set(new_for_compare.keys()))
            schema_only_new = sorted(set(new_for_compare.keys()) - set(old_for_compare.keys()))
            if schema_only_old or schema_only_new:
                print(
                    "  [WARN] Resume fingerprint schema differs between recorded and current payload; "
                    "comparing common keys only."
                )
                if schema_only_old:
                    print("    - Keys only in recorded payload: " + ", ".join(schema_only_old))
                if schema_only_new:
                    print("    - Keys only in current payload: " + ", ".join(schema_only_new))
            old_for_compare = {k: old_for_compare[k] for k in compare_keys}
            new_for_compare = {k: new_for_compare[k] for k in compare_keys}
            runtime_only_diffs = []
            if ignore_runtime:
                runtime_only_diffs = [
                    key for key in ignored_keys
                    if old_payload.get(key) != stage_state["fingerprint_payload"].get(key)
                ]
                if runtime_only_diffs:
                    print(
                        "  [WARN] Resume fingerprint ignoring runtime-only differences: "
                        + ", ".join(runtime_only_diffs)
                    )
            diffs = _diff_payloads(old_for_compare, new_for_compare)
            if diffs:
                print("  [ERROR] Resume denied: fingerprint mismatch.")
                for diff in diffs:
                    print(f"    - {diff}")
                if ignore_runtime:
                    print("  Runtime-only fingerprint keys were ignored during this comparison.")
                print("  Use --force to restart cleanly.")
                return False

            if ctx.manifest:
                ctx.manifest.set("resume_checkpoint", {
                    "path": str(md_cpt),
                    "method": "cpi",
                })
            thermostat_audit = _collect_thermostat_audit(ctx, "gmx_prod", mdp_path, resume_ndx)

            mdrun_cmd_resume = build_mdrun_command(
                ctx, "md",
                extra_args=["-cpi", str(md_cpt), "-append"],
            )

            try:
                run_gromacs_command(ctx, mdrun_cmd_resume, output_dir, self.name, stream_logs=True, log_prefix="mdrun")
            except RuntimeError:
                return False

            if not ctx.dry_run:
                done_sentinel.write_text("MD complete\n")

            _record_stage_manifest(ctx, "gmx_prod", {
                "decision": "resume",
                "reason": (
                    "checkpoint-repair + fingerprint match"
                    if repair_performed else
                    "fingerprint match and checkpoint available"
                ),
                "stage_dir": str(output_dir),
                "mdp_path": str(mdp_path),
                "mdp_sha256": stage_state["mdp_sha256"],
                "mdp_original_values": {},
                "mdp_patched_values": {},
                "inputs": {
                    "gro": {"path": str(resume_gro), "sha256": stage_state["gro_sha256"]},
                    "top": {"path": str(resume_top), "sha256": stage_state["top_sha256"]},
                    "ndx": {
                        "path": str(resume_ndx) if resume_ndx else None,
                        "sha256": stage_state["ndx_sha256"],
                    },
                    "checkpoint": {
                        "path": str(md_cpt),
                        "required": True,
                        "allow_velocity_reset": ctx.allow_velocity_reset,
                    },
                },
                "grompp_maxwarn": ctx.grompp_maxwarn,
                "grompp_command": stage_state["grompp_command"],
                "mdrun_command": " ".join(mdrun_cmd_resume),
                "stage_state_path": str(stage_state_path),
                "resume_ignore_runtime": ignore_runtime,
                "ignored_runtime_differences": runtime_only_diffs,
                "fingerprint_path_fallbacks": resume_fallback_warnings,
                "resume_repair": repair_details,
                "posres": old_state.get("posres"),
                "velocity_reset": None,
                "thermostat": thermostat_audit,
            })
            _record_thermostat_manifest(ctx, "gmx_prod", thermostat_audit)

            repro = _build_repro_json(
                ctx=ctx,
                stage="md",
                stage_dir=output_dir,
                mdp_path=mdp_path,
                mdp_orig={},
                mdp_patched={},
                gro_path=resume_gro,
                top_path=resume_top,
                ndx_path=resume_ndx,
                cpt_path=md_cpt,
                grompp_cmd=stage_state["grompp_command"],
                mdrun_cmd=" ".join(mdrun_cmd_resume),
                resume_allowed=True,
            )
            
            # Cleanup BEFORE writing repro.json so it reflects final on-disk state
            deleted_files = cleanup_large_files(output_dir, ctx, stage_tag="md")
            if deleted_files:
                repro["cleanup"] = {
                    "trr_deleted": deleted_files,
                    "enabled": True,
                }
            
            _write_repro_json(ctx, output_dir, repro)
            _write_run_report(ctx)

            print(f"  [OK] Production MD resumed and completed")
            return True

        # Fresh start (no md.cpt/md.tpr)
        # Use shared decision function for checkpoint/velocity handling
        template_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "mdp")
        template_path = template_dir / "md.mdp"
        
        cpt_available, extra_overrides, velocity_mode, policy_warnings, policy_metadata = decide_prev_state_policy(
            stage_name="md",
            prev_cpt_path=npt_cpt,
            strict_mode=_is_restart_policy_strict(ctx),
            allow_velocity_reset=ctx.allow_velocity_reset,
            allow_gro_velocity_restart=getattr(ctx, "allow_gro_velocity_restart", False),
            template_mdp_path=template_path,
            ctx=ctx,
            prev_gro_path=npt_gro,
        )
        
        # Handle error case
        if velocity_mode == "error":
            for msg in policy_warnings:
                print(f"  [ERROR] {msg}")
            return False
        
        # Print warnings
        for w in policy_warnings:
            print(f"  [WARN] {w}")
        
        # Build velocity reset info for auditing
        velocity_reset_info = {
            "velocity_mode": velocity_mode,
            "warnings": policy_warnings,
            "policy_metadata": policy_metadata,
        }
        if extra_overrides:
            velocity_reset_info["overrides"] = extra_overrides

        try:
            ndx_path = find_ndx_file(ctx)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        ctx.active_ndx_path = str(ndx_path) if ndx_path else None
        ctx.input_gro_has_velocities = _gro_has_velocities(npt_gro)

        try:
            stage_overrides, extra_overrides_resolved, effective_overrides = resolve_mdp_overrides(
                ctx, "md", extra_overrides=extra_overrides
            )
            mdp_path, orig_vals, patched_vals = prepare_stage_mdp(
                ctx,
                "md",
                output_dir,
                effective_overrides=effective_overrides,
                checkpoint_available=cpt_available,
            )
            print(f"    MDP: {mdp_path.name}")
            if patched_vals:
                print(f"    Patched: {patched_vals}")
            if velocity_mode != "checkpoint":
                print(f"    Velocity mode: {velocity_mode}")
                if extra_overrides:
                    print(f"    Velocity overrides: {extra_overrides}")
        except MDPPatchError as e:
            print(f"  [ERROR] {e}")
            return False
        # Only validate velocity reset MDP settings when in velocity_reset mode
        if velocity_mode == "velocity_reset" and extra_overrides:
            try:
                expected_gen_temp = extra_overrides.get("gen_temp")
                _assert_velocity_reset_mdp(mdp_path, "md", expected_gen_temp)
            except MDPPatchError as e:
                print(f"  [ERROR] {e}")
                return False
        try:
            _assert_ndx_for_mdp(mdp_path, ndx_path)
        except ValueError as e:
            print(f"  [ERROR] {e}")
            return False
        thermostat_audit = _collect_thermostat_audit(ctx, "gmx_prod", mdp_path, ndx_path)

        # POSRES reference handling
        try:
            posres_ref, posres_info = _resolve_posres_info(ctx, "md", npt_gro, mdp_path)
        except (ValueError, FileNotFoundError) as e:
            print(f"  [ERROR] {e}")
            return False

        grompp_cmd = build_grompp_command(
            ctx=ctx,
            mdp_path=str(mdp_path),
            structure_path=str(npt_gro),
            topology_path=str(top_path),
            output_path=str(md_tpr),
            checkpoint_path=str(npt_cpt) if npt_cpt.exists() else None,
            index_path=str(ndx_path) if ndx_path else None,
            reference_structure_path=str(posres_ref) if posres_ref else None,
            maxwarn=ctx.grompp_maxwarn,
        )
        mdrun_cmd_base = build_mdrun_command(ctx, "md")
        stage_state = _build_stage_state(
            ctx,
            "md",
            mdp_path,
            npt_gro,
            top_path,
            ndx_path,
            grompp_cmd,
            mdrun_cmd_base,
            ctx.grompp_maxwarn,
            posres_info=posres_info,
            stage_overrides=stage_overrides,
            extra_overrides=extra_overrides_resolved,
            effective_overrides=effective_overrides,
            velocity_mode=velocity_mode,
        )

        # Fresh start: record resume checkpoint (if available)
        if ctx.manifest and npt_cpt.exists():
            ctx.manifest.set("resume_checkpoint", {
                "path": str(npt_cpt),
                "method": "t",
            })

        try:
            run_gromacs_command(ctx, grompp_cmd, output_dir, self.name)
        except RuntimeError:
            return False

        if not ctx.dry_run:
            stage_state = _build_stage_state(
                ctx,
                "md",
                mdp_path,
                npt_gro,
                top_path,
                ndx_path,
                grompp_cmd,
                mdrun_cmd_base,
                ctx.grompp_maxwarn,
                posres_info=posres_info,
                stage_overrides=stage_overrides,
                extra_overrides=extra_overrides_resolved,
                effective_overrides=effective_overrides,
                velocity_mode=velocity_mode,
            )
            _write_stage_state(stage_state_path, stage_state, ctx)

        try:
            run_gromacs_command(ctx, mdrun_cmd_base, output_dir, self.name, stream_logs=True, log_prefix="mdrun")
        except RuntimeError:
            return False

        if not ctx.dry_run:
            done_sentinel.write_text("MD complete\n")
            stage_state = _build_stage_state(
                ctx,
                "md",
                mdp_path,
                npt_gro,
                top_path,
                ndx_path,
                grompp_cmd,
                mdrun_cmd_base,
                ctx.grompp_maxwarn,
                posres_info=posres_info,
                stage_overrides=stage_overrides,
                extra_overrides=extra_overrides_resolved,
                effective_overrides=effective_overrides,
                velocity_mode=velocity_mode,
            )
            _write_stage_state(stage_state_path, stage_state, ctx)

        if ctx.manifest:
            ctx.manifest.set("resources", get_resource_summary(ctx))

        _record_stage_manifest(ctx, "gmx_prod", {
            "decision": "run",
            "reason": "fresh start",
            "stage_dir": str(output_dir),
            "mdp_path": str(mdp_path),
            "mdp_sha256": stage_state["mdp_sha256"],
            "mdp_original_values": orig_vals,
            "mdp_patched_values": patched_vals,
            "inputs": {
                "gro": {"path": str(npt_gro), "sha256": stage_state["gro_sha256"]},
                "top": {"path": str(top_path), "sha256": stage_state["top_sha256"]},
                "ndx": {
                    "path": str(ndx_path) if ndx_path else None,
                    "sha256": stage_state["ndx_sha256"],
                },
                "checkpoint": {
                    "path": str(npt_cpt) if npt_cpt.exists() else None,
                    "required": True,
                    "allow_velocity_reset": ctx.allow_velocity_reset,
                },
            },
            "grompp_maxwarn": ctx.grompp_maxwarn,
            "grompp_command": stage_state["grompp_command"],
            "mdrun_command": " ".join(mdrun_cmd_base),
            "stage_state_path": str(stage_state_path),
            "posres": posres_info,
            "velocity_reset": velocity_reset_info,
            "thermostat": thermostat_audit,
        })
        _record_thermostat_manifest(ctx, "gmx_prod", thermostat_audit)

        repro = _build_repro_json(
            ctx=ctx,
            stage="md",
            stage_dir=output_dir,
            mdp_path=mdp_path,
            mdp_orig=orig_vals,
            mdp_patched=patched_vals,
            gro_path=npt_gro,
            top_path=top_path,
            ndx_path=ndx_path,
            cpt_path=npt_cpt if cpt_available else None,
            grompp_cmd=stage_state["grompp_command"],
            mdrun_cmd=" ".join(mdrun_cmd_base),
            resume_allowed=False,
        )
        
        # Cleanup BEFORE writing repro.json so it reflects final on-disk state
        deleted_files = cleanup_large_files(output_dir, ctx, stage_tag="md")
        if deleted_files:
            repro["cleanup"] = {
                "trr_deleted": deleted_files,
                "enabled": True,
            }
        
        _write_repro_json(ctx, output_dir, repro)
        _write_run_report(ctx)

        print(f"  [OK] Production MD complete")
        return True
