"""
Resource scheduler for GROMACS command args/env and HTPolyNet config injection.

Builds deterministic `mdrun` arguments plus environment variables from context:
- --gmx-nt -> mdrun -nt (thread-MPI builds only)
- --gmx-gpu-id -> mdrun -gpu_id
- --omp-num-threads -> mdrun -ntomp and OMP_NUM_THREADS
"""

import os
import shlex
from typing import Dict, List, Optional, TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .context import PipelineContext
from .gromacs_cmd import resolve_gmx_command


def _parse_positive_int(name: str, value: Optional[Any], min_value: int = 1) -> Optional[int]:
    """Parse an integer-like value and enforce a minimum value."""
    if value is None:
        return None

    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer >= {min_value}, not bool (got {value!r})")

    parsed: Optional[int] = None
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            raise ValueError(
                f"{name} must be an integer >= {min_value}; empty string is not valid."
            )
        sign_stripped = stripped[1:] if stripped.startswith(("+", "-")) else stripped
        if not sign_stripped.isdigit():
            raise ValueError(
                f"{name} must be a base-10 integer >= {min_value} (got {value!r})."
            )
        parsed = int(stripped, 10)
    else:
        raise ValueError(
            f"{name} must be an integer >= {min_value} or string integer (got {value!r})."
        )

    if parsed < min_value:
        raise ValueError(f"{name} must be an integer >= {min_value} (got {parsed})")
    return parsed


def _parse_bool_flag(name: str, value: Optional[Any], default: bool = False) -> bool:
    """Parse bool / bool-like strings for config-driven flags."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "on"}:
            return True
        if lowered in {"0", "false", "no", "off", ""}:
            return False
    raise ValueError(
        f"{name} must be a boolean or one of "
        f"'true/false/1/0/yes/no/on/off' (got {value!r})."
    )


def _validate_gpu_id(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    if not isinstance(value, str):
        value = str(value)
    if value == "" or value.strip() == "":
        raise ValueError("--gmx-gpu-id must be a non-empty string (no spaces)")
    if any(ch.isspace() for ch in value):
        raise ValueError(
            f"--gmx-gpu-id contains whitespace: {value!r}. "
            "Use comma-separated IDs with no spaces (see README)."
        )
    return value


def _resources_requested(ctx: "PipelineContext") -> bool:
    return (
        ctx.gmx_nt is not None
        or ctx.gmx_gpu_id is not None
        or ctx.omp_num_threads is not None
    )


def _get_gmx_basename(gmx_cmd: str) -> str:
    """Return basename of resolved gmx command token."""
    raw = str(gmx_cmd).strip()
    if not raw:
        return ""
    try:
        tokens = shlex.split(raw)
    except ValueError:
        tokens = raw.split()
    if not tokens:
        return ""
    return os.path.basename(tokens[0])


def _is_mpi_gromacs_command(gmx_cmd: str) -> bool:
    """
    Conservative external-MPI detection from resolved gmx command basename.

    We treat names containing 'gmx_mpi' or ending in '_mpi' as MPI builds.
    """
    base = _get_gmx_basename(gmx_cmd).lower()
    return "gmx_mpi" in base or base.endswith("_mpi")


def _validate_mpi_threading_compatibility(gmx_cmd: str, gmx_nt: Optional[int]) -> None:
    """Fail-fast when thread-MPI flags are used with external MPI GROMACS."""
    if gmx_nt is None:
        return
    if _is_mpi_gromacs_command(gmx_cmd):
        raise ValueError(
            f"--gmx-nt={gmx_nt} is incompatible with MPI GROMACS command '{gmx_cmd}'. "
            "External MPI builds use mpirun/srun to set MPI ranks, and mdrun -nt "
            "is not valid in that mode. Use --omp-num-threads to control threads "
            "per rank (-ntomp and OMP_NUM_THREADS)."
        )


def get_mdrun_args(ctx: "PipelineContext", gmx_cmd: Optional[str] = None) -> List[str]:
    """
    Build mdrun arguments from context.
    
    Args:
        ctx: Pipeline context with resource options
        
    Returns:
        List of additional mdrun arguments
    """
    args = []

    if gmx_cmd is None:
        gmx_cmd = resolve_gmx_command(ctx)

    gmx_nt = _parse_positive_int("gmx_nt", ctx.gmx_nt)
    omp_num_threads = _parse_positive_int("omp_num_threads", ctx.omp_num_threads)
    gmx_gpu_id = _validate_gpu_id(ctx.gmx_gpu_id)
    _validate_mpi_threading_compatibility(gmx_cmd, gmx_nt)

    if gmx_nt is not None:
        args.extend(["-nt", str(gmx_nt)])

    if omp_num_threads is not None:
        args.extend(["-ntomp", str(omp_num_threads)])

    if gmx_gpu_id is not None:
        args.extend(["-gpu_id", gmx_gpu_id])
    
    return args


def get_env_vars(ctx: "PipelineContext", env_fallback: bool = False) -> Dict[str, str]:
    """
    Get environment variables for GROMACS execution.
    
    Args:
        ctx: Pipeline context
        
    Returns:
        Dictionary of environment variables to set
    """
    env = {}

    _parse_positive_int("gmx_nt", ctx.gmx_nt)
    omp_num_threads = _parse_positive_int("omp_num_threads", ctx.omp_num_threads)
    gmx_gpu_id = _validate_gpu_id(ctx.gmx_gpu_id)

    if omp_num_threads is not None:
        env["OMP_NUM_THREADS"] = str(omp_num_threads)

    if env_fallback and gmx_gpu_id is not None:
        env["GMX_GPU_ID"] = gmx_gpu_id
    
    return env


def get_execution_env(ctx: "PipelineContext", env_fallback: bool = False) -> Dict[str, str]:
    """
    Get complete execution environment merging current env with overrides.
    
    Args:
        ctx: Pipeline context
        
    Returns:
        Complete environment dict for subprocess execution
    """
    env = dict(os.environ)
    env.update(get_env_vars(ctx, env_fallback=env_fallback))
    return env


def inject_into_htpolynet_config(
    config_dict: dict,
    ctx: "PipelineContext",
) -> str:
    """
    Inject resource settings into HTPolyNet config if supported.
    
    HTPolyNet may support passing mdrun arguments via its config.
    This function attempts to add them if the config structure allows.
    
    Args:
        config_dict: HTPolyNet config dictionary
        ctx: Pipeline context
        
    Returns:
        Injection method used: 'config', 'env', or 'none'
    """
    method = "none"

    _parse_positive_int("gmx_nt", ctx.gmx_nt)
    omp_num_threads = _parse_positive_int("omp_num_threads", ctx.omp_num_threads)
    gmx_gpu_id = _validate_gpu_id(ctx.gmx_gpu_id)
    add_missing = _parse_bool_flag(
        "htpolynet_inject_add_missing",
        getattr(ctx, "htpolynet_inject_add_missing", False),
        default=False,
    )

    # Backward-safe strict mode: only inject existing keys.
    # Opt-in mode: create gromacs section and add ntomp/gpu_id if missing.
    gromacs_section = None
    if "gromacs" in config_dict and isinstance(config_dict["gromacs"], dict):
        gromacs_section = config_dict["gromacs"]
    elif add_missing:
        config_dict["gromacs"] = {}
        gromacs_section = config_dict["gromacs"]

    if isinstance(gromacs_section, dict):
        if add_missing:
            if omp_num_threads is not None:
                gromacs_section["ntomp"] = omp_num_threads
                method = "config"
            if gmx_gpu_id is not None:
                gromacs_section["gpu_id"] = gmx_gpu_id
                method = "config"
        else:
            if "ntomp" in gromacs_section and omp_num_threads is not None:
                gromacs_section["ntomp"] = omp_num_threads
                method = "config"

            if "gpu_id" in gromacs_section and gmx_gpu_id is not None:
                gromacs_section["gpu_id"] = gmx_gpu_id
                method = "config"

    # If we couldn't inject via config, we'll use environment
    if method == "none" and _resources_requested(ctx):
        method = "env"

    setattr(ctx, "htpolynet_injection_method", method)
    
    return method


def get_resource_summary(ctx: "PipelineContext") -> Dict[str, Any]:
    """
    Get summary of resource settings for manifest recording.
    
    Args:
        ctx: Pipeline context
        
    Returns:
        Dictionary suitable for manifest recording
    """
    gmx_cmd = resolve_gmx_command(ctx)
    gmx_nt = _parse_positive_int("gmx_nt", ctx.gmx_nt)
    omp_num_threads = _parse_positive_int("omp_num_threads", ctx.omp_num_threads)
    gmx_gpu_id = _validate_gpu_id(ctx.gmx_gpu_id)
    gmx_is_mpi = _is_mpi_gromacs_command(gmx_cmd)
    _validate_mpi_threading_compatibility(gmx_cmd, gmx_nt)

    add_missing = _parse_bool_flag(
        "htpolynet_inject_add_missing",
        getattr(ctx, "htpolynet_inject_add_missing", False),
        default=False,
    )
    injection_method = getattr(ctx, "htpolynet_injection_method", None)
    if injection_method not in {"config", "env", "none"}:
        if not _resources_requested(ctx):
            injection_method = "none"
        elif add_missing and (omp_num_threads is not None or gmx_gpu_id is not None):
            injection_method = "config"
        else:
            injection_method = "env"

    return {
        "gmx_nt": str(gmx_nt) if gmx_nt is not None else None,
        "gmx_gpu_id": gmx_gpu_id,
        "omp_num_threads": str(omp_num_threads) if omp_num_threads is not None else None,
        "gmx_is_mpi": gmx_is_mpi,
        "htpolynet_injection_method": injection_method,
        "htpolynet_inject_add_missing": add_missing,
    }


def build_mdrun_command(
    ctx: "PipelineContext",
    deffnm: str,
    extra_args: Optional[List[str]] = None,
) -> List[str]:
    """
    Build complete mdrun command with resource arguments.
    
    Args:
        ctx: Pipeline context
        deffnm: Default filename prefix for mdrun
        extra_args: Additional arguments (e.g., -cpi for checkpoint)
        
    Returns:
        Complete command as list of strings
    """
    gmx_cmd = resolve_gmx_command(ctx)
    cmd = [gmx_cmd, "mdrun", "-deffnm", deffnm]
    
    # Add resource args
    cmd.extend(get_mdrun_args(ctx, gmx_cmd=gmx_cmd))
    
    # Add extra args
    if extra_args:
        cmd.extend(extra_args)
    
    return cmd


def build_grompp_command(
    ctx: "PipelineContext",
    mdp_path: str,
    structure_path: str,
    topology_path: str,
    output_path: str,
    checkpoint_path: Optional[str] = None,
    index_path: Optional[str] = None,
    reference_structure_path: Optional[str] = None,
    maxwarn: int = 0,
) -> List[str]:
    """
    Build grompp command.
    
    Args:
        ctx: Pipeline context
        mdp_path: Path to MDP file
        structure_path: Path to input structure (.gro)
        topology_path: Path to topology (.top)
        output_path: Path for output .tpr
        checkpoint_path: Optional checkpoint for -t (velocities)
        index_path: Optional index file
        reference_structure_path: Optional reference structure for -r (POSRES)
        maxwarn: Maximum warnings to allow
        
    Returns:
        Complete grompp command as list
    """
    gmx_cmd = resolve_gmx_command(ctx)
    cmd = [
        gmx_cmd, "grompp",
        "-f", mdp_path,
        "-c", structure_path,
        "-p", topology_path,
        "-o", output_path,
    ]
    
    if checkpoint_path:
        cmd.extend(["-t", checkpoint_path])
    
    if index_path:
        cmd.extend(["-n", index_path])
    
    if reference_structure_path:
        cmd.extend(["-r", reference_structure_path])
    
    if maxwarn > 0:
        cmd.extend(["-maxwarn", str(maxwarn)])
    
    return cmd
