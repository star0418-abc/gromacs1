"""
GROMACS executable resolver.

Priority:
1) Context override (if present)
2) GMX_CMD environment variable
3) Default: "gmx"
"""

from __future__ import annotations

import os
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .context import PipelineContext


_CTX_GMX_FIELDS = (
    "gmx_cmd",
    "gromacs_cmd",
    "gmx_executable",
    "gmx_bin",
    "gromacs_bin",
)


def resolve_gmx_command(ctx: Optional["PipelineContext"] = None) -> str:
    """Resolve the GROMACS executable to use."""
    if ctx is not None:
        for field in _CTX_GMX_FIELDS:
            value = getattr(ctx, field, None)
            if value:
                return str(value)

    env_cmd = os.environ.get("GMX_CMD")
    if env_cmd:
        env_cmd = env_cmd.strip()
        if env_cmd:
            return env_cmd

    return "gmx"
