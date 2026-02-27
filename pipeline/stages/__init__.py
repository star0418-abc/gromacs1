"""
Stage registry and base classes.
"""

from .base import BaseStage
from .packmol import PackmolStage
from .htpolynet import HTPolyNetStage
from .sanitizer import SanitizerStage
from .gromacs import GmxEmStage, GmxEqStage, GmxProdStage
from .analysis import AnalysisStage


# Stage execution order
STAGE_ORDER = [
    "packmol",
    "htpolynet",
    "sanitizer",  # NEW: ITP sanitizer (after staging, before grompp)
    "gmx_em",
    "gmx_eq",
    "gmx_prod",
    "analysis",
]

# Stage name to class mapping
STAGE_REGISTRY = {
    "packmol": PackmolStage,
    "htpolynet": HTPolyNetStage,
    "sanitizer": SanitizerStage,
    "gmx_em": GmxEmStage,
    "gmx_eq": GmxEqStage,
    "gmx_prod": GmxProdStage,
    "analysis": AnalysisStage,
}


def get_stage(name: str) -> BaseStage:
    """
    Get a stage instance by name.
    
    Args:
        name: Stage name
        
    Returns:
        Stage instance
        
    Raises:
        KeyError: If stage name not found
    """
    if name not in STAGE_REGISTRY:
        raise KeyError(f"Unknown stage: {name}. Valid stages: {list(STAGE_REGISTRY.keys())}")
    return STAGE_REGISTRY[name]()


def get_all_stages():
    """Get all stages in execution order."""
    return [STAGE_REGISTRY[name]() for name in STAGE_ORDER]


__all__ = [
    "STAGE_ORDER",
    "STAGE_REGISTRY", 
    "get_stage",
    "get_all_stages",
    "BaseStage",
    "PackmolStage",
    "HTPolyNetStage",
    "SanitizerStage",
    "GmxEmStage",
    "GmxEqStage",
    "GmxProdStage",
    "AnalysisStage",
]

