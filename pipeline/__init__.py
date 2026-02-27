"""
Pipeline package for PACKMOL → HTPOLYNET → GROMACS → RDF/CN automation.
"""

__version__ = "0.3.0"
__author__ = "Pipeline Automation"

from .path_manager import PathManager, PathViolationError
from .manifest import ManifestWriter
from .context import PipelineContext
from .dispatcher import StageDispatcher
from .composition import (
    MoleculeSpec,
    CompositionConfig,
    CompositionResult,
    BoxConfig,
    get_default_recipe,
    wt_to_counts,
    compute_box_size,
    compute_effective_density,
    get_density_candidates,
    compute_full_composition,
)
from .packmol_input import (
    resolve_molecule_template,
    generate_packmol_input,
    write_packmol_input,
    run_packmol,
    run_packmol_with_retry,
    PackmolError,
)
from .htpolynet_config import (
    HTPolyNetConfig,
    generate_system_config,
    write_system_config,
    resolve_rules_bundle,
)

__all__ = [
    # Core
    "PathManager",
    "PathViolationError", 
    "ManifestWriter",
    "PipelineContext",
    "StageDispatcher",
    # Composition
    "MoleculeSpec",
    "CompositionConfig",
    "CompositionResult",
    "BoxConfig",
    "get_default_recipe",
    "wt_to_counts",
    "compute_box_size",
    "compute_effective_density",
    "get_density_candidates",
    "compute_full_composition",
    # PACKMOL
    "resolve_molecule_template",
    "generate_packmol_input",
    "write_packmol_input",
    "run_packmol",
    "run_packmol_with_retry",
    "PackmolError",
    # HTPolyNet
    "HTPolyNetConfig",
    "generate_system_config",
    "write_system_config",
    "resolve_rules_bundle",
]
