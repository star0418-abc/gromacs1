"""
Composition model: wt% -> integer molecule counts conversion.

Handles:
- wt% specification for gel polymer electrolyte components
- Representative MW selection for polydisperse materials (OEGMA, PEGDMA)
- Integer rounding with minimum count enforcement
- Auto-scaling when constraints are violated
- Box sizing from total mass and target density
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import math
import warnings


# =============================================================================
# Constants
# =============================================================================

AVOGADRO = 6.02214076e23  # mol^-1
RECOMMENDED_MOLECULE_ROLES = {"monomer", "crosslinker", "solvent", "salt", "additive"}
COMMON_POLYDISPERSE_PATTERNS = ("OEGMA", "PEGDMA")


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class MoleculeSpec:
    """
    Specification for a molecule in the composition.

    Recommended roles: monomer, crosslinker, solvent, salt, additive.
    Freeform role strings are accepted for backward compatibility.
    """

    name: str
    mw_g_mol: float  # Molecular weight in g/mol
    wt_pct: float  # Weight percent (0-100)
    min_count: int = 1  # Minimum molecule count threshold
    role: Optional[str] = None
    mw_note: Optional[str] = None  # e.g., "representative MW for polydisperse material"
    pdi: Optional[float] = None  # Record-only polydispersity index

    def __post_init__(self):
        if self.wt_pct < 0 or self.wt_pct > 100:
            raise ValueError(f"wt_pct must be in [0, 100], got {self.wt_pct}")
        if self.mw_g_mol <= 0:
            raise ValueError(f"mw_g_mol must be positive, got {self.mw_g_mol}")
        if not isinstance(self.min_count, int) or isinstance(self.min_count, bool):
            raise ValueError(f"min_count must be an integer >= 0, got {self.min_count!r}")
        if self.min_count < 0:
            raise ValueError(f"min_count must be >= 0, got {self.min_count}")
        if self.role is not None and not isinstance(self.role, str):
            raise ValueError(f"role must be a string when provided, got {type(self.role).__name__}")
        if self.mw_note is not None and not isinstance(self.mw_note, str):
            raise ValueError(f"mw_note must be a string when provided, got {type(self.mw_note).__name__}")
        if self.pdi is not None and self.pdi <= 0:
            raise ValueError(f"pdi must be > 0 when provided, got {self.pdi}")


@dataclass
class CompositionConfig:
    """
    Full composition configuration for a system.
    
    Attributes:
        molecules: List of MoleculeSpec for each component
        anchor_name: Name of the anchor molecule for scaling (optional)
        anchor_count: Target count for anchor molecule (optional)
        target_total_molecules: Alternative: target total molecule count
        wt_deviation_threshold: Max allowed wt% deviation after rounding (%)
        rho0_g_cm3: Initial density guess for box sizing (g/cm³)
        rho0_min_g_cm3: Minimum density for retry strategy
        rho0_max_g_cm3: Maximum density for retry strategy
        gel_min_crosslinker_count: Minimum expected crosslinker count (0 disables)
        gel_policy: How to handle gel guardrail violations (off|warn|error)
        auto_scale_to_gel_min_crosslinker: If True, scale up to satisfy gel minimum
        polymerization_shrinkage_vol_frac: Expected volumetric shrinkage fraction
    """
    molecules: List[MoleculeSpec]
    anchor_name: Optional[str] = None
    anchor_count: Optional[int] = None
    target_total_molecules: Optional[int] = 500
    wt_deviation_threshold: float = 2.0  # 2% max deviation
    rho0_g_cm3: float = 1.0  # g/cm³
    rho0_min_g_cm3: float = 0.95
    rho0_max_g_cm3: float = 1.05
    gel_min_crosslinker_count: int = 0
    gel_policy: Optional[str] = None  # off | warn | error ; None => auto default
    auto_scale_to_gel_min_crosslinker: bool = False
    polymerization_shrinkage_vol_frac: float = 0.0

    def __post_init__(self):
        # Validate wt% sums to ~100
        total_wt = sum(m.wt_pct for m in self.molecules)
        if abs(total_wt - 100.0) > 0.1:
            raise ValueError(f"wt_pct must sum to 100, got {total_wt:.2f}")
        if not self.molecules:
            raise ValueError("At least one molecule must be provided")

        # Validate unique molecule names (no silent dict overwrite later).
        _validate_unique_molecule_names(self.molecules)

        if self.target_total_molecules is not None and self.target_total_molecules <= 0:
            raise ValueError(
                f"target_total_molecules must be None or > 0, got {self.target_total_molecules}"
            )
        if (self.anchor_name is None) != (self.anchor_count is None):
            raise ValueError("anchor_name and anchor_count must be provided together")
        if self.anchor_count is not None and self.anchor_count <= 0:
            raise ValueError(f"anchor_count must be > 0, got {self.anchor_count}")

        if self.rho0_g_cm3 <= 0:
            raise ValueError(f"rho0_g_cm3 must be > 0, got {self.rho0_g_cm3}")
        if self.rho0_min_g_cm3 <= 0 or self.rho0_max_g_cm3 <= 0:
            raise ValueError(
                f"rho0_min_g_cm3 and rho0_max_g_cm3 must be > 0, got "
                f"{self.rho0_min_g_cm3}, {self.rho0_max_g_cm3}"
            )
        if self.rho0_min_g_cm3 > self.rho0_max_g_cm3:
            raise ValueError(
                f"rho0_min_g_cm3 must be <= rho0_max_g_cm3, got "
                f"{self.rho0_min_g_cm3} > {self.rho0_max_g_cm3}"
            )

        if not isinstance(self.gel_min_crosslinker_count, int) or isinstance(
            self.gel_min_crosslinker_count,
            bool,
        ):
            raise ValueError(
                "gel_min_crosslinker_count must be an integer >= 0, "
                f"got {self.gel_min_crosslinker_count!r}"
            )
        if self.gel_min_crosslinker_count < 0:
            raise ValueError(
                f"gel_min_crosslinker_count must be >= 0, got {self.gel_min_crosslinker_count}"
            )

        allowed_policies = {"off", "warn", "error"}
        if self.gel_policy is None:
            self.gel_policy = "warn" if self.gel_min_crosslinker_count > 0 else "off"
        else:
            policy = str(self.gel_policy).strip().lower()
            if policy not in allowed_policies:
                raise ValueError(f"gel_policy must be one of {sorted(allowed_policies)}, got {self.gel_policy!r}")
            self.gel_policy = policy

        self.polymerization_shrinkage_vol_frac = _validate_shrinkage_vol_frac(
            self.polymerization_shrinkage_vol_frac
        )


@dataclass
class CompositionResult:
    """Result of wt% → count conversion."""
    counts: Dict[str, int]  # Molecule name → count
    mw_used: Dict[str, float]  # Molecule name → MW used
    target_wt_pct: Dict[str, float]  # Original target wt%
    actual_wt_pct: Dict[str, float]  # Actual wt% after rounding
    max_deviation_pct: float  # Max deviation from target
    total_mass_g: float  # Total system mass in grams
    total_molecules: int  # Total molecule count
    scale_factor: float  # Final scale factor applied
    warnings: Dict[str, List[str]] = field(default_factory=dict)
    gel_guardrail: Dict[str, Any] = field(default_factory=dict)
    auto_scaling: Dict[str, Any] = field(default_factory=dict)
    molecule_metadata: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Convert to dictionary for manifest."""
        data = {
            "counts": self.counts,
            "mw_used": self.mw_used,
            "target_wt_pct": self.target_wt_pct,
            "actual_wt_pct": self.actual_wt_pct,
            "max_deviation_pct": self.max_deviation_pct,
            "total_mass_g": self.total_mass_g,
            "total_molecules": self.total_molecules,
            "scale_factor": self.scale_factor,
        }
        if self.warnings:
            data["warnings"] = self.warnings
        if self.gel_guardrail:
            data["gel_guardrail"] = self.gel_guardrail
        if self.auto_scaling:
            data["auto_scaling"] = self.auto_scaling
        if self.molecule_metadata:
            data["molecule_metadata"] = self.molecule_metadata
        return data


@dataclass
class BoxConfig:
    """Box sizing configuration and result."""
    length_nm: float  # Cubic box side length in nm
    rho0_g_cm3: float  # Density used
    total_mass_g: float  # Total mass in grams
    volume_nm3: float  # Box volume in nm³
    requested_rho_g_cm3: Optional[float] = None
    polymerization_shrinkage_vol_frac: float = 0.0
    rho_effective_g_cm3: Optional[float] = None
    density_candidates_g_cm3: List[float] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to dictionary for manifest."""
        data = {
            "length_nm": self.length_nm,
            "rho0_g_cm3": self.rho0_g_cm3,
            "total_mass_g": self.total_mass_g,
            "volume_nm3": self.volume_nm3,
        }
        if self.requested_rho_g_cm3 is not None:
            data["requested_rho_g_cm3"] = self.requested_rho_g_cm3
        data["polymerization_shrinkage_vol_frac"] = self.polymerization_shrinkage_vol_frac
        data["rho_effective_g_cm3"] = (
            self.rho_effective_g_cm3
            if self.rho_effective_g_cm3 is not None
            else self.rho0_g_cm3
        )
        if self.density_candidates_g_cm3:
            data["density_candidates_g_cm3"] = self.density_candidates_g_cm3
        return data


# =============================================================================
# Default SYS_001 Recipe
# =============================================================================

def get_default_recipe() -> CompositionConfig:
    """
    Get the default recipe for SYS_001 gel polymer electrolyte.
    
    Components:
    - LiTFSI: Lithium salt
    - PC: Propylene carbonate (plasticizer)
    - tetraglyme: Co-solvent
    - OEGMA: Oligo(ethylene glycol) methacrylate monomer (polydisperse, use MW=475)
    - UEMA: Urethane ethyl methacrylate monomer
    - MPB: Monomer/initiator
    - PEGDMA: Poly(ethylene glycol) dimethacrylate crosslinker (polydisperse, use MW=750)
    - TPO-L: Photoinitiator
    
    This is a representative formulation. Adjust wt% for specific experiments.
    """
    return CompositionConfig(
        molecules=[
            MoleculeSpec(name="LiTFSI", mw_g_mol=287.09, wt_pct=15.0, min_count=10, role="salt"),
            MoleculeSpec(name="PC", mw_g_mol=102.09, wt_pct=25.0, min_count=20, role="solvent"),
            MoleculeSpec(name="tetraglyme", mw_g_mol=222.28, wt_pct=15.0, min_count=10, role="solvent"),
            MoleculeSpec(
                name="OEGMA",
                mw_g_mol=475.0,
                wt_pct=20.0,
                min_count=5,
                role="monomer",
                mw_note="representative MW for polydisperse material",
            ),
            MoleculeSpec(name="UEMA", mw_g_mol=201.22, wt_pct=10.0, min_count=5, role="monomer"),
            MoleculeSpec(name="MPB", mw_g_mol=178.19, wt_pct=5.0, min_count=3, role="additive"),
            MoleculeSpec(
                name="PEGDMA",
                mw_g_mol=750.0,
                wt_pct=8.0,
                min_count=3,
                role="crosslinker",
                mw_note="representative MW for polydisperse material",
            ),
            MoleculeSpec(name="TPO-L", mw_g_mol=316.37, wt_pct=2.0, min_count=2, role="additive"),
        ],
        target_total_molecules=500,
        wt_deviation_threshold=2.0,
        rho0_g_cm3=1.0,
        rho0_min_g_cm3=0.95,
        rho0_max_g_cm3=1.05,
    )


# =============================================================================
# Core Algorithms
# =============================================================================

def wt_to_counts(
    config: CompositionConfig,
    max_scale_iterations: int = 10,
) -> CompositionResult:
    """
    Convert weight percent composition to integer molecule counts.
    
    Algorithm:
    1. Compute mole fractions from wt%/MW
    2. Scale to target total molecules (or anchor count)
    3. Round using largest remainder method
    4. Enforce minimum counts, auto-scale if violated
    5. Recompute actual wt%, check deviation threshold
    
    Args:
        config: Composition configuration
        max_scale_iterations: Max iterations for min_count/deviation fixes
        
    Returns:
        CompositionResult with final counts and metadata
        
    Raises:
        ValueError: If composition cannot be satisfied within iterations
    """
    molecules = config.molecules
    _validate_unique_molecule_names(molecules)

    if config.target_total_molecules is not None and config.target_total_molecules <= 0:
        raise ValueError(
            f"target_total_molecules must be None or > 0, got {config.target_total_molecules}"
        )

    warning_messages: Dict[str, List[str]] = {}
    auto_scaling_events: List[Dict[str, Any]] = []

    # Step 1: Compute mole ratios (proportional to n_i = wt%_i / MW_i)
    mole_ratios = []
    for mol in molecules:
        if not isinstance(mol.min_count, int) or isinstance(mol.min_count, bool):
            raise ValueError(f"min_count must be an integer >= 0 for '{mol.name}', got {mol.min_count!r}")
        if mol.min_count < 0:
            raise ValueError(f"min_count must be >= 0 for '{mol.name}', got {mol.min_count}")
        if mol.wt_pct == 0.0:
            if mol.min_count > 0:
                raise ValueError(
                    f"Impossible constraint for '{mol.name}': wt_pct=0 gives mole_fraction=0, "
                    f"but min_count={mol.min_count} > 0."
                )
            ratio = 0.0
        else:
            ratio = mol.wt_pct / mol.mw_g_mol
        mole_ratios.append(ratio)

    total_ratio = sum(mole_ratios)
    if total_ratio <= 0:
        raise ValueError("Invalid composition: total wt_pct/MW ratio is zero; cannot compute mole fractions.")
    mole_fractions = [r / total_ratio for r in mole_ratios]

    crosslinker_indices = [
        i for i, mol in enumerate(molecules) if _normalize_role(mol.role) == "crosslinker"
    ]
    gel_policy = config.gel_policy or ("warn" if config.gel_min_crosslinker_count > 0 else "off")

    # Step 2: Determine initial scale factor
    if config.anchor_name is not None and config.anchor_count is not None:
        # Scale based on anchor molecule
        anchor_idx = next(
            (i for i, m in enumerate(molecules) if m.name == config.anchor_name),
            None
        )
        if anchor_idx is None:
            raise ValueError(f"Anchor molecule '{config.anchor_name}' not found")
        if mole_fractions[anchor_idx] <= 0:
            raise ValueError(
                f"Anchor molecule '{config.anchor_name}' has zero mole fraction; "
                "cannot compute scale factor from anchor_count."
            )
        scale_factor = config.anchor_count / mole_fractions[anchor_idx]
    else:
        # Scale based on target total
        scale_factor = (
            float(config.target_total_molecules)
            if config.target_total_molecules is not None
            else 500.0
        )

    # Iteratively scale and round until constraints are satisfied
    for iteration in range(max_scale_iterations):
        # Step 3: Compute raw counts and round with largest remainder
        raw_counts = [f * scale_factor for f in mole_fractions]
        counts = _largest_remainder_round(raw_counts)

        # Step 4: Check minimum count constraints
        min_violated = False
        for i, mol in enumerate(molecules):
            if counts[i] < mol.min_count:
                if mole_fractions[i] <= 0:
                    raise ValueError(
                        f"Impossible constraint for '{mol.name}': mole_fraction=0 and "
                        f"min_count={mol.min_count} > 0."
                    )
                # Scale up so this component meets minimum
                required_scale = mol.min_count / mole_fractions[i]
                if required_scale > scale_factor:
                    previous_scale = scale_factor
                    scale_factor = required_scale * 1.1  # 10% buffer
                    auto_scaling_events.append(
                        {
                            "iteration": iteration,
                            "reason": "min_count",
                            "molecule": mol.name,
                            "required_scale": required_scale,
                            "previous_scale_factor": previous_scale,
                            "new_scale_factor": scale_factor,
                        }
                    )
                    min_violated = True
                    break

        if (
            not min_violated
            and config.auto_scale_to_gel_min_crosslinker
            and config.gel_min_crosslinker_count > 0
            and crosslinker_indices
        ):
            gel_required_scale = scale_factor
            limiting_idx: Optional[int] = None
            limiting_expected = None
            for idx in crosslinker_indices:
                fraction = mole_fractions[idx]
                if fraction <= 0:
                    raise ValueError(
                        f"Impossible gel constraint for '{molecules[idx].name}': mole_fraction=0 "
                        f"cannot satisfy gel_min_crosslinker_count={config.gel_min_crosslinker_count}."
                    )
                expected_count = fraction * scale_factor
                required_scale = config.gel_min_crosslinker_count / fraction
                if required_scale > gel_required_scale:
                    gel_required_scale = required_scale
                    limiting_idx = idx
                    limiting_expected = expected_count
            if gel_required_scale > scale_factor:
                previous_scale = scale_factor
                scale_factor = gel_required_scale * 1.1  # same buffer style as min_count logic
                if limiting_idx is not None:
                    auto_scaling_events.append(
                        {
                            "iteration": iteration,
                            "reason": "gel_min_crosslinker",
                            "policy": "auto_scale_to_gel_min_crosslinker",
                            "molecule": molecules[limiting_idx].name,
                            "gel_min_crosslinker_count": config.gel_min_crosslinker_count,
                            "expected_count_before_scale": limiting_expected,
                            "required_scale": gel_required_scale,
                            "previous_scale_factor": previous_scale,
                            "new_scale_factor": scale_factor,
                        }
                    )
                min_violated = True

        if min_violated:
            continue  # Re-compute with new scale

        # Step 5: Compute actual wt% and check deviation
        actual_wt = _compute_actual_wt(molecules, counts)
        max_dev = max(
            abs(actual_wt[m.name] - m.wt_pct)
            for m in molecules
        )

        if max_dev <= config.wt_deviation_threshold:
            # Success! Build result
            total_mass = sum(
                counts[i] * mol.mw_g_mol / AVOGADRO
                for i, mol in enumerate(molecules)
            )

            gel_guardrail = _evaluate_gel_guardrail(
                molecules=molecules,
                mole_fractions=mole_fractions,
                scale_factor=scale_factor,
                min_crosslinker_count=config.gel_min_crosslinker_count,
                gel_policy=gel_policy,
                warning_messages=warning_messages,
            )
            molecule_metadata = _collect_molecule_metadata(molecules)
            _evaluate_polydispersity_warnings(
                molecules=molecules,
                warning_messages=warning_messages,
            )

            auto_scaling: Dict[str, Any] = {}
            if auto_scaling_events:
                auto_scaling = {
                    "enabled": bool(config.auto_scale_to_gel_min_crosslinker),
                    "events": auto_scaling_events,
                }

            return CompositionResult(
                counts={mol.name: counts[i] for i, mol in enumerate(molecules)},
                mw_used={mol.name: mol.mw_g_mol for mol in molecules},
                target_wt_pct={mol.name: mol.wt_pct for mol in molecules},
                actual_wt_pct=actual_wt,
                max_deviation_pct=max_dev,
                total_mass_g=total_mass,
                total_molecules=sum(counts),
                scale_factor=scale_factor,
                warnings=warning_messages,
                gel_guardrail=gel_guardrail,
                auto_scaling=auto_scaling,
                molecule_metadata=molecule_metadata,
            )

        # Deviation too high, scale up to get better resolution
        scale_factor *= 1.5

    raise ValueError(
        f"Could not satisfy composition constraints after {max_scale_iterations} iterations. "
        f"Try increasing target_total_molecules or relaxing wt_deviation_threshold."
    )


def _largest_remainder_round(values: List[float]) -> List[int]:
    """
    Round a list of floats to integers preserving the sum.
    
    Uses the largest remainder method (Hamilton's method).
    """
    floors = [int(math.floor(v)) for v in values]
    remainders = [v - f for v, f in zip(values, floors)]

    # How many units to distribute
    total_floor = sum(floors)
    target_sum = round(sum(values))
    to_distribute = target_sum - total_floor

    # Sort deterministically (value, then index)
    indices = list(range(len(values)))
    result = floors.copy()
    if not indices:
        return result

    if to_distribute > 0:
        # Give extra units to indices with largest remainders
        indices.sort(key=lambda i: (-remainders[i], i))
        for i in range(to_distribute):
            result[indices[i % len(indices)]] += 1
    elif to_distribute < 0:
        # Remove units from indices with smallest remainders
        indices.sort(key=lambda i: (remainders[i], i))
        for i in range(-to_distribute):
            result[indices[i % len(indices)]] -= 1

    return result


def _compute_actual_wt(
    molecules: List[MoleculeSpec], 
    counts: List[int]
) -> Dict[str, float]:
    """Compute actual weight percentages from counts."""
    masses = [c * mol.mw_g_mol for c, mol in zip(counts, molecules)]
    total_mass = sum(masses)
    
    if total_mass == 0:
        return {mol.name: 0.0 for mol in molecules}
    
    return {
        mol.name: (masses[i] / total_mass) * 100.0
        for i, mol in enumerate(molecules)
    }


def compute_box_size(
    total_mass_g: float,
    rho_g_cm3: float,
    polymerization_shrinkage_vol_frac: float = 0.0,
    density_candidates_g_cm3: Optional[List[float]] = None,
) -> BoxConfig:
    """
    Compute cubic box size from total mass and target density.
    
    Args:
        total_mass_g: Total system mass in grams
        rho_g_cm3: Target density in g/cm³
        polymerization_shrinkage_vol_frac: Expected post-polymerization volumetric
            shrinkage (0.0 means disabled)
        density_candidates_g_cm3: Optional density candidates for retry workflows
        
    Returns:
        BoxConfig with computed dimensions
    """
    if total_mass_g <= 0:
        raise ValueError(f"total_mass_g must be > 0, got {total_mass_g}")
    rho_effective = compute_effective_density(
        rho_g_cm3,
        polymerization_shrinkage_vol_frac=polymerization_shrinkage_vol_frac,
    )

    # Volume in cm³
    volume_cm3 = total_mass_g / rho_effective

    # Convert to nm³ (1 cm = 10^7 nm, so 1 cm³ = 10^21 nm³)
    volume_nm3 = volume_cm3 * 1e21

    # Cubic box side length
    length_nm = volume_nm3 ** (1/3)

    return BoxConfig(
        length_nm=length_nm,
        rho0_g_cm3=rho_effective,
        total_mass_g=total_mass_g,
        volume_nm3=volume_nm3,
        requested_rho_g_cm3=rho_g_cm3,
        polymerization_shrinkage_vol_frac=polymerization_shrinkage_vol_frac,
        rho_effective_g_cm3=rho_effective,
        density_candidates_g_cm3=list(density_candidates_g_cm3 or []),
    )


# =============================================================================
# Convenience Functions
# =============================================================================

def compute_full_composition(
    config: Optional[CompositionConfig] = None,
) -> Tuple[CompositionResult, BoxConfig]:
    """
    Compute both molecule counts and box size for a composition.
    
    Args:
        config: Composition config (uses default SYS_001 if None)
        
    Returns:
        Tuple of (CompositionResult, BoxConfig)
    """
    if config is None:
        config = get_default_recipe()

    comp_result = wt_to_counts(config)
    density_candidates = get_density_candidates(config)
    box_result = compute_box_size(
        comp_result.total_mass_g,
        config.rho0_g_cm3,
        polymerization_shrinkage_vol_frac=config.polymerization_shrinkage_vol_frac,
        density_candidates_g_cm3=density_candidates,
    )

    return comp_result, box_result


def print_composition_summary(
    comp: CompositionResult,
    box: BoxConfig,
) -> None:
    """Print a human-readable summary of the composition."""
    print("\n" + "=" * 60)
    print("COMPOSITION SUMMARY")
    print("=" * 60)
    
    print(f"\n{'Molecule':<12} {'Count':>8} {'MW (g/mol)':>12} {'Target wt%':>12} {'Actual wt%':>12}")
    print("-" * 60)
    
    for name, count in comp.counts.items():
        mw = comp.mw_used[name]
        target = comp.target_wt_pct[name]
        actual = comp.actual_wt_pct[name]
        print(f"{name:<12} {count:>8} {mw:>12.2f} {target:>12.2f} {actual:>12.2f}")
    
    print("-" * 60)
    print(f"{'Total':<12} {comp.total_molecules:>8}")
    print(f"\nMax wt% deviation: {comp.max_deviation_pct:.3f}%")
    print(f"Total mass: {comp.total_mass_g:.6e} g")
    print(f"Scale factor: {comp.scale_factor:.2f}")

    if comp.warnings:
        print("\nWarnings:")
        for category in sorted(comp.warnings.keys()):
            for msg in comp.warnings[category]:
                print(f"  - [{category}] {msg}")

    print(f"\n{'BOX SIZING':=^60}")
    requested_rho = box.requested_rho_g_cm3 if box.requested_rho_g_cm3 is not None else box.rho0_g_cm3
    if abs(box.rho0_g_cm3 - requested_rho) > 1e-12:
        print(
            f"Requested density: {requested_rho:.3f} g/cm³ "
            f"(shrinkage={box.polymerization_shrinkage_vol_frac:.1%})"
        )
        print(f"Effective density used: {box.rho0_g_cm3:.3f} g/cm³")
    else:
        print(f"Target density: {box.rho0_g_cm3:.3f} g/cm³")
    print(f"Box volume: {box.volume_nm3:.2f} nm³")
    print(f"Box length: {box.length_nm:.3f} nm")
    print("=" * 60 + "\n")


def _validate_unique_molecule_names(molecules: List[MoleculeSpec]) -> None:
    """Fail on duplicate molecule names to avoid silent dict overwrite."""
    seen = set()
    duplicates = set()
    for mol in molecules:
        if mol.name in seen:
            duplicates.add(mol.name)
        seen.add(mol.name)
    if duplicates:
        dup_list = sorted(duplicates)
        raise ValueError(f"Molecule names must be unique, duplicates: {dup_list}")


def _normalize_role(role: Optional[str]) -> Optional[str]:
    if role is None:
        return None
    role_norm = role.strip().lower()
    return role_norm if role_norm else None


def _add_warning(warning_messages: Dict[str, List[str]], key: str, message: str) -> None:
    warning_messages.setdefault(key, []).append(message)
    warnings.warn(message, RuntimeWarning, stacklevel=3)


def _evaluate_gel_guardrail(
    molecules: List[MoleculeSpec],
    mole_fractions: List[float],
    scale_factor: float,
    min_crosslinker_count: int,
    gel_policy: str,
    warning_messages: Dict[str, List[str]],
) -> Dict[str, Any]:
    crosslinker_entries = []
    for idx, mol in enumerate(molecules):
        if _normalize_role(mol.role) != "crosslinker":
            continue
        expected_count = mole_fractions[idx] * scale_factor
        crosslinker_entries.append(
            {
                "name": mol.name,
                "expected_count": expected_count,
            }
        )

    if not crosslinker_entries:
        return {}

    guardrail = {
        "policy": gel_policy,
        "gel_min_crosslinker_count": min_crosslinker_count,
        "crosslinkers": crosslinker_entries,
        "risk_detected": False,
    }

    if min_crosslinker_count <= 0 or gel_policy == "off":
        return guardrail

    risk_entries = []
    for entry in crosslinker_entries:
        expected = entry["expected_count"]
        if expected >= min_crosslinker_count:
            continue
        guardrail["risk_detected"] = True
        msg = (
            f"Gel-point risk for crosslinker '{entry['name']}': expected count "
            f"{expected:.2f} < gel_min_crosslinker_count={min_crosslinker_count}. "
            f"This can produce non-percolating microgels under PBC."
        )
        if gel_policy == "error":
            raise ValueError(msg)
        if gel_policy == "warn":
            _add_warning(warning_messages, "gel_point_risk", msg)
        risk_entries.append({"name": entry["name"], "expected_count": expected})

    if risk_entries:
        guardrail["risk_entries"] = risk_entries
    return guardrail


def _collect_molecule_metadata(molecules: List[MoleculeSpec]) -> Dict[str, Dict[str, Any]]:
    metadata: Dict[str, Dict[str, Any]] = {}
    for mol in molecules:
        entry: Dict[str, Any] = {}
        if mol.role is not None:
            entry["role"] = mol.role
            role_norm = _normalize_role(mol.role)
            if role_norm in RECOMMENDED_MOLECULE_ROLES:
                entry["role_normalized"] = role_norm
        if mol.mw_note:
            entry["mw_note"] = mol.mw_note
        if mol.pdi is not None:
            entry["pdi"] = mol.pdi
        if entry:
            metadata[mol.name] = entry
    return metadata


def _evaluate_polydispersity_warnings(
    molecules: List[MoleculeSpec],
    warning_messages: Dict[str, List[str]],
) -> None:
    for mol in molecules:
        role = _normalize_role(mol.role)
        if role not in {"monomer", "crosslinker"}:
            continue
        name_upper = mol.name.upper()
        is_common_polydisperse = any(pattern in name_upper for pattern in COMMON_POLYDISPERSE_PATTERNS)
        if not is_common_polydisperse:
            continue
        if mol.mw_note:
            continue
        msg = (
            f"Polydispersity approximation warning for '{mol.name}': fixed MW={mol.mw_g_mol} g/mol "
            "is being used without mw_note; provide mw_note (e.g., representative MW) and optional pdi."
        )
        _add_warning(warning_messages, "polydispersity_approximation", msg)


def _validate_shrinkage_vol_frac(shrinkage_vol_frac: float) -> float:
    value = float(shrinkage_vol_frac)
    if value < 0.0 or value >= 0.3:
        raise ValueError(
            "polymerization_shrinkage_vol_frac must satisfy 0 <= value < 0.3, "
            f"got {shrinkage_vol_frac}"
        )
    return value


def compute_effective_density(
    rho_g_cm3: float,
    polymerization_shrinkage_vol_frac: float = 0.0,
) -> float:
    """Compute shrinkage-adjusted effective packing density."""
    rho = float(rho_g_cm3)
    if rho <= 0:
        raise ValueError(f"rho_g_cm3 must be > 0, got {rho_g_cm3}")
    shrinkage = _validate_shrinkage_vol_frac(polymerization_shrinkage_vol_frac)
    return rho / (1.0 - shrinkage)


def get_density_candidates(config: CompositionConfig) -> List[float]:
    """
    Return deterministic effective-density candidates for callers to try.

    Order: target, min, max (duplicates removed while preserving order).
    """
    raw_candidates = [
        float(config.rho0_g_cm3),
        float(config.rho0_min_g_cm3),
        float(config.rho0_max_g_cm3),
    ]
    shrinkage = float(getattr(config, "polymerization_shrinkage_vol_frac", 0.0))
    candidates: List[float] = []
    for rho in raw_candidates:
        effective_rho = compute_effective_density(rho, shrinkage)
        if not any(abs(existing - effective_rho) < 1e-12 for existing in candidates):
            candidates.append(effective_rho)
    return candidates


if __name__ == "__main__":
    # Quick test
    comp, box = compute_full_composition()
    print_composition_summary(comp, box)
