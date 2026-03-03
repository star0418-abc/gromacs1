"""
PackmolStage: Runs PACKMOL to generate initial box configuration.

Implements:
- wt% → integer molecule count conversion
- Box sizing from mass and density
- PACKMOL input generation
- PACKMOL execution with retry and density floor
- Atomic publish to IN/systems/<SYSTEM_ID>/packmol/
- Unit validation (nm vs Å)
- Versioned publishing with content hash
- Per-system composition loading from YAML
"""

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Dict, Any, Tuple, List
import hashlib
import math
import os
import shutil
import uuid
from datetime import datetime

from .base import BaseStage

if TYPE_CHECKING:
    from ..context import PipelineContext
    from ..packmol_input import RetryAttemptInfo


# =============================================================================
# Constants
# =============================================================================

# Polymer-likeness detection
POLYMER_NAME_TOKENS = {
    "PEO", "PEG", "PVDF", "PAM", "PEGDMA", "OEGMA", "PMMA", "PVA",
    "POLYMER", "CHAIN", "NETWORK", "CROSSLINK", "MONOMER",
}
POLYMER_MW_THRESHOLD = 1000.0  # g/mol - above this, assume polymer-like

LI_NAME_TOKENS = {"LI", "LITHIUM", "LIPLUS"}
TFSI_NAME_TOKENS = {"TFSI", "FSI", "PF6", "BF4"}
SOLVENT_NAME_TOKENS = {"PC", "EC", "DMC", "EMC", "DEC", "G4", "GLYME", "THF", "GBL", "ACN", "DMF"}


# =============================================================================
# Helper Functions
# =============================================================================

def _is_polymer_like(mol_name: str, mw_g_mol: float) -> bool:
    """
    Detect if a molecule is polymer-like based on name or molecular weight.
    
    Args:
        mol_name: Molecule name
        mw_g_mol: Molecular weight in g/mol
        
    Returns:
        True if molecule appears polymer-like
    """
    if mw_g_mol >= POLYMER_MW_THRESHOLD:
        return True
    name_upper = mol_name.upper()
    return any(token in name_upper for token in POLYMER_NAME_TOKENS)


def _normalize_name(value: str) -> str:
    """Normalize names for resilient matching across files/CLI/YAML."""
    return "".join(ch for ch in value.upper() if ch.isalnum())


def _guess_element(atom_name: str, element_field: str = "") -> str:
    """Best-effort element guess from PDB atom name/element columns."""
    token = (element_field or "").strip().upper()
    if token:
        return token
    letters = "".join(ch for ch in atom_name.upper() if ch.isalpha())
    if not letters:
        return ""
    if letters.startswith("LI"):
        return "LI"
    return letters[0]


def parse_pdb_atoms(pdb_path: Path) -> List[Dict[str, Any]]:
    """
    Parse PDB atoms into lightweight records used by diagnostics and metrics.

    Returns:
        List of atom records with coordinates and basic identity fields.
    """
    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            atom_name = line[12:16].strip() if len(line) >= 16 else ""
            resname = line[17:20].strip() if len(line) >= 20 else ""
            chain_id = line[21:22].strip() if len(line) >= 22 else ""
            resid_raw = line[22:26].strip() if len(line) >= 26 else ""
            element = _guess_element(atom_name, line[76:78] if len(line) >= 78 else "")
            try:
                resid = int(resid_raw)
            except ValueError:
                resid = 0
            atoms.append(
                {
                    "x": x,
                    "y": y,
                    "z": z,
                    "atom_name": atom_name,
                    "resname": resname,
                    "chain_id": chain_id,
                    "resid": resid,
                    "element": element,
                }
            )
    if not atoms:
        raise ValueError(f"No atoms found in PDB file: {pdb_path}")
    return atoms


def _compute_extent_from_coords(coords: List[Tuple[float, float, float]]) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]:
    """Compute min/max/extent for (x, y, z) coordinates."""
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    min_xyz = (min(xs), min(ys), min(zs))
    max_xyz = (max(xs), max(ys), max(zs))
    extent_xyz = (
        max_xyz[0] - min_xyz[0],
        max_xyz[1] - min_xyz[1],
        max_xyz[2] - min_xyz[2],
    )
    return min_xyz, max_xyz, extent_xyz


def _percentile(values: List[float], q: float) -> float:
    """
    Linear percentile with deterministic interpolation.

    Args:
        values: Non-empty list of values
        q: Percentile in [0, 100]
    """
    if not values:
        raise ValueError("Cannot compute percentile of empty list")
    if len(values) == 1:
        return values[0]
    q_clamped = max(0.0, min(100.0, q))
    sorted_vals = sorted(values)
    pos = (len(sorted_vals) - 1) * q_clamped / 100.0
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def _compute_winsorized_extent(coords: List[Tuple[float, float, float]], lower_pct: float = 0.5, upper_pct: float = 99.5) -> Tuple[float, float, float]:
    """Compute axis extents from winsorized coordinates (cheap span-robust metric)."""
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    lo_x, hi_x = _percentile(xs, lower_pct), _percentile(xs, upper_pct)
    lo_y, hi_y = _percentile(ys, lower_pct), _percentile(ys, upper_pct)
    lo_z, hi_z = _percentile(zs, lower_pct), _percentile(zs, upper_pct)
    return (hi_x - lo_x, hi_y - lo_y, hi_z - lo_z)


def _sample_local_nn_distances(atoms: List[Dict[str, Any]], max_pairs: int = 4000) -> Dict[str, Any]:
    """
    Residue-local nearest-neighbor distance sampling (deterministic, bounded).

    Strategy:
    1) Group atoms by (chain_id, resname, resid).
    2) For a bounded subset of residues and atoms, compute nearest-neighbor distance
       inside the same residue.
    3) Use the median of these distances as local unit signal.
    """
    residue_map: Dict[Tuple[str, str, int], List[Tuple[float, float, float]]] = {}
    for atom in atoms:
        key = (
            str(atom.get("chain_id", "")),
            str(atom.get("resname", "")),
            int(atom.get("resid", 0)),
        )
        residue_map.setdefault(key, []).append((atom["x"], atom["y"], atom["z"]))

    distances: List[float] = []
    max_residues = 256
    max_atoms_per_residue = 24
    residue_keys = sorted(residue_map.keys(), key=lambda x: (x[0], x[1], x[2]))[:max_residues]
    residues_used = 0

    for key in residue_keys:
        coords = residue_map.get(key, [])
        if len(coords) < 2:
            continue
        residues_used += 1
        if len(coords) > max_atoms_per_residue:
            step = max(1, len(coords) // max_atoms_per_residue)
            sample_indices = list(range(0, len(coords), step))[:max_atoms_per_residue]
        else:
            sample_indices = list(range(len(coords)))

        for idx in sample_indices:
            x0, y0, z0 = coords[idx]
            best2 = float("inf")
            for j, (x1, y1, z1) in enumerate(coords):
                if j == idx:
                    continue
                dx = x0 - x1
                dy = y0 - y1
                dz = z0 - z1
                d2 = dx * dx + dy * dy + dz * dz
                if d2 < best2:
                    best2 = d2
            if best2 < float("inf") and best2 > 1e-12:
                distances.append(math.sqrt(best2))
            if len(distances) >= max_pairs:
                break
        if len(distances) >= max_pairs:
            break

    if not distances:
        return {
            "count": 0,
            "min": None,
            "median": None,
            "max": None,
            "units": "raw",
            "sampling": "none",
            "residues_used": residues_used,
        }
    sorted_d = sorted(distances)
    mid = len(sorted_d) // 2
    median = (
        sorted_d[mid]
        if len(sorted_d) % 2 == 1
        else 0.5 * (sorted_d[mid - 1] + sorted_d[mid])
    )
    return {
        "count": len(sorted_d),
        "min": sorted_d[0],
        "median": median,
        "max": sorted_d[-1],
        "units": "raw",
        "sampling": "residue_local_nearest_neighbor",
        "residues_used": residues_used,
    }


def _classify_extent_signal(ratio_extent: float) -> Tuple[Optional[str], str]:
    """Classify extent/box ratio as A, nm, or unknown."""
    ratio_a_min = 5.0
    ratio_a_max = 15.0
    ratio_nm_min = 0.3
    ratio_nm_max = 3.0
    if ratio_a_min <= ratio_extent <= ratio_a_max:
        return "A", f"winsor_extent/box={ratio_extent:.3f} in [{ratio_a_min}, {ratio_a_max}]"
    if ratio_nm_min <= ratio_extent <= ratio_nm_max:
        return "nm", f"winsor_extent/box={ratio_extent:.3f} in [{ratio_nm_min}, {ratio_nm_max}]"
    return None, f"winsor_extent/box={ratio_extent:.3f} outside heuristic windows"


def _classify_distance_signal(nn_stats: Dict[str, Any]) -> Tuple[Optional[str], str]:
    """Classify sampled local-distance median as A, nm, or unknown."""
    median = nn_stats.get("median")
    if median is None:
        return None, "no local-distance sample"
    # Typical bonded distances in raw units:
    # - raw=Å: ~1.0-1.6
    # - raw=nm: ~0.10-0.16
    if 0.70 <= median <= 2.20:
        return "A", f"nn_median={median:.3f} suggests angstrom-scale bonding"
    if 0.07 <= median <= 0.22:
        return "nm", f"nn_median={median:.3f} suggests nanometer-scale bonding"
    return None, f"nn_median={median:.3f} outside expected bonded windows"


class PackmolUnitInferenceError(ValueError):
    """Unit inference failure with attached diagnostics."""

    def __init__(self, message: str, diagnostics: Dict[str, Any]):
        super().__init__(message)
        self.diagnostics = diagnostics


class PackmolConversionError(RuntimeError):
    """Conversion failure with attached diagnostics."""

    def __init__(self, message: str, diagnostics: Dict[str, Any]):
        super().__init__(message)
        self.diagnostics = diagnostics


def infer_scale_factor(
    box_length_nm: float,
    raw_extent_xyz: Tuple[float, float, float],
    winsor_extent_xyz: Tuple[float, float, float],
    nn_distance_stats: Dict[str, Any],
    scale_override: Optional[float] = None,
    strict_mode: bool = True,
) -> tuple:
    """
    Infer PDB-to-nm scale factor using extent + local-distance signals.

    Extent uses a winsorized span to reduce false unit positives from chain/PBC spanning.
    Distance uses sampled local nearest-neighbor distances that are insensitive to global span.
    """
    if box_length_nm <= 0:
        raise ValueError(f"Invalid box_length_nm: {box_length_nm}")

    raw_extent_max = max(raw_extent_xyz)
    winsor_extent_max = max(winsor_extent_xyz)
    ratio_extent = winsor_extent_max / box_length_nm

    diagnostics: Dict[str, Any] = {
        "box_length_nm": box_length_nm,
        "raw_extent_units": "template_raw",
        "raw_extent_A": list(raw_extent_xyz),
        "winsor_extent_A": list(winsor_extent_xyz),
        "raw_extent_max_A": raw_extent_max,
        "winsor_extent_max_A": winsor_extent_max,
        "ratio_extent": ratio_extent,
        "sampled_nn_distance_stats": nn_distance_stats,
        "decision_source": None,
        "unit_inferred": None,
        "ambiguity_reason": None,
    }

    if scale_override is not None and scale_override > 0:
        unit_label = "Å" if scale_override < 0.5 else "nm"
        diagnostics["decision_source"] = "explicit_override"
        diagnostics["unit_inferred"] = "override"
        diagnostics["scale_override"] = scale_override
        return scale_override, unit_label, diagnostics

    extent_signal, extent_reason = _classify_extent_signal(ratio_extent)
    distance_signal, distance_reason = _classify_distance_signal(nn_distance_stats)
    diagnostics["extent_signal"] = extent_signal or "unknown"
    diagnostics["distance_signal"] = distance_signal or "unknown"
    diagnostics["extent_signal_reason"] = extent_reason
    diagnostics["distance_signal_reason"] = distance_reason

    chosen_signal: Optional[str] = None
    decision_source: Optional[str] = None
    ambiguity_reason: Optional[str] = None

    if extent_signal and distance_signal:
        if extent_signal == distance_signal:
            chosen_signal = extent_signal
            decision_source = "extent+distance"
        else:
            ambiguity_reason = (
                f"signal_disagreement: extent={extent_signal} ({extent_reason}); "
                f"distance={distance_signal} ({distance_reason})"
            )
    elif extent_signal:
        chosen_signal = extent_signal
        decision_source = "extent"
    elif distance_signal:
        chosen_signal = distance_signal
        decision_source = "distance"
    else:
        ambiguity_reason = (
            f"insufficient_signal: extent={extent_reason}; distance={distance_reason}"
        )

    if chosen_signal == "A":
        diagnostics["decision_source"] = decision_source
        diagnostics["unit_inferred"] = "A"
        return 0.1, "Å", diagnostics
    if chosen_signal == "nm":
        diagnostics["decision_source"] = decision_source
        diagnostics["unit_inferred"] = "nm"
        return 1.0, "nm", diagnostics

    diagnostics["unit_inferred"] = "ambiguous"
    diagnostics["ambiguity_reason"] = ambiguity_reason or "unknown_ambiguity"
    diagnostics["decision_source"] = "ambiguous_error"
    strict_note = (
        "strict_packmol_units=False does not permit ambiguous scale selection."
        if not strict_mode
        else "strict_packmol_units=True enforces fail-fast unit disambiguation."
    )
    raise PackmolUnitInferenceError(
        "Ambiguous PDB unit inference from winsorized extent and residue-local distance signals.\n"
        f"- Extent signal: {extent_reason}\n"
        f"- Distance signal: {distance_reason}\n"
        f"- Reason: {diagnostics['ambiguity_reason']}\n"
        f"- Policy: {strict_note}\n"
        "Set --packmol-pdb-scale explicitly (0.1 for Å→nm, 1.0 for nm→nm), "
        "or provide a pre-scaled input template.",
        diagnostics,
    )


class PackmolStage(BaseStage):
    """
    Stage 1: PACKMOL
    
    Reads:
    - IN/molecules/<NAME>/pdb/ - molecular structure templates
    - IN/systems/<SYSTEM_ID>/composition.yaml (if exists) or default recipe
    
    Writes:
    - OUT_GMX/<RUN_ID>/01_packmol/packmol.inp - PACKMOL input file
    - OUT_GMX/<RUN_ID>/01_packmol/initial.pdb - packed system
    - OUT_GMX/<RUN_ID>/01_packmol/initial.gro - converted GRO with box vectors
    
    Publishes (on success):
    - IN/systems/<SYSTEM_ID>/packmol/pdb/initial.pdb
    - IN/systems/<SYSTEM_ID>/packmol/gro/initial.gro
    - IN/systems/<SYSTEM_ID>/packmol/current.sha256
    """
    
    @property
    def name(self) -> str:
        return "packmol"
    
    @property
    def output_subdir(self) -> str:
        return "01_packmol"
    
    def run(self, ctx: "PipelineContext") -> bool:
        """
        Execute PACKMOL stage.
        
        1. Load composition config (per-system YAML or default)
        2. Compute box size from mass/density
        3. Generate PACKMOL input
        4. Run PACKMOL with retry and density floor
        5. Validate output (unit check)
        6. Convert PDB→GRO with explicit box vectors
        7. Publish result atomically with versioning
        """
        from ..composition import (
            wt_to_counts,
            compute_box_size,
            compute_effective_density,
            get_density_candidates,
            print_composition_summary,
        )
        from ..packmol_input import (
            write_packmol_input,
            run_packmol_with_retry,
            PackmolError,
            validate_pdb_box_extent,
            summarize_template_formats,
        )
        
        print(f"  Running PACKMOL stage...")
        print(f"  - Forcefield: {ctx.ff}")
        print(f"  - Charge method: {ctx.charge}")
        print(f"  - System: {ctx.system_id}")
        print(f"  - Strict publish: {ctx.strict_publish}")
        
        # Get paths
        output_dir = self.get_output_dir(ctx)
        molecules_dir = ctx.get_input_path("molecules")
        system_packmol_dir = ctx.get_input_path("systems", ctx.system_id, "packmol")
        
        print(f"  - Output dir: {output_dir}")
        print(f"  - Molecules dir: {molecules_dir}")
        
        # === Step 1: Load composition (per-system or default) ===
        print(f"\n  Loading composition...")
        try:
            config, composition_source, composition_extras = self._load_composition_config(ctx)
        except (FileNotFoundError, ValueError, TypeError) as e:
            print(f"\n  [ERROR] {e}")
            self.mark_failed(ctx, str(e))
            return False
        print(f"  - Composition source: {composition_source}")
        
        try:
            composition = wt_to_counts(config)
        except ValueError as e:
            msg = f"Invalid composition configuration: {e}"
            print(f"\n  [ERROR] {msg}")
            self.mark_failed(ctx, msg)
            return False
        
        # === Step 1b: Polymer percolation guardrails ===
        polymer_check_result = self._check_polymer_percolation(
            ctx, config, composition, composition_source
        )
        if not polymer_check_result:
            return False
        
        # === Step 2: Compute box size ===
        try:
            shrinkage_vol_frac = float(getattr(config, "polymerization_shrinkage_vol_frac", 0.0))
            target_density = compute_effective_density(config.rho0_g_cm3, shrinkage_vol_frac)
            rho_min_effective = compute_effective_density(config.rho0_min_g_cm3, shrinkage_vol_frac)
            rho_max_effective = compute_effective_density(config.rho0_max_g_cm3, shrinkage_vol_frac)
            density_candidates = get_density_candidates(config)
            box = compute_box_size(
                composition.total_mass_g,
                config.rho0_g_cm3,
                polymerization_shrinkage_vol_frac=shrinkage_vol_frac,
                density_candidates_g_cm3=density_candidates,
            )
        except ValueError as e:
            msg = f"Invalid density/shrinkage configuration: {e}"
            print(f"\n  [ERROR] {msg}")
            self.mark_failed(ctx, msg)
            return False
        target_box_length_nm = box.length_nm
        preassembly_settings = self._resolve_preassembly_settings(
            ctx=ctx,
            config=config,
            composition=composition,
            composition_extras=composition_extras,
        )
        
        # Print summary
        print_composition_summary(composition, box)
        if shrinkage_vol_frac > 0:
            print(
                "  - Shrinkage-aware density: "
                f"rho_input={config.rho0_g_cm3:.4f}, rho_effective={target_density:.4f} g/cm³, "
                f"shrinkage={shrinkage_vol_frac:.1%}"
            )
        if preassembly_settings["mode"] != "none":
            print(
                "  - Preassembly mode: "
                f"{preassembly_settings['mode']} (repair retries={preassembly_settings['repair_retry_count']})"
            )
        
        # === Step 3 & 4: Generate input and run PACKMOL ===
        if ctx.dry_run:
            print(f"\n  [DRY-RUN] Would generate PACKMOL input and execute")
            print(f"  [DRY-RUN] Skipping actual PACKMOL execution")
            
            # Validate templates exist in dry-run
            try:
                inp_path, templates = write_packmol_input(
                    composition=composition,
                    box=box,
                    output_dir=output_dir,
                    molecules_dir=molecules_dir,
                    placement_hints=self._build_packmol_placement_hints(
                        preassembly_settings,
                        tightening=0,
                    ),
                )
                template_format, _ = summarize_template_formats(templates)
                if template_format == "mixed":
                    print(f"  [DRY-RUN] FAILED: Mixed template formats (pdb + gro) are not supported")
                    self._record_manifest(
                        ctx, composition, box, config,
                        target_density=target_density,
                        achieved_density=target_density,
                        composition_source=composition_source,
                        dry_run=True,
                        preview_only=True,
                        publish_status="failed",
                    )
                    return False
                print(f"  [DRY-RUN] Generated: {inp_path.name}")
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=target_density,
                    composition_source=composition_source,
                    dry_run=True,
                    preview_only=True,
                )
                return True
            except FileNotFoundError as e:
                # In dry-run, missing templates is a FAILURE with actionable message
                print(f"  [DRY-RUN] FAILED: Missing templates")
                print(f"  [DRY-RUN] {e}")
                print(f"  [DRY-RUN] Please ensure all molecule templates exist before running.")
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=target_density,
                    composition_source=composition_source,
                    dry_run=True,
                    preview_only=True,
                    publish_status="failed",
                )
                return False
            except (PackmolError, ValueError) as e:
                print(f"  [DRY-RUN] FAILED: Invalid PACKMOL input settings")
                print(f"  [DRY-RUN] {e}")
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=target_density,
                    composition_source=composition_source,
                    dry_run=True,
                    preview_only=True,
                    publish_status="failed",
                )
                return False
        
        # Actual execution with density floor enforcement
        retry_attempts: List["RetryAttemptInfo"] = []
        retry_meta: Dict[str, Any] = {}
        preassembly_reports: List[Dict[str, Any]] = []
        preassembly_mode = preassembly_settings["mode"]
        requested_repairs = preassembly_settings["repair_retry_count"]
        max_repair_rounds = max(0, requested_repairs) if preassembly_mode != "none" else 0

        output_pdb: Optional[Path] = None
        templates: Dict[str, Tuple[Path, str]] = {}
        template_format = "unknown"
        gro_path: Optional[Path] = None
        final_box_nm = box.length_nm
        unit_diagnostics: Dict[str, Any] = {}
        li_contact_diag: Dict[str, Any] = {}
        retries_used = 0
        preassembly_tolerance_ang = 2.0
        if preassembly_mode != "none":
            inner_nm = float(preassembly_settings.get("li_inner_exclusion_nm", 0.19))
            preassembly_tolerance_ang = max(preassembly_tolerance_ang, inner_nm * 10.0)
            preassembly_settings["effective_packmol_tolerance_ang"] = preassembly_tolerance_ang
            if preassembly_tolerance_ang > 2.0 + 1e-9:
                print(
                    f"  [Preassembly Safety] Using PACKMOL tolerance={preassembly_tolerance_ang:.2f} Å "
                    f"(from inner exclusion {inner_nm:.3f} nm)."
                )

        for repair_round in range(max_repair_rounds + 1):
            if preassembly_mode != "none":
                print(f"\n  Preassembly attempt {repair_round + 1}/{max_repair_rounds + 1}...")
            placement_hints = self._build_packmol_placement_hints(
                preassembly_settings,
                tightening=repair_round,
            )
            packmol_seed = self._compute_packmol_seed(ctx, repair_round, preassembly_mode)
            max_packmol_retries = 3 if bool(getattr(ctx, "allow_density_reduction", False)) else 0
            try:
                output_pdb, final_box, templates, retries_used, retry_attempts, retry_meta = run_packmol_with_retry(
                    composition=composition,
                    initial_rho=target_density,
                    output_dir=output_dir,
                    molecules_dir=molecules_dir,
                    max_retries=max_packmol_retries,
                    rho_decrease_factor=0.95,
                    rho_min=rho_min_effective,
                    rho_max=rho_max_effective,
                    edge_margin_nm=ctx.packmol_edge_margin_nm,
                    box_margin_nm=ctx.packmol_edge_margin_nm,
                    tolerance=preassembly_tolerance_ang,
                    density_floor_fraction=ctx.density_floor_fraction,
                    allow_density_reduction=ctx.allow_density_reduction,
                    max_density_reduction_fraction=getattr(ctx, "packmol_max_density_reduction_fraction", 0.10),
                    pdb_scale_override=ctx.packmol_pdb_scale,
                    seed=packmol_seed,
                    placement_hints=placement_hints,
                    allow_random_seed=bool(getattr(ctx, "packmol_random_seed", False)),
                    seed_context=f"{ctx.run_id}|{ctx.system_id}|{preassembly_mode}|repair_round={repair_round}",
                )
            except PackmolError as e:
                print(f"\n  [ERROR] PACKMOL failed: {e}")
                retry_attempts = getattr(e, "retry_attempts", retry_attempts)
                retry_meta = getattr(e, "retry_meta", retry_meta)
                if retry_meta.get("template_fit"):
                    print("  [ERROR] Template-fit diagnostics were recorded in manifest (composition.retry_policy.template_fit).")
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=target_density,
                    retries_used=len(retry_attempts),
                    publish_status="failed",
                    composition_source=composition_source,
                    retry_attempts=retry_attempts,
                    retry_meta=retry_meta,
                    box_length_target_nm=target_box_length_nm,
                    box_length_final_nm=target_box_length_nm,
                    preassembly_info=(
                        {
                            "mode": preassembly_mode,
                            "settings": preassembly_settings,
                            "round_reports": preassembly_reports,
                            "repair_attempted": bool(max_repair_rounds > 0),
                            "repair_retries_used": max(0, len(preassembly_reports) - 1),
                        }
                        if preassembly_mode != "none"
                        else None
                    ),
                )
                ctx.log_command(
                    "packmol < packmol.inp",
                    self.name,
                    exit_code=1,
                    stderr=str(e),
                )
                self.mark_failed(ctx, str(e))
                return False
            except FileNotFoundError as e:
                print(f"\n  [ERROR] Missing templates: {e}")
                self.mark_failed(ctx, str(e))
                return False

            template_format, _ = summarize_template_formats(templates)
            if template_format == "mixed":
                msg = "Mixed template formats (pdb + gro) are not supported. Use a single format."
                print(f"\n  [ERROR] {msg}")
                self.mark_failed(ctx, msg)
                return False

            # Update box with final retry result
            box = final_box

            # === Step 5: Validate unit consistency ===
            print(f"\n  Validating output...")
            expected_unit = "nm" if template_format == "gro" else "A"
            valid, extent_msg = validate_pdb_box_extent(
                output_pdb,
                box.length_nm,
                expected_unit=expected_unit,
            )
            if not valid:
                print(f"  [ERROR] {extent_msg}")
                self.mark_failed(ctx, f"Unit validation failed: {extent_msg}")
                return False
            print(f"  - Unit validation: PASSED")

            # === Step 6: Convert PDB to GRO with explicit box and unit handling ===
            try:
                gro_path, final_box_nm, unit_diagnostics = self._convert_pdb_to_gro(
                    ctx, output_pdb, output_dir, box.length_nm, template_format
                )
            except (RuntimeError, ValueError) as exc:
                print(f"  [ERROR] {exc}")
                unit_diag_error = getattr(exc, "diagnostics", unit_diagnostics or {})
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=target_density,
                    retries_used=retries_used,
                    publish_status="failed",
                    composition_source=composition_source,
                    unit_diagnostics=unit_diag_error,
                    retry_attempts=retry_attempts,
                    retry_meta=retry_meta,
                    box_length_target_nm=target_box_length_nm,
                    box_length_final_nm=target_box_length_nm,
                )
                self.mark_failed(ctx, f"GRO conversion/unit inference failed: {exc}")
                return False

            if preassembly_mode == "none":
                break
            preassembly_eval = self._evaluate_preassembly_quality(
                pdb_path=output_pdb,
                templates=templates,
                preassembly_settings=preassembly_settings,
                unit_diagnostics=unit_diagnostics,
            )
            preassembly_eval["repair_round"] = repair_round
            preassembly_eval["packmol_seed"] = packmol_seed
            preassembly_eval["placement_hints"] = placement_hints
            preassembly_eval["retries_used"] = retries_used
            preassembly_reports.append(preassembly_eval)
            if preassembly_eval.get("passes_thresholds", True):
                print("  - Preassembly quality check: PASSED")
                break
            if repair_round < max_repair_rounds:
                print(
                    "  [WARNING] Preassembly quality thresholds not met; "
                    "retrying with adjusted seed/constraints."
                )
                continue
            print("  [WARNING] Preassembly quality thresholds not met after retries; proceeding with diagnostics.")
            break

        # Log PACKMOL command after accepted attempt
        ctx.log_command(
            "packmol < packmol.inp",
            self.name,
            exit_code=0,
            stdout=(
                f"Generated {output_pdb.name if output_pdb else 'initial.pdb'} after {retries_used} retry iterations; "
                f"preassembly_rounds={len(preassembly_reports)}"
            ),
        )

        if output_pdb is None:
            self.mark_failed(ctx, "PACKMOL produced no output PDB.")
            return False

        unit_diagnostics["strict_packmol_units"] = bool(getattr(ctx, "strict_packmol_units", True))
        if retry_meta:
            unit_diagnostics["retry_policy"] = retry_meta

        li_contact_diag = self._evaluate_li_contact_safety(
            pdb_path=output_pdb,
            unit_diagnostics=unit_diagnostics,
            min_contact_nm=float(getattr(ctx, "packmol_min_contact_nm", 0.18)),
            sample_cap=int(preassembly_settings.get("sample_cap", 4000)),
        )
        if li_contact_diag:
            li_contact_diag["enforcement_scope"] = "preassembly_bias_modes_only"
            li_contact_diag["enforced"] = preassembly_mode != "none"
            unit_diagnostics["li_heavy_contact"] = li_contact_diag
            li_min = li_contact_diag.get("min_distance_nm")
            if (
                preassembly_mode != "none"
                and li_min is not None
                and li_min < float(getattr(ctx, "packmol_min_contact_nm", 0.18))
            ):
                msg = (
                    f"Dangerous Li-heavy close contact detected after PACKMOL: "
                    f"min={li_min:.4f} nm < threshold {float(getattr(ctx, 'packmol_min_contact_nm', 0.18)):.4f} nm. "
                    "Adjust preassembly inner exclusion/tolerance or disable bias mode."
                )
                print(f"\n  [ERROR] {msg}")
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=target_density,
                    retries_used=retries_used,
                    publish_status="failed",
                    composition_source=composition_source,
                    unit_diagnostics=unit_diagnostics,
                    retry_attempts=retry_attempts,
                    retry_meta=retry_meta,
                    box_length_target_nm=target_box_length_nm,
                    box_length_final_nm=final_box_nm,
                    preassembly_info=(
                        {
                            "mode": preassembly_mode,
                            "settings": preassembly_settings,
                            "round_reports": preassembly_reports,
                            "final_report": preassembly_reports[-1] if preassembly_reports else None,
                            "repair_attempted": bool(max_repair_rounds > 0),
                            "repair_retries_used": max(0, len(preassembly_reports) - 1),
                            "contact_safety": li_contact_diag,
                        }
                        if preassembly_mode != "none"
                        else None
                    ),
                )
                self.mark_failed(ctx, msg)
                return False
            if (
                preassembly_mode == "none"
                and li_min is not None
                and li_min < float(getattr(ctx, "packmol_min_contact_nm", 0.18))
            ):
                print(
                    "  [WARNING] Dangerous Li-heavy close contact detected, but hard-fail "
                    "is preassembly-bias-only by design (mode=none)."
                )

        # === Recompute achieved density from FINAL box ===
        unit_diagnostics["packmol_box_nm"] = box.length_nm
        unit_diagnostics["final_box_nm"] = final_box_nm
        final_volume_nm3 = final_box_nm ** 3
        final_volume_cm3 = final_volume_nm3 * 1e-21
        achieved_density = composition.total_mass_g / final_volume_cm3
        unit_diagnostics["achieved_density_g_cm3"] = achieved_density

        box = type(box)(
            length_nm=final_box_nm,
            rho0_g_cm3=achieved_density,
            total_mass_g=composition.total_mass_g,
            volume_nm3=final_volume_nm3,
            requested_rho_g_cm3=config.rho0_g_cm3,
            polymerization_shrinkage_vol_frac=shrinkage_vol_frac,
            rho_effective_g_cm3=target_density,
            density_candidates_g_cm3=density_candidates,
        )

        if retries_used > 0 and abs(achieved_density - target_density) > 1e-9:
            print(
                f"\n  [INFO] Density adjusted from {target_density:.4f} to "
                f"{achieved_density:.4f} g/cm³ after {retries_used} retries"
            )

        # === Density floor enforcement (single definition) ===
        density_ratio = achieved_density / target_density if target_density > 0 else None
        relative_floor_density = ctx.density_floor_fraction * target_density
        absolute_floor_density = rho_min_effective
        effective_floor_density = max(relative_floor_density, absolute_floor_density)
        unit_diagnostics["density_floor"] = {
            "relative_floor_density_g_cm3": relative_floor_density,
            "absolute_floor_density_g_cm3": absolute_floor_density,
            "effective_floor_density_g_cm3": effective_floor_density,
            "achieved_density_ratio": density_ratio,
        }
        if achieved_density < effective_floor_density - 1e-12:
            if getattr(ctx, "packmol_allow_density_floor_violation", False):
                print(
                    "\n  [WARNING] OVERRIDE ACTIVE: accepting density below floor. "
                    "This can produce under-packed/porous initial states."
                )
                unit_diagnostics["density_floor_override"] = True
            else:
                msg = (
                    f"Achieved density {achieved_density:.4f} g/cm³ is below required floor "
                    f"{effective_floor_density:.4f} g/cm³ "
                    f"(relative floor={relative_floor_density:.4f}, absolute floor={absolute_floor_density:.4f})."
                )
                print(f"\n  [ERROR] {msg}")
                self._record_manifest(
                    ctx, composition, box, config,
                    target_density=target_density,
                    achieved_density=achieved_density,
                    retries_used=retries_used,
                    publish_status="failed",
                    composition_source=composition_source,
                    unit_diagnostics=unit_diagnostics,
                    retry_attempts=retry_attempts,
                    retry_meta=retry_meta,
                    box_length_target_nm=target_box_length_nm,
                    box_length_final_nm=final_box_nm,
                )
                self.mark_failed(ctx, msg)
                return False

        # Prevent publishing mixed freshness states (new PDB + stale/missing GRO).
        gro_publish_ready = bool(
            gro_path
            and gro_path.exists()
            and unit_diagnostics.get("gro_conversion_status") == "success"
        )
        unit_diagnostics["gro_publish_ready"] = gro_publish_ready
        unit_diagnostics["strict_gro_conversion"] = bool(
            getattr(ctx, "strict_gro_conversion", False)
        )
        if not gro_publish_ready:
            msg = (
                "GRO conversion did not produce a fresh publishable initial.gro; "
                "skipping IN publish to prevent PDB/GRO drift."
            )
            print(f"\n  [ERROR] {msg}")
            self._record_manifest(
                ctx, composition, box, config,
                target_density=target_density,
                achieved_density=achieved_density,
                retries_used=retries_used,
                publish_status="incomplete_missing_gro",
                packmol_publish_ok=False,
                composition_source=composition_source,
                unit_diagnostics=unit_diagnostics,
                retry_attempts=retry_attempts,
                retry_meta=retry_meta,
                box_length_target_nm=target_box_length_nm,
                box_length_final_nm=final_box_nm,
                publish_details={
                    "pdb_published": False,
                    "gro_published": False,
                    "both_published": False,
                    "gro_conversion_status": unit_diagnostics.get("gro_conversion_status"),
                    "publish_guard": "blocked_missing_gro",
                },
                preassembly_info=(
                    {
                        "mode": preassembly_mode,
                        "settings": preassembly_settings,
                        "round_reports": preassembly_reports,
                        "final_report": preassembly_reports[-1] if preassembly_reports else None,
                        "repair_attempted": bool(max_repair_rounds > 0),
                        "repair_retries_used": max(0, len(preassembly_reports) - 1),
                        "contact_safety": li_contact_diag if li_contact_diag else None,
                    }
                    if preassembly_mode != "none"
                    else None
                ),
            )
            self.mark_failed(ctx, msg)
            return False
        
        # === Step 7: Atomic publish with versioning ===
        print(f"\n  Publishing output to IN/systems/{ctx.system_id}/packmol/...")
        
        published_files = []
        publish_errors = []
        publish_status = "failed"
        published_pdb = False
        published_gro = False
        
        try:
            # Archive previous versions before overwriting
            self._archive_previous_versions(system_packmol_dir)
            
            # Publish GRO first so a missing/failed GRO never leaves new PDB + old GRO drift.
            dest_gro = system_packmol_dir / "gro" / "initial.gro"
            gro_hash = self._atomic_publish(gro_path, dest_gro)
            published_files.append({"path": str(dest_gro), "sha256": gro_hash})
            published_gro = True
            print(f"  - Published: {dest_gro}")

            # Publish PDB second once GRO update is guaranteed.
            dest_pdb = system_packmol_dir / "pdb" / "initial.pdb"
            pdb_hash = self._atomic_publish(output_pdb, dest_pdb)
            published_files.append({"path": str(dest_pdb), "sha256": pdb_hash})
            published_pdb = True
            print(f"  - Published: {dest_pdb}")
            
            # Write current.sha256 for traceability
            self._write_hash_manifest(system_packmol_dir, published_files)
            
            publish_status = "complete" if (published_pdb and published_gro) else "partial"
                
        except Exception as e:
            publish_errors.append(str(e))
            print(f"  [ERROR] Publish failed: {e}")
            
            if ctx.strict_publish:
                self.mark_failed(ctx, f"Publish failed (strict mode): {e}")
                return False
            else:
                print(
                    "  [WARNING] Continuing despite publish failure (non-strict mode). "
                    "packmol_publish_ok=false; downstream stages may read stale IN/packmol assets."
                )
                publish_status = "partial" if published_files else "failed"

        publish_details = {
            "pdb_published": bool(published_pdb),
            "gro_published": bool(published_gro),
            "both_published": bool(published_pdb and published_gro),
            "gro_conversion_status": unit_diagnostics.get("gro_conversion_status"),
        }
        
        # === Step 8: Record in manifest ===
        self._record_manifest(
            ctx, composition, box, config,
            target_density=target_density,
            achieved_density=achieved_density,
            retries_used=retries_used,
            published_files=published_files,
            publish_status=publish_status,
            packmol_publish_ok=(published_pdb and published_gro),
            composition_source=composition_source,
            unit_diagnostics=unit_diagnostics,
            retry_attempts=retry_attempts,  # NEW: Include all attempt diagnostics
            retry_meta=retry_meta,
            box_length_target_nm=target_box_length_nm,
            box_length_final_nm=final_box_nm,
            publish_details=publish_details,
            preassembly_info=(
                {
                    "mode": preassembly_mode,
                    "settings": preassembly_settings,
                    "round_reports": preassembly_reports,
                    "final_report": preassembly_reports[-1] if preassembly_reports else None,
                    "repair_attempted": bool(max_repair_rounds > 0),
                    "repair_retries_used": max(0, len(preassembly_reports) - 1),
                    "contact_safety": li_contact_diag if li_contact_diag else None,
                }
                if preassembly_mode != "none"
                else None
            ),
        )

        if not (published_pdb and published_gro):
            msg = (
                f"PACKMOL publish did not update both PDB and GRO (status={publish_status}). "
                "Stage marked failed to avoid false-success state."
            )
            self.mark_failed(ctx, msg)
            return False
        
        print(f"\n  PACKMOL stage completed successfully!")
        return True
    
    def _load_composition_config(
        self,
        ctx: "PipelineContext",
    ) -> Tuple[Any, str, Dict[str, Any]]:
        """
        Load composition config from per-system YAML or fallback to default.
        
        Returns:
            Tuple of (CompositionConfig, source_description, extra_settings)
        """
        from ..composition import get_default_recipe, CompositionConfig, MoleculeSpec
        import yaml
        
        # Check for per-system composition.yaml
        system_composition_path = ctx.get_input_path(
            "systems", ctx.system_id, "composition.yaml"
        )
        
        if system_composition_path.exists():
            print(f"  - Found per-system composition: {system_composition_path}")
            
            # Compute hash for provenance
            content_hash = self._compute_file_hash(system_composition_path)
            
            # Load YAML
            try:
                with open(system_composition_path, "r", encoding="utf-8") as f:
                    data = yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(
                    f"Malformed YAML in composition file '{system_composition_path}': {e}"
                ) from e

            data = self._validate_composition_yaml_schema(data, system_composition_path)
            
            # Parse into CompositionConfig
            molecules = []
            try:
                for mol_data in data["molecules"]:
                    molecules.append(MoleculeSpec(
                        name=mol_data["name"],
                        mw_g_mol=mol_data["mw_g_mol"],
                        wt_pct=mol_data["wt_pct"],
                        min_count=mol_data.get("min_count", 1),
                        role=mol_data.get("role"),
                        mw_note=mol_data.get("mw_note"),
                        pdi=mol_data.get("pdi"),
                    ))

                config = CompositionConfig(
                    molecules=molecules,
                    anchor_name=data.get("anchor_name"),
                    anchor_count=data.get("anchor_count"),
                    target_total_molecules=data.get("target_total_molecules", 500),
                    wt_deviation_threshold=data.get("wt_deviation_threshold", 2.0),
                    rho0_g_cm3=data.get("rho0_g_cm3", 1.0),
                    rho0_min_g_cm3=data.get("rho0_min_g_cm3", 0.95),
                    rho0_max_g_cm3=data.get("rho0_max_g_cm3", 1.05),
                    gel_min_crosslinker_count=data.get("gel_min_crosslinker_count", 0),
                    gel_policy=data.get("gel_policy"),
                    auto_scale_to_gel_min_crosslinker=data.get(
                        "auto_scale_to_gel_min_crosslinker",
                        False,
                    ),
                    polymerization_shrinkage_vol_frac=data.get(
                        "polymerization_shrinkage_vol_frac",
                        0.0,
                    ),
                )
            except (TypeError, ValueError) as e:
                raise ValueError(
                    f"Invalid composition schema/values in '{system_composition_path}': {e}"
                ) from e

            preassembly_cfg = data.get("packmol_preassembly", {})
            if not isinstance(preassembly_cfg, dict):
                preassembly_cfg = {}
            if "packmol_preassembly_mode" in data and "mode" not in preassembly_cfg:
                preassembly_cfg["mode"] = data.get("packmol_preassembly_mode")
            
            # Record in manifest
            if ctx.manifest:
                ctx.manifest.add_input_asset(
                    system_composition_path,
                    content_hash,
                    category="composition_config",
                )
            
            return (
                config,
                f"yaml:{system_composition_path}:sha256:{content_hash[:12]}",
                {"packmol_preassembly": preassembly_cfg},
            )
        
        else:
            # Use default recipe (demo/fallback mode)
            print(f"  - composition.yaml not found: IN/systems/{ctx.system_id}/composition.yaml")
            
            # GUARDRAIL: Default recipe is only allowed in dry-run or with explicit flag
            if not ctx.dry_run and not ctx.allow_default_recipe:
                msg = (
                    f"Missing composition.yaml for system '{ctx.system_id}'.\n"
                    f"  Expected: IN/systems/{ctx.system_id}/composition.yaml\n"
                    f"  For production runs, create a composition.yaml with explicit molecule counts.\n"
                    f"  For demo/testing only, use --allow-default-recipe to permit fallback."
                )
                raise FileNotFoundError(msg)
            
            if ctx.allow_default_recipe:
                print(f"  [WARNING] Using demo default recipe (--allow-default-recipe set)")
            else:
                print(f"  [DRY-RUN] Using default recipe for template generation")
            
            return (
                get_default_recipe(),
                "default_recipe:get_default_recipe()",
                {"packmol_preassembly": {}},
            )

    @staticmethod
    def _is_numeric(value: Any) -> bool:
        """Return True for int/float-like scalars except booleans."""
        return isinstance(value, (int, float)) and not isinstance(value, bool)

    def _validate_composition_yaml_schema(self, data: Any, source_path: Path) -> Dict[str, Any]:
        """Validate composition.yaml structure before constructing CompositionConfig."""
        if not isinstance(data, dict):
            raise ValueError(
                f"composition.yaml at '{source_path}' must be a mapping with a 'molecules' list."
            )

        molecules = data.get("molecules")
        if not isinstance(molecules, list) or not molecules:
            raise ValueError(
                f"composition.yaml at '{source_path}' must define a non-empty 'molecules' list."
            )

        for idx, mol in enumerate(molecules):
            prefix = f"molecules[{idx}]"
            if not isinstance(mol, dict):
                raise ValueError(f"{prefix} must be a mapping.")
            for key in ("name", "mw_g_mol", "wt_pct"):
                if key not in mol:
                    raise ValueError(f"{prefix} is missing required key '{key}'.")
            if not isinstance(mol.get("name"), str) or not mol.get("name", "").strip():
                raise ValueError(f"{prefix}.name must be a non-empty string.")
            if not self._is_numeric(mol.get("mw_g_mol")):
                raise ValueError(f"{prefix}.mw_g_mol must be numeric.")
            if not self._is_numeric(mol.get("wt_pct")):
                raise ValueError(f"{prefix}.wt_pct must be numeric.")
            if "min_count" in mol:
                min_count = mol.get("min_count")
                if not isinstance(min_count, int) or isinstance(min_count, bool):
                    raise ValueError(f"{prefix}.min_count must be an integer.")
                if min_count < 0:
                    raise ValueError(f"{prefix}.min_count must be >= 0.")
            for opt_key in ("role", "mw_note"):
                if opt_key in mol and mol.get(opt_key) is not None and not isinstance(mol.get(opt_key), str):
                    raise ValueError(f"{prefix}.{opt_key} must be a string when provided.")
            if "pdi" in mol and mol.get("pdi") is not None:
                pdi = mol.get("pdi")
                if not self._is_numeric(pdi):
                    raise ValueError(f"{prefix}.pdi must be numeric when provided.")
                if float(pdi) <= 0:
                    raise ValueError(f"{prefix}.pdi must be > 0 when provided.")

        numeric_fields = (
            "target_total_molecules",
            "wt_deviation_threshold",
            "rho0_g_cm3",
            "rho0_min_g_cm3",
            "rho0_max_g_cm3",
            "gel_min_crosslinker_count",
            "anchor_count",
            "polymerization_shrinkage_vol_frac",
        )
        for field in numeric_fields:
            if field in data and data.get(field) is not None and not self._is_numeric(data.get(field)):
                raise ValueError(f"{field} must be numeric when provided.")

        if "anchor_name" in data and data.get("anchor_name") is not None:
            if not isinstance(data.get("anchor_name"), str) or not data.get("anchor_name", "").strip():
                raise ValueError("anchor_name must be a non-empty string when provided.")

        if "gel_policy" in data and data.get("gel_policy") is not None and not isinstance(data.get("gel_policy"), str):
            raise ValueError("gel_policy must be a string when provided.")

        if (
            "auto_scale_to_gel_min_crosslinker" in data
            and data.get("auto_scale_to_gel_min_crosslinker") is not None
            and not isinstance(data.get("auto_scale_to_gel_min_crosslinker"), bool)
        ):
            raise ValueError("auto_scale_to_gel_min_crosslinker must be boolean when provided.")

        preassembly_cfg = data.get("packmol_preassembly")
        if preassembly_cfg is not None and not isinstance(preassembly_cfg, dict):
            raise ValueError("packmol_preassembly must be a mapping when provided.")

        return data

    def _resolve_preassembly_settings(
        self,
        ctx: "PipelineContext",
        config: Any,
        composition: Any,
        composition_extras: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Resolve preassembly configuration from ctx + composition YAML."""
        yaml_cfg = (composition_extras or {}).get("packmol_preassembly", {}) or {}
        if not isinstance(yaml_cfg, dict):
            yaml_cfg = {}

        valid_modes = {"none", "li_near_polymer", "li_solvent_shell"}
        ctx_mode = str(getattr(ctx, "packmol_preassembly_mode", "none") or "none").strip().lower()
        yaml_mode = str(yaml_cfg.get("mode", "none") or "none").strip().lower()
        mode = ctx_mode
        if mode == "none" and yaml_mode in valid_modes:
            mode = yaml_mode
        if mode not in valid_modes:
            mode = "none"

        polymer_names: List[str] = []
        li_names: List[str] = []
        tfsi_names: List[str] = []
        solvent_names: List[str] = []
        for mol in config.molecules:
            name = mol.name
            norm = _normalize_name(name)
            if _is_polymer_like(name, mol.mw_g_mol):
                polymer_names.append(name)
                continue
            is_tfsi = any(token in norm for token in TFSI_NAME_TOKENS)
            is_li = (norm in LI_NAME_TOKENS) or (norm.startswith("LI") and not is_tfsi)
            if is_li:
                li_names.append(name)
                continue
            if is_tfsi:
                tfsi_names.append(name)
                continue
            if any(token in norm for token in SOLVENT_NAME_TOKENS) or mol.mw_g_mol <= 260.0:
                solvent_names.append(name)

        # YAML overrides ctx defaults when ctx remains at default values.
        def _pick(
            field: str,
            default_value: Any,
        ) -> Any:
            ctx_val = getattr(ctx, field, default_value)
            yaml_val = yaml_cfg.get(field)
            if yaml_val is not None and ctx_val == default_value:
                return yaml_val
            return ctx_val

        settings: Dict[str, Any] = {
            "mode": mode,
            "li_names": sorted(li_names),
            "polymer_names": sorted(polymer_names),
            "solvent_names": sorted(solvent_names),
            "tfsi_names": sorted(tfsi_names),
            "li_fraction": float(_pick("packmol_preassembly_li_fraction", 1.0)),
            "repair_retry_count": int(_pick("packmol_preassembly_retry_count", 2)),
            "sample_cap": int(_pick("packmol_preassembly_sample_cap", 4000)),
            "li_inner_exclusion_nm": float(_pick("packmol_li_inner_exclusion_nm", 0.19)),
            "li_polymer_cutoff_nm": float(_pick("packmol_li_polymer_cutoff_nm", 0.35)),
            "li_solvent_cutoff_nm": float(_pick("packmol_li_solvent_cutoff_nm", 0.35)),
            "li_tfsi_cutoff_nm": float(_pick("packmol_li_tfsi_cutoff_nm", 0.40)),
            "target_li_polymer_fraction": float(_pick("packmol_target_li_polymer_fraction", 0.40)),
            "target_li_solvent_fraction": float(_pick("packmol_target_li_solvent_fraction", 0.50)),
            "max_li_tfsi_close_fraction": float(_pick("packmol_max_li_tfsi_close_fraction", 0.30)),
        }
        settings["li_fraction"] = max(0.0, min(1.0, settings["li_fraction"]))
        settings["repair_retry_count"] = max(0, settings["repair_retry_count"])
        settings["sample_cap"] = max(100, settings["sample_cap"])
        settings["li_inner_exclusion_nm"] = max(0.05, settings["li_inner_exclusion_nm"])
        return settings

    def _build_packmol_placement_hints(
        self,
        preassembly_settings: Dict[str, Any],
        tightening: int,
    ) -> Dict[str, Any]:
        """Build hint payload consumed by packmol_input.generate_packmol_input()."""
        mode = preassembly_settings.get("mode", "none")
        if mode == "none":
            return {"mode": "none"}
        return {
            "mode": mode,
            "li_names": list(preassembly_settings.get("li_names", [])),
            "polymer_names": list(preassembly_settings.get("polymer_names", [])),
            "solvent_names": list(preassembly_settings.get("solvent_names", [])),
            "tfsi_names": list(preassembly_settings.get("tfsi_names", [])),
            "li_fraction": float(preassembly_settings.get("li_fraction", 1.0)),
            "li_inner_exclusion_nm": float(preassembly_settings.get("li_inner_exclusion_nm", 0.19)),
            "tightening": int(tightening),
        }

    def _compute_packmol_seed(self, ctx: "PipelineContext", repair_round: int, mode: str) -> int:
        """Deterministic PACKMOL seed unless explicit random mode is requested."""
        if bool(getattr(ctx, "packmol_random_seed", False)):
            return -1
        seed_key = f"{ctx.run_id}|{ctx.system_id}|{mode}|{repair_round}"
        digest = hashlib.sha256(seed_key.encode("utf-8")).hexdigest()
        # PACKMOL accepts integer seed; keep within signed 31-bit range.
        return 1 + (int(digest[:8], 16) % 2_000_000_000)

    def _template_resnames(self, template_path: Path, template_fmt: str) -> List[str]:
        """Extract resnames from template files (PDB or GRO)."""
        resnames = set()
        try:
            with open(template_path, "r") as handle:
                if template_fmt == "gro":
                    lines = handle.readlines()
                    if len(lines) < 3:
                        return []
                    try:
                        atom_count = int(lines[1].strip())
                    except (ValueError, IndexError):
                        atom_count = max(0, len(lines) - 3)
                    atom_lines = lines[2:2 + atom_count]
                    for line in atom_lines:
                        if len(line) < 15:
                            continue
                        resname = line[5:10].strip()
                        if resname:
                            resnames.add(_normalize_name(resname))
                else:
                    for line in handle:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            resname = line[17:20].strip() if len(line) >= 20 else ""
                            if resname:
                                resnames.add(_normalize_name(resname))
        except OSError:
            return []
        return sorted(r for r in resnames if r)

    @staticmethod
    def _downsample_coords(coords: List[Tuple[float, float, float]], cap: int) -> List[Tuple[float, float, float]]:
        """Deterministic downsampling to cap list size while preserving spread."""
        if len(coords) <= cap:
            return coords
        step = max(1, len(coords) // cap)
        sampled = coords[::step]
        if len(sampled) > cap:
            sampled = sampled[:cap]
        return sampled

    def _fraction_within_cutoff(
        self,
        query_coords: List[Tuple[float, float, float]],
        target_coords: List[Tuple[float, float, float]],
        cutoff_nm: float,
        sample_cap: int,
    ) -> Dict[str, Any]:
        """Compute fraction of query points within cutoff of any target point."""
        if not query_coords:
            return {"fraction": None, "count_query": 0, "count_target": len(target_coords), "distance_min_nm": None, "distance_median_nm": None, "distance_max_nm": None}
        if not target_coords:
            return {"fraction": None, "count_query": len(query_coords), "count_target": 0, "distance_min_nm": None, "distance_median_nm": None, "distance_max_nm": None}
        q = self._downsample_coords(query_coords, cap=sample_cap)
        t = self._downsample_coords(target_coords, cap=sample_cap)
        cutoff2 = cutoff_nm * cutoff_nm
        within = 0
        min_distances: List[float] = []
        for qx, qy, qz in q:
            best2 = float("inf")
            for tx, ty, tz in t:
                dx = qx - tx
                dy = qy - ty
                dz = qz - tz
                d2 = dx * dx + dy * dy + dz * dz
                if d2 < best2:
                    best2 = d2
            if best2 <= cutoff2:
                within += 1
            min_distances.append(math.sqrt(best2))
        sorted_d = sorted(min_distances)
        mid = len(sorted_d) // 2
        median = sorted_d[mid] if len(sorted_d) % 2 == 1 else 0.5 * (sorted_d[mid - 1] + sorted_d[mid])
        return {
            "fraction": within / len(q),
            "count_query": len(q),
            "count_target": len(t),
            "distance_min_nm": sorted_d[0],
            "distance_median_nm": median,
            "distance_max_nm": sorted_d[-1],
        }

    def _evaluate_preassembly_quality(
        self,
        pdb_path: Path,
        templates: Dict[str, Tuple[Path, str]],
        preassembly_settings: Dict[str, Any],
        unit_diagnostics: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Evaluate Li-centered preassembly quality metrics on packed PDB output.

        Metrics:
        - fraction(Li within r of polymer O)
        - fraction(Li within r of solvent O)
        - fraction(Li within r of TFSI atoms) as close-pair proxy
        """
        mode = preassembly_settings.get("mode", "none")
        if mode == "none":
            return {"mode": "none", "status": "disabled", "passes_thresholds": True, "metrics": {}}

        atoms = parse_pdb_atoms(pdb_path)
        scale_factor = unit_diagnostics.get("scale_factor_used")
        if scale_factor is None:
            template_format = str(unit_diagnostics.get("template_format", "pdb"))
            scale_factor = 1.0 if template_format == "gro" else 0.1
        sample_cap = int(preassembly_settings.get("sample_cap", 4000))

        polymer_resnames: set = set()
        solvent_resnames: set = set()
        tfsi_resnames: set = set()
        for mol_name, (template_path, template_fmt) in templates.items():
            resnames = set(self._template_resnames(template_path, template_fmt))
            norm = _normalize_name(mol_name)
            if mol_name in preassembly_settings.get("polymer_names", []):
                polymer_resnames.update(resnames)
            if mol_name in preassembly_settings.get("solvent_names", []):
                solvent_resnames.update(resnames)
            if mol_name in preassembly_settings.get("tfsi_names", []):
                tfsi_resnames.update(resnames)
            if not resnames:
                if any(token in norm for token in TFSI_NAME_TOKENS):
                    tfsi_resnames.add(norm[:4])

        li_coords: List[Tuple[float, float, float]] = []
        polymer_o_coords: List[Tuple[float, float, float]] = []
        solvent_o_coords: List[Tuple[float, float, float]] = []
        tfsi_coords: List[Tuple[float, float, float]] = []
        fallback_methods: Dict[str, str] = {
            "polymer_oxygen": "template_resname_match",
            "solvent_oxygen": "template_resname_match",
        }
        degraded_reasons: List[str] = []

        for atom in atoms:
            x_nm = atom["x"] * scale_factor
            y_nm = atom["y"] * scale_factor
            z_nm = atom["z"] * scale_factor
            coord = (x_nm, y_nm, z_nm)
            res_norm = _normalize_name(atom.get("resname", ""))
            atom_norm = _normalize_name(atom.get("atom_name", ""))
            element = str(atom.get("element", "")).upper()

            is_li = element == "LI" or atom_norm.startswith("LI")
            if is_li:
                li_coords.append(coord)
            if res_norm in tfsi_resnames:
                tfsi_coords.append(coord)
            if element == "O":
                if res_norm in polymer_resnames:
                    polymer_o_coords.append(coord)
                if res_norm in solvent_resnames:
                    solvent_o_coords.append(coord)

        # Fallbacks when template resnames are sparse (e.g., generic templates).
        oxygen_pool_non_tfsi: List[Tuple[float, float, float, str]] = []
        for atom in atoms:
            if str(atom.get("element", "")).upper() != "O":
                continue
            res_norm = _normalize_name(atom.get("resname", ""))
            if res_norm in tfsi_resnames:
                continue
            oxygen_pool_non_tfsi.append(
                (atom["x"] * scale_factor, atom["y"] * scale_factor, atom["z"] * scale_factor, res_norm)
            )

        if not polymer_o_coords and preassembly_settings.get("polymer_names"):
            if solvent_resnames:
                polymer_o_coords = [
                    (x, y, z)
                    for x, y, z, res_norm in oxygen_pool_non_tfsi
                    if res_norm not in solvent_resnames
                ]
                if polymer_o_coords:
                    fallback_methods["polymer_oxygen"] = "fallback_non_tfsi_excluding_solvent_resnames"
                else:
                    fallback_methods["polymer_oxygen"] = "unresolved"
                    degraded_reasons.append(
                        "polymer oxygen fallback unresolved: solvent/non-solvent separation unavailable."
                    )
            else:
                fallback_methods["polymer_oxygen"] = "unresolved"
                degraded_reasons.append(
                    "polymer oxygen fallback unresolved: missing solvent resname metadata."
                )
        if not solvent_o_coords and preassembly_settings.get("solvent_names"):
            if polymer_resnames:
                solvent_o_coords = [
                    (x, y, z)
                    for x, y, z, res_norm in oxygen_pool_non_tfsi
                    if res_norm not in polymer_resnames
                ]
                if solvent_o_coords:
                    fallback_methods["solvent_oxygen"] = "fallback_non_tfsi_excluding_polymer_resnames"
                else:
                    fallback_methods["solvent_oxygen"] = "unresolved"
                    degraded_reasons.append(
                        "solvent oxygen fallback unresolved: polymer/non-polymer separation unavailable."
                    )
            else:
                fallback_methods["solvent_oxygen"] = "unresolved"
                degraded_reasons.append(
                    "solvent oxygen fallback unresolved: missing polymer resname metadata."
                )

        if polymer_o_coords and solvent_o_coords:
            if set(polymer_o_coords) == set(solvent_o_coords):
                solvent_o_coords = []
                fallback_methods["solvent_oxygen"] = "disabled_identical_pool"
                degraded_reasons.append(
                    "polymer and solvent oxygen fallback pools became identical; solvent metric disabled."
                )

        li_polymer = self._fraction_within_cutoff(
            li_coords,
            polymer_o_coords,
            cutoff_nm=float(preassembly_settings.get("li_polymer_cutoff_nm", 0.35)),
            sample_cap=sample_cap,
        )
        li_solvent = self._fraction_within_cutoff(
            li_coords,
            solvent_o_coords,
            cutoff_nm=float(preassembly_settings.get("li_solvent_cutoff_nm", 0.35)),
            sample_cap=sample_cap,
        )
        li_tfsi_close = self._fraction_within_cutoff(
            li_coords,
            tfsi_coords,
            cutoff_nm=float(preassembly_settings.get("li_tfsi_cutoff_nm", 0.40)),
            sample_cap=sample_cap,
        )

        metrics = {
            "li_count": len(li_coords),
            "polymer_o_count": len(polymer_o_coords),
            "solvent_o_count": len(solvent_o_coords),
            "tfsi_atom_count": len(tfsi_coords),
            "li_within_polymer_o": li_polymer,
            "li_within_solvent_o": li_solvent,
            "li_within_tfsi_close": li_tfsi_close,
            "cutoffs_nm": {
                "li_polymer": float(preassembly_settings.get("li_polymer_cutoff_nm", 0.35)),
                "li_solvent": float(preassembly_settings.get("li_solvent_cutoff_nm", 0.35)),
                "li_tfsi_close": float(preassembly_settings.get("li_tfsi_cutoff_nm", 0.40)),
            },
            "classification": {
                "polymer_resnames": sorted(polymer_resnames),
                "solvent_resnames": sorted(solvent_resnames),
                "tfsi_resnames": sorted(tfsi_resnames),
                "method": "template_resname_plus_atom_element",
                "fallback_methods": fallback_methods,
                "degraded": bool(degraded_reasons),
                "degraded_reasons": degraded_reasons,
            },
        }

        if len(li_coords) == 0:
            return {
                "mode": mode,
                "status": "no_li_detected",
                "passes_thresholds": True,
                "metrics": metrics,
                "thresholds": {},
                "violations": [],
            }

        thresholds = {
            "target_li_polymer_fraction": float(preassembly_settings.get("target_li_polymer_fraction", 0.40)),
            "target_li_solvent_fraction": float(preassembly_settings.get("target_li_solvent_fraction", 0.50)),
            "max_li_tfsi_close_fraction": float(preassembly_settings.get("max_li_tfsi_close_fraction", 0.30)),
        }
        violations: List[str] = []
        li_polymer_fraction = li_polymer.get("fraction")
        li_solvent_fraction = li_solvent.get("fraction")
        li_tfsi_fraction = li_tfsi_close.get("fraction")

        if mode == "li_near_polymer":
            if li_polymer_fraction is not None and li_polymer_fraction < thresholds["target_li_polymer_fraction"]:
                violations.append(
                    f"li_polymer_fraction={li_polymer_fraction:.3f} < {thresholds['target_li_polymer_fraction']:.3f}"
                )
            if li_tfsi_fraction is not None and li_tfsi_fraction > thresholds["max_li_tfsi_close_fraction"]:
                violations.append(
                    f"li_tfsi_close_fraction={li_tfsi_fraction:.3f} > {thresholds['max_li_tfsi_close_fraction']:.3f}"
                )
        elif mode == "li_solvent_shell":
            if li_solvent_fraction is not None and li_solvent_fraction < thresholds["target_li_solvent_fraction"]:
                violations.append(
                    f"li_solvent_fraction={li_solvent_fraction:.3f} < {thresholds['target_li_solvent_fraction']:.3f}"
                )
            if li_tfsi_fraction is not None and li_tfsi_fraction > thresholds["max_li_tfsi_close_fraction"]:
                violations.append(
                    f"li_tfsi_close_fraction={li_tfsi_fraction:.3f} > {thresholds['max_li_tfsi_close_fraction']:.3f}"
                )

        return {
            "mode": mode,
            "status": "ok" if not violations else "threshold_violated",
            "passes_thresholds": not violations,
            "thresholds": thresholds,
            "violations": violations,
            "metrics": metrics,
        }

    def _evaluate_li_contact_safety(
        self,
        pdb_path: Path,
        unit_diagnostics: Dict[str, Any],
        min_contact_nm: float,
        sample_cap: int,
    ) -> Dict[str, Any]:
        """Compute Li-to-nearest-heavy-atom distances and check min-contact threshold."""
        atoms = parse_pdb_atoms(pdb_path)
        scale_factor = unit_diagnostics.get("scale_factor_used")
        if scale_factor is None:
            template_format = str(unit_diagnostics.get("template_format", "pdb"))
            scale_factor = 1.0 if template_format == "gro" else 0.1

        li_coords: List[Tuple[float, float, float]] = []
        heavy_coords: List[Tuple[float, float, float]] = []
        for atom in atoms:
            x_nm = atom["x"] * scale_factor
            y_nm = atom["y"] * scale_factor
            z_nm = atom["z"] * scale_factor
            coord = (x_nm, y_nm, z_nm)
            atom_norm = _normalize_name(atom.get("atom_name", ""))
            element = str(atom.get("element", "")).upper()
            is_li = element == "LI" or atom_norm.startswith("LI")
            if is_li:
                li_coords.append(coord)
                continue
            if element == "H":
                continue
            heavy_coords.append(coord)

        if not li_coords:
            return {
                "li_count": 0,
                "heavy_count": len(heavy_coords),
                "sampled_li_count": 0,
                "sampled_heavy_count": min(len(heavy_coords), sample_cap),
                "min_contact_threshold_nm": min_contact_nm,
                "min_distance_nm": None,
                "median_distance_nm": None,
                "max_distance_nm": None,
                "passes_min_contact": True,
                "status": "no_li_detected",
            }
        if not heavy_coords:
            return {
                "li_count": len(li_coords),
                "heavy_count": 0,
                "sampled_li_count": min(len(li_coords), sample_cap),
                "sampled_heavy_count": 0,
                "min_contact_threshold_nm": min_contact_nm,
                "min_distance_nm": None,
                "median_distance_nm": None,
                "max_distance_nm": None,
                "passes_min_contact": True,
                "status": "no_heavy_atoms_detected",
            }

        li_sample = self._downsample_coords(li_coords, cap=max(1, sample_cap))
        heavy_sample = self._downsample_coords(heavy_coords, cap=max(1, sample_cap))
        min_distances: List[float] = []
        for qx, qy, qz in li_sample:
            best2 = float("inf")
            for tx, ty, tz in heavy_sample:
                dx = qx - tx
                dy = qy - ty
                dz = qz - tz
                d2 = dx * dx + dy * dy + dz * dz
                if d2 < best2:
                    best2 = d2
            if best2 < float("inf"):
                min_distances.append(math.sqrt(best2))

        if not min_distances:
            return {
                "li_count": len(li_coords),
                "heavy_count": len(heavy_coords),
                "sampled_li_count": len(li_sample),
                "sampled_heavy_count": len(heavy_sample),
                "min_contact_threshold_nm": min_contact_nm,
                "min_distance_nm": None,
                "median_distance_nm": None,
                "max_distance_nm": None,
                "passes_min_contact": True,
                "status": "insufficient_samples",
            }

        sorted_d = sorted(min_distances)
        mid = len(sorted_d) // 2
        median = sorted_d[mid] if len(sorted_d) % 2 == 1 else 0.5 * (sorted_d[mid - 1] + sorted_d[mid])
        min_distance = sorted_d[0]
        return {
            "li_count": len(li_coords),
            "heavy_count": len(heavy_coords),
            "sampled_li_count": len(li_sample),
            "sampled_heavy_count": len(heavy_sample),
            "min_contact_threshold_nm": min_contact_nm,
            "min_distance_nm": min_distance,
            "median_distance_nm": median,
            "max_distance_nm": sorted_d[-1],
            "passes_min_contact": min_distance >= min_contact_nm,
            "status": "ok" if min_distance >= min_contact_nm else "dangerous_close_contact",
        }
    
    def _check_polymer_percolation(
        self,
        ctx: "PipelineContext",
        config,
        composition,
        composition_source: str,
    ) -> bool:
        """
        Check polymer percolation guardrails.
        
        For polymer-like systems, enforce minimum molecule counts to ensure
        HTPolyNet can achieve network percolation.
        
        Returns:
            True if checks pass (or warn-only mode), False if stage should fail
        """
        # Identify polymer-like molecules
        polymer_species = []
        polymer_chain_count = 0
        
        for mol in config.molecules:
            if _is_polymer_like(mol.name, mol.mw_g_mol):
                polymer_species.append(mol.name)
                polymer_chain_count += composition.counts.get(mol.name, 0)
        
        is_polymer_system = len(polymer_species) > 0
        
        if not is_polymer_system:
            print(f"  - Polymer check: No polymer-like species detected (non-polymer system)")
            return True
        
        print(f"  - Polymer check: Detected polymer-like species: {', '.join(polymer_species)}")
        print(f"  - Polymer chains: {polymer_chain_count}, Total molecules: {composition.total_molecules}")
        
        # Check if target_total_molecules was explicitly set or defaulted
        if composition_source.startswith("yaml:"):
            # Warn if using default 500
            if config.target_total_molecules == 500:
                print(f"  [WARNING] target_total_molecules=500 (default) may be too low for polymer percolation.")
                print(f"            Consider setting target_total_molecules explicitly in composition.yaml.")
        
        # Check minimum thresholds
        violations = []
        
        if composition.total_molecules < ctx.min_total_molecules_polymer:
            violations.append(
                f"total_molecules={composition.total_molecules} < min={ctx.min_total_molecules_polymer}"
            )
        
        if polymer_chain_count < ctx.min_polymer_chains:
            violations.append(
                f"polymer_chains={polymer_chain_count} < min={ctx.min_polymer_chains}"
            )
        
        if not violations:
            print(f"  - Polymer percolation check: PASSED")
            return True
        
        # Build error/warning message
        msg = (
            f"Polymer percolation risk: {'; '.join(violations)}.\\n"
            f"  For HTPolyNet to form a percolating network, you typically need:\\n"
            f"  - At least {ctx.min_total_molecules_polymer} total molecules\\n"
            f"  - At least {ctx.min_polymer_chains} polymer chains\\n"
            f"  Fix: Increase target_total_molecules in composition.yaml, or adjust thresholds:\\n"
            f"       --min-total-molecules-polymer <N>  --min-polymer-chains <N>\\n"
            f"  To convert to warning-only: --no-strict-polymer-check"
        )
        
        if ctx.strict_polymer_check and not ctx.dry_run:
            print(f"\\n  [ERROR] {msg}")
            self.mark_failed(ctx, msg)
            return False
        else:
            print(f"\\n  [WARNING] {msg}")
            if ctx.dry_run:
                print(f"  [DRY-RUN] Continuing despite polymer threshold violation")
            else:
                print(f"  [INFO] Continuing due to --no-strict-polymer-check")
            return True
    
    def _record_manifest(
        self,
        ctx: "PipelineContext",
        composition,
        box,
        config,
        target_density: float = None,
        achieved_density: float = None,
        retries_used: int = 0,
        published_files: Optional[list] = None,
        publish_status: str = "complete",
        packmol_publish_ok: Optional[bool] = None,
        composition_source: str = "unknown",
        dry_run: bool = False,
        preview_only: bool = False,
        unit_diagnostics: Optional[Dict[str, Any]] = None,
        polymer_info: Optional[Dict[str, Any]] = None,
        retry_attempts: Optional[List["RetryAttemptInfo"]] = None,  # NEW: List of RetryAttemptInfo
        retry_meta: Optional[Dict[str, Any]] = None,
        box_length_target_nm: Optional[float] = None,
        box_length_final_nm: Optional[float] = None,
        density_reduction_used: Optional[bool] = None,
        density_reduction_reason: Optional[str] = None,
        publish_details: Optional[Dict[str, Any]] = None,
        preassembly_info: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Record composition and box data in manifest with canonical keys."""
        if ctx.manifest is None:
            return
        
        target_density_value = target_density if target_density is not None else box.rho0_g_cm3
        achieved_density_value = achieved_density if achieved_density is not None else box.rho0_g_cm3
        requested_rho_value = float(
            getattr(
                box,
                "requested_rho_g_cm3",
                getattr(config, "rho0_g_cm3", target_density_value),
            )
        )
        shrinkage_vol_frac = float(
            getattr(
                box,
                "polymerization_shrinkage_vol_frac",
                getattr(config, "polymerization_shrinkage_vol_frac", 0.0),
            )
        )
        rho_effective_value = float(
            getattr(box, "rho_effective_g_cm3", target_density_value)
        )
        density_candidates = list(getattr(box, "density_candidates_g_cm3", []) or [])
        if not density_candidates:
            try:
                from ..composition import get_density_candidates, compute_effective_density

                density_candidates = get_density_candidates(config)
                rho_effective_value = float(
                    getattr(
                        box,
                        "rho_effective_g_cm3",
                        compute_effective_density(requested_rho_value, shrinkage_vol_frac),
                    )
                )
            except Exception:
                pass
        achieved_over_target_ratio = (
            (achieved_density_value / target_density_value)
            if target_density_value and target_density_value > 0
            else None
        )
        target_box_nm = (
            float(box_length_target_nm)
            if box_length_target_nm is not None
            else float(box.length_nm)
        )
        final_box_nm = (
            float(box_length_final_nm)
            if box_length_final_nm is not None
            else float(box.length_nm)
        )
        box_expand_fraction = (
            (final_box_nm - target_box_nm) / target_box_nm
            if target_box_nm > 0
            else None
        )
        if retry_meta and density_reduction_used is None:
            density_reduction_used = bool(retry_meta.get("density_reduction_used", False))
        if retry_meta and density_reduction_reason is None:
            density_reduction_reason = retry_meta.get("density_reduction_reason")
        if density_reduction_reason is None:
            density_reduction_reason = "disabled_by_default"

        # Extended composition data with CANONICAL keys per README
        comp_data = {
            # Canonical keys (per README)
            "counts": composition.counts,
            "molecule_counts": composition.counts,  # backward compat alias
            "box_size_nm": box.length_nm,
            "box_volume_nm3": box.volume_nm3,
            "density_g_cm3": box.rho0_g_cm3,
            "total_mass_g": composition.total_mass_g,
            "total_molecules": composition.total_molecules,
            "rho0_input_g_cm3": requested_rho_value,
            "rho_effective_g_cm3": rho_effective_value,
            "polymerization_shrinkage_vol_frac": shrinkage_vol_frac,
            
            # Density tracking - NEW: include final_box_nm explicitly
            "target_density_g_cm3": target_density_value,
            "achieved_density_g_cm3": achieved_density_value,
            "achieved_over_target_ratio": achieved_over_target_ratio,
            "box_length_target_nm": target_box_nm,
            "box_length_final_nm": final_box_nm,
            "box_expand_fraction": box_expand_fraction,
            "final_box_nm": final_box_nm,  # backward-compatible alias
            "density_reduction_used": bool(density_reduction_used),
            "density_reduction_reason": density_reduction_reason,
            "density_retries": retries_used,
            
            # Additional metadata
            "mw_used": composition.mw_used,
            "target_wt_pct": composition.target_wt_pct,
            "actual_wt_pct": composition.actual_wt_pct,
            "max_deviation_pct": composition.max_deviation_pct,
            
            # Provenance
            "composition_source": composition_source,
            "dry_run": dry_run,
            "preview_only": preview_only,
            "publish_status": publish_status,
            "packmol_publish_ok": (
                bool(packmol_publish_ok)
                if packmol_publish_ok is not None
                else publish_status == "complete"
            ),
        }
        if density_candidates:
            comp_data["density_candidates_g_cm3"] = density_candidates
        comp_warnings = getattr(composition, "warnings", None)
        if comp_warnings:
            comp_data["warnings"] = comp_warnings
        comp_gel_guardrail = getattr(composition, "gel_guardrail", None)
        if comp_gel_guardrail:
            comp_data["gel_guardrail"] = comp_gel_guardrail
        comp_auto_scaling = getattr(composition, "auto_scaling", None)
        if comp_auto_scaling:
            comp_data["auto_scaling"] = comp_auto_scaling
        comp_molecule_metadata = getattr(composition, "molecule_metadata", None)
        if comp_molecule_metadata:
            comp_data["molecule_metadata"] = comp_molecule_metadata
        
        if published_files:
            comp_data["published_files"] = published_files
        if publish_details:
            comp_data["publish_details"] = publish_details
        
        # Add unit handling diagnostics (Part B)
        if unit_diagnostics:
            comp_data["unit_handling"] = unit_diagnostics
        
        # Add polymer percolation info (Part A)
        if polymer_info:
            comp_data["polymer_percolation"] = polymer_info

        if preassembly_info:
            comp_data["preassembly"] = preassembly_info

        if retry_meta:
            comp_data["retry_policy"] = retry_meta
        
        # NEW: Add retry attempts diagnostics (A1, A2 requirements)
        if retry_attempts:
            # Convert RetryAttemptInfo objects to dicts for JSON serialization
            comp_data["retry_attempts"] = []
            for attempt in retry_attempts:
                attempt_dict = {
                    "attempt_index": attempt.attempt_index,
                    "rho_target": attempt.rho_target,
                    "box_length_nm": attempt.box_length_nm,
                    "tolerance": attempt.tolerance,
                    "seed": attempt.seed,
                    "seed_source": getattr(attempt, "seed_source", "unknown"),
                    "composition_changed": attempt.composition_changed,
                    "success": attempt.success,
                    "inp_file": attempt.inp_file,
                    "pdb_file": attempt.pdb_file,
                    "error_reason": attempt.error_reason,
                }
                comp_data["retry_attempts"].append(attempt_dict)
            
            # Record which attempt succeeded
            successful_attempts = [a for a in retry_attempts if a.success]
            if successful_attempts:
                comp_data["successful_attempt_index"] = successful_attempts[0].attempt_index
                comp_data["successful_attempt_files"] = {
                    "inp": successful_attempts[0].inp_file,
                    "pdb": successful_attempts[0].pdb_file,
                }
            
            # NEW: density_retry_summary for auditing (A2 requirement)
            summary = {
                "total_attempts": len(retry_attempts),
                "successful_attempt": successful_attempts[0].attempt_index if successful_attempts else None,
                "final_density_ratio": None,
                "final_attempt_target_density_ratio": None,
                "floor_enforced": bool(retry_meta and retry_meta.get("effective_density_floor") is not None),
                "all_failed": not successful_attempts,
                "density_reduction_used": bool(density_reduction_used),
                "density_reduction_reason": density_reduction_reason,
            }
            # Correct semantic: final_density_ratio = achieved / target
            target_rho = target_density or box.rho0_g_cm3
            achieved_rho = achieved_density or box.rho0_g_cm3
            if target_rho and target_rho > 0:
                summary["final_density_ratio"] = achieved_rho / target_rho
            # Preserve prior meaning under explicit key for auditing.
            if successful_attempts:
                sa = successful_attempts[0]
                first_attempt = retry_attempts[0]
                if first_attempt.rho_target > 0:
                    summary["final_attempt_target_density_ratio"] = sa.rho_target / first_attempt.rho_target
            comp_data["density_retry_summary"] = summary
        
        ctx.manifest._data["composition"] = comp_data
        ctx.manifest.save()
    
    def _atomic_publish(self, src: Path, dest: Path) -> str:
        """
        Atomically publish a file with content hash.
        
        Uses os.replace for atomic overwrite semantics on all platforms.
        
        Returns:
            SHA256 hash of the published file
        """
        dest.parent.mkdir(parents=True, exist_ok=True)
        
        # Compute hash before publishing
        content_hash = self._compute_file_hash(src)
        
        # Write to temp file in same directory (for atomic replace)
        temp_path = dest.parent / f".{dest.name}.tmp.{os.getpid()}_{uuid.uuid4().hex[:8]}"
        try:
            shutil.copy2(src, temp_path)
            # Use os.replace for atomic overwrite (works on Windows and Linux)
            os.replace(temp_path, dest)
        finally:
            # Clean up temp file if replace failed
            if temp_path.exists():
                try:
                    temp_path.unlink()
                except OSError:
                    pass
        
        return content_hash
    
    def _compute_file_hash(self, file_path: Path) -> str:
        """Compute SHA256 hash of a file."""
        sha256 = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        return sha256.hexdigest()
    
    def _archive_previous_versions(self, system_packmol_dir: Path) -> None:
        """Archive previous published files before overwriting."""
        # Check if there are files to archive
        pdb_file = system_packmol_dir / "pdb" / "initial.pdb"
        gro_file = system_packmol_dir / "gro" / "initial.gro"
        
        files_to_archive = [f for f in [pdb_file, gro_file] if f.exists()]
        
        if not files_to_archive:
            return
        
        # Create archive directory with collision-safe suffixes.
        archive_root = system_packmol_dir / "_archive"
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        archive_dir = archive_root / timestamp
        suffix = 1
        while archive_dir.exists():
            archive_dir = archive_root / f"{timestamp}-{suffix:02d}"
            suffix += 1
        archive_dir.mkdir(parents=True, exist_ok=False)
        
        # Copy files to archive
        for src_file in files_to_archive:
            dest_file = archive_dir / src_file.name
            shutil.copy2(src_file, dest_file)
        
        print(f"  - Archived previous files to: {archive_dir}")
    
    def _write_hash_manifest(self, system_packmol_dir: Path, published_files: list) -> None:
        """Write current.sha256 with content hashes for traceability."""
        hash_file = system_packmol_dir / "current.sha256"
        
        lines = [
            f"# PACKMOL published outputs - generated {datetime.now().isoformat()}",
            "# Format: sha256  filename",
        ]
        
        for pf in published_files:
            path = Path(pf["path"])
            lines.append(f"{pf['sha256']}  {path.name}")
        
        self._atomic_write_text(hash_file, "\n".join(lines) + "\n")

    def _atomic_write_text(self, path: Path, content: str, encoding: str = "utf-8") -> None:
        """Write text atomically via temp-file + replace."""
        path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = path.parent / f".{path.name}.tmp.{os.getpid()}_{uuid.uuid4().hex[:8]}"
        try:
            temp_path.write_text(content, encoding=encoding)
            os.replace(temp_path, path)
        finally:
            if temp_path.exists():
                try:
                    temp_path.unlink()
                except OSError:
                    pass

    def _write_scaled_pdb(self, src_path: Path, dest_path: Path, scale_factor: float) -> None:
        """
        Write a scaled PDB (coordinates only) preserving column widths.
        Only ATOM/HETATM coordinate columns are scaled; other lines are unchanged.
        """
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(src_path, "r") as src, open(dest_path, "w") as dst:
            for line in src:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if len(line) < 54:
                        dst.write(line)
                        continue
                    try:
                        x = float(line[30:38]) * scale_factor
                        y = float(line[38:46]) * scale_factor
                        z = float(line[46:54]) * scale_factor
                        new_line = (
                            f"{line[:30]}"
                            f"{x:8.3f}{y:8.3f}{z:8.3f}"
                            f"{line[54:]}"
                        )
                        dst.write(new_line)
                    except ValueError:
                        dst.write(line)
                else:
                    dst.write(line)
    
    def _convert_pdb_to_gro(
        self,
        ctx: "PipelineContext",
        pdb_path: Path,
        output_dir: Path,
        box_length_nm: float,
        template_format: str,
    ) -> tuple:
        """
        Convert PDB to GRO using gmx editconf with robust unit handling.
        
        Density-preserving policy:
        - Treat the intended PACKMOL box as authoritative for density.
        - Center within the existing box (equivalent to gmx editconf -c).
        - Apply a tiny epsilon translation to avoid boundary artifacts.
        - Only expand the box if the structure cannot fit inside the intended box.
        
        Uses ratio-based unit detection (extent/box) with two-signal detection.
        
        Args:
            ctx: Pipeline context with config options
            pdb_path: Input PDB file
            output_dir: Output directory
            box_length_nm: Cubic box side length in nm
            
        Returns:
            Tuple of (gro_path or None, final_box_nm, unit_diagnostics dict)
        """
        import shlex
        import subprocess
        from ..gromacs_cmd import resolve_gmx_command
        
        margin_nm_requested = ctx.packmol_edge_margin_nm  # Default 0.2 nm
        max_margin_fraction = 0.02  # Adaptive cap (2% of smallest box dimension)
        max_margin_nm = max_margin_fraction * box_length_nm
        margin_nm_effective = min(margin_nm_requested, max_margin_nm)
        max_expand_fraction = max(0.0, float(getattr(ctx, "packmol_max_box_expand_fraction", 0.03)))
        epsilon_nm = 1e-3  # Small translation (<< 0.01 nm) to avoid boundary artifacts
        strict_gro = bool(getattr(ctx, "strict_gro_conversion", False))
        
        unit_diagnostics = {
            "box_length_nm": box_length_nm,
            "scale_factor_used": None,
            "unit_inferred": None,
            "decision_source": None,
            "raw_extent_units": "template_raw",
            "raw_extent_A": None,
            "winsor_extent_A": None,
            "raw_extent_max_A": None,
            "winsor_extent_max_A": None,
            "sampled_nn_distance_stats": None,
            "extent_nm": None,
            "min_coords_raw": None,
            "max_coords_raw": None,
            "required_box_nm": box_length_nm,
            "box_adjusted": False,
            "box_change_reason": "none",  # enum: "none", "extent_exceeds_box"
            "translation_nm": [0.0, 0.0, 0.0],
            "clearance_low_nm": None,
            "clearance_high_nm": None,
            "margin_nm": margin_nm_requested,
            "margin_nm_effective": margin_nm_effective,
            "max_margin_fraction": max_margin_fraction,
            "box_expand_cap_fraction": max_expand_fraction,
            "epsilon_nm": epsilon_nm,
            "intended_box_nm": box_length_nm,
            "box_delta_nm": 0.0,
            "box_delta_fraction": 0.0,
            "expected_density_drop_fraction": 0.0,
            "density_change_expected": False,
            "centering_strategy": "center_in_box_via_translate",
            "input_units_inferred": None,
            "scale_override_source": None,
            "template_format": template_format,
            "gro_conversion_status": "pending",  # enum: "success", "skipped_no_gmx", "skipped_strict_fail", "failed"
        }
        
        gmx_cmd_used = resolve_gmx_command(ctx)
        unit_diagnostics["gmx_cmd_used"] = gmx_cmd_used
        try:
            gmx_tokens = shlex.split(str(gmx_cmd_used))
        except ValueError as e:
            unit_diagnostics["gro_conversion_status"] = "skipped_strict_fail" if strict_gro else "failed"
            raise PackmolConversionError(
                f"Invalid GROMACS command in configuration: {gmx_cmd_used!r} ({e})",
                unit_diagnostics,
            ) from e
        if not gmx_tokens:
            unit_diagnostics["gro_conversion_status"] = "skipped_strict_fail" if strict_gro else "failed"
            raise PackmolConversionError(
                "Resolved GROMACS command is empty; cannot run editconf.",
                unit_diagnostics,
            )
        gmx_frontend = gmx_tokens[0]
        gmx_exe = shutil.which(gmx_frontend)
        unit_diagnostics["gmx_cmd_tokens"] = gmx_tokens
        unit_diagnostics["gmx_frontend_resolved"] = gmx_exe
        if not gmx_exe:
            # Fail-fast in strict mode when GRO conversion is required
            if strict_gro:
                unit_diagnostics["gro_conversion_status"] = "skipped_strict_fail"
                unit_diagnostics["input_units_inferred"] = "unknown"
                unit_diagnostics["scale_override_source"] = "none"
                raise PackmolConversionError(
                    f"PACKMOL requires GRO conversion but '{gmx_frontend}' (from '{gmx_cmd_used}') is not found in PATH.\n"
                    "This is a fail-fast error because --strict-gro-conversion is enabled.\n"
                    "Solutions:\n"
                    f"  1. Install GROMACS and ensure '{gmx_frontend}' is in your PATH\n"
                    "  2. Use --no-strict-gro-conversion to skip GRO conversion (pipeline will use PDB)\n"
                    "Note: Downstream GROMACS stages require initial.gro; skipping may cause failures.",
                    unit_diagnostics,
                )
            # Non-strict: warn loudly and continue
            print(f"  [WARNING] {gmx_frontend} not found, skipping GRO conversion.")
            print(f"  [WARNING] This may cause issues if downstream stages require initial.gro.")
            print(f"  [WARNING] Use --strict-gro-conversion to fail-fast in this situation.")
            unit_diagnostics["gro_conversion_status"] = "skipped_no_gmx"
            unit_diagnostics["input_units_inferred"] = "unknown"
            unit_diagnostics["scale_override_source"] = "none"
            return None, box_length_nm, unit_diagnostics
        
        gro_path = output_dir / "initial.gro"
        
        # === Step 1: Parse PDB extent ===
        try:
            atoms = parse_pdb_atoms(pdb_path)
            coords = [(a["x"], a["y"], a["z"]) for a in atoms]
            min_xyz, max_xyz, extent_xyz = _compute_extent_from_coords(coords)
            winsor_extent_xyz = _compute_winsorized_extent(coords, lower_pct=0.5, upper_pct=99.5)
            nn_stats = _sample_local_nn_distances(atoms, max_pairs=4000)
            extent_max = max(extent_xyz)
            winsor_extent_max = max(winsor_extent_xyz)
            unit_diagnostics["raw_extent_A"] = list(extent_xyz)
            unit_diagnostics["winsor_extent_A"] = list(winsor_extent_xyz)
            unit_diagnostics["raw_extent_max_A"] = extent_max
            unit_diagnostics["winsor_extent_max_A"] = winsor_extent_max
            unit_diagnostics["sampled_nn_distance_stats"] = nn_stats
            unit_diagnostics["min_coords_raw"] = list(min_xyz)
            unit_diagnostics["max_coords_raw"] = list(max_xyz)
            print(f"  - PDB extent (raw): {extent_max:.2f} (from {len(atoms)} atoms)")
            print(f"  - PDB extent (winsorized 0.5-99.5%): {winsor_extent_max:.2f}")
            if nn_stats.get("median") is not None:
                print(
                    "  - Local NN sample (raw units): "
                    f"min={nn_stats['min']:.3f}, median={nn_stats['median']:.3f}, max={nn_stats['max']:.3f}, n={nn_stats['count']}"
                )
            print(f"  - PDB bounds: min={min_xyz}, max={max_xyz}")
        except ValueError as e:
            msg = f"Failed to parse PDB extent: {e}"
            print(f"  [ERROR] {msg}")
            raise PackmolConversionError(msg, unit_diagnostics)
        
        # === Step 2: Determine scale factor ===
        try:
            if template_format == "gro":
                if ctx.packmol_pdb_scale is not None and abs(ctx.packmol_pdb_scale - 1.0) > 1e-6:
                    raise PackmolUnitInferenceError(
                        "GRO template detected. packmol-pdb-scale must be 1.0 for nm inputs. "
                        "Convert templates to PDB in Å or remove the override.",
                        unit_diagnostics,
                    )
                scale_factor = 1.0
                unit_label = "nm"
                unit_diagnostics["decision_source"] = "template_gro"
                unit_diagnostics["unit_inferred"] = "nm"
                unit_diagnostics["scale_override_source"] = "template_gro"
            else:
                scale_factor, unit_label, diag = infer_scale_factor(
                    box_length_nm=box_length_nm,
                    raw_extent_xyz=extent_xyz,
                    winsor_extent_xyz=winsor_extent_xyz,
                    nn_distance_stats=nn_stats,
                    scale_override=ctx.packmol_pdb_scale,
                    strict_mode=getattr(ctx, "strict_packmol_units", True),
                )
                unit_diagnostics.update(diag)
                unit_diagnostics["scale_override_source"] = (
                    "override" if ctx.packmol_pdb_scale is not None else "inferred"
                )
            unit_diagnostics["input_units_inferred"] = unit_label
            unit_diagnostics["scale_factor_used"] = scale_factor
            print(f"  - Unit detection: {unit_diagnostics.get('decision_source', 'unknown')} → {unit_label} (scale={scale_factor})")
        except PackmolUnitInferenceError as e:
            unit_diagnostics.update(getattr(e, "diagnostics", {}) or {})
            unit_diagnostics["gro_conversion_status"] = "skipped_strict_fail"
            raise PackmolUnitInferenceError(str(e), unit_diagnostics)
        except ValueError as e:
            unit_diagnostics["gro_conversion_status"] = "skipped_strict_fail"
            raise PackmolUnitInferenceError(str(e), unit_diagnostics)

        # === Step 3: Convert bounds to nm and compute centering translation ===
        # After scaling, bounds will be:
        min_nm = tuple(c * scale_factor for c in min_xyz)
        max_nm = tuple(c * scale_factor for c in max_xyz)
        extent_nm = tuple(c * scale_factor for c in extent_xyz)
        extent_nm_max = max(extent_nm)
        
        unit_diagnostics["extent_nm"] = extent_nm_max
        
        # Density-preserving policy:
        # 1) Keep the intended PACKMOL box unless the structure cannot fit.
        # 2) Center within the box (equivalent to gmx editconf -c).
        # 3) Apply a tiny epsilon clearance; expand only if needed to fit.
        
        axis_extents = [max_nm[i] - min_nm[i] for i in range(3)]
        max_extent_nm = max(axis_extents)
        # Enforce epsilon clearance from both sides for all axes.
        required_box_nm = max(box_length_nm, max_extent_nm + 2.0 * epsilon_nm)

        if required_box_nm > box_length_nm + 1e-12:
            unit_diagnostics["box_adjusted"] = True
            unit_diagnostics["box_change_reason"] = "extent_or_clearance_exceeds_box"
            unit_diagnostics["box_delta_nm"] = required_box_nm - box_length_nm
            unit_diagnostics["box_delta_fraction"] = (required_box_nm - box_length_nm) / box_length_nm
            unit_diagnostics["density_change_expected"] = True
            density_ratio = (box_length_nm / required_box_nm) ** 3 if required_box_nm > 0 else 1.0
            unit_diagnostics["expected_density_drop_fraction"] = max(0.0, 1.0 - density_ratio)
            if unit_diagnostics["box_delta_fraction"] > max_expand_fraction + 1e-12:
                msg = (
                    f"Required conversion box expansion ({unit_diagnostics['box_delta_fraction']*100:.2f}%) "
                    f"exceeds cap ({max_expand_fraction*100:.2f}%). Template span is too large for the intended box.\n"
                    "Provide a compact/pre-relaxed conformer, or intentionally increase box size via density settings."
                )
                raise PackmolConversionError(msg, unit_diagnostics)

        translation = [0.0, 0.0, 0.0]
        for i in range(3):
            axis_center = 0.5 * (min_nm[i] + max_nm[i])
            translation[i] = 0.5 * required_box_nm - axis_center
        
        unit_diagnostics["translation_nm"] = translation
        unit_diagnostics["required_box_nm"] = required_box_nm
        
        # Compute final clearances after translation
        final_min = [min_nm[i] + translation[i] for i in range(3)]
        final_max = [max_nm[i] + translation[i] for i in range(3)]
        clearance_low = min(final_min)
        clearance_high = min(required_box_nm - fm for fm in final_max)

        # Logical closure: ensure epsilon clearance is truly satisfied.
        if clearance_low < epsilon_nm - 1e-12 or clearance_high < epsilon_nm - 1e-12:
            required_box_nm = max(
                required_box_nm,
                max_extent_nm + 2.0 * epsilon_nm + 1e-6,
            )
            unit_diagnostics["required_box_nm"] = required_box_nm
            unit_diagnostics["box_adjusted"] = True
            unit_diagnostics["box_change_reason"] = "epsilon_clearance_enforced"
            unit_diagnostics["box_delta_nm"] = required_box_nm - box_length_nm
            unit_diagnostics["box_delta_fraction"] = (required_box_nm - box_length_nm) / box_length_nm
            unit_diagnostics["density_change_expected"] = True
            density_ratio = (box_length_nm / required_box_nm) ** 3 if required_box_nm > 0 else 1.0
            unit_diagnostics["expected_density_drop_fraction"] = max(0.0, 1.0 - density_ratio)
            if unit_diagnostics["box_delta_fraction"] > max_expand_fraction + 1e-12:
                msg = (
                    f"Required conversion box expansion ({unit_diagnostics['box_delta_fraction']*100:.2f}%) "
                    f"exceeds cap ({max_expand_fraction*100:.2f}%) after epsilon-clearance enforcement."
                )
                raise PackmolConversionError(msg, unit_diagnostics)
            for i in range(3):
                axis_center = 0.5 * (min_nm[i] + max_nm[i])
                translation[i] = 0.5 * required_box_nm - axis_center
            final_min = [min_nm[i] + translation[i] for i in range(3)]
            final_max = [max_nm[i] + translation[i] for i in range(3)]
            clearance_low = min(final_min)
            clearance_high = min(required_box_nm - fm for fm in final_max)

        unit_diagnostics["clearance_low_nm"] = clearance_low
        unit_diagnostics["clearance_high_nm"] = clearance_high
        unit_diagnostics["margin_satisfied"] = (
            clearance_low >= margin_nm_effective and clearance_high >= margin_nm_effective
        )
        
        if unit_diagnostics["box_adjusted"]:
            print(
                f"  [WARNING] Conversion box adjusted ({unit_diagnostics['box_change_reason']}): "
                f"{box_length_nm:.4f} -> {required_box_nm:.4f} nm "
                f"(+{unit_diagnostics['box_delta_fraction']*100:.2f}%)."
            )
            print(
                "            Expected density drop from conversion box change: "
                f"{unit_diagnostics.get('expected_density_drop_fraction', 0.0)*100:.2f}%"
            )
        elif not unit_diagnostics["margin_satisfied"] and margin_nm_effective > 0.0:
            print(
                f"  [INFO] Requested margin {margin_nm_requested:.4f} nm "
                f"(effective cap {margin_nm_effective:.4f} nm) cannot be fully met in box "
                f"{box_length_nm:.4f} nm; using centering within the intended box"
            )
        
        print(f"  - Translation (centered): ({translation[0]:.4f}, {translation[1]:.4f}, {translation[2]:.4f}) nm")
        print(
            f"  - Clearances: low={clearance_low:.4f} nm, high={clearance_high:.4f} nm "
            f"(requested margin={margin_nm_requested:.4f} nm, effective cap={margin_nm_effective:.4f} nm)"
        )
        
        # === Step 4: Write scaled PDB (nm) if needed ===
        scaled_pdb = pdb_path
        if abs(scale_factor - 1.0) > 1e-9:
            scaled_pdb = output_dir / "initial_scaled.pdb"
            try:
                self._write_scaled_pdb(pdb_path, scaled_pdb, scale_factor)
                unit_diagnostics["scaled_pdb_path"] = str(scaled_pdb)
            except Exception as e:
                unit_diagnostics["gro_conversion_status"] = "failed"
                if strict_gro:
                    raise PackmolConversionError(
                        f"Failed to write scaled PDB for GRO conversion: {e}",
                        unit_diagnostics,
                    ) from e
                print(f"  [ERROR] Failed to write scaled PDB: {e}")
                return None, box_length_nm, unit_diagnostics

        # === Step 5: Run gmx editconf with centering translation (NO -scale flag) ===
        cmd = [
            *gmx_tokens, "editconf",
            "-f", str(scaled_pdb),
            "-o", str(gro_path),
            "-box", str(required_box_nm), str(required_box_nm), str(required_box_nm),
            "-bt", "cubic",
            # Use explicit -translate (centered) to preserve the intended box
            "-translate", str(translation[0]), str(translation[1]), str(translation[2]),
        ]

        print(f"  - Running: {' '.join(str(c) for c in cmd[:8])}...")  # Truncated for readability
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60,
            )
            
            if result.returncode == 0 and gro_path.exists():
                print(f"  - Converted to GRO with box: {required_box_nm:.4f} nm")
                unit_diagnostics["gro_conversion_status"] = "success"
                return gro_path, required_box_nm, unit_diagnostics
            else:
                unit_diagnostics["gro_conversion_status"] = "failed"
                stderr_preview = (result.stderr or result.stdout or "").strip()[:300]
                msg = (
                    "GRO conversion command failed: "
                    f"exit_code={result.returncode}; stderr='{stderr_preview}'"
                )
                if strict_gro:
                    raise PackmolConversionError(msg, unit_diagnostics)
                print(f"  [WARNING] {msg}")
                return None, box_length_nm, unit_diagnostics
                
        except PackmolConversionError:
            raise
        except Exception as e:
            unit_diagnostics["gro_conversion_status"] = "failed"
            if strict_gro:
                raise PackmolConversionError(f"GRO conversion failed: {e}", unit_diagnostics) from e
            print(f"  [WARNING] GRO conversion failed: {e}")
            return None, box_length_nm, unit_diagnostics

