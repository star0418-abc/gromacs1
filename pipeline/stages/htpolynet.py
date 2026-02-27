"""
HTPolyNetStage: Runs HTPolyNet for polymer network generation.

Implements:
- Sandbox workdir enforcement (OUT_GMX/<RUN_ID>/02_htpolynet/workdir/)
- Dynamic system_config.yaml generation from composition
- Forcefield-coupled reaction rules selection
- PACKMOL output integration
- Atomic publish to IN/systems/<SYSTEM_ID>/htpolynet/ with versioning
- Staging to IN/systems/<SYSTEM_ID>/gromacs/ for GROMACS consumption
"""

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Dict, List, Tuple, Any, Union
import shutil
import subprocess
import os
import datetime
import uuid
import re
import yaml

from .base import BaseStage
from ..utils import parse_bool

if TYPE_CHECKING:
    from ..context import PipelineContext

def _normalize_name(value: str) -> str:
    """Normalize molecule names for resilient key matching."""
    return "".join(ch for ch in value.upper() if ch.isalnum())


def _load_yaml_file(path: Path) -> Any:
    """Load a YAML file into Python objects."""
    with open(path, "r") as handle:
        data = yaml.safe_load(handle)
    return data if data is not None else {}


def _compact_dict(d: Dict[str, Any]) -> Dict[str, Any]:
    """
    Return dict with None values removed.
    
    Used to avoid overwriting existing manifest keys with None values
    during dict merge operations.
    """
    return {k: v for k, v in d.items() if v is not None}


def _safe_unlink(path: Path) -> None:
    """Best-effort unlink for temp files/symlinks."""
    try:
        if path.exists() or path.is_symlink():
            path.unlink()
    except OSError:
        pass


class HTPolyNetError(Exception):
    """Error during HTPolyNet execution."""
    pass


class HTPolyNetStage(BaseStage):
    """
    Stage 2: HTPolyNet
    
    Reads:
    - IN/systems/<SYSTEM_ID>/htpolynet/rules/<FF>/ - reaction rules
    - OUT_GMX/<RUN_ID>/01_packmol/initial.pdb - packed system from PACKMOL
    - Manifest composition data (N_i, box_length_nm)
    
    Writes:
    - OUT_GMX/<RUN_ID>/02_htpolynet/system_config.yaml - generated config
    - OUT_GMX/<RUN_ID>/02_htpolynet/workdir/ - sandbox for intermediates
    - OUT_GMX/<RUN_ID>/02_htpolynet/system.{gro,itp,top} - final outputs
    
    Publishes (on success):
    - IN/systems/<SYSTEM_ID>/htpolynet/gro/system_<RUN_ID>.gro + system_current.gro
    - IN/systems/<SYSTEM_ID>/htpolynet/itp/system_<RUN_ID>.itp + system_current.itp
    - IN/systems/<SYSTEM_ID>/htpolynet/top/system_<RUN_ID>.top + system_current.top
    
    Stages (on success):
    - IN/systems/<SYSTEM_ID>/gromacs/gro/system.gro
    - IN/systems/<SYSTEM_ID>/gromacs/itp/system.itp
    - IN/systems/<SYSTEM_ID>/gromacs/top/system.top
    """
    
    @property
    def name(self) -> str:
        return "htpolynet"
    
    @property
    def output_subdir(self) -> str:
        return "02_htpolynet"

    def _coerce_reactive_site_count(self, raw: Any, source_desc: str, molecule_name: str) -> int:
        """Validate and coerce reactive_site_count to a non-negative integer."""
        try:
            value = int(raw)
        except (TypeError, ValueError) as exc:
            raise HTPolyNetError(
                f"Invalid reactive_site_count for molecule '{molecule_name}' from {source_desc}: {raw!r}"
            ) from exc
        if value < 0:
            raise HTPolyNetError(
                f"Invalid reactive_site_count for molecule '{molecule_name}' from {source_desc}: {value} (<0)"
            )
        return value

    def _resolve_reactive_site_counts(
        self,
        ctx: "PipelineContext",
        molecule_counts: Dict[str, int],
    ) -> Tuple[Dict[str, int], Dict[str, str]]:
        """
        Resolve per-species reactive_site_count metadata.

        Source precedence (highest to lowest):
        1) IN/systems/<SYSTEM_ID>/htpolynet/reactive_site_counts.{yaml,yml}
        2) IN/systems/<SYSTEM_ID>/composition.yaml (molecules[].reactive_site_count)
        3) IN/molecules/<NAME>/meta.yaml
        """
        in_dir = ctx.path_manager.in_dir
        values_by_norm: Dict[str, int] = {}
        source_by_norm: Dict[str, str] = {}
        name_by_norm: Dict[str, str] = {}

        def _record(name: str, raw_value: Any, source_desc: str) -> None:
            norm = _normalize_name(name)
            if not norm:
                return
            prior_name = name_by_norm.get(norm)
            if prior_name is not None and prior_name != name:
                prior_source = source_by_norm.get(norm, "<unknown>")
                raise HTPolyNetError(
                    "Reactive-site normalization collision detected:\n"
                    f"  - '{prior_name}' from {prior_source}\n"
                    f"  - '{name}' from {source_desc}\n"
                    f"Both normalize to key '{norm}'. Rename molecules or disambiguate names."
                )
            coerced = self._coerce_reactive_site_count(raw_value, source_desc, name)
            values_by_norm[norm] = coerced
            source_by_norm[norm] = source_desc
            name_by_norm[norm] = name

        # Lowest precedence: per-molecule meta.yaml
        for mol_name in sorted(molecule_counts.keys()):
            meta_path = in_dir / "molecules" / mol_name / "meta.yaml"
            if not meta_path.exists():
                continue
            meta_data = _load_yaml_file(meta_path)
            if isinstance(meta_data, dict) and "reactive_site_count" in meta_data:
                _record(mol_name, meta_data.get("reactive_site_count"), str(meta_path))

        # Middle precedence: composition.yaml molecule entries
        composition_path = in_dir / "systems" / ctx.system_id / "composition.yaml"
        if composition_path.exists():
            composition_data = _load_yaml_file(composition_path)
            if isinstance(composition_data, dict):
                molecules = composition_data.get("molecules", [])
                if isinstance(molecules, list):
                    for item in molecules:
                        if not isinstance(item, dict):
                            continue
                        name = item.get("name")
                        if not name:
                            continue
                        if "reactive_site_count" in item:
                            _record(str(name), item.get("reactive_site_count"), str(composition_path))

        # Highest precedence: system-level reactive_site_counts.yaml
        reactive_count_paths = [
            in_dir / "systems" / ctx.system_id / "htpolynet" / "reactive_site_counts.yaml",
            in_dir / "systems" / ctx.system_id / "htpolynet" / "reactive_site_counts.yml",
        ]
        for reactive_path in reactive_count_paths:
            if not reactive_path.exists():
                continue
            reactive_data = _load_yaml_file(reactive_path)
            if isinstance(reactive_data, dict):
                # Accept either flat mapping {NAME: count} or {molecules: [{name, reactive_site_count}]}
                if "molecules" in reactive_data and isinstance(reactive_data.get("molecules"), list):
                    for item in reactive_data.get("molecules", []):
                        if not isinstance(item, dict):
                            continue
                        name = item.get("name")
                        if not name:
                            continue
                        if "reactive_site_count" not in item:
                            raise HTPolyNetError(
                                f"Missing reactive_site_count entry for molecule '{name}' in {reactive_path}"
                            )
                        _record(str(name), item.get("reactive_site_count"), str(reactive_path))
                else:
                    for name, raw in reactive_data.items():
                        if isinstance(raw, dict) and "reactive_site_count" in raw:
                            _record(str(name), raw.get("reactive_site_count"), str(reactive_path))
                        else:
                            _record(str(name), raw, str(reactive_path))
            else:
                raise HTPolyNetError(
                    f"Invalid reactive-site metadata format in {reactive_path}: expected mapping."
                )

        resolved: Dict[str, int] = {}
        resolved_sources: Dict[str, str] = {}
        missing: List[str] = []
        for mol_name in sorted(molecule_counts.keys()):
            if int(molecule_counts.get(mol_name, 0)) <= 0:
                continue
            norm = _normalize_name(mol_name)
            if norm not in values_by_norm:
                missing.append(mol_name)
                continue
            resolved[mol_name] = values_by_norm[norm]
            resolved_sources[mol_name] = source_by_norm[norm]

        if missing:
            msg = (
                "Missing reactive_site_count metadata for molecules used in HTPolyNet: "
                + ", ".join(missing)
                + ".\nProvide reactive_site_count in one of:\n"
                + f"  - IN/systems/{ctx.system_id}/htpolynet/reactive_site_counts.yaml\n"
                + f"  - IN/systems/{ctx.system_id}/composition.yaml (molecules[].reactive_site_count)\n"
                + "  - IN/molecules/<NAME>/meta.yaml\n"
                + "The pipeline does not infer crosslinker behavior from molecule names."
            )
            raise HTPolyNetError(msg)

        return resolved, resolved_sources
    
    def run(self, ctx: "PipelineContext") -> bool:
        """
        Execute HTPolyNet stage.
        
        1. Load composition from manifest/PACKMOL
        2. Resolve rules file for selected FF
        3. Generate system_config.yaml
        4. Run HTPolyNet in sandbox workdir
        5. Publish outputs atomically with versioning
        6. Stage to GROMACS input pool
        """
        from ..htpolynet_config import (
            HTPolyNetConfig,
            write_system_config,
            resolve_rules_bundle,
        )
        
        print(f"  Running HTPolyNet stage...")
        print(f"  - Forcefield: {ctx.ff}")
        print(f"  - System: {ctx.system_id}")
        
        # Get paths
        output_dir = self.get_output_dir(ctx)
        workdir = output_dir / "workdir"
        workdir.mkdir(parents=True, exist_ok=True)
        
        in_dir = ctx.path_manager.in_dir
        htpolynet_publish_dir = in_dir / "systems" / ctx.system_id / "htpolynet"
        gromacs_stage_dir = in_dir / "systems" / ctx.system_id / "gromacs"
        
        print(f"  - Output dir: {output_dir}")
        print(f"  - Sandbox workdir: {workdir}")
        
        # === Step 1: Load composition from manifest (fail-fast) ===
        composition_data = self._get_composition_from_manifest(ctx)
        # Priority: initial_counts -> counts -> molecule_counts
        molecule_counts = {}
        counts_key_used = None
        if composition_data.get("initial_counts"):
            molecule_counts = composition_data.get("initial_counts", {})
            counts_key_used = "initial_counts"
        elif composition_data.get("counts"):
            molecule_counts = composition_data.get("counts", {})
            counts_key_used = "counts"
        elif composition_data.get("molecule_counts"):
            molecule_counts = composition_data.get("molecule_counts", {})
            counts_key_used = "molecule_counts"

        # Prefer 'box_size_nm', fallback to 'box_length_nm'
        box_length_nm = composition_data.get("box_size_nm") or composition_data.get("box_length_nm")
        
        # Fail-fast: require valid composition data
        if not molecule_counts:
            print(f"  [ERROR] Missing molecule counts in manifest composition")
            print(f"         Run packmol stage first to compute composition.")
            return False
        if counts_key_used:
            print(f"  - Composition source: composition.{counts_key_used}")
        if box_length_nm is None:
            print(f"  [ERROR] Missing box_size_nm in manifest composition")
            return False
        
        print(f"  - Molecules: {sum(molecule_counts.values())} total")
        print(f"  - Box size: {box_length_nm:.3f} nm")
        
        # Get PACKMOL output (initial structure) - prefer IN/ if available
        packmol_output_dir = ctx.get_output_path("01_packmol")
        initial_gro, initial_source = self._find_initial_structure(packmol_output_dir, ctx)
        
        if not initial_gro:
            print(f"  [ERROR] No initial structure found from PACKMOL stage")
            return False
        
        print(f"  - Initial structure: {initial_gro.name} (from {initial_source})")
        
        # Read box from .gro file if available, otherwise use manifest box
        packmol_box_vectors = None
        packmol_box_is_triclinic = False
        packmol_box_warnings: List[str] = []
        box_source = None
        box_scalar_policy = None  # Set when triclinic approximation is used
        original_box_matrix_nm = None  # Fix D: store original 3x3 if triclinic

        allow_triclinic_unsafe_diag = parse_bool(
            getattr(ctx, "allow_triclinic_unsafe_diagonal_approx", None),
            False,
        )
        if initial_gro.suffix.lower() == ".gro":
            try:
                packmol_box_vectors, packmol_box_is_triclinic, packmol_box_warnings = self._read_box_from_gro(
                    initial_gro
                )
            except HTPolyNetError as e:
                # Fix F: Catch box parsing errors and return False
                print(f"  [ERROR] Failed to parse box from GRO file: {e}")
                return False
            
            # Print any box warnings
            for warn in packmol_box_warnings:
                print(f"  [WARN] {warn}")
            
            if packmol_box_vectors:
                # Use box from file (prefer file over manifest for accuracy)
                if packmol_box_is_triclinic:
                    # Fix D: Store original 3x3 matrix for audit
                    original_box_matrix_nm = packmol_box_vectors
                    if not allow_triclinic_unsafe_diag:
                        print(f"  [ERROR] Triclinic box detected in {initial_gro} (9-float GRO box line).")
                        print(f"         HTPolyNet config here assumes orthorhombic/cubic box vectors.")
                        print(f"         Diagonal collapsing is unsafe: it can distort volume/density and PBC geometry.")
                        print(f"         Re-box to orthorhombic/cubic first, then rerun HTPolyNet.")
                        print(f"         Example: gmx editconf -f input.gro -o reboxed_ortho.gro -box Lx Ly Lz -angles 90 90 90")
                        print(f"         Unsafe override (not recommended): --allow-triclinic-unsafe-diagonal-approx")
                        return False
                    # Explicit unsafe escape hatch
                    diag = [packmol_box_vectors[i][i] for i in range(3)]
                    box_length_nm = max(diag)
                    packmol_box_vectors = (diag[0], diag[1], diag[2])
                    box_source = "gro_triclinic_unsafe_diagonal_approx"
                    box_scalar_policy = "UNSAFE_DIAGONAL_APPROX"
                    print(f"  [WARN] ===============================================================")
                    print(f"  [WARN] UNSAFE TRICLINIC DIAGONAL APPROXIMATION ENABLED")
                    print(f"  [WARN] Original 3x3 box matrix will be recorded in manifest for audit.")
                    print(f"  [WARN] Approximated box used for config: {diag[0]:.4f} x {diag[1]:.4f} x {diag[2]:.4f} nm")
                    print(f"  [WARN] This may break periodic geometry and invalidate density-dependent behavior.")
                    print(f"  [WARN] ===============================================================")
                else:
                    # Orthorhombic: keep full vector info, use max for scalar
                    box_length_nm = max(packmol_box_vectors)
                    box_source = "gro"
                    print(f"  - Box from .gro: {packmol_box_vectors[0]:.4f} x {packmol_box_vectors[1]:.4f} x {packmol_box_vectors[2]:.4f} nm")
            elif box_length_nm is None:
                print(f"  [ERROR] Could not read box from .gro file and no box_size_nm in manifest")
                return False
            else:
                box_source = "manifest"
        else:
            # .pdb file - require explicit box from manifest
            if box_length_nm is None:
                ctx_box_length = getattr(ctx, "htpolynet_box_length_nm", None)
                if ctx_box_length is None:
                    ctx_box_length = getattr(ctx, "box_length_nm", None)
                if ctx_box_length is None:
                    ctx_box_length = getattr(ctx, "box_size_nm", None)
                if ctx_box_length is not None:
                    box_length_nm = ctx_box_length
                    box_source = "context"
                    print(f"  - Box size from context: {box_length_nm:.3f} nm")
            if box_length_nm is None:
                print(f"  [ERROR] Using .pdb file but no box_size_nm in manifest")
                print(f"         PDB files do not contain box information.")
                print(f"         Provide box_size_nm in composition or pass a box length via context/CLI.")
                return False
            # Create box_vectors for manifest recording (assume cubic from manifest)
            packmol_box_vectors = (box_length_nm, box_length_nm, box_length_nm)
            if box_source is None:
                box_source = "manifest"
        
        # Record input source and PACKMOL box in manifest
        if ctx.manifest:
            htpolynet_data = ctx.manifest._data.setdefault("htpolynet", {})
            htpolynet_data["initial_structure_source"] = initial_source
            htpolynet_data["initial_structure_file"] = str(initial_gro)
            htpolynet_data["box_source"] = box_source
            htpolynet_data["triclinic_detected"] = packmol_box_is_triclinic
            if box_scalar_policy:
                htpolynet_data["box_scalar_policy"] = box_scalar_policy
            
            # Store PACKMOL box (pre-HTPolyNet) in box section
            box_data = ctx.manifest._data.setdefault("box", {})
            if packmol_box_vectors:
                if packmol_box_is_triclinic:
                    # Fix D: Store BOTH original 3x3 AND approximated diagonal
                    if original_box_matrix_nm:
                        box_data["original_box_matrix_nm"] = [list(row) for row in original_box_matrix_nm]
                    box_data["packmol_box_diag_nm"] = list(packmol_box_vectors)
                    box_data["packmol_box_is_triclinic"] = True
                    if box_scalar_policy:
                        box_data["box_scalar_policy"] = box_scalar_policy
                else:
                    box_data["packmol_box_vectors_nm"] = list(packmol_box_vectors)
                    box_data["packmol_box_is_triclinic"] = False
                box_data["packmol_box_source"] = str(initial_gro)
            
            if packmol_box_warnings:
                htpolynet_data["packmol_box_warnings"] = packmol_box_warnings
        
        # === Step 2: Resolve rules file ===
        try:
            rules_dir, rules_file, rules_source = resolve_rules_bundle(ctx.system_id, ctx.ff, in_dir)
        except FileNotFoundError as e:
            print(f"  [ERROR] HTPolyNet rules resolution failed:\n{e}")
            return False
        
        print(f"  - Rules dir ({rules_source}): {rules_dir}")
        print(f"  - Rules file: {rules_file}")
        if ctx.manifest:
            ctx.manifest.add_input_asset(
                rules_file,
                ctx.path_manager.hash_file(rules_file),
                category="htpolynet_rules",
            )
        
        # === Step 3: Generate system_config.yaml ===
        curing_steps = getattr(ctx, "htpolynet_curing_steps", None)
        if curing_steps is not None and not self._validate_curing_steps(curing_steps):
            print(f"  [ERROR] Invalid htpolynet_curing_steps format.")
            print(f"         Expected list of dicts with temperature_K and duration_ps.")
            return False

        try:
            reactive_site_counts, reactive_site_sources = self._resolve_reactive_site_counts(
                ctx=ctx,
                molecule_counts=molecule_counts,
            )
        except HTPolyNetError as e:
            print(f"  [ERROR] {e}")
            return False

        crosslinkers: Dict[str, Dict[str, Any]] = {}
        chain_stoppers: Dict[str, Dict[str, Any]] = {}
        for mol_name in sorted(molecule_counts.keys()):
            count = int(molecule_counts[mol_name])
            if count <= 0:
                continue
            site_count = int(reactive_site_counts[mol_name])
            if site_count >= 2:
                crosslinkers[mol_name] = {
                    "count": count,
                    "reactive_site_count": site_count,
                    "classification": (
                        "multifunctional_crosslinker" if site_count >= 3 else "crosslinker"
                    ),
                }
            elif site_count == 1:
                chain_stoppers[mol_name] = {
                    "count": count,
                    "reactive_site_count": site_count,
                }

        print(f"  - Reactive-site metadata resolved for {len(reactive_site_counts)} species")
        print(f"  - Crosslinkers (reactive_site_count >=2): {len(crosslinkers)}")
        print(f"  - Chain-stoppers (reactive_site_count ==1): {len(chain_stoppers)}")
        if ctx.manifest:
            ht_data = ctx.manifest._data.setdefault("htpolynet", {})
            ht_data["reactive_site_counts"] = dict(sorted(reactive_site_counts.items()))
            ht_data["reactive_site_count_sources"] = dict(sorted(reactive_site_sources.items()))
            ht_data["crosslinkers"] = dict(sorted(crosslinkers.items()))
            ht_data["chain_stoppers"] = dict(sorted(chain_stoppers.items()))

        gelation_ok, gelation_summary = self._run_gelation_precheck(
            ctx=ctx,
            molecule_counts=molecule_counts,
            reactive_site_counts=reactive_site_counts,
        )
        if ctx.manifest:
            ctx.manifest._data.setdefault("htpolynet", {})["gelation_precheck"] = gelation_summary
            ctx.manifest.save()
        if not gelation_ok:
            for reason in gelation_summary.get("fail_reasons", []):
                print(f"  [ERROR] Gelation precheck: {reason}")
            print(f"         Override options:")
            print(f"           --skip-gelation-precheck")
            print(f"           --allow-gelation-precheck-fail")
            return False
        if gelation_summary.get("status") == "failed_overridden":
            for reason in gelation_summary.get("fail_reasons", []):
                print(f"  [WARN] Gelation precheck override: {reason}")

        config = HTPolyNetConfig(
            project_name=f"{ctx.system_id}_{ctx.run_id}",
            forcefield=ctx.ff,
            initial_gro=initial_gro,
            workdir=workdir,
            output_dir=output_dir,
            rules_file=rules_file,
            molecule_counts=molecule_counts,
            box_length_nm=box_length_nm,
            box_vectors_nm=packmol_box_vectors if packmol_box_vectors else None,
            reactive_site_counts=reactive_site_counts,
            crosslinkers=crosslinkers,
            chain_stoppers=chain_stoppers,
            curing_steps=curing_steps,
        )
        
        config_path = output_dir / "system_config.yaml"
        write_system_config(config, config_path)
        print(f"  - Generated config: {config_path.name}")

        # === Step 4: Run HTPolyNet ===
        placeholder_used = False
        placeholder_propagated = False
        placeholder_poison_pill_enabled = False
        placeholder_poison_pill_include: Optional[str] = None
        allow_placeholder_stage_to_gromacs = parse_bool(
            getattr(ctx, "allow_placeholder_stage_to_gromacs", None),
            False,
        )
        allow_placeholder_gromacs_compile = parse_bool(
            getattr(ctx, "allow_placeholder_gromacs_compile", None),
            False,
        )

        outputs: Dict[str, Path] = {}
        output_dir_selected: Optional[Path] = None
        ndx_path: Optional[Path] = None

        if ctx.dry_run:
            print(f"\n  [DRY-RUN] Would execute htpolynet run -f {config_path}")
            print(f"  [DRY-RUN] Skipping actual HTPolyNet execution")
            self._record_manifest(
                ctx,
                config,
                rules_file,
                rules_dir=rules_dir,
                rules_source=rules_source,
                molecule_counts_pre=molecule_counts,
                dry_run=True,
                composition_key_used=counts_key_used,
                box_source=box_source,
            )
            return True

        try:
            success, outputs, output_dir_selected, ndx_path = self._run_htpolynet(config_path, workdir, ctx)
            if not success:
                print(f"  [ERROR] HTPolyNet execution failed")
                return False
        except HTPolyNetError as e:
            print(f"  [ERROR] HTPolyNet error: {e}")
            ctx.log_command(
                f"htpolynet run -f {config_path}",
                self.name,
                exit_code=1,
                stderr=str(e),
            )
            return False
        except FileNotFoundError as e:
            print(f"  [ERROR] HTPolyNet not found: {e}")
            if not ctx.allow_placeholder:
                print(f"  Pipeline must fail-fast when htpolynet is unavailable.")
                print(f"  Use --allow-placeholder for demo/testing only.")
                return False

            placeholder_used = True
            placeholder_propagated = parse_bool(getattr(ctx, "allow_placeholder_propagate", None), False)
            print(f"  [WARN] --allow-placeholder enabled, creating NON-PHYSICAL outputs")
            print(f"         These outputs are NOT valid for production simulations.")
            outputs = self._create_placeholder_outputs(output_dir, ctx)
            output_dir_selected = output_dir

            if placeholder_propagated:
                print(f"  [WARN] ===============================================================")
                print(f"  [WARN] DEMO MODE: Placeholder outputs will be published to IN/htpolynet")
                print(f"  [WARN] These files are not physically meaningful.")
                print(f"  [WARN] Manifest flag: htpolynet.placeholder_propagated_downstream_unsafe=true")
                print(f"  [WARN] ===============================================================")
                if allow_placeholder_gromacs_compile:
                    print(f"  [WARN] --allow-placeholder-gromacs-compile enabled; poison pill will NOT be injected.")
                else:
                    placeholder_poison_pill_include = self._apply_placeholder_poison_pill(
                        outputs.get("top"),
                        ctx.run_id,
                    )
                    if placeholder_poison_pill_include:
                        placeholder_poison_pill_enabled = True
                        print(f"  [WARN] Poison pill enabled: system.top includes missing '{placeholder_poison_pill_include}'.")
                        print(f"  [WARN] This intentionally prevents accidental successful grompp.")
            else:
                print(f"  [WARN] Placeholders remain in OUT_GMX only (not published/staged).")
                print(f"         Use --allow-placeholder-propagate for explicit demo publish behavior.")
                self._record_manifest(
                    ctx,
                    config,
                    rules_file,
                    rules_dir=rules_dir,
                    rules_source=rules_source,
                    molecule_counts_pre=molecule_counts,
                    placeholder=placeholder_used,
                    placeholder_propagated=placeholder_propagated,
                    placeholder_poison_pill_enabled=placeholder_poison_pill_enabled,
                    placeholder_poison_pill_include=placeholder_poison_pill_include,
                    published_files=[],
                    staged_files=[],
                    composition_key_used=counts_key_used,
                    selected_output_dir=output_dir_selected,
                    ndx_path=ndx_path,
                    box_source=box_source,
                )
                return False

        # === Step 5: Publish outputs atomically (real or explicit demo placeholders) ===
        print(f"\n  Publishing to IN/systems/{ctx.system_id}/htpolynet/...")
        published_files: List[str] = []
        try:
            published_files = self._publish_outputs(outputs, htpolynet_publish_dir, ctx.run_id)
            for pf in published_files:
                print(f"  - Published: {Path(pf).name}")
        except Exception as e:
            print(f"  [ERROR] Publish failed: {e}")
            return False

        # === Step 6: Stage to GROMACS input pool (policy-aware for placeholders) ===
        staged_files: List[str] = []
        should_stage_to_gromacs = True
        if placeholder_used and placeholder_propagated and not allow_placeholder_stage_to_gromacs:
            should_stage_to_gromacs = False
            print(f"\n  [WARN] Placeholder propagation active: staging to IN/.../gromacs is disabled by default.")
            print(f"         To stage demo placeholders downstream, rerun with:")
            print(f"           --allow-placeholder-stage-to-gromacs")
            print(f"         Downstream grompp remains blocked by poison pill unless:")
            print(f"           --allow-placeholder-gromacs-compile")

        if should_stage_to_gromacs:
            print(f"\n  Staging to IN/systems/{ctx.system_id}/gromacs/...")
            try:
                staged_files = self._stage_to_gromacs(
                    outputs,
                    gromacs_stage_dir,
                    run_id=ctx.run_id,
                    ndx_path=ndx_path,
                )
                for sf in staged_files:
                    print(f"  - Staged: {Path(sf).name}")
            except Exception as e:
                print(f"  [ERROR] Staging failed: {e}")
                print(f"  Cannot proceed without staged files for GROMACS.")
                return False
        else:
            print(f"  - Staging skipped by placeholder safety policy.")

        # === Step 7: Post-publish completeness check ===
        staged_gro = gromacs_stage_dir / "gro" / "system.gro"
        staged_itp = gromacs_stage_dir / "itp" / "system.itp"
        staged_top = gromacs_stage_dir / "top" / "system.top"
        post_htpolynet_gro: Optional[Path] = None

        if should_stage_to_gromacs:
            expected_staged = [staged_gro, staged_itp, staged_top]
            missing_staged = [str(p) for p in expected_staged if not p.exists()]
            if missing_staged:
                print(f"  [ERROR] Post-publish completeness check FAILED!")
                print(f"  Missing expected staged files:")
                for m in missing_staged:
                    print(f"    - {m}")
                print(f"  HTPolyNet may have produced partial outputs or staging failed silently.")
                return False
            print(f"  - Completeness check: all expected staged files present")
            post_htpolynet_gro = staged_gro
        else:
            post_htpolynet_gro = outputs.get("gro")

        # === Step 8: Parse composition_after from staged/published system.top ===
        composition_after = None
        top_for_composition = staged_top if staged_top.exists() else outputs.get("top")
        if top_for_composition and top_for_composition.exists():
            composition_after = self._parse_molecules_from_top(top_for_composition)
            if composition_after:
                print(
                    f"  - Parsed [molecules] from {top_for_composition.name}: "
                    f"{sum(composition_after.values())} molecules"
                )

        # === Step 9: Record in manifest (single source of truth for placeholder state) ===
        self._record_manifest(
            ctx,
            config,
            rules_file,
            rules_dir=rules_dir,
            rules_source=rules_source,
            molecule_counts_pre=molecule_counts,
            molecule_counts_post=composition_after,
            published_files=published_files,
            staged_files=staged_files,
            dry_run=False,
            placeholder=placeholder_used,
            placeholder_propagated=placeholder_propagated,
            placeholder_poison_pill_enabled=placeholder_poison_pill_enabled,
            placeholder_poison_pill_include=placeholder_poison_pill_include,
            post_htpolynet_gro=post_htpolynet_gro,
            composition_key_used=counts_key_used,
            selected_output_dir=output_dir_selected,
            ndx_path=ndx_path,
            box_source=box_source,
        )

        print(f"\n  HTPolyNet stage completed successfully!")
        return True
    
    def _get_composition_from_manifest(self, ctx: "PipelineContext") -> Dict:
        """Load composition data from manifest."""
        if ctx.manifest is None:
            return {}
        
        return ctx.manifest._data.get("composition", {})
    
    def _read_box_from_gro(
        self,
        gro_path: Path,
    ) -> Tuple[Optional[Union[Tuple[float, float, float], Tuple[Tuple[float, ...], ...]]], bool, List[str]]:
        """
        Read box vectors from the last NON-EMPTY line of a GRO file.
        
        GRO format: The last line contains box vectors in nm:
        - Cubic/orthorhombic: v1x v2y v3z (3 values)
        - Triclinic: v1x v2y v3z v1y v1z v2x v2z v3x v3y (9 values)
        
        Fix B: Robust parsing that handles:
        - Trailing blank lines
        - Validates ALL tokens parse as floats
        - Accepts ONLY 3 or 9 float tokens (otherwise treats as corrupted)
        
        Returns:
            Tuple of:
            - Box vectors (3-tuple for ortho, 3x3 for triclinic), or None on failure
            - is_triclinic flag
            - List of warning messages
        
        Raises:
            HTPolyNetError: If file is corrupted (invalid box line format)
        """
        warnings: List[str] = []
        
        if not gro_path.exists():
            return None, False, warnings
        
        try:
            with open(gro_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 3:
                raise HTPolyNetError(
                    f"GRO file {gro_path} has fewer than 3 lines. "
                    f"Expected: title, atom count, coordinates, box. "
                    f"Possible causes: truncated file, partial write, disk full, interrupted transfer."
                )
            
            # Fix B: Find last NON-EMPTY line as box candidate
            box_line = None
            for i in range(len(lines) - 1, -1, -1):
                stripped = lines[i].strip()
                if stripped:
                    box_line = stripped
                    break
            
            if box_line is None:
                raise HTPolyNetError(
                    f"GRO file {gro_path} contains no non-empty lines. "
                    f"Possible causes: file is empty, corrupted, or only whitespace."
                )
            
            parts = box_line.split()
            
            # Fix B: Validate ALL tokens are floats before accepting
            parsed_floats: List[float] = []
            for token in parts:
                try:
                    parsed_floats.append(float(token))
                except ValueError:
                    # Not a pure float box line - likely a corrupted atom line
                    raise HTPolyNetError(
                        f"GRO box line in {gro_path} contains non-numeric token: '{token}'\n"
                        f"Full line: '{box_line}'\n"
                        f"Possible causes: truncated file, partial write, disk full, interrupted transfer.\n"
                        f"Expected format: 3 floats (orthorhombic) or 9 floats (triclinic)."
                    )
            
            # Fix B: Accept ONLY exactly 3 or 9 float tokens
            float_count = len(parsed_floats)
            if float_count not in (3, 9):
                raise HTPolyNetError(
                    f"GRO box line in {gro_path} has {float_count} numeric values, expected 3 or 9.\n"
                    f"Full line: '{box_line}'\n"
                    f"Possible causes: truncated file (missing box line), corrupted coordinates, "
                    f"atom line with velocities being misread as box.\n"
                    f"Verify file integrity and ensure it ends with a valid box vector line."
                )
            
            if float_count == 9:
                # Triclinic box: v1x v2y v3z v1y v1z v2x v2z v3x v3y
                # Store as 3x3 matrix: [[v1x, v1y, v1z], [v2x, v2y, v2z], [v3x, v3y, v3z]]
                v1x, v2y, v3z = parsed_floats[0], parsed_floats[1], parsed_floats[2]
                v1y, v1z = parsed_floats[3], parsed_floats[4]
                v2x, v2z = parsed_floats[5], parsed_floats[6]
                v3x, v3y = parsed_floats[7], parsed_floats[8]
                triclinic_matrix = (
                    (v1x, v1y, v1z),
                    (v2x, v2y, v2z),
                    (v3x, v3y, v3z),
                )
                warnings.append(f"Triclinic box detected: storing full 3x3 matrix")
                return triclinic_matrix, True, warnings
            else:  # float_count == 3
                box_x, box_y, box_z = parsed_floats[0], parsed_floats[1], parsed_floats[2]
                
                # Anisotropy check
                max_dim = max(box_x, box_y, box_z)
                if max_dim > 0:
                    anisotropy_xy = abs(box_x - box_y) / max_dim
                    anisotropy_yz = abs(box_y - box_z) / max_dim
                    anisotropy_xz = abs(box_x - box_z) / max_dim
                    max_anisotropy = max(anisotropy_xy, anisotropy_yz, anisotropy_xz)
                    if max_anisotropy > 0.02:  # 2% threshold
                        warnings.append(
                            f"Box anisotropy detected: Lx={box_x:.4f}, Ly={box_y:.4f}, Lz={box_z:.4f} nm. "
                            f"Max anisotropy: {max_anisotropy*100:.1f}%. "
                            f"This may affect periodic boundary conditions and analysis."
                        )
                return (box_x, box_y, box_z), False, warnings
            
        except HTPolyNetError:
            # Re-raise our own errors
            raise
        except (IOError, OSError) as e:
            raise HTPolyNetError(
                f"Failed to read GRO file {gro_path}: {e}\n"
                f"Check file permissions and disk health."
            )
    
    def _count_atoms_in_gro(self, gro_path: Path) -> Optional[int]:
        """
        Count atoms in a GRO file from the second line.
        
        Returns:
            Number of atoms, or None if parsing fails
        """
        if not gro_path.exists():
            return None
        
        try:
            with open(gro_path, 'r') as f:
                _title = f.readline()  # First line: title
                count_line = f.readline().strip()  # Second line: atom count
                return int(count_line)
        except (ValueError, IOError):
            return None
    
    def _find_initial_structure(
        self,
        packmol_dir: Path,
        ctx: "PipelineContext",
    ) -> Tuple[Optional[Path], str]:
        """
        Find initial structure with strict provenance enforcement.
        
        Priority order:
        1. IN/systems/<SYSTEM_ID>/packmol/gro/ (.gro preferred over .pdb)
        2. IN/systems/<SYSTEM_ID>/packmol/pdb/
        3. OUT fallback ONLY if ctx.unsafe_allow_out_fallback is True
        
        Strict mode: Inputs must come from IN/published. OUT fallback is
        disabled by default to prevent "stale artifact" contamination.
        
        Returns:
            Tuple of (path, source_description)
        """
        # Try IN/systems/<SYSTEM_ID>/packmol/ first (published outputs)
        # Prefer .gro over .pdb for better coordinate precision and box info
        published_dir = ctx.get_input_path("systems", ctx.system_id, "packmol")
        
        # Priority: gro subdir first, then pdb subdir
        for subdir in ["gro", "pdb"]:
            subdir_path = published_dir / subdir
            if not subdir_path.exists():
                continue
            # Prefer .gro files, then .pdb
            for ext in [".gro", ".pdb"]:
                for prefix in ["initial", "system_current", "system"]:
                    candidate = subdir_path / f"{prefix}{ext}"
                    if candidate.exists():
                        return candidate, f"IN/published/{subdir}"
        
        # OUT fallback: only if explicitly allowed (unsafe mode)
        unsafe_fallback = getattr(ctx, 'unsafe_allow_out_fallback', False)
        if unsafe_fallback:
            for prefix in ["initial", "packed_system", "system"]:
                for ext in [".gro", ".pdb"]:
                    candidate = packmol_dir / f"{prefix}{ext}"
                    if candidate.exists():
                        print(f"  [WARN] ============================================")
                        print(f"  [WARN] Using OUT fallback for initial structure!")
                        print(f"  [WARN] Path: {candidate}")
                        print(f"  [WARN] This violates strict input provenance.")
                        print(f"  [WARN] Consider publishing PACKMOL outputs to IN/ first.")
                        print(f"  [WARN] ============================================")
                        return candidate, "OUT/packmol (UNSAFE FALLBACK)"
        else:
            # Check if OUT files exist but are inaccessible due to strict mode
            out_candidates = []
            for prefix in ["initial", "packed_system", "system"]:
                for ext in [".gro", ".pdb"]:
                    candidate = packmol_dir / f"{prefix}{ext}"
                    if candidate.exists():
                        out_candidates.append(str(candidate))
            if out_candidates:
                print(f"  [INFO] Found {len(out_candidates)} candidate(s) in OUT/, but OUT fallback is disabled.")
                print(f"         Candidates: {out_candidates[:3]}..." if len(out_candidates) > 3 else f"         Candidates: {out_candidates}")
                print(f"         To use OUT fallback, set --unsafe-allow-out-fallback (not recommended).")
                print(f"         Preferred: Run PACKMOL stage to publish outputs to IN/.")
        
        return None, "not_found"
    
    def _run_htpolynet(
        self,
        config_path: Path,
        workdir: Path,
        ctx: "PipelineContext",
    ) -> Tuple[bool, Dict[str, Path], Optional[Path], Optional[Path]]:
        """
        Execute HTPolyNet.
        
        Returns:
            Tuple of (success, output_files_dict, selected_output_dir, optional_ndx_path)
        """
        htpolynet_exe = shutil.which("htpolynet")
        if not htpolynet_exe:
            raise FileNotFoundError("htpolynet executable not found in PATH")
        
        # Build command - run in workdir to sandbox all intermediates
        cmd = [
            htpolynet_exe, "run",
            "-f", str(config_path),
            "-d", str(workdir),
        ]
        
        print(f"  - Executing: {' '.join(cmd)}")
        
        # Bug #5: Respect ctx.omp_num_threads, do NOT hard-code
        env = self._build_env(ctx)
        if ctx.omp_num_threads is not None:
            print(f"  - OMP_NUM_THREADS: {ctx.omp_num_threads}")
        
        # Configurable timeout (default: 43200s = 12 hours)
        # None means no timeout (explicitly allowed but warned)
        timeout_sec = getattr(ctx, 'htpolynet_timeout_sec', 43200)
        if timeout_sec is None:
            print(f"  [WARN] HTPolyNet timeout is disabled (timeout=None). Job may hang forever.")
            print(f"         Consider setting a finite timeout via --htpolynet-timeout.")
        else:
            print(f"  - Timeout: {timeout_sec}s ({timeout_sec/3600:.1f}h)")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=workdir,
                timeout=timeout_sec,
                env=env,
            )
        except subprocess.TimeoutExpired:
            raise HTPolyNetError(
                f"HTPolyNet timed out after {timeout_sec}s ({timeout_sec/3600:.1f}h). "
                f"For large/high-crosslink systems, increase timeout via:\n"
                f"  --htpolynet-timeout <seconds>  (e.g., 86400 for 24h)\n"
                f"  or set ctx.htpolynet_timeout_sec in code/config."
            )
        
        ctx.log_command(
            " ".join(cmd),
            self.name,
            exit_code=result.returncode,
            stdout=result.stdout[-500:] if result.stdout else None,
            stderr=result.stderr[-500:] if result.stderr else None,
        )

        def _record_qc_metrics(metrics: Dict[str, Any]) -> None:
            if ctx.manifest:
                ctx.manifest.set_htpolynet_metrics(
                    conversion=metrics.get("conversion"),
                    gel_fraction=metrics.get("gel_fraction"),
                    qc_status=metrics.get("status", "unknown"),
                    warnings=metrics.get("warnings", []),
                )

        qc_metrics = self._parse_htpolynet_qc(result.stdout or "", result.stderr or "")

        if result.returncode != 0:
            _record_qc_metrics(qc_metrics)
            # Print last 50 lines of output for debugging
            if result.stdout:
                lines = result.stdout.strip().split("\n")
                print(f"\n  Last 50 lines of stdout:")
                for line in lines[-50:]:
                    print(f"    {line}")
            if result.stderr:
                print(f"\n  Stderr: {result.stderr[-500:]}")
            return False, {}, None, None

        # Strict QC parse contract for successful real HTPolyNet runs.
        allow_missing_qc = parse_bool(getattr(ctx, "allow_missing_htpolynet_qc", None), False)
        if qc_metrics.get("conversion") is None and qc_metrics.get("gel_fraction") is None:
            qc_msg = (
                "Could not parse HTPolyNet QC labels ('conversion' or 'gel fraction') from stdout/stderr."
            )
            qc_metrics.setdefault("warnings", []).append(qc_msg)
            if allow_missing_qc:
                qc_metrics["status"] = "missing_allowed"
                print(f"  [WARN] {qc_msg}")
                print(f"         Proceeding only because --allow-missing-htpolynet-qc is set.")
            else:
                qc_metrics["status"] = "failed_missing_qc"
                print(f"  [ERROR] {qc_msg}")
                print(f"         Refusing to continue to avoid silent non-crosslinked outputs downstream.")
                print(f"         Override (unsafe): --allow-missing-htpolynet-qc")
                _record_qc_metrics(qc_metrics)
                return False, {}, None, None

        # Keep minimum conversion behavior.
        min_conversion = getattr(ctx, 'min_conversion', None)
        conversion_parsed = qc_metrics.get("conversion")
        conversion_status = qc_metrics.get("status", "unknown")

        if min_conversion is not None:
            if conversion_parsed is not None:
                # Check against threshold
                if conversion_parsed < min_conversion:
                    print(f"  [ERROR] Conversion {conversion_parsed:.1%} below threshold {min_conversion:.1%}")
                    qc_metrics["status"] = "failed_min_conversion"
                    _record_qc_metrics(qc_metrics)
                    return False, {}, None, None
            else:
                # Conversion not parsed but threshold was set - fail to prevent silent bypass
                print(f"  [ERROR] min_conversion={min_conversion:.1%} set but conversion could not be parsed from HTPolyNet output")
                print(f"         qc_status={conversion_status}")
                print(f"         This prevents silent bypass of QC threshold.")
                print(f"         To proceed without conversion checking, remove --min-conversion flag.")
                qc_metrics["status"] = "failed_missing_conversion_for_threshold"
                _record_qc_metrics(qc_metrics)
                return False, {}, None, None

        # Find output files (fail-fast on partial outputs), then mark completion.
        _record_qc_metrics(qc_metrics)
        outputs, selected_dir, ndx_path = self._find_htpolynet_outputs(workdir, config_path.parent)
        self._write_htpolynet_completion_marker(
            workdir=workdir,
            config_path=config_path,
            selected_output_dir=selected_dir,
            run_id=ctx.run_id,
        )
        return True, outputs, selected_dir, ndx_path
    
    def _build_env(self, ctx: "PipelineContext") -> Dict[str, str]:
        """
        Build environment for HTPolyNet subprocess.
        
        Bug #5: Respect ctx.omp_num_threads instead of hard-coding.
        """
        env = os.environ.copy()
        if ctx.omp_num_threads is not None:
            env["OMP_NUM_THREADS"] = str(ctx.omp_num_threads)
        return env

    def _make_unique_temp_path(self, directory: Path, label: str, run_id: str) -> Path:
        """Build a unique temp path safe for concurrent workers."""
        safe_label = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)
        token = uuid.uuid4().hex
        return directory / f".{safe_label}.{run_id}.pid{os.getpid()}.{token}.tmp"

    def _write_htpolynet_completion_marker(
        self,
        workdir: Path,
        config_path: Path,
        selected_output_dir: Path,
        run_id: str,
    ) -> None:
        """
        Write a completion marker only after validated final outputs are located.
        """
        marker_path = workdir / "htpolynet.completed"
        tmp_path = self._make_unique_temp_path(workdir, "htpolynet.completed", run_id)
        content = (
            f"completed_at={datetime.datetime.now().isoformat()}\n"
            f"config_path={config_path}\n"
            f"selected_output_dir={selected_output_dir}\n"
        )
        try:
            tmp_path.write_text(content)
            os.replace(tmp_path, marker_path)
        finally:
            _safe_unlink(tmp_path)
    
    def _parse_htpolynet_qc(self, stdout: str, stderr: str) -> Dict[str, Any]:
        """
        Parse HTPolyNet output for crosslinking QC metrics using anchored labels only.

        Accepted labels (case-insensitive):
          - ^(total )?conversion[:=] <value> (%|fraction)
          - ^gel fraction[:=] <value> (%|fraction)
        """
        qc: Dict[str, Any] = {
            "status": "missing",
            "conversion": None,
            "gel_fraction": None,
            "warnings": [],
        }

        number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
        conversion_pat = re.compile(
            rf"^\s*(?:total\s+)?conversion\s*[:=]\s*(?P<value>{number})\s*(?P<unit>%|fraction)\b.*$",
            flags=re.IGNORECASE,
        )
        gel_fraction_pat = re.compile(
            rf"^\s*gel\s*fraction\s*[:=]\s*(?P<value>{number})\s*(?P<unit>%|fraction)\b.*$",
            flags=re.IGNORECASE,
        )

        def _convert_metric(value_str: str, unit: str, raw_line: str, label: str) -> Optional[float]:
            try:
                value = float(value_str)
            except ValueError:
                qc["warnings"].append(f"Failed to parse numeric {label} from line: {raw_line.strip()}")
                return None
            unit_norm = unit.lower()
            if unit_norm == "%":
                converted = value / 100.0
            else:
                converted = value
            if converted < 0.0 or converted > 1.0:
                qc["warnings"].append(
                    f"Ignoring out-of-range {label}={converted:.6g} parsed from line: {raw_line.strip()}"
                )
                return None
            return converted

        conversion_matches: List[Dict[str, Any]] = []
        gel_matches: List[Dict[str, Any]] = []
        for stream_name, content in (("stdout", stdout), ("stderr", stderr)):
            for line in content.splitlines():
                conv_match = conversion_pat.match(line)
                if conv_match:
                    converted = _convert_metric(
                        conv_match.group("value"),
                        conv_match.group("unit"),
                        line,
                        "conversion",
                    )
                    if converted is not None:
                        conversion_matches.append(
                            {
                                "value": converted,
                                "line": line.strip(),
                                "stream": stream_name,
                            }
                        )
                    continue

                gel_match = gel_fraction_pat.match(line)
                if gel_match:
                    converted = _convert_metric(
                        gel_match.group("value"),
                        gel_match.group("unit"),
                        line,
                        "gel_fraction",
                    )
                    if converted is not None:
                        gel_matches.append(
                            {
                                "value": converted,
                                "line": line.strip(),
                                "stream": stream_name,
                            }
                        )
                    continue

                line_lower = line.lower()
                if "warning" in line_lower or "warn" in line_lower:
                    qc["warnings"].append(line.strip())

        def _finalize_metric(
            matches: List[Dict[str, Any]],
            metric_key: str,
            metric_label: str,
        ) -> None:
            if not matches:
                return
            # Deterministic rule: use last anchored match.
            qc[metric_key] = matches[-1]["value"]
            unique_vals = sorted({round(float(m["value"]), 8) for m in matches})
            if len(matches) > 1 and len(unique_vals) > 1:
                qc["warnings"].append(
                    f"Conflicting anchored {metric_label} matches detected; using last match."
                )
                for m in matches:
                    qc["warnings"].append(
                        f"{metric_label}_match[{m['stream']}]: {m['line']}"
                    )

        _finalize_metric(conversion_matches, "conversion", "conversion")
        _finalize_metric(gel_matches, "gel_fraction", "gel_fraction")

        if qc["conversion"] is not None or qc["gel_fraction"] is not None:
            qc["status"] = "parsed"

        return qc

    def _validate_curing_steps(self, curing_steps: Any) -> bool:
        """
        Validate optional multi-step curing protocol structure.
        Expected: list of dicts with temperature_K and duration_ps (numbers).
        """
        if not isinstance(curing_steps, list):
            return False
        for step in curing_steps:
            if not isinstance(step, dict):
                return False
            if "temperature_K" not in step or "duration_ps" not in step:
                return False
            try:
                float(step["temperature_K"])
                float(step["duration_ps"])
            except (TypeError, ValueError):
                return False
            if "ramp_ps" in step:
                try:
                    float(step["ramp_ps"])
                except (TypeError, ValueError):
                    return False
        return True

    def _run_gelation_precheck(
        self,
        ctx: "PipelineContext",
        molecule_counts: Dict[str, int],
        reactive_site_counts: Dict[str, int],
    ) -> Tuple[bool, Dict[str, Any]]:
        """
        Conservative gelation feasibility precheck.

        Uses degree moments over reactive species (f_i > 0), weighted by molecule counts.
        """
        skip_precheck = parse_bool(getattr(ctx, "skip_gelation_precheck", None), False)
        allow_fail = parse_bool(getattr(ctx, "allow_gelation_precheck_fail", None), False)
        eps = 1.0e-9

        reactive_species: List[Dict[str, Any]] = []
        for name in sorted(molecule_counts.keys()):
            count = int(molecule_counts.get(name, 0))
            f_i = int(reactive_site_counts.get(name, 0))
            if count > 0 and f_i > 0:
                reactive_species.append(
                    {
                        "name": name,
                        "count": count,
                        "reactive_site_count": f_i,
                    }
                )

        summary: Dict[str, Any] = {
            "skipped": skip_precheck,
            "allow_fail_override": allow_fail,
            "override_used": None,
            "eps": eps,
            "reactive_species": reactive_species,
            "total_reactive_molecules": 0,
            "has_species_f_ge_3": False,
            "mean_f": None,
            "mean_f2": None,
            "mean_f2_minus_mean_f": None,
            "p_c_est": None,
            "fail_reasons": [],
            "pass": False,
            "status": "not_run",
        }

        if skip_precheck:
            summary["pass"] = True
            summary["status"] = "skipped"
            summary["override_used"] = "skip_gelation_precheck"
            print(f"  [WARN] Gelation precheck skipped (--skip-gelation-precheck).")
            return True, summary

        total_n = sum(int(spec["count"]) for spec in reactive_species)
        summary["total_reactive_molecules"] = total_n
        summary["has_species_f_ge_3"] = any(int(spec["reactive_site_count"]) >= 3 for spec in reactive_species)

        if total_n > 0:
            weighted_f = sum(int(spec["count"]) * int(spec["reactive_site_count"]) for spec in reactive_species)
            weighted_f2 = sum(
                int(spec["count"]) * (int(spec["reactive_site_count"]) ** 2)
                for spec in reactive_species
            )
            mean_f = weighted_f / total_n
            mean_f2 = weighted_f2 / total_n
            denom = mean_f2 - mean_f
            summary["mean_f"] = mean_f
            summary["mean_f2"] = mean_f2
            summary["mean_f2_minus_mean_f"] = denom
            if denom > 0.0:
                summary["p_c_est"] = mean_f / denom

        fail_reasons: List[str] = []
        if not summary["has_species_f_ge_3"]:
            fail_reasons.append(
                "No species with reactive_site_count >= 3; network gelation/percolation is infeasible."
            )

        denom_val = summary["mean_f2_minus_mean_f"]
        if denom_val is None:
            fail_reasons.append("No reactive species (reactive_site_count > 0) available for gelation precheck.")
        elif denom_val <= 0.0:
            fail_reasons.append(
                f"mean_f2 - mean_f = {denom_val:.6g} <= 0; gel point criterion cannot be satisfied."
            )

        p_c_est = summary.get("p_c_est")
        if p_c_est is not None and p_c_est > (1.0 + eps):
            fail_reasons.append(
                f"Estimated gel threshold p_c_est={p_c_est:.6g} > 1.0; recipe is conservatively non-gelling."
            )

        summary["fail_reasons"] = fail_reasons
        if fail_reasons:
            if allow_fail:
                summary["pass"] = True
                summary["status"] = "failed_overridden"
                summary["override_used"] = "allow_gelation_precheck_fail"
                print(f"  [WARN] Gelation precheck failed but override is enabled (--allow-gelation-precheck-fail).")
            else:
                summary["pass"] = False
                summary["status"] = "failed"
                return False, summary
        else:
            summary["pass"] = True
            summary["status"] = "passed"

        if summary["mean_f"] is not None and summary["mean_f2"] is not None:
            print(
                f"  - Gelation precheck: mean_f={summary['mean_f']}, "
                f"mean_f2={summary['mean_f2']}, p_c_est={summary['p_c_est']}"
            )
        else:
            print(f"  - Gelation precheck: insufficient reactive species for moment calculation.")
        return True, summary

    def _apply_placeholder_poison_pill(self, top_path: Optional[Path], run_id: str) -> Optional[str]:
        """
        Insert an intentional missing include into placeholder system.top.
        This keeps placeholder files human-readable while mechanically unsafe by default.
        """
        if top_path is None or not top_path.exists():
            return None
        missing_include = f"PLACEHOLDER_POISON_PILL_{run_id}.itp"
        marker = "PLACEHOLDER_POISON_PILL=1"
        content = top_path.read_text()
        if marker in content:
            return missing_include
        poison_header = (
            ";;; PLACEHOLDER_POISON_PILL=1\n"
            ";;; Intentional missing include to block accidental grompp.\n"
            f'#include "{missing_include}"\n\n'
        )
        top_path.write_text(poison_header + content)
        return missing_include
    
    def _find_htpolynet_outputs(
        self,
        workdir: Path,
        output_dir: Path,
    ) -> Tuple[Dict[str, Path], Optional[Path], Optional[Path]]:
        """
        Find HTPolyNet output files with deterministic selection.

        Fix H: Output directory selection robustness:
          - Prefer workdir outputs first (only search output_dir if not found in workdir)
          - De-duplicate identical candidates via resolved realpath
          - Only treat as ambiguous if two DIFFERENT directories both contain valid trio

        Fix C: Zero-byte file checks - require st_size > 0 for all outputs

        Selection order (highest priority first):
          1) workdir/systems/final-results
          2) workdir/systems/final
          3) workdir/output/final
          4) output_dir/ subdirs (only if not found in workdir)
          5) other directories as deterministic fallback

        A directory is valid only if it contains all required files:
        system.gro, system.itp, system.top (index.ndx optional) AND all are non-empty.
        """
        required_names = {
            "gro": "system.gro",
            "itp": "system.itp",
            "top": "system.top",
        }
        optional_ndx_names = ["index.ndx", "system.ndx"]

        priority_relpaths = [
            ("systems/final-results", 0),
            ("systems/final", 1),
            ("output/final", 2),
        ]
        
        # Fix H: Search workdir first, only fall back to output_dir
        search_roots = [workdir, output_dir]
        seen_realpaths: set = set()  # Fix H: Track seen directories by realpath

        def is_valid_output_dir(cand_dir: Path) -> bool:
            """Check if directory has all required non-empty output files."""
            for name in required_names.values():
                fpath = cand_dir / name
                if not fpath.exists():
                    return False
                # Fix C: Require non-zero size
                try:
                    if fpath.stat().st_size == 0:
                        return False
                except OSError:
                    return False
            return True

        def get_output_dict(cand_dir: Path) -> Dict[str, Path]:
            return {k: cand_dir / v for k, v in required_names.items()}

        candidates: List[Tuple[int, Path, Dict[str, Path]]] = []
        for root_idx, root in enumerate(search_roots):
            if not root.exists():
                continue
            for relpath, base_rank in priority_relpaths:
                cand_dir = root / relpath
                if not cand_dir.exists():
                    continue
                # Fix H: Deduplicate by realpath
                try:
                    real = cand_dir.resolve()
                    if real in seen_realpaths:
                        continue
                    seen_realpaths.add(real)
                except OSError:
                    pass
                
                if is_valid_output_dir(cand_dir):
                    # Fix H: Add root index to rank so workdir wins over output_dir
                    rank = base_rank + (root_idx * 10)  # 0-2 for workdir, 10-12 for output_dir
                    outputs = get_output_dict(cand_dir)
                    candidates.append((rank, cand_dir, outputs))

        if candidates:
            min_rank = min(c[0] for c in candidates)
            best = [c for c in candidates if c[0] == min_rank]
            
            # Fix H: Deduplicate best by realpath (shouldn't be needed but safety check)
            unique_best: List[Tuple[int, Path, Dict[str, Path]]] = []
            best_realpaths: set = set()
            for c in best:
                try:
                    real = c[1].resolve()
                    if real not in best_realpaths:
                        best_realpaths.add(real)
                        unique_best.append(c)
                except OSError:
                    unique_best.append(c)
            best = unique_best
            
            if len(best) > 1:
                details = "\n".join(
                    f"  - {c[1]} ({', '.join(sorted(p.name for p in c[2].values()))})"
                    for c in best
                )
                raise HTPolyNetError(
                    f"Multiple valid output directories found at priority {min_rank}:\n{details}\n"
                    f"Clean workdir or remove stale outputs to disambiguate."
                )
            _, selected_dir, outputs = best[0]
        else:
            other_candidates: List[Tuple[Path, Dict[str, Path]]] = []
            for root in search_roots:
                if not root.exists():
                    continue
                for dirpath, _, filenames in os.walk(root):
                    dirpath = Path(dirpath)
                    if not dirpath.is_dir():
                        continue
                    # Fix H: Deduplicate by realpath
                    try:
                        real = dirpath.resolve()
                        if real in seen_realpaths:
                            continue
                        seen_realpaths.add(real)
                    except OSError:
                        pass
                    
                    # Fix C: Check both existence AND non-zero size
                    if is_valid_output_dir(dirpath):
                        outputs = get_output_dict(dirpath)
                        other_candidates.append((dirpath, outputs))
                        
            if not other_candidates:
                raise HTPolyNetError(
                    f"Missing required outputs (system.gro/.itp/.top) under {workdir}.\n"
                    f"HTPolyNet may have failed or produced partial/empty outputs.\n"
                    f"All output files must exist and be non-empty (st_size > 0)."
                )
            if len(other_candidates) > 1:
                details = "\n".join(
                    f"  - {d[0]} ({', '.join(sorted(d[1][k].name for k in d[1]))})"
                    for d in other_candidates
                )
                raise HTPolyNetError(
                    f"Ambiguous output selection: found multiple directories with complete outputs:\n{details}\n"
                    f"Clean the workdir and re-run HTPolyNet."
                )
            selected_dir, outputs = other_candidates[0]

        ndx_path = None
        for name in optional_ndx_names:
            candidate = selected_dir / name
            if candidate.exists():
                ndx_path = candidate
                break

        print(f"  - Selected output dir: {selected_dir}")
        for ext, path in outputs.items():
            file_size = path.stat().st_size
            print(f"  - Found {ext}: {path} ({file_size} bytes)")
        if ndx_path:
            print(f"  - Found ndx: {ndx_path}")
        else:
            print(f"  - No index.ndx found alongside outputs")

        return outputs, selected_dir, ndx_path
    
    def _create_placeholder_outputs(
        self,
        output_dir: Path,
        ctx: "PipelineContext",
    ) -> Dict[str, Path]:
        """
        Create placeholder outputs when htpolynet is not available.
        
        SAFETY: Placeholders are clearly marked as non-physical and include
        FF-aware comb-rule to avoid misleading defaults.
        """
        outputs = {}
        
        # Determine comb-rule based on forcefield (GAFF/AMBER→2, OPLS→3)
        ff_upper = ctx.ff.upper() if hasattr(ctx, 'ff') else 'GAFF2'
        if 'OPLS' in ff_upper:
            comb_rule = 3
            ff_label = 'OPLS-AA'
        else:
            comb_rule = 2
            ff_label = 'GAFF2/AMBER'
        
        # Timestamp for sentinel (ISO 8601)
        timestamp = datetime.datetime.now().isoformat()
        
        # Determine if placeholders will be propagated (affects wording)
        allow_propagate = parse_bool(getattr(ctx, "allow_placeholder_propagate", None), False)
        
        # Placeholder banner with explicit sentinel for downstream detection
        placeholder_banner = (
            f";;; PLACEHOLDER=1\n"
            f";;; run_id={ctx.run_id}\n"
            f";;; timestamp={timestamp}\n"
            f";;; ================================================================\n"
            f";;; WARNING: THIS IS A PLACEHOLDER FILE - NOT PHYSICALLY VALID!\n"
            f";;; ================================================================\n"
            f";;; Generated because HTPolyNet was not available.\n"
            f";;; DO NOT use for production simulations.\n"
        )
        if allow_propagate:
            placeholder_banner += (
                f";;; DEMO MODE: This placeholder WAS propagated to IN/.\n"
                f";;; Results are for pipeline wiring/testing ONLY!\n"
            )
        else:
            placeholder_banner += (
                f";;; This placeholder stays in OUT_GMX only (NOT published to IN/).\n"
            )
        placeholder_banner += (
            f";;; ================================================================\n"
        )
        
        # GRO placeholder with explicit sentinel in title
        gro_path = output_dir / "system.gro"
        gro_path.write_text(
            f"PLACEHOLDER=1 run_id={ctx.run_id} ts={timestamp} ; NOT PHYSICAL - {ctx.system_id}\n"
            "    0\n"
            "   5.00000   5.00000   5.00000\n"
        )
        outputs["gro"] = gro_path
        
        # ITP placeholder
        itp_path = output_dir / "system.itp"
        itp_path.write_text(
            placeholder_banner +
            f"; HTPolyNet crosslinked topology - {ctx.system_id}\n"
            "; Generated by pipeline placeholder\n"
            "\n"
            "[ moleculetype ]\n"
            "; name  nrexcl\n"
            f"{ctx.system_id}  3\n"
            "\n"
            "[ atoms ]\n"
            "; nr  type  resnr  residu  atom  cgnr  charge  mass\n"
        )
        outputs["itp"] = itp_path
        
        # TOP placeholder with FF-aware comb-rule
        top_path = output_dir / "system.top"
        top_path.write_text(
            placeholder_banner +
            f"; HTPolyNet system topology - {ctx.system_id}\n"
            f"; Forcefield: {ff_label} (placeholder comb-rule={comb_rule})\n"
            ";\n"
            "\n"
            "[ defaults ]\n"
            "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n"
            f"1         {comb_rule}          yes        0.5      0.8333\n"
            "\n"
            f'#include "system.itp"\n'
            "\n"
            "[ system ]\n"
            f"; {ctx.system_id} PLACEHOLDER (not crosslinked)\n"
            f"{ctx.system_id}\n"
            "\n"
            "[ molecules ]\n"
            "; name  count\n"
            f"{ctx.system_id}  1\n"
        )
        outputs["top"] = top_path
        
        print(f"  - Created placeholder: system.gro (NON-PHYSICAL)")
        print(f"  - Created placeholder: system.itp (NON-PHYSICAL)")
        print(f"  - Created placeholder: system.top (comb-rule={comb_rule} for {ff_label})")
        
        return outputs
    
    def _publish_outputs(
        self,
        outputs: Dict[str, Path],
        publish_dir: Path,
        run_id: str,
    ) -> List[str]:
        """
        Publish outputs atomically with versioning.
        
        Creates:
        - system_<RUN_ID>.{gro,itp,top} (versioned)
        - system_current.{gro,itp,top} (symlink to versioned)
        """
        published = []
        
        for ext, src_path in outputs.items():
            if not src_path or not src_path.exists():
                continue
            
            dest_dir = publish_dir / ext
            dest_dir.mkdir(parents=True, exist_ok=True)
            
            # Versioned filename
            versioned_name = f"system_{run_id}.{ext}"
            versioned_path = dest_dir / versioned_name
            
            # Atomic copy: write to temp, then rename
            temp_path = self._make_unique_temp_path(dest_dir, versioned_name, run_id)
            try:
                shutil.copy2(src_path, temp_path)
                os.replace(temp_path, versioned_path)
                published.append(str(versioned_path))
            finally:
                _safe_unlink(temp_path)
            
            # Fix G: Atomic symlink update using temp + os.replace pattern
            # This avoids the brief window where the file doesn't exist
            current_name = f"system_current.{ext}"
            current_path = dest_dir / current_name
            temp_current_path = self._make_unique_temp_path(dest_dir, current_name, run_id)
            
            # Create new symlink/copy at temp location first
            try:
                # Try symlink first
                temp_current_path.symlink_to(versioned_name)
                # Atomic replace
                os.replace(temp_current_path, current_path)
                published.append(str(current_path))
            except OSError:
                # Symlinks may not work on all systems, use atomic copy instead
                _safe_unlink(temp_current_path)
                temp_current_copy = self._make_unique_temp_path(
                    dest_dir, f"{current_name}.copy", run_id
                )
                try:
                    shutil.copy2(versioned_path, temp_current_copy)
                    os.replace(temp_current_copy, current_path)
                    published.append(str(current_path))
                finally:
                    _safe_unlink(temp_current_copy)
            finally:
                _safe_unlink(temp_current_path)
        
        return published
    
    def _stage_to_gromacs(
        self,
        outputs: Dict[str, Path],
        gromacs_dir: Path,
        run_id: str,
        ndx_path: Optional[Path] = None,
    ) -> List[str]:
        """
        Stage outputs to GROMACS input pool.
        
        Creates:
        - system.{gro,itp,top} in gromacs/{gro,itp,top}/ directories
        - index.ndx in gromacs/ndx/ if available
        """
        staged = []
        
        for ext, src_path in outputs.items():
            if not src_path or not src_path.exists():
                continue
            
            dest_dir = gromacs_dir / ext
            dest_dir.mkdir(parents=True, exist_ok=True)
            
            # Simple name for GROMACS stage
            dest_name = f"system.{ext}"
            dest_path = dest_dir / dest_name
            
            # Atomic copy
            temp_path = self._make_unique_temp_path(dest_dir, dest_name, run_id)
            try:
                shutil.copy2(src_path, temp_path)
                os.replace(temp_path, dest_path)
                staged.append(str(dest_path))
            finally:
                _safe_unlink(temp_path)

        if ndx_path and ndx_path.exists():
            ndx_dir = gromacs_dir / "ndx"
            ndx_dir.mkdir(parents=True, exist_ok=True)
            dest_name = "index.ndx"
            dest_path = ndx_dir / dest_name
            temp_path = self._make_unique_temp_path(ndx_dir, dest_name, run_id)
            try:
                shutil.copy2(ndx_path, temp_path)
                os.replace(temp_path, dest_path)
                staged.append(str(dest_path))
                print(f"  - Staged index: {dest_path.name}")
            finally:
                _safe_unlink(temp_path)
        else:
            print(f"  - No index.ndx to stage")
        
        return staged
    
    def _parse_molecules_from_top(self, top_path: Path) -> Optional[Dict[str, int]]:
        """
        Parse [molecules] section from a system.top file.
        
        Returns:
            Dict of molecule_name -> count, or None if parsing fails
        """
        if not top_path.exists():
            return None
        
        try:
            content = top_path.read_text()
        except Exception:
            return None
        
        molecules: Dict[str, int] = {}
        in_molecules = False
        
        for line in content.splitlines():
            stripped = line.strip()
            
            # Detect section headers
            if stripped.startswith("[") and stripped.endswith("]"):
                section_name = stripped[1:-1].strip().lower()
                in_molecules = (section_name == "molecules")
                continue
            
            if in_molecules:
                # Skip comments and empty lines
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                
                parts = data.split()
                if len(parts) >= 2:
                    try:
                        mol_name = parts[0]
                        mol_count = int(parts[1])
                        molecules[mol_name] = molecules.get(mol_name, 0) + mol_count
                    except ValueError:
                        continue
        
        return molecules if molecules else None
    
    def _record_manifest(
        self,
        ctx: "PipelineContext",
        config: "HTPolyNetConfig",
        rules_file: Optional[Path],
        rules_dir: Optional[Path] = None,
        rules_source: Optional[str] = None,
        molecule_counts_pre: Optional[Dict[str, int]] = None,
        molecule_counts_post: Optional[Dict[str, int]] = None,
        published_files: Optional[List[str]] = None,
        staged_files: Optional[List[str]] = None,
        dry_run: bool = False,
        placeholder: bool = False,
        placeholder_propagated: bool = False,  # Fix A: Track propagation
        placeholder_poison_pill_enabled: bool = False,
        placeholder_poison_pill_include: Optional[str] = None,
        post_htpolynet_gro: Optional[Path] = None,
        composition_key_used: Optional[str] = None,
        selected_output_dir: Optional[Path] = None,
        ndx_path: Optional[Path] = None,
        box_source: Optional[str] = None,
    ) -> None:
        """
        Record HTPolyNet execution in manifest.
        
        MANIFEST IMMUTABILITY POLICY:
        - composition.initial_counts: IMMUTABLE after PACKMOL stage sets it
        - composition.counts: ALIAS for initial_counts (backward compat shim)
        - composition.post_counts: Post-crosslink state from HTPolyNet [molecules]
        
        We NEVER overwrite initial_counts; we only ADD post_counts.
        """
        if ctx.manifest is None:
            return
        
        htpolynet_data = {
            "project_name": config.project_name,
            "forcefield": config.forcefield,
            "initial_structure": str(config.initial_gro),
            "workdir": str(config.workdir),
            "config_file": str(config.output_dir / "system_config.yaml"),
            "rules_file": str(rules_file) if rules_file else None,
            "rules_dir": str(rules_dir) if rules_dir else None,
            "rules_source": rules_source,
            "reactive_site_counts": dict(sorted(config.reactive_site_counts.items())),
            "crosslinkers": dict(sorted(config.crosslinkers.items())),
            "chain_stoppers": dict(sorted(config.chain_stoppers.items())),
            "max_conversion": config.max_conversion,
            "dry_run": dry_run,
            "placeholder": placeholder,
            "placeholder_propagated": bool(placeholder and placeholder_propagated),
            "placeholder_propagated_downstream_unsafe": bool(placeholder and placeholder_propagated),
            "placeholder_poison_pill_enabled": bool(placeholder and placeholder_poison_pill_enabled),
        }
        if config.curing_steps:
            htpolynet_data["curing_steps"] = config.curing_steps
        
        # Fix A: Add placeholder and propagation tracking
        if placeholder:
            htpolynet_data["placeholder_warning"] = (
                "This is a PLACEHOLDER run. Outputs are NOT physically valid. "
                "Do not use for production simulations or analysis."
            )
            if placeholder_poison_pill_include:
                htpolynet_data["placeholder_poison_pill_include"] = placeholder_poison_pill_include
            if placeholder_propagated:
                htpolynet_data["placeholder_propagation_warning"] = (
                    "DEMO ONLY: Placeholders were propagated to IN/. "
                    "Results are NOT physically meaningful!"
                )
        if box_source:
            htpolynet_data["box_source"] = box_source
        if composition_key_used:
            htpolynet_data["composition_key_used"] = composition_key_used
        if selected_output_dir:
            htpolynet_data["selected_output_dir"] = str(selected_output_dir)
        if ndx_path:
            htpolynet_data["ndx_path"] = str(ndx_path)
        
        # === MANIFEST IMMUTABILITY: initial_counts preserved, post_counts added ===
        composition = ctx.manifest._data.setdefault("composition", {})
        
        # Migrate legacy manifests: if only 'counts' exists, treat as initial_counts
        if "counts" in composition and "initial_counts" not in composition:
            composition["initial_counts"] = composition["counts"]
        
        # Record pre-reaction composition (should already be set by PACKMOL)
        if molecule_counts_pre:
            htpolynet_data["composition_pre_reaction"] = molecule_counts_pre
            # Ensure initial_counts is set (immutable after first set)
            if "initial_counts" not in composition:
                composition["initial_counts"] = molecule_counts_pre
                composition["counts"] = molecule_counts_pre  # backward compat shim
                composition["molecule_counts"] = molecule_counts_pre  # legacy alias
        
        # Record post-reaction composition in SEPARATE field (do NOT overwrite counts)
        if molecule_counts_post:
            htpolynet_data["composition_after_reaction"] = molecule_counts_post
            # Store in composition.post_counts (new field)
            composition["post_counts"] = molecule_counts_post
            composition["post_counts_source"] = "htpolynet_system_top"
            composition["post_counts_timestamp"] = datetime.datetime.now().isoformat()
            # DO NOT update composition.counts - it stays as initial_counts
        
        if published_files is not None:
            htpolynet_data["published_files"] = published_files
        if staged_files is not None:
            htpolynet_data["staged_files"] = staged_files
        
        # === Record post-HTPolyNet box vectors if available ===
        if post_htpolynet_gro and post_htpolynet_gro.exists():
            post_box, post_triclinic, post_warnings = self._read_box_from_gro(post_htpolynet_gro)
            if post_box:
                box_data = ctx.manifest._data.setdefault("box", {})
                if post_triclinic:
                    box_data["post_htpolynet_box_vectors_nm"] = [list(row) for row in post_box]
                    box_data["post_htpolynet_box_is_triclinic"] = True
                else:
                    box_data["post_htpolynet_box_vectors_nm"] = list(post_box)
                    box_data["post_htpolynet_box_is_triclinic"] = False
                box_data["post_htpolynet_box_source"] = str(post_htpolynet_gro)
            for warn in post_warnings:
                htpolynet_data.setdefault("post_htpolynet_box_warnings", []).append(warn)
        
        # Merge with existing htpolynet data (preserve existing keys, don't overwrite with None)
        existing_htpolynet = ctx.manifest._data.get("htpolynet", {})
        for key in (
            "placeholder_warning",
            "placeholder_propagation_warning",
            "placeholder_poison_pill_include",
            "placeholder_poison_pill_enabled",
            "placeholder_propagated",
            "placeholder_propagated_downstream_unsafe",
        ):
            existing_htpolynet.pop(key, None)
        existing_htpolynet.update(_compact_dict(htpolynet_data))
        ctx.manifest._data["htpolynet"] = existing_htpolynet
        ctx.manifest.save()
