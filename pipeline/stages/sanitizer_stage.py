"""
Sanitizer stage orchestrator extracted from sanitizer.py.
"""

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple
import os

from .base import BaseStage
from ..itp_sanitizer import IncludeResolution
from ..gromacs_cmd import resolve_gmx_command
from .topology_sanitizer import (
    SanitizerError,
    TopologySanitizerMixin,
    parse_section_name,
    normalize_include_priority,
    DEFAULT_INCLUDE_PRIORITY,
    DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE,
)
from .spatial_checker import SpatialCheckerMixin

if TYPE_CHECKING:
    from ..context import PipelineContext


class SanitizerStage(TopologySanitizerMixin, SpatialCheckerMixin, BaseStage):
    """
    Stage 2.5: ITP Sanitizer (between htpolynet and gromacs)

    Reads:
    - IN/systems/<SYSTEM_ID>/gromacs/itp/*.itp (staged ITPs)
    - IN/systems/<SYSTEM_ID>/htpolynet/itp/*.itp (htpolynet ITPs)
    - IN/molecules/<NAME>/itp/*.itp (referenced molecule ITPs)
    - IN/forcefield/<FF>/gromacs/*.itp (forcefield atomtypes)

    Writes:
    - IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_<RUN_ID>.itp
    - IN/systems/<SYSTEM_ID>/gromacs/combined_atomtypes_current.itp (link)
    - IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_<RUN_ID>/
    - IN/systems/<SYSTEM_ID>/gromacs/itp_sanitized_current/ (link)
    - IN/systems/<SYSTEM_ID>/gromacs/top/system.top (updated)
    """
    @property
    def name(self) -> str:
        return "sanitizer"
    
    @property
    def output_subdir(self) -> str:
        return "02_5_sanitizer"

    @staticmethod
    def _path_identity(path: Path) -> str:
        """Stable path identity for dedup/map keys."""
        try:
            return str(path.resolve(strict=False))
        except OSError:
            return str(path.absolute())

    @classmethod
    def _manifest_safe_value(cls, value: Any) -> Any:
        """Convert values to manifest-safe JSON-like structures."""
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {
                str(k): cls._manifest_safe_value(v)
                for k, v in value.items()
            }
        if isinstance(value, (list, tuple, set)):
            return [cls._manifest_safe_value(v) for v in value]
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        return str(value)

    @classmethod
    def _serialize_manifest_entry(cls, obj: Any) -> Dict[str, Any]:
        """Safely serialize include-resolution/diagnostic entries."""
        if is_dataclass(obj):
            raw = asdict(obj)
        elif isinstance(obj, dict):
            raw = obj
        elif hasattr(obj, "__dict__"):
            raw = {
                key: value
                for key, value in vars(obj).items()
                if not key.startswith("_")
            }
        else:
            return {"value": cls._manifest_safe_value(obj)}
        if not isinstance(raw, dict):
            return {"value": cls._manifest_safe_value(raw)}
        return {
            str(key): cls._manifest_safe_value(value)
            for key, value in raw.items()
        }
    
    def run(self, ctx: "PipelineContext") -> bool:
        """
        Execute ITP sanitization and charge neutrality preflight.
        
        Returns:
            True if successful, False otherwise
        """
        from ..itp_sanitizer import (
            ItpSanitizer,
            ItpSanitizerError,
            generate_system_top,
        )
        from ..charge_neutrality import (
            ChargeNeutralityChecker,
            ChargeNeutralityError,
        )
        
        print("  Running ITP Sanitizer stage...")
        self.diagnostics: List[IncludeResolution] = []
        # Reset per-run cached paths for safe instance reuse.
        self._forcefield_dir = None
        self._sanitized_current_dir = None
        
        # Get paths
        system_gromacs_dir = ctx.get_input_path(
            "systems", ctx.system_id, "gromacs"
        )
        system_htpolynet_dir = ctx.get_input_path(
            "systems", ctx.system_id, "htpolynet"
        )
        molecules_dir = ctx.get_input_path("molecules")
        
        # Task F: Flexible forcefield directory resolution
        # Try both dashed and non-dashed directory names
        ff_base = ctx.ff.lower()
        ff_candidates = [
            ctx.get_input_path("forcefield", ff_base.replace("-", "")),
            ctx.get_input_path("forcefield", ff_base),
        ]
        forcefield_dir = None
        for ff_candidate in ff_candidates:
            if ff_candidate.exists():
                forcefield_dir = ff_candidate
                break
        if forcefield_dir is None:
            raise SanitizerError(
                f"Forcefield directory not found. Tried:\n"
                + "\n".join(f"  - {c}" for c in ff_candidates)
            )
        self._forcefield_dir = forcefield_dir
        
        print(f"  - System GROMACS dir: {system_gromacs_dir}")
        print(f"  - Forcefield: {ctx.ff} (dir: {forcefield_dir.name})")
        
        # Determine molecule counts early (ordered if possible)
        top_path = system_gromacs_dir / "top" / "system.top"
        top_parse = self._get_ordered_molecules_from_top(top_path)
        ordered_molecules = top_parse.ordered
        top_exists = top_path.exists()
        output_dir = self.get_output_dir(ctx)
        molecule_count_source: Dict[str, Any] = {
            "raw_top_path": str(top_path),
            "raw_top_exists": bool(top_exists),
            "raw_top_uncertain": bool(top_parse.uncertain),
            "raw_top_uncertainty_reasons": list(top_parse.uncertainty_reasons),
            "source": None,
            "fallback_used": False,
            "trusted": False,
        }

        if top_exists and (top_parse.uncertain or not ordered_molecules):
            if top_parse.uncertain:
                print(
                    "  [WARN] Raw [molecules] parse is uncertain due to conditional directives. "
                    "Requiring grompp -pp preprocessed topology for trusted counts."
                )
            else:
                print(
                    "  [WARN] system.top exists but [molecules] could not be parsed. "
                    "Requiring grompp -pp preprocessed topology for trusted counts."
                )
            gro_for_pp = self._find_coords_gro(system_gromacs_dir)
            if gro_for_pp is None:
                raise SanitizerError(
                    "Cannot derive trusted [molecules] counts: system.top is present but raw parsing "
                    "is not trustworthy and no coordinate GRO is available for grompp -pp preprocessing."
                )
            try:
                self._expected_atoms_from_grompp(
                    gmx_exe=resolve_gmx_command(ctx),
                    system_top_path=top_path,
                    gro_path=gro_for_pp,
                    output_dir=output_dir,
                    ctx=ctx,
                )
            except SanitizerError as exc:
                print(
                    "  [WARN] Unable to preprocess topology for trusted [molecules] extraction: "
                    f"{exc}"
                )
            pp_path = self._grompp_preprocessed_top_path(output_dir)
            preprocessed_molecules = self._parse_preprocessed_molecules(pp_path)
            molecule_count_source["preprocessed_top_path"] = str(pp_path)
            if preprocessed_molecules:
                ordered_molecules = preprocessed_molecules
                molecule_counts = {name: count for name, count in ordered_molecules}
                molecule_count_source["source"] = "grompp_preprocessed_top"
                molecule_count_source["raw_top_used"] = False
                molecule_count_source["trusted"] = True
                print("  Using [molecules] from grompp -pp preprocessed topology")
            else:
                raise SanitizerError(
                    "Cannot derive trusted [molecules] counts for staged topology.\n"
                    f"  system.top: {top_path}\n"
                    f"  preprocessed_top: {pp_path}\n"
                    "Raw [molecules] parsing is uncertain or empty, and grompp -pp did not yield "
                    "trusted counts. Refusing fallback to manifest pre-reaction counts."
                )
        elif ordered_molecules:
            molecule_counts = {name: count for name, count in ordered_molecules}
            molecule_count_source["source"] = "raw_system_top"
            molecule_count_source["raw_top_used"] = True
            molecule_count_source["trusted"] = True
            print("  Using [molecules] from HTPolyNet-staged system.top")
        else:
            molecule_counts = self._get_molecule_counts(ctx)
            ordered_molecules = list(molecule_counts.items())
            molecule_count_source["source"] = "manifest"
            molecule_count_source["fallback_used"] = True
            molecule_count_source["raw_top_used"] = False
            molecule_count_source["trusted"] = bool(ordered_molecules)
            print("  Using [molecules] from manifest (pre-reaction counts)")
        self._molecule_count_source = molecule_count_source
        
        needed_molecule_names = {
            name for name, count in ordered_molecules if count > 0
        }
        
        # Collect all ITP files with proper classification
        itp_paths: List[Path] = []
        ff_map: Dict[str, str] = {}  # file_path -> forcefield identifier
        class_map: Dict[str, str] = {}  # file_path -> file class
        
        # 1. Forcefield atomtypes (if exists) - highest priority
        ff_gromacs_dir = forcefield_dir / "gromacs"
        if ff_gromacs_dir.exists():
            ff_itps = self._sorted_paths_casefold(ff_gromacs_dir.glob("*.itp"))
            for itp in ff_itps:
                itp_key = self._path_identity(itp)
                if itp_key in ff_map:
                    continue
                itp_paths.append(itp)
                ff_map[itp_key] = ctx.ff  # Only FF files get ctx.ff
                class_map[itp_key] = "forcefield"
                print(f"    + FF: {itp.name}")
        
        # 2. Staged system ITPs from gromacs/itp/
        staged_itp_dir = system_gromacs_dir / "itp"
        if staged_itp_dir.exists():
            staged_itps = self._sorted_paths_casefold(staged_itp_dir.glob("*.itp"))
            for itp in staged_itps:
                itp_key = self._path_identity(itp)
                if itp_key not in ff_map:  # Avoid duplicates
                    itp_paths.append(itp)
                    ff_map[itp_key] = "UNKNOWN"  # Not from forcefield
                    class_map[itp_key] = "staged"
                    print(f"    + Staged: {itp.name}")
        
        # 3. HTPolyNet produced ITPs
        htpolynet_itp_dir = system_htpolynet_dir / "itp"
        if htpolynet_itp_dir.exists():
            htpolynet_itps = self._sorted_paths_casefold(htpolynet_itp_dir.glob("*.itp"))
            for itp in htpolynet_itps:
                itp_key = self._path_identity(itp)
                if itp_key not in ff_map:  # Avoid duplicates
                    itp_paths.append(itp)
                    ff_map[itp_key] = "UNKNOWN"  # From HTPolyNet, not FF library
                    class_map[itp_key] = "htpolynet"
                    print(f"    + HTPolyNet: {itp.name}")
        
        # 4. Molecule ITPs (molecules library)
        needed_name_map = {n.casefold(): n for n in needed_molecule_names}
        if molecules_dir.exists():
            mol_dirs = sorted(molecules_dir.iterdir(), key=lambda p: p.name.casefold())
            for mol_dir in mol_dirs:
                if not mol_dir.is_dir():
                    continue
                mol_itp_dir = mol_dir / "itp"
                if mol_itp_dir.exists():
                    mol_itps = self._sorted_paths_casefold(mol_itp_dir.glob("*.itp"))
                    for itp in mol_itps:
                        itp_key = self._path_identity(itp)
                        if itp_key in ff_map:
                            continue
                        mol_names = self._get_moleculetype_names(itp, ctx=ctx)
                        matched_names = sorted(
                            {
                                needed_name_map[m.casefold()]
                                for m in mol_names
                                if m.casefold() in needed_name_map
                            },
                            key=lambda n: n.casefold(),
                        )
                        if not matched_names:
                            continue
                        if itp_key not in ff_map:
                            itp_paths.append(itp)
                            ff_map[itp_key] = "UNKNOWN"
                            class_map[itp_key] = "molecule"
                            print(
                                f"    + Molecules-lib {','.join(matched_names)}: {itp.name}"
                            )
        
        # Hardening v7: Expand include-closure for ALL ITPs (not just forcefield)
        # This ensures moleculetypes in included files are discovered
        print("  Expanding include-closure for all ITPs...")
        itp_paths, class_map, ff_map = self._expand_itp_include_closure(
            itp_paths, ctx, class_map, ff_map
        )
        
        # Build exact-name moleculetype candidates from all sources
        all_candidates: Dict[str, List[Path]] = {}
        casefold_index: Dict[str, Set[str]] = {}  # casefold -> set of exact names
        
        for itp in itp_paths:
            for mol_name in self._get_moleculetype_names(itp, ctx=ctx):
                all_candidates.setdefault(mol_name, []).append(itp)
                casefold_index.setdefault(mol_name.casefold(), set()).add(mol_name)
        
        # Fail-fast on case collisions (different exact names that casefold the same)
        for cf_key, exact_names in casefold_index.items():
            if len(exact_names) > 1:
                raise SanitizerError(
                    f"Case-insensitive moleculetype collision will break grompp:\n"
                    f"  Names: {sorted(exact_names)}\n"
                    f"  Rename one of these moleculetypes to avoid OS-dependent behavior."
                )

        # Select one deterministic winner per exact moleculetype name.
        all_moleculetypes: Dict[str, Path] = {}
        for mol_name in sorted(all_candidates, key=lambda n: (n.casefold(), n)):
            winner = self._select_moleculetype_winner(
                mol_name=mol_name,
                candidate_paths=all_candidates[mol_name],
                class_map=class_map,
                ctx=ctx,
            )
            all_moleculetypes[mol_name] = winner
        
        # Map needed molecules from the single winner map above (no contradictory paths).
        molecule_itp_paths: Dict[str, Path] = {}
        if needed_molecule_names:
            still_missing = []
            for n in sorted(needed_molecule_names, key=lambda x: x.casefold()):
                exact_names = sorted(
                    casefold_index.get(n.casefold(), set()),
                    key=lambda x: (x.casefold(), x),
                )
                if not exact_names:
                    still_missing.append(n)
                    continue
                winner = all_moleculetypes.get(exact_names[0])
                if winner is None:
                    still_missing.append(n)
                    continue
                molecule_itp_paths[n] = winner
            
            if still_missing:
                missing_sorted = sorted(still_missing, key=lambda x: x.casefold())
                # Count ITPs per source category
                ff_count = sum(1 for p, c in class_map.items() if c == "forcefield")
                staged_count = sum(1 for p, c in class_map.items() if c == "staged")
                htpoly_count = sum(1 for p, c in class_map.items() if c == "htpolynet")
                mol_count = sum(1 for p, c in class_map.items() if c == "molecule")
                
                raise SanitizerError(
                    f"Missing moleculetype(s) for [molecules]: {missing_sorted}\n"
                    f"Searched in:\n"
                    f"  - Forcefield: {forcefield_dir}/gromacs/ ({ff_count} files)\n"
                    f"  - Staged: {system_gromacs_dir}/itp/ ({staged_count} files)\n"
                    f"  - HTPolyNet: {system_htpolynet_dir}/itp/ ({htpoly_count} files)\n"
                    f"  - Molecules: {molecules_dir}/ ({mol_count} files)\n"
                    f"Suggestions:\n"
                    f"  - Verify HTPolyNet stage outputs\n"
                    f"  - Check #include paths in topology\n"
                    f"  - Run `grompp -pp` to see full preprocessed topology"
                )
        
        if not itp_paths:
            print("  [WARN] No ITP files found to sanitize")
            print("  Creating empty combined atomtypes file...")
            
            # Create minimal outputs
            self._create_minimal_outputs(ctx, system_gromacs_dir, allow_default_defaults=ctx.allow_default_defaults)
            return True
        
        print(f"  Found {len(itp_paths)} ITP files to process")
        
        # Get options from context
        allow_override = getattr(ctx, 'allow_override', False)
        atomtype_prefix = getattr(ctx, 'atomtype_prefix', None)
        strict_charge = getattr(ctx, 'strict_charge', False)
        
        # Initialize sanitizer
        sanitizer = ItpSanitizer(
            source_ff=ctx.ff,
            allow_override=allow_override,
            atomtype_prefix=atomtype_prefix,
            strict_charge=strict_charge,
        )
        mixed_defaults_warnings: List[str] = []
        
        try:
            # Run sanitization
            include_roots = [top_path] if top_path.exists() else []
            result = sanitizer.run(
                itp_paths=itp_paths,
                output_dir=system_gromacs_dir,
                run_id=ctx.run_id,
                ff_map=ff_map,
                class_map=class_map,
                ctx=ctx,
                include_roots=include_roots,
                allow_default_defaults=ctx.allow_default_defaults,
            )
            self._sanitized_current_dir = result.sanitized_current_dir
            
            print(f"  [OK] Sanitization complete:")
            print(f"       Atomtypes: {result.atomtype_count}")
            print(f"       Files sanitized: {len(result.sanitized_files)}")
            print(f"       Molecule ITPs: {len(result.molecule_itp_files)}")
            
            if result.conflicts_overridden:
                print(f"  [WARN] {len(result.conflicts_overridden)} conflicts overridden")
                for conflict in result.conflicts_overridden:
                    winner_file = conflict.winner.source_files[0] if conflict.winner and conflict.winner.source_files else "?"
                    print(f"         - {conflict.atomtype_name}: winner = {Path(winner_file).name}")
            
            if result.charge_warnings:
                print(f"  [WARN] {len(result.charge_warnings)} charge warnings (non-fatal)")
            
            if result.defaults_used:
                print(f"  [ defaults ] comb-rule: {result.defaults_used.comb_rule}")
            if getattr(result, "mixed_defaults_detected", False):
                defaults_report = getattr(result, "mixed_defaults_report", {})
                mixed_defaults_warnings.append(
                    "Mixed [defaults] tuples detected; this is an unsafe topology mode."
                )
                mixed_defaults_warnings.append(
                    "Mixed defaults classification="
                    f"{defaults_report.get('classification', 'unknown')}, "
                    f"differing_fields={defaults_report.get('differing_fields', [])}"
                )
                mixed_defaults_warnings.append(
                    "Mixed defaults chosen_policy="
                    f"{defaults_report.get('chosen_policy', 'unknown')}, "
                    f"primary_signature={defaults_report.get('primary_defaults_signature', 'n/a')}"
                )
                for warning in mixed_defaults_warnings:
                    print(f"  [WARN] {warning}")
                if getattr(ctx, "allow_mixed_defaults", False):
                    print("  [WARN] Continuing due to explicit unsafe override (--allow-mixed-defaults).")
                    print(
                        "  [WARN] Override reason: "
                        f"{getattr(ctx, 'allow_mixed_defaults_reason', '')}"
                    )
            if getattr(result, "nonbond_params_secondary_current_path", None):
                print(
                    "  [WARN] Added sanitizer-generated within-group preserve include: "
                    f"{result.nonbond_params_secondary_current_path}"
                )
            if getattr(result, "nonbond_params_cross_group_current_path", None):
                print(
                    "  [WARN] Added sanitizer-generated cross-group include: "
                    f"{result.nonbond_params_cross_group_current_path}"
                )
            prefix_summary = getattr(result, "prefix_policy_summary", {})
            if prefix_summary:
                print(
                    "  Prefix policy: "
                    f"{prefix_summary.get('policy')} "
                    f"-> {prefix_summary.get('applied_strategy')}"
                )
            if getattr(result, "prefix_injected_types_current_path", None):
                print(
                    "  [WARN] Added sanitizer-generated prefixed [*types] include: "
                    f"{result.prefix_injected_types_current_path}"
                )
            
        except ItpSanitizerError as e:
            print(f"  [FAIL] Sanitization failed: {e}")
            raise SanitizerError(str(e))

        if result.defaults_used is None:
            raise SanitizerError(
                "No [ defaults ] block found in inputs. Refusing to invent defaults.\n"
                "Provide explicit defaults in input ITPs or enable --allow-default-defaults."
            )
        if getattr(result, "defaults_guessed", False):
            print("  [WARN] Using fallback [ defaults ] per ctx.allow_default_defaults")
        if bool(getattr(ctx, "allow_mixed_defaults", False)) and bool(
            getattr(result, "mixed_defaults_detected", False)
        ):
            cross_group_summary = getattr(result, "nonbond_params_cross_group_summary", {})
            policy = (
                cross_group_summary.get("policy")
                if isinstance(cross_group_summary, dict)
                else None
            )
            status = (
                cross_group_summary.get("status")
                if isinstance(cross_group_summary, dict)
                else None
            )
            cross_group_current = getattr(
                result, "nonbond_params_cross_group_current_path", None
            )
            if policy == "generate":
                if not cross_group_current:
                    raise SanitizerError(
                        "Mixed [defaults] override requested cross-group remediation (policy=generate) "
                        "but no nonbond_params cross-group include was produced."
                    )
                cross_group_path = Path(cross_group_current)
                if not cross_group_path.exists():
                    raise SanitizerError(
                        "Mixed [defaults] cross-group remediation include path is missing on disk.\n"
                        f"  expected: {cross_group_path}\n"
                        f"  status: {status}"
                    )
        
        cross_ff_conflicts = [
            c for c in result.conflicts_detected
            if c.conflict_type == "cross_ff_collision"
        ]
        if cross_ff_conflicts:
            msg = (
                f"Detected {len(cross_ff_conflicts)} cross-forcefield atomtype conflicts. "
                "This suggests mixed forcefield profiles."
            )
            if ctx.strict_forcefield_consistency:
                raise SanitizerError(msg)
            print(f"  [WARN] {msg}")
        
        # v4: Pass comb_rule for comb-rule-aware LJ outlier detection
        comb_rule = result.defaults_used.comb_rule if result.defaults_used else None
        outlier_warnings = self._warn_atomtype_outliers(result.combined_atomtypes_current_path, ctx, comb_rule)
        if outlier_warnings:
            print(f"  [WARN] LJ outlier warnings: {len(outlier_warnings)}")
        
        # v4 Issue 3: Invoke gen-pairs/[pairs] consistency check
        # Scan all source ITPs (not just sanitized_files) to catch FF/staged/htpolynet-derived pairs
        if result.defaults_used:
            self._warn_pairs_consistency(
                result.defaults_used,
                itp_paths,
                source_defaults_map=getattr(result, "source_defaults_map", None),
                ctx=ctx,
            )
        
        # === Charge Neutrality Check ===
        
        # Fail fast if no molecule counts available
        if not molecule_counts:
            print("  [ERROR] No molecule counts available from HTPolyNet or manifest")
            print("         Run HTPolyNet stage first or ensure manifest has composition data")
            raise SanitizerError("No molecule counts available for charge neutrality check")
        
        sanitized_molecule_itp_paths = self._map_sanitized_itps_to_molecules(
            result.sanitized_current_dir,
            needed_molecule_names,
            ctx=ctx,
        )
        active_molecule_names = sorted(needed_molecule_names, key=str.casefold)
        if active_molecule_names and not sanitized_molecule_itp_paths:
            if not result.sanitized_current_dir.exists():
                raise SanitizerError(
                    "Charge neutrality check cannot run: active molecules are present but "
                    "sanitized ITP output directory is missing.\n"
                    f"  active_molecules: {active_molecule_names}\n"
                    f"  expected_dir: {result.sanitized_current_dir}"
                )
            raise SanitizerError(
                "Charge neutrality check cannot run: active molecules are present but "
                "sanitized molecule-to-ITP mapping is empty.\n"
                f"  active_molecules: {active_molecule_names}\n"
                f"  sanitized_dir: {result.sanitized_current_dir}\n"
                "Ensure sanitizer outputs include moleculetype definitions for all active molecules."
            )
        
        if sanitized_molecule_itp_paths:
            print("  Running charge neutrality preflight...")
            
            # Parse target allowlist from comma-separated string
            target_allowlist = self._parse_csv_casefold_set(ctx.charge_fix_target_allowlist)
            
            # Part 2: Parse new comma-separated settings
            allowed_atomnames = None
            if ctx.charge_fix_target_atomnames:
                allowed_atomnames = set(
                    s.strip() for s in ctx.charge_fix_target_atomnames.split(",") if s.strip()
                )
            
            disallow_atomtypes = None
            if ctx.charge_fix_disallowed_atomtypes:
                disallow_atomtypes = set(
                    s.strip() for s in ctx.charge_fix_disallowed_atomtypes.split(",") if s.strip()
                )
            
            # Hardening v5 - Issue G: Pre-classify ionic moleculetypes
            ion_classifications, ion_classification_audit = self._classify_ionic_moleculetypes(
                sanitized_molecule_itp_paths, ctx
            )
            protected_classes = {
                "protected_ion",
                "ionic",
                "monatomic_ion",
                "protected_solvent",
                "protected_polymer",
                "protected_configured",
            }
            protected_detected = any(
                cls in protected_classes for cls in ion_classifications.values()
            )
            if ctx.manifest:
                ctx.manifest.set_sanitizer_output(
                    "charge_classification",
                    ion_classification_audit,
                )
            effective_target_allowlist = self._derive_charge_fix_checker_allowlist(
                ion_classifications,
                ion_classification_audit,
                target_allowlist,
            )
            ion_classification_audit["checker_target_allowlist"] = sorted(
                effective_target_allowlist,
                key=str.casefold,
            )

            # Optional explicit correction target(s) from CLI. Checker supports one target;
            # we resolve deterministically and fail-fast on ambiguous strict requests.
            requested_target: Optional[str] = None
            requested_targets_raw = self._parse_csv_casefold_set(
                getattr(ctx, "charge_fix_target_molecules", None)
            )
            if requested_targets_raw:
                available_by_casefold = {
                    name.casefold(): name for name in sanitized_molecule_itp_paths.keys()
                }
                resolved_targets: List[str] = []
                unknown_targets: List[str] = []
                for token_cf in sorted(requested_targets_raw):
                    resolved = available_by_casefold.get(token_cf)
                    if resolved is None:
                        unknown_targets.append(token_cf)
                        continue
                    resolved_targets.append(resolved)
                if unknown_targets:
                    msg = (
                        "Unknown --charge-fix-target-molecules entries: "
                        f"{unknown_targets}. Available molecules: "
                        f"{sorted(sanitized_molecule_itp_paths.keys(), key=str.casefold)}"
                    )
                    if ctx.strict_charge_neutrality:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                if len(resolved_targets) > 1:
                    msg = (
                        "Multiple explicit charge-fix targets provided; only one target is supported "
                        f"per correction run: {resolved_targets}"
                    )
                    if ctx.strict_charge_neutrality:
                        raise SanitizerError(msg)
                    resolved_targets = sorted(set(resolved_targets), key=str.casefold)
                    requested_target = resolved_targets[0]
                    print(
                        f"  [WARN] {msg}. Using deterministic first target: '{requested_target}'."
                    )
                elif len(resolved_targets) == 1:
                    requested_target = resolved_targets[0]

            # Strict safety hardening: when protected species exist, require explicit target allowlist.
            if ctx.strict_charge_neutrality and protected_detected and target_allowlist is None:
                raise SanitizerError(
                    "Strict charge-neutrality mode requires an explicit "
                    "--charge-fix-target-allowlist when protected ions/solvents/polymers are present."
                )
            
            # Task C: Coalesce None to safe defaults at point of use
            max_delta_limit = (
                ctx.charge_fix_max_delta_per_atom
                if ctx.charge_fix_max_delta_per_atom is not None
                else 1e-5
            )
            max_total_limit = (
                ctx.charge_fix_max_total
                if ctx.charge_fix_max_total is not None
                else 1e-4
            )
            dipole_limit = (
                ctx.charge_fix_max_dipole_shift_debye
                if ctx.charge_fix_max_dipole_shift_debye is not None
                else 0.01
            )
            
            try:
                # Part 2: Use new context fields for threshold and correction limits
                checker = ChargeNeutralityChecker(
                    threshold=ctx.charge_neutrality_correct,  # Max |Q| for correction
                    max_delta_per_atom=max_delta_limit,
                    max_total_delta=max_total_limit,
                    max_dipole_shift_debye=dipole_limit,
                    strict_charge_physics=ctx.strict_charge_physics,
                    # Safe-subset config
                    correction_method=getattr(ctx, "charge_fix_method", "safe_subset"),
                    target_allowlist=effective_target_allowlist,
                    allowed_atomnames=allowed_atomnames,
                    disallow_atomtypes=disallow_atomtypes,
                    # Part 2: New strict mode and threshold options
                    strict_mode=ctx.strict_charge_neutrality,
                    allow_ions=ctx.charge_fix_allow_ions,
                    warn_threshold=ctx.charge_neutrality_warn,
                    neutrality_tol=ctx.charge_neutrality_tol,
                    protected_polymers={
                        name for name, cls in ion_classifications.items() if cls == "protected_polymer"
                    },
                    polymer_net_charge_tol=getattr(ctx, "polymer_net_charge_tol", 1e-3),
                    polymer_correction_method=getattr(ctx, "charge_fix_polymer_method", "skip_if_small"),
                    polymer_exclusion_bonds=getattr(ctx, "charge_fix_polymer_exclusion_bonds", 2),
                    rounding_moleculetype_tol=getattr(ctx, "charge_fix_moleculetype_rounding_tol", 1e-4),
                    allow_non_rounding_correction=getattr(ctx, "charge_fix_allow_non_rounding", False),
                    allow_unparsed_atoms_lines=getattr(ctx, "allow_unparsed_atoms_lines", False),
                )
                correction_result = checker.run(
                    molecule_itp_paths=sanitized_molecule_itp_paths,
                    molecule_counts=molecule_counts,
                    output_dir=result.sanitized_current_dir,
                    run_id=ctx.run_id,
                    target_molecule=requested_target,
                )
                
                charge_summary = {
                    "computed_q": correction_result.original_q,
                    "corrected": correction_result.corrected,
                    "method": correction_result.correction_method,
                    "effective_method": correction_result.effective_method,
                    "method_desc": correction_result.method_desc,
                    "patched_itp": str(correction_result.patched_itp_path)
                        if correction_result.patched_itp_path else None,
                    "max_delta_applied": correction_result.max_delta_applied,
                    "total_delta": correction_result.total_delta,
                    "atoms_adjusted": correction_result.atoms_adjusted,
                    "target_molecule": correction_result.target_molecule,
                    "target_itp_moleculetype_name": correction_result.target_itp_moleculetype_name,
                    "patched_atoms": correction_result.patched_atoms,
                    # New fields for enhanced manifest logging (Issues A/B)
                    "safe_subset_used": correction_result.safe_subset_used,
                    "patched_atom_details": correction_result.patched_atom_details,
                    "patched_atom_details_truncated": correction_result.patched_atom_details_truncated,
                    "correction_refused": correction_result.correction_refused,
                    "refuse_reason": correction_result.refuse_reason,
                    "unparsed_atoms_lines": correction_result.unparsed_atoms_lines,
                    "unparsed_atoms_lines_truncated": correction_result.unparsed_atoms_lines_truncated,
                    "polymer_policy": correction_result.polymer_policy,
                    "rounding_only_passed": correction_result.rounding_only_passed,
                    "rounding_moleculetype_tol": correction_result.rounding_moleculetype_tol,
                    "per_moleculetype_report": correction_result.per_moleculetype_report,
                    "dominant_charge_contributors": correction_result.dominant_charge_contributors,
                    "hetero_atoms_modified": correction_result.hetero_atoms_modified,
                    "hetero_atoms_modified_count": correction_result.hetero_atoms_modified_count,
                    "hetero_atoms_modified_details": correction_result.hetero_atoms_modified_details,
                    # Record actual limits used (Task C)
                    "limits_used": {
                        "max_delta_per_atom": max_delta_limit,
                        "max_total": max_total_limit,
                        "max_dipole_shift_debye": dipole_limit,
                    },
                    "strict_mode_enforced": bool(ctx.strict_charge_neutrality),
                    "classification_audit": ion_classification_audit,
                }
                electrolyte_warning = self._warn_electrolyte_fixed_charge_limitations(
                    correction_result=correction_result,
                    classifications=ion_classifications,
                    ctx=ctx,
                )
                if electrolyte_warning:
                    charge_summary["electrolyte_limitations"] = electrolyte_warning
                
                if correction_result.corrected:
                    # Task C: Use coalesced limits for comparison
                    if correction_result.max_delta_applied is not None and (
                        correction_result.max_delta_applied > max_delta_limit
                    ):
                        raise SanitizerError(
                            f"Charge correction exceeds per-atom limit: |Δq|={correction_result.max_delta_applied:.3e} "
                            f"> {max_delta_limit:.3e} (limit used)"
                        )
                    # v4 Issue 5: Use abs(total_delta) and guard for None
                    if (correction_result.total_delta is not None 
                        and max_total_limit is not None
                        and abs(correction_result.total_delta) > max_total_limit
                    ):
                        raise SanitizerError(
                            f"Charge correction exceeds total limit: |ΔQ|={abs(correction_result.total_delta):.3e} "
                            f"> {max_total_limit:.3e} (limit used)"
                        )
                    
                    # v4 Issue 5: Warn for uniform_all correction method
                    if (
                        getattr(correction_result, "effective_method", None) == "uniform_all"
                        or str(correction_result.correction_method or "").startswith("uniform_all:")
                    ):
                        print(
                            f"  [WARN] uniform_all correction applied (total_delta={correction_result.total_delta:.3e}). "
                            "Verify protected patterns/atomtypes were not affected."
                        )
                    
                    dipole_info = self._compute_dipole_shift(
                        system_gromacs_dir=system_gromacs_dir,
                        ordered_molecules=ordered_molecules,
                        molecule_itp_paths=sanitized_molecule_itp_paths,
                        correction_result=correction_result,
                        ctx=ctx,  # Hardening v5: explicit ctx passing
                    )
                    if dipole_info is not None:
                        charge_summary["dipole"] = {
                            "computed": dipole_info.computed,
                            "dipole_before_debye": dipole_info.dipole_before_debye,
                            "dipole_after_debye": dipole_info.dipole_after_debye,
                            "delta_dipole_debye": dipole_info.delta_dipole_debye,
                            "avg_delta_dipole_debye": dipole_info.delta_dipole_debye,
                            "max_delta_dipole_debye": dipole_info.max_delta_dipole_debye,
                            "reliable": bool(getattr(dipole_info, "reliable", False)),
                            "reliability_reasons": list(getattr(dipole_info, "reliability_reasons", [])),
                            "coordinates_file": dipole_info.coordinates_file,
                            "skip_reasons": list(getattr(dipole_info, "skip_reasons", [])),
                        }
                        if dipole_info.computed and dipole_info.delta_dipole_debye is not None:
                            max_delta_val = (
                                dipole_info.max_delta_dipole_debye
                                if dipole_info.max_delta_dipole_debye is not None
                                else dipole_info.delta_dipole_debye
                            )
                            print(
                                "  Dipole shift summary: "
                                f"avg |Δμ|={dipole_info.delta_dipole_debye:.4f} D, "
                                f"max |Δμ|={max_delta_val:.4f} D"
                            )
                            # v4: Dipole check is heuristic by default
                            # Only hard-fail if strict_dipole_check=True AND no gating reasons
                            if dipole_info.delta_dipole_debye > dipole_limit:
                                msg = (
                                    f"Dipole shift |Δμ|={dipole_info.delta_dipole_debye:.4f} D exceeds "
                                    f"limit {dipole_limit:.4f} D"
                                )
                                # Check gating: skip enforcement for giant/percolated molecules
                                skip_reasons = getattr(dipole_info, 'skip_reasons', [])
                                can_enforce = bool(getattr(dipole_info, "reliable", False)) and not skip_reasons
                                if ctx.strict_dipole_check and can_enforce:
                                    raise SanitizerError(msg)
                                else:
                                    reason_str = f"; skip reasons: {skip_reasons}" if skip_reasons else " (heuristic mode)"
                                    print(f"  [WARN] {msg}{reason_str}")
                    
                    backup_itp = self._sync_corrected_itp(
                        correction_result=correction_result,
                        molecule_itp_paths=sanitized_molecule_itp_paths,
                        sanitized_dir=result.sanitized_dir,
                        sanitized_current_dir=result.sanitized_current_dir,
                    )
                    if backup_itp is not None:
                        charge_summary["pre_patch_backup_itp"] = str(backup_itp)
                else:
                    # v4: Use ctx.charge_neutrality_tol instead of hard-coded 1e-10
                    allowed_skip_statuses = {"polymer_acceptable_drift", "polymer_skipped_non_strict"}
                    if (
                        abs(correction_result.original_q) > ctx.charge_neutrality_tol
                        and correction_result.neutrality_status not in allowed_skip_statuses
                    ):
                        raise SanitizerError(
                            f"Charge neutrality correction failed: {correction_result.correction_method}"
                        )
                
                # Hardening v5 - Issue G: Post-correction ion protection check
                protection_audit = self._verify_ion_protection(
                    correction_result,
                    ion_classifications,
                    ctx,
                    classification_audit=ion_classification_audit,
                )
                violations = protection_audit.get("violations", [])
                charge_summary["target_policy"] = protection_audit
                if violations:
                    for v in violations:
                        print(f"  [ERROR] {v}")
                    raise SanitizerError(
                        f"Charge protection policy violations: {violations}. "
                        "Use --charge-fix-allow-ions/--charge-fix-allow-solvents, and for polymer "
                        "targets provide --charge-fix-target-allowlist with "
                        "--charge-fix-polymer-method spread_safe."
                    )
                
                # Record in manifest
                if ctx.manifest:
                    ctx.manifest.set_sanitizer_output("charge_neutrality", charge_summary)
                    
            except ChargeNeutralityError as e:
                print(f"  [FAIL] {e}")
                raise SanitizerError(str(e))
        else:
            print("  [SKIP] Charge neutrality check (no active molecules in [molecules])")
        
        # === Validate GRO/TOP consistency ===
        # Fail fast if coordinate file doesn't match topology molecule counts
        gro_dir = system_gromacs_dir / "gro"
        gro_path = gro_dir / "system.gro"
        if gro_path.exists() and molecule_counts and sanitized_molecule_itp_paths:
            if not self._validate_gro_top_consistency(
                gro_path,
                system_top_path=top_path if top_path.exists() else None,
                molecule_counts=molecule_counts,
                molecule_itp_paths=sanitized_molecule_itp_paths,
                ctx=ctx,
                output_dir=self.get_output_dir(ctx),
            ):
                raise SanitizerError(
                    f"GRO/TOP atom count mismatch. The coordinate file does not match "
                    f"the topology. This can happen if HTPolyNet crosslinking changed "
                    f"the molecule composition but the GRO was not updated."
                )
        
        # === Update system.top ===
        
        top_dir = system_gromacs_dir / "top"
        top_dir.mkdir(parents=True, exist_ok=True)
        
        system_top_path = top_dir / "system.top"
        molecule_itp_include_paths = self._ordered_molecule_itp_paths(
            ordered_molecules,
            sanitized_molecule_itp_paths,
            result.sanitized_current_dir,
        )
        extra_includes: List[str] = []
        prefix_injected_current = getattr(result, "prefix_injected_types_current_path", None)
        if prefix_injected_current and Path(prefix_injected_current).exists():
            extra_includes.append(
                Path(os.path.relpath(prefix_injected_current, top_dir)).as_posix()
            )
        cross_group_nonbond_current = getattr(result, "nonbond_params_cross_group_current_path", None)
        if cross_group_nonbond_current and Path(cross_group_nonbond_current).exists():
            extra_includes.append(
                Path(os.path.relpath(cross_group_nonbond_current, top_dir)).as_posix()
            )
        secondary_nonbond_current = getattr(result, "nonbond_params_secondary_current_path", None)
        if secondary_nonbond_current and Path(secondary_nonbond_current).exists():
            extra_includes.append(
                Path(os.path.relpath(secondary_nonbond_current, top_dir)).as_posix()
            )
        
        existing_top_before = system_top_path.exists()
        molecule_source_used = molecule_count_source.get("source")
        preserve_existing_molecules = (
            existing_top_before and molecule_source_used == "raw_system_top"
        )
        update_strategy = "generated_from_trusted_counts"
        strategy_details: Dict[str, Any] = {}
        strategy_degraded = False
        if preserve_existing_molecules:
            print(f"  Updating existing system.top (preserving [molecules] section)...")
            top_content = self._update_system_top_includes(
                existing_top=system_top_path,
                combined_atomtypes_path=result.combined_atomtypes_current_path,
                molecule_itp_files=molecule_itp_include_paths,
                defaults=result.defaults_used,
                top_dir=top_dir,
                allow_default_defaults=ctx.allow_default_defaults,
                ctx=ctx,
                extra_includes=extra_includes,
            )
            strategy_details = dict(getattr(self, "_system_top_update_status", {}))
            update_strategy = str(strategy_details.get("mode", "updated_existing_block"))
            strategy_degraded = bool(strategy_details.get("degraded", False))
        else:
            if existing_top_before and molecule_source_used != "raw_system_top":
                print(
                    "  Rebuilding existing system.top from trusted molecule-count source "
                    f"'{molecule_source_used}' to avoid stale [molecules] state."
                )
            # Generate new system.top with proper include order
            top_content = generate_system_top(
                combined_atomtypes_path=result.combined_atomtypes_current_path,
                sanitized_itp_paths=result.sanitized_files,
                system_name=ctx.system_id,
                molecule_counts=molecule_counts,
                defaults=result.defaults_used,
                molecule_itp_files=molecule_itp_include_paths,
                top_dir=top_dir,
                allow_guess_defaults=ctx.allow_default_defaults,
                extra_includes=extra_includes,
            )

        defaults_sections_in_top = sum(
            1 for line in top_content.splitlines() if parse_section_name(line) == "defaults"
        )
        if defaults_sections_in_top > 1:
            raise SanitizerError(
                "Generated system.top contains multiple [ defaults ] sections. "
                "Refusing to proceed because GROMACS accepts only one global defaults tuple."
            )
        
        previous_top = ""
        if system_top_path.exists():
            previous_top = self.safe_read_text(system_top_path, strict=True)
        top_changed = previous_top != top_content
        if top_changed:
            self._atomic_write_text(system_top_path, top_content, encoding="utf-8")
            print(f"  Wrote: {system_top_path}")
        else:
            print(f"  system.top unchanged: {system_top_path}")
        self._validate_system_top_includes(system_top_path, ctx)
        self._system_top_update_status = {
            "mode": "updated_existing" if existing_top_before else "generated_new",
            "strategy": update_strategy,
            "degraded": bool(strategy_degraded),
            "changed": bool(top_changed),
            "molecule_count_source": molecule_source_used,
            "preserved_existing_molecules": bool(preserve_existing_molecules),
            "strategy_details": strategy_details,
        }
        
        # === Record in manifest ===
        
        if ctx.manifest:
            include_priority_raw = getattr(
                ctx,
                "itp_include_dirs_priority",
                DEFAULT_INCLUDE_PRIORITY,
            )
            include_priority_norm = normalize_include_priority(include_priority_raw)
            ctx.manifest.set_sanitizer_output(
                "include_resolution_policy",
                {
                    "itp_include_dirs_priority": include_priority_raw,
                    "normalized_priority": include_priority_norm,
                    "allow_include_shadowing": bool(getattr(ctx, "allow_include_shadowing", False)),
                    "allow_unsafe_include_escape": bool(
                        getattr(ctx, "allow_unsafe_include_escape", DEFAULT_ALLOW_UNSAFE_INCLUDE_ESCAPE)
                    ),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "system_top_update",
                {
                    "path": str(system_top_path),
                    "status": getattr(self, "_system_top_update_status", {}),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "molecule_count_source",
                getattr(self, "_molecule_count_source", {}),
            )
            ctx.manifest.set_sanitizer_output(
                "combined_atomtypes",
                {
                    "versioned": str(result.combined_atomtypes_path),
                    "current": str(result.combined_atomtypes_current_path),
                    "atomtype_count": result.atomtype_count,
                }
            )
            ctx.manifest.set_sanitizer_output(
                "sanitized_itps",
                {
                    "versioned_dir": str(result.sanitized_dir),
                    "current_dir": str(result.sanitized_current_dir),
                    "file_count": len(result.sanitized_files),
                    "molecule_itp_count": len(result.molecule_itp_files),
                    "files": [str(f) for f in result.sanitized_files],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "conflicts",
                {
                    "detected": len(result.conflicts_detected),
                    "overridden": len(result.conflicts_overridden),
                    "allow_override_flag": allow_override,
                    "atomtype_prefix": atomtype_prefix,
                    "overridden_details": [
                        {
                            "atomtype": c.atomtype_name,
                            "winner": Path(c.winner.source_files[0]).name if c.winner and c.winner.source_files else None,
                        }
                        for c in result.conflicts_overridden
                    ],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "defaults",
                {
                    "used": result.defaults_used is not None,
                    "nbfunc": result.defaults_used.nbfunc if result.defaults_used else None,
                    "comb_rule": result.defaults_used.comb_rule if result.defaults_used else None,
                    "source_file": Path(result.defaults_used.source_file).name if result.defaults_used else None,
                    "guessed": getattr(result, "defaults_guessed", False),
                    "allow_mixed_defaults": bool(getattr(ctx, "allow_mixed_defaults", False)),
                    "allow_mixed_defaults_reason": getattr(ctx, "allow_mixed_defaults_reason", None),
                    "mixed_defaults_cross_group_policy": getattr(
                        ctx, "mixed_defaults_cross_group_policy", None
                    ),
                    "mixed_defaults_cross_group_rule": getattr(
                        ctx, "mixed_defaults_cross_group_rule", None
                    ),
                    "mixed_defaults_cross_group_max_pairs": getattr(
                        ctx, "mixed_defaults_cross_group_max_pairs", None
                    ),
                    "mixed_defaults_cross_group_reason": getattr(
                        ctx, "mixed_defaults_cross_group_reason", None
                    ),
                    "mixed_defaults_preserve_within_group": bool(
                        getattr(ctx, "mixed_defaults_preserve_within_group", False)
                    ),
                    "unsafe_mixed_defaults_override_used": bool(
                        getattr(ctx, "allow_mixed_defaults", False)
                        and getattr(result, "mixed_defaults_detected", False)
                    ),
                    "primary_defaults_signature": getattr(result, "mixed_defaults_report", {}).get(
                        "primary_defaults_signature"
                    ),
                    "secondary_signatures": getattr(result, "mixed_defaults_report", {}).get(
                        "secondary_signatures", []
                    ),
                    "chosen_policy": getattr(result, "mixed_defaults_report", {}).get("chosen_policy"),
                    "classification": getattr(result, "mixed_defaults_report", {}).get("classification"),
                    "differing_fields": getattr(result, "mixed_defaults_report", {}).get(
                        "differing_fields", []
                    ),
                    "comb_rules_present": getattr(result, "mixed_defaults_report", {}).get(
                        "comb_rules_present", []
                    ),
                    "grouped_report": getattr(result, "mixed_defaults_report", {}).get("groups", []),
                    "source_defaults": [
                        {
                            "source_file": str(path),
                            "nbfunc": entry.nbfunc,
                            "comb_rule": entry.comb_rule,
                            "gen_pairs": entry.gen_pairs,
                            "fudge_lj": entry.fudge_lj,
                            "fudge_qq": entry.fudge_qq,
                        }
                        for path, entry in sorted(
                            getattr(result, "source_defaults_map", {}).items(),
                            key=lambda item: item[0].as_posix().casefold(),
                        )
                    ],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "mixed_defaults_within_group_preservation",
                {
                    "versioned": str(result.nonbond_params_secondary_path)
                    if getattr(result, "nonbond_params_secondary_path", None)
                    else None,
                    "current": str(result.nonbond_params_secondary_current_path)
                    if getattr(result, "nonbond_params_secondary_current_path", None)
                    else None,
                    "summary": getattr(result, "nonbond_params_secondary_summary", {}),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "mixed_defaults_cross_group_override",
                {
                    "versioned": str(result.nonbond_params_cross_group_path)
                    if getattr(result, "nonbond_params_cross_group_path", None)
                    else None,
                    "current": str(result.nonbond_params_cross_group_current_path)
                    if getattr(result, "nonbond_params_cross_group_current_path", None)
                    else None,
                    "summary": getattr(result, "nonbond_params_cross_group_summary", {}),
                },
            )
            ctx.manifest.set_sanitizer_output(
                "prefix_implicit_topology",
                {
                    "policy_summary": getattr(result, "prefix_policy_summary", {}),
                    "injected_types_versioned": str(result.prefix_injected_types_path)
                    if getattr(result, "prefix_injected_types_path", None)
                    else None,
                    "injected_types_current": str(result.prefix_injected_types_current_path)
                    if getattr(result, "prefix_injected_types_current_path", None)
                    else None,
                },
            )
            if bool(getattr(ctx, "allow_mixed_defaults", False)) and bool(
                getattr(result, "mixed_defaults_detected", False)
            ):
                unsafe_overrides = ctx.manifest.get("unsafe_overrides", [])
                if not isinstance(unsafe_overrides, list):
                    unsafe_overrides = []
                unsafe_overrides.append(
                    {
                        "kind": "mixed_defaults_override",
                        "unsafe_override": True,
                        "reason": getattr(ctx, "allow_mixed_defaults_reason", None),
                        "primary_defaults_signature": getattr(
                            result, "mixed_defaults_report", {}
                        ).get("primary_defaults_signature"),
                        "secondary_signatures": getattr(
                            result, "mixed_defaults_report", {}
                        ).get("secondary_signatures", []),
                        "chosen_policy": getattr(result, "mixed_defaults_report", {}).get(
                            "chosen_policy"
                        ),
                    }
                )
                cross_group_summary = getattr(result, "nonbond_params_cross_group_summary", {})
                if (
                    isinstance(cross_group_summary, dict)
                    and cross_group_summary.get("policy") == "generate"
                ):
                    unsafe_overrides.append(
                        {
                            "kind": "mixed_defaults_cross_group_generate",
                            "unsafe_override": True,
                            "reason": cross_group_summary.get("reason"),
                            "rule": cross_group_summary.get("rule"),
                            "pair_count": cross_group_summary.get("pair_count"),
                            "file": cross_group_summary.get("current_path"),
                        }
                    )
                ctx.manifest.set("unsafe_overrides", unsafe_overrides)
            ctx.manifest.set_sanitizer_output(
                "charge_warnings",
                {
                    "count": len(result.charge_warnings),
                    "strict_mode": strict_charge,
                    "warnings": [str(w) for w in result.charge_warnings],
                }
            )
            ctx.manifest.set_sanitizer_output(
                "lj_outlier_warnings",
                {
                    "count": len(outlier_warnings),
                    "warnings": outlier_warnings,
                },
            )
            ctx.manifest.set_sanitizer_output(
                "mixed_defaults_warnings",
                {
                    "count": len(mixed_defaults_warnings),
                    "warnings": mixed_defaults_warnings,
                },
            )
            
            # v5 Hardening: Include resolution report (Issue 1)
            ctx.manifest.set_sanitizer_output(
                "include_resolution",
                [self._serialize_manifest_entry(d) for d in self.diagnostics]
            )
            ctx.manifest.save()
        
        # Log output for stage directory
        output_dir = self.get_output_dir(ctx)
        log_path = output_dir / "sanitizer.log"
        
        log_lines = [
            "ITP Sanitizer completed",
            f"Atomtypes: {result.atomtype_count}",
            f"Files: {len(result.sanitized_files)}",
            f"Molecule ITPs: {len(result.molecule_itp_files)}",
            f"Conflicts: {len(result.conflicts_detected)}",
            f"Overridden: {len(result.conflicts_overridden)}",
            f"Charge warnings: {len(result.charge_warnings)}",
            f"LJ outlier warnings: {len(outlier_warnings)}",
        ]
        
        if result.defaults_used:
            log_lines.append(f"Defaults: nbfunc={result.defaults_used.nbfunc}, comb-rule={result.defaults_used.comb_rule}")
        defaults_report = getattr(result, "mixed_defaults_report", {})
        if defaults_report:
            log_lines.append(
                "Defaults policy: "
                f"classification={defaults_report.get('classification')}, "
                f"chosen_policy={defaults_report.get('chosen_policy')}, "
                f"primary={defaults_report.get('primary_defaults_signature')}"
            )
        nonbond_summary = getattr(result, "nonbond_params_secondary_summary", {})
        if nonbond_summary:
            log_lines.append(
                "Secondary nonbond_params preserve: "
                f"requested={nonbond_summary.get('requested')}, "
                f"applied={nonbond_summary.get('applied')}, "
                f"pairs={nonbond_summary.get('pair_count')}, "
                f"reason={nonbond_summary.get('reason')}"
            )
        cross_group_summary = getattr(result, "nonbond_params_cross_group_summary", {})
        if cross_group_summary:
            log_lines.append(
                "Cross-group nonbond_params: "
                f"policy={cross_group_summary.get('policy')}, "
                f"rule={cross_group_summary.get('rule')}, "
                f"applied={cross_group_summary.get('applied')}, "
                f"pairs={cross_group_summary.get('pair_count')}, "
                f"status={cross_group_summary.get('status')}"
            )
        prefix_summary = getattr(result, "prefix_policy_summary", {})
        if prefix_summary:
            log_lines.append(
                "Prefix implicit-topology policy: "
                f"policy={prefix_summary.get('policy')}, "
                f"strategy={prefix_summary.get('applied_strategy')}, "
                f"effective_prefix={prefix_summary.get('effective_prefix')}"
            )
        if mixed_defaults_warnings:
            log_lines.append("Mixed defaults safety warnings:")
            for warning in mixed_defaults_warnings:
                log_lines.append(f"  - {warning}")
        if outlier_warnings:
            log_lines.append("LJ outlier details:")
            for warning in outlier_warnings:
                log_lines.append(f"  - {warning}")
        
        if result.conflicts_overridden:
            log_lines.append("\nOverridden conflicts (winner selected by priority):")
            for c in result.conflicts_overridden:
                winner_file = Path(c.winner.source_files[0]).name if c.winner and c.winner.source_files else "?"
                log_lines.append(f"  - {c.atomtype_name}: {winner_file}")
        
        self._atomic_write_text(log_path, "\n".join(log_lines) + "\n", encoding="utf-8")
        
        print("  [OK] ITP Sanitizer stage complete")
        return True
    

