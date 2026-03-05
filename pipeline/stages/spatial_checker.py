"""
Spatial consistency helpers extracted from sanitizer.py.

Contains GRO/PBC/unwrap/dipole/grompp-based consistency checks.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple
import math
import os
import re
import shlex
import subprocess

from .topology_sanitizer import (
    SanitizerError,
    DEFAULT_GRO_COORD_RESOLUTION_NM,
    GRO_ABSURD_ATOMS,
    FLOAT_FINDER_RE,
    GPE_UNWRAP_LOOP_TOL_NM_ENV,
    UNWRAP_LOOP_TOL_FACTOR,
    DEFAULT_GROMPP_PP_TIMEOUT_S,
    GROMPP_TIMEOUT_FLOOR_S,
    GROMPP_TIMEOUT_PER_1K_LINES,
    GROMPP_TIMEOUT_PER_10_INCLUDES,
    GROMPP_TIMEOUT_RETRY_MULTIPLIER,
    parse_section_name,
    strip_gmx_comment,
    _extract_floats_from_text,
    _decimal_places_from_numeric_token,
)

if TYPE_CHECKING:
    from ..context import PipelineContext


@dataclass
class GROAtom:
    """Atom data from a GRO file line."""
    resnr: int
    resname: str
    atomname: str
    x: float
    y: float
    z: float


@dataclass
class GROParseResult:
    """Result of parsing a GRO file with box vectors."""
    atoms: List[GROAtom]
    box_x: float
    box_y: float
    box_z: float
    is_triclinic: bool
    coord_resolution_nm: float = DEFAULT_GRO_COORD_RESOLUTION_NM
    title: str = ""


class SpatialCheckerMixin:
    """GRO/PBC/dipole/grompp consistency mixin."""
    def _find_coords_gro(self, system_gromacs_dir: Path) -> Optional[Path]:
        """Locate a coordinate file for dipole estimation."""
        coords_path = system_gromacs_dir / "coords" / "system.gro"
        if coords_path.exists():
            return coords_path
        gro_path = system_gromacs_dir / "gro" / "system.gro"
        if gro_path.exists():
            return gro_path
        return None

    def _read_gro_coords(
        self,
        gro_path: Path,
        retain_resnames: Optional[Set[str]] = None,
        allow_full_parse: bool = False,
    ) -> Optional[List[Tuple[float, float, float]]]:
        """Read GRO coordinates with safe-by-default selective retention."""
        if retain_resnames is None and not allow_full_parse:
            # Compatibility mode: no residue filter means "do not materialize atoms".
            # Keep this path O(1) memory by only validating the GRO header count line.
            return [] if self._count_atoms_in_gro(gro_path) is not None else None

        retain_filter: Optional[Set[str]]
        if retain_resnames is None:
            declared_atoms = self._count_atoms_in_gro(gro_path)
            if declared_atoms is None:
                return None
            if declared_atoms > GRO_ABSURD_ATOMS:
                raise SanitizerError(
                    f"_read_gro_coords() refused full GRO parse for '{gro_path}' "
                    f"(declared atoms={declared_atoms}, cap={GRO_ABSURD_ATOMS}). "
                    "Use retain_resnames for selective reads, or call "
                    "_parse_gro_file(..., retain_resnames=...) with a targeted filter."
                )
            retain_filter = None
        else:
            retain_filter = {name.upper() for name in retain_resnames}
        result = self._parse_gro_file(gro_path, retain_resnames=retain_filter)
        if result is None:
            return None
        return [(a.x, a.y, a.z) for a in result.atoms]

    def _parse_gro_file(
        self,
        gro_path: Path,
        strict: bool = False,
        max_atoms: Optional[int] = None,
        retain_resnames: Optional[Set[str]] = None,
    ) -> Optional[GROParseResult]:
        """
        Parse a GRO file returning atoms with residue info and box vectors.
        
        Hardening v6: Robust parsing with truncation/corruption detection.
        Hardening v7: Streamed parsing (avoid loading huge files into memory).
        Hardening v8: Optional selective retention by residue name for dipole checks.
        
        GRO format (fixed columns):
        - Line 1: title
        - Line 2: natoms (integer count)
        - Lines 3 to 2+natoms: atom lines
        - Last line: box vectors (3 or 9 floats)
        
        Atom line format (fixed columns):
        - Columns 1-5: residue number
        - Columns 6-10: residue name  
        - Columns 11-15: atom name
        - Columns 16-20: atom number
        - Columns 21-28: x coordinate (nm)
        - Columns 29-36: y coordinate (nm)
        - Columns 37-44: z coordinate (nm)
        
        Args:
            gro_path: Path to GRO file
            strict: If True, fail on any parsing issue; else return None
            max_atoms: Optional safety cap for declared atoms (non-strict early exit)
            retain_resnames: Optional uppercase residue names to retain in memory.
                             When set, all lines are still parsed/validated but only
                             matching atoms are materialized.
            
        Returns:
            GROParseResult with atoms, box dimensions, and triclinic flag
            
        Raises:
            SanitizerError: In strict mode on truncation or format errors
        """
        try:
            with open(gro_path, "r") as handle:
                title_line = handle.readline()
                count_line = handle.readline()
                if not title_line or not count_line:
                    if strict:
                        raise SanitizerError(
                            f"GRO file {gro_path} is too short (missing title/count line)"
                        )
                    return None
                title = title_line.strip()
                try:
                    natoms_declared = int(count_line.strip())
                except ValueError:
                    if strict:
                        raise SanitizerError(
                            f"GRO file {gro_path}: line 2 is not a valid atom count: {count_line!r}"
                        )
                    return None

                if natoms_declared < 0:
                    if strict:
                        raise SanitizerError(
                            f"GRO file {gro_path}: negative atom count declared: {natoms_declared}"
                        )
                    return None

                # Safety cap for absurdly large files (non-strict only)
                cap = max_atoms if max_atoms is not None else GRO_ABSURD_ATOMS
                if cap is not None and natoms_declared > cap:
                    msg = (
                        f"GRO file {gro_path}: declared atom count {natoms_declared} exceeds safety cap {cap}. "
                        "This parser keeps only selected residues for dipole checks; "
                        "reduce system size for sanitizer dipole checks or raise the cap intentionally."
                    )
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                    return None

                # Parse atom lines (line numbers start at 3 for first atom)
                atoms: List[GROAtom] = []
                parse_errors: List[str] = []
                atoms_parsed_total = 0
                max_coord_decimals = 0
                for i in range(natoms_declared):
                    line = handle.readline()
                    if not line:
                        parse_errors.append(f"Line {i + 3}: EOF before atom data")
                        break
                    if len(line) < 20:
                        parse_errors.append(
                            f"Line {i + 3}: too short for atom ({len(line)} chars)"
                        )
                        continue
                    try:
                        resnr = int(line[0:5].strip())
                        resname = line[5:10].strip()
                        atomname = line[10:15].strip()
                        x_token = line[20:28].strip()
                        y_token = line[28:36].strip()
                        z_token = line[36:44].strip()
                        x = float(x_token)
                        y = float(y_token)
                        z = float(z_token)
                        for token in (x_token, y_token, z_token):
                            decimals = _decimal_places_from_numeric_token(token)
                            if decimals is not None:
                                max_coord_decimals = max(max_coord_decimals, decimals)
                        atoms_parsed_total += 1
                        if retain_resnames is None or resname.upper() in retain_resnames:
                            atoms.append(GROAtom(resnr, resname, atomname, x, y, z))
                    except (ValueError, IndexError) as e:
                        parse_errors.append(f"Line {i + 3}: {e}")
                        continue

                # Hardening v6: Verify natoms matches parsed atoms
                if atoms_parsed_total != natoms_declared:
                    msg = (
                        f"GRO file {gro_path}: atom count mismatch!\n"
                        f"  Declared: {natoms_declared}, Parsed: {atoms_parsed_total}\n"
                        f"  File may be truncated or corrupted."
                    )
                    if parse_errors:
                        msg += f"\n  Parse errors: {parse_errors[:5]}"
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                    if atoms_parsed_total == 0:
                        return None

                # Find the box line robustly by scanning remaining lines
                box_floats: List[float] = []
                last_nonempty: Optional[str] = None
                for line in handle:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    last_nonempty = stripped
                    parts = stripped.split()
                    try:
                        floats = [float(x) for x in parts]
                        if len(floats) == 3 or len(floats) == 9:
                            box_floats = floats
                            break
                    except ValueError:
                        continue

                if not box_floats and last_nonempty:
                    parts = last_nonempty.split()
                    try:
                        floats = [float(x) for x in parts]
                        if len(floats) == 3 or len(floats) == 9:
                            box_floats = floats
                    except ValueError:
                        pass

                if not box_floats:
                    msg = f"GRO file {gro_path}: no valid box vector line found"
                    if strict:
                        raise SanitizerError(msg)
                    print(f"  [WARN] {msg}")
                    box_floats = [0.0, 0.0, 0.0]
        except Exception as e:
            if strict:
                raise SanitizerError(f"Cannot read GRO file {gro_path}: {e}")
            return None
        
        is_triclinic = len(box_floats) >= 9
        box_x = box_floats[0] if len(box_floats) >= 1 else 0.0
        box_y = box_floats[1] if len(box_floats) >= 2 else 0.0
        box_z = box_floats[2] if len(box_floats) >= 3 else 0.0
        coord_resolution_nm = (
            10.0 ** (-max_coord_decimals)
            if max_coord_decimals > 0
            else DEFAULT_GRO_COORD_RESOLUTION_NM
        )
        
        return GROParseResult(
            atoms=atoms,
            box_x=box_x,
            box_y=box_y,
            box_z=box_z,
            is_triclinic=is_triclinic,
            coord_resolution_nm=coord_resolution_nm,
            title=title,
        )

    def _get_molecule_charges(self, itp_path: Path, mol_name: str) -> Optional[List[float]]:
        """Extract charges list for a specific moleculetype."""
        from ..charge_neutrality import ChargeParser
        per_mol = ChargeParser.parse_atoms_charges_per_molecule(itp_path)
        if mol_name in per_mol:
            return per_mol[mol_name][0]
        key = mol_name.casefold()
        matches = [name for name in per_mol if name.casefold() == key]
        if len(matches) == 1:
            return per_mol[matches[0]][0]
        return None

    def _compute_dipole_shift(
        self,
        system_gromacs_dir: Path,
        ordered_molecules: List[Tuple[str, int]],
        molecule_itp_paths: Dict[str, Path],
        correction_result,
        ctx: Optional["PipelineContext"] = None,
    ):
        """
        Compute approximate dipole shift for corrected molecule(s).
        
        PBC-aware implementation using residue grouping and topology-aware
        bond-graph unwrapping when bonds are available.
        
        v4: Adds gating logic for heuristic mode:
        - Giant molecule detection (>dipole_unwrap_max_atoms)
        - Net charge guard (skip enforcement if molecule is charged)
        - Percolation check via bond-loop image inconsistency (preferred)
          with span-based fallback only when bond graph is unavailable
        - Triclinic fallback to diagonal approximation
        """
        from ..charge_neutrality import ChargeParser, DipoleResult

        skip_reasons: List[str] = []
        reliability_reasons: List[str] = []

        def _add_reason(target: List[str], reason: str) -> None:
            if reason not in target:
                target.append(reason)
        
        if not correction_result.patched_itp_path:
            return None
        
        target = correction_result.target_molecule
        if not target and correction_result.patched_itp_path:
            name = correction_result.patched_itp_path.stem
            if name.endswith("_corrected"):
                target = name[: -len("_corrected")]
        if not target:
            return DipoleResult(computed=False)
        target_casefold_map = {name.casefold(): name for name in molecule_itp_paths}
        target_resolved = target_casefold_map.get(target.casefold())
        if not target_resolved:
            return DipoleResult(computed=False)
        target = target_resolved
        target_itp_name = (
            getattr(correction_result, "target_itp_moleculetype_name", None) or target
        )
        # Handle GRO 5-char truncation: try target and target[:5]
        target_resnames = {target.upper(), target[:5].upper()}

        charges_before = self._get_molecule_charges(molecule_itp_paths[target], target_itp_name)
        charges_after = self._get_molecule_charges(correction_result.patched_itp_path, target_itp_name)
        if not charges_before or not charges_after or len(charges_before) != len(charges_after):
            return DipoleResult(computed=False, reliable=False)
        bond_source = correction_result.patched_itp_path or molecule_itp_paths[target]
        bond_graph = ChargeParser.parse_molecule_bond_graph(bond_source, target_itp_name)
        has_bond_graph = bool(bond_graph)
        if not has_bond_graph:
            _add_reason(reliability_reasons, "bond_graph_missing")
            _add_reason(reliability_reasons, "dipole_unwrap_requires_bond_graph")
            _add_reason(skip_reasons, "bond_graph_missing")
            msg = (
                f"Dipole reliability failure for '{target_itp_name}': "
                "bond graph is unavailable; skipping dipole check."
            )
            if ctx and ctx.strict_dipole_check:
                raise SanitizerError(f"{msg} Disable strict dipole checks or provide bond topology.")
            print(f"  [WARN] {msg} (non-strict mode marks result unreliable)")
            gro_path = self._find_coords_gro(system_gromacs_dir)
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                skip_reasons=skip_reasons,
                coordinates_file=str(gro_path) if gro_path else None,
            )
        
        # v4: Giant molecule gating (skip enforcement for large molecules)
        expected_natoms = len(charges_before)
        # Hardening v5: Use explicit ctx parameter (no more getattr(self, '_ctx'))
        dipole_max_atoms = ctx.dipole_unwrap_max_atoms if ctx else 512
        if expected_natoms > dipole_max_atoms:
            reason = f"giant_molecule ({expected_natoms} > {dipole_max_atoms} atoms)"
            skip_reasons.append(reason)
            _add_reason(reliability_reasons, reason)
            print(f"  [INFO] Dipole enforcement gated: {reason}")
            gro_path = self._find_coords_gro(system_gromacs_dir)
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                skip_reasons=skip_reasons,
                coordinates_file=str(gro_path) if gro_path else None,
            )

        gro_path = self._find_coords_gro(system_gromacs_dir)
        if not gro_path:
            _add_reason(reliability_reasons, "coordinates_missing")
            return DipoleResult(computed=False, reliable=False, reliability_reasons=reliability_reasons)

        # Use selective parser: retain only target residue atoms for dipole checks.
        gro_data = self._parse_gro_file(gro_path, retain_resnames=target_resnames)
        if gro_data is None or not gro_data.atoms:
            _add_reason(reliability_reasons, "target_residue_coordinates_missing")
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                coordinates_file=str(gro_path),
            )

        # Skip if box dimensions are zero/invalid
        if gro_data.box_x <= 0 or gro_data.box_y <= 0 or gro_data.box_z <= 0:
            print("  [INFO] Dipole shift check skipped (reason: invalid box dimensions)")
            _add_reason(reliability_reasons, "invalid_box_dimensions")
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                coordinates_file=str(gro_path),
            )
        
        # v4: Net charge guard (dipole is origin-dependent for charged molecules)
        charge_tol = ctx.charge_neutrality_tol if ctx else 1e-6
        net_before = sum(charges_before)
        net_after = sum(charges_after)
        if abs(net_before) > charge_tol or abs(net_after) > charge_tol:
            reason = f"charged_molecule (net_q={net_before:.4e}/{net_after:.4e})"
            skip_reasons.append(reason)
            _add_reason(reliability_reasons, "charged_molecule_origin_dependent")
            print(f"  [INFO] Dipole enforcement gated: {reason}")
        
        # Hardening v6 - Issue C: Triclinic box policy
        # CRITICAL: Never apply orthorhombic MIC to triclinic boxes - it gives wrong results
        # Proper triclinic MIC requires box matrix inversion which is NOT implemented
        triclinic_policy_used: Optional[str] = None
        if gro_data.is_triclinic:
            # Determine effective policy
            policy = getattr(ctx, "dipole_triclinic_policy", "skip") if ctx else "skip"
            triclinic_policy_used = policy
            
            # If strict mode and no explicit policy, treat as "error"
            if ctx and ctx.strict_dipole_check and policy == "skip":
                # Check if user explicitly set the policy (vs default)
                # Default is "skip", so if strict we escalate to error unless user overrode
                if not hasattr(ctx, "_dipole_triclinic_policy_explicit"):
                    policy = "error"
                    triclinic_policy_used = "error (escalated from skip in strict mode)"
            
            if policy == "error":
                raise SanitizerError(
                    "Cannot compute accurate dipole shift for triclinic box.\n"
                    "Orthorhombic minimum-image convention (MIC) is mathematically invalid "
                    "for tilted boxes - it will give incorrect results.\n\n"
                    "Options:\n"
                    "  1. Re-box to orthorhombic: gmx editconf -bt cubic\n"
                    "  2. Skip check: --dipole-triclinic-policy skip\n"
                    "  3. Disable strict mode: --no-strict-dipole-check"
                )
            elif policy == "mic":
                # Hardening v6: REMOVED misleading "correct MIC" claim
                # Proper triclinic MIC requires matrix inversion - NOT IMPLEMENTED
                # Fail clearly rather than silently giving wrong results
                raise SanitizerError(
                    "Triclinic MIC policy='mic' is not supported.\n"
                    "Proper triclinic minimum-image convention requires box matrix inversion "
                    "which is not implemented in this version.\n\n"
                    "Options:\n"
                    "  1. Re-box to orthorhombic: gmx editconf -bt cubic\n"
                    "  2. Skip check: --dipole-triclinic-policy skip\n"
                    "  3. Use GROMACS 'gmx dipoles' for accurate triclinic dipole analysis"
                )
            else:  # policy == "skip"
                skip_reasons.append("triclinic_box (orthorhombic MIC invalid for tilted boxes)")
                _add_reason(reliability_reasons, "triclinic_box")
                print("  [INFO] Triclinic box detected; skipping dipole check (policy=skip)")
                return DipoleResult(
                    computed=False,
                    reliable=False,
                    reliability_reasons=reliability_reasons,
                    skip_reasons=skip_reasons,
                    coordinates_file=str(gro_path),
                )
        
        # Build residue groups: (resnr, resname) -> list of atoms
        residue_groups: Dict[Tuple[int, str], List[GROAtom]] = {}
        for atom in gro_data.atoms:
            key = (atom.resnr, atom.resname.upper())
            if key not in residue_groups:
                residue_groups[key] = []
            residue_groups[key].append(atom)
        
        # Find groups matching target
        target_groups: List[List[GROAtom]] = []
        for (resnr, resname), atoms in residue_groups.items():
            if resname in target_resnames:
                target_groups.append(atoms)
        
        if not target_groups:
            # No matching residues found
            _add_reason(reliability_reasons, "target_residue_not_found")
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                coordinates_file=str(gro_path),
            )
        
        group_max_atoms = (
            getattr(ctx, "dipole_group_max_atoms", dipole_max_atoms) if ctx else dipole_max_atoms
        )
        valid_groups: List[List[GROAtom]] = []
        for group in target_groups:
            group_size = len(group)
            if group_size > group_max_atoms:
                reason = f"group_too_large ({group_size} > {group_max_atoms} atoms)"
                if reason not in skip_reasons:
                    skip_reasons.append(reason)
                    _add_reason(reliability_reasons, reason)
                    print(f"  [INFO] Dipole enforcement gated: {reason}")
                continue
            if group_size == expected_natoms:
                valid_groups.append(group)
        
        if not valid_groups:
            # Atom count mismatch in all groups or all groups gated
            _add_reason(reliability_reasons, "no_valid_target_groups")
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                skip_reasons=skip_reasons,
                coordinates_file=str(gro_path),
            )
        
        # Compute dipole for each valid group with PBC unwrapping.
        # Prefer topology-aware bond graph unwrap; fallback to span-gated MIC unwrap.
        # v4: Issue 7 - Per-group vector delta: ||mu_after_vec - mu_before_vec||
        box = (gro_data.box_x, gro_data.box_y, gro_data.box_z)
        delta_magnitudes: List[float] = []  # Per-group vector delta norms
        mu_before_mags: List[float] = []  # For reporting
        mu_after_mags: List[float] = []
        percolation_threshold = (
            ctx.dipole_percolation_threshold if ctx else 0.45
        )
        loop_tol_nm = self._resolve_unwrap_loop_tol_nm(gro_data.coord_resolution_nm)
        
        for group in valid_groups:
            # Fast wrapped-span pre-gate: conservative early skip for obvious spanning groups.
            wrapped_span = self._wrapped_span(group)
            wrapped_span_ratio = max(
                wrapped_span[0] / box[0] if box[0] > 0 else 0,
                wrapped_span[1] / box[1] if box[1] > 0 else 0,
                wrapped_span[2] / box[2] if box[2] > 0 else 0,
            )
            if wrapped_span_ratio > percolation_threshold:
                print(
                    "  [INFO] Dipole wrapped-span suspicion: "
                    f"ratio={wrapped_span_ratio:.2f} > {percolation_threshold:.2f}. "
                    "Continuing with unwrap/percolation checks."
                )

            coords_list, topo_percolated, topo_msg = self._unwrap_molecule_bond_graph(
                group,
                box,
                bond_graph,
                tol_nm=loop_tol_nm,
            )
            if topo_percolated:
                reason = "percolated_topology"
                if reason not in skip_reasons:
                    skip_reasons.append(reason)
                    _add_reason(reliability_reasons, reason)
                    detail = f" ({topo_msg})" if topo_msg else ""
                    print(f"  [INFO] Dipole enforcement gated: {reason}{detail}")
                continue
            if not coords_list:
                continue
            
            # Compute vector dipoles (returns (mx, my, mz) in Debye)
            mu_before_vec = self._dipole_vector_unwrapped(coords_list, charges_before)
            mu_after_vec = self._dipole_vector_unwrapped(coords_list, charges_after)
            
            # v4 Issue 7: Per-group vector delta norm
            delta_vec = (
                mu_after_vec[0] - mu_before_vec[0],
                mu_after_vec[1] - mu_before_vec[1],
                mu_after_vec[2] - mu_before_vec[2],
            )
            delta_mag = math.sqrt(delta_vec[0]**2 + delta_vec[1]**2 + delta_vec[2]**2)
            delta_magnitudes.append(delta_mag)
            
            # Also store magnitudes for reporting
            mu_before_mags.append(math.sqrt(sum(x**2 for x in mu_before_vec)))
            mu_after_mags.append(math.sqrt(sum(x**2 for x in mu_after_vec)))
        
        if not delta_magnitudes:
            return DipoleResult(
                computed=False,
                reliable=False,
                reliability_reasons=reliability_reasons,
                skip_reasons=skip_reasons,
                coordinates_file=str(gro_path),
            )
        
        # Report average and max delta
        avg_before = sum(mu_before_mags) / len(mu_before_mags)
        avg_after = sum(mu_after_mags) / len(mu_after_mags)
        avg_delta = sum(delta_magnitudes) / len(delta_magnitudes)
        max_delta = max(delta_magnitudes)
        reliable = not skip_reasons and not reliability_reasons
        
        # Use average delta as the primary metric
        return DipoleResult(
            dipole_before_debye=avg_before,
            dipole_after_debye=avg_after,
            delta_dipole_debye=avg_delta,
            max_delta_dipole_debye=max_delta,
            coordinates_file=str(gro_path),
            computed=True,
            reliable=reliable,
            reliability_reasons=reliability_reasons,
            skip_reasons=skip_reasons,
        )

    def _unwrap_molecule_pbc(
        self, atoms: List[GROAtom], box: Tuple[float, float, float]
    ) -> List[GROAtom]:
        """
        Apply minimum-image unwrapping to a molecule's atoms (orthorhombic).
        
        Task A: Uses first atom as reference, unwraps others relative to it.
        This ensures a molecule straddling a PBC boundary is contiguous.
        """
        if not atoms:
            return atoms
        
        Lx, Ly, Lz = box
        ref = atoms[0]
        unwrapped = [ref]
        
        for atom in atoms[1:]:
            dx = atom.x - ref.x
            dy = atom.y - ref.y
            dz = atom.z - ref.z
            
            # Minimum image convention
            dx -= round(dx / Lx) * Lx
            dy -= round(dy / Ly) * Ly
            dz -= round(dz / Lz) * Lz
            
            unwrapped.append(GROAtom(
                resnr=atom.resnr,
                resname=atom.resname,
                atomname=atom.atomname,
                x=ref.x + dx,
                y=ref.y + dy,
                z=ref.z + dz,
            ))
        
        return unwrapped
    
    def _unwrap_molecule_pbc_with_span(
        self,
        atoms: List[GROAtom],
        box: Tuple[float, float, float],
        early_exit_ratio: Optional[float] = None,
    ) -> Tuple[Optional[List[Tuple[float, float, float]]], Tuple[float, float, float], bool]:
        """
        Apply minimum-image unwrapping with span tracking.

        Returns compact coordinates (x, y, z) tuples instead of GROAtom objects.
        Supports optional early exit when span ratio exceeds a threshold.
        """
        if not atoms:
            return [], (0.0, 0.0, 0.0), False
        
        Lx, Ly, Lz = box
        ref = atoms[0]
        
        # Initialize min/max tracking with reference atom
        min_x = max_x = ref.x
        min_y = max_y = ref.y
        min_z = max_z = ref.z
        
        coords: List[Tuple[float, float, float]] = [(ref.x, ref.y, ref.z)]
        
        def _span_ratio() -> float:
            span_x = max_x - min_x
            span_y = max_y - min_y
            span_z = max_z - min_z
            return max(
                span_x / Lx if Lx > 0 else 0.0,
                span_y / Ly if Ly > 0 else 0.0,
                span_z / Lz if Lz > 0 else 0.0,
            )
        
        for atom in atoms[1:]:
            dx = atom.x - ref.x
            dy = atom.y - ref.y
            dz = atom.z - ref.z
            
            # Minimum image convention
            dx -= round(dx / Lx) * Lx
            dy -= round(dy / Ly) * Ly
            dz -= round(dz / Lz) * Lz
            
            new_x = ref.x + dx
            new_y = ref.y + dy
            new_z = ref.z + dz
            
            # Stream min/max update
            min_x = min(min_x, new_x)
            max_x = max(max_x, new_x)
            min_y = min(min_y, new_y)
            max_y = max(max_y, new_y)
            min_z = min(min_z, new_z)
            max_z = max(max_z, new_z)
            
            coords.append((new_x, new_y, new_z))

            if early_exit_ratio is not None and _span_ratio() > early_exit_ratio:
                span = (max_x - min_x, max_y - min_y, max_z - min_z)
                return None, span, True
        
        span = (max_x - min_x, max_y - min_y, max_z - min_z)
        return coords, span, False

    def _resolve_unwrap_loop_tol_nm(self, base_resolution_nm: Optional[float]) -> float:
        """
        Resolve unwrap loop-consistency tolerance in nm.

        Precedence:
        1) Env override: GPE_UNWRAP_LOOP_TOL_NM
        2) Derived default: max(base_resolution_nm * UNWRAP_LOOP_TOL_FACTOR, 1e-6)
        """
        raw = os.environ.get(GPE_UNWRAP_LOOP_TOL_NM_ENV, "").strip()
        if raw:
            try:
                value = float(raw)
                if math.isfinite(value) and value > 0:
                    return value
            except ValueError:
                pass
            print(
                f"  [WARN] Invalid {GPE_UNWRAP_LOOP_TOL_NM_ENV}={raw!r}; "
                "falling back to derived tolerance."
            )

        resolution = (
            base_resolution_nm
            if base_resolution_nm is not None and base_resolution_nm > 0
            else DEFAULT_GRO_COORD_RESOLUTION_NM
        )
        return max(resolution * UNWRAP_LOOP_TOL_FACTOR, 1e-6)

    def _unwrap_molecule_bond_graph(
        self,
        atoms: List[GROAtom],
        box: Tuple[float, float, float],
        bond_graph: Dict[int, Set[int]],
        tol_nm: Optional[float] = None,
    ) -> Tuple[Optional[List[Tuple[float, float, float]]], bool, Optional[str]]:
        """
        Unwrap coordinates via bond-graph traversal and detect topology percolation.

        For revisits through alternate graph paths, inconsistent implied image vectors
        indicate percolation/topological wrapping and the dipole check is skipped.
        """
        if not atoms:
            return [], False, None

        natoms = len(atoms)
        Lx, Ly, Lz = box
        adjacency: Dict[int, Set[int]] = {i: set() for i in range(1, natoms + 1)}
        for i, neighbors in bond_graph.items():
            if i < 1 or i > natoms:
                continue
            for j in neighbors:
                if j < 1 or j > natoms or i == j:
                    continue
                adjacency[i].add(j)
                adjacency[j].add(i)

        if not any(adjacency.values()):
            coords, _, _ = self._unwrap_molecule_pbc_with_span(atoms, box, early_exit_ratio=None)
            return coords, False, None

        unwrapped: Dict[int, Tuple[float, float, float]] = {}
        images: Dict[int, Tuple[int, int, int]] = {}
        tol = tol_nm if tol_nm is not None and tol_nm > 0 else 1e-6

        for root in range(1, natoms + 1):
            if root in unwrapped:
                continue
            root_atom = atoms[root - 1]
            unwrapped[root] = (root_atom.x, root_atom.y, root_atom.z)
            images[root] = (0, 0, 0)
            queue = deque([root])

            while queue:
                i = queue.popleft()
                xi, yi, zi = unwrapped[i]
                img_i = images[i]
                atom_i = atoms[i - 1]

                for j in sorted(adjacency.get(i, set())):
                    atom_j = atoms[j - 1]
                    dx_w = atom_j.x - atom_i.x
                    dy_w = atom_j.y - atom_i.y
                    dz_w = atom_j.z - atom_i.z

                    sx = int(round(dx_w / Lx)) if Lx > 0 else 0
                    sy = int(round(dy_w / Ly)) if Ly > 0 else 0
                    sz = int(round(dz_w / Lz)) if Lz > 0 else 0

                    dx = dx_w - sx * Lx
                    dy = dy_w - sy * Ly
                    dz = dz_w - sz * Lz

                    implied_pos = (xi + dx, yi + dy, zi + dz)
                    implied_img = (img_i[0] - sx, img_i[1] - sy, img_i[2] - sz)

                    if j not in unwrapped:
                        unwrapped[j] = implied_pos
                        images[j] = implied_img
                        queue.append(j)
                        continue

                    existing_pos = unwrapped[j]
                    existing_img = images[j]
                    if implied_img != existing_img:
                        diff_img = (
                            implied_img[0] - existing_img[0],
                            implied_img[1] - existing_img[1],
                            implied_img[2] - existing_img[2],
                        )
                        if diff_img != (0, 0, 0):
                            return None, True, f"image_loop_shift={diff_img}"

                    dx_pos = implied_pos[0] - existing_pos[0]
                    dy_pos = implied_pos[1] - existing_pos[1]
                    dz_pos = implied_pos[2] - existing_pos[2]
                    if abs(dx_pos) <= tol and abs(dy_pos) <= tol and abs(dz_pos) <= tol:
                        continue

                    kx = int(round(dx_pos / Lx)) if Lx > 0 else 0
                    ky = int(round(dy_pos / Ly)) if Ly > 0 else 0
                    kz = int(round(dz_pos / Lz)) if Lz > 0 else 0
                    residual_ok = (
                        abs(dx_pos - kx * Lx) <= tol
                        and abs(dy_pos - ky * Ly) <= tol
                        and abs(dz_pos - kz * Lz) <= tol
                    )
                    if residual_ok and (kx, ky, kz) != (0, 0, 0):
                        return None, True, f"loop_lattice_shift={(kx, ky, kz)}"
                    return None, True, "inconsistent_bond_loop"

        coords: List[Tuple[float, float, float]] = []
        for idx in range(1, natoms + 1):
            coords.append(unwrapped.get(idx, (atoms[idx - 1].x, atoms[idx - 1].y, atoms[idx - 1].z)))
        return coords, False, None

    def _wrapped_span(self, atoms: List[GROAtom]) -> Tuple[float, float, float]:
        """Compute axis-aligned span from wrapped coordinates without allocations."""
        if not atoms:
            return (0.0, 0.0, 0.0)
        min_x = max_x = atoms[0].x
        min_y = max_y = atoms[0].y
        min_z = max_z = atoms[0].z
        for atom in atoms[1:]:
            x = atom.x
            y = atom.y
            z = atom.z
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)
            min_z = min(min_z, z)
            max_z = max(max_z, z)
        return (max_x - min_x, max_y - min_y, max_z - min_z)
    
    def _dipole_vector_unwrapped(
        self, coords: List[Tuple[float, float, float]], charges: List[float]
    ) -> Tuple[float, float, float]:
        """
        Compute dipole vector (Debye) using geometric center.
        
        v4 Issue 7: Returns vector components for proper delta calculation.
        
        Args:
            coords: List of (x, y, z) tuples in nm (should be PBC-unwrapped)
            charges: List of partial charges (e)
            
        Returns:
            Dipole vector (mx, my, mz) in Debye
        """
        if not coords or len(coords) != len(charges):
            return (0.0, 0.0, 0.0)
        
        # Geometric center for origin
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        
        # Compute dipole components relative to geometric center
        mux = muy = muz = 0.0
        for (x, y, z), q in zip(coords, charges):
            mux += q * (x - cx)
            muy += q * (y - cy)
            muz += q * (z - cz)
        
        # Convert from e*nm to Debye: 1 e*nm = 48.0321 Debye
        return (mux * 48.0321, muy * 48.0321, muz * 48.0321)

    def _dipole_magnitude(self, coords: List[Tuple[float, float, float]], charges: List[float]) -> float:
        """Compute dipole magnitude (Debye) using molecule-centered coordinates. Legacy."""
        return self._dipole_magnitude_unwrapped(coords, charges)

    def _dipole_magnitude_unwrapped(
        self, coords: List[Tuple[float, float, float]], charges: List[float]
    ) -> float:
        """
        Compute dipole magnitude (Debye) using geometric center.
        
        Physics Note (Issue E hardening):
        ---------------------------------
        The dipole moment 渭 = 危 q_i * r_i is ORIGIN-DEPENDENT for charged molecules.
        For a molecule with net charge Q 鈮?0, shifting the origin by vector R changes
        the dipole by 螖Q * R.
        
        For NEUTRAL molecules (which we're correcting toward), the dipole is
        origin-independent, so any choice of center gives the same result.
        
        We use GEOMETRIC CENTER because:
        1. It minimizes numerical error from large |r| values
        2. It provides consistent comparison before/after charge correction
        3. For the PURPOSE of checking dipole SHIFT (螖渭 = 渭_after - 渭_before),
           the origin choice CANCELS OUT as long as it's consistent
        
        Alternative: Charge-weighted center could be used if comparing absolute
        dipole magnitudes to experimental values, but for SHIFT calculations
        (our use case), geometric center is equivalent and simpler.
        
        Units:
        - Input coords: nm (GROMACS native)
        - Input charges: elementary charges (e)
        - Output: Debye (1 e路nm = 48.0321 Debye)
        
        Args:
            coords: List of (x, y, z) tuples in nm (should be PBC-unwrapped)
            charges: List of partial charges (e)
            
        Returns:
            Dipole magnitude in Debye
        """
        if not coords or len(coords) != len(charges):
            return 0.0
        
        # Geometric center for origin
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        
        # Compute dipole components relative to geometric center
        mux = muy = muz = 0.0
        for (x, y, z), q in zip(coords, charges):
            mux += q * (x - cx)
            muy += q * (y - cy)
            muz += q * (z - cz)
        
        mu_nm = math.sqrt(mux * mux + muy * muy + muz * muz)
        # Convert from e*nm to Debye: 1 e*nm = 48.0321 Debye
        return mu_nm * 48.0321

    def _grompp_preprocess_workdir(self, output_dir: Path) -> Path:
        """Dedicated debug workdir for sanitizer grompp -pp calls."""
        return output_dir / "grompp_preprocess"

    def _grompp_preprocessed_top_path(self, output_dir: Path) -> Path:
        """Path of preprocessed topology emitted by sanitizer grompp -pp helper."""
        return self._grompp_preprocess_workdir(output_dir) / "system.preprocessed.top"

    def _tail_text_file(self, path: Path, max_lines: int = 50) -> List[str]:
        """Read a bounded tail from a text file without loading it entirely."""
        if not path.exists():
            return []
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                return list(deque(handle, maxlen=max_lines))
        except OSError:
            return []

    def _expected_atoms_from_grompp(
        self,
        gmx_exe: str,
        system_top_path: Path,
        gro_path: Path,
        output_dir: Path,
        ctx: "PipelineContext",
    ) -> Optional[int]:
        """
        Run grompp -pp to preprocess topology and parse expected atom count.
        
        Hardening v6: Complexity-based timeout estimation with retry logic.
        Hardening v8: Use streamed log files (no capture_output in memory).
        """
        workdir = self._grompp_preprocess_workdir(output_dir)
        workdir.mkdir(parents=True, exist_ok=True)
        mdp_path = self._select_mdp_for_preprocess(ctx, workdir)
        pp_path = self._grompp_preprocessed_top_path(output_dir)
        tpr_path = workdir / "preprocess.tpr"
        resolved_gmx_cmd = str(gmx_exe).strip() if gmx_exe is not None else ""
        if not resolved_gmx_cmd:
            resolved_gmx_cmd = resolve_gmx_command(ctx)
        gmx_tokens = shlex.split(resolved_gmx_cmd) if resolved_gmx_cmd else []
        if not gmx_tokens:
            print("  [WARN] Empty GROMACS command for grompp preprocess")
            return None
        cmd = gmx_tokens + [
            "grompp",
            "-f",
            str(mdp_path),
            "-c",
            str(gro_path),
            "-p",
            str(system_top_path),
            "-pp",
            str(pp_path),
            "-o",
            str(tpr_path),
            "-maxwarn",
            str(ctx.grompp_maxwarn),
        ]
        
        # Hardening v6: Complexity-based timeout estimation
        complexity_metrics = self._estimate_topology_complexity(system_top_path, ctx)
        
        # Base timeout from ctx > env var > default
        base_timeout_s = getattr(ctx, "grompp_preprocess_timeout_s", None)
        if base_timeout_s is None:
            env_timeout = os.environ.get("GMX_GROMPP_PP_TIMEOUT_S", "")
            if env_timeout:
                try:
                    base_timeout_s = int(env_timeout)
                except ValueError:
                    base_timeout_s = DEFAULT_GROMPP_PP_TIMEOUT_S
            else:
                base_timeout_s = DEFAULT_GROMPP_PP_TIMEOUT_S
        
        # Compute complexity-adjusted timeout
        timeout_s = GROMPP_TIMEOUT_FLOOR_S
        timeout_s += complexity_metrics.get("lines", 0) // 1000 * GROMPP_TIMEOUT_PER_1K_LINES
        timeout_s += complexity_metrics.get("includes", 0) // 10 * GROMPP_TIMEOUT_PER_10_INCLUDES
        
        # Add GRO atom-based scaling
        gro_atoms = self._count_atoms_in_gro(gro_path)
        if gro_atoms is not None and gro_atoms > 10000:
            extra_time = ((gro_atoms - 10000) // 10000) * 60
            if gro_atoms > 100000:
                extra_time = int(extra_time * 1.5) + 600
            timeout_s += extra_time
        
        # Never reduce below user/env-set timeout or floor
        timeout_s = max(timeout_s, base_timeout_s, GROMPP_TIMEOUT_FLOOR_S)
        
        # Hardening v6: Run with retry logic
        retry_count = 0
        max_retries = 1  # Try once, then retry with longer timeout
        current_timeout = timeout_s
        result: Optional[subprocess.CompletedProcess] = None
        attempt_logs: List[Dict[str, str]] = []
        
        while retry_count <= max_retries:
            attempt_idx = retry_count + 1
            stdout_log = workdir / f"grompp_pp_attempt{attempt_idx}.stdout.log"
            stderr_log = workdir / f"grompp_pp_attempt{attempt_idx}.stderr.log"
            attempt_logs.append(
                {"stdout": str(stdout_log), "stderr": str(stderr_log)}
            )
            try:
                with open(stdout_log, "w", encoding="utf-8", errors="replace") as out_h, open(
                    stderr_log, "w", encoding="utf-8", errors="replace"
                ) as err_h:
                    result = subprocess.run(
                        cmd,
                        cwd=str(workdir),
                        stdout=out_h,
                        stderr=err_h,
                        text=False,
                        timeout=current_timeout,
                    )
                break  # Success
            except subprocess.TimeoutExpired:
                retry_count += 1
                stderr_tail = "".join(self._tail_text_file(stderr_log, max_lines=50)).strip()
                if retry_count <= max_retries:
                    # Retry with longer timeout
                    current_timeout = int(current_timeout * GROMPP_TIMEOUT_RETRY_MULTIPLIER)
                    print(
                        f"  [WARN] grompp preprocessing timed out after {timeout_s}s, "
                        f"retrying with {current_timeout}s timeout..."
                    )
                    continue
                else:
                    # Final failure
                    msg = (
                        f"grompp preprocessing timed out after {current_timeout}s (including retry).\n"
                        f"  Command: {' '.join(cmd)}\n"
                        f"  Complexity metrics: {complexity_metrics}\n"
                        f"  GRO atoms: {gro_atoms}\n"
                        f"  Stderr log: {stderr_log}\n"
                        f"  Stderr tail:\n{stderr_tail or '  (no stderr tail)'}\n"
                        f"  To increase timeout:\n"
                        f"    - Set environment variable: GMX_GROMPP_PP_TIMEOUT_S={int(current_timeout * 2)}\n"
                        f"    - Or add grompp_preprocess_timeout_s to context\n"
                        f"  Possible causes:\n"
                        f"    - Very large system (>100k atoms)\n"
                        f"    - Complex topology with many includes\n"
                        f"    - Slow filesystem or network-mounted storage"
                    )
                    print(f"  [ERROR] {msg}")
                    if ctx.strict_gro_top_check:
                        raise SanitizerError(
                            f"grompp preprocessing timed out (strict mode). {msg}"
                        )
                    return None
            except Exception as e:
                stderr_tail = "".join(self._tail_text_file(stderr_log, max_lines=50)).strip()
                print(
                    f"  [WARN] grompp preprocess failed: {e}\n"
                    f"         stderr log: {stderr_log}\n"
                    f"         stderr tail: {stderr_tail or '(no stderr tail)'}"
                )
                if ctx.strict_gro_top_check:
                    raise SanitizerError(f"grompp preprocess failed for GRO/TOP validation: {e}")
                return None
        
        if result is None:
            return None

        if result.returncode != 0 or not pp_path.exists():
            last_logs = attempt_logs[-1] if attempt_logs else {}
            stderr_path = Path(last_logs.get("stderr", "")) if last_logs.get("stderr") else None
            stderr_tail = (
                "".join(self._tail_text_file(stderr_path, max_lines=50)).strip()
                if stderr_path is not None
                else ""
            )
            print(
                "  [WARN] grompp preprocess failed "
                f"(exit={result.returncode}, pp_exists={pp_path.exists()}).\n"
                f"         stderr log: {stderr_path}\n"
                f"         stderr tail:\n{stderr_tail or '  (no stderr tail)'}"
            )
            if ctx.strict_gro_top_check:
                raise SanitizerError(
                    "grompp preprocess failed for GRO/TOP validation "
                    f"(exit={result.returncode})."
                )
            return None
        
        expected_atoms = self._parse_preprocessed_top(pp_path)
        if expected_atoms is None:
            print("  [WARN] Could not parse expected atoms from preprocessed topology.")
        
        # Hardening v6: Record complexity metrics in manifest
        if ctx.manifest:
            ctx.manifest.set_sanitizer_output("grompp_preprocess", {
                "timeout_used_s": current_timeout,
                "retry_count": retry_count,
                "complexity_metrics": complexity_metrics,
                "gro_atoms": gro_atoms,
                "expected_atoms": expected_atoms,
                "gmx_cmd": resolved_gmx_cmd,
                "debug_logs": attempt_logs,
                "debug_logs_excluded_from_fingerprints": True,
            })

        return expected_atoms
    
    def _estimate_topology_complexity(
        self, top_path: Path, ctx: "PipelineContext"
    ) -> Dict[str, int]:
        """
        Estimate topology complexity for timeout calculation.
        
        Hardening v6: Counts lines, includes, and other complexity indicators.
        
        Returns:
            Dict with 'lines', 'includes', 'moleculetypes' counts
        """
        metrics = {"lines": 0, "includes": 0, "moleculetypes": 0}
        
        if not top_path.exists():
            return metrics

        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return metrics
        lines = content.splitlines()
        metrics["lines"] = len(lines)

        for line in lines:
            if line.strip().startswith("#include"):
                metrics["includes"] += 1
            section = parse_section_name(line)
            if section == "moleculetype":
                metrics["moleculetypes"] += 1
        
        return metrics

    def _select_mdp_for_preprocess(self, ctx: "PipelineContext", workdir: Path) -> Path:
        """Select an existing MDP for preprocessing or create a minimal one."""
        mdp_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs", "mdp")
        if mdp_dir.exists():
            mdp_files = self._sorted_paths_casefold(mdp_dir.glob("*.mdp"))
            for preferred in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]:
                for mdp in mdp_files:
                    if mdp.name.lower() == preferred:
                        return mdp
            if mdp_files:
                return mdp_files[0]
        minimal_mdp = workdir / "minimal.mdp"
        self._atomic_write_text(
            minimal_mdp,
            "integrator = md\n"
            "nsteps = 0\n"
            "cutoff-scheme = Verlet\n"
            "coulombtype = PME\n"
            "rcoulomb = 1.0\n"
            "vdwtype = Cut-off\n"
            "rvdw = 1.0\n"
            "pbc = xyz\n"
        )
        return minimal_mdp

    def _parse_preprocessed_top(self, top_path: Path) -> Optional[int]:
        """Parse a preprocessed topology to compute expected atom count."""
        if not top_path.exists():
            return None
        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return None
        current_section = ""
        current_mol = None
        mol_atom_counts: Dict[str, int] = {}
        molecules_counts: Dict[str, int] = {}
        for line in content.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if current_section == "moleculetype":
                    current_mol = None
                continue
            if current_section == "moleculetype":
                data = stripped.split(";")[0].strip()
                if data:
                    parts = data.split()
                    if parts:
                        current_mol = parts[0]
                        mol_atom_counts.setdefault(current_mol, 0)
                continue
            if current_section == "atoms" and current_mol:
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                parts = data.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    mol_atom_counts[current_mol] += 1
                continue
            if current_section == "molecules":
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                parts = data.split()
                if len(parts) >= 2:
                    try:
                        molecules_counts[parts[0]] = molecules_counts.get(parts[0], 0) + int(parts[1])
                    except ValueError:
                        continue
        if not molecules_counts or not mol_atom_counts:
            return None
        total = 0
        for mol_name, count in molecules_counts.items():
            atom_count = mol_atom_counts.get(mol_name)
            if atom_count is None:
                return None
            total += atom_count * count
        return total

    def _parse_preprocessed_molecules(self, top_path: Path) -> List[Tuple[str, int]]:
        """Parse ordered [molecules] counts from a preprocessed topology."""
        if not top_path.exists():
            return []
        content = self.safe_read_text(top_path, strict=False)
        if not content:
            return []

        ordered: List[Tuple[str, int]] = []
        seen_indices: Dict[str, int] = {}
        in_molecules = False
        for line in content.splitlines():
            stripped = line.strip()
            section = parse_section_name(line)
            if section is not None:
                in_molecules = (section == "molecules")
                continue
            if not in_molecules:
                continue
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            data = strip_gmx_comment(stripped).strip()
            if not data:
                continue
            parts = data.split()
            if len(parts) < 2:
                continue
            try:
                name = parts[0]
                count = int(parts[1])
            except ValueError:
                continue
            if name in seen_indices:
                idx = seen_indices[name]
                old_name, old_count = ordered[idx]
                ordered[idx] = (old_name, old_count + count)
            else:
                seen_indices[name] = len(ordered)
                ordered.append((name, count))
        return ordered
    

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
    
    def _validate_gro_top_consistency(
        self,
        gro_path: Path,
        system_top_path: Optional[Path],
        molecule_counts: Dict[str, int],
        molecule_itp_paths: Dict[str, Path],
        ctx: "PipelineContext",
        output_dir: Path,
    ) -> bool:
        """
        Validate that GRO file atom count matches expected from topology [molecules].
        
        This prevents downstream grompp failures due to coordinate/topology mismatch.
        
        Args:
            gro_path: Path to GRO file
            molecule_counts: Dict of molecule_name -> count from [molecules]
            molecule_itp_paths: Dict of molecule_name -> itp_path
            ctx: Pipeline context for logging
            
        Returns:
            True if consistent, False if mismatch detected (caller should fail)
        """
        gro_atoms = self._count_atoms_in_gro(gro_path)
        if gro_atoms is None:
            print(f"  [WARN] Could not read atom count from {gro_path}")
            return True  # Can't validate, allow to proceed
        
        # Option A: use grompp -pp if available to preprocess topology
        expected_atoms = None
        gmx_exe = resolve_gmx_command(ctx)
        if system_top_path:
            expected_atoms = self._expected_atoms_from_grompp(
                gmx_exe=gmx_exe,
                system_top_path=system_top_path,
                gro_path=gro_path,
                output_dir=output_dir,
                ctx=ctx,
            )
        
        # Option B: fallback parser
        if expected_atoms is None:
            expected_atoms, uncertain = self._expected_atoms_from_itps(
                molecule_counts, molecule_itp_paths
            )
            if uncertain:
                msg = "Preprocessor directives detected in [ atoms ]; GRO/TOP check is uncertain."
                if ctx.strict_gro_top_check:
                    raise SanitizerError(msg)
                print(f"  [WARN] {msg} Skipping GRO/TOP check.")
                return True
            if expected_atoms is None:
                print("  [WARN] Could not determine expected atom count from ITPs.")
                return True
        
        if gro_atoms != expected_atoms:
            print(f"  [ERROR] GRO/TOP atom count mismatch!")
            print(f"         GRO file: {gro_atoms} atoms")
            print(f"         Expected from [molecules]: {expected_atoms} atoms")
            print(f"         Molecule counts: {molecule_counts}")
            print(f"         This usually means the coordinate file doesn't match the topology.")
            return False
        
        print(f"  [OK] GRO/TOP consistency check: {gro_atoms} atoms match")
        return True
    
    def _expected_atoms_from_itps(
        self,
        molecule_counts: Dict[str, int],
        molecule_itp_paths: Dict[str, Path],
    ) -> Tuple[Optional[int], bool]:
        """Compute expected atom count from ITPs (fallback, preprocessor-blind)."""
        expected_atoms = 0
        uncertain = False
        unknown = []
        for mol_name, mol_count in molecule_counts.items():
            if mol_count <= 0:
                continue
            itp_path = molecule_itp_paths.get(mol_name)
            if not itp_path:
                unknown.append(mol_name)
                continue
            count, has_preproc = self._count_atoms_in_itp(itp_path, mol_name)
            if has_preproc:
                uncertain = True
            if count is None:
                unknown.append(mol_name)
                continue
            expected_atoms += count * mol_count
        if unknown:
            print(f"  [WARN] Cannot validate atom counts for molecules: {unknown}")
            return None, uncertain
        return expected_atoms, uncertain

    def _count_atoms_in_itp(self, itp_path: Path, mol_name: str) -> Tuple[Optional[int], bool]:
        """
        Count atoms in a molecule ITP file from the [atoms] section.
        
        Returns:
            Tuple of (atom_count, has_preprocessor_directives)
        """
        if not itp_path.exists():
            return None, False
        content = self.safe_read_text(itp_path, strict=False)
        if not content:
            return None, False
        current_section = ""
        current_mol = None
        atom_count = 0
        has_preproc = False
        for line in content.splitlines():
            stripped = line.strip()
            if stripped.startswith("#"):
                if current_section == "atoms":
                    has_preproc = True
                continue
            section = parse_section_name(line)
            if section is not None:
                current_section = section
                if current_section == "moleculetype":
                    current_mol = None
                continue
            if current_section == "moleculetype":
                data = stripped.split(";")[0].strip()
                if data and not data.startswith(";"):
                    parts = data.split()
                    if parts:
                        current_mol = parts[0]
                continue
            if current_section == "atoms":
                if current_mol != mol_name:
                    continue
                data = stripped.split(";")[0].strip()
                if not data:
                    continue
                parts = data.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    atom_count += 1
        return (atom_count if atom_count > 0 else None), has_preproc

__all__ = ["SanitizerStage", "SanitizerError", "GROAtom", "GROParseResult"]

