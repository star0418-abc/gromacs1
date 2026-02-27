"""
PACKMOL input file generator.

Generates packmol.inp from composition and box configuration.
Resolves molecule templates from IN/molecules/<NAME>/pdb/ (primary) or gro/ (fallback).
"""

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING
import shutil

if TYPE_CHECKING:
    from .context import PipelineContext
    from .composition import CompositionResult, BoxConfig


# =============================================================================
# Retry Diagnostics
# =============================================================================

@dataclass
class RetryAttemptInfo:
    """Diagnostic info for a single PACKMOL attempt."""
    attempt_index: int
    rho_target: float  # Density used for this attempt (g/cm³)
    box_length_nm: float  # Computed box length for this attempt
    tolerance: float  # PACKMOL tolerance (Å)
    seed: int  # Random seed used
    composition_changed: bool  # Whether molecule counts were changed
    success: bool  # Whether this attempt succeeded
    inp_file: str  # Path to input file for this attempt
    pdb_file: str  # Path to output PDB for this attempt (if generated)
    error_reason: Optional[str] = None  # Error reason if failed
    seed_source: str = "provided"  # How the seed was resolved


def _normalize_name(value: str) -> str:
    """Normalize molecule-like names for robust matching."""
    return "".join(ch for ch in value.upper() if ch.isalnum())


POLYMER_NAME_TOKENS = {
    "PEO", "PEG", "PVDF", "PAM", "PEGDMA", "OEGMA", "PMMA", "PVA",
    "POLYMER", "CHAIN", "NETWORK", "CROSSLINK", "MONOMER",
}
POLYMER_MW_THRESHOLD = 1000.0
PACKMOL_SUCCESS_MARKERS = (
    "Success!",
    "PACKMOL Ended Successfully",
)
MIN_INSIDE_BOX_SPAN_ANG = 1e-3


def _is_polymer_like(mol_name: str, mw_g_mol: float) -> bool:
    """Best-effort polymer-like classification for template-fit diagnostics."""
    if mw_g_mol >= POLYMER_MW_THRESHOLD:
        return True
    return any(token in mol_name.upper() for token in POLYMER_NAME_TOKENS)


def _contains_packmol_success(stdout: str, stderr: str) -> bool:
    """Return True when either stream contains a PACKMOL success marker."""
    return any(marker in stdout or marker in stderr for marker in PACKMOL_SUCCESS_MARKERS)


def _derive_deterministic_packmol_seed(seed_context: str, attempt_index: int) -> int:
    """Derive a deterministic PACKMOL seed in [1, 2_000_000_000]."""
    seed_key = f"{seed_context}|attempt={attempt_index}"
    digest = hashlib.sha256(seed_key.encode("utf-8")).hexdigest()
    return 1 + (int(digest[:8], 16) % 2_000_000_000)


def _validate_inside_box_geometry(box_length_ang: float, margin: float) -> float:
    """Validate that margin leaves a non-degenerate PACKMOL inside-box span."""
    usable_span = box_length_ang - 2.0 * margin
    if usable_span < MIN_INSIDE_BOX_SPAN_ANG:
        raise ValueError(
            "Invalid PACKMOL inside-box region: usable span is too small after margins. "
            f"box_length_ang={box_length_ang:.6f}, margin={margin:.6f}, "
            f"usable_span={usable_span:.6f}, required_min_span={MIN_INSIDE_BOX_SPAN_ANG:.6f}"
        )
    return usable_span


def _safe_inside_box(
    lo: Tuple[float, float, float],
    hi: Tuple[float, float, float],
    margin: float,
    box_length_ang: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Clamp a requested inside-box region to a valid PACKMOL region.

    Ensures the region stays within [margin, box_length_ang - margin] and has
    strictly positive span on all axes.
    """
    _validate_inside_box_geometry(box_length_ang=box_length_ang, margin=margin)
    lower_bound = margin
    upper_bound = box_length_ang - margin

    x0 = max(lower_bound, min(lo[0], upper_bound - MIN_INSIDE_BOX_SPAN_ANG))
    y0 = max(lower_bound, min(lo[1], upper_bound - MIN_INSIDE_BOX_SPAN_ANG))
    z0 = max(lower_bound, min(lo[2], upper_bound - MIN_INSIDE_BOX_SPAN_ANG))
    x1 = min(upper_bound, max(hi[0], x0 + MIN_INSIDE_BOX_SPAN_ANG))
    y1 = min(upper_bound, max(hi[1], y0 + MIN_INSIDE_BOX_SPAN_ANG))
    z1 = min(upper_bound, max(hi[2], z0 + MIN_INSIDE_BOX_SPAN_ANG))
    return (x0, y0, z0, x1, y1, z1)


def _full_box_region(box_length_ang: float, margin: float) -> Tuple[float, float, float, float, float, float]:
    """Return the default full-box region with edge margin."""
    _validate_inside_box_geometry(box_length_ang=box_length_ang, margin=margin)
    box_max = box_length_ang - margin
    return (margin, margin, margin, box_max, box_max, box_max)


def _central_box_region(
    box_length_ang: float,
    margin: float,
    fraction: float,
) -> Tuple[float, float, float, float, float, float]:
    """Return a centered sub-box occupying `fraction` of the usable span."""
    usable = _validate_inside_box_geometry(box_length_ang=box_length_ang, margin=margin)
    frac = max(0.15, min(1.0, fraction))
    half_span = 0.5 * usable * frac
    center = 0.5 * box_length_ang
    lo = (center - half_span, center - half_span, center - half_span)
    hi = (center + half_span, center + half_span, center + half_span)
    return _safe_inside_box(lo, hi, margin=margin, box_length_ang=box_length_ang)


def _z_slab_region(
    box_length_ang: float,
    margin: float,
    z_start_fraction: float,
) -> Tuple[float, float, float, float, float, float]:
    """Return a slab that spans the upper z part of the box."""
    full = _full_box_region(box_length_ang, margin)
    z0_full = full[2]
    z1_full = full[5]
    span = max(1e-3, z1_full - z0_full)
    z_start = z0_full + max(0.0, min(0.95, z_start_fraction)) * span
    lo = (full[0], full[1], z_start)
    hi = (full[3], full[4], full[5])
    return _safe_inside_box(lo, hi, margin=margin, box_length_ang=box_length_ang)


def _build_structure_blocks(
    counts: Dict[str, int],
    box_length_ang: float,
    margin: float,
    placement_hints: Optional[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Build per-structure placement blocks for PACKMOL, optionally with hints.

    Returns:
        (blocks, hint_diagnostics)
    """
    mode = str((placement_hints or {}).get("mode", "none")).strip() or "none"
    li_names = {_normalize_name(v) for v in (placement_hints or {}).get("li_names", [])}
    polymer_names = {_normalize_name(v) for v in (placement_hints or {}).get("polymer_names", [])}
    solvent_names = {_normalize_name(v) for v in (placement_hints or {}).get("solvent_names", [])}
    tfsi_names = {_normalize_name(v) for v in (placement_hints or {}).get("tfsi_names", [])}

    li_fraction = float((placement_hints or {}).get("li_fraction", 1.0))
    li_fraction = max(0.0, min(1.0, li_fraction))
    tightening = int((placement_hints or {}).get("tightening", 0))

    full_region = _full_box_region(box_length_ang, margin)
    blocks: List[Dict[str, Any]] = []
    hint_diag: Dict[str, Any] = {
        "mode": mode,
        "li_fraction": li_fraction,
        "tightening": tightening,
        "applied": False,
        "regions": [],
    }

    # Tightening provides a deterministic way to nudge retry constraints.
    li_core_fraction = max(0.45, 0.72 - 0.04 * tightening)
    polymer_core_fraction = max(0.58, 0.82 - 0.03 * tightening)
    solvent_core_fraction = max(0.62, 0.88 - 0.03 * tightening)
    tfsi_slab_start = min(0.82, 0.55 + 0.05 * tightening)

    for mol_name in sorted(counts.keys()):
        count = int(counts[mol_name])
        if count <= 0:
            continue

        mol_norm = _normalize_name(mol_name)
        entries: List[Tuple[int, Tuple[float, float, float, float, float, float], str]] = [
            (count, full_region, "default"),
        ]

        if mode == "li_near_polymer":
            if mol_norm in polymer_names:
                entries = [
                    (count, _central_box_region(box_length_ang, margin, polymer_core_fraction), "polymer_core"),
                ]
            elif mol_norm in li_names:
                li_core_count = int(round(count * li_fraction))
                li_core_count = max(0, min(count, li_core_count))
                li_rest = count - li_core_count
                entries = []
                if li_core_count > 0:
                    entries.append(
                        (
                            li_core_count,
                            _central_box_region(box_length_ang, margin, li_core_fraction),
                            "li_core",
                        )
                    )
                if li_rest > 0:
                    entries.append((li_rest, full_region, "li_bulk"))
        elif mode == "li_solvent_shell":
            if mol_norm in li_names:
                entries = [
                    (count, _central_box_region(box_length_ang, margin, li_core_fraction), "li_core"),
                ]
            elif mol_norm in solvent_names:
                entries = [
                    (count, _central_box_region(box_length_ang, margin, solvent_core_fraction), "solvent_core"),
                ]
            elif mol_norm in tfsi_names:
                entries = [
                    (count, _z_slab_region(box_length_ang, margin, tfsi_slab_start), "tfsi_slab"),
                ]

        for block_count, region, tag in entries:
            if block_count <= 0:
                continue
            blocks.append(
                {
                    "molecule": mol_name,
                    "count": block_count,
                    "region": region,
                    "tag": tag,
                }
            )
            hint_diag["regions"].append(
                {
                    "molecule": mol_name,
                    "count": block_count,
                    "tag": tag,
                    "inside_box_ang": [round(v, 4) for v in region],
                }
            )

    hint_diag["applied"] = any(item["tag"] != "default" for item in hint_diag["regions"])
    return blocks, hint_diag


def _parse_template_extent_nm(
    template_path: Path,
    template_format: str,
    pdb_scale_nm: float,
) -> Tuple[float, int]:
    """
    Compute max axis span in nm for a template.

    Returns:
        (span_nm, atom_count)
    """
    coords: List[Tuple[float, float, float]] = []
    fmt = template_format.lower()
    if fmt == "gro":
        with open(template_path, "r") as handle:
            lines = handle.readlines()
        for line in lines[2:]:
            if len(line) < 44:
                continue
            try:
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
            except ValueError:
                continue
            coords.append((x, y, z))
    else:
        with open(template_path, "r") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                try:
                    x = float(line[30:38]) * pdb_scale_nm
                    y = float(line[38:46]) * pdb_scale_nm
                    z = float(line[46:54]) * pdb_scale_nm
                except ValueError:
                    continue
                coords.append((x, y, z))
    if not coords:
        return 0.0, 0
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    span_nm = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    return span_nm, len(coords)


def _compute_template_fit_diagnostics(
    composition: "CompositionResult",
    templates: Dict[str, Tuple[Path, str]],
    box_length_nm: float,
    margin_nm: float,
    pdb_scale_override: Optional[float] = None,
) -> Dict[str, Any]:
    """Compute per-template fit diagnostics against the intended box."""
    rows: List[Dict[str, Any]] = []
    limiting_rows: List[Dict[str, Any]] = []
    limiting_polymer_rows: List[Dict[str, Any]] = []

    usable_box_nm = max(0.0, box_length_nm - 2.0 * margin_nm)
    pdb_scale_nm = pdb_scale_override if pdb_scale_override and pdb_scale_override > 0 else 0.1

    for mol_name in sorted(templates.keys()):
        template_path, template_format = templates[mol_name]
        span_nm, atom_count = _parse_template_extent_nm(
            template_path=template_path,
            template_format=template_format,
            pdb_scale_nm=pdb_scale_nm,
        )
        mw = float(composition.mw_used.get(mol_name, 0.0))
        polymer_like = _is_polymer_like(mol_name, mw)
        fit_ok = span_nm <= usable_box_nm + 1e-6 and span_nm <= box_length_nm + 1e-6
        row = {
            "name": mol_name,
            "polymer_like": polymer_like,
            "span_nm": span_nm,
            "box_nm": box_length_nm,
            "margin_nm": margin_nm,
            "usable_box_nm": usable_box_nm,
            "fit_ok": fit_ok,
            "template_format": template_format,
            "atom_count": atom_count,
            "template_path": str(template_path),
            "unit_scale_assumed_nm_per_raw": 1.0 if template_format == "gro" else pdb_scale_nm,
        }
        rows.append(row)
        if not fit_ok:
            limiting_rows.append(row)
            if polymer_like:
                limiting_polymer_rows.append(row)

    return {
        "box_nm": box_length_nm,
        "margin_nm": margin_nm,
        "usable_box_nm": usable_box_nm,
        "rows": rows,
        "limiting_rows": limiting_rows,
        "limiting_polymer_rows": limiting_polymer_rows,
        "limiting_detected": bool(limiting_rows),
        "polymer_limiting_detected": bool(limiting_polymer_rows),
    }


def _print_template_fit_table(fit_diag: Dict[str, Any]) -> None:
    """Print concise template-fit diagnostics to console."""
    rows = fit_diag.get("rows", [])
    if not rows:
        print("  [Template Fit] no templates available")
        return
    box_nm = float(fit_diag.get("box_nm", 0.0))
    margin_nm = float(fit_diag.get("margin_nm", 0.0))
    usable_nm = float(fit_diag.get("usable_box_nm", 0.0))
    print(
        "  [Template Fit] "
        f"box={box_nm:.4f} nm, margin={margin_nm:.4f} nm, usable={usable_nm:.4f} nm"
    )
    print("    name              polymer  span_nm  box_nm  margin_nm  fit_ok")
    for row in rows:
        print(
            "    "
            f"{row['name'][:16]:<16}  "
            f"{'yes' if row['polymer_like'] else 'no ':<7}  "
            f"{row['span_nm']:>7.3f}  "
            f"{row['box_nm']:>6.3f}  "
            f"{row['margin_nm']:>9.3f}  "
            f"{'yes' if row['fit_ok'] else 'no'}"
        )


# =============================================================================
# Template Resolution
# =============================================================================

def resolve_molecule_template(
    molecule_name: str,
    molecules_dir: Path,
    prefer_format: str = "pdb",
) -> Tuple[Path, str]:
    """
    Resolve the path to a molecule's structure template.
    
    Resolution order:
    1. IN/molecules/<NAME>/pdb/<NAME>.pdb (preferred)
    2. IN/molecules/<NAME>/gro/<NAME>.gro (fallback)
    
    Args:
        molecule_name: Name of the molecule (e.g., "LiTFSI")
        molecules_dir: Path to IN/molecules/ directory
        prefer_format: Preferred format ("pdb" or "gro")
        
    Returns:
        Tuple of (resolved_path, format)
        
    Raises:
        FileNotFoundError: If neither PDB nor GRO template exists
    """
    mol_dir = molecules_dir / molecule_name
    
    # Try preferred format first
    formats = ["pdb", "gro"] if prefer_format == "pdb" else ["gro", "pdb"]
    
    for fmt in formats:
        template_dir = mol_dir / fmt
        # Try exact name match first
        template_file = template_dir / f"{molecule_name}.{fmt}"
        if template_file.exists():
            return template_file, fmt
        
        # Try any file with the format extension (if only one exists)
        if template_dir.exists():
            candidates = list(template_dir.glob(f"*.{fmt}"))
            if len(candidates) == 1:
                return candidates[0], fmt
            elif len(candidates) > 1:
                raise FileNotFoundError(
                    f"Ambiguous templates for {molecule_name}: found {len(candidates)} "
                    f".{fmt} files in {template_dir}. Please ensure exactly one template."
                )
    
    raise FileNotFoundError(
        f"No template found for molecule '{molecule_name}'. "
        f"Expected: {mol_dir}/pdb/{molecule_name}.pdb or {mol_dir}/gro/{molecule_name}.gro"
    )


def resolve_all_templates(
    molecule_names: List[str],
    molecules_dir: Path,
    prefer_format: str = "pdb",
) -> Dict[str, Tuple[Path, str]]:
    """
    Resolve templates for all molecules in composition.
    
    Args:
        molecule_names: List of molecule names
        molecules_dir: Path to IN/molecules/ directory
        prefer_format: Preferred format ("pdb" or "gro")
        
    Returns:
        Dict mapping molecule name to (path, format) tuple
        
    Raises:
        FileNotFoundError: If any template is missing
    """
    templates = {}
    missing = []
    
    for name in molecule_names:
        try:
            templates[name] = resolve_molecule_template(name, molecules_dir, prefer_format)
        except FileNotFoundError as e:
            missing.append(str(e))
    
    if missing:
        raise FileNotFoundError(
            "Missing molecule templates:\n" + "\n".join(f"  - {m}" for m in missing)
        )
    
    return templates


def summarize_template_formats(
    templates: Dict[str, Tuple[Path, str]],
) -> Tuple[str, Dict[str, int]]:
    """
    Summarize template formats and counts.

    Returns:
        Tuple of (summary_format, counts_by_format)
        summary_format is one of: "pdb", "gro", "mixed", "unknown"
    """
    counts: Dict[str, int] = {}
    for _, (_, fmt) in templates.items():
        counts[fmt] = counts.get(fmt, 0) + 1

    if not counts:
        return "unknown", counts
    if len(counts) == 1:
        return next(iter(counts.keys())), counts
    return "mixed", counts


# =============================================================================
# PACKMOL Input Generation
# =============================================================================

def generate_packmol_input(
    counts: Dict[str, int],
    templates: Dict[str, Tuple[Path, str]],
    box_length_nm: float,
    output_file: Path,
    tolerance: float = 2.0,
    seed: int = -1,
    filetype: str = "pdb",
    edge_margin_nm: Optional[float] = None,
    placement_hints: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Generate PACKMOL input file content.
    
    Args:
        counts: Dict mapping molecule name to count
        templates: Dict mapping molecule name to (path, format) tuple
        box_length_nm: Cubic box side length in nanometers
        output_file: Path where PACKMOL should write the output structure
        tolerance: Minimum distance between atoms (Angstrom)
        seed: Random seed (-1 for random)
        filetype: Output file type ("pdb" or "xyz")
        placement_hints: Optional placement-bias hints (Li/polymer/solvent modes)
        
    Returns:
        PACKMOL input file content as string
    """
    # Convert nm to Angstrom for PACKMOL
    box_length_ang = box_length_nm * 10.0
    output_target = output_file.resolve().as_posix()
    
    lines = [
        f"# PACKMOL input file",
        f"# Generated for MD simulation pipeline",
        f"# Box size: {box_length_nm:.3f} nm = {box_length_ang:.2f} Å",
        f"",
        f"tolerance {tolerance}",
        f"filetype {filetype}",
        f"output {output_target}",
        f"seed {seed}",
        f"",
    ]
    
    # Use configurable margin from box edges (default: tolerance)
    # edge_margin_nm allows explicit control over boundary safety
    if edge_margin_nm is not None:
        margin = edge_margin_nm * 10.0  # Convert nm to Å for PACKMOL
    else:
        margin = tolerance

    blocks, hint_diag = _build_structure_blocks(
        counts=counts,
        box_length_ang=box_length_ang,
        margin=margin,
        placement_hints=placement_hints,
    )
    if hint_diag.get("applied"):
        lines.append(
            f"# Placement hints enabled: mode={hint_diag.get('mode')} "
            f"(tightening={hint_diag.get('tightening', 0)})"
        )
        lines.append("")

    for block in blocks:
        mol_name = block["molecule"]
        count = int(block["count"])
        x0, y0, z0, x1, y1, z1 = block["region"]
        tag = block.get("tag", "default")
        template_path, _template_fmt = templates[mol_name]
        lines.extend([
            f"# {mol_name}: {count} molecules ({tag})",
            f"structure {template_path}",
            f"  number {count}",
            f"  inside box {x0:.3f} {y0:.3f} {z0:.3f} {x1:.3f} {y1:.3f} {z1:.3f}",
            f"end structure",
            f"",
        ])
    
    return "\n".join(lines)


def write_packmol_input(
    composition: "CompositionResult",
    box: "BoxConfig", 
    output_dir: Path,
    molecules_dir: Path,
    prefer_format: str = "pdb",
    placement_hints: Optional[Dict[str, Any]] = None,
) -> Tuple[Path, Dict[str, Tuple[Path, str]]]:
    """
    Generate and write PACKMOL input file.
    
    Args:
        composition: CompositionResult with molecule counts
        box: BoxConfig with box dimensions
        output_dir: Directory to write packmol.inp
        molecules_dir: Path to IN/molecules/ directory
        prefer_format: Preferred template format
        placement_hints: Optional placement-bias hints for PACKMOL regions
        
    Returns:
        Tuple of (packmol_inp_path, templates_dict)
    """
    # Resolve templates
    templates = resolve_all_templates(
        list(composition.counts.keys()),
        molecules_dir,
        prefer_format,
    )
    
    # Generate input content
    try:
        content = generate_packmol_input(
            counts=composition.counts,
            templates=templates,
            box_length_nm=box.length_nm,
            output_file=output_dir / "initial.pdb",
            placement_hints=placement_hints,
        )
    except ValueError as exc:
        raise PackmolError(f"Invalid PACKMOL box/margin settings: {exc}") from exc
    
    # Write input file
    output_dir.mkdir(parents=True, exist_ok=True)
    inp_path = output_dir / "packmol.inp"
    inp_path.write_text(content)
    
    return inp_path, templates


# =============================================================================
# PACKMOL Execution
# =============================================================================

class PackmolError(Exception):
    """Error during PACKMOL execution."""
    pass


def run_packmol(
    inp_file: Path,
    working_dir: Path,
    timeout: int = 300,
) -> Tuple[bool, str, str]:
    """
    Execute PACKMOL.
    
    Args:
        inp_file: Path to packmol.inp
        working_dir: Working directory for execution
        timeout: Timeout in seconds
        
    Returns:
        Tuple of (success, stdout, stderr)
    """
    import subprocess
    
    # Check PACKMOL is available
    packmol_exe = shutil.which("packmol")
    if not packmol_exe:
        raise PackmolError("PACKMOL executable not found in PATH")
    
    try:
        with open(inp_file, "r") as stdin_handle:
            result = subprocess.run(
                [packmol_exe],
                stdin=stdin_handle,
                capture_output=True,
                text=True,
                cwd=working_dir,
                timeout=timeout,
            )
        
        stdout = result.stdout
        stderr = result.stderr
        
        # PACKMOL can emit completion markers to either stream.
        success = _contains_packmol_success(stdout, stderr)
        
        return success, stdout, stderr
        
    except subprocess.TimeoutExpired:
        raise PackmolError(f"PACKMOL timed out after {timeout} seconds")
    except Exception as e:
        raise PackmolError(f"PACKMOL execution failed: {e}")


def check_packmol_output(
    stdout: str,
    stderr: str,
) -> Tuple[bool, Optional[str]]:
    """
    Check PACKMOL output for success or specific error conditions.
    
    Args:
        stdout: PACKMOL stdout
        stderr: PACKMOL stderr
        
    Returns:
        Tuple of (success, error_reason)
    """
    # Success indicators
    if _contains_packmol_success(stdout, stderr):
        return True, None

    stdout_l = stdout.lower()
    stderr_l = stderr.lower()
    combined_l = f"{stdout_l}\n{stderr_l}"

    # Specific failure classes first (more actionable than generic STOP/ERROR).
    specific_error_indicators = [
        ("no such file", "file_missing"),
        ("cannot open file", "file_missing"),
        ("can't open file", "file_missing"),
        ("error opening file", "file_missing"),
        ("output file", "output_file_missing"),
        ("could not read", "template_read_error"),
        ("error reading", "template_read_error"),
        ("failed to read", "template_read_error"),
        ("structure file", "template_read_error"),
        ("inside box", "invalid_inside_box"),
        ("invalid box", "invalid_inside_box"),
        ("box limits", "invalid_inside_box"),
        ("must be greater", "invalid_inside_box"),
        ("syntax error", "input_parse_error"),
        ("unknown keyword", "input_parse_error"),
        ("floating invalid", "numeric_parse_error"),
    ]
    for indicator, reason in specific_error_indicators:
        if indicator in combined_l:
            return False, reason
    
    # Error indicators
    error_indicators = [
        ("Could not fit", "overlap_failure"),
        ("STOP", "stop_error"),
        ("ERROR", "general_error"),
        ("impossible", "impossible_packing"),
    ]
    
    for indicator, reason in error_indicators:
        indicator_l = indicator.lower()
        if indicator_l in stdout_l or indicator_l in stderr_l:
            return False, reason
    
    # Unknown failure
    return False, "unknown_failure"


# =============================================================================
# Unit Validation
# =============================================================================

def validate_pdb_box_extent(
    pdb_path: Path,
    expected_box_nm: float,
    tolerance_factor: float = 0.5,
    expected_unit: str = "A",
    min_fraction: Optional[float] = None,
) -> Tuple[bool, str]:
    """
    Validate that PDB coordinate extent matches expected box size.
    
    Detects nm/Å unit mismatches by checking if coordinates are ~10x off.
    
    Args:
        pdb_path: Path to PDB file
        expected_box_nm: Expected cubic box length in nm
        tolerance_factor: Allowed deviation factor (0.5 = 50%)
        expected_unit: "A" for Å (default) or "nm" if coordinates are expected in nm
        min_fraction: Optional explicit lower threshold as a fraction of expected box.
            When None, derived from tolerance_factor.
        
    Returns:
        Tuple of (is_valid, message)
    """
    max_coord = 0.0
    min_coord = float('inf')
    atom_count = 0
    
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    max_coord = max(max_coord, x, y, z)
                    min_coord = min(min_coord, x, y, z)
                    atom_count += 1
                except (ValueError, IndexError):
                    continue
    
    if atom_count == 0:
        return False, f"No atoms found in PDB file: {pdb_path}"
    
    extent_raw = max_coord - min_coord
    expected_unit = expected_unit.lower()
    tol = max(0.0, float(tolerance_factor))
    lower_fraction = 0.1 * tol if min_fraction is None else float(min_fraction)
    lower_fraction = max(0.0, min(0.99, lower_fraction))
    upper_fraction = max(lower_fraction, 1.0 + 0.2 * tol)
    threshold_msg = (
        f"(thresholds: lower={lower_fraction:.3f}x, upper={upper_fraction:.3f}x, "
        f"tolerance_factor={tol:.3f})"
    )

    if expected_unit == "nm":
        extent_nm = extent_raw
        if extent_nm > expected_box_nm * upper_fraction:
            return False, (
                f"PDB extent ({extent_nm:.3f} nm) exceeds expected box "
                f"({expected_box_nm:.3f} nm). Check for unit mismatch. {threshold_msg}"
            )
        if extent_nm < expected_box_nm * lower_fraction:
            return False, (
                f"PDB extent ({extent_nm:.3f} nm) is much smaller than expected box "
                f"({expected_box_nm:.3f} nm). Coordinates may be in Å. {threshold_msg}"
            )
        return True, (
            f"PDB extent {extent_nm:.3f} nm within expected box "
            f"{expected_box_nm:.3f} nm {threshold_msg}"
        )

    # Default: Å
    extent_ang = extent_raw
    expected_box_ang = expected_box_nm * 10.0

    if extent_ang > expected_box_ang * upper_fraction:
        return False, (
            f"PDB extent ({extent_ang:.1f} Å) exceeds expected box "
            f"({expected_box_ang:.1f} Å). Check for unit mismatch. {threshold_msg}"
        )

    if extent_ang < expected_box_ang * lower_fraction:
        return False, (
            f"PDB extent ({extent_ang:.1f} Å) is much smaller than expected box "
            f"({expected_box_ang:.1f} Å). Coordinates may already be in nm. {threshold_msg}"
        )

    return True, (
        f"PDB extent {extent_ang:.1f} Å within expected box "
        f"{expected_box_ang:.1f} Å {threshold_msg}"
    )


# =============================================================================
# Retry Strategy
# =============================================================================

@dataclass
class RetryResult:
    """Complete result from PACKMOL retry execution."""
    output_pdb: Path  # Final successful PDB output
    final_box: "BoxConfig"  # Box configuration used in successful attempt
    templates: Dict[str, Tuple[Path, str]]  # Resolved templates
    successful_attempt_index: int  # 0-based index of successful attempt
    attempts: List[RetryAttemptInfo]  # All attempt diagnostics
    rho0: float  # Initial target density (immutable reference)
    effective_rho_min: float  # Effective minimum density used
    composition_changed: bool  # Whether any attempt changed composition


def run_packmol_with_retry(
    composition: "CompositionResult",
    initial_rho: float,
    output_dir: Path,
    molecules_dir: Path,
    max_retries: int = 3,
    rho_decrease_factor: float = 0.95,
    prefer_format: str = "pdb",
    rho_min: Optional[float] = None,  # Absolute minimum density from config
    rho_max: Optional[float] = None,  # Maximum density from config
    edge_margin_nm: Optional[float] = None,  # Edge margin for molecule placement
    box_margin_nm: Optional[float] = None,  # Extra margin added to box size (nm)
    tolerance: float = 2.0,  # PACKMOL tolerance in Å
    seed: int = -1,  # Random seed
    density_floor_fraction: float = 0.85,  # Floor as fraction of initial density
    allow_density_reduction: bool = False,
    max_density_reduction_fraction: float = 0.10,
    pdb_scale_override: Optional[float] = None,
    placement_hints: Optional[Dict[str, Any]] = None,
    allow_random_seed: bool = False,
    seed_context: Optional[str] = None,
) -> Tuple[
    Path,
    "BoxConfig",
    Dict[str, Tuple[Path, str]],
    int,
    List[RetryAttemptInfo],
    Dict[str, Any],
]:
    """
    Run PACKMOL with retry strategy on failure.
    
    SAFEGUARDS IMPLEMENTED:
    1) 85% floor is based on INITIAL rho0, not per-attempt rho_target
    2) Each attempt gets versioned outputs: packmol_try{N}.inp AND packmol_try{N}.pdb
    3) Density schedule only when enabled: rho_i = max(rho0 * (0.95 ** attempt), effective_rho_min)
    4) Molecule counts stay CONSTANT (composition_changed=False)
    
    Default behavior keeps target density fixed across retries.
    Density backoff (larger box) is only used when allow_density_reduction=True.
    Any backoff is clamped by both absolute and relative density floors.
    
    Args:
        composition: CompositionResult with molecule counts
        initial_rho: Initial density guess (g/cm³) - this is rho0
        output_dir: Directory for PACKMOL output
        molecules_dir: Path to IN/molecules/ directory
        max_retries: Maximum retry attempts
        rho_decrease_factor: Factor to multiply density on retry (<1.0)
        prefer_format: Preferred template format
        rho_min: Absolute minimum density (g/cm³) from config
        rho_max: Maximum density (g/cm³) from config
        edge_margin_nm: Edge margin for molecule placement (nm)
        box_margin_nm: Additional margin applied to box length (nm)
        tolerance: PACKMOL tolerance in Angstroms
        seed: Base seed for PACKMOL. If >0, attempt index offsets are deterministic.
            If -1, deterministic derivation is used unless allow_random_seed=True.
        density_floor_fraction: Floor as fraction of initial density (default 0.85)
        allow_density_reduction: Allow retries to lower density/expand box
        max_density_reduction_fraction: Hard cap on total density reduction from rho0
        pdb_scale_override: Optional explicit scale for PDB template diagnostics
        placement_hints: Optional placement-bias hints for structure regions
        allow_random_seed: Keep seed=-1 random behavior (explicit opt-in, non-reproducible)
        seed_context: Optional context string for deterministic seed derivation when seed=-1
        
    Returns:
        Tuple of (output_pdb_path, final_box_config, templates_dict, retries_used, attempt_list, retry_meta)
        
    Raises:
        PackmolError: If all retries fail
    """
    from .composition import compute_box_size
    
    # === B1: Compute effective minimum density based on INITIAL rho0 ===
    rho0 = initial_rho  # Immutable reference - this is the original target
    relative_floor = rho0 * density_floor_fraction  # Relative floor from initial rho0
    reduction_cap_fraction = max(0.0, min(float(max_density_reduction_fraction), 0.95))
    cap_floor = rho0 * (1.0 - reduction_cap_fraction)
    config_floor = rho_min if rho_min is not None else 0.0
    effective_rho_min = max(config_floor, relative_floor, cap_floor)
    
    # Also respect rho_max if provided
    effective_rho_max = rho_max if rho_max is not None else float('inf')
    
    print(
        f"  [Retry Config] rho0={rho0:.4f}, floor_fraction={density_floor_fraction:.0%}, "
        f"effective_rho_min={effective_rho_min:.4f}, allow_density_reduction={allow_density_reduction}"
    )
    print(
        "                 density backoff cap: "
        f"max_reduction={reduction_cap_fraction:.1%} "
        f"(floor_from_cap={cap_floor:.4f} g/cm³)"
    )
    if rho_min is not None:
        print(
            f"                 config rho_min={rho_min:.4f}, "
            f"chosen floor={effective_rho_min:.4f} (max of absolute+relative)"
        )
    if not allow_density_reduction:
        print("                 retries keep density fixed; only seed/hints vary.")

    default_seed_context = (
        f"{output_dir.resolve().as_posix()}|"
        f"rho0={rho0:.8f}|"
        f"counts={','.join(f'{k}:{composition.counts[k]}' for k in sorted(composition.counts.keys()))}"
    )
    seed_context_key = str(seed_context).strip() if seed_context else default_seed_context
    if seed == -1 and allow_random_seed:
        print(
            "  [Retry Seed] WARNING: seed=-1 keeps PACKMOL internal randomness; "
            "attempts are not reproducible."
        )
        seed_mode = "random_non_reproducible"
    elif seed == -1:
        print(
            "  [Retry Seed] seed=-1 received; deriving deterministic attempt seeds "
            f"from context='{seed_context_key}'."
        )
        seed_mode = "derived_from_context"
    elif seed > 0:
        print(f"  [Retry Seed] base seed={seed} (deterministic increment per attempt).")
        seed_mode = "provided_base_seed"
    else:
        seed_mode = "provided_literal_seed"

    def _seed_meta() -> Dict[str, Any]:
        return {
            "seed_input": seed,
            "seed_mode": seed_mode,
            "allow_random_seed": allow_random_seed,
            "seed_context": seed_context_key if (seed == -1 and not allow_random_seed) else None,
        }

    templates = None
    attempts: List[RetryAttemptInfo] = []
    output_dir.mkdir(parents=True, exist_ok=True)
    template_fit_diag: Optional[Dict[str, Any]] = None
    reduction_blocked_by_template_fit = False
    final_rho_i = rho0
    final_box_length_nm = 0.0

    for attempt_idx in range(max_retries + 1):
        attempt_num = attempt_idx + 1
        attempt_label = f"Attempt {attempt_num}/{max_retries + 1}"

        # === B1: Compute density for this attempt ===
        if allow_density_reduction:
            # rho_i = max(rho0 * (factor ** attempt_index), effective_rho_min)
            rho_unclamped = rho0 * (rho_decrease_factor ** attempt_idx)
            rho_i = max(rho_unclamped, effective_rho_min)
        else:
            rho_unclamped = rho0
            rho_i = rho0
        rho_i = min(rho_i, effective_rho_max)  # Also clamp to max

        # Compute base box size for this density
        box = compute_box_size(composition.total_mass_g, rho_i)
        final_rho_i = rho_i

        # Expand box to accommodate margin (reduce chance of later enlargement)
        if box_margin_nm and box_margin_nm > 0:
            adjusted_length = box.length_nm + 2.0 * box_margin_nm
            volume_nm3 = adjusted_length ** 3
            volume_cm3 = volume_nm3 * 1e-21
            adjusted_rho = composition.total_mass_g / volume_cm3
            box = type(box)(
                length_nm=adjusted_length,
                rho0_g_cm3=adjusted_rho,
                total_mass_g=composition.total_mass_g,
                volume_nm3=volume_nm3,
                requested_rho_g_cm3=getattr(box, "requested_rho_g_cm3", rho_i),
                polymerization_shrinkage_vol_frac=getattr(box, "polymerization_shrinkage_vol_frac", 0.0),
                rho_effective_g_cm3=getattr(box, "rho_effective_g_cm3", rho_i),
                density_candidates_g_cm3=list(getattr(box, "density_candidates_g_cm3", []) or []),
            )
        final_box_length_nm = box.length_nm

        print(
            f"  [{attempt_label}] rho_i={rho_i:.4f} g/cm³ (unclamped={rho_unclamped:.4f}), "
            f"Box={box.length_nm:.3f} nm"
        )

        # === B2: Version attempt inputs AND outputs ===
        attempt_inp_path = output_dir / f"packmol_try{attempt_num}.inp"
        attempt_pdb_name = f"packmol_try{attempt_num}.pdb"
        attempt_pdb_path = output_dir / attempt_pdb_name
        
        # Resolve templates (only once)
        if templates is None:
            templates = resolve_all_templates(
                list(composition.counts.keys()),
                molecules_dir,
                prefer_format,
            )
            diagnostic_margin_nm = edge_margin_nm if edge_margin_nm is not None else (tolerance / 10.0)
            template_fit_diag = _compute_template_fit_diagnostics(
                composition=composition,
                templates=templates,
                box_length_nm=box.length_nm,
                margin_nm=diagnostic_margin_nm,
                pdb_scale_override=pdb_scale_override,
            )
            if template_fit_diag.get("polymer_limiting_detected", False):
                reduction_blocked_by_template_fit = True
                _print_template_fit_table(template_fit_diag)
                msg = (
                    "Template-fit failure: at least one polymer-like template span exceeds the intended "
                    "PACKMOL box (including margin). Provide a compact/pre-relaxed conformer instead of "
                    "relying on density backoff."
                )
                err = PackmolError(msg)
                setattr(err, "retry_attempts", attempts)
                setattr(
                    err,
                    "retry_meta",
                    {
                        "rho0": rho0,
                        "allow_density_reduction": allow_density_reduction,
                        "relative_density_floor": relative_floor,
                        "absolute_density_floor": config_floor if rho_min is not None else None,
                        "cap_density_floor": cap_floor,
                        "effective_density_floor": effective_rho_min,
                        "density_floor_fraction": density_floor_fraction,
                        "max_density_reduction_fraction": reduction_cap_fraction,
                        **_seed_meta(),
                        "template_fit": template_fit_diag,
                        "template_fit_blocked_density_reduction": True,
                    },
                )
                raise err

        if seed > 0:
            attempt_seed = 1 + ((seed - 1 + attempt_idx) % 2_000_000_000)
            attempt_seed_source = "provided_base_seed+attempt_index"
        elif seed == -1:
            if allow_random_seed:
                attempt_seed = -1
                attempt_seed_source = "packmol_internal_random"
            else:
                attempt_seed = _derive_deterministic_packmol_seed(seed_context_key, attempt_idx)
                attempt_seed_source = "derived_from_context+attempt_index"
        else:
            attempt_seed = seed
            attempt_seed_source = "provided_literal_seed"

        # Generate input content with attempt-specific output name
        try:
            content = generate_packmol_input(
                counts=composition.counts,
                templates=templates,
                box_length_nm=box.length_nm,
                output_file=attempt_pdb_path,
                tolerance=tolerance,
                seed=attempt_seed,
                edge_margin_nm=edge_margin_nm,
                placement_hints=placement_hints,
            )
        except ValueError as exc:
            raise PackmolError(f"Invalid PACKMOL box/margin settings: {exc}") from exc
        
        # Write versioned input file
        attempt_inp_path.write_text(content)
        
        # Build attempt info (will update success/error after run)
        attempt_info = RetryAttemptInfo(
            attempt_index=attempt_idx,
            rho_target=rho_i,
            box_length_nm=box.length_nm,
            tolerance=tolerance,
            seed=attempt_seed,
            composition_changed=False,  # We never change composition by default
            success=False,
            inp_file=str(attempt_inp_path),
            pdb_file=str(attempt_pdb_path),
            error_reason=None,
            seed_source=attempt_seed_source,
        )
        print(f"  [{attempt_label}] seed={attempt_seed} ({attempt_seed_source})")
        
        # Run PACKMOL
        try:
            success, stdout, stderr = run_packmol(attempt_inp_path, output_dir)
            
            if success and attempt_pdb_path.exists():
                print(f"  [{attempt_label}] PACKMOL succeeded!")
                attempt_info.success = True
                attempts.append(attempt_info)
                
                # === B2: Copy successful attempt to canonical names ===
                canonical_inp = output_dir / "packmol.inp"
                canonical_pdb = output_dir / "initial.pdb"
                shutil.copy2(attempt_inp_path, canonical_inp)
                shutil.copy2(attempt_pdb_path, canonical_pdb)
                print(f"  [{attempt_label}] Copied to canonical: packmol.inp, initial.pdb")
                
                # Return with canonical PDB path for downstream compatibility
                reduction_used = any(a.rho_target < rho0 - 1e-12 for a in attempts)
                if not allow_density_reduction:
                    reduction_reason = "disabled_by_default"
                elif reduction_used:
                    reduction_reason = "packmol_retry_backoff"
                else:
                    reduction_reason = "not_needed_success_first_attempt"
                max_box_expand_fraction = (
                    (rho0 / effective_rho_min) ** (1.0 / 3.0) - 1.0
                    if (effective_rho_min > 0 and rho0 > 0)
                    else 0.0
                )
                return (
                    canonical_pdb,
                    box,
                    templates,
                    attempt_idx,  # retries_used = attempt_idx (0 means first attempt succeeded)
                    attempts,
                    {
                        "rho0": rho0,
                        "allow_density_reduction": allow_density_reduction,
                        "relative_density_floor": relative_floor,
                        "absolute_density_floor": config_floor if rho_min is not None else None,
                        "cap_density_floor": cap_floor,
                        "effective_density_floor": effective_rho_min,
                        "density_floor_fraction": density_floor_fraction,
                        "max_density_reduction_fraction": reduction_cap_fraction,
                        "max_box_expand_fraction_from_density": max_box_expand_fraction,
                        "density_reduction_used": reduction_used,
                        "density_reduction_reason": reduction_reason,
                        **_seed_meta(),
                        "template_fit": template_fit_diag,
                        "template_fit_blocked_density_reduction": reduction_blocked_by_template_fit,
                    },
                )
            else:
                if not attempt_pdb_path.exists():
                    reason = "output_file_missing"
                    print(f"  [{attempt_label}] PACKMOL reported success but output file missing")
                else:
                    _, reason = check_packmol_output(stdout, stderr)
                    print(f"  [{attempt_label}] PACKMOL failed: {reason}")
                
                attempt_info.error_reason = reason
                
                # Save failure log
                log_file = output_dir / f"packmol_try{attempt_num}.log"
                log_file.write_text(
                    "\n".join(
                        [
                            "ATTEMPT_METADATA:",
                            f"attempt_index={attempt_idx}",
                            f"rho_target={rho_i:.8f}",
                            f"box_length_nm={box.length_nm:.8f}",
                            f"tolerance_ang={tolerance:.6f}",
                            f"seed={attempt_seed}",
                            f"seed_source={attempt_seed_source}",
                            f"allow_random_seed={allow_random_seed}",
                            (
                                f"seed_context={seed_context_key}"
                                if seed == -1 and not allow_random_seed
                                else "seed_context=None"
                            ),
                            f"error_reason={reason}",
                            "",
                            "STDOUT:",
                            stdout,
                            "",
                            "STDERR:",
                            stderr,
                            "",
                        ]
                    )
                )
                
        except PackmolError as e:
            attempt_info.error_reason = str(e)
            print(f"  [{attempt_label}] PACKMOL error: {e}")
            log_file = output_dir / f"packmol_try{attempt_num}.log"
            if not log_file.exists():
                log_file.write_text(
                    "\n".join(
                        [
                            "ATTEMPT_METADATA:",
                            f"attempt_index={attempt_idx}",
                            f"rho_target={rho_i:.8f}",
                            f"box_length_nm={box.length_nm:.8f}",
                            f"tolerance_ang={tolerance:.6f}",
                            f"seed={attempt_seed}",
                            f"seed_source={attempt_seed_source}",
                            f"allow_random_seed={allow_random_seed}",
                            (
                                f"seed_context={seed_context_key}"
                                if seed == -1 and not allow_random_seed
                                else "seed_context=None"
                            ),
                            f"error_reason={attempt_info.error_reason}",
                            "",
                            "PACKMOL_EXCEPTION:",
                            str(e),
                            "",
                        ]
                    )
                )

        attempts.append(attempt_info)

        if attempt_idx < max_retries and allow_density_reduction:
            next_rho_unclamped = rho0 * (rho_decrease_factor ** (attempt_idx + 1))
            next_rho = min(max(next_rho_unclamped, effective_rho_min), effective_rho_max)
            if next_rho < rho_i - 1e-12:
                if template_fit_diag:
                    # Recompute against current (pre-reduction) box for a precise decision point.
                    diagnostic_margin_nm = edge_margin_nm if edge_margin_nm is not None else (tolerance / 10.0)
                    template_fit_diag = _compute_template_fit_diagnostics(
                        composition=composition,
                        templates=templates,
                        box_length_nm=box.length_nm,
                        margin_nm=diagnostic_margin_nm,
                        pdb_scale_override=pdb_scale_override,
                    )
                    _print_template_fit_table(template_fit_diag)
                    if template_fit_diag.get("polymer_limiting_detected", False):
                        reduction_blocked_by_template_fit = True
                        msg = (
                            "PACKMOL failure appears limited by oversized polymer-like template span at the "
                            "current target density/box. Provide a compact (coiled/pre-relaxed) polymer "
                            "conformer, or intentionally increase the target box via density settings."
                        )
                        err = PackmolError(msg)
                        setattr(err, "retry_attempts", attempts)
                        setattr(
                            err,
                            "retry_meta",
                            {
                                "rho0": rho0,
                                "allow_density_reduction": allow_density_reduction,
                                "relative_density_floor": relative_floor,
                                "absolute_density_floor": config_floor if rho_min is not None else None,
                                "cap_density_floor": cap_floor,
                                "effective_density_floor": effective_rho_min,
                                "density_floor_fraction": density_floor_fraction,
                                "max_density_reduction_fraction": reduction_cap_fraction,
                                "max_box_expand_fraction_from_density": (
                                    (rho0 / effective_rho_min) ** (1.0 / 3.0) - 1.0
                                    if (effective_rho_min > 0 and rho0 > 0)
                                    else 0.0
                                ),
                                "density_reduction_used": any(a.rho_target < rho0 - 1e-12 for a in attempts),
                                "density_reduction_reason": "template_fit_blocked",
                                **_seed_meta(),
                                "template_fit": template_fit_diag,
                                "template_fit_blocked_density_reduction": True,
                            },
                        )
                        raise err

    if template_fit_diag:
        _print_template_fit_table(template_fit_diag)
    reduction_used = any(a.rho_target < rho0 - 1e-12 for a in attempts)
    if not allow_density_reduction:
        reduction_reason = "disabled_by_default"
    elif reduction_used:
        reduction_reason = "packmol_retry_backoff"
    else:
        reduction_reason = "attempts_exhausted_without_backoff"
    max_box_expand_fraction = (
        (rho0 / effective_rho_min) ** (1.0 / 3.0) - 1.0
        if (effective_rho_min > 0 and rho0 > 0)
        else 0.0
    )
    retry_meta = {
        "rho0": rho0,
        "allow_density_reduction": allow_density_reduction,
        "relative_density_floor": relative_floor,
        "absolute_density_floor": config_floor if rho_min is not None else None,
        "cap_density_floor": cap_floor,
        "effective_density_floor": effective_rho_min,
        "density_floor_fraction": density_floor_fraction,
        "max_density_reduction_fraction": reduction_cap_fraction,
        "max_box_expand_fraction_from_density": max_box_expand_fraction,
        "density_reduction_used": reduction_used,
        "density_reduction_reason": reduction_reason,
        **_seed_meta(),
        "template_fit": template_fit_diag,
        "template_fit_blocked_density_reduction": reduction_blocked_by_template_fit,
    }
    err = PackmolError(
        f"PACKMOL failed after {len(attempts)} attempts. "
        f"rho0={rho0:.4f}, effective_rho_min={effective_rho_min:.4f}, "
        f"final rho={final_rho_i:.4f} g/cm³, Box={final_box_length_nm:.3f} nm. "
        f"Check versioned logs (packmol_try*.log) in {output_dir} for details."
    )
    setattr(err, "retry_attempts", attempts)
    setattr(err, "retry_meta", retry_meta)
    raise err

