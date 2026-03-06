from pathlib import Path
from types import SimpleNamespace
import math
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.stages.topology_sanitizer import (
    ATOMS_ROW_STATUS_INVALID,
    ATOMS_ROW_STATUS_PARSED,
    ATOMS_ROW_STATUS_UNSUPPORTED,
    DEFAULT_MOLECULE_SIG_QUANT_STEP,
    SanitizerError,
    TOP_MOLECULES_ENTRY_STATUS_PARSED,
    TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED,
    TopMoleculesParseResult,
    TopologySanitizerMixin,
    parse_atoms_row,
    parse_top_molecules_entry,
)


class DummySanitizer(TopologySanitizerMixin):
    def __init__(self) -> None:
        self.diagnostics = []
        self._forcefield_dir = None
        self._sanitized_current_dir = None


class DummyContext:
    def __init__(self, root: Path, **overrides) -> None:
        self.root = root
        self.system_id = "SYS"
        self.project_root = str(root)
        self.itp_include_dirs = ""
        self.itp_include_dirs_priority = "forcefield_first"
        self.strict_include_resolution = False
        self.allow_include_shadowing = True
        self.allow_unsafe_include_escape = False
        self.strict_forcefield_consistency = False
        self.strict_charge_neutrality = False
        self.charge_neutrality_tol = 1e-6
        self.charge_fix_protect_resnames = None
        self.charge_fix_target_allowlist = None
        self.charge_fix_allow_ions = False
        self.charge_fix_allow_solvents = False
        self.charge_fix_max_delta_per_atom = 1e-5
        self.charge_fix_polymer_method = "skip_if_small"
        self.polymer_net_charge_tol = 1e-3
        self.manifest = None
        self.lj_outlier_policy = "warn"
        self.strict_lj_validation = False
        self.lj_bounds_profile = "aa"
        self.lj_outlier_thresholds = None
        for key, value in overrides.items():
            setattr(self, key, value)

    def get_input_path(self, *parts) -> Path:
        return self.root / "IN" / Path(*parts)


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _simple_itp(name: str, charge: float = 0.0) -> str:
    return (
        "[ moleculetype ]\n"
        f"{name} 3\n"
        "[ atoms ]\n"
        f"1 CT 1 {name} C1 1 {charge:.6f} 12.011\n"
        f"2 HC 1 {name} H1 1 0.000000 1.008\n"
    )


def test_recursive_include_resolution_rebuilds_current_file_search_stack(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(
        tmp_path,
        itp_include_dirs="user_include",
        itp_include_dirs_priority="sanitized_first",
    )

    forcefield_dir = tmp_path / "IN/forcefield/oplsaa"
    sanitizer._forcefield_dir = forcefield_dir

    root_itp = forcefield_dir / "gromacs/root.itp"
    child_itp = tmp_path / "IN/systems/SYS/gromacs/itp/child.itp"
    user_leaf = tmp_path / "user_include/leaf.itp"
    ff_leaf = forcefield_dir / "gromacs/leaf.itp"

    _write(root_itp, '#include "child.itp"\n[ moleculetype ]\nROOT 3\n')
    _write(child_itp, '#include "leaf.itp"\n')
    _write(user_leaf, "[ moleculetype ]\nUSER_LEAF 3\n")
    _write(ff_leaf, "[ moleculetype ]\nFF_LEAF 3\n")

    names = sanitizer._get_moleculetype_names(root_itp, ctx=ctx)

    assert "USER_LEAF" in names
    assert "FF_LEAF" not in names


def test_uncertain_same_name_conflict_fails_closed_in_strict_mode(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_charge_neutrality=True)
    itp1 = tmp_path / "a.itp"
    itp2 = tmp_path / "b.itp"
    content = (
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "#ifdef SOMEFLAG\n"
        "#endif\n"
        "1 CT 1 MOL C1 1 0.000000 12.011\n"
    )
    _write(itp1, content)
    _write(itp2, content)

    with pytest.raises(SanitizerError, match="Uncertain same-name moleculetype comparison"):
        sanitizer._moleculetype_signatures_differ(itp1, itp2, "MOL", ctx=ctx)


def test_uncertain_same_name_conflict_is_not_treated_as_equal_in_non_strict_mode(
    tmp_path,
    capsys,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    itp1 = tmp_path / "a.itp"
    itp2 = tmp_path / "b.itp"
    content = (
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "#ifdef SOMEFLAG\n"
        "#endif\n"
        "1 CT 1 MOL C1 1 0.000000 12.011\n"
    )
    _write(itp1, content)
    _write(itp2, content)

    assert sanitizer._moleculetype_signatures_differ(itp1, itp2, "MOL", ctx=ctx) is True
    assert "Treating definitions as conflicting" in capsys.readouterr().out


def test_charge_fix_classification_strict_fails_on_active_parse_uncertainty(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_charge_neutrality=True)
    broken = tmp_path / "broken.itp"
    _write(
        broken,
        "[ moleculetype ]\n"
        "BROKEN 3\n"
        "[ atoms ]\n"
        "#ifdef SOMEFLAG\n"
        "1 CT 1 BROKEN C1 1 0.000000 12.011\n"
        "#endif\n",
    )

    with pytest.raises(SanitizerError, match="Charge-fix classification for active molecule 'BROKEN' is degraded"):
        sanitizer._classify_ionic_moleculetypes({"BROKEN": broken}, ctx=ctx)


def test_degraded_and_protected_polymer_classification_refuses_unallowlisted_targets(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    broken = tmp_path / "broken.itp"
    polymer = tmp_path / "polymer.itp"
    _write(
        broken,
        "[ moleculetype ]\n"
        "BROKEN 3\n"
        "[ atoms ]\n"
        "#ifdef SOMEFLAG\n"
        "1 CT 1 BROKEN C1 1 0.000000 12.011\n"
        "#endif\n",
    )
    _write(
        polymer,
        "[ moleculetype ]\n"
        "POLYMER_CHAIN 3\n"
        "[ atoms ]\n"
        "#ifdef SOMEFLAG\n"
        "1 CT 1 POLYMER_CHAIN C1 1 0.000000 12.011\n"
        "#endif\n",
    )

    classifications, audit = sanitizer._classify_ionic_moleculetypes(
        {"BROKEN": broken, "POLYMER_CHAIN": polymer},
        ctx=ctx,
    )

    assert classifications["BROKEN"] == "unknown"
    assert classifications["POLYMER_CHAIN"] == "protected_polymer"
    assert audit["degraded_molecules"]["BROKEN"]["status"] == "degraded"
    assert audit["degraded_molecules"]["POLYMER_CHAIN"]["classification"] == "protected_polymer"

    correction_result = SimpleNamespace(
        target_molecule="BROKEN",
        corrected=True,
        correction_refused=False,
        max_delta_applied=1e-6,
        total_delta=1e-6,
        correction_method="safe_subset",
        atoms_adjusted=1,
        effective_method="safe_subset",
        neutrality_status="corrected",
        polymer_policy=None,
    )
    protection_audit = sanitizer._verify_ion_protection(
        correction_result,
        classifications,
        ctx,
        classification_audit=audit,
    )
    assert any("BROKEN" in item and "classification" in item for item in protection_audit["violations"])

    allow_ctx = DummyContext(tmp_path, charge_fix_target_allowlist="BROKEN")
    allow_audit = sanitizer._verify_ion_protection(
        correction_result,
        classifications,
        allow_ctx,
        classification_audit=audit,
    )
    assert allow_audit["violations"] == []
    assert allow_audit["target_from_allowlist"] is True


def test_charge_fix_protect_resnames_extends_built_in_protection(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, charge_fix_protect_resnames="G4_ALT")
    itp = tmp_path / "g4.itp"
    _write(itp, _simple_itp("G4_ALT"))

    classifications, audit = sanitizer._classify_ionic_moleculetypes({"G4_ALT": itp}, ctx=ctx)

    assert classifications["G4_ALT"] == "protected_configured"
    assert audit["details"]["G4_ALT"]["source"] == "charge_fix_protect_resnames"


def test_moleculetype_signature_quantization_is_deterministic_at_boundaries():
    sanitizer = DummySanitizer()
    step = DEFAULT_MOLECULE_SIG_QUANT_STEP
    half = step / 2.0
    below = math.nextafter(half, 0.0)
    above = math.nextafter(half, 1.0)

    assert sanitizer._quantize_moleculetype_charge(half, step) == 0.0
    assert sanitizer._quantize_moleculetype_charge(below, step) == 0.0
    assert sanitizer._quantize_moleculetype_charge(above, step) == step
    assert sanitizer._quantize_moleculetype_charge(-above, step) == -step


def test_parse_atoms_row_parses_standard_8_column_record():
    result = parse_atoms_row("1 CT 2 MOL C1 7 -0.125000 12.011")

    assert result.status == ATOMS_ROW_STATUS_PARSED
    assert result.record is not None
    assert result.record.atom_id == 1
    assert result.record.atom_type == "CT"
    assert result.record.resnr == 2
    assert result.record.residue == "MOL"
    assert result.record.atom_name == "C1"
    assert result.record.cgnr == 7
    assert result.record.charge == pytest.approx(-0.125)
    assert result.record.mass == pytest.approx(12.011)


@pytest.mark.parametrize(
    ("line", "expected_cgnr", "expected_charge"),
    [
        ("1 HC 1 MOL H1 0.125000 1.008", None, 0.125),
        ("1 HC 1 MOL H1 3-0.125000 1.008", 3, -0.125),
    ],
)
def test_parse_atoms_row_preserves_supported_7_column_variants(
    line,
    expected_cgnr,
    expected_charge,
):
    result = parse_atoms_row(line)

    assert result.status == ATOMS_ROW_STATUS_PARSED
    assert result.record is not None
    assert result.record.cgnr == expected_cgnr
    assert result.record.charge == pytest.approx(expected_charge)
    assert result.record.mass == pytest.approx(1.008)


def test_parse_atoms_row_returns_explicit_failure_states():
    invalid = parse_atoms_row("1 CT 1 MOL C1")
    unsupported = parse_atoms_row("1 CT 1 MOL C1 1 0.0 12.0 2 0.0 12.0")

    assert invalid.status == ATOMS_ROW_STATUS_INVALID
    assert invalid.record is None
    assert unsupported.status == ATOMS_ROW_STATUS_UNSUPPORTED
    assert unsupported.record is None


def test_parse_atoms_row_rejects_ambiguous_raw_7_column_integer_token():
    result = parse_atoms_row("1 CT 1 MOL C1 1 -0.125000")

    assert result.status == ATOMS_ROW_STATUS_UNSUPPORTED
    assert result.record is None


def test_parse_atoms_row_rejects_raw_6_column_glued_charge_mass_form():
    result = parse_atoms_row("1 CT 1 MOL C1 0.125000-12.011")

    assert result.status == ATOMS_ROW_STATUS_INVALID
    assert result.record is None


def test_parse_top_molecules_entry_parses_integer_count():
    result = parse_top_molecules_entry("SOL 42")

    assert result.status == TOP_MOLECULES_ENTRY_STATUS_PARSED
    assert result.entry is not None
    assert result.entry.name == "SOL"
    assert result.entry.count_token == "42"
    assert result.entry.count_value == 42


def test_parse_top_molecules_entry_surfaces_unresolved_macro_count():
    result = parse_top_molecules_entry("POLY N_POLY")

    assert result.status == TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED
    assert result.entry is not None
    assert result.entry.name == "POLY"
    assert result.entry.count_token == "N_POLY"
    assert result.entry.count_value is None


def test_extract_moleculetype_atom_charges_remains_compatible_with_stage1_ir(tmp_path):
    sanitizer = DummySanitizer()
    itp = tmp_path / "mol.itp"
    _write(
        itp,
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "1 CT 1 MOL C1 1 -0.125000 12.011\n"
        "2 HC 1 MOL H1 0.125000 1.008\n",
    )

    charges, uncertain = sanitizer._extract_moleculetype_atom_charges(itp, "MOL")

    assert charges == pytest.approx([-0.125, 0.125])
    assert uncertain is False


def test_get_ordered_molecules_from_top_preserves_ordered_counts_and_raw_entries(tmp_path):
    sanitizer = DummySanitizer()
    top = tmp_path / "system.top"
    _write(
        top,
        "[ molecules ]\n"
        "A 1\n"
        "B N_B\n"
        "A 2\n",
    )

    result = sanitizer._get_ordered_molecules_from_top(top)

    assert result.ordered == [("A", 3)]
    assert [(entry.name, entry.count_token, entry.count_value) for entry in result.entries] == [
        ("A", "1", 1),
        ("B", "N_B", None),
        ("A", "2", 2),
    ]
    assert result.uncertain is True
    assert any("N_B" in reason for reason in result.uncertainty_reasons)


def test_top_molecules_parse_result_positional_init_remains_compatible():
    result = TopMoleculesParseResult([("A", 1)], True, ["raw-uncertain"])

    assert result.ordered == [("A", 1)]
    assert result.uncertain is True
    assert result.uncertainty_reasons == ["raw-uncertain"]
    assert result.entries == []


def test_negative_sigma_policy_warns_by_default_and_errors_when_requested(tmp_path):
    sanitizer = DummySanitizer()
    combined = tmp_path / "combined_atomtypes.itp"
    _write(
        combined,
        "[ atomtypes ]\n"
        "NEG 12.011 0.000 A -1.000000e-01 2.500000e-01\n",
    )

    warn_ctx = DummyContext(tmp_path)
    warnings = sanitizer._warn_atomtype_outliers(combined, warn_ctx, comb_rule=2)
    assert any("Negative sigma" in warning for warning in warnings)

    error_ctx = DummyContext(tmp_path, lj_outlier_policy="error")
    with pytest.raises(SanitizerError, match="LJ outlier policy=error triggered"):
        sanitizer._warn_atomtype_outliers(combined, error_ctx, comb_rule=2)
