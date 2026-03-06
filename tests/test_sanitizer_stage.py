from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.itp_sanitizer import AtomtypeEntry, DefaultsEntry, ItpSanitizer, generate_system_top
from pipeline.stages.sanitizer_stage import SanitizerError, SanitizerStage


class FakeManifest:
    def __init__(self) -> None:
        self.data = {}
        self.sanitizer_outputs = {}

    def set_sanitizer_output(self, key, value) -> None:
        self.sanitizer_outputs[key] = value

    def set(self, key, value) -> None:
        self.data[key] = value

    def get(self, key, default=None):
        return self.data.get(key, default)

    def save(self) -> None:
        return None


class FakeContext:
    def __init__(self, root: Path, **overrides) -> None:
        self.root = root
        self.system_id = "SYS"
        self.ff = "OPLS-AA"
        self.run_id = "RUN123"
        self.manifest = FakeManifest()
        self.allow_default_defaults = True
        self.strict_include_resolution = False
        self.allow_include_shadowing = True
        self.allow_unsafe_include_escape = False
        self.itp_include_dirs_priority = "forcefield_first"
        self.itp_include_dirs = ""
        self.project_root = str(root)
        self.allow_override = False
        self.atomtype_prefix = None
        self.strict_charge = False
        self.strict_forcefield_consistency = False
        self.charge_fix_target_allowlist = None
        self.charge_fix_target_atomnames = None
        self.charge_fix_disallowed_atomtypes = None
        self.strict_charge_neutrality = False
        self.charge_fix_max_delta_per_atom = None
        self.charge_fix_max_total = None
        self.charge_fix_max_dipole_shift_debye = None
        self.charge_neutrality_correct = 1e-4
        self.strict_charge_physics = False
        self.charge_fix_method = "safe_subset"
        self.charge_fix_allow_ions = False
        self.charge_fix_allow_solvents = False
        self.charge_fix_protect_resnames = None
        self.charge_neutrality_warn = 1e-5
        self.charge_neutrality_tol = 1e-6
        self.polymer_net_charge_tol = 1e-3
        self.charge_fix_polymer_method = "skip_if_small"
        self.charge_fix_polymer_exclusion_bonds = 2
        self.charge_fix_moleculetype_rounding_tol = 1e-4
        self.charge_fix_allow_non_rounding = False
        self.allow_unparsed_atoms_lines = False
        self.strict_dipole_check = False
        self.strict_gro_top_check = False
        self.grompp_maxwarn = 0
        self.strict_top_update = False
        self.allow_mixed_defaults = False
        self.allow_mixed_defaults_reason = None
        self.mixed_defaults_cross_group_policy = None
        self.mixed_defaults_cross_group_rule = None
        self.mixed_defaults_cross_group_max_pairs = None
        self.mixed_defaults_cross_group_reason = None
        self.mixed_defaults_preserve_within_group = False
        self.gmx_cmd = "gmx"
        for key, value in overrides.items():
            setattr(self, key, value)

    def get_input_path(self, *parts) -> Path:
        return self.root / "IN" / Path(*parts)

    def ensure_output_dir(self, subdir: str) -> Path:
        path = self.root / "OUT_GMX" / self.run_id / subdir
        path.mkdir(parents=True, exist_ok=True)
        return path


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _prepare_base_tree(tmp_path: Path, top_content: str) -> FakeContext:
    ctx = FakeContext(tmp_path)
    _write(
        tmp_path / "IN/forcefield/oplsaa/gromacs/forcefield.itp",
        "[ atomtypes ]\n; none\n",
    )
    _write(
        tmp_path / "IN/systems/SYS/gromacs/top/system.top",
        top_content,
    )
    (tmp_path / "IN/systems/SYS/gromacs/itp").mkdir(parents=True, exist_ok=True)
    (tmp_path / "IN/systems/SYS/htpolynet/itp").mkdir(parents=True, exist_ok=True)
    (tmp_path / "IN/molecules").mkdir(parents=True, exist_ok=True)
    return ctx


def _make_stub_result(system_gromacs_dir: Path, run_id: str, mixed_defaults=False):
    combined_versioned = system_gromacs_dir / f"combined_atomtypes_{run_id}.itp"
    combined_current = system_gromacs_dir / "combined_atomtypes_current.itp"
    sanitized_dir = system_gromacs_dir / f"itp_sanitized_{run_id}"
    sanitized_current = system_gromacs_dir / "itp_sanitized_current"
    sanitized_dir.mkdir(parents=True, exist_ok=True)
    sanitized_current.mkdir(parents=True, exist_ok=True)
    _write(combined_versioned, "[ atomtypes ]\n")
    _write(combined_current, "[ atomtypes ]\n")
    dummy_itp = sanitized_current / "dummy.itp"
    _write(dummy_itp, "[ moleculetype ]\nDUMMY 3\n")
    return SimpleNamespace(
        atomtype_count=1,
        sanitized_files=[dummy_itp],
        molecule_itp_files=[dummy_itp],
        conflicts_overridden=[],
        charge_warnings=[],
        defaults_used=SimpleNamespace(
            nbfunc=1,
            comb_rule=2,
            gen_pairs="yes",
            fudge_lj=0.5,
            fudge_qq=0.8333,
            source_file=str(combined_current),
        ),
        mixed_defaults_detected=mixed_defaults,
        mixed_defaults_report={
            "classification": "comb-rule-only",
            "chosen_policy": "unsafe_override",
            "primary_defaults_signature": "sig-primary",
            "differing_fields": ["comb_rule"],
            "secondary_signatures": ["sig-secondary"],
        }
        if mixed_defaults
        else {},
        nonbond_params_secondary_current_path=None,
        nonbond_params_cross_group_current_path=None,
        nonbond_params_secondary_path=None,
        nonbond_params_cross_group_path=None,
        nonbond_params_secondary_summary={},
        nonbond_params_cross_group_summary={},
        prefix_policy_summary={},
        prefix_injected_types_current_path=None,
        prefix_injected_types_path=None,
        conflicts_detected=[],
        combined_atomtypes_current_path=combined_current,
        combined_atomtypes_path=combined_versioned,
        sanitized_dir=sanitized_dir,
        sanitized_current_dir=sanitized_current,
        source_defaults_map={},
        defaults_guessed=False,
    )


def _patch_sanitizer_run(monkeypatch, ctx: FakeContext, result=None):
    from pipeline import itp_sanitizer

    class StubItpSanitizer:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def run(self, **kwargs):
            system_gromacs_dir = ctx.get_input_path("systems", ctx.system_id, "gromacs")
            return result or _make_stub_result(system_gromacs_dir, ctx.run_id)

    monkeypatch.setattr(itp_sanitizer, "ItpSanitizer", StubItpSanitizer)


def test_typing_imports_present_for_optional_and_set():
    content = Path("pipeline/stages/sanitizer_stage.py").read_text(encoding="utf-8")
    assert "Optional" in content
    assert "Set" in content


def test_serialize_manifest_entry_handles_non_dataclass_objects(tmp_path):
    stage = SanitizerStage()

    class Dummy:
        def __init__(self, path: Path):
            self.path = path
            self.values = {"x", "y"}

    payload = stage._serialize_manifest_entry(Dummy(tmp_path / "a.itp"))
    assert payload["path"] == str(tmp_path / "a.itp")
    assert sorted(payload["values"]) == ["x", "y"]


def test_uncertain_molecule_counts_fail_without_trusted_preprocessed_source(tmp_path):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\n#ifdef TEST\nA 1\n#endif\n",
    )
    stage = SanitizerStage()
    with pytest.raises(SanitizerError, match="Cannot derive trusted \\[molecules\\] counts"):
        stage.run(ctx)


def test_existing_system_top_reports_updated_existing_when_regenerated(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\n#ifdef TEST\nA 1\n#endif\n",
    )
    _write(tmp_path / "IN/systems/SYS/gromacs/gro/system.gro", "stub\n0\n0 0 0\n")
    _patch_sanitizer_run(monkeypatch, ctx)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_expected_atoms_from_grompp", lambda **kwargs: 0)
    monkeypatch.setattr(stage, "_parse_preprocessed_molecules", lambda path: [("A", 0)])
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    assert stage.run(ctx) is True
    status = ctx.manifest.sanitizer_outputs["system_top_update"]["status"]
    truth_source = ctx.manifest.sanitizer_outputs["topology_truth_source"]
    assert status["mode"] == "updated_existing"
    assert status["preserved_existing_molecules"] is False
    assert truth_source["source"] == "grompp_preprocessed_topology"


def test_stage_prefers_existing_preprocessed_topology_truth_source(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 1\n",
    )
    _patch_sanitizer_run(monkeypatch, ctx)

    stage = SanitizerStage()
    pp_path = stage._grompp_preprocessed_top_path(ctx.ensure_output_dir(stage.output_subdir))
    _write(pp_path, "[ molecules ]\nA 0\n")
    monkeypatch.setattr(stage, "_expected_atoms_from_grompp", lambda **kwargs: (_ for _ in ()).throw(AssertionError("grompp -pp should not run when preferred preprocessed topology already exists")))
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    assert stage.run(ctx) is True
    truth_source = ctx.manifest.sanitizer_outputs["topology_truth_source"]
    count_source = ctx.manifest.sanitizer_outputs["molecule_count_source"]
    assert truth_source["source"] == "preprocessed_topology"
    assert truth_source["is_preprocessed"] is True
    assert count_source["source"] == "preferred_preprocessed_topology"


def test_stage_records_python_fallback_truth_source_when_preprocessed_unavailable(
    tmp_path,
    monkeypatch,
):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 0\n",
    )
    _patch_sanitizer_run(monkeypatch, ctx)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    assert stage.run(ctx) is True
    truth_source = ctx.manifest.sanitizer_outputs["topology_truth_source"]
    count_source = ctx.manifest.sanitizer_outputs["molecule_count_source"]
    assert truth_source["source"] == "python_fallback_raw_topology"
    assert truth_source["fallback_used"] is True
    assert count_source["topology_truth_source"]["source"] == "python_fallback_raw_topology"


def test_active_molecules_require_sanitized_mapping(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 1\n",
    )
    _write(
        tmp_path / "IN/systems/SYS/gromacs/itp/a.itp",
        "[ moleculetype ]\nA 3\n[ atoms ]\n1 C 1 A C1 1 0.0 12.0\n",
    )
    _patch_sanitizer_run(monkeypatch, ctx)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)
    monkeypatch.setattr(stage, "_map_sanitized_itps_to_molecules", lambda *args, **kwargs: {})

    with pytest.raises(SanitizerError, match="mapping is empty"):
        stage.run(ctx)


def test_stage_does_not_use_forcefield_only_include_prepass(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 0\n",
    )
    _patch_sanitizer_run(monkeypatch, ctx)

    from pipeline import itp_sanitizer

    def _unexpected(*args, **kwargs):
        raise AssertionError("ItpParser.resolve_includes should not be called by sanitizer_stage")

    monkeypatch.setattr(itp_sanitizer.ItpParser, "resolve_includes", _unexpected)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    assert stage.run(ctx) is True


def test_path_identity_normalizes_equivalent_paths(tmp_path):
    stage = SanitizerStage()
    p = tmp_path / "x" / "y.itp"
    _write(p, "")
    alt = tmp_path / "x" / ".." / "x" / "y.itp"
    assert stage._path_identity(p) == stage._path_identity(alt)


def test_garbled_symbols_removed_from_stage_messages():
    content = Path("pipeline/stages/sanitizer_stage.py").read_text(encoding="utf-8")
    assert "螖" not in content
    assert "|Δq|" in content
    assert "|ΔQ|" in content
    assert "|Δμ|" in content


def test_mixed_defaults_generate_requires_cross_group_artifact(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 0\n",
    )
    ctx.allow_mixed_defaults = True
    ctx.allow_mixed_defaults_reason = "test"
    result = _make_stub_result(
        ctx.get_input_path("systems", ctx.system_id, "gromacs"),
        ctx.run_id,
        mixed_defaults=True,
    )
    result.nonbond_params_cross_group_summary = {"policy": "generate", "status": "generated"}
    result.nonbond_params_cross_group_current_path = (
        ctx.get_input_path("systems", ctx.system_id, "gromacs")
        / "itp_sanitized_current"
        / "missing_cross_group.itp"
    )
    _patch_sanitizer_run(monkeypatch, ctx, result=result)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    with pytest.raises(SanitizerError, match="cross-group remediation include path is missing"):
        stage.run(ctx)


def test_mixed_defaults_generate_includes_cross_group_artifact_in_system_top(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 0\n",
    )
    ctx.allow_mixed_defaults = True
    ctx.allow_mixed_defaults_reason = "test"
    result = _make_stub_result(
        ctx.get_input_path("systems", ctx.system_id, "gromacs"),
        ctx.run_id,
        mixed_defaults=True,
    )
    cross_group = (
        ctx.get_input_path("systems", ctx.system_id, "gromacs")
        / "itp_sanitized_current"
        / "nonbond_params_cross_group_current.itp"
    )
    _write(cross_group, "[ nonbond_params ]\n")
    result.nonbond_params_cross_group_summary = {"policy": "generate", "status": "generated"}
    result.nonbond_params_cross_group_current_path = cross_group
    _patch_sanitizer_run(monkeypatch, ctx, result=result)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    assert stage.run(ctx) is True
    system_top = ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")
    content = system_top.read_text(encoding="utf-8")
    assert 'nonbond_params_cross_group_current.itp' in content


def test_stage_orders_sanitizer_generated_extra_includes_deterministically(tmp_path, monkeypatch):
    ctx = _prepare_base_tree(
        tmp_path,
        "[ molecules ]\nA 0\n",
    )
    ctx.allow_mixed_defaults = True
    ctx.allow_mixed_defaults_reason = "test"
    result = _make_stub_result(
        ctx.get_input_path("systems", ctx.system_id, "gromacs"),
        ctx.run_id,
        mixed_defaults=True,
    )
    prefix_include = (
        ctx.get_input_path("systems", ctx.system_id, "gromacs")
        / "itp_sanitized_current"
        / "prefix_injected_types_current.itp"
    )
    cross_group = (
        ctx.get_input_path("systems", ctx.system_id, "gromacs")
        / "itp_sanitized_current"
        / "nonbond_params_cross_group_current.itp"
    )
    secondary_include = (
        ctx.get_input_path("systems", ctx.system_id, "gromacs")
        / "itp_sanitized_current"
        / "nonbond_params_secondary_current.itp"
    )
    _write(prefix_include, "[ bondtypes ]\n")
    _write(cross_group, "[ nonbond_params ]\n")
    _write(secondary_include, "[ nonbond_params ]\n")
    result.prefix_injected_types_current_path = prefix_include
    result.nonbond_params_cross_group_current_path = cross_group
    result.nonbond_params_cross_group_summary = {"policy": "generate", "status": "generated"}
    result.nonbond_params_secondary_current_path = secondary_include
    _patch_sanitizer_run(monkeypatch, ctx, result=result)

    stage = SanitizerStage()
    monkeypatch.setattr(stage, "_warn_atomtype_outliers", lambda *args, **kwargs: [])
    monkeypatch.setattr(stage, "_warn_pairs_consistency", lambda *args, **kwargs: None)

    assert stage.run(ctx) is True
    system_top = ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")
    content = system_top.read_text(encoding="utf-8")
    assert content.count('prefix_injected_types_current.itp') == 1
    assert content.count('nonbond_params_cross_group_current.itp') == 1
    assert content.count('nonbond_params_secondary_current.itp') == 1
    assert content.index('prefix_injected_types_current.itp') < content.index('nonbond_params_cross_group_current.itp')
    assert content.index('nonbond_params_cross_group_current.itp') < content.index('nonbond_params_secondary_current.itp')


def test_generate_cross_group_nonbond_params_produces_current_artifact(tmp_path):
    ctx = FakeContext(
        tmp_path,
        allow_mixed_defaults=True,
        allow_mixed_defaults_reason="test",
        mixed_defaults_cross_group_policy="generate",
        mixed_defaults_cross_group_reason="test",
    )
    sanitizer = ItpSanitizer(source_ff="OPLS-AA")
    primary_source = tmp_path / "primary.itp"
    secondary_source = tmp_path / "secondary.itp"
    primary_defaults = DefaultsEntry(
        nbfunc=1,
        comb_rule=2,
        gen_pairs="yes",
        fudge_lj=0.5,
        fudge_qq=0.8333,
        source_file=str(primary_source),
    )
    secondary_defaults = DefaultsEntry(
        nbfunc=1,
        comb_rule=3,
        gen_pairs="yes",
        fudge_lj=0.5,
        fudge_qq=0.8333,
        source_file=str(secondary_source),
    )
    sanitizer.canonical_defaults = primary_defaults
    sanitizer.registry = {
        "AA": AtomtypeEntry(
            name="AA",
            mass=12.0,
            charge=0.0,
            ptype="A",
            sigma=0.30,
            epsilon=0.20,
            source_ff="OPLS-AA",
            source_files=[str(primary_source)],
        ),
        "BB": AtomtypeEntry(
            name="BB",
            mass=12.0,
            charge=0.0,
            ptype="A",
            sigma=0.40,
            epsilon=0.10,
            source_ff="GAFF2",
            source_files=[str(secondary_source)],
        ),
    }

    versioned_path, current_path, summary = sanitizer._generate_cross_group_nonbond_params(
        output_dir=tmp_path,
        run_id="RUN123",
        source_defaults_map={
            primary_source: primary_defaults,
            secondary_source: secondary_defaults,
        },
        defaults_report={
            "mixed_defaults_detected": True,
            "classification": "comb-rule-only",
            "primary_defaults_signature": primary_defaults.signature(),
            "secondary_signatures": [secondary_defaults.signature()],
        },
        scan_paths=[],
        ctx=ctx,
    )

    assert summary["status"] == "generated"
    assert summary["pair_count"] == 1
    assert versioned_path == tmp_path / "nonbond_params_cross_group_RUN123.itp"
    assert current_path == tmp_path / "nonbond_params_cross_group_current.itp"
    assert versioned_path.exists()
    assert current_path.exists()
    content = current_path.read_text(encoding="utf-8")
    assert "[ nonbond_params ]" in content
    assert "AA" in content
    assert "BB" in content
    assert not list(tmp_path.glob(".nonbond_params_cross_group_current.itp.tmp.*"))


def test_generate_system_top_dedupes_and_orders_extra_and_molecule_includes(tmp_path):
    top_dir = tmp_path / "gromacs/top"
    combined_atomtypes = top_dir / "combined_atomtypes_current.itp"
    molecule_a = tmp_path / "gromacs/itp/A.itp"
    molecule_b = tmp_path / "gromacs/itp/B.itp"
    _write(combined_atomtypes, "[ atomtypes ]\n")
    _write(molecule_a, "[ moleculetype ]\nA 3\n")
    _write(molecule_b, "[ moleculetype ]\nB 3\n")

    content = generate_system_top(
        combined_atomtypes_path=combined_atomtypes,
        sanitized_itp_paths=[molecule_a, molecule_b, molecule_a],
        molecule_itp_files=[molecule_b, molecule_a, molecule_a],
        system_name="SYS",
        molecule_counts={"B": 1, "A": 2},
        defaults=DefaultsEntry(
            nbfunc=1,
            comb_rule=2,
            gen_pairs="yes",
            fudge_lj=0.5,
            fudge_qq=0.8333,
            source_file=str(combined_atomtypes),
        ),
        top_dir=top_dir,
        extra_includes=[
            "prefix_injected_types_current.itp",
            "nonbond_params_cross_group_current.itp",
            "prefix_injected_types_current.itp",
            "nonbond_params_secondary_current.itp",
        ],
    )

    assert content.count('prefix_injected_types_current.itp') == 1
    assert content.count('nonbond_params_cross_group_current.itp') == 1
    assert content.count('nonbond_params_secondary_current.itp') == 1
    assert content.count('#include "../itp/B.itp"') == 1
    assert content.count('#include "../itp/A.itp"') == 1
    assert content.index('combined_atomtypes_current.itp') < content.index('prefix_injected_types_current.itp')
    assert content.index('prefix_injected_types_current.itp') < content.index('nonbond_params_cross_group_current.itp')
    assert content.index('nonbond_params_cross_group_current.itp') < content.index('nonbond_params_secondary_current.itp')
    assert content.index('nonbond_params_secondary_current.itp') < content.index('../itp/B.itp')
    assert content.index('../itp/B.itp') < content.index('../itp/A.itp')
