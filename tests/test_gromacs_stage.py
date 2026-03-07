import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.stages import gromacs
from pipeline.stages.gromacs import (
    GmxEmStage,
    GmxEqStage,
    GmxProdStage,
    _build_stage_state,
    _get_gromacs_provenance,
    _write_gmx_eq_metrics,
)


class FakeManifest:
    def __init__(self) -> None:
        self.data = {}
        self.commands = []
        self.stage_status = {}

    def set(self, key, value) -> None:
        self.data[key] = value

    def get(self, key, default=None):
        return self.data.get(key, default)

    def log_command(self, command, stage, **kwargs) -> None:
        self.commands.append({"command": command, "stage": stage, **kwargs})

    def set_stage_status(self, stage, status, error=None) -> None:
        self.stage_status[stage] = {"status": status, "error": error}


class FakeContext:
    def __init__(self, root: Path, **overrides) -> None:
        self.root = root
        self.project_root = str(root)
        self.run_id = "RUN123"
        self.system_id = "SYS"
        self.resume = True
        self.force = False
        self.dry_run = False
        self.grompp_maxwarn = 0
        self.allow_velocity_reset = False
        self.allow_gro_velocity_restart = False
        self.allow_resume_repair = False
        self.resume_ignore_runtime = False
        self.strict_mdp_validation = False
        self.strict_restart_policy = None
        self.temperature = 300.0
        self.pressure = None
        self.nsteps_em = None
        self.nsteps_nvt = None
        self.nsteps_npt = None
        self.nsteps_md = None
        self.tc_grps_mode = "system"
        self.comm_grps_policy = "System"
        self.force_gen_vel = False
        self.allow_continuation_override = False
        self.gen_seed = None
        self.npt_barostat = None
        self.gmx_eq_metrics_window_ps = 20.0
        self.manifest = FakeManifest()
        self.gmx_cmd = "gmx"
        self.active_ndx_path = None
        self.input_gro_has_velocities = None
        self.tc_grps_auto_decision = None
        self.cleanup_success = False
        for key, value in overrides.items():
            setattr(self, key, value)

    def get_input_path(self, *parts: str) -> Path:
        return self.root / "IN" / Path(*parts)

    def get_output_path(self, *parts: str) -> Path:
        return self.root / "OUT_GMX" / self.run_id / Path(*parts)

    def ensure_output_dir(self, *parts: str) -> Path:
        path = self.root / "OUT_GMX" / self.run_id / Path(*parts)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def log_command(self, command, stage: str, **kwargs) -> None:
        self.manifest.log_command(command, stage, **kwargs)


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_gro(path: Path, *, with_velocities: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "Test GRO",
        "2",
    ]
    if with_velocities:
        lines.append(
            f"{1:5d}{'MOL':<5}{'A1':>5}{1:5d}"
            f"{0.000:8.3f}{0.000:8.3f}{0.000:8.3f}"
            f"{0.1000:8.4f}{0.2000:8.4f}{0.3000:8.4f}"
        )
        lines.append(
            f"{1:5d}{'MOL':<5}{'A2':>5}{2:5d}"
            f"{0.100:8.3f}{0.100:8.3f}{0.100:8.3f}"
            f"{0.4000:8.4f}{0.5000:8.4f}{0.6000:8.4f}"
        )
    else:
        lines.append(
            f"{1:5d}{'MOL':<5}{'A1':>5}{1:5d}"
            f"{0.000:8.3f}{0.000:8.3f}{0.000:8.3f}"
        )
        lines.append(
            f"{1:5d}{'MOL':<5}{'A2':>5}{2:5d}"
            f"{0.100:8.3f}{0.100:8.3f}{0.100:8.3f}"
        )
    lines.append("   1.00000   1.00000   1.00000")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_topology(path: Path) -> None:
    _write(path, "[ system ]\nTest system\n[ molecules ]\nMOL 1\n")


def _write_mdp(path: Path) -> None:
    _write(
        path,
        "integrator = md\n"
        "nsteps = 1000\n"
        "ref_t = 300\n"
        "gen_temp = 300\n"
        "tc-grps = System\n",
    )


def _write_valid_done_state(stage_dir: Path, stage_name: str, prefix: str, *, extra_outputs=None) -> None:
    extra_outputs = extra_outputs or []
    _write(stage_dir / f"{prefix}.mdp", "stub mdp\n")
    _write(stage_dir / f"{prefix}.tpr", "stub tpr\n")
    _write(stage_dir / f"{prefix}.gro", "stub gro\n")
    for relpath in extra_outputs:
        _write(stage_dir / relpath, "stub output\n")
    _write_json(
        stage_dir / "stage_state.json",
        {"stage": stage_name, "fingerprint_payload": {"stage": stage_name}},
    )
    _write_json(stage_dir / "repro.json", {"stage": stage_name})
    _write(stage_dir / "done.ok", "completed\n")


def _install_common_stage_stubs(monkeypatch):
    monkeypatch.setattr(
        gromacs,
        "_get_git_info",
        lambda project_root: {"commit": "abc1234", "dirty": "false"},
    )
    monkeypatch.setattr(
        gromacs,
        "_get_gromacs_provenance",
        lambda **kwargs: {
            "gmx_command_tokens": ["gmx"],
            "gmx_launcher_token": "gmx",
            "gmx_launcher_resolved": "/usr/bin/gmx",
            "gmx_cmd_used": "gmx",
            "gmx_bin_resolved": "/usr/bin/gmx",
            "gmx_binary_token": "gmx",
            "gmx_binary_resolved": "/usr/bin/gmx",
            "gromacs_version": "GROMACS version: stub",
        },
    )
    monkeypatch.setattr(gromacs, "cleanup_large_files", lambda *args, **kwargs: [])
    monkeypatch.setattr(gromacs, "get_resource_summary", lambda ctx: {"gmx_nt": None})
    monkeypatch.setattr(gromacs, "find_ndx_file", lambda ctx: None)
    monkeypatch.setattr(
        gromacs,
        "_resolve_posres_info",
        lambda ctx, stage, current_structure, mdp_path: (
            None,
            {"active": False, "reference_path": None, "reason": "not_active"},
        ),
    )
    monkeypatch.setattr(
        gromacs,
        "resolve_mdp_overrides",
        lambda ctx, stage, extra_overrides=None: ({}, extra_overrides or {}, extra_overrides or {}),
    )

    def fake_prepare_stage_mdp(ctx, stage, output_dir, effective_overrides=None, checkpoint_available=None):
        mdp_path = output_dir / f"{stage}.mdp"
        _write_mdp(mdp_path)
        return mdp_path, {}, effective_overrides or {}

    monkeypatch.setattr(gromacs, "prepare_stage_mdp", fake_prepare_stage_mdp)

    def fake_build_grompp_command(
        *,
        ctx,
        mdp_path,
        structure_path,
        topology_path,
        output_path,
        checkpoint_path=None,
        index_path=None,
        reference_structure_path=None,
        maxwarn=0,
    ):
        cmd = ["gmx", "grompp", "-f", str(mdp_path), "-c", str(structure_path), "-p", str(topology_path), "-o", str(output_path)]
        if checkpoint_path:
            cmd.extend(["-t", str(checkpoint_path)])
        if index_path:
            cmd.extend(["-n", str(index_path)])
        if reference_structure_path:
            cmd.extend(["-r", str(reference_structure_path)])
        return cmd

    monkeypatch.setattr(gromacs, "build_grompp_command", fake_build_grompp_command)

    def fake_build_mdrun_command(ctx, stage, extra_args=None):
        deffnm = stage
        cmd = ["gmx", "mdrun", "-deffnm", deffnm]
        if extra_args:
            cmd.extend(extra_args)
        return cmd

    monkeypatch.setattr(gromacs, "build_mdrun_command", fake_build_mdrun_command)

    calls = []

    def fake_run_gromacs_command(ctx, cmd, workdir, stage_name, stream_logs=False, log_prefix=None):
        tokens = [str(token) for token in cmd]
        calls.append(tokens)
        workdir = Path(workdir)
        if "grompp" in tokens:
            output_path = Path(tokens[tokens.index("-o") + 1])
            _write(output_path, "compiled tpr\n")
            return
        deffnm = tokens[tokens.index("-deffnm") + 1] if "-deffnm" in tokens else workdir.name
        _write(workdir / f"{deffnm}.gro", "gro output\n")
        _write(workdir / f"{deffnm}.cpt", "checkpoint\n")
        _write(workdir / f"{deffnm}.xtc", "trajectory\n")
        _write(workdir / f"{deffnm}.log", "mdrun log\n")
        logs_dir = workdir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        if log_prefix:
            _write(logs_dir / f"{log_prefix}.stdout.log", "streamed stdout\n")
        else:
            _write(logs_dir / "mdrun.stdout.log", "streamed stdout\n")

    monkeypatch.setattr(gromacs, "run_gromacs_command", fake_run_gromacs_command)
    return calls


def _build_md_state_for_resume(ctx: FakeContext, monkeypatch, *, velocity_mode=None):
    _install_common_stage_stubs(monkeypatch)
    output_dir = ctx.ensure_output_dir("03_gromacs", "md")
    npt_dir = ctx.ensure_output_dir("03_gromacs", "npt")
    top_path = ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")
    _write_topology(top_path)
    npt_gro = npt_dir / "npt.gro"
    _write_gro(npt_gro)
    mdp_path = output_dir / "md.mdp"
    _write_mdp(mdp_path)
    _write(output_dir / "md.tpr", "md tpr\n")
    _write(output_dir / "md.cpt", "md checkpoint\n")
    stage_state = _build_stage_state(
        ctx,
        "md",
        mdp_path,
        npt_gro,
        top_path,
        None,
        "gmx grompp -f md.mdp -o md.tpr",
        ["gmx", "mdrun", "-deffnm", "md"],
        ctx.grompp_maxwarn,
        posres_info={"active": False, "reference_path": None, "reason": "not_active"},
        stage_overrides={},
        extra_overrides={},
        effective_overrides={},
        velocity_mode=velocity_mode,
    )
    _write_json(output_dir / "stage_state.json", stage_state)
    return output_dir, npt_dir, top_path, npt_gro


def test_gmx_em_run_writes_done_ok_after_success(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path)
    calls = _install_common_stage_stubs(monkeypatch)
    gro_path = ctx.get_input_path("systems", ctx.system_id, "gromacs", "gro", "system.gro")
    top_path = ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")
    _write_gro(gro_path)
    _write_topology(top_path)
    monkeypatch.setattr(gromacs, "verify_staged_inputs", lambda ctx: (gro_path, top_path))

    stage = GmxEmStage()
    assert stage.run(ctx) is True

    output_dir = ctx.get_output_path("03_gromacs", "em")
    assert calls
    assert (output_dir / "stage_state.json").exists()
    assert (output_dir / "repro.json").exists()
    assert (output_dir / "done.ok").exists()


def test_gmx_em_skips_only_when_done_ok_exists(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path)
    output_dir = ctx.ensure_output_dir("03_gromacs", "em")
    _write_valid_done_state(output_dir, "em", "em")

    monkeypatch.setattr(gromacs, "verify_staged_inputs", lambda ctx: pytest.fail("should skip before input checks"))

    stage = GmxEmStage()
    assert stage.run(ctx) is True


def test_gmx_prod_fresh_run_writes_done_ok_after_success(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path)
    _install_common_stage_stubs(monkeypatch)

    npt_dir = ctx.ensure_output_dir("03_gromacs", "npt")
    npt_gro = npt_dir / "npt.gro"
    _write_gro(npt_gro)
    _write(npt_dir / "npt.cpt", "npt checkpoint\n")

    top_path = ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")
    _write_topology(top_path)
    monkeypatch.setattr(gromacs, "verify_staged_inputs", lambda ctx: (npt_gro, top_path))

    stage = GmxProdStage()
    assert stage.run(ctx) is True

    output_dir = ctx.get_output_path("03_gromacs", "md")
    assert (output_dir / "stage_state.json").exists()
    assert (output_dir / "repro.json").exists()
    assert (output_dir / "done.ok").exists()


def test_gmx_prod_resume_run_writes_done_ok_when_resume_matches(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path)
    _build_md_state_for_resume(ctx, monkeypatch, velocity_mode=None)
    monkeypatch.setattr(gromacs, "verify_staged_inputs", lambda ctx: (None, ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")))

    stage = GmxProdStage()
    assert stage.run(ctx) is True

    output_dir = ctx.get_output_path("03_gromacs", "md")
    assert (output_dir / "done.ok").exists()


def test_gmx_prod_same_stage_resume_keeps_velocity_mode_in_fingerprint(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path)
    output_dir, _, _, _ = _build_md_state_for_resume(ctx, monkeypatch, velocity_mode="checkpoint")
    monkeypatch.setattr(gromacs, "verify_staged_inputs", lambda ctx: (None, ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")))

    stage = GmxProdStage()
    assert stage.run(ctx) is True

    state = json.loads((output_dir / "stage_state.json").read_text(encoding="utf-8"))
    assert state["velocity_mode"] == "checkpoint"
    assert state["fingerprint_payload"]["velocity_mode"] == "checkpoint"


def test_gmx_eq_npt_mutually_exclusive_restart_flags_fail_cleanly(tmp_path, monkeypatch, capsys):
    ctx = FakeContext(
        tmp_path,
        allow_velocity_reset=True,
        allow_gro_velocity_restart=True,
    )
    nvt_dir = ctx.ensure_output_dir("03_gromacs", "nvt")
    _write_gro(nvt_dir / "nvt.gro")

    stage = GmxEqStage()
    assert stage._run_npt(ctx) is False

    captured = capsys.readouterr()
    assert "mutually exclusive" in captured.out


def test_gmx_prod_mutually_exclusive_restart_flags_fail_cleanly(tmp_path, monkeypatch, capsys):
    ctx = FakeContext(
        tmp_path,
        allow_velocity_reset=True,
        allow_gro_velocity_restart=True,
    )
    npt_dir = ctx.ensure_output_dir("03_gromacs", "npt")
    _write_gro(npt_dir / "npt.gro")

    top_path = ctx.get_input_path("systems", ctx.system_id, "gromacs", "top", "system.top")
    _write_topology(top_path)
    monkeypatch.setattr(gromacs, "verify_staged_inputs", lambda ctx: (None, top_path))

    stage = GmxProdStage()
    assert stage.run(ctx) is False

    captured = capsys.readouterr()
    assert "mutually exclusive" in captured.out


def test_get_gromacs_provenance_separates_launcher_from_binary(monkeypatch):
    which_map = {
        "mpirun": "/usr/bin/mpirun",
        "gmx_mpi": "/opt/gromacs/bin/gmx_mpi",
    }
    monkeypatch.setattr(gromacs.shutil, "which", lambda token: which_map.get(token))
    monkeypatch.setattr(gromacs, "_get_gromacs_version", lambda token: f"GROMACS version via {token}")

    provenance = _get_gromacs_provenance(
        ctx=None,
        grompp_cmd="mpirun -np 4 gmx_mpi grompp -f md.mdp -o md.tpr",
        mdrun_cmd_base=None,
    )

    assert provenance["gmx_command_tokens"][:4] == ["mpirun", "-np", "4", "gmx_mpi"]
    assert provenance["gmx_launcher_token"] == "mpirun"
    assert provenance["gmx_launcher_resolved"] == "/usr/bin/mpirun"
    assert provenance["gmx_cmd_used"] == "gmx_mpi"
    assert provenance["gmx_bin_resolved"] == "/opt/gromacs/bin/gmx_mpi"
    assert provenance["gmx_binary_token"] == "gmx_mpi"
    assert provenance["gmx_binary_resolved"] == "/opt/gromacs/bin/gmx_mpi"
    assert provenance["gromacs_version"] == "GROMACS version via gmx_mpi"


def test_write_gmx_eq_metrics_defaults_to_20ps_and_marks_diagnostic_window(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path, gmx_eq_metrics_window_ps=0.0)
    npt_dir = ctx.ensure_output_dir("03_gromacs", "npt")
    captured = {}

    def fake_build_metrics(ctx_arg, npt_dir_arg, window_ps):
        captured["window_ps"] = window_ps
        return {
            "stage": "gmx_eq",
            "temperature_K": {"avg": 300.0, "std": 1.0},
            "pressure_bar": {"avg": 1.0, "std": 0.1},
            "density_kg_m3": {"avg": 1000.0, "std": 10.0},
            "density_g_cm3": {"avg": 1.0, "std": 0.01},
            "volume_nm3": {"avg": 1.0, "std": 0.01},
            "last_window_ps": 20.0,
            "notes": [],
        }

    monkeypatch.setattr(gromacs, "_build_gmx_eq_metrics", fake_build_metrics)

    alias_path = _write_gmx_eq_metrics(ctx, npt_dir)
    primary_metrics = json.loads((npt_dir / "stage_metrics.json").read_text(encoding="utf-8"))
    alias_metrics = json.loads(alias_path.read_text(encoding="utf-8"))

    assert captured["window_ps"] == 20.0
    assert primary_metrics["window_ps_configured"] == 20.0
    assert primary_metrics["window_semantics"] == "diagnostic_tail_window"
    assert alias_metrics == primary_metrics


def test_write_gmx_eq_metrics_preserves_actual_window_with_configured_request(tmp_path, monkeypatch):
    ctx = FakeContext(tmp_path, gmx_eq_metrics_window_ps=75.0)
    npt_dir = ctx.ensure_output_dir("03_gromacs", "npt")
    captured = {}

    def fake_build_metrics(ctx_arg, npt_dir_arg, window_ps):
        captured["window_ps"] = window_ps
        return {
            "stage": "gmx_eq",
            "temperature_K": {"avg": 300.0, "std": 1.0},
            "pressure_bar": {"avg": 1.0, "std": 0.1},
            "density_kg_m3": {"avg": 1000.0, "std": 10.0},
            "density_g_cm3": {"avg": 1.0, "std": 0.01},
            "volume_nm3": {"avg": 1.0, "std": 0.01},
            "last_window_ps": 12.5,
            "notes": [],
        }

    monkeypatch.setattr(gromacs, "_build_gmx_eq_metrics", fake_build_metrics)

    _write_gmx_eq_metrics(ctx, npt_dir)
    metrics = json.loads((npt_dir / "stage_metrics.json").read_text(encoding="utf-8"))

    assert captured["window_ps"] == 75.0
    assert metrics["window_ps_configured"] == 75.0
    assert metrics["last_window_ps"] == 12.5
    assert metrics["window_semantics"] == "diagnostic_tail_window"
