"""
Unit tests for decide_prev_state_policy restart behavior.

Tests acceptance criteria:
A) prev_cpt_path exists → checkpoint restart
B) prev_cpt_path missing + both flags False + strict → error
C) prev_cpt_path missing + allow_velocity_reset=True → velocity_reset mode
D) prev_cpt_path missing + allow_gro_velocity_restart=True + GRO has velocities → allowed
E) prev_cpt_path missing + both flags True → ValueError (mutual exclusivity)
F) prev_cpt_path missing + allow_gro_velocity_restart=True + GRO lacks velocities → error
G) prev_cpt_path missing + no flags + non-strict requires GRO velocity validation
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.stages.gromacs import decide_prev_state_policy, _gro_has_velocities


class MockContext:
    """Minimal mock of PipelineContext for testing."""
    def __init__(self, **kwargs):
        self.temperature = 300.0
        self.allow_velocity_reset = False
        self.allow_gro_velocity_restart = False
        self.strict_mdp_validation = False
        for k, v in kwargs.items():
            setattr(self, k, v)


def create_gro_file(with_velocities: bool) -> Path:
    """Create a temporary GRO file with or without velocities."""
    fd, path = tempfile.mkstemp(suffix=".gro")
    with open(path, "w") as f:
        f.write("Test GRO file\n")
        f.write("2\n")  # 2 atoms
        if with_velocities:
            f.write(
                f"{1:5d}{'MOL':<5}{'A1':>5}{1:5d}"
                f"{0.000:8.3f}{0.000:8.3f}{0.000:8.3f}"
                f"{0.1000:8.4f}{0.2000:8.4f}{0.3000:8.4f}\n"
            )
            f.write(
                f"{1:5d}{'MOL':<5}{'A2':>5}{2:5d}"
                f"{0.100:8.3f}{0.100:8.3f}{0.100:8.3f}"
                f"{0.1000:8.4f}{0.2000:8.4f}{0.3000:8.4f}\n"
            )
        else:
            f.write(
                f"{1:5d}{'MOL':<5}{'A1':>5}{1:5d}"
                f"{0.000:8.3f}{0.000:8.3f}{0.000:8.3f}\n"
            )
            f.write(
                f"{1:5d}{'MOL':<5}{'A2':>5}{2:5d}"
                f"{0.100:8.3f}{0.100:8.3f}{0.100:8.3f}\n"
            )
        f.write("   1.00000   1.00000   1.00000\n")  # Box
    return Path(path)


def create_cpt_file() -> Path:
    """Create a dummy checkpoint file."""
    fd, path = tempfile.mkstemp(suffix=".cpt")
    with open(path, "wb") as f:
        f.write(b"dummy checkpoint")
    return Path(path)


def create_mdp_file() -> Path:
    """Create a minimal MDP template."""
    fd, path = tempfile.mkstemp(suffix=".mdp")
    with open(path, "w") as f:
        f.write("integrator = md\n")
        f.write("nsteps = 1000\n")
        f.write("ref_t = 300\n")
        f.write("gen_temp = 300\n")
    return Path(path)


class TestGROVelocityDetection:
    """Test _gro_has_velocities helper."""
    
    def test_gro_with_velocities(self):
        """GRO file with velocities should return True."""
        gro = create_gro_file(with_velocities=True)
        try:
            assert _gro_has_velocities(gro) == True
        finally:
            gro.unlink()
    
    def test_gro_without_velocities(self):
        """GRO file without velocities should return False."""
        gro = create_gro_file(with_velocities=False)
        try:
            assert _gro_has_velocities(gro) == False
        finally:
            gro.unlink()
    
    def test_nonexistent_file(self):
        """Non-existent file should return None (unreadable)."""
        assert _gro_has_velocities(Path("/nonexistent/file.gro")) is None


class TestRestartPolicy:
    """Test decide_prev_state_policy behavior."""
    
    def test_a_checkpoint_exists(self):
        """A) Checkpoint exists → checkpoint restart."""
        cpt = create_cpt_file()
        try:
            cpt_available, _, velocity_mode, _, _ = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=True,
                allow_velocity_reset=False,
                allow_gro_velocity_restart=False,
            )
            assert cpt_available == True
            assert velocity_mode == "checkpoint"
        finally:
            cpt.unlink()
    
    def test_b_no_cpt_both_false_strict(self):
        """B) No checkpoint + both flags False + strict → error."""
        cpt = Path("/nonexistent.cpt")
        cpt_available, _, velocity_mode, warnings, _ = decide_prev_state_policy(
            stage_name="test",
            prev_cpt_path=cpt,
            strict_mode=True,
            allow_velocity_reset=False,
            allow_gro_velocity_restart=False,
        )
        assert cpt_available == False
        assert velocity_mode == "error"
        assert len(warnings) > 0
    
    def test_c_no_cpt_velocity_reset(self):
        """C) No checkpoint + allow_velocity_reset=True → velocity_reset mode."""
        cpt = Path("/nonexistent.cpt")
        mdp = create_mdp_file()
        ctx = MockContext(temperature=300.0)
        try:
            cpt_available, extra, velocity_mode, _, _ = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=True,
                allow_velocity_reset=True,
                allow_gro_velocity_restart=False,
                template_mdp_path=mdp,
                ctx=ctx,
            )
            assert cpt_available == False
            assert velocity_mode == "velocity_reset"
            assert extra is not None
            assert extra.get("gen_vel") == "yes"
        finally:
            mdp.unlink()
    
    def test_d_no_cpt_gro_restart_with_velocities(self):
        """D) No checkpoint + allow_gro_velocity_restart=True + GRO has velocities → allowed."""
        cpt = Path("/nonexistent.cpt")
        gro = create_gro_file(with_velocities=True)
        try:
            cpt_available, extra, velocity_mode, _, meta = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=True,
                allow_velocity_reset=False,
                allow_gro_velocity_restart=True,
                prev_gro_path=gro,
            )
            assert cpt_available == False
            assert velocity_mode == "gro_velocity_restart"
            assert extra is not None
            assert extra.get("gen_vel") == "no"
            assert meta["gro_velocity_validation_performed"] is True
            assert meta["gro_velocity_validation_outcome"] == "has_velocities"
        finally:
            gro.unlink()
    
    def test_e_mutual_exclusivity(self):
        """E) Both flags True → ValueError."""
        cpt = Path("/nonexistent.cpt")
        try:
            decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=True,
                allow_velocity_reset=True,
                allow_gro_velocity_restart=True,
            )
            assert False, "Expected ValueError"
        except ValueError as e:
            assert "mutually exclusive" in str(e)
    
    def test_f_gro_restart_without_velocities(self):
        """F) No checkpoint + allow_gro_velocity_restart=True + GRO lacks velocities → error."""
        cpt = Path("/nonexistent.cpt")
        gro = create_gro_file(with_velocities=False)
        try:
            cpt_available, _, velocity_mode, warnings, meta = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=True,
                allow_velocity_reset=False,
                allow_gro_velocity_restart=True,
                prev_gro_path=gro,
            )
            assert cpt_available == False
            assert velocity_mode == "error"
            assert len(warnings) > 0
            assert "lacks velocities" in warnings[0]
            assert meta["gro_velocity_validation_outcome"] == "missing_velocities"
        finally:
            gro.unlink()

    def test_g_no_cpt_non_strict_no_flags_requires_validation(self):
        """G) Non-strict/no-flags path must validate .gro velocities."""
        cpt = Path("/nonexistent.cpt")
        gro_missing_vel = create_gro_file(with_velocities=False)
        gro_with_vel = create_gro_file(with_velocities=True)
        try:
            cpt_available, _, velocity_mode, warnings, meta = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=False,
                allow_velocity_reset=False,
                allow_gro_velocity_restart=False,
                prev_gro_path=gro_missing_vel,
            )
            assert cpt_available is False
            assert velocity_mode == "error"
            assert "lacks velocities" in warnings[0]
            assert meta["gro_velocity_validation_performed"] is True
            assert meta["gro_velocity_validation_outcome"] == "missing_velocities"

            cpt_available, _, velocity_mode, warnings, meta = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=False,
                allow_velocity_reset=False,
                allow_gro_velocity_restart=False,
                prev_gro_path=gro_with_vel,
            )
            assert cpt_available is False
            assert velocity_mode == "gro_velocity_restart"
            assert "auto-fallback" in warnings[0]
            assert meta["gro_velocity_validation_performed"] is True
            assert meta["gro_velocity_validation_outcome"] == "has_velocities"

            cpt_available, _, velocity_mode, warnings, _ = decide_prev_state_policy(
                stage_name="test",
                prev_cpt_path=cpt,
                strict_mode=False,
                allow_velocity_reset=False,
                allow_gro_velocity_restart=False,
                prev_gro_path=None,
            )
            assert cpt_available is False
            assert velocity_mode == "error"
            assert "no prev_gro_path" in warnings[0].lower()
        finally:
            gro_missing_vel.unlink()
            gro_with_vel.unlink()


def run_tests():
    """Run all tests."""
    passed = 0
    failed = 0
    
    print("=== GRO Velocity Detection Tests ===")
    t = TestGROVelocityDetection()
    for name in ["test_gro_with_velocities", "test_gro_without_velocities", "test_nonexistent_file"]:
        try:
            getattr(t, name)()
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print("\n=== Restart Policy Tests ===")
    t2 = TestRestartPolicy()
    for name in ["test_a_checkpoint_exists", "test_b_no_cpt_both_false_strict",
                 "test_c_no_cpt_velocity_reset", "test_d_no_cpt_gro_restart_with_velocities",
                 "test_e_mutual_exclusivity", "test_f_gro_restart_without_velocities",
                 "test_g_no_cpt_non_strict_no_flags_requires_validation"]:
        try:
            getattr(t2, name)()
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
        except Exception as e:
            print(f"  FAIL: {name} - {type(e).__name__}: {e}")
            failed += 1
    
    print(f"\n=== Results: {passed} passed, {failed} failed ===")
    return failed == 0


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
