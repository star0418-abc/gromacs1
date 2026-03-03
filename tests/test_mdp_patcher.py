"""
Unit tests for MDP Patcher validation logic.

Tests P0/P1 safety improvements:
- P0.1: Stage policy for continuation/gen_vel
- P0.2: Expert restart guard with checkpoint-context guard
- P0.3: gen_temp validation when gen_vel=yes
- P0.B: Continuation input requirement messaging
- P1.5: Validation coverage consistency
- P1.7: Numeric parse error handling
"""
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.mdp_patcher import (
    _validate_gen_vel_continuation_consistency,
    _validate_gen_vel_consistency,
    _validate_vector_lengths,
    copy_and_patch_mdp,
    get_mdp_overrides_for_stage,
    parse_ndx_group_sizes,
    parse_ndx_groups,
    validate_mdp_consistency,
    MDPPatchError,
)


# =============================================================================
# Mock PipelineContext for testing
# =============================================================================
class MockContext:
    """Minimal mock of PipelineContext for testing."""
    def __init__(self, **kwargs):
        # Set defaults
        self.temperature = None
        self.pressure = None
        self.nsteps_em = None
        self.nsteps_nvt = None
        self.nsteps_npt = None
        self.nsteps_md = None
        self.tc_grps_mode = "system"
        self.comm_grps_policy = "match_tc"
        self.force_gen_vel = False
        self.allow_continuation_override = False
        self.allow_gro_velocity_restart = False
        self.input_gro_has_velocities = None
        self.active_ndx_path = None
        self.strict_mdp_validation = False
        self.comm_mode = "Linear"
        self.gen_seed = None
        self.npt_barostat = "Berendsen"
        self.allow_comm_mode_angular_with_pbc = False
        self.tc_grps_min_atoms = 50
        self.allow_small_tc_grps = False
        self.system_type = None
        # Override with provided values
        for k, v in kwargs.items():
            setattr(self, k, v)


# =============================================================================
# P0.1: Stage policy for continuation/gen_vel
# =============================================================================
class TestP01StagePolicyLoop:
    """Test that stage policy closes the loop for both continuation and gen_vel."""
    
    def test_nvt_defaults_fresh_start(self):
        """NVT should produce continuation=no, gen_vel=yes."""
        ctx = MockContext()
        overrides = get_mdp_overrides_for_stage(ctx, "nvt")
        assert overrides.get("continuation") == "no"
        assert overrides.get("gen_vel") == "yes"
    
    def test_npt_defaults_continuation(self):
        """NPT should produce continuation=yes, gen_vel=no."""
        ctx = MockContext()
        overrides = get_mdp_overrides_for_stage(ctx, "npt")
        assert overrides.get("continuation") == "yes"
        assert overrides.get("gen_vel") == "no"
    
    def test_md_defaults_continuation(self):
        """MD (production) should produce continuation=yes, gen_vel=no."""
        ctx = MockContext()
        overrides = get_mdp_overrides_for_stage(ctx, "md")
        assert overrides.get("continuation") == "yes"
        assert overrides.get("gen_vel") == "no"

    def test_allow_gro_velocity_restart_prefers_continuation_no(self):
        """Expert gro restart should default to continuation=no in NPT/MD stage policy."""
        ctx = MockContext(allow_gro_velocity_restart=True)
        npt_overrides = get_mdp_overrides_for_stage(ctx, "npt")
        md_overrides = get_mdp_overrides_for_stage(ctx, "md")
        assert npt_overrides.get("continuation") == "no"
        assert npt_overrides.get("gen_vel") == "no"
        assert md_overrides.get("continuation") == "no"
        assert md_overrides.get("gen_vel") == "no"
    
    def test_force_gen_vel_overrides_npt(self):
        """force_gen_vel=True should override gen_vel for NPT and set continuation=no."""
        ctx = MockContext(force_gen_vel=True)
        overrides = get_mdp_overrides_for_stage(ctx, "npt")
        assert overrides.get("gen_vel") == "yes"
        # force_gen_vel MUST set continuation=no to avoid contradiction (gen_vel=yes + continuation=yes is invalid)
        assert overrides.get("continuation") == "no"

    def test_force_gen_vel_propagates_gen_seed_for_npt_md(self):
        """When force_gen_vel is active, gen_seed should be set for NPT/MD as well."""
        ctx = MockContext(force_gen_vel=True, gen_seed=4242)
        npt_overrides = get_mdp_overrides_for_stage(ctx, "npt")
        md_overrides = get_mdp_overrides_for_stage(ctx, "md")
        assert npt_overrides.get("gen_vel") == "yes"
        assert md_overrides.get("gen_vel") == "yes"
        assert npt_overrides.get("gen_seed") == 4242
        assert md_overrides.get("gen_seed") == 4242
    
    def test_allow_continuation_override_skips_continuation(self):
        """allow_continuation_override=True should not set continuation."""
        ctx = MockContext(allow_continuation_override=True)
        overrides = get_mdp_overrides_for_stage(ctx, "npt")
        assert "continuation" not in overrides
        # gen_vel should still be set
        assert overrides.get("gen_vel") == "no"


# =============================================================================
# P0.2: Expert restart guard
# =============================================================================
class TestP02ExpertRestartGuard:
    """Test expert restart guard with allow_gro_velocity_restart and checkpoint-context."""
    
    def test_gen_vel_no_continuation_no_warning_nonstrict(self):
        """gen_vel=no, continuation=no without flag should produce WARNING in non-strict."""
        params = {"gen_vel": "no", "continuation": "no"}
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "npt", strict=False, ctx=None, checkpoint_available=False
        )
        assert len(errors) == 0, f"Expected no errors, got: {errors}"
        assert len(warnings) == 1
        assert "RISK" in warnings[0]
        assert "zero velocities" in warnings[0].lower()
    
    def test_gen_vel_no_continuation_no_error_strict(self):
        """gen_vel=no, continuation=no without flag should produce ERROR in strict mode."""
        params = {"gen_vel": "no", "continuation": "no"}
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "npt", strict=True, ctx=None, checkpoint_available=False
        )
        assert len(errors) == 1
        assert "RISK" in errors[0]
    
    def test_expert_restart_allowed_with_flag(self):
        """With allow_gro_velocity_restart and known velocities, expert restart is allowed."""
        ctx = MockContext(allow_gro_velocity_restart=True, input_gro_has_velocities=True)
        params = {"gen_vel": "no", "continuation": "no"}
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "npt", strict=True, ctx=ctx, checkpoint_available=False
        )
        # Should be allowed (logged as acknowledgment)
        assert len(errors) == 0, f"Expected no errors, got: {errors}"
        assert any("expert restart acknowledged" in w.lower() for w in warnings)

    def test_expert_restart_errors_when_gro_velocities_missing(self):
        """allow_gro_velocity_restart must hard-fail when .gro velocities are known missing."""
        ctx = MockContext(allow_gro_velocity_restart=True, input_gro_has_velocities=False)
        params = {"gen_vel": "no", "continuation": "no"}
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "md", strict=False, ctx=ctx, checkpoint_available=False
        )
        assert len(errors) == 1
        assert "hotspots/constraint failures" in errors[0]
        assert "gen_vel=yes with gen_temp" in errors[0]
    
    def test_expert_restart_blocked_when_checkpoint_available(self):
        """Even with allow_gro_velocity_restart, ERROR when checkpoint is available."""
        ctx = MockContext(allow_gro_velocity_restart=True)
        params = {"gen_vel": "no", "continuation": "no"}
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "npt", strict=True, ctx=ctx, checkpoint_available=True
        )
        # Should still be error because checkpoint exists
        assert len(errors) == 1
        assert "checkpoint" in errors[0].lower()
    
    def test_contradictory_combo_always_error(self):
        """gen_vel=yes + continuation=yes should ALWAYS be an error."""
        params = {"gen_vel": "yes", "continuation": "yes"}
        # Even in non-strict mode
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "nvt", strict=False, ctx=None, checkpoint_available=False
        )
        assert len(errors) == 1
        assert "contradictory" in errors[0].lower()


# =============================================================================
# P0.3: gen_temp validation when gen_vel=yes
# =============================================================================
class TestP03GenTempValidation:
    """Test gen_temp presence and numeric validation."""
    
    def test_gen_temp_missing_warning_nonstrict(self):
        """gen_vel=yes with missing gen_temp should warn in non-strict mode."""
        params = {"gen_vel": "yes", "ref_t": "300"}
        errors, warnings = _validate_gen_vel_consistency(params, "nvt", strict=False)
        # Should have warning about missing gen_temp
        assert any("gen_temp is not set" in w for w in warnings)
    
    def test_gen_temp_missing_error_strict(self):
        """gen_vel=yes with missing gen_temp should error in strict mode."""
        params = {"gen_vel": "yes", "ref_t": "300"}
        errors, warnings = _validate_gen_vel_consistency(params, "nvt", strict=True)
        assert any("gen_temp is not set" in e for e in errors)
    
    def test_gen_temp_non_numeric_caught(self):
        """gen_temp='300K' (non-numeric) should be caught."""
        params = {"gen_vel": "yes", "gen_temp": "300K", "ref_t": "300"}
        errors, warnings = _validate_gen_vel_consistency(params, "nvt", strict=True)
        assert any("not a valid number" in e for e in errors)
    
    def test_ref_t_non_numeric_caught(self):
        """ref_t='300K' (non-numeric) should be caught."""
        params = {"gen_vel": "yes", "gen_temp": "300", "ref_t": "300K"}
        errors, warnings = _validate_gen_vel_consistency(params, "nvt", strict=True)
        assert any("non-numeric" in e for e in errors)


# =============================================================================
# P0.B: Continuation input requirement messaging
# =============================================================================
class TestP0BContinuationMessaging:
    """Test that continuation=yes produces INFO message about input requirements."""
    
    def test_continuation_yes_produces_info_message(self):
        """continuation=yes should produce INFO about input requirements."""
        params = {"gen_vel": "no", "continuation": "yes"}
        errors, warnings = _validate_gen_vel_continuation_consistency(
            params, "npt", strict=False, ctx=None, checkpoint_available=False
        )
        # Should have INFO message
        assert any("continuation=yes" in w and "Checkpoint file" in w for w in warnings)


# =============================================================================
# Task A: Semantic key override collision fix
# =============================================================================
class TestSemanticKeyOverrideCollision:
    """Test that overrides update existing variant instead of creating duplicates."""
    
    def test_override_underscore_when_template_has_hyphen(self, tmp_path):
        """Override gen_temp should update gen-temp if template uses hyphen."""
        from pipeline.mdp_patcher import copy_and_patch_mdp
        
        # Create template with gen-temp (hyphenated)
        template = tmp_path / "test.mdp"
        template.write_text("gen-temp = 0\ngen_vel = yes\n")
        output = tmp_path / "output.mdp"
        
        # Override with underscore variant
        orig, patched = copy_and_patch_mdp(
            template, output, {"gen_temp": 300}, validate=False
        )
        
        # Read output and verify
        content = output.read_text()
        assert "gen-temp = 300" in content, f"Expected gen-temp=300, got: {content}"
        # Should NOT have both variants
        lines = [l.strip() for l in content.split("\n") if "gen" in l.lower() and "temp" in l.lower()]
        assert len([l for l in lines if l.startswith("gen") and "temp" in l]) == 1, \
            f"Should have exactly one gen*temp line, got: {lines}"
    
    def test_override_hyphen_when_template_has_underscore(self, tmp_path):
        """Override gen-temp should update gen_temp if template uses underscore."""
        from pipeline.mdp_patcher import copy_and_patch_mdp
        
        template = tmp_path / "test.mdp"
        template.write_text("gen_temp = 0\ngen_vel = yes\n")
        output = tmp_path / "output.mdp"
        
        orig, patched = copy_and_patch_mdp(
            template, output, {"gen-temp": 300}, validate=False
        )
        
        content = output.read_text()
        assert "gen_temp = 300" in content, f"Expected gen_temp=300, got: {content}"


# =============================================================================
# Task C: Temperature injection for missing keys
# =============================================================================
class TestTemperatureInjection:
    """Test that temperature overrides inject missing ref_t/gen_temp."""
    
    def test_inject_ref_t_when_missing(self, tmp_path):
        """Temperature override should inject ref_t when template lacks it."""
        from pipeline.mdp_patcher import copy_and_patch_mdp
        
        template = tmp_path / "test.mdp"
        template.write_text("nsteps = 1000\n")
        output = tmp_path / "output.mdp"
        
        orig, patched = copy_and_patch_mdp(
            template, output, {"temperature": 300}, validate=False
        )
        
        content = output.read_text()
        # Should have ref-t or ref_t = 300
        assert "ref-t = 300" in content or "ref_t = 300" in content, \
            f"Expected ref_t/ref-t=300 to be injected, got: {content}"
    
    def test_inject_gen_temp_when_gen_vel_yes(self, tmp_path):
        """When gen_vel=yes and gen_temp missing, temperature override injects it."""
        from pipeline.mdp_patcher import copy_and_patch_mdp
        
        template = tmp_path / "test.mdp"
        template.write_text("nsteps = 1000\ngen_vel = yes\n")
        output = tmp_path / "output.mdp"
        
        orig, patched = copy_and_patch_mdp(
            template, output, {"temperature": 300}, validate=False
        )
        
        content = output.read_text()
        # Should have gen-temp or gen_temp = 300
        assert "gen-temp = 300" in content or "gen_temp = 300" in content, \
            f"Expected gen_temp/gen-temp=300 to be injected, got: {content}"


class TestSplitAndCommPolicy:
    """Sanity checks for override outcome status and comm-grps policy hardening."""

    def test_split_skip_in_non_strict_is_warning_not_error(self, tmp_path):
        template = tmp_path / "test.mdp"
        template.write_text(
            "tc-grps = A B\n"
            "tau_t = 0.1 0.2\n"
            "ref_t = 300 300\n"
            "tcoupl = v-rescale\n"
        )
        output = tmp_path / "output.mdp"
        ctx = MockContext(tc_grps_groups="A B C", strict_mdp_validation=False)
        orig, patched, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"tc_grps_mode": "split"},
            validate=True,
            ctx=ctx,
            stage="nvt",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("split skipped" in w.lower() for w in warns)
        content = output.read_text()
        assert "tc-grps = A B" in content

    def test_comm_policy_auto_chooses_system_for_polymer_risk(self, tmp_path):
        template = tmp_path / "test.mdp"
        template.write_text(
            "tc-grps = Polymer Solvent\n"
            "comm-grps = Polymer Solvent\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "output.mdp"
        ctx = MockContext(strict_mdp_validation=False)
        orig, patched, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"comm_grps_policy": "auto"},
            validate=True,
            ctx=ctx,
            stage="md",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("using system-wide com removal" in w.lower() for w in warns)
        content = output.read_text()
        assert "comm-grps = System" in content

    def test_comm_policy_match_tc_warns_strongly_for_polymer_risk(self, tmp_path):
        template = tmp_path / "test.mdp"
        template.write_text(
            "tc-grps = Polymer Solvent\n"
            "comm-grps = System\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "output.mdp"
        ndx = tmp_path / "index.ndx"
        ndx.write_text("[ Polymer ]\n1\n[ Solvent ]\n2\n")
        ctx = MockContext(active_ndx_path=str(ndx), strict_mdp_validation=False)
        orig, patched, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"comm_grps_policy": "match_tc"},
            validate=True,
            ctx=ctx,
            stage="md",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("strong warning" in w.lower() for w in warns)

    def test_comm_policy_missing_ndx_group_strict_errors(self, tmp_path):
        template = tmp_path / "test.mdp"
        template.write_text(
            "tc-grps = GroupA GroupB\n"
            "comm-grps = System\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "output.mdp"
        ndx = tmp_path / "index.ndx"
        ndx.write_text("[ GroupA ]\n1\n")
        ctx = MockContext(active_ndx_path=str(ndx), strict_mdp_validation=True)
        try:
            copy_and_patch_mdp(
                template,
                output,
                {"comm_grps_policy": "match_tc"},
                validate=True,
                ctx=ctx,
                stage="npt",
            )
        except MDPPatchError as e:
            assert "missing in index file" in str(e)
        else:
            raise AssertionError("Expected MDPPatchError for missing strict comm-grps ndx group")


class TestValidateConsistencyCtxPassThrough:
    """Ensure final validation honors context flags and .gro velocity detection."""

    def test_validate_uses_ctx_for_velocity_restart_guard(self, tmp_path):
        mdp_path = tmp_path / "test.mdp"
        mdp_path.write_text("gen_vel = no\ncontinuation = no\n")
        ctx = MockContext(allow_gro_velocity_restart=True, input_gro_has_velocities=False)
        errors, warns = validate_mdp_consistency(
            mdp_path,
            expected_values={},
            stage="npt",
            strict=False,
            ctx=ctx,
            checkpoint_available=False,
        )
        assert len(errors) == 1
        assert "hotspots/constraint failures" in errors[0]


class TestCompressibilitySafety:
    def test_pr_gel_injects_small_compressibility_with_dimensionality(self, tmp_path):
        template = tmp_path / "npt.mdp"
        template.write_text(
            "pcoupl = Berendsen\n"
            "pcoupltype = semiisotropic\n"
            "ref_p = 1.0 1.0\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(system_type="gel")
        copy_and_patch_mdp(
            template,
            output,
            {"npt_barostat": "Parrinello-Rahman"},
            validate=True,
            ctx=ctx,
            stage="npt",
        )
        content = output.read_text()
        assert "compressibility = 1.0e-6 1.0e-6" in content

    def test_pr_gel_water_like_compressibility_warns_without_override(self, tmp_path):
        template = tmp_path / "npt.mdp"
        template.write_text(
            "pcoupl = Berendsen\n"
            "pcoupltype = isotropic\n"
            "compressibility = 4.5e-5\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(system_type="gel")
        _, _, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"npt_barostat": "Parrinello-Rahman"},
            validate=True,
            ctx=ctx,
            stage="npt",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("water-like compressibility" in w.lower() for w in warns)
        content = output.read_text()
        assert "compressibility = 4.5e-5" in content


class TestCommModeAngularPbcSafety:
    def test_angular_with_pbc_overrides_to_linear_by_default(self, tmp_path):
        template = tmp_path / "md.mdp"
        template.write_text(
            "pbc = xyz\n"
            "tc-grps = System\n"
            "comm-grps = System\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(comm_mode="Angular")
        _, _, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"comm_grps_policy": "system"},
            validate=True,
            ctx=ctx,
            stage="md",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("overriding comm-mode to linear" in w.lower() for w in warns)
        assert "comm-mode = Linear" in output.read_text()

    def test_angular_with_pbc_strict_fails_fast(self, tmp_path):
        template = tmp_path / "md.mdp"
        template.write_text(
            "pbc = xyz\n"
            "tc-grps = System\n"
            "comm-grps = System\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(comm_mode="Angular", strict_mdp_validation=True)
        try:
            copy_and_patch_mdp(
                template,
                output,
                {"comm_grps_policy": "system"},
                validate=True,
                ctx=ctx,
                stage="md",
            )
        except MDPPatchError as e:
            assert "unsafe/ill-defined" in str(e)
        else:
            raise AssertionError("Expected MDPPatchError in strict mode for Angular+PBC")

    def test_angular_with_pbc_expert_opt_in_keeps_setting(self, tmp_path):
        template = tmp_path / "md.mdp"
        template.write_text(
            "pbc = xyz\n"
            "tc-grps = System\n"
            "comm-grps = System\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(
            comm_mode="Angular",
            allow_comm_mode_angular_with_pbc=True,
        )
        _, _, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"comm_grps_policy": "system"},
            validate=True,
            ctx=ctx,
            stage="md",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("high visibility warning" in w.lower() for w in warns)
        assert "comm-mode = Angular" in output.read_text()


class TestTinyTcGrpsGuard:
    def test_auto_detect_tiny_group_falls_back_to_system(self, tmp_path):
        ndx = tmp_path / "index.ndx"
        ndx.write_text(
            "[ Tiny ]\n" + " ".join(str(i) for i in range(1, 11)) + "\n"
            "[ Bulk ]\n" + " ".join(str(i) for i in range(11, 211)) + "\n"
        )
        template = tmp_path / "nvt.mdp"
        template.write_text(
            "tc-grps = System\n"
            "tau_t = 0.1\n"
            "ref_t = 300\n"
            "tcoupl = v-rescale\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(tc_grps_ndx_path=str(ndx), tc_grps_min_atoms=50)
        _, _, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"tc_grps_mode": "split"},
            validate=True,
            ctx=ctx,
            stage="nvt",
            return_diagnostics=True,
        )
        assert errors == []
        assert any("tiny thermostat group" in w.lower() for w in warns)
        assert "tc-grps = System" in output.read_text()

    def test_explicit_tiny_group_strict_errors(self, tmp_path):
        ndx = tmp_path / "index.ndx"
        ndx.write_text(
            "[ Li ]\n1 2 3 4\n"
            "[ Polymer ]\n" + " ".join(str(i) for i in range(5, 205)) + "\n"
        )
        template = tmp_path / "nvt.mdp"
        template.write_text(
            "tc-grps = System\n"
            "tau_t = 0.1 0.1\n"
            "ref_t = 300 300\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(
            tc_grps_groups="Li Polymer",
            tc_grps_ndx_path=str(ndx),
            strict_mdp_validation=True,
            tc_grps_min_atoms=50,
        )
        try:
            copy_and_patch_mdp(
                template,
                output,
                {"tc_grps_mode": "split"},
                validate=True,
                ctx=ctx,
                stage="nvt",
            )
        except MDPPatchError as e:
            assert "tiny thermostat group" in str(e).lower()
            assert "allow-small-tc-grps" in str(e)
        else:
            raise AssertionError("Expected strict error for explicit tiny tc-grps")


class TestSystemModeRewriteSafety:
    def test_system_mode_rewrites_split_template_state(self, tmp_path):
        template = tmp_path / "npt.mdp"
        template.write_text(
            "tc-grps = Polymer NonPolymer\n"
            "tau_t = 0.1 0.1\n"
            "ref_t = 300 300\n"
            "comm-grps = Polymer NonPolymer\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(tc_grps_mode="system", comm_mode="Linear")
        copy_and_patch_mdp(
            template,
            output,
            {"tc_grps_mode": "system"},
            validate=True,
            ctx=ctx,
            stage="npt",
        )
        content = output.read_text()
        assert "tc-grps = System" in content
        assert "tau_t = 0.1" in content
        assert "ref_t = 300" in content
        assert "comm-grps = System" in content

    def test_system_fallback_refuses_nonuniform_vector_collapse(self, tmp_path):
        ndx = tmp_path / "index.ndx"
        ndx.write_text(
            "[ Tiny ]\n" + " ".join(str(i) for i in range(1, 11)) + "\n"
            "[ Bulk ]\n" + " ".join(str(i) for i in range(11, 211)) + "\n"
            "[ System ]\n" + " ".join(str(i) for i in range(1, 211)) + "\n"
        )
        template = tmp_path / "nvt.mdp"
        template.write_text(
            "tc-grps = Tiny Bulk\n"
            "tau_t = 0.5 0.1\n"
            "ref_t = 300 305\n"
            "tcoupl = v-rescale\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext(
            tc_grps_mode="auto",
            tc_grps_auto_split_groups="Tiny Bulk",
            tc_grps_ndx_path=str(ndx),
            tc_grps_min_atoms=50,
        )
        try:
            copy_and_patch_mdp(
                template,
                output,
                {"tc_grps_mode": "split"},
                validate=True,
                ctx=ctx,
                stage="nvt",
            )
        except MDPPatchError as e:
            assert "Cannot safely collapse tau_t" in str(e)
        else:
            raise AssertionError("Expected MDPPatchError for non-uniform split->system collapse")


class TestAutoSplitSourceOfTruth:
    def test_auto_decision_groups_are_applied_without_first_two_fallback(self, tmp_path):
        ndx = tmp_path / "index.ndx"
        ndx.write_text(
            "[ GroupA ]\n" + " ".join(str(i) for i in range(1, 81)) + "\n"
            "[ GroupB ]\n" + " ".join(str(i) for i in range(81, 161)) + "\n"
            "[ Polymer ]\n" + " ".join(str(i) for i in range(1, 401)) + "\n"
            "[ NonPolymer ]\n" + " ".join(str(i) for i in range(401, 801)) + "\n"
            "[ System ]\n" + " ".join(str(i) for i in range(1, 801)) + "\n"
        )
        ctx = MockContext(
            tc_grps_mode="auto",
            active_ndx_path=str(ndx),
            tc_grps_min_atoms=50,
            tc_grps_min_fraction=0.01,
            comm_grps_policy="system",
        )
        overrides = get_mdp_overrides_for_stage(ctx, "nvt")
        assert overrides.get("tc_grps_mode") == "split"
        assert getattr(ctx, "tc_grps_auto_split_groups", None) == "Polymer NonPolymer"

        template = tmp_path / "nvt.mdp"
        template.write_text(
            "tc-grps = System\n"
            "tau_t = 0.1 0.1\n"
            "ref_t = 300 300\n"
            "comm-grps = System\n"
            "comm-mode = Linear\n"
            "tcoupl = v-rescale\n"
        )
        output = tmp_path / "out.mdp"
        copy_and_patch_mdp(template, output, overrides, validate=True, ctx=ctx, stage="nvt")
        content = output.read_text()
        assert "tc-grps = Polymer NonPolymer" in content
        assert "tc-grps = GroupA GroupB" not in content


class TestNdxInlineCommentParsing:
    def test_ndx_parser_ignores_inline_comments_for_groups_and_sizes(self, tmp_path):
        ndx = tmp_path / "index.ndx"
        ndx.write_text(
            "[ Polymer ] ; polymer atoms\n"
            "1 2 3 ; tail comment\n"
            "4 5\n"
            "[ NonPolymer ] ; solvent\n"
            "6 7 ; ions\n"
        )
        groups = parse_ndx_groups(ndx)
        sizes = parse_ndx_group_sizes(ndx)
        assert groups == ["Polymer", "NonPolymer"]
        assert sizes["Polymer"] == 5
        assert sizes["NonPolymer"] == 2


class TestPolicyNoneAndNstepsRobustness:
    def test_comm_policy_none_skips_comm_mismatch_validation(self, tmp_path):
        template = tmp_path / "md.mdp"
        template.write_text(
            "tc-grps = A B\n"
            "tau_t = 1.0 1.0\n"
            "ref_t = 300 300\n"
            "comm-grps = A\n"
            "comm-mode = Linear\n"
        )
        output = tmp_path / "out.mdp"
        ctx = MockContext()
        _, _, errors, warns = copy_and_patch_mdp(
            template,
            output,
            {"comm_grps_policy": "none"},
            validate=True,
            ctx=ctx,
            stage="md",
            return_diagnostics=True,
        )
        assert errors == []
        assert not any("comm-grps has" in w for w in warns)
        params = {
            "tc-grps": "A B",
            "tau_t": "1.0 1.0",
            "ref_t": "300 300",
            "comm-grps": "A",
            "comm-mode": "None",
        }
        vec_errors, vec_warns = _validate_vector_lengths(params, "md", strict=False)
        assert vec_errors == []
        assert vec_warns == []

    def test_invalid_nsteps_raises_helpful_mdppatcherror(self, tmp_path):
        template = tmp_path / "nvt.mdp"
        template.write_text("nsteps = 1000\n")
        output = tmp_path / "out.mdp"
        try:
            copy_and_patch_mdp(
                template,
                output,
                {"nsteps": "abc"},
                validate=True,
                ctx=MockContext(),
                stage="nvt",
            )
        except MDPPatchError as e:
            text = str(e)
            assert "Stage 'nvt'" in text
            assert "nsteps" in text
            assert "'abc'" in text
        else:
            raise AssertionError("Expected MDPPatchError for non-integer nsteps")


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
