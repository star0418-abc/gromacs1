#!/usr/bin/env python3
"""Simple test runner for MDP patcher tests."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from tests.test_mdp_patcher import (
    TestP01StagePolicyLoop,
    TestP02ExpertRestartGuard,
    TestP03GenTempValidation,
    TestP0BContinuationMessaging,
    TestSemanticKeyOverrideCollision,
    TestTemperatureInjection,
)

def run_tests():
    passed = 0
    failed = 0
    
    print("=== P0.1 Tests: Stage Policy Loop ===")
    t = TestP01StagePolicyLoop()
    for name in ["test_nvt_defaults_fresh_start", "test_npt_defaults_continuation", 
                 "test_md_defaults_continuation", "test_force_gen_vel_overrides_npt",
                 "test_allow_continuation_override_skips_continuation"]:
        try:
            getattr(t, name)()
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print("\n=== P0.2 Tests: Expert Restart Guard ===")
    t2 = TestP02ExpertRestartGuard()
    for name in ["test_gen_vel_no_continuation_no_warning_nonstrict",
                 "test_gen_vel_no_continuation_no_error_strict",
                 "test_expert_restart_allowed_with_flag",
                 "test_expert_restart_blocked_when_checkpoint_available",
                 "test_contradictory_combo_always_error"]:
        try:
            getattr(t2, name)()
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print("\n=== P0.3 Tests: gen_temp Validation ===")
    t3 = TestP03GenTempValidation()
    for name in ["test_gen_temp_missing_warning_nonstrict",
                 "test_gen_temp_missing_error_strict",
                 "test_gen_temp_non_numeric_caught",
                 "test_ref_t_non_numeric_caught"]:
        try:
            getattr(t3, name)()
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print("\n=== P0.B Tests: Continuation Messaging ===")
    t4 = TestP0BContinuationMessaging()
    for name in ["test_continuation_yes_produces_info_message"]:
        try:
            getattr(t4, name)()
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print("\n=== Task A Tests: Semantic Key Override Collision ===")
    import tempfile
    from pathlib import Path
    t5 = TestSemanticKeyOverrideCollision()
    for name in ["test_override_underscore_when_template_has_hyphen",
                 "test_override_hyphen_when_template_has_underscore"]:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                getattr(t5, name)(Path(tmp))
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print("\n=== Task C Tests: Temperature Injection ===")
    t6 = TestTemperatureInjection()
    for name in ["test_inject_ref_t_when_missing",
                 "test_inject_gen_temp_when_gen_vel_yes"]:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                getattr(t6, name)(Path(tmp))
            print(f"  PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {name} - {e}")
            failed += 1
    
    print(f"\n=== Results: {passed} passed, {failed} failed ===")
    return failed == 0

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
