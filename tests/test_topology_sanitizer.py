from pathlib import Path
from types import SimpleNamespace
import math
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.charge_neutrality import (
    ChargeNeutralityChecker,
    ChargeNeutralityError,
    ChargeParser,
)
from pipeline.stages.topology_sanitizer import (
    ATOMS_ROW_STATUS_INVALID,
    ATOMS_ROW_STATUS_PARSED,
    ATOMS_ROW_STATUS_UNSUPPORTED,
    DEFAULT_MOLECULE_SIG_QUANT_STEP,
    MOLECULETYPE_COMPARE_UNCERTAIN,
    SANITIZER_BLOCK_BEGIN,
    SANITIZER_BLOCK_END,
    SanitizerError,
    TOP_MOLECULES_ENTRY_STATUS_PARSED,
    TOP_MOLECULES_ENTRY_STATUS_UNRESOLVED,
    TopMoleculesParseResult,
    TopologySanitizerMixin,
    _matches_protected_pattern,
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


def test_nested_include_resolution_uses_current_included_file_directory(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    forcefield_dir = tmp_path / "IN/forcefield/oplsaa"
    sanitizer._forcefield_dir = forcefield_dir

    root_itp = forcefield_dir / "gromacs/root.itp"
    child_itp = tmp_path / "IN/systems/SYS/gromacs/itp/nested/child.itp"
    leaf_itp = tmp_path / "IN/systems/SYS/gromacs/itp/nested/deeper/leaf.itp"

    _write(root_itp, '#include "nested/child.itp"\n[ moleculetype ]\nROOT 3\n')
    _write(child_itp, '#include "deeper/leaf.itp"\n[ moleculetype ]\nCHILD 3\n')
    _write(leaf_itp, "[ moleculetype ]\nLEAF 3\n")

    names = sanitizer._get_moleculetype_names(root_itp, ctx=ctx)

    assert names == ["LEAF", "CHILD", "ROOT"]


@pytest.mark.parametrize(
    ("priority", "expected_first"),
    [
        ("sanitized_first", "user_include"),
        ("forcefield_first", "IN/systems/SYS/gromacs/top"),
    ],
)
def test_build_include_search_paths_uses_single_priority_source(
    tmp_path,
    priority,
    expected_first,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(
        tmp_path,
        itp_include_dirs="user_include",
        itp_include_dirs_priority=priority,
    )
    sanitizer._forcefield_dir = tmp_path / "IN/forcefield/oplsaa"
    sanitizer._sanitized_current_dir = tmp_path / "IN/systems/SYS/gromacs/itp_sanitized_current"

    base_itp = tmp_path / "IN/systems/SYS/gromacs/itp/base.itp"
    for path in [
        base_itp.parent,
        tmp_path / "user_include",
        tmp_path / "IN/systems/SYS/gromacs/top",
        tmp_path / "IN/systems/SYS/gromacs/itp",
        tmp_path / "IN/systems/SYS/htpolynet/itp",
        tmp_path / "IN/forcefield/oplsaa/gromacs",
        tmp_path / "IN/systems/SYS/gromacs/itp_sanitized_current",
    ]:
        path.mkdir(parents=True, exist_ok=True)
    _write(base_itp, "")

    search_paths = sanitizer._build_include_search_paths(ctx, base_itp)
    search_strings = [str(path.relative_to(tmp_path)).replace("\\", "/") for path in search_paths]

    assert search_strings[0] == "IN/systems/SYS/gromacs/itp"
    assert search_strings[1] == expected_first


def test_resolve_include_strict_fails_on_shadowing_when_not_allowed(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(
        tmp_path,
        strict_include_resolution=True,
        allow_include_shadowing=False,
        itp_include_dirs="user_include",
        itp_include_dirs_priority="sanitized_first",
    )
    forcefield_dir = tmp_path / "IN/forcefield/oplsaa"
    sanitizer._forcefield_dir = forcefield_dir

    root_itp = forcefield_dir / "gromacs/root.itp"
    user_leaf = tmp_path / "user_include/leaf.itp"
    ff_leaf = forcefield_dir / "gromacs/leaf.itp"
    _write(root_itp, "#include \"leaf.itp\"\n")
    _write(user_leaf, "[ moleculetype ]\nUSER 3\n")
    _write(ff_leaf, "[ moleculetype ]\nFF 3\n")

    with pytest.raises(SanitizerError, match="Include shadowing violation"):
        sanitizer._resolve_include_in_dirs(
            "leaf.itp",
            root_itp,
            [],
            strict=True,
            check_shadowing=True,
            allow_shadowing=False,
            ctx=ctx,
        )


def test_resolve_include_non_strict_reports_shadowing_warning(tmp_path, capsys):
    sanitizer = DummySanitizer()
    ctx = DummyContext(
        tmp_path,
        allow_include_shadowing=False,
        itp_include_dirs="user_include",
        itp_include_dirs_priority="sanitized_first",
    )
    forcefield_dir = tmp_path / "IN/forcefield/oplsaa"
    sanitizer._forcefield_dir = forcefield_dir

    root_itp = forcefield_dir / "gromacs/root.itp"
    user_leaf = tmp_path / "user_include/leaf.itp"
    ff_leaf = forcefield_dir / "gromacs/leaf.itp"
    _write(root_itp, "#include \"leaf.itp\"\n")
    _write(user_leaf, "[ moleculetype ]\nUSER 3\n")
    _write(ff_leaf, "[ moleculetype ]\nFF 3\n")

    resolved = sanitizer._resolve_include_in_dirs(
        "leaf.itp",
        root_itp,
        [],
        strict=False,
        check_shadowing=True,
        allow_shadowing=False,
        ctx=ctx,
    )

    assert resolved == ff_leaf.resolve()
    assert "Include shadowing detected" in capsys.readouterr().out


def test_get_moleculetype_names_strict_fails_closed_on_unresolved_conditionals(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_include_resolution=True)
    itp = tmp_path / "macro.itp"
    _write(
        itp,
        "#ifdef USE_ALT\n"
        "[ moleculetype ]\n"
        "ALT 3\n"
        "#endif\n",
    )

    with pytest.raises(SanitizerError, match="Cannot safely moleculetype name extraction"):
        sanitizer._get_moleculetype_names(itp, ctx=ctx)


def test_get_moleculetype_names_strict_fails_closed_on_macro_like_include(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_include_resolution=True)
    itp = tmp_path / "macro_include.itp"
    _write(
        itp,
        "#include SOME_MACRO\n"
        "[ moleculetype ]\n"
        "BASE 3\n",
    )

    with pytest.raises(SanitizerError, match="#include\\(macro-like\\)"):
        sanitizer._get_moleculetype_names(itp, ctx=ctx)


def test_get_moleculetype_names_non_strict_reports_degraded_macro_fallback(
    tmp_path,
    capsys,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    itp = tmp_path / "macro.itp"
    _write(
        itp,
        "#ifdef USE_ALT\n"
        "[ moleculetype ]\n"
        "ALT 3\n"
        "#endif\n"
        "[ moleculetype ]\n"
        "BASE 3\n",
    )

    names = sanitizer._get_moleculetype_names(itp, ctx=ctx)

    assert names == ["BASE"]
    assert "Conservative Python fallback used for moleculetype name extraction" in capsys.readouterr().out


def test_expand_include_closure_reports_degraded_conditional_includes_in_non_strict_mode(
    tmp_path,
    capsys,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    root_itp = tmp_path / "IN/systems/SYS/gromacs/itp/root.itp"
    leaf_itp = tmp_path / "IN/systems/SYS/gromacs/itp/leaf.itp"
    _write(
        root_itp,
        "#ifdef USE_LEAF\n"
        "#include \"leaf.itp\"\n"
        "#endif\n"
        "[ moleculetype ]\n"
        "ROOT 3\n",
    )
    _write(leaf_itp, "[ moleculetype ]\nLEAF 3\n")

    expanded, class_map, ff_map = sanitizer._expand_itp_include_closure(
        [root_itp],
        ctx,
        {str(root_itp): "staged"},
        {str(root_itp): "UNKNOWN"},
    )

    assert expanded == [root_itp.resolve()]
    assert class_map[str(root_itp.resolve())] == "staged"
    assert ff_map[str(root_itp.resolve())] == "UNKNOWN"
    assert "Conservative Python fallback used for include closure discovery" in capsys.readouterr().out


def test_validate_system_top_includes_strict_fails_closed_on_macro_like_include(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_include_resolution=True)
    top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    _write(top, "#include SOME_MACRO\n[molecules]\nA 1\n")

    with pytest.raises(SanitizerError, match="#include\\(macro-like\\)"):
        sanitizer._validate_system_top_includes(top, ctx)


def test_validate_system_top_includes_non_strict_reports_degraded_conditionals(
    tmp_path,
    capsys,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    leaf = tmp_path / "IN/systems/SYS/gromacs/top/leaf.itp"
    _write(
        top,
        "#ifdef USE_LEAF\n"
        "#include \"leaf.itp\"\n"
        "#endif\n"
        "[molecules]\n"
        "A 1\n",
    )
    _write(leaf, "")

    sanitizer._validate_system_top_includes(top, ctx)

    assert "Conservative Python fallback used for system.top include validation" in capsys.readouterr().out


def test_update_system_top_includes_strict_fails_closed_on_preprocessor_ambiguity(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_include_resolution=True)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "#ifdef USE_FF\n"
        "#include \"forcefield.itp\"\n"
        "#endif\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    with pytest.raises(SanitizerError, match="Cannot safely system.top managed-block update"):
        sanitizer._update_system_top_includes(
            existing_top,
            combined_atomtypes,
            [molecule_itp],
            top_dir=existing_top.parent,
            ctx=ctx,
        )


def test_update_system_top_includes_inserts_new_managed_block_after_forcefield_include(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "; Banner\n"
        "#define USE_SANITIZED\n"
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        "; User comment\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    updated = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_itp],
        top_dir=existing_top.parent,
        ctx=ctx,
    )

    lines = updated.splitlines()
    assert lines.count(SANITIZER_BLOCK_BEGIN) == 1
    assert lines.count(SANITIZER_BLOCK_END) == 1
    assert "[ defaults ]" not in updated
    assert lines.index(SANITIZER_BLOCK_BEGIN) > lines.index('#include "oplsaa.ff/forcefield.itp"')
    assert lines.index(SANITIZER_BLOCK_BEGIN) < lines.index("; User comment")
    assert '#include "combined_atomtypes_current.itp"' in updated
    assert '#include "../itp/A.itp"' in updated


def test_update_system_top_includes_uses_safe_anchor_before_non_forcefield_includes(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "; Banner\n"
        "#define KEEP_THIS\n"
        "#include \"user_custom_before_block.itp\"\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    updated = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_itp],
        defaults=SimpleNamespace(
            nbfunc=1,
            comb_rule=2,
            gen_pairs="yes",
            fudge_lj=0.5,
            fudge_qq=0.8333,
        ),
        top_dir=existing_top.parent,
        ctx=ctx,
    )

    lines = updated.splitlines()
    assert lines.index(SANITIZER_BLOCK_BEGIN) > lines.index("#define KEEP_THIS")
    assert lines.index(SANITIZER_BLOCK_BEGIN) < lines.index('#include "user_custom_before_block.itp"')
    assert "[ defaults ]" in updated


def test_update_system_top_includes_replaces_existing_block_without_duplicate_includes(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "; Combined atomtypes\n"
        "#include \"old_combined.itp\"\n"
        "#include \"extra.itp\"\n"
        "\n"
        "; Molecule definitions (sanitized)\n"
        "#include \"../itp/OLD.itp\"\n"
        "\n"
        f"{SANITIZER_BLOCK_END}\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    updated = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_itp],
        top_dir=existing_top.parent,
        ctx=ctx,
        extra_includes=["extra.itp", "extra.itp"],
    )

    assert updated.count(SANITIZER_BLOCK_BEGIN) == 1
    assert updated.count(SANITIZER_BLOCK_END) == 1
    assert "old_combined.itp" not in updated
    assert "../itp/OLD.itp" not in updated
    assert updated.count('#include "extra.itp"') == 1
    assert updated.count('#include "../itp/A.itp"') == 1
    assert "[ system ]\nDemo" in updated


def test_update_system_top_includes_is_idempotent_on_rerun(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_a = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    molecule_b = tmp_path / "IN/systems/SYS/gromacs/itp/B.itp"
    _write(
        existing_top,
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "B 1\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_a, "[ moleculetype ]\nA 3\n")
    _write(molecule_b, "[ moleculetype ]\nB 3\n")

    first = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_b, molecule_a, molecule_a],
        top_dir=existing_top.parent,
        ctx=ctx,
        extra_includes=["extra.itp", "extra.itp"],
    )
    _write(existing_top, first)
    second = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_b, molecule_a, molecule_a],
        top_dir=existing_top.parent,
        ctx=ctx,
        extra_includes=["extra.itp", "extra.itp"],
    )

    assert second == first
    assert second.count(SANITIZER_BLOCK_BEGIN) == 1
    assert second.count(SANITIZER_BLOCK_END) == 1
    assert second.count('#include "extra.itp"') == 1
    assert second.count('#include "../itp/B.itp"') == 1
    assert second.count('#include "../itp/A.itp"') == 1


def test_update_system_top_includes_strict_rejects_multiple_managed_blocks(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, strict_top_update=True)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "#include \"old_a.itp\"\n"
        f"{SANITIZER_BLOCK_END}\n"
        "; keep me\n"
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "#include \"old_b.itp\"\n"
        f"{SANITIZER_BLOCK_END}\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    with pytest.raises(SanitizerError, match="Multiple sanitizer managed blocks detected"):
        sanitizer._update_system_top_includes(
            existing_top,
            combined_atomtypes,
            [molecule_itp],
            top_dir=existing_top.parent,
            ctx=ctx,
        )


def test_update_system_top_includes_non_strict_collapses_multiple_blocks_without_losing_user_content(
    tmp_path,
    capsys,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "#include \"old_a.itp\"\n"
        f"{SANITIZER_BLOCK_END}\n"
        "#include \"user_custom.itp\"\n"
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "#include \"old_b.itp\"\n"
        f"{SANITIZER_BLOCK_END}\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    updated = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_itp],
        top_dir=existing_top.parent,
        ctx=ctx,
    )

    assert "Multiple sanitizer managed blocks detected" in capsys.readouterr().out
    assert updated.count(SANITIZER_BLOCK_BEGIN) == 1
    assert updated.count(SANITIZER_BLOCK_END) == 1
    assert '#include "user_custom.itp"' in updated
    assert "old_a.itp" not in updated
    assert "old_b.itp" not in updated


def test_update_system_top_includes_non_strict_rebuilds_safe_malformed_managed_fragments(
    tmp_path,
    capsys,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "#include \"old_combined.itp\"\n"
        f"{SANITIZER_BLOCK_END}\n"
        f"{SANITIZER_BLOCK_END}\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    updated = sanitizer._update_system_top_includes(
        existing_top,
        combined_atomtypes,
        [molecule_itp],
        top_dir=existing_top.parent,
        ctx=ctx,
    )

    assert "Malformed sanitizer managed markers in system.top" in capsys.readouterr().out
    assert updated.count(SANITIZER_BLOCK_BEGIN) == 1
    assert updated.count(SANITIZER_BLOCK_END) == 1
    assert "old_combined.itp" not in updated
    assert '#include "combined_atomtypes_current.itp"' in updated
    assert '#include "../itp/A.itp"' in updated
    assert "[ system ]\nDemo" in updated


def test_update_system_top_includes_non_strict_fails_closed_on_unsafe_malformed_managed_fragments(
    tmp_path,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    existing_top = tmp_path / "IN/systems/SYS/gromacs/top/system.top"
    combined_atomtypes = tmp_path / "IN/systems/SYS/gromacs/top/combined_atomtypes_current.itp"
    molecule_itp = tmp_path / "IN/systems/SYS/gromacs/itp/A.itp"
    _write(
        existing_top,
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        f"{SANITIZER_BLOCK_BEGIN}\n"
        "#include \"old_combined.itp\"\n"
        "[ system ]\n"
        "Demo\n"
        "[ molecules ]\n"
        "A 1\n",
    )
    _write(combined_atomtypes, "")
    _write(molecule_itp, "[ moleculetype ]\nA 3\n")

    with pytest.raises(
        SanitizerError,
        match="Cannot safely rebuild non-strict system.top because malformed managed-block fragments overlap non-managed topology content",
    ):
        sanitizer._update_system_top_includes(
            existing_top,
            combined_atomtypes,
            [molecule_itp],
            top_dir=existing_top.parent,
            ctx=ctx,
        )


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

    comparison = sanitizer._compare_moleculetype_ir(itp1, itp2, "MOL", ctx=ctx)

    assert comparison.status == MOLECULETYPE_COMPARE_UNCERTAIN
    assert sanitizer._moleculetype_signatures_differ(itp1, itp2, "MOL", ctx=ctx) is True
    assert "Treating definitions as conflicting" in capsys.readouterr().out


def test_same_name_moleculetype_comparison_uses_ir_for_valid_supported_rows(tmp_path):
    sanitizer = DummySanitizer()
    itp1 = tmp_path / "a.itp"
    itp2 = tmp_path / "b.itp"
    _write(
        itp1,
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "1 CT 1 MOL C1 1 -0.1250000 12.0110000\n"
        "2 HC 1 MOL H1 1 0.1250000 1.0080000\n",
    )
    _write(
        itp2,
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "1 CT 1 MOL C1 1 -0.1250004 12.0110004\n"
        "2 HC 1 MOL H1 1 0.1250004 1.0080004\n",
    )

    assert sanitizer._moleculetype_signatures_differ(itp1, itp2, "MOL") is False


def test_same_name_moleculetype_comparison_detects_local_atomtype_lj_conflict(tmp_path):
    sanitizer = DummySanitizer()
    itp1 = tmp_path / "a.itp"
    itp2 = tmp_path / "b.itp"
    base = (
        "[ atomtypes ]\n"
        "{atomtype_row}\n"
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "1 CT 1 MOL C1 1 0.000000 12.011\n"
    )
    _write(
        itp1,
        base.format(atomtype_row="CT 12.011 0.000 A 0.300000 0.100000"),
    )
    _write(
        itp2,
        base.format(atomtype_row="CT 12.011 0.000 A 0.300000 0.200000"),
    )

    assert sanitizer._moleculetype_signatures_differ(itp1, itp2, "MOL") is True


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


@pytest.mark.parametrize(
    ("name", "patterns"),
    [
        ("PEG-400", frozenset({"PEG_400"})),
        ("PEG_400", frozenset({"PEG-400"})),
        ("TFSI-1", frozenset({"TFSI"})),
        ("TFSI_1", frozenset({"TFSI"})),
        ("TFSI-alpha", frozenset({"TFSI"})),
    ],
)
def test_matches_protected_pattern_supports_separator_bearing_names(name, patterns):
    assert _matches_protected_pattern(name, patterns) is True


@pytest.mark.parametrize(
    ("name", "pattern"),
    [
        ("LIPID_GPC", "PC"),
        ("POLYMER_CHAIN", "PC"),
        ("SODIUM", "NA"),
        ("LIPID", "LI"),
    ],
)
def test_matches_protected_pattern_keeps_short_tokens_boundary_safe(name, pattern):
    assert _matches_protected_pattern(name, frozenset({pattern})) is False


@pytest.mark.parametrize("mol_name", ["PEG-400", "PEG_400"])
def test_charge_fix_protect_resnames_supports_separator_names_in_classification(
    tmp_path,
    mol_name,
):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path, charge_fix_protect_resnames="PEG-400,PEG_400")
    itp = tmp_path / f"{mol_name}.itp"
    _write(itp, _simple_itp(mol_name))

    classifications, audit = sanitizer._classify_ionic_moleculetypes({mol_name: itp}, ctx=ctx)

    assert classifications[mol_name] == "protected_configured"
    assert audit["details"][mol_name]["source"] == "charge_fix_protect_resnames"


def test_degraded_or_protected_targets_are_removed_from_default_checker_allowlist(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    broken = tmp_path / "broken.itp"
    solvent = tmp_path / "pc.itp"
    _write(
        broken,
        "[ moleculetype ]\n"
        "BROKEN 3\n"
        "[ atoms ]\n"
        "#ifdef SOMEFLAG\n"
        "1 CT 1 BROKEN C1 1 0.000000 12.011\n"
        "#endif\n",
    )
    _write(solvent, _simple_itp("PC"))

    classifications, audit = sanitizer._classify_ionic_moleculetypes(
        {"BROKEN": broken, "PC": solvent},
        ctx=ctx,
    )
    allowlist = sanitizer._derive_charge_fix_checker_allowlist(
        classifications,
        audit,
        explicit_allowlist=None,
    )

    assert "BROKEN" not in allowlist
    assert "PC" not in allowlist
    assert allowlist == set()


@pytest.mark.parametrize(
    "atoms_rows",
    [
        (
            "1 CT 1 MOL C1 -0.125000 12.011\n"
            "2 HC 1 MOL H1 0.125000 1.008\n"
        ),
        (
            "1 CT 1 MOL C1 3-0.125000 12.011\n"
            "2 HC 1 MOL H1 4+0.125000 1.008\n"
        ),
    ],
)
def test_charge_fix_classification_accepts_supported_7_column_variants(
    tmp_path,
    atoms_rows,
):
    sanitizer = DummySanitizer()
    itp = tmp_path / "mol.itp"
    _write(
        itp,
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        f"{atoms_rows}",
    )

    classifications, audit = sanitizer._classify_ionic_moleculetypes({"MOL": itp})

    assert classifications["MOL"] == "neutral"
    assert audit["details"]["MOL"]["status"] == "exact"
    assert audit["details"]["MOL"]["total_charge"] == pytest.approx(0.0)


@pytest.mark.parametrize(
    "atoms_rows",
    [
        (
            "1 CT 1 MOL C1 -0.125000 12.011\n"
            "2 HC 1 MOL H1 0.125000 1.008\n"
        ),
        (
            "1 CT 1 MOL C1 3-0.125000 12.011\n"
            "2 HC 1 MOL H1 4+0.125000 1.008\n"
        ),
    ],
)
def test_charge_parser_remains_compatible_with_supported_7_column_variants(
    tmp_path,
    atoms_rows,
):
    itp = tmp_path / "mol.itp"
    _write(
        itp,
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        f"{atoms_rows}",
    )

    result = ChargeParser.get_molecule_charge(itp, "MOL", include_details=True)

    assert result.total_charge == pytest.approx(0.0)
    assert result.atom_count == 2
    assert result.unparsed_atoms_lines == []


@pytest.mark.parametrize(
    "atoms_row",
    [
        "1 CT 1 MOL C1 1 -0.125000",
        "1 CT 1 MOL C1 0.125000-12.011",
        "1 CT 1 MOL C1 1 0.125000 -12.011",
    ],
)
def test_charge_parser_fails_closed_for_stage1_unsupported_atoms_rows(
    tmp_path,
    atoms_row,
):
    itp = tmp_path / "mol.itp"
    _write(
        itp,
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        f"{atoms_row}\n",
    )

    with pytest.raises(ChargeNeutralityError, match="Unparseable \\[ atoms \\] data line"):
        ChargeParser.get_molecule_charge(itp, "MOL", include_details=True)


def test_allowlisted_polymer_requires_spread_safe_before_any_correction(tmp_path):
    itp = tmp_path / "polymer.itp"
    _write(
        itp,
        "[ moleculetype ]\n"
        "POLYMER_CHAIN 3\n"
        "[ atoms ]\n"
        "1 CT 1 POLYMER_CHAIN C1 1 0.500000 12.011\n"
        "2 HC 1 POLYMER_CHAIN H1 1 0.500000 1.008\n",
    )

    checker = ChargeNeutralityChecker(
        threshold=2.0,
        target_allowlist={"POLYMER_CHAIN"},
        strict_mode=False,
        protected_polymers={"POLYMER_CHAIN"},
        polymer_correction_method="skip_if_small",
        allow_non_rounding_correction=True,
    )

    result = checker.run(
        molecule_itp_paths={"POLYMER_CHAIN": itp},
        molecule_counts={"POLYMER_CHAIN": 1},
        output_dir=tmp_path / "out",
        run_id="TEST",
        target_molecule="POLYMER_CHAIN",
    )

    assert result.corrected is False
    assert result.neutrality_status == "polymer_skipped_non_strict"
    assert "spread_safe" in (result.refuse_reason or "")
    assert list((tmp_path / "out").glob("*corrected*.itp")) == []


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


def test_map_sanitized_itps_to_molecules_maps_all_expected_names_deterministically(tmp_path):
    sanitizer = DummySanitizer()
    ctx = DummyContext(tmp_path)
    sanitized_current_dir = tmp_path / "itp_sanitized_current"
    _write(sanitized_current_dir / "z_b.itp", _simple_itp("B"))
    _write(sanitized_current_dir / "a_a.itp", _simple_itp("A"))

    mapping = sanitizer._map_sanitized_itps_to_molecules(
        sanitized_current_dir,
        {"B", "A"},
        ctx=ctx,
    )

    assert list(mapping.keys()) == ["A", "B"]
    assert mapping["A"].name == "a_a.itp"
    assert mapping["B"].name == "z_b.itp"


def test_ordered_molecule_itp_paths_preserves_first_seen_order_and_dedupes_files(tmp_path):
    sanitizer = DummySanitizer()
    a_itp = tmp_path / "a.itp"
    b_itp = tmp_path / "b.itp"
    _write(a_itp, "")
    _write(b_itp, "")

    result = sanitizer._ordered_molecule_itp_paths(
        [("B", 1), ("A", 2), ("D", 1), ("B", 3), ("C", 0)],
        {
            "A": a_itp,
            "B": b_itp,
            "D": a_itp,
        },
        tmp_path,
    )

    assert result == [b_itp, a_itp]


def test_ordered_molecule_itp_paths_fails_when_expected_include_is_missing(tmp_path):
    sanitizer = DummySanitizer()
    a_itp = tmp_path / "a.itp"
    _write(a_itp, "")

    with pytest.raises(SanitizerError, match="No sanitized ITP found for molecule 'B'"):
        sanitizer._ordered_molecule_itp_paths(
            [("A", 1), ("B", 1)],
            {"A": a_itp},
            tmp_path,
        )


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
