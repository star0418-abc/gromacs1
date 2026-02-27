#!/usr/bin/env python3
"""
Export a publication-ready computational_details.md from run_report.json.
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


EXPECTED_CANONICAL_STAGES = ("em", "nvt", "npt", "md")
CANONICAL_LABELS = {
    "em": "EM",
    "nvt": "NVT",
    "npt": "NPT",
    "md": "MD",
}

STAGE_TOKEN_ALIASES = {
    "em": {
        "em",
        "min",
        "minim",
        "minimize",
        "minimization",
        "energymin",
        "energyminimization",
    },
    "nvt": {"nvt"},
    "npt": {"npt"},
    "md": {"md", "prod", "production"},
}

MDP_KEY_SPECS: Sequence[Tuple[str, str]] = (
    ("integrator", "integrator"),
    ("dt", "dt"),
    ("nsteps", "nsteps"),
    ("continuation", "continuation"),
    ("constraints", "constraints"),
    ("constraint-algorithm", "constraint-algorithm"),
    ("lincs-order", "lincs-order"),
    ("lincs-iter", "lincs-iter"),
    ("lincs-warnangle", "lincs-warnangle"),
    ("cutoff-scheme", "cutoff-scheme"),
    ("ns_type", "ns-type"),
    ("nstlist", "nstlist"),
    ("coulombtype", "coulombtype"),
    ("rcoulomb", "rcoulomb"),
    ("rvdw", "rvdw"),
    ("vdwtype", "vdwtype"),
    ("rlist", "rlist"),
    ("pme_order", "pme-order"),
    ("fourierspacing", "fourierspacing"),
    ("DispCorr", "dispcorr"),
    ("tcoupl", "tcoupl"),
    ("tc-grps", "tc-grps"),
    ("tau_t", "tau-t"),
    ("ref_t", "ref-t"),
    ("pcoupl", "pcoupl"),
    ("pcoupltype", "pcoupltype"),
    ("tau_p", "tau-p"),
    ("ref_p", "ref-p"),
    ("compressibility", "compressibility"),
    ("gen_vel", "gen-vel"),
    ("gen_temp", "gen-temp"),
    ("gen_seed", "gen-seed"),
    ("pbc", "pbc"),
)

TRUE_STRINGS = {"1", "true", "yes", "on", "y"}
CONTINUITY_TERMS = ("resume", "restart", "velocity", "checkpoint", "maxwarn")
KEY_POLICIES = (
    "allow_velocity_reset",
    "allow_gro_velocity_restart",
    "resume_allowed",
    "grompp_maxwarn",
)


def _load_json(path: Path) -> Dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RuntimeError(f"Failed to read JSON: {path} ({exc})") from exc


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _normalize_mdp_key(key: str) -> str:
    return key.strip().lower().replace("_", "-")


def _parse_mdp_file(path: Path) -> Dict[str, str]:
    params: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.split(";", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key_norm = _normalize_mdp_key(key)
            params[key_norm] = value.strip()
    return params


def _extract_absolute_mdp_values(params: Dict[str, str]) -> List[Tuple[str, str]]:
    values: List[Tuple[str, str]] = []
    for display_key, normalized_key in MDP_KEY_SPECS:
        raw_value = params.get(normalized_key)
        value = raw_value if raw_value else "not set"
        values.append((display_key, value))
    return values


def _json_scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def _expand_stage_token(token: str) -> List[str]:
    candidates = [token]
    for prefix in ("gmx", "gromacs"):
        if token.startswith(prefix) and len(token) > len(prefix):
            candidates.append(token[len(prefix) :])
        if token.endswith(prefix) and len(token) > len(prefix):
            candidates.append(token[: -len(prefix)])
    return [c for c in candidates if c]


def _canonical_stage_name(stage_name: str) -> Optional[str]:
    cleaned = re.sub(r"[^a-z0-9]+", "_", stage_name.lower()).strip("_")
    if not cleaned:
        return None
    tokens = [token for token in cleaned.split("_") if token]
    expanded = set()
    for token in tokens:
        expanded.update(_expand_stage_token(token))
    expanded.update(_expand_stage_token(cleaned))
    matches = []
    for canonical, aliases in STAGE_TOKEN_ALIASES.items():
        if aliases & expanded:
            matches.append(canonical)
    if len(matches) == 1:
        return matches[0]
    return None


def _parse_require_stages(raw: str) -> List[str]:
    if not raw.strip():
        return []
    parsed: List[str] = []
    for part in raw.split(","):
        stage = part.strip().lower()
        if not stage:
            continue
        if stage not in EXPECTED_CANONICAL_STAGES:
            raise RuntimeError(
                f"Invalid value in --require-stages: {stage}. "
                f"Allowed values: {', '.join(EXPECTED_CANONICAL_STAGES)}"
            )
        if stage not in parsed:
            parsed.append(stage)
    return parsed


def _format_duplicate_stage_error(
    duplicates: Dict[str, List[Tuple[int, str]]]
) -> RuntimeError:
    lines = []
    for canonical in sorted(duplicates):
        mapped = ", ".join(
            f"{name} (index {index})" for index, name in duplicates[canonical]
        )
        lines.append(f"- {canonical}: {mapped}")
    msg = (
        "Multiple run_report stage entries map to the same canonical stage. "
        "Cannot export deterministically:\n" + "\n".join(lines)
    )
    return RuntimeError(msg)


def _collect_stage_records(
    report: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    stages_raw = report.get("stages", [])
    if stages_raw is None:
        stages_raw = []
    if not isinstance(stages_raw, list):
        raise RuntimeError("run_report.json field 'stages' must be a list")

    records: List[Dict[str, Any]] = []
    canonical_hits: Dict[str, List[Tuple[int, str]]] = {}
    for index, entry in enumerate(stages_raw):
        if isinstance(entry, dict):
            stage_name = entry.get("stage")
            stage_name_str = (
                str(stage_name)
                if stage_name is not None and str(stage_name).strip()
                else f"<missing stage name @ index {index}>"
            )
            entry_dict = entry
        else:
            stage_name_str = f"<invalid stage entry @ index {index}>"
            entry_dict = {}
        canonical = _canonical_stage_name(stage_name_str)
        if canonical:
            canonical_hits.setdefault(canonical, []).append((index, stage_name_str))
        records.append(
            {
                "index": index,
                "stage_name": stage_name_str,
                "canonical": canonical,
                "entry": entry_dict,
            }
        )

    duplicates = {
        canonical: hits
        for canonical, hits in canonical_hits.items()
        if len(hits) > 1
    }
    if duplicates:
        raise _format_duplicate_stage_error(duplicates)

    by_canonical: Dict[str, Dict[str, Any]] = {}
    for record in records:
        canonical = record.get("canonical")
        if canonical:
            by_canonical[canonical] = record
    return records, by_canonical


def _resolve_report_path(run_root: Path, raw_path: Any) -> Optional[Path]:
    if raw_path is None:
        return None
    path_str = str(raw_path).strip()
    if not path_str:
        return None
    path = Path(path_str)
    if not path.is_absolute():
        path = run_root / path
    return path.resolve()


def _ensure_required_stages(
    by_canonical: Dict[str, Dict[str, Any]], required: Sequence[str]
) -> None:
    if not required:
        return
    missing = [stage for stage in required if stage not in by_canonical]
    if missing:
        raise RuntimeError(
            "Missing required stages in run_report.json: " + ", ".join(missing)
        )


def _index_manifest_stages(
    report: Dict[str, Any],
) -> Tuple[Dict[str, List[Tuple[str, Dict[str, Any]]]], Dict[str, Dict[str, Any]]]:
    manifest = report.get("manifest", {})
    if not isinstance(manifest, dict):
        return {}, {}
    stages = manifest.get("stages", {})
    if not isinstance(stages, dict):
        return {}, {}
    by_key: Dict[str, Dict[str, Any]] = {}
    by_canonical: Dict[str, List[Tuple[str, Dict[str, Any]]]] = {}
    for key, value in stages.items():
        if not isinstance(value, dict):
            continue
        key_str = str(key)
        by_key[key_str] = value
        canonical = _canonical_stage_name(key_str)
        if canonical:
            by_canonical.setdefault(canonical, []).append((key_str, value))
    return by_canonical, by_key


def _manifest_entries_for_stage(
    stage_record: Dict[str, Any],
    manifest_by_canonical: Dict[str, List[Tuple[str, Dict[str, Any]]]],
    manifest_by_key: Dict[str, Dict[str, Any]],
) -> List[Tuple[str, Dict[str, Any]]]:
    entries: List[Tuple[str, Dict[str, Any]]] = []
    seen = set()
    canonical = stage_record.get("canonical")
    if canonical:
        for key, value in manifest_by_canonical.get(canonical, []):
            entries.append((key, value))
            seen.add(key)
    stage_name = str(stage_record.get("stage_name", ""))
    if stage_name in manifest_by_key and stage_name not in seen:
        entries.append((stage_name, manifest_by_key[stage_name]))
    return entries


def _flatten_matching_fields(
    value: Any,
    prefix: str,
    terms: Sequence[str],
) -> List[Tuple[str, Any]]:
    fields: List[Tuple[str, Any]] = []
    if isinstance(value, dict):
        for key in sorted(value):
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            fields.extend(_flatten_matching_fields(value[key], child_prefix, terms))
        return fields
    if isinstance(value, list):
        for index, item in enumerate(value):
            child_prefix = f"{prefix}[{index}]"
            fields.extend(_flatten_matching_fields(item, child_prefix, terms))
        return fields
    if any(term in prefix.lower() for term in terms):
        fields.append((prefix, value))
    return fields


def _unique_pairs(pairs: Iterable[Tuple[str, Any]]) -> List[Tuple[str, Any]]:
    seen = set()
    unique: List[Tuple[str, Any]] = []
    for key, value in sorted(pairs, key=lambda item: item[0]):
        signature = (key, _json_scalar(value))
        if signature in seen:
            continue
        seen.add(signature)
        unique.append((key, value))
    return unique


def _scan_used_signals(value: Any, prefix: str) -> List[str]:
    signals: List[str] = []
    if isinstance(value, dict):
        for key in sorted(value):
            item = value[key]
            path = f"{prefix}.{key}" if prefix else str(key)
            key_l = str(key).lower()
            path_l = path.lower()
            if isinstance(item, (dict, list)):
                signals.extend(_scan_used_signals(item, path))
                continue
            item_str_l = str(item).strip().lower() if isinstance(item, str) else ""

            if key_l == "decision" and item_str_l == "resume":
                signals.append(f"{path}=resume")
                continue
            if key_l == "velocity_mode" and item_str_l in {
                "velocity_reset",
                "gro_velocity_restart",
            }:
                signals.append(f"{path}={item_str_l}")
                continue
            if isinstance(item, bool) and item:
                if "allow_" in path_l or key_l in {
                    "resume_allowed",
                    "resume_requested",
                    "allow_velocity_reset",
                    "allow_gro_velocity_restart",
                }:
                    continue
                if any(term in path_l for term in ("resume", "restart", "velocity_reset")):
                    signals.append(f"{path}=true")
                continue
            if isinstance(item, str) and item_str_l == "used":
                if any(term in path_l for term in ("resume", "restart", "velocity")):
                    signals.append(f"{path}=used")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            signals.extend(_scan_used_signals(item, f"{prefix}[{index}]"))
    return sorted(set(signals))


def _continuity_status(allowed: Any, used: bool) -> str:
    if used:
        return "used"
    if isinstance(allowed, bool):
        return "allowed" if allowed else "not allowed"
    return "unknown"


def _is_yes(value: Optional[str]) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() in TRUE_STRINGS


def _render_pairs_block(pairs: Sequence[Tuple[str, Any]]) -> List[str]:
    if not pairs:
        return ["(none)"]
    return [f"{key} = {_json_scalar(value)}" for key, value in pairs]


def _resolve_out_path(run_root: Path, out_raw: Optional[str]) -> Path:
    if not out_raw:
        return run_root / "computational_details.md"
    out_path = Path(out_raw)
    if not out_path.is_absolute():
        out_path = run_root / out_path
    return out_path.resolve()


def _collect_software_info(
    stage_records: Sequence[Dict[str, Any]],
) -> Tuple[str, str, str]:
    for record in stage_records:
        entry = record.get("entry", {})
        gromacs = entry.get("gromacs", {}) if isinstance(entry, dict) else {}
        if not isinstance(gromacs, dict):
            continue
        version = gromacs.get("version")
        gmx_bin = gromacs.get("gmx_bin")
        if version or gmx_bin:
            source = f"{record.get('stage_name')} (index {record.get('index')})"
            return str(version or "unknown"), str(gmx_bin or "unknown"), source
    return "unknown", "unknown", "unknown"


def _stage_sort_key(stage_record: Dict[str, Any]) -> Tuple[int, int]:
    canonical = stage_record.get("canonical")
    if canonical in EXPECTED_CANONICAL_STAGES:
        return (0, EXPECTED_CANONICAL_STAGES.index(canonical))
    return (1, int(stage_record.get("index", 0)))


def _build_stage_detail(
    run_root: Path,
    stage_record: Dict[str, Any],
    report: Dict[str, Any],
    manifest_by_canonical: Dict[str, List[Tuple[str, Dict[str, Any]]]],
    manifest_by_key: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    entry = stage_record.get("entry", {})
    mdp = entry.get("mdp", {}) if isinstance(entry.get("mdp"), dict) else {}
    mdp_path_raw = mdp.get("path")
    mdp_path_resolved = _resolve_report_path(run_root, mdp_path_raw)
    mdp_sha_report = mdp.get("sha256")
    mdp_sha_computed: Optional[str] = None
    mdp_unavailable_reason: Optional[str] = None
    mdp_params: Optional[List[Tuple[str, str]]] = None
    mdp_param_lookup: Dict[str, str] = {}

    if mdp_path_resolved is None:
        mdp_unavailable_reason = "mdp.path missing in run_report.json"
    elif not mdp_path_resolved.exists():
        mdp_unavailable_reason = f"MDP file not found: {mdp_path_resolved}"
    else:
        mdp_sha_computed = _sha256_file(mdp_path_resolved)
        if isinstance(mdp_sha_report, str) and mdp_sha_report.strip():
            if mdp_sha_computed.lower() != mdp_sha_report.strip().lower():
                raise RuntimeError(
                    "MDP SHA256 mismatch for stage "
                    f"'{stage_record.get('stage_name')}' (index {stage_record.get('index')}): "
                    f"run_report.json mdp.sha256={mdp_sha_report}, "
                    f"computed sha256={mdp_sha_computed}, file={mdp_path_resolved}. "
                    "Refusing to export non-reproducible computational details."
                )
        mdp_param_lookup = _parse_mdp_file(mdp_path_resolved)
        mdp_params = _extract_absolute_mdp_values(mdp_param_lookup)

    policies = entry.get("policies", {}) if isinstance(entry, dict) else {}
    if not isinstance(policies, dict):
        policies = {}

    manifest_entries = _manifest_entries_for_stage(
        stage_record,
        manifest_by_canonical,
        manifest_by_key,
    )

    continuity_fields: List[Tuple[str, Any]] = []
    for key in sorted(policies):
        if key in KEY_POLICIES or any(term in key.lower() for term in CONTINUITY_TERMS):
            continuity_fields.append((f"policies.{key}", policies[key]))
    continuity_fields.extend(_flatten_matching_fields(entry, "stage", CONTINUITY_TERMS))
    for manifest_stage_name, manifest_stage_value in manifest_entries:
        continuity_fields.extend(
            _flatten_matching_fields(
                manifest_stage_value,
                f"manifest.stages.{manifest_stage_name}",
                CONTINUITY_TERMS,
            )
        )
    continuity_fields = _unique_pairs(continuity_fields)

    used_signals: List[str] = []
    used_signals.extend(_scan_used_signals(entry, "stage"))
    for manifest_stage_name, manifest_stage_value in manifest_entries:
        used_signals.extend(
            _scan_used_signals(
                manifest_stage_value,
                f"manifest.stages.{manifest_stage_name}",
            )
        )
    used_signals = sorted(set(used_signals))

    used_velocity_reset = any(
        "velocity_mode=velocity_reset" in item or "velocity_reset" in item and "=true" in item
        for item in used_signals
    )
    used_gro_velocity_restart = any(
        "velocity_mode=gro_velocity_restart" in item
        or "gro_velocity_restart" in item and "=true" in item
        for item in used_signals
    )
    used_resume = any(
        item.endswith("=resume")
        or "resume" in item and ("=true" in item or "=used" in item)
        for item in used_signals
    )

    continuity_status = [
        (
            "velocity_reset",
            _continuity_status(policies.get("allow_velocity_reset"), used_velocity_reset),
        ),
        (
            "gro_velocity_restart",
            _continuity_status(
                policies.get("allow_gro_velocity_restart"),
                used_gro_velocity_restart,
            ),
        ),
        ("resume", _continuity_status(policies.get("resume_allowed"), used_resume)),
    ]

    mdp_inference: List[str] = []
    if not mdp_param_lookup:
        mdp_inference.append("unknown / not inferable from MDP")
    else:
        if _is_yes(mdp_param_lookup.get("gen-vel")):
            mdp_inference.append("velocities generated by GROMACS (gen_vel=yes)")
        if _is_yes(mdp_param_lookup.get("continuation")):
            mdp_inference.append("continuation=yes")
        if not mdp_inference:
            mdp_inference.append("unknown / not inferable from MDP")

    patch = mdp.get("patch", {}) if isinstance(mdp.get("patch"), dict) else {}
    delta = patch.get("patched", {}) if isinstance(patch.get("patched"), dict) else {}

    return {
        "stage_record": stage_record,
        "manifest_entries": manifest_entries,
        "mdp_path_raw": mdp_path_raw,
        "mdp_path_resolved": mdp_path_resolved,
        "mdp_sha_report": mdp_sha_report,
        "mdp_sha_computed": mdp_sha_computed,
        "mdp_unavailable_reason": mdp_unavailable_reason,
        "mdp_params": mdp_params,
        "policies": policies,
        "continuity_fields": continuity_fields,
        "continuity_status": continuity_status,
        "mdp_inference": mdp_inference,
        "used_signals": used_signals,
        "delta_patch": delta,
    }


def _write_markdown(
    out_path: Path,
    run_root: Path,
    report_path: Path,
    report: Dict[str, Any],
    stage_records: Sequence[Dict[str, Any]],
    by_canonical: Dict[str, Dict[str, Any]],
    stage_details: Sequence[Dict[str, Any]],
) -> None:
    gmx_version, gmx_bin, gmx_source = _collect_software_info(stage_records)
    report_timestamp = report.get("generated_at_utc", "unknown")

    lines: List[str] = [
        "# Computational Details",
        "",
        "## Software",
        f"- GROMACS version: {gmx_version}",
        f"- GROMACS binary: {gmx_bin}",
        f"- Software source stage entry: {gmx_source}",
        "",
        "## Run Provenance",
        f"- run_root: {run_root}",
        f"- run_report: {report_path}",
        f"- run_id: {report.get('run_id', 'unknown')}",
        f"- system_id: {report.get('system_id', 'unknown')}",
        f"- run_report generated_at_utc: {report_timestamp}",
        "",
        "## Stage Summary",
        "- Expected ensemble sequence: EM -> NVT -> NPT -> MD",
    ]

    for stage in EXPECTED_CANONICAL_STAGES:
        label = CANONICAL_LABELS[stage]
        record = by_canonical.get(stage)
        if record is None:
            lines.append(f"- {label}: not run")
            continue
        lines.append(
            f"- {label}: found as '{record.get('stage_name')}' "
            f"(index {record.get('index')})"
        )

    extra_stages = [rec for rec in stage_records if rec.get("canonical") is None]
    if extra_stages:
        extra_text = ", ".join(
            f"{rec.get('stage_name')} (index {rec.get('index')})"
            for rec in extra_stages
        )
        lines.append(f"- Additional non-canonical stage entries: {extra_text}")

    for detail in stage_details:
        stage_record = detail["stage_record"]
        stage_name = stage_record.get("stage_name")
        canonical = stage_record.get("canonical")
        canonical_label = CANONICAL_LABELS.get(canonical, "UNKNOWN")
        stage_index = stage_record.get("index")
        manifest_entries = detail.get("manifest_entries", [])

        lines.extend(
            [
                "",
                f"## Stage {canonical_label}",
                f"- Report stage entry: {stage_name} (index {stage_index})",
                f"- Canonical stage: {canonical_label if canonical else 'unknown'}",
            ]
        )
        if manifest_entries:
            lines.append(
                "- Matched manifest stage entries: "
                + ", ".join(name for name, _ in manifest_entries)
            )
        else:
            lines.append("- Matched manifest stage entries: none")

        lines.extend(
            [
                f"- MDP path (run_report): {_json_scalar(detail.get('mdp_path_raw'))}",
                f"- MDP path (resolved): {_json_scalar(detail.get('mdp_path_resolved'))}",
                f"- MDP sha256 (run_report): {_json_scalar(detail.get('mdp_sha_report'))}",
                f"- MDP sha256 (computed): {_json_scalar(detail.get('mdp_sha_computed'))}",
            ]
        )

        mdp_unavailable_reason = detail.get("mdp_unavailable_reason")
        if mdp_unavailable_reason:
            lines.append(f"- MDP status: unavailable ({mdp_unavailable_reason})")
        else:
            lines.append("- MDP status: available")

        lines.append("- Final MDP parameters (absolute):")
        lines.append("```text")
        if detail.get("mdp_params") is None:
            lines.append("MDP unavailable: absolute parameters unknown")
        else:
            lines.extend(
                _render_pairs_block(
                    detail.get("mdp_params", [])
                )
            )
        lines.append("```")

        policy_pairs = sorted(detail.get("policies", {}).items(), key=lambda item: item[0])
        lines.append("- Policies:")
        lines.append("```text")
        lines.extend(_render_pairs_block(policy_pairs))
        lines.append("```")

        lines.append("- Continuity notes:")
        lines.append("```text")
        for key, status in detail.get("continuity_status", []):
            lines.append(f"{key}: {status}")
        for inference in detail.get("mdp_inference", []):
            lines.append(f"mdp_inference: {inference}")
        used_signals = detail.get("used_signals", [])
        if used_signals:
            for signal in used_signals:
                lines.append(f"event_evidence: {signal}")
        else:
            lines.append("event_evidence: none found; reporting allowed/not-allowed only")
        lines.append("```")

        lines.append("- Continuity-related report fields:")
        lines.append("```text")
        lines.extend(_render_pairs_block(detail.get("continuity_fields", [])))
        lines.append("```")

    delta_stages = []
    for detail in stage_details:
        delta = detail.get("delta_patch", {})
        if not isinstance(delta, dict) or not delta:
            continue
        stage_record = detail["stage_record"]
        stage_name = stage_record.get("stage_name")
        stage_label = CANONICAL_LABELS.get(stage_record.get("canonical"), "UNKNOWN")
        delta_items = ", ".join(
            f"{key}={_json_scalar(delta[key])}" for key in sorted(delta)
        )
        delta_stages.append((stage_label, stage_name, delta_items))

    if delta_stages:
        lines.extend(
            [
                "",
                "## MDP CLI Overrides (Delta)",
                "- Source: run_report.json -> stage.mdp.patch.patched",
                "- This section is delta-only and not the final absolute MDP state.",
            ]
        )
        for stage_label, stage_name, delta_items in delta_stages:
            lines.append(f"- {stage_label} ({stage_name}): {delta_items}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Export computational_details.md from run_report.json"
    )
    parser.add_argument(
        "--run-root",
        required=True,
        help="Path to OUT_GMX/<RUN_ID> directory",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output markdown path (default: <run_root>/computational_details.md)",
    )
    parser.add_argument(
        "--require-stages",
        default="",
        help=(
            "Comma-separated canonical stages required for strict publication mode. "
            "Allowed values: em,nvt,npt,md. Default: none."
        ),
    )
    args = parser.parse_args(argv)

    run_root = Path(args.run_root).resolve()
    report_path = run_root / "run_report.json"
    if not report_path.exists():
        raise RuntimeError(f"run_report.json not found: {report_path}")

    report = _load_json(report_path)
    stage_records, by_canonical = _collect_stage_records(report)
    required_stages = _parse_require_stages(args.require_stages)
    _ensure_required_stages(by_canonical, required_stages)

    manifest_by_canonical, manifest_by_key = _index_manifest_stages(report)
    ordered_stage_records = sorted(stage_records, key=_stage_sort_key)
    stage_details = [
        _build_stage_detail(
            run_root,
            stage_record,
            report,
            manifest_by_canonical,
            manifest_by_key,
        )
        for stage_record in ordered_stage_records
    ]

    out_path = _resolve_out_path(run_root, args.out)
    _write_markdown(
        out_path,
        run_root,
        report_path.resolve(),
        report,
        ordered_stage_records,
        by_canonical,
        stage_details,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
