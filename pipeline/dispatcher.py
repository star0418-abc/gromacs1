"""
StageDispatcher: orchestrates pipeline stages, run-level reporting, QC gates, and provenance.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

from .provenance import collect_run_provenance, write_provenance_text
from .stages import STAGE_ORDER, get_stage

if TYPE_CHECKING:
    from .context import PipelineContext


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


class StageDispatcher:
    """
    Orchestrates execution of pipeline stages and writes run_report.json continuously.
    """

    _REVIEW_POINTS_ALL = ("htpolynet", "gmx_eq")
    _REPRO_STAGE_DIRS = (
        "03_gromacs/em",
        "03_gromacs/nvt",
        "03_gromacs/npt",
        "03_gromacs/md",
    )

    def __init__(self, ctx: "PipelineContext"):
        self.ctx = ctx
        self._run_root: Path = ctx.get_output_path()
        self._run_root.mkdir(parents=True, exist_ok=True)
        self._run_report_path: Path = self._run_root / "run_report.json"
        self._dispatcher_error_log_path: Path = self._run_root / "dispatcher_error.log"

        self._run_started_utc = _utc_now()
        self._run_started_mono = time.monotonic()
        self._run_status = "pending"

        self._stages_to_run = self.get_stages_to_run()
        self._stage_reports: List[Dict[str, Any]] = []
        self._warnings: List[str] = []
        self._qc_events: List[Dict[str, Any]] = []
        self._qc_review_points: List[Dict[str, Any]] = []
        self._missing_metrics_warned: set[str] = set()
        self._provenance: Dict[str, Any] = {}

        qc_policy = getattr(ctx, "qc_policy", None)
        if qc_policy is None:
            qc_policy = "warn" if ctx.stage == "all" else "off"
        self._qc_policy = qc_policy
        self.ctx.qc_policy = qc_policy
        self._qc_density_rel_tol = float(getattr(ctx, "qc_density_rel_tol", 0.05))
        raw_qc_stages = getattr(ctx, "qc_enable_after_stages", ["gmx_eq"])
        if isinstance(raw_qc_stages, str):
            parsed_stages = [tok.strip() for tok in raw_qc_stages.split(",") if tok.strip()]
        elif isinstance(raw_qc_stages, (list, tuple, set)):
            parsed_stages = [str(tok).strip() for tok in raw_qc_stages if str(tok).strip()]
        else:
            parsed_stages = ["gmx_eq"]
        self._qc_enable_after_stages = set(parsed_stages)

    def run(self) -> bool:
        """
        Run pipeline stages and keep run_report.json up to date.
        """
        print("\n" + "=" * 60)
        print("PIPELINE EXECUTION")
        print("=" * 60)
        print(f"  Run ID: {self.ctx.run_id}")
        print(f"  System: {self.ctx.system_id}")
        print(f"  Forcefield: {self.ctx.ff}")
        print(f"  Charge method: {self.ctx.charge}")
        print(f"  Stage: {self.ctx.stage}")
        print(f"  Resume: {self.ctx.resume}")
        print(f"  Force: {self.ctx.force}")
        print(f"  QC policy: {self._qc_policy}")

        self._run_status = "running"
        try:
            self._collect_provenance()
            self._write_run_report()

            for stage_name in self._stages_to_run:
                stage_ok = self._execute_stage(stage_name)
                self._write_run_report()
                if not stage_ok:
                    return False

            self._run_status = "completed"
            self._write_run_report()
            print("\n" + "=" * 60)
            print("  PIPELINE COMPLETE")
            print("=" * 60)
            return True
        except KeyboardInterrupt:
            self.ctx.run_interrupted = True
            self._run_status = "interrupted"
            summary = "KeyboardInterrupt in dispatcher."
            tb_text = traceback.format_exc()
            self._write_dispatcher_error_log(None, summary, tb_text)
            print("\n  [ABORT] Interrupted by user.")
            if self.ctx.verbose:
                print(tb_text)
            self._write_run_report()
            return False
        except Exception as exc:
            self._run_status = "failed"
            summary = f"Dispatcher error: {type(exc).__name__}: {exc}"
            tb_text = traceback.format_exc()
            self._write_dispatcher_error_log(None, summary, tb_text)
            print(f"\n  [ABORT] {summary}")
            if self.ctx.verbose:
                print(tb_text)
            self._write_run_report()
            return False

    def _collect_provenance(self) -> None:
        self._provenance = collect_run_provenance(self.ctx)
        provenance_txt = write_provenance_text(self._run_root, self._provenance)
        self._provenance["provenance_txt_path"] = str(provenance_txt)

        if self.ctx.manifest:
            tool_versions = self._provenance.get("tool_versions", {})
            for tool, version in tool_versions.items():
                self.ctx.manifest.set_tool_version(tool, str(version))

        print("\n  Tool versions:")
        tool_versions = self._provenance.get("tool_versions", {})
        for tool_name in ("packmol", "gromacs", "htpolynet"):
            version = tool_versions.get(tool_name, "unknown")
            print(f"    - {tool_name}: {version}")
        print(f"  Provenance: {provenance_txt}")

    def _execute_stage(self, stage_name: str) -> bool:
        stage_start_mono = time.monotonic()
        stage_started_utc = _utc_now()
        stage_report: Dict[str, Any] = {
            "stage": stage_name,
            "started_at_utc": stage_started_utc,
            "ended_at_utc": None,
            "duration_s": None,
            "status": "running",
            "error_summary": None,
            "traceback_path": None,
        }
        try:
            try:
                stage = get_stage(stage_name)
            except KeyError as exc:
                summary = str(exc)
                trace_path = self._write_dispatcher_error_log(stage_name, summary, traceback.format_exc())
                stage_report["status"] = "failed"
                stage_report["error_summary"] = summary
                stage_report["traceback_path"] = trace_path
                self._run_status = "failed"
                return False

            try:
                success = stage.execute(self.ctx)
            except KeyboardInterrupt:
                self.ctx.run_interrupted = True
                self._run_status = "interrupted"
                summary = f"Interrupted by user during stage '{stage_name}'."
                tb_text = traceback.format_exc()
                trace_path = self._write_dispatcher_error_log(stage_name, summary, tb_text)
                stage_report["status"] = "interrupted"
                stage_report["error_summary"] = summary
                stage_report["traceback_path"] = trace_path
                print(f"\n  [ABORT] {summary}")
                if self.ctx.verbose:
                    print(tb_text)
                return False
            except Exception as exc:
                self._run_status = "failed"
                summary = f"{type(exc).__name__}: {exc}"
                tb_text = traceback.format_exc()
                trace_path = self._write_dispatcher_error_log(stage_name, summary, tb_text)
                stage_report["status"] = "failed"
                stage_report["error_summary"] = summary
                stage_report["traceback_path"] = trace_path
                print(f"\n  [ABORT] Stage {stage_name} raised {summary}")
                if self.ctx.verbose:
                    print(tb_text)
                return False

            if not success:
                self._run_status = "failed"
                summary = f"Stage '{stage_name}' returned False."
                trace_path = self._write_dispatcher_error_log(stage_name, summary, None)
                stage_report["status"] = "failed"
                stage_report["error_summary"] = summary
                stage_report["traceback_path"] = trace_path
                print(f"\n  [ABORT] {summary}")
                return False

            manifest_status = self._manifest_stage_status(stage_name)
            stage_report["status"] = (
                manifest_status if manifest_status in {"completed", "skipped"} else "completed"
            )

            qc_stop = self._run_stage_qc(stage_name, stage_report)
            if self._should_emit_review_point(stage_name):
                self._emit_review_point(stage_name, stage_report.get("metrics_path"))

            if qc_stop:
                return False

            return True
        finally:
            stage_report["ended_at_utc"] = _utc_now()
            stage_report["duration_s"] = round(max(0.0, time.monotonic() - stage_start_mono), 6)
            self._stage_reports.append(stage_report)
        
    def _manifest_stage_status(self, stage_name: str) -> Optional[str]:
        if not self.ctx.manifest:
            return None
        stages = self.ctx.manifest.get("stages", {})
        if not isinstance(stages, dict):
            return None
        stage_entry = stages.get(stage_name)
        if isinstance(stage_entry, dict):
            status = stage_entry.get("status")
            if isinstance(status, str) and status.strip():
                return status.strip()
        return None

    def _stage_metrics_candidates(self, stage_name: str) -> List[Path]:
        run_root = self._run_root
        mapping = {
            "packmol": [run_root / "01_packmol" / "stage_metrics.json"],
            "htpolynet": [run_root / "02_htpolynet" / "stage_metrics.json"],
            "sanitizer": [run_root / "02_5_sanitizer" / "stage_metrics.json"],
            "gmx_em": [run_root / "03_gromacs" / "em" / "stage_metrics.json"],
            "gmx_eq": [
                run_root / "03_gromacs" / "gmx_eq" / "stage_metrics.json",
                run_root / "03_gromacs" / "npt" / "stage_metrics.json",
                run_root / "03_gromacs" / "nvt" / "stage_metrics.json",
            ],
            "gmx_prod": [run_root / "03_gromacs" / "md" / "stage_metrics.json"],
            "analysis": [run_root / "analysis" / "stage_metrics.json"],
        }
        return mapping.get(stage_name, [run_root / stage_name / "stage_metrics.json"])

    def _load_stage_metrics(self, stage_name: str) -> Tuple[Optional[Path], Optional[Dict[str, Any]], Optional[str]]:
        for candidate in self._stage_metrics_candidates(stage_name):
            if not candidate.exists():
                continue
            try:
                payload = json.loads(candidate.read_text(encoding="utf-8"))
                if isinstance(payload, dict):
                    return candidate, payload, None
                return candidate, None, "stage_metrics.json is not a JSON object."
            except Exception as exc:
                return candidate, None, f"Failed to parse stage metrics: {type(exc).__name__}: {exc}"
        return None, None, None

    def _float_or_none(self, value: Any) -> Optional[float]:
        if value is None:
            return None
        if isinstance(value, dict):
            value = value.get("avg")
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _evaluate_qc_metrics(
        self,
        stage_name: str,
        metrics_path: Path,
        metrics: Dict[str, Any],
    ) -> Dict[str, Any]:
        messages: List[str] = []
        checks: Dict[str, Any] = {}

        temperature_avg = self._float_or_none(metrics.get("temperature_K"))
        pressure_avg = self._float_or_none(metrics.get("pressure_bar"))
        density_g_avg = self._float_or_none(metrics.get("density_g_cm3"))
        density_kg_avg = self._float_or_none(metrics.get("density_kg_m3"))
        if density_g_avg is None and density_kg_avg is not None:
            density_g_avg = density_kg_avg / 1000.0

        expected_density = self._float_or_none(metrics.get("expected_density_g_cm3"))
        notes = metrics.get("notes")
        if not isinstance(notes, list):
            notes = []
        if any(str(note).strip() == "failed_to_extract_metrics" for note in notes):
            return {
                "stage": stage_name,
                "metrics_path": str(metrics_path),
                "policy": self._qc_policy,
                "decision": "skip",
                "messages": [],
                "checks": {},
                "notes": notes,
                "recorded_at_utc": _utc_now(),
            }

        if temperature_avg is not None:
            checks["temperature_K_avg"] = temperature_avg
            if math.isnan(temperature_avg) or temperature_avg < 50.0 or temperature_avg > 1000.0:
                messages.append(
                    f"Stage {stage_name}: average temperature is non-physical ({temperature_avg} K)."
                )

        if pressure_avg is not None:
            checks["pressure_bar_avg"] = pressure_avg
            if math.isnan(pressure_avg) or abs(pressure_avg) > 5000.0:
                messages.append(
                    f"Stage {stage_name}: average pressure is extreme ({pressure_avg} bar)."
                )

        if expected_density is not None and expected_density > 0 and density_g_avg is not None:
            rel_err = abs(density_g_avg - expected_density) / expected_density
            checks["density_rel_err"] = rel_err
            if rel_err > self._qc_density_rel_tol:
                messages.append(
                    "NPT density deviates by "
                    f"{rel_err * 100:.1f}% from expected; production sampling likely meaningless "
                    "under wrong free volume"
                )

        decision = "pass"
        if messages:
            decision = "error" if self._qc_policy == "error" else "warn"

        return {
            "stage": stage_name,
            "metrics_path": str(metrics_path),
            "policy": self._qc_policy,
            "decision": decision,
            "messages": messages,
            "checks": checks,
            "notes": notes,
            "recorded_at_utc": _utc_now(),
        }

    def _run_stage_qc(self, stage_name: str, stage_report: Dict[str, Any]) -> bool:
        metrics_path, metrics, metrics_error = self._load_stage_metrics(stage_name)
        if metrics_path is not None:
            stage_report["metrics_path"] = str(metrics_path)

        qc_enabled = self._qc_policy != "off" and stage_name in self._qc_enable_after_stages
        if not qc_enabled:
            return False

        if metrics_path is None:
            if self.ctx.verbose and stage_name not in self._missing_metrics_warned:
                msg = (
                    f"QC enabled after stage '{stage_name}', but no stage_metrics.json was found. "
                    "Skipping QC checks."
                )
                self._warnings.append(msg)
                self._missing_metrics_warned.add(stage_name)
                print(f"  [QC][WARN] {msg}")
            return False

        if metrics_error is not None:
            event = {
                "stage": stage_name,
                "metrics_path": str(metrics_path),
                "policy": self._qc_policy,
                "decision": "error" if self._qc_policy == "error" else "warn",
                "messages": [metrics_error],
                "checks": {},
                "notes": ["failed_to_parse_stage_metrics"],
                "recorded_at_utc": _utc_now(),
            }
        else:
            event = self._evaluate_qc_metrics(stage_name, metrics_path, metrics or {})

        self._qc_events.append(event)
        stage_report["qc_decision"] = event.get("decision")
        event_messages = event.get("messages", [])
        if event_messages:
            for msg in event_messages:
                self._warnings.append(msg)
                level = "ERROR" if event.get("decision") == "error" else "WARN"
                print(f"  [QC][{level}] {msg}")

        if event.get("decision") == "error":
            summary = f"QC policy=error blocked continuation after stage '{stage_name}'."
            self._run_status = "failed_qc"
            stage_report["status"] = "failed"
            stage_report["error_summary"] = summary
            stage_report["traceback_path"] = self._write_dispatcher_error_log(stage_name, summary, None)
            print(f"  [ABORT] {summary}")
            return True

        if (
            self.ctx.stage == "all"
            and bool(getattr(self.ctx, "all_qc_stop_on_warn", False))
            and event.get("decision") == "warn"
        ):
            summary = f"QC warning triggered early stop after stage '{stage_name}'."
            self._run_status = "stopped_on_qc_warning"
            stage_report["status"] = "failed"
            stage_report["error_summary"] = summary
            stage_report["traceback_path"] = self._write_dispatcher_error_log(stage_name, summary, None)
            print(f"  [ABORT] {summary}")
            return True

        return False

    def _should_emit_review_point(self, stage_name: str) -> bool:
        return (
            self.ctx.stage == "all"
            and self._qc_policy == "warn"
            and stage_name in self._REVIEW_POINTS_ALL
        )

    def _emit_review_point(self, stage_name: str, metrics_path: Optional[str]) -> None:
        review_paths: List[str] = []
        if stage_name == "htpolynet":
            review_paths.append(str(self._run_root / "02_htpolynet"))
            review_paths.append(str(self._run_root / "02_htpolynet" / "workdir"))
        elif stage_name == "gmx_eq":
            review_paths.append(str(self._run_root / "03_gromacs" / "nvt"))
            review_paths.append(str(self._run_root / "03_gromacs" / "npt"))
        else:
            review_paths.extend(str(path.parent) for path in self._stage_metrics_candidates(stage_name))

        print("\n" + "!" * 60)
        print(f"  [REVIEW POINT] Stage '{stage_name}' finished. Inspect outputs before continuing:")
        for path in review_paths:
            print(f"    - {path}")
        if metrics_path:
            print(f"    - metrics: {metrics_path}")
        print("!" * 60)

        self._qc_review_points.append(
            {
                "stage": stage_name,
                "recorded_at_utc": _utc_now(),
                "paths": review_paths,
                "metrics_path": metrics_path,
            }
        )

    def _collect_repro_stages(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
        stages_data: List[Dict[str, Any]] = []
        errors: List[Dict[str, str]] = []
        for stage_rel in self._REPRO_STAGE_DIRS:
            repro_path = self._run_root / stage_rel / "repro.json"
            if not repro_path.exists():
                continue
            try:
                stages_data.append(json.loads(repro_path.read_text(encoding="utf-8")))
            except Exception as exc:
                errors.append(
                    {
                        "path": str(repro_path),
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                )
        return stages_data, errors

    def _atomic_write_json(self, path: Path, data: Dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(data, indent=2, sort_keys=True)
        temp_path: Optional[Path] = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=str(path.parent),
                prefix=f".{path.name}.tmp.",
                suffix=".json",
                delete=False,
            ) as fh:
                fh.write(payload)
                fh.flush()
                os.fsync(fh.fileno())
                temp_path = Path(fh.name)
            os.replace(temp_path, path)
            temp_path = None
            try:
                flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
                dir_fd = os.open(str(path.parent), flags)
                try:
                    os.fsync(dir_fd)
                finally:
                    os.close(dir_fd)
            except OSError:
                pass
        finally:
            if temp_path is not None and temp_path.exists():
                try:
                    temp_path.unlink()
                except OSError:
                    pass

    def _write_dispatcher_error_log(
        self,
        stage_name: Optional[str],
        summary: str,
        tb_text: Optional[str],
    ) -> str:
        lines = [
            f"timestamp_utc: {_utc_now()}",
            f"run_id: {self.ctx.run_id}",
            f"system_id: {self.ctx.system_id}",
            f"stage: {stage_name or 'dispatcher'}",
            f"summary: {summary}",
            "",
            "traceback:",
            (tb_text or "<no traceback available; stage returned False>").rstrip(),
            "",
        ]
        self._dispatcher_error_log_path.write_text("\n".join(lines), encoding="utf-8")
        return str(self._dispatcher_error_log_path)

    def _write_run_report(self) -> None:
        ended_at = _utc_now()
        duration_s = round(max(0.0, time.monotonic() - self._run_started_mono), 6)

        repro_stages, repro_errors = self._collect_repro_stages()
        manifest_data = self.ctx.manifest.data if self.ctx.manifest else {}

        qc_summary = {
            "pass": sum(1 for event in self._qc_events if event.get("decision") == "pass"),
            "skip": sum(1 for event in self._qc_events if event.get("decision") == "skip"),
            "warn": sum(1 for event in self._qc_events if event.get("decision") == "warn"),
            "error": sum(1 for event in self._qc_events if event.get("decision") == "error"),
        }

        report = {
            "run_id": self.ctx.run_id,
            "system_id": self.ctx.system_id,
            "stage_selection": self.ctx.stage,
            "resume": self.ctx.resume,
            "force": self.ctx.force,
            "stages_to_run": self._stages_to_run,
            "run_status": self._run_status,
            "started_at_utc": self._run_started_utc,
            "ended_at_utc": ended_at,
            "duration_s": duration_s,
            "warnings": self._warnings,
            "qc": {
                "policy": self._qc_policy,
                "density_rel_tol": self._qc_density_rel_tol,
                "enable_after_stages": sorted(self._qc_enable_after_stages),
                "all_qc_stop_on_warn": bool(getattr(self.ctx, "all_qc_stop_on_warn", False)),
                "events": self._qc_events,
                "review_points": self._qc_review_points,
                "summary": qc_summary,
            },
            "provenance": self._provenance,
            "dispatcher_error_log": (
                str(self._dispatcher_error_log_path) if self._dispatcher_error_log_path.exists() else None
            ),
            "stage_reports": self._stage_reports,
            "manifest": manifest_data,
            "stages": repro_stages,
            "errors": repro_errors,
            "summary": {
                "decisions": manifest_data.get("stages", {}) if isinstance(manifest_data, dict) else {},
                "stage_count": len(self._stage_reports),
            },
        }

        self._atomic_write_json(self._run_report_path, report)

    def get_stages_to_run(self) -> List[str]:
        """Get list of stage names that will run in dispatcher order."""
        if self.ctx.stage == "all":
            return list(STAGE_ORDER)
        return [self.ctx.stage]
