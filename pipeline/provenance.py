"""
Run-level provenance collection for dispatcher reports.

Captures raw tool version probes and runtime environment decisions.
"""

from __future__ import annotations

import json
import os
import platform
import shlex
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, TYPE_CHECKING

from .gromacs_cmd import resolve_gmx_command
from .resource_scheduler import get_resource_summary

if TYPE_CHECKING:
    from .context import PipelineContext


PROBE_OUTPUT_CAP_CHARS = 4000


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _best_effort_fsync(fd: int) -> None:
    try:
        os.fsync(fd)
    except OSError:
        pass


def _best_effort_fsync_dir(directory: Path) -> None:
    dir_flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        dir_flags |= os.O_DIRECTORY
    try:
        dir_fd = os.open(str(directory), dir_flags)
    except OSError:
        return
    try:
        _best_effort_fsync(dir_fd)
    finally:
        try:
            os.close(dir_fd)
        except OSError:
            pass


def _atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.parent / f"{path.name}.tmp.{os.getpid()}"
    try:
        with open(tmp_path, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            _best_effort_fsync(handle.fileno())
        os.replace(tmp_path, path)
        _best_effort_fsync_dir(path.parent)
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass


def _cap_text(text: Optional[str], cap_chars: int = PROBE_OUTPUT_CAP_CHARS) -> tuple[str, bool, int]:
    content = str(text or "")
    original_len = len(content)
    if original_len <= cap_chars:
        return content, False, original_len
    return content[:cap_chars], True, original_len


def _set_capped_text_fields(payload: Dict[str, Any], field_name: str, text: Optional[str]) -> None:
    preview, was_truncated, original_len = _cap_text(text)
    payload[field_name] = preview
    payload[f"{field_name}_truncated"] = was_truncated
    payload[f"{field_name}_original_len"] = original_len


def _tokenize_command(command: Any) -> Dict[str, Any]:
    if command is None:
        return {"raw": "", "tokens": [], "tokenization_warning": None}
    if isinstance(command, (list, tuple)):
        tokens = [str(part) for part in command if part is not None and str(part).strip()]
        return {
            "raw": " ".join(tokens),
            "tokens": tokens,
            "tokenization_warning": None,
        }

    raw = str(command or "").strip()
    if not raw:
        return {"raw": "", "tokens": [], "tokenization_warning": None}
    try:
        return {
            "raw": raw,
            "tokens": shlex.split(raw, posix=True),
            "tokenization_warning": None,
        }
    except ValueError as exc:
        return {
            "raw": raw,
            "tokens": [raw],
            "tokenization_warning": (
                f"shlex.split(posix=True) failed: {exc}. "
                "Command was treated as a single token."
            ),
        }


def _resolve_executable_path(token: Optional[str]) -> Optional[str]:
    if not token:
        return None
    resolved = shutil.which(token)
    if resolved:
        return resolved
    candidate = Path(token).expanduser()
    if candidate.exists():
        try:
            return str(candidate.resolve())
        except Exception:
            return str(candidate)
    return None


def _run_probe(
    cmd: List[str],
    *,
    timeout_s: int = 8,
    input_text: Optional[str] = None,
    tokenization_warning: Optional[str] = None,
) -> Dict[str, Any]:
    original_tokens = [str(part) for part in cmd if part is not None and str(part).strip()]
    cmd_str = " ".join(original_tokens)
    executable = original_tokens[0] if original_tokens else None
    executable_path = _resolve_executable_path(executable)
    probe_payload: Dict[str, Any] = {
        "command": cmd_str,
        "command_tokens": list(original_tokens),
        "original_command_tokens": list(original_tokens),
        "executable": executable,
        "executable_path": executable_path,
        "returncode": None,
        "error": None,
        "captured_at_utc": _utc_now(),
    }
    if tokenization_warning:
        probe_payload["tokenization_warning"] = tokenization_warning

    if executable_path is None:
        probe_payload["status"] = "not_found"
        _set_capped_text_fields(probe_payload, "stdout", "")
        _set_capped_text_fields(probe_payload, "stderr", "")
        probe_payload["error"] = "executable_not_found"
        return probe_payload

    executed_tokens = list(original_tokens)
    if executed_tokens:
        executed_tokens[0] = executable_path
    probe_payload["command_tokens"] = list(executed_tokens)

    try:
        result = subprocess.run(
            executed_tokens,
            input=input_text,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
        probe_payload["status"] = "ok" if result.returncode == 0 else "error"
        probe_payload["returncode"] = result.returncode
        _set_capped_text_fields(probe_payload, "stdout", result.stdout)
        _set_capped_text_fields(probe_payload, "stderr", result.stderr)
        return probe_payload
    except Exception as exc:
        probe_payload["status"] = "error"
        probe_payload["error"] = f"{type(exc).__name__}: {exc}"
        _set_capped_text_fields(probe_payload, "stdout", "")
        _set_capped_text_fields(probe_payload, "stderr", "")
        return probe_payload


def _extract_version_summary(tool_name: str, probe: Dict[str, Any]) -> str:
    status = str(probe.get("status") or "")
    if status == "not_found":
        return "not found"
    if status not in {"ok", "error"}:
        return "detection failed"
    stdout = str(probe.get("stdout") or "")
    stderr = str(probe.get("stderr") or "")
    combined_lines = [line.strip() for line in f"{stdout}\n{stderr}".splitlines() if line.strip()]

    if tool_name == "gromacs":
        for line in combined_lines:
            if "GROMACS version" in line:
                return line.split(":", 1)[-1].strip() or line
        for line in combined_lines:
            if line.lower().startswith("gromacs"):
                return line
    elif tool_name == "packmol":
        for line in combined_lines:
            upper = line.upper()
            if "PACKMOL" in upper and "VERSION" in upper:
                return line
    elif tool_name == "htpolynet":
        if combined_lines:
            return combined_lines[0]

    if status == "error":
        return "detection failed"
    return "detected (version unknown)"


def _parse_gromacs_build_features(raw_text: str) -> Dict[str, Optional[str]]:
    features: Dict[str, Optional[str]] = {
        "precision": None,
        "mpi": None,
        "openmp": None,
        "simd": None,
        "gpu_support": None,
    }
    key_values: Dict[str, str] = {}
    for line in (raw_text or "").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key_norm = key.strip().lower()
        value_clean = value.strip()
        if key_norm and value_clean:
            key_values[key_norm] = value_clean

    for key, value in key_values.items():
        if features["precision"] is None and "precision" in key:
            features["precision"] = value
        if features["mpi"] is None and "mpi" in key:
            features["mpi"] = value
        if features["openmp"] is None and "openmp" in key:
            features["openmp"] = value
        if features["simd"] is None and "simd" in key:
            features["simd"] = value
        if features["gpu_support"] is None and any(tok in key for tok in ("gpu", "cuda", "hip", "sycl", "opencl")):
            features["gpu_support"] = value
    return features


def _collect_runtime_provenance() -> Dict[str, Any]:
    uname = platform.uname()
    return {
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "uname": {
            "system": uname.system,
            "node": uname.node,
            "release": uname.release,
            "version": uname.version,
            "machine": uname.machine,
            "processor": uname.processor,
        },
        "hostname": socket.gethostname(),
        "cwd": str(Path.cwd()),
    }


def _build_git_error(
    kind: str,
    message: str,
    *,
    cwd: Optional[Path],
    command: Optional[List[str]] = None,
    returncode: Optional[int] = None,
    stdout: Optional[str] = None,
    stderr: Optional[str] = None,
) -> Dict[str, Any]:
    error: Dict[str, Any] = {
        "kind": kind,
        "message": message,
        "cwd": str(cwd) if cwd is not None else None,
    }
    if command is not None:
        error["command"] = command
    if returncode is not None:
        error["returncode"] = returncode
    if stdout is not None:
        _set_capped_text_fields(error, "stdout", stdout)
    if stderr is not None:
        _set_capped_text_fields(error, "stderr", stderr)
    return error


def _run_git_command(repo_root: Path, args: List[str]) -> subprocess.CompletedProcess[str]:
    git_path = shutil.which("git")
    if not git_path:
        raise FileNotFoundError("git executable not found")
    return subprocess.run(
        [git_path, *args],
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        check=False,
    )


def _collect_code_provenance(project_root: Optional[Path]) -> Dict[str, Any]:
    code: Dict[str, Any] = {
        "git_commit": None,
        "git_branch": None,
        "git_dirty": None,
    }
    if project_root is None:
        code["error"] = _build_git_error(
            "project_root_unavailable",
            "Pipeline context did not provide project_root for git provenance.",
            cwd=None,
        )
        return code

    repo_root = Path(project_root)
    git_path = shutil.which("git")
    if not git_path:
        code["error"] = _build_git_error(
            "git_unavailable",
            "git executable not found on PATH.",
            cwd=repo_root,
        )
        return code

    repo_check = subprocess.run(
        [git_path, "rev-parse", "--is-inside-work-tree"],
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        check=False,
    )
    repo_check_stdout = (repo_check.stdout or "").strip().lower()
    if repo_check.returncode != 0 or repo_check_stdout != "true":
        code["error"] = _build_git_error(
            "not_a_git_repo",
            "Unable to confirm that project_root is inside a git work tree.",
            cwd=repo_root,
            command=[git_path, "rev-parse", "--is-inside-work-tree"],
            returncode=repo_check.returncode,
            stdout=repo_check.stdout,
            stderr=repo_check.stderr,
        )
        return code

    try:
        commit_result = _run_git_command(repo_root, ["rev-parse", "HEAD"])
        if commit_result.returncode != 0:
            code["error"] = _build_git_error(
                "git_command_failed",
                "Failed to resolve git commit.",
                cwd=repo_root,
                command=[git_path, "rev-parse", "HEAD"],
                returncode=commit_result.returncode,
                stdout=commit_result.stdout,
                stderr=commit_result.stderr,
            )
            return code

        branch_result = _run_git_command(repo_root, ["rev-parse", "--abbrev-ref", "HEAD"])
        if branch_result.returncode != 0:
            code["error"] = _build_git_error(
                "git_command_failed",
                "Failed to resolve git branch.",
                cwd=repo_root,
                command=[git_path, "rev-parse", "--abbrev-ref", "HEAD"],
                returncode=branch_result.returncode,
                stdout=branch_result.stdout,
                stderr=branch_result.stderr,
            )
            return code

        dirty_result = _run_git_command(repo_root, ["status", "--porcelain"])
        if dirty_result.returncode != 0:
            code["error"] = _build_git_error(
                "git_command_failed",
                "Failed to resolve git dirty state.",
                cwd=repo_root,
                command=[git_path, "status", "--porcelain"],
                returncode=dirty_result.returncode,
                stdout=dirty_result.stdout,
                stderr=dirty_result.stderr,
            )
            return code
    except Exception as exc:
        code["error"] = _build_git_error(
            "git_probe_error",
            f"{type(exc).__name__}: {exc}",
            cwd=repo_root,
        )
        return code

    code["git_commit"] = (commit_result.stdout or "").strip() or None
    code["git_branch"] = (branch_result.stdout or "").strip() or None
    code["git_dirty"] = bool((dirty_result.stdout or "").strip())
    return code


def _probe_packmol_version() -> Dict[str, Any]:
    primary = _run_probe(["packmol", "--version"], timeout_s=5)
    primary_summary = _extract_version_summary("packmol", primary)
    if primary.get("status") == "not_found" or primary_summary != "detection failed":
        return primary

    fallback = _run_probe(["packmol"], timeout_s=5, input_text="")
    fallback["fallback_from"] = {
        "command": primary.get("command"),
        "command_tokens": primary.get("command_tokens"),
        "original_command_tokens": primary.get("original_command_tokens"),
        "status": primary.get("status"),
        "returncode": primary.get("returncode"),
        "error": primary.get("error"),
        "stdout": primary.get("stdout"),
        "stdout_truncated": primary.get("stdout_truncated"),
        "stdout_original_len": primary.get("stdout_original_len"),
        "stderr": primary.get("stderr"),
        "stderr_truncated": primary.get("stderr_truncated"),
        "stderr_original_len": primary.get("stderr_original_len"),
    }
    return fallback


def collect_run_provenance(ctx: "PipelineContext") -> Dict[str, Any]:
    """Collect tool and environment provenance once at dispatcher startup."""
    gmx_cmd = resolve_gmx_command(ctx)
    gmx_command_info = _tokenize_command(gmx_cmd)
    gmx_tokens = list(gmx_command_info.get("tokens") or [])
    if not gmx_tokens:
        gmx_tokens = ["gmx"]
    gmx_frontend = gmx_tokens[0]

    gmx_version_probe = _run_probe(
        gmx_tokens + ["--version"],
        timeout_s=10,
        tokenization_warning=gmx_command_info.get("tokenization_warning"),
    )
    mdrun_version_probe = _run_probe(
        gmx_tokens + ["mdrun", "-version"],
        timeout_s=10,
        tokenization_warning=gmx_command_info.get("tokenization_warning"),
    )
    gmx_raw = "\n".join(
        [
            str(gmx_version_probe.get("stdout") or ""),
            str(gmx_version_probe.get("stderr") or ""),
            str(mdrun_version_probe.get("stdout") or ""),
            str(mdrun_version_probe.get("stderr") or ""),
        ]
    ).strip()
    gmx_features = _parse_gromacs_build_features(gmx_raw)

    packmol_probe = _probe_packmol_version()
    htpolynet_probe = _run_probe(["htpolynet", "--version"], timeout_s=5)

    tools = {
        "packmol": {
            "version_probe": packmol_probe,
            "version_summary": _extract_version_summary("packmol", packmol_probe),
        },
        "gromacs": {
            "gmx_cmd": gmx_cmd,
            "gmx_cmd_tokens": gmx_tokens,
            "gmx_cmd_token": gmx_frontend,
            "gmx_cmd_path": _resolve_executable_path(gmx_frontend),
            "gmx_cmd_tokenization_warning": gmx_command_info.get("tokenization_warning"),
            "version_probe": gmx_version_probe,
            "mdrun_version_probe": mdrun_version_probe,
            "version_summary": _extract_version_summary("gromacs", gmx_version_probe),
            "build_features": gmx_features,
        },
        "htpolynet": {
            "version_probe": htpolynet_probe,
            "version_summary": _extract_version_summary("htpolynet", htpolynet_probe),
        },
    }

    try:
        resource_pinning = get_resource_summary(ctx)
    except Exception as exc:
        resource_pinning = {"error": f"{type(exc).__name__}: {exc}"}

    env = {
        "GMX_CMD_env": os.environ.get("GMX_CMD"),
        "resolved_GMX_CMD": gmx_cmd,
        "resolved_GMX_CMD_tokens": gmx_tokens,
        "resolved_GMX_CMD_path": _resolve_executable_path(gmx_frontend),
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        "GMX_GPU_ID": os.environ.get("GMX_GPU_ID"),
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "ctx_resource_overrides": {
            "gmx_nt": getattr(ctx, "gmx_nt", None),
            "gmx_gpu_id": getattr(ctx, "gmx_gpu_id", None),
            "omp_num_threads": getattr(ctx, "omp_num_threads", None),
        },
        "resource_pinning": resource_pinning,
    }

    tool_versions = {
        "packmol": tools["packmol"]["version_summary"],
        "gromacs": tools["gromacs"]["version_summary"],
        "htpolynet": tools["htpolynet"]["version_summary"],
    }

    return {
        "captured_at_utc": _utc_now(),
        "tools": tools,
        "env": env,
        "runtime": _collect_runtime_provenance(),
        "code": _collect_code_provenance(getattr(ctx, "project_root", None)),
        "tool_versions": tool_versions,
    }


def _format_text_value(value: Any) -> str:
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _append_probe_lines(lines: List[str], label: str, probe: Dict[str, Any]) -> None:
    lines.append(f"{label}_command: {probe.get('command')}")
    lines.append(f"{label}_command_tokens: {_format_text_value(probe.get('command_tokens'))}")
    lines.append(
        f"{label}_original_command_tokens: {_format_text_value(probe.get('original_command_tokens'))}"
    )
    lines.append(f"{label}_status: {probe.get('status')}")
    lines.append(f"{label}_returncode: {probe.get('returncode')}")
    lines.append(f"{label}_executable: {probe.get('executable')}")
    lines.append(f"{label}_executable_path: {probe.get('executable_path')}")
    if probe.get("tokenization_warning"):
        lines.append(f"{label}_tokenization_warning: {probe.get('tokenization_warning')}")
    if probe.get("error"):
        lines.append(f"{label}_error: {probe.get('error')}")
    lines.append(f"{label}_stdout_truncated: {probe.get('stdout_truncated')}")
    lines.append(f"{label}_stdout_original_len: {probe.get('stdout_original_len')}")
    lines.append(f"{label}_stdout:")
    lines.append(str(probe.get("stdout") or "").rstrip())
    lines.append(f"{label}_stderr_truncated: {probe.get('stderr_truncated')}")
    lines.append(f"{label}_stderr_original_len: {probe.get('stderr_original_len')}")
    lines.append(f"{label}_stderr:")
    lines.append(str(probe.get("stderr") or "").rstrip())
    fallback_from = probe.get("fallback_from")
    if isinstance(fallback_from, dict):
        lines.append(f"{label}_fallback_from: {_format_text_value(fallback_from)}")


def write_provenance_text(run_root: Path, provenance: Dict[str, Any]) -> Path:
    """Write a human-readable provenance artifact with raw probe outputs."""
    run_root.mkdir(parents=True, exist_ok=True)
    out_path = run_root / "provenance.txt"

    tools = provenance.get("tools", {}) if isinstance(provenance, dict) else {}
    env = provenance.get("env", {}) if isinstance(provenance, dict) else {}
    runtime = provenance.get("runtime", {}) if isinstance(provenance, dict) else {}
    code = provenance.get("code", {}) if isinstance(provenance, dict) else {}

    lines: List[str] = []
    lines.append(f"captured_at_utc: {provenance.get('captured_at_utc')}")
    lines.append("")
    lines.append("[runtime]")
    for key in ["python_version", "platform", "hostname", "cwd"]:
        lines.append(f"{key}: {_format_text_value(runtime.get(key))}")
    lines.append(f"uname: {_format_text_value(runtime.get('uname'))}")
    lines.append("")
    lines.append("[code]")
    for key in ["git_commit", "git_branch", "git_dirty"]:
        lines.append(f"{key}: {_format_text_value(code.get(key))}")
    if "error" in code:
        lines.append(f"error: {_format_text_value(code.get('error'))}")
    lines.append("")
    lines.append("[environment]")
    for key in [
        "GMX_CMD_env",
        "resolved_GMX_CMD",
        "resolved_GMX_CMD_tokens",
        "resolved_GMX_CMD_path",
        "OMP_NUM_THREADS",
        "GMX_GPU_ID",
        "CUDA_VISIBLE_DEVICES",
    ]:
        lines.append(f"{key}: {_format_text_value(env.get(key))}")
    lines.append("ctx_resource_overrides:")
    lines.append(_format_text_value(env.get("ctx_resource_overrides")))
    lines.append("resource_pinning:")
    lines.append(_format_text_value(env.get("resource_pinning")))
    lines.append("")

    for tool_name in ("packmol", "gromacs", "htpolynet"):
        tool_payload = tools.get(tool_name, {})
        lines.append(f"[tool:{tool_name}]")
        lines.append(f"version_summary: {tool_payload.get('version_summary')}")
        if tool_name == "gromacs":
            lines.append(f"gmx_cmd: {tool_payload.get('gmx_cmd')}")
            lines.append(f"gmx_cmd_path: {tool_payload.get('gmx_cmd_path')}")
            lines.append(f"gmx_cmd_tokens: {_format_text_value(tool_payload.get('gmx_cmd_tokens'))}")
            lines.append(
                f"gmx_cmd_tokenization_warning: "
                f"{_format_text_value(tool_payload.get('gmx_cmd_tokenization_warning'))}"
            )
            lines.append(f"build_features: {_format_text_value(tool_payload.get('build_features'))}")
        version_probe = tool_payload.get("version_probe", {})
        _append_probe_lines(lines, "version_probe", version_probe)
        if tool_name == "gromacs":
            mdrun_probe = tool_payload.get("mdrun_version_probe", {})
            lines.append("")
            _append_probe_lines(lines, "mdrun_version_probe", mdrun_probe)
        lines.append("")

    _atomic_write_text(out_path, "\n".join(lines).rstrip() + "\n")
    return out_path
