"""
Crash Context Preservation.

On stage failure:
- Print failure-centric excerpts from relevant logs (to stderr)
- Print absolute workdir path
- Print exact command executed
- Write crash_report.txt atomically
"""

from dataclasses import dataclass
from datetime import datetime
import logging
import os
from pathlib import Path
import shlex
import sys
import tempfile
from typing import List, Optional, Sequence, Tuple, Union
from uuid import uuid4

logger = logging.getLogger(__name__)

_TAIL_READ_CAP_BYTES = 5 * 1024 * 1024
_SCORING_TAIL_BYTES = 256 * 1024
_TAIL_CHUNK_BYTES = 64 * 1024
_TAIL_PROGRESSIVE_STEPS = (200, 800, 2000)
_PRIMARY_WINDOW_BEFORE = 80
_PRIMARY_WINDOW_AFTER = 120
_TRAILING_TAIL_LINES = 80
_FALLBACK_LINES_NO_MARKER = 200

_FAILURE_MARKERS: Tuple[Tuple[str, str], ...] = (
    ("fatal error", "Fatal error"),
    ("error", "ERROR"),
    ("segmentation fault", "Segmentation fault"),
    ("constraint", "constraint"),
    ("lincs", "LINCS"),
    ("shake", "shake"),
    ("step", "step"),
    ("nan", "nan"),
    ("infinite", "infinite"),
    ("bond", "bond"),
    ("domain decomposition", "domain decomposition"),
    ("particle moved more than", "particle moved more than"),
    ("blowing up", "blowing up"),
    ("pressure", "pressure"),
    ("temperature", "temperature"),
)

_MARKER_SCORE_BY_LABEL = {
    "Fatal error": 260,
    "Segmentation fault": 220,
    "constraint": 180,
    "LINCS": 220,
    "shake": 160,
    "nan": 200,
    "infinite": 180,
    "domain decomposition": 160,
    "particle moved more than": 220,
    "blowing up": 180,
    "bond": 90,
    "ERROR": 40,
    "step": 35,
    "pressure": 30,
    "temperature": 30,
}

_PHYSICS_SIGNAL_MARKERS = {
    "Fatal error",
    "Segmentation fault",
    "constraint",
    "LINCS",
    "shake",
    "nan",
    "infinite",
    "domain decomposition",
    "particle moved more than",
    "blowing up",
    "bond",
}

_WRAPPER_NOISE_MARKERS: Tuple[Tuple[str, str], ...] = (
    ("mpirun has exited", "mpirun has exited"),
    ("mpi_abort", "MPI_Abort"),
    ("srun:", "srun:"),
    ("pmi", "PMI"),
)

_ARTIFACT_PATTERNS: Tuple[str, ...] = (
    "*.tpr",
    "*.cpt",
    "*.gro",
    "*.edr",
    "*.log",
    "*.mdp",
    "*.top",
    "*.itp",
)

_ENV_KEYS: Tuple[str, ...] = (
    "OMP_NUM_THREADS",
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "CUDA_VISIBLE_DEVICES",
)


@dataclass(frozen=True)
class SmartExcerpt:
    primary_title: str
    primary_text: str
    trailing_title: Optional[str]
    trailing_text: Optional[str]
    markers_found: List[str]
    tail_lines_examined: int
    truncated_to_cap: bool


@dataclass(frozen=True)
class LogEvaluation:
    path: Path
    filename_score: int
    content_score: int
    total_score: int
    mtime: float
    markers_found: List[str]
    wrapper_markers_found: List[str]
    wrapper_noise: bool
    reason: str


@dataclass(frozen=True)
class SelectedLog:
    role: str
    role_note: str
    evaluation: LogEvaluation


def _dedupe_keep_order(items: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def _eprint(message: str = "") -> None:
    print(message, file=sys.stderr)


def _read_tail_text(
    file_path: Path,
    max_bytes: int = _TAIL_READ_CAP_BYTES,
) -> Tuple[str, bool]:
    """
    Read file tail by scanning from the end in chunks.

    Returns:
        (decoded_tail_text, truncated_to_cap)
    """
    with open(file_path, "rb") as f:
        f.seek(0, 2)
        file_size = f.tell()
        if file_size == 0:
            return "", False

        pos = file_size
        bytes_read = 0
        chunks: List[bytes] = []

        while pos > 0 and bytes_read < max_bytes:
            read_size = min(_TAIL_CHUNK_BYTES, pos, max_bytes - bytes_read)
            pos -= read_size
            f.seek(pos)
            data = f.read(read_size)
            chunks.append(data)
            bytes_read += read_size

        truncated = pos > 0
        text = b"".join(reversed(chunks)).decode("utf-8", errors="replace")
        return text, truncated


def read_last_lines(file_path: Path, n: int = 50) -> str:
    """
    Read the last n lines of a file.
    """
    if not file_path.exists():
        return f"[File not found: {file_path}]"

    try:
        text, truncated = _read_tail_text(file_path, _TAIL_READ_CAP_BYTES)
        if text == "":
            return "[Empty file]"
        lines = text.splitlines()
        if len(lines) > n:
            lines = lines[-n:]
        if truncated:
            lines.append("[truncated after 5MB; showing last available lines]")
        return "\n".join(lines)
    except Exception as exc:
        return f"[Error reading file: {exc}]"


def _safe_mtime(path: Path) -> float:
    try:
        return path.stat().st_mtime
    except Exception:
        return 0.0


def _find_markers_in_text(text: str, markers: Sequence[Tuple[str, str]]) -> List[str]:
    lower = text.lower()
    found = [label for needle, label in markers if needle in lower]
    return _dedupe_keep_order(found)


def _find_marker_line_hits(lines: Sequence[str]) -> List[Tuple[int, List[str]]]:
    hits: List[Tuple[int, List[str]]] = []
    for idx, line in enumerate(lines):
        lower = line.lower()
        matched = [label for needle, label in _FAILURE_MARKERS if needle in lower]
        if matched:
            hits.append((idx, _dedupe_keep_order(matched)))
    return hits


def _build_smart_excerpt(file_path: Path) -> SmartExcerpt:
    if not file_path.exists():
        return SmartExcerpt(
            primary_title="Tail excerpt unavailable (file missing)",
            primary_text=f"[File not found: {file_path}]",
            trailing_title=None,
            trailing_text=None,
            markers_found=[],
            tail_lines_examined=0,
            truncated_to_cap=False,
        )
    if not file_path.is_file():
        return SmartExcerpt(
            primary_title="Tail excerpt unavailable (not a file)",
            primary_text=f"[Not a file: {file_path}]",
            trailing_title=None,
            trailing_text=None,
            markers_found=[],
            tail_lines_examined=0,
            truncated_to_cap=False,
        )

    try:
        tail_text, truncated = _read_tail_text(file_path, _TAIL_READ_CAP_BYTES)
    except Exception as exc:
        return SmartExcerpt(
            primary_title="Tail excerpt unavailable (read error)",
            primary_text=f"[Error reading file: {exc}]",
            trailing_title=None,
            trailing_text=None,
            markers_found=[],
            tail_lines_examined=0,
            truncated_to_cap=False,
        )

    if tail_text == "":
        return SmartExcerpt(
            primary_title="Tail excerpt unavailable (empty file)",
            primary_text="[Empty file]",
            trailing_title=None,
            trailing_text=None,
            markers_found=[],
            tail_lines_examined=0,
            truncated_to_cap=truncated,
        )

    all_lines = tail_text.splitlines()
    if not all_lines:
        return SmartExcerpt(
            primary_title="Tail excerpt unavailable (no decodable lines)",
            primary_text="[No decodable text]",
            trailing_title=None,
            trailing_text=None,
            markers_found=[],
            tail_lines_examined=0,
            truncated_to_cap=truncated,
        )

    hits: List[Tuple[int, List[str]]] = []
    lines_examined = min(_TAIL_PROGRESSIVE_STEPS[0], len(all_lines))
    search_region = all_lines[-lines_examined:]

    for step in _TAIL_PROGRESSIVE_STEPS:
        lines_examined = min(step, len(all_lines))
        search_region = all_lines[-lines_examined:]
        hits = _find_marker_line_hits(search_region)
        if hits:
            break
        if lines_examined >= len(all_lines):
            break

    if hits:
        marker_idx = hits[-1][0]
        start = max(0, marker_idx - _PRIMARY_WINDOW_BEFORE)
        end = min(len(search_region), marker_idx + _PRIMARY_WINDOW_AFTER + 1)
        primary_lines = list(search_region[start:end])
        markers_found = _dedupe_keep_order(
            [marker for _, markers in hits for marker in markers]
        )
        title = (
            "Primary failure excerpt "
            f"(markers: {', '.join(markers_found)}; "
            f"searched tail: {lines_examined} lines)"
        )
        trailing_lines = all_lines[-min(_TRAILING_TAIL_LINES, len(all_lines)):]
        if truncated:
            primary_lines.append("[tail read capped at last 5MB]")
            trailing_lines = list(trailing_lines) + ["[tail read capped at last 5MB]"]
        return SmartExcerpt(
            primary_title=title,
            primary_text="\n".join(primary_lines),
            trailing_title=f"Trailing tail (last {_TRAILING_TAIL_LINES} lines)",
            trailing_text="\n".join(trailing_lines),
            markers_found=markers_found,
            tail_lines_examined=lines_examined,
            truncated_to_cap=truncated,
        )

    fallback_n = min(_FALLBACK_LINES_NO_MARKER, len(all_lines))
    fallback_lines = list(all_lines[-fallback_n:])
    if truncated:
        fallback_lines.append("[tail read capped at last 5MB]")
    return SmartExcerpt(
        primary_title=(
            "Tail excerpt (no failure markers found in tail; "
            f"showing last {fallback_n} lines)"
        ),
        primary_text="\n".join(fallback_lines),
        trailing_title=None,
        trailing_text=None,
        markers_found=[],
        tail_lines_examined=lines_examined,
        truncated_to_cap=truncated,
    )


def _score_log(path: Path) -> int:
    """
    Filename-only hints (content scoring is applied separately).
    Higher score means better candidate for physics failure context.
    """
    name = path.name.lower()
    suffix = path.suffix.lower()
    score = 0

    if name in {"md.log", "npt.log", "nvt.log", "em.log"}:
        score += 160
    elif name == "grompp.log":
        score += 140
    elif name == "mdrun.log":
        score += 130
    elif suffix == ".log":
        score += 100

    if "stderr" in name:
        score -= 15
    if suffix == ".err":
        score -= 25
    if any(token in name for token in ("mpirun", "slurm", "srun", "pmi")):
        score -= 40

    return score


def _evaluate_log(path: Path) -> LogEvaluation:
    name = path.name.lower()
    filename_score = _score_log(path)
    reason_bits: List[str] = [f"filename_score={filename_score}"]

    markers_found: List[str] = []
    physics_markers_found: List[str] = []
    wrapper_markers_found: List[str] = []
    wrapper_noise = False
    content_score = 0

    try:
        tail_text, _ = _read_tail_text(path, _SCORING_TAIL_BYTES)
        markers_found = _find_markers_in_text(tail_text, _FAILURE_MARKERS)
        physics_markers_found = [
            marker for marker in markers_found if marker in _PHYSICS_SIGNAL_MARKERS
        ]
        wrapper_markers_found = _find_markers_in_text(tail_text, _WRAPPER_NOISE_MARKERS)

        if markers_found:
            marker_bonus = min(
                sum(_MARKER_SCORE_BY_LABEL.get(marker, 40) for marker in markers_found),
                520,
            )
            content_score += marker_bonus
            reason_bits.append(
                "content markers found: " + ", ".join(markers_found)
            )

        if wrapper_markers_found:
            if physics_markers_found:
                content_score -= 40
                reason_bits.append(
                    "wrapper markers also present: " + ", ".join(wrapper_markers_found)
                )
            else:
                wrapper_noise = True
                content_score -= 220
                reason_bits.append(
                    "wrapper-noise markers: " + ", ".join(wrapper_markers_found)
                )

        if tail_text.strip() == "":
            content_score -= 20
            reason_bits.append("tail is empty")
    except Exception as exc:
        content_score -= 80
        reason_bits.append(f"tail unreadable: {exc}")

    if (
        (("stderr" in name) or path.suffix.lower() == ".err")
        and wrapper_noise
        and not physics_markers_found
    ):
        content_score -= 60
        reason_bits.append("stderr/err demoted due to wrapper noise")

    total_score = filename_score + content_score
    reason_bits.append(f"content_score={content_score}")
    reason_bits.append(f"total_score={total_score}")

    return LogEvaluation(
        path=path,
        filename_score=filename_score,
        content_score=content_score,
        total_score=total_score,
        mtime=_safe_mtime(path),
        markers_found=markers_found,
        wrapper_markers_found=wrapper_markers_found,
        wrapper_noise=wrapper_noise,
        reason="; ".join(reason_bits),
    )


def _sort_log_candidates(log_files: List[Path]) -> List[Path]:
    return sorted(
        log_files,
        key=lambda p: (-_safe_mtime(p), p.as_posix().lower()),
    )


def _looks_like_wrapper_log(evaluation: LogEvaluation) -> bool:
    name = evaluation.path.name.lower()
    suffix = evaluation.path.suffix.lower()
    return (
        evaluation.wrapper_noise
        or "stderr" in name
        or suffix == ".err"
        or any(token in name for token in ("mpirun", "slurm", "srun", "pmi"))
    )


def _select_log_files(log_files: List[Path], max_logs: int = 2) -> List[SelectedLog]:
    evaluations = [_evaluate_log(path) for path in log_files if path.exists() and path.is_file()]
    if not evaluations:
        return []

    ranked = sorted(
        evaluations,
        key=lambda e: (-e.total_score, -e.mtime, e.path.as_posix().lower()),
    )

    primary = ranked[0]
    selected = [
        SelectedLog(
            role="Primary log",
            role_note="best physics failure signal",
            evaluation=primary,
        )
    ]

    if max_logs > 1:
        secondary_candidates = [e for e in ranked[1:] if _looks_like_wrapper_log(e)]
        if secondary_candidates:
            selected.append(
                SelectedLog(
                    role="Secondary log",
                    role_note="MPI/wrapper stderr (may be noisy)",
                    evaluation=secondary_candidates[0],
                )
            )

    return selected[:max_logs]


def _format_command(command: Union[str, List[str]]) -> str:
    if isinstance(command, list):
        try:
            return shlex.join(command)
        except AttributeError:
            return " ".join(command)
    return command


def _stage_dir_status_message(stage_dir: Path) -> Optional[str]:
    stage_abs = stage_dir.absolute()
    if not stage_dir.exists():
        return f"Stage directory is missing: {stage_abs}"
    if not stage_dir.is_dir():
        return f"Stage path exists but is not a directory: {stage_abs}"
    return None


def _resolve_selected_logs(
    stage_dir: Path,
    log_paths: Optional[List[Path]],
    max_logs: int = 2,
) -> List[SelectedLog]:
    if log_paths is None:
        candidates = find_log_files(stage_dir)
    else:
        candidates = [Path(p) for p in log_paths]
    return _select_log_files(candidates, max_logs=max_logs)


def find_log_files(stage_dir: Path) -> List[Path]:
    """
    Find relevant text log files in a stage directory.
    """
    if not stage_dir.exists() or not stage_dir.is_dir():
        return []

    allowed_suffixes = {
        ".log",
        ".out",
        ".err",
        ".txt",
        ".mdp",
        ".mdpout",
        ".top",
        ".itp",
    }
    excluded_suffixes = {".xtc", ".trr", ".tpr", ".cpt", ".gro", ".edr"}
    log_files: List[Path] = []
    seen = set()

    def _is_text_candidate(path: Path) -> bool:
        if not path.is_file():
            return False
        suffix = path.suffix.lower()
        if suffix in excluded_suffixes:
            return False
        return suffix in allowed_suffixes

    def _add_candidate(path: Path) -> None:
        key = path.as_posix()
        if key in seen:
            return
        seen.add(key)
        if _is_text_candidate(path):
            log_files.append(path)

    for path in stage_dir.iterdir():
        _add_candidate(path)
    logs_dir = stage_dir / "logs"
    if logs_dir.is_dir():
        for path in logs_dir.iterdir():
            _add_candidate(path)
    return _sort_log_candidates(log_files)


def _selection_summary_lines(selected: Sequence[SelectedLog]) -> List[str]:
    lines: List[str] = []
    if not selected:
        return lines
    lines.append("Selected Logs:")
    for selected_log in selected:
        eval_info = selected_log.evaluation
        marker_summary = (
            ", ".join(eval_info.markers_found)
            if eval_info.markers_found
            else "none"
        )
        lines.append(
            f"  - {selected_log.role}: {eval_info.path} "
            f"({selected_log.role_note})"
        )
        lines.append(f"    reason: {eval_info.reason}")
        lines.append(f"    content markers found: {marker_summary}")
    return lines


def _artifact_summary_lines(stage_dir: Path) -> List[str]:
    lines = ["Artifacts summary:"]
    if not stage_dir.exists() or not stage_dir.is_dir():
        lines.append("  [stage directory unavailable]")
        return lines

    logs_dir = stage_dir / "logs"
    for pattern in _ARTIFACT_PATTERNS:
        matches: List[Path] = []
        seen = set()
        for path in sorted(stage_dir.glob(pattern), key=lambda p: p.name.lower()):
            if path.is_file():
                key = path.as_posix()
                if key not in seen:
                    matches.append(path)
                    seen.add(key)
        if pattern == "*.log" and logs_dir.is_dir():
            for path in sorted(logs_dir.glob(pattern), key=lambda p: p.name.lower()):
                if path.is_file():
                    key = path.as_posix()
                    if key not in seen:
                        matches.append(path)
                        seen.add(key)

        lines.append(f"  - {pattern}: {'present' if matches else 'absent'} ({len(matches)} file(s))")
        for path in matches:
            try:
                size = path.stat().st_size
                rel = path.relative_to(stage_dir)
            except Exception:
                size = -1
                rel = path
            if size >= 0:
                lines.append(f"      {rel} ({size} bytes)")
            else:
                lines.append(f"      {rel} (size unavailable)")

    return lines


def _environment_snapshot_lines() -> List[str]:
    lines = ["Environment snapshot:"]
    env = os.environ

    for key in _ENV_KEYS:
        value = env.get(key, "<unset>")
        lines.append(f"  - {key}={value}")

    gmx_keys = sorted(k for k in env if k.startswith("GMX_"))
    if gmx_keys:
        for key in gmx_keys:
            lines.append(f"  - {key}={env.get(key, '')}")
    else:
        lines.append("  - GMX_* variables: <none set>")

    return lines


def _excerpt_block_lines(selected_log: SelectedLog) -> List[str]:
    eval_info = selected_log.evaluation
    excerpt = _build_smart_excerpt(eval_info.path)
    lines = [
        "-" * 60,
        f"{selected_log.role}: {eval_info.path}",
        f"Selection note: {selected_log.role_note}",
        f"Selection reason: {eval_info.reason}",
        "-" * 60,
        excerpt.primary_title,
        excerpt.primary_text,
    ]
    if excerpt.trailing_title and excerpt.trailing_text is not None:
        lines.extend(["", excerpt.trailing_title, excerpt.trailing_text])
    lines.append("")
    return lines


def print_crash_context(
    stage_dir: Path,
    command: Union[str, List[str]],
    log_paths: Optional[List[Path]] = None,
    exit_code: Optional[int] = None,
) -> None:
    """
    Print crash context to stderr.
    """
    cmd_str = _format_command(command)
    stage_dir_status = _stage_dir_status_message(stage_dir)

    _eprint("\n" + "=" * 60)
    _eprint("CRASH CONTEXT")
    _eprint("=" * 60)
    _eprint(f"\n[EXIT CODE] {exit_code if exit_code is not None else 'Unknown'}")
    _eprint(f"\n[WORKING DIRECTORY]\n{stage_dir.absolute()}")
    if stage_dir_status:
        _eprint(f"[STAGE DIRECTORY STATUS] {stage_dir_status}")
    _eprint(f"\n[COMMAND]\n{cmd_str}")

    selected = _resolve_selected_logs(stage_dir, log_paths, max_logs=2)
    if selected:
        _eprint("\n[SELECTED LOGS]")
        for line in _selection_summary_lines(selected):
            _eprint(line)
        for selected_log in selected:
            for line in _excerpt_block_lines(selected_log):
                _eprint(line)
    else:
        if stage_dir_status:
            _eprint(f"\n[No log files found: {stage_dir_status}]")
        else:
            _eprint("\n[No log files found]")

    _eprint("=" * 60)


def _atomic_write_report(
    candidate_dir: Path,
    report_name: str,
    report_text: str,
) -> Path:
    candidate_dir.mkdir(parents=True, exist_ok=True)
    report_path = candidate_dir / report_name
    tmp_path = candidate_dir / (
        f"{report_name}.tmp.{os.getpid()}.{uuid4().hex}"
    )

    try:
        with open(tmp_path, "w", encoding="utf-8") as f:
            f.write(report_text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, report_path)
        return report_path
    except Exception:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass
        raise


def write_crash_report(
    stage_dir: Path,
    command: Union[str, List[str]],
    log_paths: Optional[List[Path]] = None,
    exit_code: Optional[int] = None,
) -> Path:
    """
    Write crash_report.txt to stage directory.
    """
    cmd_str = _format_command(command)
    stage_dir_status = _stage_dir_status_message(stage_dir)
    selected = _resolve_selected_logs(stage_dir, log_paths, max_logs=2)

    lines = [
        "=" * 60,
        "CRASH REPORT",
        f"Generated: {datetime.now().isoformat()}",
        "=" * 60,
        "",
        f"Exit Code: {exit_code if exit_code is not None else 'Unknown'}",
        "",
        "Working Directory:",
        str(stage_dir.absolute()),
        "",
        "Stage Directory Status:",
        stage_dir_status if stage_dir_status else "OK",
        "",
        "Command Executed:",
        cmd_str,
        "",
    ]

    if selected:
        lines.extend(_selection_summary_lines(selected))
        lines.append("")
        for selected_log in selected:
            lines.extend(_excerpt_block_lines(selected_log))
    else:
        if stage_dir_status:
            lines.append(f"[No log files found: {stage_dir_status}]")
        else:
            lines.append("[No log files found]")
        lines.append("")

    lines.extend(_artifact_summary_lines(stage_dir))
    lines.append("")
    lines.extend(_environment_snapshot_lines())
    lines.extend(
        [
            "",
            "=" * 60,
            "To debug:",
            f"  cd {stage_dir.absolute()}",
            "  # Inspect logs and re-run command manually",
            "=" * 60,
        ]
    )
    report_text = "\n".join(lines)

    temp_dir = Path(tempfile.gettempdir())
    candidate_dirs = [stage_dir, stage_dir.parent, Path.cwd(), temp_dir]
    attempted = set()
    write_errors: List[str] = []

    for candidate_dir in candidate_dirs:
        key = candidate_dir.as_posix()
        if key in attempted:
            continue
        attempted.add(key)
        try:
            return _atomic_write_report(candidate_dir, "crash_report.txt", report_text)
        except Exception as exc:
            write_errors.append(f"{candidate_dir}: {exc}")

    fallback_name = f"crash_report.{os.getpid()}.{uuid4().hex}.txt"
    if write_errors:
        report_text = (
            report_text
            + "\n\n[Report write errors]\n"
            + "\n".join(write_errors)
            + "\n"
        )

    fallback_path = _atomic_write_report(temp_dir, fallback_name, report_text)
    logger.warning(
        "Failed to write crash_report.txt in standard locations; wrote fallback report to %s",
        fallback_path,
    )
    return fallback_path


def handle_stage_failure(
    stage_dir: Path,
    command: Union[str, List[str]],
    exit_code: Optional[int] = None,
    log_paths: Optional[List[Path]] = None,
) -> Path:
    """
    Complete crash context handling: print to stderr and write report.
    """
    print_crash_context(stage_dir, command, log_paths, exit_code)
    report_path = write_crash_report(stage_dir, command, log_paths, exit_code)
    _eprint(f"\n[CRASH REPORT] Written to: {report_path}")
    return report_path
