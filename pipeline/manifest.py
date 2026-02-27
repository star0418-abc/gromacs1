"""
ManifestWriter: Records pipeline execution metadata.

This module focuses on three guarantees:
1) Crash-safe manifest writes (atomic replace + best-effort durability)
2) Backward-compatible schema migration with explicit composition semantics
3) JSON-safe serialization without implicit default=str coercion
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import shlex
import shutil
import subprocess
import uuid
from datetime import date, datetime, time, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Union

try:  # Optional dependency; manifest normalization supports numpy scalars/arrays when present.
    import numpy as np  # type: ignore
except Exception:  # pragma: no cover - optional dependency may be unavailable
    np = None


SCHEMA_VERSION = 2
MANIFEST_VERSION = "1.1"

COMMAND_PREVIEW_CAP_CHARS = 4000
CHARGE_AFFECTED_ATOMS_PREVIEW_CAP = 200
ENV_PREVIEW_MAX_ITEMS = 32

_JSON_SAFE_PRIMITIVES = (str, int, float, bool)


def _now_utc() -> datetime:
    return datetime.now(timezone.utc)


def _now_utc_iso() -> str:
    return _now_utc().isoformat()


def _utc_stamp_for_filename() -> str:
    return _now_utc().strftime("%Y%m%dT%H%M%SZ")


def _best_effort_fsync(fd: int) -> None:
    try:
        os.fsync(fd)
    except OSError:
        # Durability is best-effort on platforms/filesystems that support fsync.
        pass


def _best_effort_fsync_dir(directory: Path) -> None:
    try:
        dir_fd = os.open(str(directory), os.O_RDONLY)
    except OSError:
        return
    try:
        _best_effort_fsync(dir_fd)
    finally:
        try:
            os.close(dir_fd)
        except OSError:
            pass


def _atomic_write_json(path: Path, data: Dict[str, Any]) -> None:
    """
    Atomically write JSON to disk using temp-file + replace semantics.

    Write steps:
    1) Write JSON to temp file in same directory
    2) Flush + fsync temp file (best-effort)
    3) os.replace(temp, final) for atomic swap
    4) fsync parent directory (best-effort) to persist rename on POSIX
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp_path = path.parent / f"{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    payload = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)

    try:
        with open(tmp_path, "w", encoding="utf-8") as handle:
            handle.write(payload)
            handle.write("\n")
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


class ManifestWriter:
    """
    Manages the manifest.json file for a pipeline run.

    The manifest records all metadata about the run including:
    - Configuration options
    - Tool versions
    - Input file hashes for reproducibility
    - Commands executed
    - Stage completion status
    """

    def __init__(self, manifest_path: Union[str, Path]):
        """
        Initialize ManifestWriter.

        Args:
            manifest_path: Path to manifest.json file
        """
        self.manifest_path = Path(manifest_path)
        self._data: Dict[str, Any] = self._load_or_init()

    # ------------------------------------------------------------------
    # Load / init / migration
    # ------------------------------------------------------------------

    def _load_or_init(self) -> Dict[str, Any]:
        """Load existing manifest or initialize a new one with migration."""
        if not self.manifest_path.exists():
            return self._migrate_manifest(self._init_manifest())

        try:
            with open(self.manifest_path, "r", encoding="utf-8") as handle:
                loaded = json.load(handle)
            if not isinstance(loaded, dict):
                raise ValueError("manifest root must be a JSON object")
            return self._migrate_manifest(loaded)
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            backup_path, backup_error = self._backup_corrupt_manifest_file()
            recovered = self._init_manifest()
            recovered["recovered_from_corrupt_manifest"] = True
            recovered["corrupt_backup_path"] = str(backup_path) if backup_path else None
            recovered["recovered_at_utc"] = _now_utc_iso()
            recovered["corruption_error"] = f"{type(exc).__name__}: {exc}"
            if backup_error:
                recovered["corrupt_backup_error"] = backup_error

            migrated = self._migrate_manifest(recovered)
            _atomic_write_json(self.manifest_path, migrated)
            return migrated

    def _backup_corrupt_manifest_file(self) -> tuple[Optional[Path], Optional[str]]:
        """
        Rename a corrupt manifest for post-mortem audit.

        Returns:
            (backup_path, error_message)
        """
        if not self.manifest_path.exists():
            return None, "manifest file disappeared before backup"

        base_name = f"{self.manifest_path.name}.corrupt.{_utc_stamp_for_filename()}"
        candidate = self.manifest_path.with_name(base_name)
        counter = 1
        while candidate.exists():
            candidate = self.manifest_path.with_name(f"{base_name}.{counter}")
            counter += 1

        try:
            self.manifest_path.replace(candidate)
            return candidate, None
        except OSError as exc:
            return None, str(exc)

    def _init_manifest(self) -> Dict[str, Any]:
        """Initialize a fresh manifest structure."""
        now = _now_utc_iso()
        return {
            "version": MANIFEST_VERSION,
            "schema_version": SCHEMA_VERSION,
            "created_at": now,
            "created_at_utc": now,
            "updated_at": now,
            "updated_at_utc": now,
            "options": {},
            "tool_versions": {},
            "input_assets": [],
            "composition": {},
            "mdp_patches": {},
            "sanitizer_outputs": {},
            "commands": [],
            "charge_corrections": [],
            "stages": {},
        }

    def _migrate_manifest(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply lightweight schema migration while preserving backward compatibility.

        Migration notes:
        - Keeps old top-level keys
        - Promotes legacy composition fields to canonical immutable semantics
        - Synthesizes box vectors from legacy scalar box_size_nm when needed
        """
        migrated = dict(data) if isinstance(data, dict) else {}
        now = _now_utc_iso()

        if "version" not in migrated:
            migrated["version"] = MANIFEST_VERSION

        schema_version = migrated.get("schema_version")
        if not isinstance(schema_version, int) or schema_version < SCHEMA_VERSION:
            migrated["schema_version"] = SCHEMA_VERSION

        created_at = migrated.get("created_at") or migrated.get("created_at_utc") or now
        updated_at = migrated.get("updated_at") or migrated.get("updated_at_utc") or created_at
        migrated["created_at"] = self._normalize_timestamp_value(created_at)
        migrated["created_at_utc"] = migrated.get("created_at_utc") or migrated["created_at"]
        migrated["updated_at"] = self._normalize_timestamp_value(updated_at)
        migrated["updated_at_utc"] = migrated.get("updated_at_utc") or migrated["updated_at"]

        self._ensure_top_level_defaults(migrated)
        self._migrate_composition_block(migrated)

        return migrated

    def _ensure_top_level_defaults(self, data: Dict[str, Any]) -> None:
        defaults: Dict[str, Any] = {
            "options": {},
            "tool_versions": {},
            "input_assets": [],
            "composition": {},
            "mdp_patches": {},
            "sanitizer_outputs": {},
            "commands": [],
            "charge_corrections": [],
            "stages": {},
        }
        for key, default in defaults.items():
            current = data.get(key)
            if isinstance(default, dict):
                if not isinstance(current, dict):
                    data[key] = {}
            elif isinstance(default, list):
                if not isinstance(current, list):
                    data[key] = []

    def _migrate_composition_block(self, data: Dict[str, Any]) -> None:
        composition = data.get("composition")
        if not isinstance(composition, dict):
            composition = {}
            data["composition"] = composition

        initial_counts = None
        if isinstance(composition.get("initial_counts"), dict):
            initial_counts = self._normalize_counts(composition["initial_counts"], "composition.initial_counts")
        elif isinstance(composition.get("counts"), dict):
            initial_counts = self._normalize_counts(composition["counts"], "composition.counts")
        elif isinstance(composition.get("molecule_counts"), dict):
            initial_counts = self._normalize_counts(composition["molecule_counts"], "composition.molecule_counts")

        if initial_counts is not None:
            composition["initial_counts"] = dict(initial_counts)
            composition["counts"] = dict(initial_counts)
            composition["molecule_counts"] = dict(initial_counts)

        if isinstance(composition.get("post_counts"), dict):
            composition["post_counts"] = self._normalize_counts(composition["post_counts"], "composition.post_counts")

        if not isinstance(composition.get("history"), list):
            composition["history"] = []
        if not isinstance(composition.get("conflicts"), list):
            composition["conflicts"] = []

        # Legacy scalar migration: synthesize vectors when only box_size_nm is present.
        normalized_box = None
        if composition.get("box") is not None:
            try:
                normalized_box = self._normalize_box(composition["box"], field_name="composition.box")
            except (TypeError, ValueError):
                normalized_box = None

        if normalized_box is None:
            if composition.get("box_size_nm") is not None:
                normalized_box = self._normalize_box(
                    composition.get("box_size_nm"),
                    field_name="composition.box_size_nm",
                )
            else:
                top_box = data.get("box")
                if isinstance(top_box, dict) and top_box.get("packmol_box_vectors_nm") is not None:
                    normalized_box = self._normalize_box(
                        {
                            "box_vectors_nm": top_box.get("packmol_box_vectors_nm"),
                            "box_type": "triclinic" if top_box.get("packmol_box_is_triclinic") else "orthorhombic",
                        },
                        field_name="box.packmol_box_vectors_nm",
                    )

        if normalized_box is not None:
            composition["box"] = normalized_box
            if normalized_box.get("box_size_nm") is not None and composition.get("box_size_nm") is None:
                composition["box_size_nm"] = normalized_box["box_size_nm"]
            self._record_box_alias(normalized_box, stage="initial", target_data=data)

    # ------------------------------------------------------------------
    # Save / serialization
    # ------------------------------------------------------------------

    def save(self) -> None:
        """Save manifest to disk using atomic replace semantics."""
        now = _now_utc_iso()
        self._data.setdefault("created_at", now)
        self._data.setdefault("created_at_utc", self._data.get("created_at"))
        self._data["updated_at"] = now
        self._data["updated_at_utc"] = now
        self._data["schema_version"] = max(
            SCHEMA_VERSION,
            int(self._data.get("schema_version", SCHEMA_VERSION)) if str(self._data.get("schema_version", "")).isdigit() else SCHEMA_VERSION,
        )

        self._ensure_top_level_defaults(self._data)
        self._enforce_composition_invariants()
        self._apply_size_caps()

        serialization_warnings: List[str] = []
        serializable = self._to_json_safe(
            self._data,
            path="$",
            warnings=serialization_warnings,
            seen=set(),
        )
        if serialization_warnings:
            existing = serializable.get("serialization_warnings", [])
            if not isinstance(existing, list):
                existing = []
            existing.extend(serialization_warnings)
            # Keep deterministic and bounded warning list.
            serializable["serialization_warnings"] = sorted(set(existing))[-200:]

        _atomic_write_json(self.manifest_path, serializable)
        self._data = serializable

    def load(self) -> Dict[str, Any]:
        """Reload manifest from disk."""
        self._data = self._load_or_init()
        return self._data

    @property
    def data(self) -> Dict[str, Any]:
        """Get the manifest data dictionary."""
        return self._data

    # ------------------------------------------------------------------
    # Generic get/set
    # ------------------------------------------------------------------

    def set(self, key: str, value: Any) -> None:
        """Set an arbitrary key-value pair in the manifest."""
        self._data[key] = value
        self.save()

    def get(self, key: str, default: Any = None) -> Any:
        """Get a value from the manifest."""
        return self._data.get(key, default)

    # ------------------------------------------------------------------
    # Options
    # ------------------------------------------------------------------

    def set_options(
        self,
        ff: str,
        charge: str,
        system_id: str,
        run_id: str,
        stage: str,
        resume: bool,
        force: bool,
    ) -> None:
        """Record pipeline options."""
        self._data["options"] = {
            "forcefield": ff,
            "charge_method": charge,
            "system_id": system_id,
            "run_id": run_id,
            "stage": stage,
            "resume": resume,
            "force": force,
        }
        self.save()

    # ------------------------------------------------------------------
    # Tool versions
    # ------------------------------------------------------------------

    def set_tool_version(self, tool_name: str, version: str) -> None:
        """Record a tool version."""
        self._data.setdefault("tool_versions", {})[tool_name] = version
        self.save()

    def detect_tool_versions(self) -> Dict[str, str]:
        """Auto-detect versions of required tools."""
        versions: Dict[str, str] = {}

        # Packmol
        if shutil.which("packmol"):
            try:
                result = subprocess.run(
                    ["packmol"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                for line in result.stdout.split("\n"):
                    if "PACKMOL" in line.upper() and "VERSION" in line.upper():
                        versions["packmol"] = line.strip()
                        break
                else:
                    versions["packmol"] = "detected (version unknown)"
            except Exception:
                versions["packmol"] = "detection failed"
        else:
            versions["packmol"] = "not found"

        # GROMACS
        if shutil.which("gmx"):
            try:
                result = subprocess.run(
                    ["gmx", "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                for line in result.stdout.split("\n"):
                    if "GROMACS version" in line:
                        versions["gromacs"] = line.split(":")[-1].strip()
                        break
                else:
                    versions["gromacs"] = "detected (version unknown)"
            except Exception:
                versions["gromacs"] = "detection failed"
        else:
            versions["gromacs"] = "not found"

        # HTPolyNet
        if shutil.which("htpolynet"):
            try:
                result = subprocess.run(
                    ["htpolynet", "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                versions["htpolynet"] = result.stdout.strip() or "detected"
            except Exception:
                versions["htpolynet"] = "detection failed"
        else:
            versions["htpolynet"] = "not found"

        self._data.setdefault("tool_versions", {}).update(versions)
        self.save()
        return versions

    # ------------------------------------------------------------------
    # Input assets
    # ------------------------------------------------------------------

    def add_input_asset(
        self,
        path: Union[str, Path],
        content_hash: str,
        category: Optional[str] = None,
        fail_on_hash_mismatch: bool = True,
    ) -> None:
        """
        Record an input asset with hash-aware dedup behavior.

        Rules:
        - Same path + same SHA256: no-op
        - Same path + different SHA256: fail-fast by default for reproducibility
        - Optional revision mode when fail_on_hash_mismatch=False
        """
        path_str = str(path)
        now = _now_utc_iso()
        assets = self._data.setdefault("input_assets", [])
        if not isinstance(assets, list):
            assets = []
            self._data["input_assets"] = assets

        existing = None
        for item in assets:
            if isinstance(item, dict) and item.get("path") == path_str:
                existing = item
                break

        if existing is not None:
            existing_hash = existing.get("sha256")
            if existing_hash == content_hash:
                return

            if fail_on_hash_mismatch:
                raise ValueError(
                    "Input asset hash mismatch for existing path "
                    f"'{path_str}': old={existing_hash}, new={content_hash}. "
                    "Refusing to overwrite provenance entry."
                )

            revision = {
                "path": path_str,
                "sha256": content_hash,
                "category": category,
                "supersedes_sha256": existing_hash,
                "detected_at_utc": now,
                "added_at": now,
                "added_at_utc": now,
            }
            assets.append(revision)
            self.save()
            return

        asset = {
            "path": path_str,
            "sha256": content_hash,
            "category": category,
            "added_at": now,
            "added_at_utc": now,
        }
        assets.append(asset)
        self.save()

    # ------------------------------------------------------------------
    # Composition / box
    # ------------------------------------------------------------------

    def set_initial_composition(
        self,
        counts: Dict[str, int],
        box: Optional[Union[float, int, Sequence[Any], Dict[str, Any]]] = None,
        density: Optional[float] = None,
        source: str = "packmol",
        strict: bool = True,
    ) -> None:
        """
        Set immutable initial composition (first writer wins).

        If initial_counts already exists and differs:
        - strict=True: raise ValueError
        - strict=False: append to composition.conflicts and do not overwrite
        """
        composition = self._data.setdefault("composition", {})
        if not isinstance(composition, dict):
            composition = {}
            self._data["composition"] = composition

        normalized_counts = self._normalize_counts(counts, "composition.initial_counts")
        existing_initial_raw = composition.get("initial_counts")
        existing_initial = (
            self._normalize_counts(existing_initial_raw, "composition.initial_counts")
            if isinstance(existing_initial_raw, dict)
            else None
        )

        if existing_initial is not None and existing_initial != normalized_counts:
            conflict = {
                "stage": source,
                "reason": "initial_counts_conflict",
                "existing_counts": existing_initial,
                "attempted_counts": normalized_counts,
                "at_utc": _now_utc_iso(),
            }
            conflicts = composition.setdefault("conflicts", [])
            if not isinstance(conflicts, list):
                conflicts = []
                composition["conflicts"] = conflicts
            conflicts.append(conflict)

            if strict:
                raise ValueError(
                    "composition.initial_counts is immutable and already set to a different value. "
                    "Use set_post_composition() for post-crosslink counts or call "
                    "set_initial_composition(..., strict=False) to record conflict without overwrite."
                )

            self._append_composition_history(
                composition,
                stage=source,
                counts=normalized_counts,
                reason="initial_conflict_recorded_no_overwrite",
            )
            self.save()
            return

        composition["initial_counts"] = normalized_counts
        composition["counts"] = dict(normalized_counts)
        composition["molecule_counts"] = dict(normalized_counts)
        composition["initial_counts_source"] = source
        composition["initial_counts_set_at_utc"] = _now_utc_iso()

        if density is not None:
            composition["density_kg_m3"] = float(density)

        if box is not None:
            normalized_box = self._normalize_box(box, field_name="composition.box")
            composition["box"] = normalized_box
            if normalized_box.get("box_size_nm") is not None:
                composition["box_size_nm"] = normalized_box["box_size_nm"]
            self._record_box_alias(normalized_box, stage="initial")

        self._append_composition_history(
            composition,
            stage=source,
            counts=normalized_counts,
            reason="set_initial_composition",
        )
        self.save()

    def set_post_composition(
        self,
        post_counts: Dict[str, int],
        box: Optional[Union[float, int, Sequence[Any], Dict[str, Any]]] = None,
        source: str = "htpolynet",
    ) -> None:
        """Record post-crosslink composition while preserving immutable initial_counts."""
        composition = self._data.setdefault("composition", {})
        if not isinstance(composition, dict):
            composition = {}
            self._data["composition"] = composition

        normalized_post = self._normalize_counts(post_counts, "composition.post_counts")
        composition["post_counts"] = normalized_post
        composition["post_counts_source"] = source
        composition["post_counts_timestamp"] = _now_utc_iso()  # backward-compatible key
        composition["post_counts_timestamp_utc"] = composition["post_counts_timestamp"]

        if box is not None:
            normalized_box = self._normalize_box(box, field_name="composition.post_box")
            composition["post_box"] = normalized_box
            self._record_box_alias(normalized_box, stage="post")

        self._append_composition_history(
            composition,
            stage=source,
            counts=normalized_post,
            reason="set_post_composition",
        )
        self.save()

    def set_composition(
        self,
        counts: Dict[str, int],
        box_size_nm: Optional[float] = None,
        density_kg_m3: Optional[float] = None,
    ) -> None:
        """
        Backward-compatible composition setter.

        Behavior:
        - If initial_counts missing: treat as set_initial_composition(...)
        - If initial_counts exists: treat as set_post_composition(...), and add compatibility note
        """
        composition = self._data.setdefault("composition", {})
        if not isinstance(composition, dict):
            composition = {}
            self._data["composition"] = composition

        if "initial_counts" not in composition:
            self.set_initial_composition(
                counts=counts,
                box=box_size_nm,
                density=density_kg_m3,
                source="legacy_set_composition",
                strict=True,
            )
            return

        notes = composition.setdefault("compatibility_notes", [])
        if not isinstance(notes, list):
            notes = []
            composition["compatibility_notes"] = notes
        notes.append(
            {
                "at_utc": _now_utc_iso(),
                "note": (
                    "set_composition() called after initial_counts was set; interpreted as post composition. "
                    "Prefer set_post_composition() for explicit stage-aware semantics."
                ),
            }
        )

        self.set_post_composition(
            post_counts=counts,
            box=box_size_nm,
            source="legacy_set_composition",
        )

        if density_kg_m3 is not None:
            composition = self._data.setdefault("composition", {})
            if isinstance(composition, dict):
                composition["density_kg_m3"] = float(density_kg_m3)
            self.save()

    def _append_composition_history(
        self,
        composition: Dict[str, Any],
        stage: str,
        counts: Dict[str, int],
        reason: str,
    ) -> None:
        history = composition.setdefault("history", [])
        if not isinstance(history, list):
            history = []
            composition["history"] = history
        history.append(
            {
                "stage": stage,
                "counts": counts,
                "reason": reason,
                "at_utc": _now_utc_iso(),
            }
        )

    def _record_box_alias(
        self,
        normalized_box: Dict[str, Any],
        stage: str,
        target_data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Keep top-level box aliases for existing readers while storing canonical vectors.
        """
        data = target_data if target_data is not None else self._data
        box_data = data.setdefault("box", {})
        if not isinstance(box_data, dict):
            box_data = {}
            data["box"] = box_data

        vectors = normalized_box.get("box_vectors_nm")
        if vectors is None:
            return

        is_triclinic = normalized_box.get("box_type") == "triclinic"
        if stage == "initial":
            box_data["packmol_box_vectors_nm"] = vectors
            box_data["packmol_box_is_triclinic"] = is_triclinic
        elif stage == "post":
            box_data["post_htpolynet_box_vectors_nm"] = vectors
            box_data["post_htpolynet_box_is_triclinic"] = is_triclinic

    def _normalize_box(
        self,
        box: Union[float, int, Sequence[Any], Dict[str, Any]],
        field_name: str,
    ) -> Dict[str, Any]:
        """
        Normalize box representation.

        Accepted inputs:
        - scalar float/int: legacy cubic box size
        - [3] numeric: orthorhombic vectors
        - [9] numeric: triclinic matrix flattened row-major
        - [[3],[3],[3]] numeric matrix: triclinic vectors
        - {box_vectors_nm: ..., box_type: ...}
        """
        if box is None:
            raise ValueError(f"{field_name} cannot be None")

        if isinstance(box, Mapping):
            raw_vectors = box.get("box_vectors_nm")
            raw_type = box.get("box_type")

            if raw_vectors is None and box.get("box_size_nm") is not None:
                raw_vectors = box.get("box_size_nm")

            if raw_vectors is None:
                raise ValueError(
                    f"{field_name} mapping must include 'box_vectors_nm' or 'box_size_nm'"
                )

            normalized = self._normalize_box(raw_vectors, field_name=f"{field_name}.box_vectors_nm")
            if raw_type is not None:
                if raw_type not in {"orthorhombic", "triclinic"}:
                    raise ValueError(
                        f"{field_name}.box_type must be 'orthorhombic' or 'triclinic', got {raw_type!r}"
                    )
                normalized["box_type"] = raw_type
            return normalized

        if self._is_number(box):
            length = self._normalize_float(box, field_name=field_name, require_positive=True)
            return {
                "box_vectors_nm": [length, length, length],
                "box_type": "orthorhombic",
                "box_size_nm": length,
            }

        if isinstance(box, Sequence) and not isinstance(box, (str, bytes, bytearray)):
            seq = list(box)

            # 3x3 matrix
            if len(seq) == 3 and all(
                isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray))
                for row in seq
            ):
                matrix = [
                    [
                        self._normalize_float(value, field_name=f"{field_name}[{row_idx}][{col_idx}]", require_positive=False)
                        for col_idx, value in enumerate(list(row))
                    ]
                    for row_idx, row in enumerate(seq)
                ]
                if any(len(row) != 3 for row in matrix):
                    raise ValueError(f"{field_name} matrix must be 3x3")
                self._validate_nonzero_box_matrix(matrix, field_name)
                return {
                    "box_vectors_nm": matrix,
                    "box_type": "triclinic",
                }

            # Orthorhombic 3-vector
            if len(seq) == 3 and all(self._is_number(value) for value in seq):
                vectors = [
                    self._normalize_float(value, field_name=f"{field_name}[{idx}]", require_positive=True)
                    for idx, value in enumerate(seq)
                ]
                normalized: Dict[str, Any] = {
                    "box_vectors_nm": vectors,
                    "box_type": "orthorhombic",
                }
                if self._is_isotropic_vector(vectors):
                    normalized["box_size_nm"] = vectors[0]
                return normalized

            # Flattened triclinic length-9
            if len(seq) == 9 and all(self._is_number(value) for value in seq):
                values = [
                    self._normalize_float(value, field_name=f"{field_name}[{idx}]", require_positive=False)
                    for idx, value in enumerate(seq)
                ]
                matrix = [values[0:3], values[3:6], values[6:9]]
                self._validate_nonzero_box_matrix(matrix, field_name)
                return {
                    "box_vectors_nm": matrix,
                    "box_type": "triclinic",
                }

        raise ValueError(
            f"Unsupported box format for {field_name}. Expected scalar, [3], [9], 3x3 matrix, "
            "or mapping with box_vectors_nm/box_size_nm."
        )

    @staticmethod
    def _is_number(value: Any) -> bool:
        return isinstance(value, (int, float)) and not isinstance(value, bool)

    @staticmethod
    def _normalize_float(value: Any, field_name: str, require_positive: bool) -> float:
        if isinstance(value, bool):
            raise TypeError(f"{field_name} must be numeric, got bool")
        if not isinstance(value, (int, float)):
            raise TypeError(f"{field_name} must be numeric, got {type(value).__name__}")
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f"{field_name} must be finite, got {value!r}")
        if require_positive and numeric <= 0.0:
            raise ValueError(f"{field_name} must be > 0, got {numeric}")
        return numeric

    @staticmethod
    def _validate_nonzero_box_matrix(matrix: List[List[float]], field_name: str) -> None:
        for row_idx, row in enumerate(matrix):
            if not any(abs(component) > 0.0 for component in row):
                raise ValueError(f"{field_name} row {row_idx} cannot be all zeros")

    @staticmethod
    def _is_isotropic_vector(vectors: Sequence[float], tol: float = 1e-12) -> bool:
        if len(vectors) != 3:
            return False
        return abs(vectors[0] - vectors[1]) <= tol and abs(vectors[1] - vectors[2]) <= tol

    def _normalize_counts(self, counts: Mapping[str, Any], field_name: str) -> Dict[str, int]:
        if not isinstance(counts, Mapping):
            raise TypeError(f"{field_name} must be a mapping of molecule->count")

        normalized: Dict[str, int] = {}
        for molecule, raw_count in counts.items():
            mol_key = str(molecule)
            if not mol_key:
                raise ValueError(f"{field_name} contains empty molecule key")
            if isinstance(raw_count, bool):
                raise TypeError(f"{field_name}.{mol_key} must be integer-like, got bool")
            try:
                count_float = float(raw_count)
            except (TypeError, ValueError):
                raise TypeError(
                    f"{field_name}.{mol_key} must be integer-like, got {raw_count!r}"
                ) from None
            if not math.isfinite(count_float) or not count_float.is_integer():
                raise ValueError(
                    f"{field_name}.{mol_key} must be an integer value, got {raw_count!r}"
                )
            count_int = int(count_float)
            if count_int < 0:
                raise ValueError(f"{field_name}.{mol_key} must be >= 0, got {count_int}")
            normalized[mol_key] = count_int

        return dict(sorted(normalized.items(), key=lambda item: item[0]))

    def _enforce_composition_invariants(self) -> None:
        """
        Keep canonical composition semantics even if callers mutate _data directly.
        """
        composition = self._data.get("composition")
        if not isinstance(composition, dict):
            composition = {}
            self._data["composition"] = composition

        initial_counts = None
        if isinstance(composition.get("initial_counts"), dict):
            initial_counts = self._normalize_counts(composition["initial_counts"], "composition.initial_counts")
        elif isinstance(composition.get("counts"), dict):
            initial_counts = self._normalize_counts(composition["counts"], "composition.counts")
        elif isinstance(composition.get("molecule_counts"), dict):
            initial_counts = self._normalize_counts(composition["molecule_counts"], "composition.molecule_counts")

        if initial_counts is not None:
            composition["initial_counts"] = dict(initial_counts)
            composition["counts"] = dict(initial_counts)
            composition["molecule_counts"] = dict(initial_counts)

        if isinstance(composition.get("post_counts"), dict):
            composition["post_counts"] = self._normalize_counts(composition["post_counts"], "composition.post_counts")

        if not isinstance(composition.get("history"), list):
            composition["history"] = []
        if not isinstance(composition.get("conflicts"), list):
            composition["conflicts"] = []

        # Keep box migration/invariants active on every save.
        normalized_box = None
        if composition.get("box") is not None:
            normalized_box = self._normalize_box(composition["box"], field_name="composition.box")
        elif composition.get("box_size_nm") is not None:
            normalized_box = self._normalize_box(
                composition["box_size_nm"],
                field_name="composition.box_size_nm",
            )

        if normalized_box is not None:
            composition["box"] = normalized_box
            if normalized_box.get("box_size_nm") is not None:
                composition["box_size_nm"] = normalized_box["box_size_nm"]
            self._record_box_alias(normalized_box, stage="initial")

    # ------------------------------------------------------------------
    # HTPolyNet metrics
    # ------------------------------------------------------------------

    def set_htpolynet_metrics(
        self,
        conversion: Optional[float] = None,
        gel_fraction: Optional[float] = None,
        qc_status: str = "unknown",
        warnings: Optional[List[str]] = None,
    ) -> None:
        """
        Record HTPolyNet crosslinking QC metrics.

        Args:
            conversion: Crosslinking conversion rate (0.0-1.0)
            gel_fraction: Gel fraction achieved
            qc_status: 'parsed', 'unknown', or 'failed'
            warnings: List of warning messages from HTPolyNet
        """
        if "htpolynet" not in self._data or not isinstance(self._data.get("htpolynet"), dict):
            self._data["htpolynet"] = {}

        recorded_at = _now_utc_iso()
        self._data["htpolynet"]["qc_metrics"] = {
            "conversion": conversion,
            "gel_fraction": gel_fraction,
            "qc_status": qc_status,
            "warnings": warnings or [],
            "recorded_at": recorded_at,
            "recorded_at_utc": recorded_at,
        }
        self.save()

    # ------------------------------------------------------------------
    # MDP patches
    # ------------------------------------------------------------------

    def set_mdp_patch(
        self,
        stage: str,
        original_values: Dict[str, Any],
        patched_values: Dict[str, Any],
    ) -> None:
        """Record MDP file patches for a stage."""
        patched_at = _now_utc_iso()
        self._data.setdefault("mdp_patches", {})[stage] = {
            "original": original_values,
            "patched": patched_values,
            "patched_at": patched_at,
            "patched_at_utc": patched_at,
        }
        self.save()

    # ------------------------------------------------------------------
    # Sanitizer outputs
    # ------------------------------------------------------------------

    def set_sanitizer_output(
        self,
        output_type: str,
        data: Any,
    ) -> None:
        """
        Record sanitizer outputs.

        Args:
            output_type: Type of output (combined_atomtypes, sanitized_itps, etc.)
            data: Output data
        """
        generated_at = _now_utc_iso()
        self._data.setdefault("sanitizer_outputs", {})[output_type] = {
            "data": data,
            "generated_at": generated_at,
            "generated_at_utc": generated_at,
        }
        self.save()

    # ------------------------------------------------------------------
    # Commands
    # ------------------------------------------------------------------

    def log_command(
        self,
        command: Union[str, Sequence[Any]],
        stage: str,
        exit_code: Optional[int] = None,
        stdout: Optional[str] = None,
        stderr: Optional[str] = None,
        returncode: Optional[int] = None,
        duration_s: Optional[float] = None,
        cwd: Optional[Union[str, Path]] = None,
        env: Optional[Mapping[str, Any]] = None,
        argv: Optional[Sequence[Any]] = None,
    ) -> None:
        """
        Log a command that was executed.

        Stores both machine-usable argv and human-readable command string.
        """
        now = _now_utc_iso()

        argv_list: List[str]
        if argv is not None:
            argv_list = [str(item) for item in argv]
            cmd_string = shlex.join(argv_list)
        elif isinstance(command, Sequence) and not isinstance(command, (str, bytes, bytearray)):
            argv_list = [str(item) for item in command]
            cmd_string = shlex.join(argv_list)
        else:
            cmd_string = str(command)
            try:
                argv_list = shlex.split(cmd_string)
            except ValueError:
                argv_list = [cmd_string]

        rc = returncode if returncode is not None else exit_code
        cmd_entry: Dict[str, Any] = {
            "argv": argv_list,
            "cmd": cmd_string,
            "command": cmd_string,  # backward-compatible alias
            "stage": stage,
            "returncode": rc,
            "exit_code": rc,  # backward-compatible alias
            "executed_at": now,
            "executed_at_utc": now,
        }

        if duration_s is not None:
            cmd_entry["duration_s"] = float(duration_s)
        if cwd is not None:
            cmd_entry["cwd"] = str(cwd)

        if env is not None:
            env_items = sorted((str(k), str(v)) for k, v in env.items())
            truncated = len(env_items) > ENV_PREVIEW_MAX_ITEMS
            preview_items = env_items[:ENV_PREVIEW_MAX_ITEMS]
            cmd_entry["env"] = {k: v for k, v in preview_items}
            if truncated:
                cmd_entry["env_truncated"] = True
                cmd_entry["env_total_items"] = len(env_items)

        if stdout is not None:
            preview, was_truncated, original_length = self._cap_text(str(stdout), COMMAND_PREVIEW_CAP_CHARS)
            cmd_entry["stdout_preview"] = preview
            if was_truncated:
                cmd_entry["stdout_preview_truncated"] = True
                cmd_entry["stdout_preview_original_len"] = original_length

        if stderr is not None:
            preview, was_truncated, original_length = self._cap_text(str(stderr), COMMAND_PREVIEW_CAP_CHARS)
            cmd_entry["stderr_preview"] = preview
            if was_truncated:
                cmd_entry["stderr_preview_truncated"] = True
                cmd_entry["stderr_preview_original_len"] = original_length

        self._data.setdefault("commands", []).append(cmd_entry)
        self.save()

    # ------------------------------------------------------------------
    # Charge corrections
    # ------------------------------------------------------------------

    def add_charge_correction(
        self,
        molecule: str,
        original_charge: float,
        corrected_charge: float,
        method: str,
        target_file: Optional[Union[str, Path]] = None,
        target_section: Optional[str] = None,
        affected_atoms: Optional[List[Dict[str, Any]]] = None,
        distribution_rule: Optional[str] = None,
        total_delta: Optional[float] = None,
        validation: Optional[Dict[str, Any]] = None,
        **extra_fields: Any,
    ) -> None:
        """
        Record a charge correction with optional atom-level audit detail.

        Heavy atom-level data is capped to keep manifest sizes bounded.
        """
        applied_at = _now_utc_iso()
        original = float(original_charge)
        corrected = float(corrected_charge)
        expected_delta = corrected - original if total_delta is None else float(total_delta)

        correction: Dict[str, Any] = {
            "molecule": molecule,
            "original_charge": original,
            "corrected_charge": corrected,
            "method": method,
            "total_delta": expected_delta,
            "applied_at": applied_at,
            "applied_at_utc": applied_at,
        }

        if target_file is not None:
            correction["target_file"] = str(target_file)
        if target_section is not None:
            correction["target_section"] = target_section
        if distribution_rule is not None:
            correction["distribution_rule"] = distribution_rule

        if affected_atoms is not None:
            normalized_atoms = self._normalize_affected_atoms(affected_atoms)
            atom_count = len(normalized_atoms)
            preview = normalized_atoms[:CHARGE_AFFECTED_ATOMS_PREVIEW_CAP]
            truncated = atom_count > CHARGE_AFFECTED_ATOMS_PREVIEW_CAP

            correction["affected_atoms"] = preview
            correction["affected_atoms_count"] = atom_count
            if truncated:
                correction["affected_atoms_truncated"] = True

            checksum_src = json.dumps(normalized_atoms, sort_keys=True, separators=(",", ":"))
            correction["affected_atoms_checksum_sha256"] = hashlib.sha256(
                checksum_src.encode("utf-8")
            ).hexdigest()

            sum_delta = 0.0
            for atom in normalized_atoms:
                atom_delta = atom.get("delta")
                if isinstance(atom_delta, (int, float)):
                    sum_delta += float(atom_delta)

            correction_validation = {
                "sum_affected_atom_deltas": sum_delta,
                "expected_total_delta": expected_delta,
                "corrected_minus_original": corrected - original,
                "sum_matches_total_delta": abs(sum_delta - expected_delta) <= 1e-6,
            }
            if isinstance(validation, dict):
                correction_validation.update(validation)
            correction["validation"] = correction_validation
        elif isinstance(validation, dict):
            correction["validation"] = validation

        if extra_fields:
            correction.update(extra_fields)

        self._data.setdefault("charge_corrections", []).append(correction)
        self.save()

    def _normalize_affected_atoms(self, atoms: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        normalized: List[Dict[str, Any]] = []
        for idx, atom in enumerate(atoms):
            if not isinstance(atom, dict):
                raise TypeError(f"affected_atoms[{idx}] must be a dict")

            entry: Dict[str, Any] = {}
            if "atom_index" in atom:
                entry["atom_index"] = int(atom["atom_index"])
            if "atom_name" in atom:
                entry["atom_name"] = str(atom["atom_name"])
            if "residue" in atom:
                entry["residue"] = str(atom["residue"])
            if "before" in atom:
                entry["before"] = float(atom["before"])
            if "after" in atom:
                entry["after"] = float(atom["after"])
            if "delta" in atom:
                entry["delta"] = float(atom["delta"])
            elif "before" in entry and "after" in entry:
                entry["delta"] = entry["after"] - entry["before"]

            # Keep any additional primitive fields for audit.
            for key, value in atom.items():
                if key in entry:
                    continue
                if value is None or isinstance(value, _JSON_SAFE_PRIMITIVES):
                    entry[str(key)] = value

            normalized.append(entry)

        return normalized

    # ------------------------------------------------------------------
    # Stage status
    # ------------------------------------------------------------------

    def set_stage_status(
        self,
        stage: str,
        status: str,
        message: Optional[str] = None,
    ) -> None:
        """
        Record stage completion status.

        Args:
            stage: Stage name
            status: Status (pending, running, completed, failed, skipped)
            message: Optional status message
        """
        now = _now_utc_iso()
        self._data.setdefault("stages", {})[stage] = {
            "status": status,
            "message": message,
            "updated_at": now,
            "updated_at_utc": now,
        }
        self.save()

    # ------------------------------------------------------------------
    # Helpers: size caps + serialization normalization
    # ------------------------------------------------------------------

    def _apply_size_caps(self) -> None:
        """
        Enforce size caps on hot-path verbose fields.

        This keeps frequent manifest writes reasonably efficient.
        """
        commands = self._data.get("commands")
        if isinstance(commands, list):
            for entry in commands:
                if not isinstance(entry, dict):
                    continue
                for key in ("stdout_preview", "stderr_preview"):
                    value = entry.get(key)
                    if isinstance(value, str):
                        preview, truncated, original_length = self._cap_text(
                            value,
                            COMMAND_PREVIEW_CAP_CHARS,
                        )
                        entry[key] = preview
                        if truncated:
                            entry[f"{key}_truncated"] = True
                            entry[f"{key}_original_len"] = original_length

        charge_corrections = self._data.get("charge_corrections")
        if isinstance(charge_corrections, list):
            for correction in charge_corrections:
                if not isinstance(correction, dict):
                    continue
                atoms = correction.get("affected_atoms")
                if isinstance(atoms, list) and len(atoms) > CHARGE_AFFECTED_ATOMS_PREVIEW_CAP:
                    correction["affected_atoms"] = atoms[:CHARGE_AFFECTED_ATOMS_PREVIEW_CAP]
                    correction["affected_atoms_truncated"] = True
                    correction["affected_atoms_count"] = len(atoms)

    @staticmethod
    def _cap_text(text: str, cap: int) -> tuple[str, bool, int]:
        if len(text) <= cap:
            return text, False, len(text)
        return text[:cap], True, len(text)

    def _to_json_safe(
        self,
        value: Any,
        path: str,
        warnings: List[str],
        seen: Set[int],
    ) -> Any:
        if value is None:
            return None
        if isinstance(value, bool):
            return value
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            if not math.isfinite(value):
                raise ValueError(f"Non-finite float at {path}: {value!r}")
            return value
        if isinstance(value, str):
            return value

        if isinstance(value, Path):
            return str(value)

        if isinstance(value, datetime):
            dt = value
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            return dt.astimezone(timezone.utc).isoformat()
        if isinstance(value, date):
            return value.isoformat()
        if isinstance(value, time):
            return value.isoformat()

        if np is not None:
            if isinstance(value, np.generic):
                return self._to_json_safe(value.item(), path, warnings, seen)
            if isinstance(value, np.ndarray):
                return self._to_json_safe(value.tolist(), path, warnings, seen)

        if isinstance(value, Mapping):
            obj_id = id(value)
            if obj_id in seen:
                raise ValueError(f"Cyclic reference detected at {path}")
            seen.add(obj_id)
            out: Dict[str, Any] = {}
            for raw_key, raw_val in value.items():
                key = raw_key if isinstance(raw_key, str) else str(raw_key)
                out[key] = self._to_json_safe(raw_val, f"{path}.{key}", warnings, seen)
            seen.remove(obj_id)
            return out

        if isinstance(value, (list, tuple, set, frozenset)):
            obj_id = id(value)
            if obj_id in seen:
                raise ValueError(f"Cyclic reference detected at {path}")
            seen.add(obj_id)
            out_list = [
                self._to_json_safe(item, f"{path}[{idx}]", warnings, seen)
                for idx, item in enumerate(value)
            ]
            seen.remove(obj_id)
            return out_list

        # Explicit fallback: no silent stringification of unknown objects.
        warning = (
            f"Unsupported type at {path}: {type(value).__module__}.{type(value).__qualname__}; "
            "stored repr fallback"
        )
        warnings.append(warning)
        return {
            "warning": "unsupported_type_repr_fallback",
            "type": f"{type(value).__module__}.{type(value).__qualname__}",
            "repr": repr(value),
        }

    @staticmethod
    def _normalize_timestamp_value(value: Any) -> str:
        if isinstance(value, datetime):
            dt = value if value.tzinfo is not None else value.replace(tzinfo=timezone.utc)
            return dt.astimezone(timezone.utc).isoformat()
        if isinstance(value, str):
            return value
        return _now_utc_iso()
