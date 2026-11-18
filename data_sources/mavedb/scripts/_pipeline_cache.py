"""Durable, checked stage results independent of Nextflow task work directories.

Callers decide stage dependencies and whether execution is authorized. This module
stores completed results in ordinary output folders, rejects incompatible reuse, and tracks
unfinished attempts. A missing result is different from a damaged saved result.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import shutil
import tempfile
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

MANIFEST_NAME = "pipeline_manifest.json"
SCHEMA_VERSION = 2
PUBLIC_STAGE_PATHS = {
    "ingest": "ingest/current",
    "taxonomy": "references/taxonomy",
    "process": "processed",
    "tables": "tables",
    "gmm": "gmm-fit",
    "merge": "gmm",
    "qc": "qc",
    "reports": "reports",
}


class CacheError(ValueError):
    """A saved result cannot safely satisfy a requested stage."""


def file_identity(path: Path) -> dict[str, int | str]:
    """Identify file content without using its location or modification time."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {"bytes": Path(path).stat().st_size, "sha256": digest.hexdigest()}


def input_identity(paths: dict[str, Path]) -> dict[str, dict[str, int | str]]:
    return {name: file_identity(path) for name, path in sorted(paths.items())}


def _write_json(path: Path, value: dict[str, Any]) -> None:
    """Replace a JSON record atomically on its destination filesystem."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CacheError(
            f"Cannot read saved pipeline record {path}: {error}"
        ) from error
    if not isinstance(value, dict):
        raise CacheError(f"Saved pipeline record is not an object: {path}")
    return value


def _stage_name(stage: str) -> str:
    if not isinstance(stage, str) or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.:-]*", stage
    ):
        raise CacheError(f"Invalid stage name: {stage!r}")
    return stage


def _safe_relative(value: Any) -> Path:
    if not isinstance(value, str) or not value:
        raise CacheError("Saved pipeline record has an invalid output path")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path == Path("."):
        raise CacheError(f"Saved pipeline path escapes its result: {value!r}")
    return path


def stage_path(stage: str) -> Path:
    """Return a stage's ordinary output folder, independent of attempt IDs."""
    _stage_name(stage)
    if stage in PUBLIC_STAGE_PATHS:
        return Path(PUBLIC_STAGE_PATHS[stage])
    if re.fullmatch(r"gmm-plan-[a-f0-9]{32}", stage):
        return Path("gmm-plan")
    batch = re.fullmatch(r"gmm-batch-[a-f0-9]{32}-([1-9][0-9]*)", stage)
    if batch:
        return Path("gmm-batches") / f"batch{batch[1]}"
    raise CacheError(f"Unknown pipeline stage: {stage!r}")


class ResultStore:
    """Checked stage outputs in real directories below an analysis output root.

    The small root-level index records completion and pending attempts. Each
    stage owns one stable directory with a completion manifest. ``lock`` guards
    a complete run; individual state operations use a separate worker lock.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.directory = self.root
        self.state_path = self.root / "pipeline_state.json"

    @contextmanager
    def lock(self) -> Iterator[None]:
        self.root.mkdir(parents=True, exist_ok=True)
        with (self.root / ".pipeline-run.lock").open("a+") as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise CacheError(
                    f"Another pipeline is writing to {self.root}; use a separate --outdir"
                ) from error
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    @contextmanager
    def _state_lock(self) -> Iterator[None]:
        self.root.mkdir(parents=True, exist_ok=True)
        with (self.root / ".pipeline-state.lock").open("a+") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    def _migration_error(self) -> CacheError:
        return CacheError(
            f"Legacy pipeline links or .pipeline storage found in {self.root}. "
            "Migrate these outputs before running the pipeline: "
            ".venv/bin/python data_sources/mavedb/scripts/migrate_pipeline_outputs.py "
            f"--outdir {str(self.root)!r}"
        )

    def _state(self) -> dict[str, Any]:
        legacy = self.root / ".pipeline"
        if not self.state_path.exists():
            if legacy.exists() or legacy.is_symlink():
                raise self._migration_error()
            return {"schema_version": SCHEMA_VERSION, "results": {}, "pending": {}}
        if self.state_path.is_symlink():
            raise self._migration_error()
        state = _read_json(self.state_path)
        if (
            state.get("schema_version") != SCHEMA_VERSION
            or not isinstance(state.get("results"), dict)
            or not isinstance(state.get("pending"), dict)
        ):
            raise CacheError(
                f"Unsupported or damaged pipeline state: {self.state_path}"
            )
        return state

    def _destination(self, stage: str) -> Path:
        relative = stage_path(stage)
        destination = self.root / relative
        current = destination
        while current != self.root:
            if current.is_symlink():
                raise self._migration_error()
            current = current.parent
        return destination

    def begin(self, stage: str, rerun: bool = False) -> str:
        """Return a persisted unfinished attempt ID; rerun starts a fresh one."""
        stage_path(stage)
        with self._state_lock():
            state = self._state()
            self._destination(stage)
            attempt = state["pending"].get(stage)
            if stage in state["pending"] and (
                not isinstance(attempt, str)
                or not re.fullmatch(r"[a-f0-9]{32}", attempt)
            ):
                raise CacheError(f"Invalid unfinished attempt for {stage}")
            if rerun or attempt is None:
                attempt = uuid.uuid4().hex
                state["pending"][stage] = attempt
                _write_json(self.state_path, state)
            return attempt

    def pending(self, stage: str) -> str | None:
        stage_path(stage)
        with self._state_lock():
            pending = self._state()["pending"]
            attempt = pending.get(stage)
            if stage in pending and (
                not isinstance(attempt, str)
                or not re.fullmatch(r"[a-f0-9]{32}", attempt)
            ):
                raise CacheError(f"Invalid unfinished attempt for {stage}")
        return attempt

    def _record(self, stage: str, state: dict[str, Any]) -> dict[str, Any] | None:
        destination = self._destination(stage)
        results = state["results"]
        if stage not in results:
            return None
        entry = results[stage]
        advice = f"Restore the saved result or explicitly use --rerun {stage}"
        if not isinstance(entry, dict):
            raise CacheError(f"Damaged saved {stage} record. {advice}")
        if not isinstance(entry.get("generation"), str) or not re.fullmatch(
            r"[a-f0-9]{32}", entry["generation"]
        ):
            raise CacheError(f"Invalid saved {stage} attempt identity. {advice}")
        if _safe_relative(entry.get("path")) != stage_path(stage):
            raise CacheError(f"Invalid saved {stage} result location. {advice}")
        if not destination.is_dir():
            raise CacheError(f"Saved {stage} result directory is missing. {advice}")
        manifest_path = destination / MANIFEST_NAME
        try:
            if file_identity(manifest_path) != entry.get("manifest"):
                raise CacheError(f"Saved {stage} completion record differs. {advice}")
        except OSError as error:
            raise CacheError(
                f"Cannot read saved {stage} completion record. {advice}"
            ) from error
        manifest = _read_json(manifest_path)
        if (
            manifest.get("schema_version") != SCHEMA_VERSION
            or manifest.get("stage") != stage
            or manifest.get("generation") != entry.get("generation")
            or manifest.get("result_path") != stage_path(stage).as_posix()
            or manifest.get("status") != "completed"
            or not isinstance(manifest.get("settings"), dict)
            or not isinstance(manifest.get("inputs"), dict)
            or not isinstance(manifest.get("outputs"), dict)
        ):
            raise CacheError(f"Invalid saved {stage} completion record. {advice}")
        return manifest

    def record(self, stage: str) -> dict[str, Any] | None:
        """Read completion metadata cheaply; use get to verify all outputs."""
        stage_path(stage)
        with self._state_lock():
            return self._record(stage, self._state())

    def get(
        self,
        stage: str,
        settings: dict[str, Any] | None = None,
        inputs: dict[str, Any] | None = None,
    ) -> Path | None:
        """Return checked outputs or None only when this stage was never saved."""
        stage_path(stage)
        with self._state_lock():
            manifest = self._record(stage, self._state())
            if manifest is None:
                return None
            result = self._destination(stage)
            advice = f"Restore the saved result or explicitly use --rerun {stage}"
            try:
                if any(path.is_symlink() for path in result.rglob("*")):
                    raise CacheError(
                        f"Saved {stage} result contains a symbolic link. {advice}"
                    )
                actual = {
                    path.relative_to(result).as_posix()
                    for path in result.rglob("*")
                    if path.is_file() and path != result / MANIFEST_NAME
                }
                if actual != set(manifest["outputs"]):
                    raise CacheError(
                        f"Saved {stage} output inventory differs. {advice}"
                    )
                for name, identity in manifest["outputs"].items():
                    if file_identity(result / _safe_relative(name)) != identity:
                        raise CacheError(
                            f"Saved {stage} output differs: {name}. {advice}"
                        )
            except OSError as error:
                raise CacheError(
                    f"Cannot validate saved {stage} outputs: {error}. {advice}"
                ) from error
            if settings is not None and manifest["settings"] != settings:
                raise CacheError(
                    f"Settings differ from saved {stage} result; use --rerun {stage} or a new --outdir"
                )
            if inputs is not None and manifest["inputs"] != inputs:
                raise CacheError(
                    f"Inputs differ from saved {stage} result; use --rerun {stage}"
                )
            return result

    def commit(
        self,
        stage: str,
        source: Path,
        settings: dict[str, Any],
        inputs: dict[str, Any],
    ) -> Path:
        """Replace one real output directory and its index, with rollback.

        Stage producers validate their own contracts before registration. The
        temporary backup is removed after success; no historical output copies
        or versioned directories are retained.
        """
        relative = stage_path(stage)
        source = Path(source)
        if not source.is_dir():
            raise CacheError(f"Stage {stage} output directory does not exist: {source}")
        if not isinstance(settings, dict) or not isinstance(inputs, dict):
            raise CacheError("Stage settings and input identities must be objects")
        json.dumps({"settings": settings, "inputs": inputs}, allow_nan=False)
        with self._state_lock():
            state = self._state()
            destination = self._destination(stage)
            if destination.exists() and not destination.is_dir():
                raise CacheError(f"Stage output path is not a directory: {destination}")
            if destination.parent.resolve().is_relative_to(source.resolve()):
                raise CacheError(
                    "Stage output source cannot contain its destination parent"
                )
            generation = state["pending"].get(stage, uuid.uuid4().hex)
            if not isinstance(generation, str) or not re.fullmatch(
                r"[a-f0-9]{32}", generation
            ):
                raise CacheError(f"Invalid unfinished attempt for {stage}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            workspace = Path(
                tempfile.mkdtemp(
                    prefix=f".{destination.name}-staging-", dir=destination.parent
                )
            )
            cleanup_workspace = True
            try:
                staged = workspace / "result"
                backup = workspace / "previous"
                previous_state = workspace / "pipeline_state.previous.json"
                had_state = self.state_path.exists()
                if had_state:
                    shutil.copyfile(self.state_path, previous_state)
                shutil.copytree(source, staged)
                for name in (MANIFEST_NAME, ".pipeline-result.json"):
                    (staged / name).unlink(missing_ok=True)
                outputs = {
                    path.relative_to(staged).as_posix(): file_identity(path)
                    for path in sorted(staged.rglob("*"))
                    if path.is_file()
                }
                if not outputs:
                    raise CacheError(f"Stage {stage} produced no output files")
                manifest = {
                    "schema_version": SCHEMA_VERSION,
                    "status": "completed",
                    "stage": stage,
                    "generation": generation,
                    "result_path": relative.as_posix(),
                    "created_at": datetime.now(timezone.utc).isoformat(),
                    "settings": settings,
                    "inputs": inputs,
                    "implementation": {"registry": file_identity(Path(__file__))},
                    "outputs": outputs,
                }
                _write_json(staged / MANIFEST_NAME, manifest)
                promoted = False
                try:
                    if destination.exists():
                        destination.rename(backup)
                    staged.rename(destination)
                    promoted = True
                    # Attempt-specific fitting keys share stable output slots.
                    # Remove superseded keys so no index entry points at another
                    # attempt's completion manifest.
                    state["results"] = {
                        name: entry
                        for name, entry in state["results"].items()
                        if name == stage
                        or not isinstance(entry, dict)
                        or entry.get("path") != relative.as_posix()
                    }
                    state["results"][stage] = {
                        "path": relative.as_posix(),
                        "generation": generation,
                        "manifest": file_identity(destination / MANIFEST_NAME),
                    }
                    state["pending"].pop(stage, None)
                    _write_json(self.state_path, state)
                except (OSError, ValueError, KeyboardInterrupt, SystemExit) as error:
                    try:
                        if promoted:
                            shutil.rmtree(destination)
                        if backup.exists():
                            backup.rename(destination)
                        if had_state:
                            os.replace(previous_state, self.state_path)
                        else:
                            self.state_path.unlink(missing_ok=True)
                    except (
                        OSError,
                        ValueError,
                        KeyboardInterrupt,
                        SystemExit,
                    ) as rollback_error:
                        cleanup_workspace = False
                        raise CacheError(
                            f"Could not register {stage} or restore its previous output: {rollback_error}. "
                            f"Recovery files are preserved at {workspace}. "
                            f"Restore any 'previous' directory to {destination} and any "
                            f"'pipeline_state.previous.json' file to {self.state_path} before retrying. "
                            f"Original registration failure: {error}"
                        ) from rollback_error
                    raise
            finally:
                if cleanup_workspace:
                    shutil.rmtree(workspace)
            return destination
