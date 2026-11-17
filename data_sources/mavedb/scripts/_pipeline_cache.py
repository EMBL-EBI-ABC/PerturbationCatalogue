"""Durable, checked stage results independent of Nextflow task work directories.

Callers decide stage dependencies and whether execution is authorized. This module
only stores immutable completed results, rejects incompatible reuse, and tracks
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

MANIFEST_NAME = ".pipeline-result.json"
SCHEMA_VERSION = 1


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


class ResultStore:
    """Completed stage generations below an analysis output directory.

    ``lock`` is a nonblocking whole-run writer lock. Individual state operations
    use a different blocking lock, allowing independent worker processes to
    register completed batches without losing each other's updates.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.directory = self.root / ".pipeline"
        self.state_path = self.directory / "state.json"

    @contextmanager
    def lock(self) -> Iterator[None]:
        self.directory.mkdir(parents=True, exist_ok=True)
        with (self.directory / "run.lock").open("a+") as handle:
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
        self.directory.mkdir(parents=True, exist_ok=True)
        with (self.directory / "state.lock").open("a+") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    def _state(self) -> dict[str, Any]:
        if not self.state_path.exists():
            return {"schema_version": SCHEMA_VERSION, "results": {}, "pending": {}}
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

    def begin(self, stage: str, rerun: bool = False) -> str:
        """Return an unfinished attempt ID, creating a fresh one for reruns."""
        _stage_name(stage)
        with self._state_lock():
            state = self._state()
            attempt = state["pending"].get(stage)
            if attempt is not None and (
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
        """Read an unfinished attempt without starting one."""
        _stage_name(stage)
        with self._state_lock():
            attempt = self._state()["pending"].get(stage)
        if attempt is not None and (
            not isinstance(attempt, str) or not re.fullmatch(r"[a-f0-9]{32}", attempt)
        ):
            raise CacheError(f"Invalid unfinished attempt for {stage}")
        return attempt

    def record(self, stage: str) -> dict[str, Any] | None:
        """Read completion metadata cheaply; use get to verify output content."""
        _stage_name(stage)
        with self._state_lock():
            results = self._state()["results"]
            if stage not in results:
                return None
            entry = results[stage]
        if not isinstance(entry, dict):
            raise CacheError(f"Damaged saved {stage} record; use --rerun {stage}")
        relative = _safe_relative(entry.get("path"))
        if relative.parent != Path(".pipeline") / "results" / stage:
            raise CacheError(f"Invalid saved {stage} result location")
        manifest_path = self.root / relative / MANIFEST_NAME
        try:
            if file_identity(manifest_path) != entry.get("manifest"):
                raise CacheError(
                    f"Saved {stage} completion record differs; use --rerun {stage}"
                )
        except OSError as error:
            raise CacheError(
                f"Cannot read saved {stage} completion record; use --rerun {stage}"
            ) from error
        return _read_json(manifest_path)

    def get(
        self,
        stage: str,
        settings: dict[str, Any] | None = None,
        inputs: dict[str, Any] | None = None,
    ) -> Path | None:
        """Return checked output paths or None only when never registered."""
        _stage_name(stage)
        with self._state_lock():
            results = self._state()["results"]
            if stage not in results:
                return None
            entry = results[stage]
        advice = f"Restore the saved result or explicitly use --rerun {stage}"
        if not isinstance(entry, dict):
            raise CacheError(f"Damaged saved {stage} record. {advice}")
        relative = _safe_relative(entry.get("path"))
        expected_parent = Path(".pipeline") / "results" / stage
        if relative.parent != expected_parent:
            raise CacheError(f"Invalid saved {stage} result location. {advice}")
        result = self.root / relative
        manifest_path = result / MANIFEST_NAME
        try:
            if not result.is_dir() or result.is_symlink():
                raise CacheError(f"Saved {stage} result directory is missing. {advice}")
            if file_identity(manifest_path) != entry.get("manifest"):
                raise CacheError(f"Saved {stage} completion record differs. {advice}")
            manifest = _read_json(manifest_path)
            if (
                manifest.get("schema_version") != SCHEMA_VERSION
                or manifest.get("stage") != stage
                or manifest.get("generation") != relative.name
                or manifest.get("status") != "completed"
                or not isinstance(manifest.get("settings"), dict)
                or not isinstance(manifest.get("inputs"), dict)
                or not isinstance(manifest.get("outputs"), dict)
            ):
                raise CacheError(f"Invalid saved {stage} completion record. {advice}")
            if any(path.is_symlink() for path in result.rglob("*")):
                raise CacheError(
                    f"Saved {stage} result contains a symbolic link. {advice}"
                )
            actual = {
                path.relative_to(result).as_posix()
                for path in result.rglob("*")
                if path.is_file() and path != manifest_path
            }
            if actual != set(manifest["outputs"]):
                raise CacheError(f"Saved {stage} output inventory differs. {advice}")
            for name, identity in manifest["outputs"].items():
                output = result / _safe_relative(name)
                if output.is_symlink() or file_identity(output) != identity:
                    raise CacheError(f"Saved {stage} output differs: {name}. {advice}")
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
        """Copy a verified producer output and promote it only after completion.

        Callers validate the stage's scientific/file contract before committing.
        Source symlinks are dereferenced so results survive removal of task work.
        Previous generations are retained if copying or state promotion fails.
        """
        _stage_name(stage)
        source = Path(source)
        if not source.is_dir():
            raise CacheError(f"Stage {stage} output directory does not exist: {source}")
        # Check serializability before expensive copying or any persistent mutation.
        json.dumps({"settings": settings, "inputs": inputs}, allow_nan=False)
        if not isinstance(settings, dict) or not isinstance(inputs, dict):
            raise CacheError("Stage settings and input identities must be objects")
        with self._state_lock():
            state = self._state()
            generation = state["pending"].get(stage) or uuid.uuid4().hex
            if not isinstance(generation, str) or not re.fullmatch(
                r"[a-f0-9]{32}", generation
            ):
                raise CacheError(f"Invalid unfinished attempt for {stage}")
            parent = self.directory / "results" / stage
            parent.mkdir(parents=True, exist_ok=True)
            destination = parent / generation
            # A checkpoint may have survived a failed state promotion. Never
            # overwrite it; another immutable generation is cheap and unambiguous.
            if destination.exists():
                generation = uuid.uuid4().hex
                destination = parent / generation
            with tempfile.TemporaryDirectory(
                prefix=".staging-", dir=parent
            ) as temporary:
                staged = Path(temporary) / "result"
                shutil.copytree(source, staged)
                (staged / MANIFEST_NAME).unlink(missing_ok=True)
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
                    "created_at": datetime.now(timezone.utc).isoformat(),
                    "settings": settings,
                    "inputs": inputs,
                    "implementation": {"registry": file_identity(Path(__file__))},
                    "outputs": outputs,
                }
                _write_json(staged / MANIFEST_NAME, manifest)
                staged.rename(destination)
            state["results"][stage] = {
                "path": destination.relative_to(self.root).as_posix(),
                "manifest": file_identity(destination / MANIFEST_NAME),
            }
            state["pending"].pop(stage, None)
            _write_json(self.state_path, state)
        return destination
