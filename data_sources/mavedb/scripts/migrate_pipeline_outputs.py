#!/usr/bin/env python3
"""Move legacy hidden pipeline results into verified ordinary output directories.

The default invocation is a read-only migration preview (apart from writer lock
files). ``--apply`` copies active results and resumable checkpoints without
removing ``.pipeline``. ``--remove-legacy`` is a separate, fully verified cleanup.
Run cleanup before running the pipeline again; a changed migrated state requires
another explicit review instead of deleting its historical source automatically.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import re
import shutil
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.scripts._pipeline_cache import CacheError, file_identity

LEGACY_MANIFEST = ".pipeline-result.json"
MANIFEST = "pipeline_manifest.json"
STATE = "pipeline_state.json"
RECEIPT = "pipeline_migration.json"
PUBLICATIONS = {
    "ingest": ["ingest/current"],
    "taxonomy": ["references/taxonomy"],
    "process": ["processed"],
    "tables": ["tables", "compatibility"],
    "gmm": ["gmm-fit"],
    "merge": ["gmm"],
    "qc": ["qc"],
    "reports": ["reports"],
}
ATTEMPT = re.compile(r"gmm-plan-([a-f0-9]{32})\Z")
BATCH = re.compile(r"gmm-batch-([a-f0-9]{32})-([1-9][0-9]*)\Z")


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CacheError(f"Cannot read migration record {path}: {error}") from error
    if not isinstance(value, dict):
        raise CacheError(f"Migration record must be a JSON object: {path}")
    return value


def _write(path: Path, value: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _relative(value: Any) -> Path:
    if not isinstance(value, str) or not value:
        raise CacheError("Invalid migration result path")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path == Path("."):
        raise CacheError(f"Migration result path escapes its output: {value!r}")
    return path


def _stage_path(stage: str) -> Path:
    # Import lazily so --help is usable while a checkout is being upgraded.
    from data_sources.mavedb.scripts._pipeline_cache import stage_path

    return Path(stage_path(stage))


def _present(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def _check_target(root: Path, relative: Path) -> None:
    for parent in relative.parents:
        if parent != Path(".") and (root / parent).is_symlink():
            raise CacheError(
                f"Publication parent is a symbolic link: {root / parent}. "
                "Materialize that parent before migration."
            )


def _inside(path: Path, root: Path) -> bool:
    try:
        return path.resolve(strict=True).is_relative_to(root.resolve(strict=True))
    except (OSError, RuntimeError):
        return False


def _inventory(
    path: Path,
    root: Path,
    *,
    exclude: set[str] | None = None,
    allow_links: bool = False,
) -> dict[str, Any]:
    """Hash dereferenced content, rejecting escaping links and directory cycles."""
    excluded = exclude or set()
    found: dict[str, Any] = {}

    def visit(source: Path, relative: Path, parents: set[Path]) -> None:
        if not _inside(source, root):
            raise CacheError(f"Missing or escaping migration source: {source}")
        if source.is_symlink() and not allow_links:
            raise CacheError(f"Expected ordinary files, found symbolic link: {source}")
        resolved = source.resolve(strict=True)
        if source.is_dir():
            if resolved in parents:
                raise CacheError(f"Cyclic migration source: {source}")
            for child in sorted(source.iterdir()):
                child_relative = relative / child.name
                if child_relative.as_posix() not in excluded:
                    visit(child, child_relative, parents | {resolved})
        elif source.is_file():
            found[relative.as_posix() if relative != Path(".") else ""] = file_identity(
                source
            )
        else:
            raise CacheError(f"Unsupported migration source file type: {source}")

    visit(path, Path("."), set())
    return found


def _legacy_result(
    root: Path, stage: str, entry: Any, *, outputs: bool = True
) -> tuple[Path, dict[str, Any]]:
    if not isinstance(entry, dict):
        raise CacheError(f"Invalid legacy result entry: {stage}")
    relative = _relative(entry.get("path"))
    if relative.parent != Path(".pipeline/results") / stage or not re.fullmatch(
        r"[a-f0-9]{32}", relative.name
    ):
        raise CacheError(f"Invalid legacy result location for {stage}: {relative}")
    source = root / relative
    if source.is_symlink() or not _inside(source, root / ".pipeline"):
        raise CacheError(f"Missing or escaping legacy result: {source}")
    manifest_path = source / LEGACY_MANIFEST
    if file_identity(manifest_path) != entry.get("manifest"):
        raise CacheError(f"Legacy completion record hash differs: {stage}")
    record = _json(manifest_path)
    if (
        record.get("schema_version") != 1
        or record.get("status") != "completed"
        or record.get("stage") != stage
        or record.get("generation") != relative.name
        or not isinstance(record.get("outputs"), dict)
        or not record["outputs"]
        or any(
            not isinstance(record.get(name), dict)
            for name in ("settings", "inputs", "implementation")
        )
    ):
        raise CacheError(f"Invalid legacy completion metadata: {stage}")
    for name in record["outputs"]:
        _relative(name)
    if MANIFEST in record["outputs"]:
        raise CacheError(f"Legacy scientific output conflicts with {MANIFEST}: {stage}")
    if (
        outputs
        and _inventory(source, root, exclude={LEGACY_MANIFEST}) != record["outputs"]
    ):
        raise CacheError(f"Legacy output inventory or hashes differ: {stage}")
    return source, record


def _select(
    root: Path, state: dict[str, Any]
) -> tuple[dict[str, tuple[Path, dict[str, Any]]], list[str]]:
    results = state["results"]
    selected = {
        stage: _legacy_result(root, stage, entry)
        for stage, entry in results.items()
        if stage in PUBLICATIONS
    }
    warnings = []
    plans = {}
    for stage, entry in results.items():
        if ATTEMPT.fullmatch(stage):
            try:
                plans[stage] = _legacy_result(root, stage, entry, outputs=False)
            except (OSError, CacheError) as error:
                if ATTEMPT.fullmatch(stage)[1] == state["pending"].get("gmm"):
                    raise
                warnings.append(f"Skipped unusable historical plan {stage}: {error}")
    attempt = state["pending"].get("gmm")
    if attempt is None and plans:
        collected = selected.get("gmm")
        matching = {}
        if collected is not None:
            identity = file_identity(collected[0] / "batch_plan.json")
            matching = {
                stage: value
                for stage, value in plans.items()
                if file_identity(value[0] / "batch_plan.json") == identity
            }
        choices = matching or plans
        chosen = max(
            choices, key=lambda stage: (choices[stage][1].get("created_at", ""), stage)
        )
        attempt = ATTEMPT.fullmatch(chosen)[1]
    if attempt is not None:
        plan_stage = f"gmm-plan-{attempt}"
        batch_stages = [
            stage
            for stage in results
            if BATCH.fullmatch(stage) and BATCH.fullmatch(stage)[1] == attempt
        ]
        if plan_stage in results:
            selected[plan_stage] = _legacy_result(root, plan_stage, results[plan_stage])
            expected_plan = file_identity(selected[plan_stage][0] / "batch_plan.json")
            for stage in batch_stages:
                selected[stage] = _legacy_result(root, stage, results[stage])
                summary = _json(selected[stage][0] / "summary.json")
                if summary.get("plan") != expected_plan or summary.get(
                    "batch_index"
                ) != int(BATCH.fullmatch(stage)[2]):
                    raise CacheError(
                        f"Saved batch differs from selected GMM plan: {stage}"
                    )
        elif batch_stages:
            raise CacheError(
                "Pending GMM batches have no saved plan; restore the plan first"
            )
    return selected, warnings


def _covered(relative: Path, occupied: set[Path]) -> bool:
    return any(relative == target or target in relative.parents for target in occupied)


def _legacy_links(root: Path, excluded: set[Path] | None = None) -> list[Path]:
    """Find visible links whose targets depend on the removable legacy tree."""
    found = []
    excluded = (excluded or set()) | {Path(".pipeline")}
    legacy = root / ".pipeline"
    for parent, directories, files in os.walk(root, followlinks=False):
        base = Path(parent)
        directories[:] = [
            name
            for name in directories
            if not _covered((base / name).relative_to(root), excluded)
        ]
        for name in directories + files:
            path = base / name
            if path.is_symlink():
                resolved = path.resolve(strict=False)
                if resolved.is_relative_to(legacy):
                    found.append(path.relative_to(root))
    return sorted(found)


def plan_migration(root: Path) -> tuple[dict[str, Any], dict[Path, dict[str, Any]]]:
    """Validate legacy sources and build a transaction without changing outputs."""
    root = Path(root).resolve()
    legacy = root / ".pipeline"
    if legacy.is_symlink() or not legacy.is_dir():
        raise CacheError(f"No ordinary legacy .pipeline directory found in {root}")
    state = _json(legacy / "state.json")
    if (
        state.get("schema_version") != 1
        or not isinstance(state.get("results"), dict)
        or not isinstance(state.get("pending"), dict)
    ):
        raise CacheError("Expected a schema-1 legacy pipeline state")
    for stage, attempt in state["pending"].items():
        if (
            not isinstance(stage, str)
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]*", stage)
            or not isinstance(attempt, str)
            or not re.fullmatch(r"[a-f0-9]{32}", attempt)
        ):
            raise CacheError("Invalid pending attempt in legacy pipeline state")
    selected, warnings = _select(root, state)
    operations: dict[Path, dict[str, Any]] = {}
    for stage, (source, record) in selected.items():
        target = _relative(_stage_path(stage).as_posix())
        _check_target(root, target)
        if target in operations:
            raise CacheError(f"Selected stages collide at {target}")
        destination = root / target
        if _present(destination) and not destination.is_symlink():
            existing = _inventory(destination, root, exclude={LEGACY_MANIFEST})
            if not destination.is_dir() or existing != record["outputs"]:
                raise CacheError(
                    f"Ordinary publication conflicts with registered active result: {destination}. "
                    "Preserve or move that publication before applying migration."
                )
        elif destination.is_symlink() and not _inside(destination, root):
            raise CacheError(
                f"Publication link escapes the output directory: {destination}"
            )
        if destination.is_symlink() and destination.resolve() != source.resolve():
            warnings.append(f"Stale published snapshot {target} retained in .pipeline")
        operations[target] = {"source": source, "stage": stage, "record": record}

    for stage, candidates in PUBLICATIONS.items():
        target = _relative(_stage_path(stage).as_posix())
        for candidate in map(Path, candidates):
            source = root / candidate
            if not _present(source):
                continue
            _check_target(root, target)
            if target not in operations:
                if not source.is_dir():
                    raise CacheError(f"Published result is not a directory: {source}")
                inventory = _inventory(source, root, allow_links=True)
                operations[target] = {"source": source, "inventory": inventory}
            elif candidate != target:
                existing = operations[target]
                if not source.is_symlink():
                    expected = existing.get(
                        "inventory", existing.get("record", {}).get("outputs")
                    )
                    if _inventory(source, root, exclude={LEGACY_MANIFEST}) != expected:
                        raise CacheError(
                            f"Ordinary publication alias conflicts with active result: {source}"
                        )
                warnings.append(
                    f"Obsolete publication alias {candidate} will be removed"
                )
            if candidate != target:
                operations[candidate] = {"remove": True}

    occupied = set(operations)
    for relative in _legacy_links(root, occupied):
        _check_target(root, relative)
        source = root / relative
        operations[relative] = {
            "source": source,
            "inventory": _inventory(source, root, allow_links=True),
        }
    selected_paths = {str(path.relative_to(root)) for path, _ in selected.values()}
    history = [
        path.relative_to(root).as_posix()
        for path in sorted((legacy / "results").glob("*/*"))
        if path.is_dir() and path.relative_to(root).as_posix() not in selected_paths
    ]
    report = {
        "status": "dry-run",
        "selected": {
            stage: _stage_path(stage).as_posix() for stage in sorted(selected)
        },
        "unregistered_publications": [
            path.as_posix()
            for path, operation in operations.items()
            if "inventory" in operation
        ],
        "skipped_keys": sorted(set(state["results"]) - set(selected)),
        "retained_legacy_history": history,
        "pending": state["pending"],
        "warnings": warnings,
        "legacy_state": file_identity(legacy / "state.json"),
        "bytes_to_copy": sum(
            identity["bytes"]
            for operation in operations.values()
            for identity in operation.get(
                "inventory", operation.get("record", {}).get("outputs", {})
            ).values()
        ),
    }
    return report, operations


def _remove(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _apply(root: Path, report: dict[str, Any], operations: dict[Path, dict]) -> None:
    if _present(root / STATE) or _present(root / RECEIPT):
        raise CacheError(
            "Direct state or migration receipt already exists; verify it before continuing"
        )
    staging = Path(tempfile.mkdtemp(prefix=".pipeline-migrate-", dir=root))
    preserve_recovery = False
    try:
        (staging / "new").mkdir()
        (staging / "backups").mkdir()
        state = {"schema_version": 2, "results": {}, "pending": report["pending"]}
        unregistered = {}
        prepared = []
        for index, (relative, operation) in enumerate(operations.items()):
            destination = staging / "new" / str(index)
            if not operation.get("remove"):
                source = operation["source"]
                if source.is_dir():
                    shutil.copytree(source, destination)
                else:
                    shutil.copyfile(source, destination)
                if "record" in operation:
                    (destination / LEGACY_MANIFEST).unlink()
                    record = operation["record"] | {
                        "schema_version": 2,
                        "result_path": relative.as_posix(),
                    }
                    if _inventory(destination, root) != record["outputs"]:
                        raise CacheError(f"Copied output identity differs: {relative}")
                    _write(destination / MANIFEST, record)
                    state["results"][operation["stage"]] = {
                        "path": relative.as_posix(),
                        "generation": record["generation"],
                        "manifest": file_identity(destination / MANIFEST),
                    }
                else:
                    inventory = _inventory(destination, root)
                    if inventory != operation["inventory"]:
                        raise CacheError(
                            f"Copied unregistered output differs: {relative}"
                        )
                    unregistered[relative.as_posix()] = {
                        "kind": "directory" if destination.is_dir() else "file",
                        "outputs": inventory,
                    }
            prepared.append(
                (relative, None if operation.get("remove") else destination)
            )
        _write(staging / "state.json", state)
        receipt = {
            "schema_version": 1,
            "status": "applied",
            "state_identity": file_identity(staging / "state.json"),
            "unregistered": unregistered,
            "report": report | {"status": "applied"},
        }
        _write(staging / "receipt.json", receipt)
        if file_identity(root / ".pipeline/state.json") != report["legacy_state"]:
            raise CacheError(
                "Legacy state changed while staging; retry with stable inputs"
            )
        prepared.extend(
            [
                (Path(STATE), staging / "state.json"),
                (Path(RECEIPT), staging / "receipt.json"),
            ]
        )
        changed = []
        try:
            for index, (relative, source) in enumerate(prepared):
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                backup = staging / "backups" / str(index)
                change = {"destination": destination, "backup": None, "placed": False}
                changed.append(change)
                if _present(destination):
                    os.replace(destination, backup)
                    change["backup"] = backup
                if source is not None:
                    os.replace(source, destination)
                    change["placed"] = True
            verify_migration(root, excluded={staging.relative_to(root)})
        except (OSError, ValueError, RuntimeError, KeyboardInterrupt) as error:
            recovery_errors = []
            for change in reversed(changed):
                try:
                    if change["placed"]:
                        _remove(change["destination"])
                    if change["backup"] is not None:
                        os.replace(change["backup"], change["destination"])
                except (OSError, ValueError, RuntimeError) as recovery_error:
                    recovery_errors.append(str(recovery_error))
            if recovery_errors:
                preserve_recovery = True
                raise CacheError(
                    f"Migration failed ({error}) and rollback could not restore every path. "
                    f"Recovery files are preserved at {staging}; restore its backups before retrying. "
                    + "; ".join(recovery_errors)
                ) from error
            raise
    finally:
        if not preserve_recovery:
            shutil.rmtree(staging)


def verify_migration(
    root: Path, *, excluded: set[Path] | None = None
) -> dict[str, Any]:
    """Recheck every migrated stage and unregistered publication before cleanup."""
    root = Path(root).resolve()
    receipt = _json(root / RECEIPT)
    if receipt.get("schema_version") != 1 or receipt.get("status") != "applied":
        raise CacheError("Missing or invalid migration receipt")
    legacy = root / ".pipeline"
    if _present(legacy):
        if legacy.is_symlink() or not legacy.is_dir():
            raise CacheError("Legacy pipeline storage is not an ordinary directory")
        expected = receipt.get("report", {}).get("legacy_state")
        try:
            unchanged = file_identity(legacy / "state.json") == expected
        except OSError as error:
            raise CacheError(
                "Cannot verify unchanged legacy pipeline state; preserve .pipeline and review its results"
            ) from error
        if not unchanged:
            raise CacheError(
                "Legacy pipeline state changed since migration; preserve .pipeline and review its newer results"
            )
    if file_identity(root / STATE) != receipt.get("state_identity"):
        raise CacheError(
            "Migrated state changed; cleanup must precede further pipeline runs"
        )
    state = _json(root / STATE)
    if state.get("schema_version") != 2 or not isinstance(state.get("results"), dict):
        raise CacheError("Invalid direct pipeline state")
    for stage, entry in state["results"].items():
        relative = _relative(entry.get("path"))
        if relative != _stage_path(stage):
            raise CacheError(f"Unexpected direct result location: {stage}")
        source = root / relative
        manifest_path = source / MANIFEST
        if file_identity(manifest_path) != entry.get("manifest"):
            raise CacheError(f"Direct completion record hash differs: {stage}")
        record = _json(manifest_path)
        if (
            record.get("schema_version") != 2
            or record.get("status") != "completed"
            or record.get("stage") != stage
            or record.get("result_path") != relative.as_posix()
            or record.get("generation") != entry.get("generation")
            or _inventory(source, root, exclude={MANIFEST}) != record.get("outputs")
        ):
            raise CacheError(f"Direct output verification failed: {stage}")
    unregistered = receipt.get("unregistered")
    if not isinstance(unregistered, dict):
        raise CacheError(
            "Migration receipt is missing unregistered publication records"
        )
    for relative, expected in unregistered.items():
        source = root / _relative(relative)
        if (expected.get("kind") == "directory") != source.is_dir() or _inventory(
            source, root
        ) != expected.get("outputs"):
            raise CacheError(
                f"Unregistered publication verification failed: {relative}"
            )
    links = _legacy_links(root, excluded)
    if links:
        raise CacheError(
            "Visible links still depend on .pipeline: " + ", ".join(map(str, links))
        )
    return {
        "status": "verified",
        "registered": sorted(state["results"]),
        "unregistered": sorted(unregistered),
    }


@contextmanager
def migration_lock(root: Path) -> Iterator[None]:
    """Interoperate with both legacy and direct-layout whole-run writer locks."""
    handles = []
    try:
        for path in (root / ".pipeline-run.lock", root / ".pipeline/run.lock"):
            if not path.parent.is_dir():
                continue
            handle = path.open("a+")
            handles.append(handle)
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise CacheError(
                    f"Another workflow or migration is writing to {root}"
                ) from error
        yield
    finally:
        for handle in reversed(handles):
            fcntl.flock(handle, fcntl.LOCK_UN)
            handle.close()


def migrate_outputs(
    root: Path,
    *,
    apply: bool = False,
    verify: bool = False,
    remove_legacy: bool = False,
) -> dict[str, Any]:
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise CacheError(f"Output directory does not exist: {root}")
    with migration_lock(root):
        recovery = sorted(root.glob(".pipeline-migrate-*"))
        if recovery:
            raise CacheError(
                "Unfinished migration workspace requires recovery before continuing: "
                + ", ".join(map(str, recovery))
            )
        if _present(root / STATE):
            report = verify_migration(root)
        else:
            if (verify or remove_legacy) and not apply:
                raise CacheError(
                    "Apply the migration before verification or legacy removal"
                )
            report, operations = plan_migration(root)
            if apply:
                _apply(root, report, operations)
                report["status"] = "applied"
        if verify or remove_legacy:
            verified = verify_migration(root)
            report["verification"] = verified
        if remove_legacy:
            legacy = root / ".pipeline"
            if legacy.is_symlink():
                raise CacheError("Refusing to remove a symbolic-link legacy directory")
            if legacy.exists():
                shutil.rmtree(legacy)
            report["legacy_removed"] = True
        return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Copy and promote direct outputs; retain legacy data",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify every migrated output without regenerating data",
    )
    parser.add_argument(
        "--remove-legacy",
        action="store_true",
        help="Verify and remove .pipeline after a successful migration",
    )
    args = parser.parse_args(argv)
    try:
        report = migrate_outputs(
            args.outdir,
            apply=args.apply,
            verify=args.verify,
            remove_legacy=args.remove_legacy,
        )
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"MaveDB output migration: {error}\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
