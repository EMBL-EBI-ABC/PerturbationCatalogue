"""Durable SQLite persistence for one self-contained curation run.

The store deliberately owns all run state.  JSON is used only as a canonical
SQLite column encoding for flexible Pydantic payloads; no workflow method reads
or writes an intermediate JSON file.
"""

from __future__ import annotations

import contextlib
import json
import os
import sqlite3
import tempfile
import unicodedata
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

DATABASE_FILENAME = "curation_run.sqlite3"
SCHEMA_VERSION = 1
RUN_NAME_MAX_LENGTH = 64
RUN_STATES = {"active", "failed", "completed"}
STEP_STATES = {"running", "completed", "failed"}
ARTIFACT_KINDS = {"evidence", "normalized", "backfilled", "final"}
REVIEW_STATUSES = {"Pending", "Approved", "Rejected", "Accepted"}


class CurationRunError(RuntimeError):
    """Raised when an operation violates the curation-run contract."""


class CurationRunConflictError(CurationRunError):
    """Raised when review changes race with a newer persisted decision."""

    def __init__(self, conflicts: list[dict[str, Any]]):
        self.conflicts = conflicts
        super().__init__(
            "Review changed in another session: " + _canonical_json(conflicts)
        )


@dataclass(frozen=True)
class RunSummary:
    run_name: str
    state: str
    revision: int


@dataclass(frozen=True)
class StepExecution:
    execution_id: str
    step: str
    state: str
    selected_item_ids: tuple[str, ...]
    retry_of: str | None


@dataclass(frozen=True)
class Artifact:
    artifact_id: str
    item_id: str
    kind: str
    payload: dict[str, Any]
    metadata: dict[str, Any]
    execution_id: str
    parent_artifact_id: str | None
    revision_id: str | None


@dataclass(frozen=True)
class ReviewChange:
    candidate_id: str
    scope_type: str
    scope_key: str = ""
    status: str = "Pending"
    term: str | None = None
    is_null: bool = False
    rationale: str | None = None
    expected_event_id: int | None = None


@dataclass(frozen=True)
class SchemaOperation:
    operation_id: str
    status: str
    old_hash: str
    new_hash: str


def utc_now() -> str:
    """Return a sortable UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def validate_run_name(value: str) -> str:
    """Validate a portable run directory name."""
    if not isinstance(value, str) or not value or len(value) > RUN_NAME_MAX_LENGTH:
        raise ValueError(
            "Run name must contain 1-64 letters, digits, underscores, or hyphens"
        )
    if value in {".", ".."} or not value[0].isalnum():
        raise ValueError("Run name must begin with a letter or digit")
    if any(
        character
        not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
        for character in value
    ):
        raise ValueError(
            "Run name must contain only letters, digits, underscores, or hyphens"
        )
    return value


def candidate_id_for(field_name: str, original_term: str) -> str:
    """Return a stable candidate identity across discovery retries."""
    identity = "\0".join(
        (
            unicodedata.normalize("NFKC", str(field_name).casefold()),
            unicodedata.normalize("NFKC", str(original_term).casefold()),
        )
    )
    return sha256(identity.encode("utf-8")).hexdigest()


def _canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _json_hash(value: object) -> str:
    return sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _text_hash(value: str) -> str:
    return sha256(value.encode("utf-8")).hexdigest()


def _source_label(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Source labels must be non-empty strings")
    return Path(value).name


class CurationRunStore:
    """Persistence module for a fresh, self-contained curation run."""

    def __init__(self, run_dir: str | Path):
        self.run_dir = Path(run_dir).resolve()
        self.database_path = self.run_dir / DATABASE_FILENAME

    @classmethod
    def create(
        cls,
        runs_root: str | Path,
        run_name: str,
        configuration: Mapping[str, Any],
    ) -> "CurationRunStore":
        """Create a new empty run; existing directories are never reused."""
        name = validate_run_name(run_name)
        run_dir = Path(runs_root).resolve() / name
        if run_dir.exists():
            raise FileExistsError(f"Curation run already exists: {run_dir}")
        run_dir.mkdir(parents=True, exist_ok=False)
        store = cls(run_dir)
        try:
            with store._connection() as connection:
                store._initialize_schema(connection)
                now = utc_now()
                values = {
                    "run_name": name,
                    "state": "active",
                    "revision": "0",
                    "configuration_json": _canonical_json(dict(configuration)),
                    "created_at": now,
                    "updated_at": now,
                }
                connection.executemany(
                    "INSERT INTO run_metadata(key, value) VALUES (?, ?)", values.items()
                )
                connection.commit()
        except Exception:
            with contextlib.suppress(OSError):
                store.database_path.unlink()
            with contextlib.suppress(OSError):
                run_dir.rmdir()
            raise
        return store

    @classmethod
    def open(cls, runs_root: str | Path, run_name: str) -> "CurationRunStore":
        """Open a native run and refuse legacy run directories."""
        name = validate_run_name(run_name)
        run_dir = Path(runs_root).resolve() / name
        legacy_markers = (
            run_dir / "review_state.sqlite3",
            run_dir / "pipeline_manifest.json",
            run_dir / "step3_ontology_candidates",
        )
        if any(marker.exists() for marker in legacy_markers):
            raise CurationRunError(
                f"Legacy curation run is unsupported: {run_dir}. Create a new run instead."
            )
        store = cls(run_dir)
        if not store.database_path.is_file():
            raise FileNotFoundError(
                f"Curation run database not found: {store.database_path}"
            )
        with store._connection() as connection:
            store._initialize_schema(connection)
            saved_name = store._metadata_value(connection, "run_name")
            if saved_name != name:
                raise CurationRunError(
                    "Run directory and database run name do not match"
                )
        return store

    @contextlib.contextmanager
    def _connection(self) -> Iterator[sqlite3.Connection]:
        connection = sqlite3.connect(self.database_path, timeout=10)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        connection.execute("PRAGMA journal_mode = WAL")
        connection.execute("PRAGMA busy_timeout = 10000")
        try:
            yield connection
        finally:
            connection.close()

    @staticmethod
    def _initialize_schema(connection: sqlite3.Connection) -> None:
        connection.execute(
            "CREATE TABLE IF NOT EXISTS run_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        version_row = connection.execute(
            "SELECT value FROM run_metadata WHERE key = 'schema_version'"
        ).fetchone()
        if version_row is None:
            connection.execute(
                "INSERT INTO run_metadata(key, value) VALUES ('schema_version', ?)",
                (str(SCHEMA_VERSION),),
            )
        elif int(version_row["value"]) != SCHEMA_VERSION:
            raise CurationRunError(
                f"Unsupported curation-run schema version: {version_row['value']}"
            )

        connection.executescript("""
            CREATE TABLE IF NOT EXISTS source_documents (
                document_id TEXT PRIMARY KEY,
                source_label TEXT NOT NULL UNIQUE,
                doi TEXT,
                original_path TEXT NOT NULL,
                content_text TEXT NOT NULL,
                content_hash TEXT NOT NULL,
                ingested_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS mavedb_snapshots (
                snapshot_id TEXT PRIMARY KEY,
                urn TEXT NOT NULL UNIQUE,
                source_label TEXT NOT NULL,
                original_path TEXT NOT NULL,
                payload_json TEXT NOT NULL,
                prompt_metadata_json TEXT NOT NULL,
                payload_hash TEXT NOT NULL,
                ingested_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS curation_items (
                item_id TEXT PRIMARY KEY,
                item_label TEXT NOT NULL UNIQUE,
                primary_dataset_id TEXT,
                prompt_context_json TEXT NOT NULL,
                context_hash TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS item_documents (
                item_id TEXT NOT NULL REFERENCES curation_items(item_id) ON DELETE CASCADE,
                document_id TEXT NOT NULL REFERENCES source_documents(document_id) ON DELETE RESTRICT,
                PRIMARY KEY(item_id, document_id)
            );
            CREATE TABLE IF NOT EXISTS item_mavedb_snapshots (
                item_id TEXT NOT NULL REFERENCES curation_items(item_id) ON DELETE CASCADE,
                snapshot_id TEXT NOT NULL REFERENCES mavedb_snapshots(snapshot_id) ON DELETE RESTRICT,
                PRIMARY KEY(item_id, snapshot_id)
            );
            CREATE TABLE IF NOT EXISTS step_executions (
                execution_id TEXT PRIMARY KEY,
                step TEXT NOT NULL,
                state TEXT NOT NULL CHECK(state IN ('running', 'completed', 'failed')),
                selected_item_ids_json TEXT NOT NULL,
                configuration_json TEXT NOT NULL,
                configuration_hash TEXT NOT NULL,
                retry_of TEXT REFERENCES step_executions(execution_id),
                started_at TEXT NOT NULL,
                completed_at TEXT,
                error_summary TEXT
            );
            CREATE INDEX IF NOT EXISTS step_executions_step_idx
            ON step_executions(step, started_at);
            CREATE TABLE IF NOT EXISTS artifact_versions (
                artifact_seq INTEGER PRIMARY KEY AUTOINCREMENT,
                artifact_id TEXT NOT NULL UNIQUE,
                item_id TEXT NOT NULL REFERENCES curation_items(item_id) ON DELETE RESTRICT,
                kind TEXT NOT NULL CHECK(kind IN ('evidence', 'normalized', 'backfilled', 'final')),
                execution_id TEXT NOT NULL REFERENCES step_executions(execution_id) ON DELETE RESTRICT,
                parent_artifact_id TEXT REFERENCES artifact_versions(artifact_id),
                revision_id TEXT,
                payload_json TEXT NOT NULL,
                payload_hash TEXT NOT NULL,
                metadata_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS artifact_versions_latest_idx
            ON artifact_versions(item_id, kind, artifact_seq DESC);
            CREATE TABLE IF NOT EXISTS discoveries (
                discovery_id TEXT PRIMARY KEY,
                execution_id TEXT NOT NULL REFERENCES step_executions(execution_id),
                input_artifact_ids_json TEXT NOT NULL,
                input_hash TEXT NOT NULL,
                status TEXT NOT NULL CHECK(status IN ('completed', 'failed')),
                error_summary TEXT,
                created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS candidates (
                candidate_id TEXT PRIMARY KEY,
                discovery_id TEXT NOT NULL REFERENCES discoveries(discovery_id) ON DELETE RESTRICT,
                field_name TEXT NOT NULL,
                original_term TEXT NOT NULL,
                payload_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS candidate_sources (
                candidate_id TEXT NOT NULL REFERENCES candidates(candidate_id) ON DELETE CASCADE,
                item_id TEXT NOT NULL REFERENCES curation_items(item_id) ON DELETE RESTRICT,
                evidence_json TEXT NOT NULL,
                evidence_hash TEXT NOT NULL,
                PRIMARY KEY(candidate_id, item_id, evidence_hash)
            );
            CREATE TABLE IF NOT EXISTS review_events (
                event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                candidate_id TEXT NOT NULL REFERENCES candidates(candidate_id) ON DELETE RESTRICT,
                scope_type TEXT NOT NULL CHECK(scope_type IN ('global', 'item')),
                scope_key TEXT NOT NULL DEFAULT '',
                status TEXT NOT NULL CHECK(status IN ('Pending', 'Approved', 'Rejected', 'Accepted')),
                term TEXT,
                is_null INTEGER NOT NULL DEFAULT 0 CHECK(is_null IN (0, 1)),
                rationale TEXT,
                session_id TEXT NOT NULL,
                created_at TEXT NOT NULL,
                previous_event_id INTEGER REFERENCES review_events(event_id)
            );
            CREATE INDEX IF NOT EXISTS review_events_latest_idx
            ON review_events(candidate_id, scope_type, scope_key, event_id DESC);
            CREATE TABLE IF NOT EXISTS backfill_revisions (
                revision_id TEXT PRIMARY KEY,
                execution_id TEXT NOT NULL REFERENCES step_executions(execution_id),
                review_revision INTEGER NOT NULL,
                source_artifact_ids_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS final_edit_events (
                event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                artifact_id TEXT NOT NULL REFERENCES artifact_versions(artifact_id) ON DELETE RESTRICT,
                field_name TEXT NOT NULL,
                old_value_json TEXT,
                new_value_json TEXT NOT NULL,
                session_id TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS schema_operations (
                operation_id TEXT PRIMARY KEY,
                status TEXT NOT NULL CHECK(status IN ('prepared', 'applied', 'failed')),
                schema_file_path TEXT NOT NULL,
                old_source_text TEXT NOT NULL,
                new_source_text TEXT NOT NULL,
                old_hash TEXT NOT NULL,
                new_hash TEXT NOT NULL,
                diff_text TEXT NOT NULL,
                approved_terms_json TEXT NOT NULL,
                prepared_at TEXT NOT NULL,
                finalized_at TEXT,
                error_summary TEXT
            );
            CREATE TABLE IF NOT EXISTS run_events (
                event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                execution_id TEXT REFERENCES step_executions(execution_id),
                item_id TEXT REFERENCES curation_items(item_id),
                level TEXT NOT NULL,
                message TEXT NOT NULL,
                details_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS run_events_time_idx ON run_events(event_id DESC);
            CREATE TABLE IF NOT EXISTS exports (
                export_id TEXT PRIMARY KEY,
                output_dir TEXT NOT NULL,
                format TEXT NOT NULL,
                content_hash TEXT NOT NULL,
                final_artifact_ids_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            """)
        connection.commit()

    @staticmethod
    def _metadata_value(connection: sqlite3.Connection, key: str) -> str | None:
        row = connection.execute(
            "SELECT value FROM run_metadata WHERE key = ?", (key,)
        ).fetchone()
        return str(row["value"]) if row else None

    def _assert_active(self, connection: sqlite3.Connection) -> None:
        state = self._metadata_value(connection, "state")
        if state != "active":
            raise CurationRunError(
                f"Run is {state}; completed and failed runs are read-only"
            )

    def _bump_revision(self, connection: sqlite3.Connection) -> int:
        current = int(self._metadata_value(connection, "revision") or "0") + 1
        now = utc_now()
        connection.execute(
            "UPDATE run_metadata SET value = ? WHERE key = 'revision'", (str(current),)
        )
        connection.execute(
            "UPDATE run_metadata SET value = ? WHERE key = 'updated_at'", (now,)
        )
        return current

    def summary(self) -> RunSummary:
        with self._connection() as connection:
            return RunSummary(
                run_name=self._metadata_value(connection, "run_name") or "",
                state=self._metadata_value(connection, "state") or "",
                revision=int(self._metadata_value(connection, "revision") or "0"),
            )

    def configuration(self) -> dict[str, Any]:
        with self._connection() as connection:
            payload = self._metadata_value(connection, "configuration_json") or "{}"
        return json.loads(payload)

    def snapshot_document(
        self,
        source_label: str,
        content_text: str,
        original_path: str | Path,
        doi: str | None = None,
    ) -> str:
        """Persist an immutable publication-text snapshot and return its ID."""
        label = _source_label(source_label)
        if not content_text:
            raise ValueError(f"Publication text is empty: {label}")
        document_id = str(uuid.uuid4())
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                existing = connection.execute(
                    "SELECT document_id, content_hash FROM source_documents WHERE source_label = ?",
                    (label,),
                ).fetchone()
                content_hash = _text_hash(content_text)
                if existing:
                    if existing["content_hash"] != content_hash:
                        raise CurationRunError(
                            f"Source label already has a different snapshot: {label}"
                        )
                    connection.rollback()
                    return str(existing["document_id"])
                connection.execute(
                    """INSERT INTO source_documents(
                    document_id, source_label, doi, original_path, content_text, content_hash, ingested_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        document_id,
                        label,
                        doi,
                        str(Path(original_path).resolve()),
                        content_text,
                        content_hash,
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return document_id

    def snapshot_mavedb_entry(
        self,
        urn: str,
        source_label: str,
        payload: Mapping[str, Any],
        prompt_metadata: Mapping[str, Any],
        original_path: str | Path,
    ) -> str:
        """Persist an immutable MaveDB metadata snapshot and return its ID."""
        if not urn:
            raise ValueError("MaveDB URN is required")
        snapshot_id = str(uuid.uuid4())
        payload_json = _canonical_json(dict(payload))
        payload_hash = _text_hash(payload_json)
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                existing = connection.execute(
                    "SELECT snapshot_id, payload_hash FROM mavedb_snapshots WHERE urn = ?",
                    (urn,),
                ).fetchone()
                if existing:
                    if existing["payload_hash"] != payload_hash:
                        raise CurationRunError(
                            f"MaveDB URN already has a different snapshot: {urn}"
                        )
                    connection.rollback()
                    return str(existing["snapshot_id"])
                connection.execute(
                    """INSERT INTO mavedb_snapshots(
                    snapshot_id, urn, source_label, original_path, payload_json, prompt_metadata_json,
                    payload_hash, ingested_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        snapshot_id,
                        urn,
                        _source_label(source_label),
                        str(Path(original_path).resolve()),
                        payload_json,
                        _canonical_json(dict(prompt_metadata)),
                        payload_hash,
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return snapshot_id

    def create_item(
        self,
        item_label: str,
        document_ids: Sequence[str],
        snapshot_ids: Sequence[str],
        prompt_context: Mapping[str, Any],
        primary_dataset_id: str | None,
    ) -> str:
        """Create one immutable prompt-context work item."""
        if not document_ids:
            raise ValueError("A curation item requires at least one source document")
        if not snapshot_ids:
            raise ValueError("A curation item requires at least one MaveDB snapshot")
        item_id = str(uuid.uuid4())
        context_json = _canonical_json(dict(prompt_context))
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """INSERT INTO curation_items(
                    item_id, item_label, primary_dataset_id, prompt_context_json, context_hash, created_at
                    ) VALUES (?, ?, ?, ?, ?, ?)""",
                    (
                        item_id,
                        _source_label(item_label),
                        primary_dataset_id,
                        context_json,
                        _text_hash(context_json),
                        utc_now(),
                    ),
                )
                connection.executemany(
                    "INSERT INTO item_documents(item_id, document_id) VALUES (?, ?)",
                    [(item_id, document_id) for document_id in document_ids],
                )
                connection.executemany(
                    "INSERT INTO item_mavedb_snapshots(item_id, snapshot_id) VALUES (?, ?)",
                    [(item_id, snapshot_id) for snapshot_id in snapshot_ids],
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return item_id

    def list_items(self) -> list[dict[str, Any]]:
        with self._connection() as connection:
            rows = connection.execute(
                "SELECT item_id, item_label, primary_dataset_id, prompt_context_json FROM curation_items ORDER BY item_label"
            ).fetchall()
            return [
                {
                    "item_id": row["item_id"],
                    "item_label": row["item_label"],
                    "primary_dataset_id": row["primary_dataset_id"],
                    "prompt_context": json.loads(row["prompt_context_json"]),
                }
                for row in rows
            ]

    def load_item(self, item_id: str) -> dict[str, Any]:
        with self._connection() as connection:
            item = connection.execute(
                "SELECT * FROM curation_items WHERE item_id = ?", (item_id,)
            ).fetchone()
            if item is None:
                raise KeyError(f"Unknown curation item: {item_id}")
            documents = connection.execute(
                """SELECT d.* FROM source_documents d
                JOIN item_documents link ON link.document_id = d.document_id
                WHERE link.item_id = ? ORDER BY d.source_label""",
                (item_id,),
            ).fetchall()
            snapshots = connection.execute(
                """SELECT m.* FROM mavedb_snapshots m
                JOIN item_mavedb_snapshots link ON link.snapshot_id = m.snapshot_id
                WHERE link.item_id = ? ORDER BY m.urn""",
                (item_id,),
            ).fetchall()
        return {
            "item_id": item["item_id"],
            "item_label": item["item_label"],
            "primary_dataset_id": item["primary_dataset_id"],
            "prompt_context": json.loads(item["prompt_context_json"]),
            "documents": [dict(row) for row in documents],
            "mavedb_snapshots": [
                {
                    **dict(row),
                    "payload": json.loads(row["payload_json"]),
                    "prompt_metadata": json.loads(row["prompt_metadata_json"]),
                }
                for row in snapshots
            ],
        }

    def start_step(
        self,
        step: str,
        selected_item_ids: Sequence[str],
        configuration: Mapping[str, Any],
        retry_of: str | None = None,
    ) -> StepExecution:
        """Create one running execution in the active run."""
        if not step:
            raise ValueError("Step name is required")
        selected = tuple(sorted(set(selected_item_ids)))
        execution_id = str(uuid.uuid4())
        configuration_json = _canonical_json(dict(configuration))
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                known_items = {
                    row["item_id"]
                    for row in connection.execute("SELECT item_id FROM curation_items")
                }
                unknown = set(selected) - known_items
                if unknown:
                    raise KeyError(f"Unknown curation item IDs: {sorted(unknown)}")
                connection.execute(
                    """INSERT INTO step_executions(
                    execution_id, step, state, selected_item_ids_json, configuration_json,
                    configuration_hash, retry_of, started_at
                    ) VALUES (?, ?, 'running', ?, ?, ?, ?, ?)""",
                    (
                        execution_id,
                        step,
                        _canonical_json(selected),
                        configuration_json,
                        _text_hash(configuration_json),
                        retry_of,
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return StepExecution(execution_id, step, "running", selected, retry_of)

    def finish_step(self, execution_id: str) -> None:
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                row = connection.execute(
                    "SELECT state FROM step_executions WHERE execution_id = ?",
                    (execution_id,),
                ).fetchone()
                if row is None:
                    raise KeyError(f"Unknown execution: {execution_id}")
                if row["state"] != "running":
                    raise CurationRunError(f"Execution is already {row['state']}")
                connection.execute(
                    "UPDATE step_executions SET state = 'completed', completed_at = ? WHERE execution_id = ?",
                    (utc_now(), execution_id),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise

    def fail_step(self, execution_id: str, error_summary: str) -> None:
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """UPDATE step_executions SET state = 'failed', completed_at = ?, error_summary = ?
                    WHERE execution_id = ? AND state = 'running'""",
                    (utc_now(), error_summary, execution_id),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise

    def latest_failed_execution(self, step: str) -> str | None:
        with self._connection() as connection:
            row = connection.execute(
                """SELECT execution_id FROM step_executions
                WHERE step = ? AND state = 'failed' ORDER BY started_at DESC LIMIT 1""",
                (step,),
            ).fetchone()
        return str(row["execution_id"]) if row else None

    def step_history(self) -> list[dict[str, Any]]:
        """Return execution status without exposing database details to callers."""
        with self._connection() as connection:
            rows = connection.execute(
                """SELECT execution_id, step, state, selected_item_ids_json, retry_of,
                started_at, completed_at, error_summary
                FROM step_executions ORDER BY started_at"""
            ).fetchall()
        return [
            {
                **dict(row),
                "selected_item_ids": tuple(json.loads(row["selected_item_ids_json"])),
            }
            for row in rows
        ]

    def put_artifact(
        self,
        item_id: str,
        kind: str,
        payload: Mapping[str, Any],
        metadata: Mapping[str, Any],
        execution_id: str,
        parent_artifact_id: str | None = None,
        revision_id: str | None = None,
    ) -> Artifact:
        """Append a validated immutable artifact produced by a running execution."""
        if kind not in ARTIFACT_KINDS:
            raise ValueError(f"Unsupported artifact kind: {kind}")
        artifact_id = str(uuid.uuid4())
        payload_dict = dict(payload)
        metadata_dict = dict(metadata)
        payload_json = _canonical_json(payload_dict)
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                execution = connection.execute(
                    "SELECT state FROM step_executions WHERE execution_id = ?",
                    (execution_id,),
                ).fetchone()
                if execution is None or execution["state"] != "running":
                    raise CurationRunError(
                        "Artifacts can only be written by a running execution"
                    )
                connection.execute(
                    """INSERT INTO artifact_versions(
                    artifact_id, item_id, kind, execution_id, parent_artifact_id, revision_id,
                    payload_json, payload_hash, metadata_json, created_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        artifact_id,
                        item_id,
                        kind,
                        execution_id,
                        parent_artifact_id,
                        revision_id,
                        payload_json,
                        _text_hash(payload_json),
                        _canonical_json(metadata_dict),
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return Artifact(
            artifact_id,
            item_id,
            kind,
            payload_dict,
            metadata_dict,
            execution_id,
            parent_artifact_id,
            revision_id,
        )

    def latest_artifact(self, item_id: str, kind: str) -> Artifact | None:
        with self._connection() as connection:
            row = connection.execute(
                """SELECT * FROM artifact_versions WHERE item_id = ? AND kind = ?
                ORDER BY artifact_seq DESC LIMIT 1""",
                (item_id, kind),
            ).fetchone()
        return self._artifact_from_row(row) if row else None

    def latest_artifacts(
        self, kind: str, item_ids: Sequence[str] | None = None
    ) -> list[Artifact]:
        with self._connection() as connection:
            parameters: list[object] = [kind]
            filter_sql = ""
            if item_ids is not None:
                selected = tuple(sorted(set(item_ids)))
                if not selected:
                    return []
                placeholders = ",".join("?" for _ in selected)
                filter_sql = f" AND item_id IN ({placeholders})"
                parameters.extend(selected)
            rows = connection.execute(
                f"""SELECT artifact.* FROM artifact_versions artifact
                JOIN (
                    SELECT item_id, MAX(artifact_seq) AS max_seq
                    FROM artifact_versions WHERE kind = ?{filter_sql} GROUP BY item_id
                ) latest ON latest.max_seq = artifact.artifact_seq
                ORDER BY artifact.item_id""",
                parameters,
            ).fetchall()
        return [self._artifact_from_row(row) for row in rows]

    @staticmethod
    def _artifact_from_row(row: sqlite3.Row) -> Artifact:
        return Artifact(
            artifact_id=str(row["artifact_id"]),
            item_id=str(row["item_id"]),
            kind=str(row["kind"]),
            payload=json.loads(row["payload_json"]),
            metadata=json.loads(row["metadata_json"]),
            execution_id=str(row["execution_id"]),
            parent_artifact_id=row["parent_artifact_id"],
            revision_id=row["revision_id"],
        )

    def record_discovery(
        self,
        execution_id: str,
        input_artifact_ids: Sequence[str],
        candidates: Mapping[str, Sequence[Mapping[str, Any]]],
        status: str = "completed",
        error_summary: str | None = None,
    ) -> str:
        """Persist a Step 3 discovery and its immutable candidate evidence."""
        if status not in {"completed", "failed"}:
            raise ValueError("Discovery status must be completed or failed")
        discovery_id = str(uuid.uuid4())
        input_ids = tuple(sorted(set(input_artifact_ids)))
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """INSERT INTO discoveries(
                    discovery_id, execution_id, input_artifact_ids_json, input_hash, status, error_summary, created_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        discovery_id,
                        execution_id,
                        _canonical_json(input_ids),
                        _json_hash(input_ids),
                        status,
                        error_summary,
                        utc_now(),
                    ),
                )
                for field_name, entries in candidates.items():
                    for entry in entries:
                        original_term = str(
                            entry.get("original_term")
                            or entry.get("proposed_new_term")
                            or "Other"
                        ).strip()
                        candidate_id = candidate_id_for(field_name, original_term)
                        candidate_payload = dict(entry)
                        candidate_payload.setdefault("original_term", original_term)
                        candidate_payload.setdefault("term", original_term)
                        candidate_payload.setdefault("status", "Pending")
                        connection.execute(
                            """INSERT INTO candidates(
                            candidate_id, discovery_id, field_name, original_term, payload_json, created_at
                            ) VALUES (?, ?, ?, ?, ?, ?)
                            ON CONFLICT(candidate_id) DO UPDATE SET
                                discovery_id = excluded.discovery_id,
                                field_name = excluded.field_name,
                                original_term = excluded.original_term,
                                payload_json = excluded.payload_json""",
                            (
                                candidate_id,
                                discovery_id,
                                field_name,
                                original_term,
                                _canonical_json(candidate_payload),
                                utc_now(),
                            ),
                        )
                        for evidence in candidate_payload.get(
                            "supporting_evidence", []
                        ):
                            if not isinstance(evidence, Mapping):
                                continue
                            item_id = str(evidence.get("item_id") or "")
                            if not item_id:
                                continue
                            evidence_json = _canonical_json(dict(evidence))
                            connection.execute(
                                """INSERT OR IGNORE INTO candidate_sources(
                                candidate_id, item_id, evidence_json, evidence_hash
                                ) VALUES (?, ?, ?, ?)""",
                                (
                                    candidate_id,
                                    item_id,
                                    evidence_json,
                                    _text_hash(evidence_json),
                                ),
                            )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return discovery_id

    def load_review_rows(self) -> list[dict[str, Any]]:
        """Return effective global and item-scoped decisions for the latest discovery."""
        with self._connection() as connection:
            rows = connection.execute("""SELECT c.* FROM candidates c
                JOIN discoveries d ON d.discovery_id = c.discovery_id
                WHERE d.status = 'completed'
                ORDER BY d.created_at DESC, c.field_name, c.original_term""").fetchall()
            result: list[dict[str, Any]] = []
            seen: set[str] = set()
            for candidate in rows:
                candidate_id = str(candidate["candidate_id"])
                if candidate_id in seen:
                    continue
                seen.add(candidate_id)
                payload = json.loads(candidate["payload_json"])
                global_event = self._latest_review_event(
                    connection, candidate_id, "global", ""
                )
                sources = connection.execute(
                    "SELECT item_id, evidence_json FROM candidate_sources WHERE candidate_id = ? ORDER BY item_id",
                    (candidate_id,),
                ).fetchall()
                item_decisions: dict[str, dict[str, Any]] = {}
                for source in sources:
                    item_id = str(source["item_id"])
                    event = self._latest_review_event(
                        connection, candidate_id, "item", item_id
                    )
                    if event:
                        item_decisions[item_id] = self._event_to_decision(event)
                decision = (
                    self._event_to_decision(global_event)
                    if global_event
                    else {
                        "event_id": None,
                        "status": "Pending",
                        "term": payload.get("term") or candidate["original_term"],
                        "is_null": False,
                        "rationale": payload.get("rationale"),
                    }
                )
                result.append(
                    {
                        "candidate_id": candidate_id,
                        "field_name": candidate["field_name"],
                        "original_term": candidate["original_term"],
                        "term": decision["term"],
                        "status": decision["status"],
                        "is_null": decision["is_null"],
                        "rationale": decision["rationale"],
                        "event_id": decision["event_id"],
                        "supporting_evidence": [
                            json.loads(source["evidence_json"]) for source in sources
                        ],
                        "item_decisions": item_decisions,
                    }
                )
        return result

    @staticmethod
    def _latest_review_event(
        connection: sqlite3.Connection,
        candidate_id: str,
        scope_type: str,
        scope_key: str,
    ) -> sqlite3.Row | None:
        return connection.execute(
            """SELECT * FROM review_events WHERE candidate_id = ? AND scope_type = ? AND scope_key = ?
            ORDER BY event_id DESC LIMIT 1""",
            (candidate_id, scope_type, scope_key),
        ).fetchone()

    @staticmethod
    def _event_to_decision(event: sqlite3.Row) -> dict[str, Any]:
        return {
            "event_id": int(event["event_id"]),
            "status": str(event["status"]),
            "term": event["term"],
            "is_null": bool(event["is_null"]),
            "rationale": event["rationale"],
        }

    def commit_review_changes(
        self, changes: Iterable[ReviewChange], session_id: str
    ) -> int:
        """Atomically append conflict-checked review events and return the new revision."""
        queued = list(changes)
        if not queued:
            return self.summary().revision
        conflicts: list[dict[str, Any]] = []
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                for change in queued:
                    if change.scope_type not in {"global", "item"}:
                        raise ValueError("Review scope must be global or item")
                    if change.status not in REVIEW_STATUSES:
                        raise ValueError(f"Invalid review status: {change.status}")
                    scope_key = (
                        "" if change.scope_type == "global" else change.scope_key
                    )
                    latest = self._latest_review_event(
                        connection, change.candidate_id, change.scope_type, scope_key
                    )
                    latest_id = int(latest["event_id"]) if latest else None
                    if change.expected_event_id != latest_id:
                        conflicts.append(
                            {
                                "candidate_id": change.candidate_id,
                                "scope_type": change.scope_type,
                                "scope_key": scope_key,
                                "expected_event_id": change.expected_event_id,
                                "actual_event_id": latest_id,
                            }
                        )
                if conflicts:
                    raise CurationRunConflictError(conflicts)
                for change in queued:
                    scope_key = (
                        "" if change.scope_type == "global" else change.scope_key
                    )
                    latest = self._latest_review_event(
                        connection, change.candidate_id, change.scope_type, scope_key
                    )
                    connection.execute(
                        """INSERT INTO review_events(
                        candidate_id, scope_type, scope_key, status, term, is_null, rationale,
                        session_id, created_at, previous_event_id
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                        (
                            change.candidate_id,
                            change.scope_type,
                            scope_key,
                            change.status,
                            change.term,
                            int(change.is_null),
                            change.rationale,
                            session_id,
                            utc_now(),
                            latest["event_id"] if latest else None,
                        ),
                    )
                revision = self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return revision

    def resolve_approved_terms(self) -> dict[str, dict[str, str | None]]:
        """Resolve global terms then item overrides for Step 4 preview and commit."""
        updates: dict[str, dict[str, str | None]] = {}
        for row in self.load_review_rows():
            status = row["status"]
            if status in {"Approved", "Accepted"}:
                value = None if row["is_null"] else str(row["term"] or "").strip()
                if value is None or value:
                    for evidence in row["supporting_evidence"]:
                        item_id = evidence.get("item_id")
                        if item_id:
                            updates.setdefault(str(item_id), {})[
                                row["field_name"]
                            ] = value
            for item_id, decision in row["item_decisions"].items():
                if decision["status"] not in {"Approved", "Accepted"}:
                    continue
                value = (
                    None if decision["is_null"] else str(decision["term"] or "").strip()
                )
                if value is None or value:
                    updates.setdefault(item_id, {})[row["field_name"]] = value
        return updates

    def approved_schema_terms(self) -> dict[str, list[str]]:
        result: dict[str, list[str]] = {}
        for field_name, updates in self.resolve_approved_terms().items():
            del field_name
            for field, value in updates.items():
                if value and value not in result.setdefault(field, []):
                    result[field].append(value)
        return result

    def create_backfill_revision(
        self, execution_id: str, source_artifact_ids: Sequence[str]
    ) -> str:
        revision_id = str(uuid.uuid4())
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """INSERT INTO backfill_revisions(
                    revision_id, execution_id, review_revision, source_artifact_ids_json, created_at
                    ) VALUES (?, ?, ?, ?, ?)""",
                    (
                        revision_id,
                        execution_id,
                        int(self._metadata_value(connection, "revision") or "0"),
                        _canonical_json(tuple(sorted(source_artifact_ids))),
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return revision_id

    def save_final_edit(
        self, artifact_id: str, field_name: str, new_value: object, session_id: str
    ) -> int:
        """Append a final-payload edit instead of mutating its source artifact."""
        artifact = self.artifact_by_id(artifact_id)
        old_value = self.effective_final_payload(artifact_id).get(field_name)
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                cursor = connection.execute(
                    """INSERT INTO final_edit_events(
                    artifact_id, field_name, old_value_json, new_value_json, session_id, created_at
                    ) VALUES (?, ?, ?, ?, ?, ?)""",
                    (
                        artifact.artifact_id,
                        field_name,
                        _canonical_json(old_value),
                        _canonical_json(new_value),
                        session_id,
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
                return int(cursor.lastrowid)
            except Exception:
                connection.rollback()
                raise

    def artifact_by_id(self, artifact_id: str) -> Artifact:
        with self._connection() as connection:
            row = connection.execute(
                "SELECT * FROM artifact_versions WHERE artifact_id = ?", (artifact_id,)
            ).fetchone()
        if row is None:
            raise KeyError(f"Unknown artifact: {artifact_id}")
        return self._artifact_from_row(row)

    def effective_final_payload(self, artifact_id: str) -> dict[str, Any]:
        artifact = self.artifact_by_id(artifact_id)
        if artifact.kind != "final":
            raise ValueError("Final edits apply only to final artifacts")
        payload = dict(artifact.payload)
        with self._connection() as connection:
            rows = connection.execute(
                """SELECT field_name, new_value_json FROM final_edit_events
                WHERE artifact_id = ? ORDER BY event_id""",
                (artifact_id,),
            ).fetchall()
        for row in rows:
            payload[row["field_name"]] = json.loads(row["new_value_json"])
        return payload

    def begin_schema_operation(
        self,
        schema_file_path: str | Path,
        old_source_text: str,
        new_source_text: str,
        diff_text: str,
        approved_terms: Mapping[str, Sequence[str]],
    ) -> SchemaOperation:
        operation_id = str(uuid.uuid4())
        old_hash = _text_hash(old_source_text)
        new_hash = _text_hash(new_source_text)
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """INSERT INTO schema_operations(
                    operation_id, status, schema_file_path, old_source_text, new_source_text,
                    old_hash, new_hash, diff_text, approved_terms_json, prepared_at
                    ) VALUES (?, 'prepared', ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        operation_id,
                        str(Path(schema_file_path).resolve()),
                        old_source_text,
                        new_source_text,
                        old_hash,
                        new_hash,
                        diff_text,
                        _canonical_json(
                            {key: list(value) for key, value in approved_terms.items()}
                        ),
                        utc_now(),
                    ),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        return SchemaOperation(operation_id, "prepared", old_hash, new_hash)

    def finish_schema_operation(
        self, operation_id: str, status: str, error_summary: str | None = None
    ) -> None:
        if status not in {"applied", "failed"}:
            raise ValueError("Schema operation status must be applied or failed")
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """UPDATE schema_operations SET status = ?, finalized_at = ?, error_summary = ?
                    WHERE operation_id = ? AND status = 'prepared'""",
                    (status, utc_now(), error_summary, operation_id),
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise

    def log_event(
        self,
        level: str,
        message: str,
        execution_id: str | None = None,
        item_id: str | None = None,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute(
                """INSERT INTO run_events(execution_id, item_id, level, message, details_json, created_at)
                VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    execution_id,
                    item_id,
                    level,
                    message,
                    _canonical_json(dict(details or {})),
                    utc_now(),
                ),
            )
            connection.commit()

    def recent_events(self, limit: int = 100) -> list[dict[str, Any]]:
        with self._connection() as connection:
            rows = connection.execute(
                "SELECT * FROM run_events ORDER BY event_id DESC LIMIT ?", (limit,)
            ).fetchall()
        return [
            {
                **dict(row),
                "details": json.loads(row["details_json"]),
            }
            for row in reversed(rows)
        ]

    def complete_run(self) -> None:
        with self._connection() as connection:
            self._assert_active(connection)
            connection.execute("BEGIN IMMEDIATE")
            try:
                unfinished = connection.execute(
                    "SELECT step FROM step_executions WHERE state = 'running' LIMIT 1"
                ).fetchone()
                if unfinished:
                    raise CurationRunError("Cannot complete a run with a running step")
                final_artifact = connection.execute(
                    "SELECT 1 FROM artifact_versions WHERE kind = 'final' LIMIT 1"
                ).fetchone()
                if final_artifact is None:
                    raise CurationRunError(
                        "Cannot complete a run without final artifacts"
                    )
                connection.execute(
                    "UPDATE run_metadata SET value = 'completed' WHERE key = 'state'"
                )
                self._bump_revision(connection)
                connection.commit()
            except Exception:
                connection.rollback()
                raise

    def record_export(
        self,
        output_dir: str | Path,
        format_name: str,
        content_hash: str,
        final_artifact_ids: Sequence[str],
    ) -> str:
        """Record an explicit final-deliverable export; completed runs allow this write."""
        export_id = str(uuid.uuid4())
        with self._connection() as connection:
            if self._metadata_value(connection, "state") != "completed":
                raise CurationRunError(
                    "Final export is available only after run completion"
                )
            connection.execute(
                """INSERT INTO exports(
                export_id, output_dir, format, content_hash, final_artifact_ids_json, created_at
                ) VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    export_id,
                    str(Path(output_dir).resolve()),
                    format_name,
                    content_hash,
                    _canonical_json(tuple(sorted(final_artifact_ids))),
                    utc_now(),
                ),
            )
            connection.commit()
        return export_id
