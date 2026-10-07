"""Run-oriented orchestration for the SQLite curation workflow.

This module is the external seam for curation work.  Callers provide inputs at
creation time and subsequently address a run by name; intermediate workflow
state never crosses the seam as filesystem paths or JSON files.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence, Type

from pydantic import BaseModel

from curation_tools.study_curation.llm.curation_run_store import (
    Artifact,
    CurationRunError,
    CurationRunStore,
    ReviewChange,
    RunSummary,
    SchemaOperation,
)
from curation_tools.study_curation.llm.final_metadata import (
    merge_final_metadata_from_entries,
)
from curation_tools.study_curation.llm.llm_curation_schema import (
    EvidenceExtractionSchema,
    FieldCandidates,
    SpecificTermExtractionSchema,
)
from curation_tools.study_curation.paths import CURATION_RUNS_DIR
from curation_tools.study_curation.sources.mavedb import (
    extract_curated_mavedb_prompt_metadata,
)

STEP_NAMES = {"step1", "step2", "step3", "step4", "step5"}


@dataclass(frozen=True)
class SourceDocumentInput:
    """Publication text and provenance to snapshot as a run source."""

    source_label: str
    content_text: str
    original_path: str | Path
    doi: str | None = None


@dataclass(frozen=True)
class MaveDBSnapshotInput:
    """MaveDB record payload and provenance to snapshot into a run."""

    urn: str
    source_label: str
    payload: Mapping[str, Any]
    original_path: str | Path


@dataclass(frozen=True)
class CurationItemInput:
    """Connect source snapshots and prompt context for one curation item."""

    item_label: str
    document_labels: tuple[str, ...]
    source_urns: tuple[str, ...]
    prompt_context: Mapping[str, Any]
    primary_dataset_id: str | None = None


@dataclass(frozen=True)
class CreateRunRequest:
    """Inputs and settings required to create a fully snapshotted run."""

    run_name: str
    documents: tuple[SourceDocumentInput, ...]
    mavedb_snapshots: tuple[MaveDBSnapshotInput, ...]
    items: tuple[CurationItemInput, ...]
    prompt_templates: Mapping[str, str]
    model_settings: Mapping[str, Any]
    schema_path: str | Path
    configuration: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class StepSummary:
    """Outcome and item-level details for one workflow step execution."""

    run: RunSummary
    step: str
    execution_id: str | None
    state: str
    processed_item_ids: tuple[str, ...]
    skipped_item_ids: tuple[str, ...]
    failures: tuple[str, ...]


@dataclass(frozen=True)
class RunStatus:
    """Current run summary together with its items, executions, and events."""

    summary: RunSummary
    items: tuple[dict[str, Any], ...]
    executions: tuple[dict[str, Any], ...]
    events: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class ReviewState:
    """Candidate review rows at the run's current revision."""

    revision: int
    candidates: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class ReviewRevision:
    """Revision and event count produced by committing review changes."""

    revision: int
    event_count: int


@dataclass(frozen=True)
class BackfillChange:
    """Proposed field change produced by the metadata backfill preview."""

    item_id: str
    artifact_id: str
    field_name: str
    old_value: Any
    new_value: Any


@dataclass(frozen=True)
class ExportResult:
    """Paths and identifier for a completed final metadata export."""

    output_dir: Path
    json_paths: tuple[Path, ...]
    csv_path: Path
    export_id: str


@dataclass(frozen=True)
class ArtifactInspection:
    """One persisted artifact rendered for a dashboard or other read-only caller."""

    artifact_id: str
    item_id: str
    item_label: str
    dataset_id: str | None
    source_urns: tuple[str, ...]
    source_documents: tuple[str, ...]
    kind: str
    payload: Mapping[str, Any]
    metadata: Mapping[str, Any]


LLMCallable = Callable[
    [str, Type[BaseModel], Mapping[str, Any]], BaseModel | Mapping[str, Any]
]


def _canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _hash(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _payload(value: BaseModel | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(value, BaseModel):
        return value.model_dump()
    if isinstance(value, Mapping):
        return dict(value)
    raise TypeError("LLM response must be a Pydantic model or mapping")


class CurationWorkflow:
    """Deep module joining pure LLM operations with durable run persistence."""

    def __init__(
        self,
        runs_root: str | Path = CURATION_RUNS_DIR,
        llm_call: LLMCallable | None = None,
    ) -> None:
        """Configure the run store location and optional LLM implementation.

        Args:
            runs_root: Directory where curation runs are persisted.
            llm_call: Callable used for model requests; defaults to the
                Instructor-backed implementation.
        """
        self.runs_root = Path(runs_root).resolve()
        self._llm_call = llm_call or self._call_instructor

    def create_run(self, request: CreateRunRequest) -> RunSummary:
        """Create and fully snapshot a fresh run before any LLM work begins."""
        if not request.documents:
            raise ValueError(
                "A curation run requires at least one publication document"
            )
        if not request.items:
            raise ValueError("A curation run requires at least one curation item")

        configuration = self._configuration_snapshot(request)
        store = CurationRunStore.create(self.runs_root, request.run_name, configuration)
        document_ids: dict[str, str] = {}
        snapshot_ids: dict[str, str] = {}
        try:
            for document in request.documents:
                document_ids[document.source_label] = store.snapshot_document(
                    document.source_label,
                    document.content_text,
                    document.original_path,
                    document.doi,
                )
            for snapshot in request.mavedb_snapshots:
                snapshot_ids[snapshot.urn] = store.snapshot_mavedb_entry(
                    snapshot.urn,
                    snapshot.source_label,
                    snapshot.payload,
                    extract_curated_mavedb_prompt_metadata(dict(snapshot.payload)),
                    snapshot.original_path,
                )
            for item in request.items:
                missing_documents = set(item.document_labels) - set(document_ids)
                missing_snapshots = set(item.source_urns) - set(snapshot_ids)
                if missing_documents or missing_snapshots:
                    raise ValueError(
                        f"Item {item.item_label!r} references unsnapshotted inputs: "
                        f"documents={sorted(missing_documents)}, urns={sorted(missing_snapshots)}"
                    )
                store.create_item(
                    item.item_label,
                    [document_ids[label] for label in item.document_labels],
                    [snapshot_ids[urn] for urn in item.source_urns],
                    item.prompt_context,
                    item.primary_dataset_id,
                )
            store.log_event("INFO", "Run created with immutable source snapshots")
        except Exception as exc:
            store.log_event(
                "ERROR", "Run creation did not finish", details={"error": str(exc)}
            )
            raise
        return store.summary()

    def run_step(
        self,
        run_ref: str,
        step: str,
        selected_item_ids: Sequence[str] | None = None,
        retry_failed: bool = False,
    ) -> StepSummary:
        """Execute missing work for one step, preserving all existing artifacts."""
        if step not in STEP_NAMES:
            raise ValueError(f"Unknown step: {step}")
        store = self._store(run_ref)
        self._require_active(store)
        selected = tuple(
            selected_item_ids or [item["item_id"] for item in store.list_items()]
        )
        eligible, skipped = self._eligible_items(store, step, selected, retry_failed)
        if not eligible:
            return StepSummary(
                store.summary(), step, None, "skipped", (), tuple(skipped), ()
            )

        retry_of = store.latest_failed_execution(step) if retry_failed else None
        execution = store.start_step(
            step,
            eligible,
            self._step_configuration(store, step),
            retry_of=retry_of,
        )
        failures: list[str] = []
        try:
            if step == "step1":
                self._run_step1(store, execution.execution_id, eligible)
            elif step == "step2":
                self._run_step2(store, execution.execution_id, eligible)
            elif step == "step3":
                self._run_step3(store, execution.execution_id, eligible)
            elif step == "step4":
                self._run_step4(store, execution.execution_id, eligible)
            else:
                self._run_step5(store, execution.execution_id, eligible)
        except Exception as exc:
            failures.append(str(exc))
            store.log_event(
                "ERROR",
                f"{step} failed",
                execution.execution_id,
                details={"error": str(exc)},
            )
            store.fail_step(execution.execution_id, str(exc))
            return StepSummary(
                store.summary(),
                step,
                execution.execution_id,
                "failed",
                (),
                tuple(skipped),
                tuple(failures),
            )
        store.finish_step(execution.execution_id)
        store.log_event("INFO", f"{step} completed", execution.execution_id)
        return StepSummary(
            store.summary(),
            step,
            execution.execution_id,
            "completed",
            tuple(eligible),
            tuple(skipped),
            (),
        )

    def get_status(self, run_ref: str) -> RunStatus:
        """Return the run summary, item rows, step history, and recent events."""
        store = self._store(run_ref)
        return RunStatus(
            store.summary(),
            tuple(store.list_items()),
            tuple(store.step_history()),
            tuple(store.recent_events()),
        )

    def run_configuration(self, run_ref: str) -> dict[str, Any]:
        """Return the immutable configuration snapshot for display only."""
        return self._store(run_ref).configuration()

    def load_final_metadata(self, run_ref: str) -> tuple[dict[str, Any], ...]:
        """Return the effective final revision, including append-only edits."""
        store = self._store(run_ref)
        return tuple(
            {
                "artifact_id": artifact.artifact_id,
                **store.effective_final_payload(artifact.artifact_id),
            }
            for artifact in store.latest_artifacts("final")
        )

    def inspect_items(self, run_ref: str) -> tuple[dict[str, Any], ...]:
        """Return one status row per stable item for the dashboard's step tables."""
        store = self._store(run_ref)
        artifacts = {
            kind: {
                artifact.item_id: artifact for artifact in store.latest_artifacts(kind)
            }
            for kind in ("evidence", "normalized", "backfilled", "final")
        }
        execution_states: dict[str, dict[str, str]] = {
            step: {} for step in ("step1", "step2", "step4", "step5")
        }
        for execution in store.step_history():
            states = execution_states.get(execution["step"])
            if states is None:
                continue
            for item_id in execution["selected_item_ids"]:
                states[item_id] = execution["state"]

        def artifact_status(kind: str, step: str, item_id: str) -> str:
            if item_id in artifacts[kind]:
                return "ready"
            execution_state = execution_states[step].get(item_id)
            if execution_state in {"failed", "running"}:
                return execution_state
            if execution_state == "completed":
                return "incomplete"
            return "pending"

        rows: list[dict[str, Any]] = []
        for item in store.list_items():
            full_item = store.load_item(item["item_id"])
            rows.append(
                {
                    "item_id": item["item_id"],
                    "item_label": item["item_label"],
                    "dataset_id": item["primary_dataset_id"],
                    "source_urns": tuple(
                        entry["urn"] for entry in full_item["mavedb_snapshots"]
                    ),
                    "source_documents": tuple(
                        document["source_label"] for document in full_item["documents"]
                    ),
                    "evidence": artifact_status("evidence", "step1", item["item_id"]),
                    "normalized": artifact_status(
                        "normalized", "step2", item["item_id"]
                    ),
                    "backfilled": artifact_status(
                        "backfilled", "step4", item["item_id"]
                    ),
                    "final": artifact_status("final", "step5", item["item_id"]),
                }
            )
        return tuple(rows)

    def inspect_artifacts(
        self, run_ref: str, kind: str
    ) -> tuple[ArtifactInspection, ...]:
        """Load current persisted artifacts for a step-result inspection table."""
        store = self._store(run_ref)
        result: list[ArtifactInspection] = []
        for artifact in store.latest_artifacts(kind):
            item = store.load_item(artifact.item_id)
            result.append(
                ArtifactInspection(
                    artifact_id=artifact.artifact_id,
                    item_id=artifact.item_id,
                    item_label=item["item_label"],
                    dataset_id=item["primary_dataset_id"],
                    source_urns=tuple(
                        entry["urn"] for entry in item["mavedb_snapshots"]
                    ),
                    source_documents=tuple(
                        document["source_label"] for document in item["documents"]
                    ),
                    kind=artifact.kind,
                    payload=artifact.payload,
                    metadata=artifact.metadata,
                )
            )
        return tuple(result)

    def load_review(self, run_ref: str) -> ReviewState:
        """Load candidate review rows and their current revision."""
        store = self._store(run_ref)
        return ReviewState(store.summary().revision, tuple(store.load_review_rows()))

    def commit_review_changes(
        self,
        run_ref: str,
        changes: Iterable[ReviewChange | Mapping[str, Any]],
        expected_revisions: Mapping[str, int | None] | None,
        session_id: str,
    ) -> ReviewRevision:
        """Persist review edits, checking supplied revisions for conflicts.

        Changes may be ``ReviewChange`` instances or mappings. Missing expected
        event identifiers are filled from ``expected_revisions`` when available.
        """
        store = self._store(run_ref)
        queued: list[ReviewChange] = []
        for change in changes:
            if isinstance(change, ReviewChange):
                queued.append(change)
                continue
            values = dict(change)
            key = self._review_key(
                str(values["candidate_id"]),
                str(values.get("scope_type", "global")),
                str(values.get("scope_key", "")),
            )
            values.setdefault(
                "expected_event_id",
                (expected_revisions or {}).get(
                    key, (expected_revisions or {}).get(values["candidate_id"])
                ),
            )
            queued.append(ReviewChange(**values))
        revision = store.commit_review_changes(queued, session_id)
        store.log_event(
            "INFO", "Review changes committed", details={"count": len(queued)}
        )
        return ReviewRevision(revision, len(queued))

    def preview_backfill(self, run_ref: str) -> list[BackfillChange]:
        """Return proposed metadata backfills without changing persisted data."""
        store = self._store(run_ref)
        return self._backfill_changes(
            store, self._source_artifacts(store, "normalized")
        )

    def apply_schema_terms(
        self, run_ref: str, expected_schema_hash: str
    ) -> SchemaOperation:
        """Apply approved schema terms if the source still matches its hash.

        Raises:
            CurationRunError: If the schema changed after it was reviewed.
        """
        store = self._store(run_ref)
        configuration = store.configuration()
        schema_path = Path(configuration["schema_path"]).resolve()
        actual_hash = hashlib.sha256(schema_path.read_bytes()).hexdigest()
        if actual_hash != expected_schema_hash:
            raise CurationRunError(
                "Schema source changed; reload the schema before applying terms"
            )
        from curation_tools.study_curation.llm.schema_updater import apply_schema_update

        operation = apply_schema_update(
            store.approved_schema_terms(), schema_path=schema_path, run_store=store
        )
        store.log_event(
            "INFO",
            "Schema terms applied",
            details={"operation_id": operation.operation_id},
        )
        return operation

    def preview_schema_update(self, run_ref: str) -> dict[str, Any]:
        """Render the schema diff using persisted approvals without applying it."""
        store = self._store(run_ref)
        configuration = store.configuration()
        schema_path = Path(configuration["schema_path"]).resolve()
        source = schema_path.read_text(encoding="utf-8")
        from curation_tools.study_curation.llm.schema_updater import get_schema_diff

        return {
            "schema_hash": hashlib.sha256(source.encode("utf-8")).hexdigest(),
            "approved_terms": store.approved_schema_terms(),
            "diff": get_schema_diff(store.approved_schema_terms(), schema_path),
        }

    def save_final_edits(
        self,
        run_ref: str,
        edits: Iterable[Mapping[str, Any]],
        session_id: str,
    ) -> int:
        """Append final field edits for an active run and return their count."""
        store = self._store(run_ref)
        self._require_active(store)
        count = 0
        for edit in edits:
            store.save_final_edit(
                str(edit["artifact_id"]),
                str(edit["field_name"]),
                edit.get("new_value"),
                session_id,
            )
            count += 1
        if count:
            store.log_event("INFO", "Final edits committed", details={"count": count})
        return count

    def seal_run(self, run_ref: str) -> RunSummary:
        """Mark a run complete and return its final summary."""
        store = self._store(run_ref)
        store.log_event("INFO", "Run sealed")
        store.complete_run()
        return store.summary()

    def export_final(
        self,
        run_ref: str,
        output_dir: str | Path,
        overwrite: bool = False,
    ) -> ExportResult:
        """Write final per-item JSON and combined CSV deliverables.

        The run must be sealed before export. A non-empty destination is kept
        unless ``overwrite`` is true.
        """
        store = self._store(run_ref)
        if store.summary().state != "completed":
            raise CurationRunError("Seal the run before exporting final deliverables")
        destination = Path(output_dir).resolve()
        if destination.exists() and any(destination.iterdir()):
            if not overwrite:
                raise FileExistsError(
                    f"Export destination is not empty: {destination}; pass overwrite=True"
                )
            protected = {
                Path(destination.anchor).resolve(),
                self.runs_root,
                Path.cwd().resolve(),
            }
            if destination in protected:
                raise ValueError("Refusing to overwrite a broad export destination")
            shutil.rmtree(destination)
        destination.mkdir(parents=True, exist_ok=True)

        artifacts = store.latest_artifacts("final")
        rows: list[dict[str, Any]] = []
        json_paths: list[Path] = []
        digest = hashlib.sha256()
        for artifact in artifacts:
            payload = store.effective_final_payload(artifact.artifact_id)
            dataset_id = str(payload.get("dataset_id") or artifact.item_id)
            filename = self._safe_filename(dataset_id) + ".json"
            json_path = destination / filename
            rendered = (
                json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
            )
            json_path.write_text(rendered, encoding="utf-8")
            digest.update(rendered.encode("utf-8"))
            rows.append(payload)
            json_paths.append(json_path)
        csv_path = destination / "curation_final_metadata.csv"
        fieldnames = sorted({key for row in rows for key in row})
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(rows)
        digest.update(csv_path.read_bytes())
        export_id = store.record_export(
            destination,
            "json+csv",
            digest.hexdigest(),
            [artifact.artifact_id for artifact in artifacts],
        )
        return ExportResult(destination, tuple(json_paths), csv_path, export_id)

    def _store(self, run_ref: str) -> CurationRunStore:
        return CurationRunStore.open(self.runs_root, run_ref)

    @staticmethod
    def _require_active(store: CurationRunStore) -> None:
        if store.summary().state != "active":
            raise CurationRunError("Completed or failed runs are immutable")

    @staticmethod
    def _review_key(candidate_id: str, scope_type: str, scope_key: str) -> str:
        return "|".join((candidate_id, scope_type, scope_key))

    def _configuration_snapshot(self, request: CreateRunRequest) -> dict[str, Any]:
        schema_path = Path(request.schema_path).resolve()
        schema_source_text = schema_path.read_text(encoding="utf-8")
        prompt_templates = dict(request.prompt_templates)
        return {
            **dict(request.configuration or {}),
            "schema_path": str(schema_path),
            "schema_source_text": schema_source_text,
            "schema_source_hash": hashlib.sha256(
                schema_source_text.encode("utf-8")
            ).hexdigest(),
            "pydantic_schema_fingerprints": {
                "evidence": _hash(EvidenceExtractionSchema.model_json_schema()),
                "normalized": _hash(SpecificTermExtractionSchema.model_json_schema()),
                "candidates": _hash(FieldCandidates.model_json_schema()),
            },
            "prompt_templates": prompt_templates,
            "prompt_template_hashes": {
                key: _hash(value) for key, value in prompt_templates.items()
            },
            "model_settings": dict(request.model_settings),
        }

    def _step_configuration(self, store: CurationRunStore, step: str) -> dict[str, Any]:
        configuration = store.configuration()
        schema_key = {
            "step1": "evidence",
            "step2": "normalized",
            "step3": "candidates",
        }.get(step)
        return {
            "step": step,
            "model_settings": configuration["model_settings"],
            "prompt_hash": configuration["prompt_template_hashes"].get(step),
            "schema_fingerprint": configuration["pydantic_schema_fingerprints"].get(
                schema_key
            ),
            "schema_source_hash": configuration["schema_source_hash"],
        }

    def _eligible_items(
        self,
        store: CurationRunStore,
        step: str,
        selected: Sequence[str],
        retry_failed: bool,
    ) -> tuple[tuple[str, ...], list[str]]:
        all_items = {item["item_id"] for item in store.list_items()}
        unknown = set(selected) - all_items
        if unknown:
            raise KeyError(f"Unknown curation item IDs: {sorted(unknown)}")
        required_kind = {
            "step1": "evidence",
            "step2": "normalized",
            "step4": "backfilled",
            "step5": "final",
        }.get(step)
        if step == "step3":
            if retry_failed:
                failed_id = store.latest_failed_execution(step)
                if failed_id:
                    failed = next(
                        (
                            entry
                            for entry in store.step_history()
                            if entry["execution_id"] == failed_id
                        ),
                        None,
                    )
                    selected = failed["selected_item_ids"] if failed else selected
            other_items = self._step3_items_with_supported_other_values(store, selected)
            completed = (
                set()
                if retry_failed
                else {
                    item_id
                    for entry in store.step_history()
                    if entry["step"] == "step3" and entry["state"] == "completed"
                    for item_id in entry["selected_item_ids"]
                }
            )
            eligible = tuple(
                item_id
                for item_id in selected
                if item_id in other_items and item_id not in completed
            )
            skipped = [item_id for item_id in selected if item_id not in eligible]
            return eligible, skipped
        complete = {
            artifact.item_id
            for artifact in store.latest_artifacts(required_kind or "evidence")
        }
        if retry_failed:
            failed_id = store.latest_failed_execution(step)
            if failed_id:
                failed = next(
                    (
                        entry
                        for entry in store.step_history()
                        if entry["execution_id"] == failed_id
                    ),
                    None,
                )
                selected = failed["selected_item_ids"] if failed else selected
        eligible = tuple(item_id for item_id in selected if item_id not in complete)
        skipped = [item_id for item_id in selected if item_id in complete]
        return eligible, skipped

    @staticmethod
    def _step3_items_with_supported_other_values(
        store: CurationRunStore, item_ids: Sequence[str]
    ) -> set[str]:
        normalized = {
            artifact.item_id: artifact
            for artifact in store.latest_artifacts("normalized", item_ids)
        }
        evidence = {
            artifact.item_id: artifact
            for artifact in store.latest_artifacts("evidence", item_ids)
        }
        return {
            item_id
            for item_id, artifact in normalized.items()
            if (source := evidence.get(item_id))
            and any(
                value == "Other" and source.payload.get(f"{field_name}_evidence")
                for field_name, value in artifact.payload.items()
            )
        }

    def _run_step1(
        self, store: CurationRunStore, execution_id: str, item_ids: Sequence[str]
    ) -> None:
        configuration = store.configuration()
        if len(item_ids) > 1:
            self._parallel_items(
                configuration,
                item_ids,
                lambda item_id: self._run_step1(store, execution_id, (item_id,)),
            )
            return
        template = configuration["prompt_templates"]["step1"]
        evidence_schema = self._schema_model(store, "EvidenceExtractionSchema")
        normalization_schema = self._schema_model(store, "SpecificTermExtractionSchema")
        for item_id in item_ids:
            item = store.load_item(item_id)
            text = "\n\n".join(
                document["content_text"] for document in item["documents"]
            )
            context = self._prompt_context(item)
            prompt = template.format(
                publication_full_text=text,
                supplementary_metadata=_canonical_json(context),
                supplementary_mavedb_metadata=_canonical_json(context),
                controlled_vocabulary_hints=self._controlled_vocabulary_hints(
                    normalization_schema
                ),
            )
            result = _payload(
                self._llm_call(prompt, evidence_schema, configuration["model_settings"])
            )
            result.update(self._provenance(item))
            result["curation_agent_type"] = "LLM"
            result["curation_agent_name"] = str(
                configuration["model_settings"].get("model", "")
            )
            store.put_artifact(
                item_id,
                "evidence",
                result,
                self._artifact_metadata(store, item, prompt, "evidence"),
                execution_id,
            )
            store.log_event("INFO", "Step 1 item completed", execution_id, item_id)

    def _run_step2(
        self, store: CurationRunStore, execution_id: str, item_ids: Sequence[str]
    ) -> None:
        configuration = store.configuration()
        if len(item_ids) > 1:
            self._parallel_items(
                configuration,
                item_ids,
                lambda item_id: self._run_step2(store, execution_id, (item_id,)),
            )
            return
        template = configuration["prompt_templates"]["step2"]
        normalization_schema = self._schema_model(store, "SpecificTermExtractionSchema")
        evidence_by_item = {
            a.item_id: a for a in store.latest_artifacts("evidence", item_ids)
        }
        for item_id in item_ids:
            evidence = evidence_by_item.get(item_id)
            if evidence is None:
                raise CurationRunError(
                    f"Step 2 requires a Step 1 artifact for item {item_id}"
                )
            item = store.load_item(item_id)
            clean_evidence = {
                key: value
                for key, value in evidence.payload.items()
                if key.endswith("_evidence") or key in {"dataset_id", "data_modality"}
            }
            prompt = template.format(
                step1_evidence=json.dumps(clean_evidence, ensure_ascii=False, indent=2),
                supplementary_mavedb_metadata=json.dumps(
                    self._prompt_context(item), ensure_ascii=False, indent=2
                ),
            )
            result = _payload(
                self._llm_call(
                    prompt,
                    normalization_schema,
                    configuration["model_settings"],
                )
            )
            result.update(self._provenance(item))
            result["curation_agent_type"] = "LLM"
            result["curation_agent_name"] = str(
                configuration["model_settings"].get("model", "")
            )
            store.put_artifact(
                item_id,
                "normalized",
                result,
                self._artifact_metadata(store, item, prompt, "normalized"),
                execution_id,
                parent_artifact_id=evidence.artifact_id,
            )
            store.log_event("INFO", "Step 2 item completed", execution_id, item_id)

    def _run_step3(
        self, store: CurationRunStore, execution_id: str, item_ids: Sequence[str]
    ) -> None:
        configuration = store.configuration()
        template = configuration["prompt_templates"]["step3"]
        candidate_schema = self._schema_model(store, "FieldCandidates")
        normalization_schema = self._schema_model(store, "SpecificTermExtractionSchema")
        normalized = store.latest_artifacts("normalized", item_ids)
        evidence = {
            artifact.item_id: artifact
            for artifact in store.latest_artifacts("evidence", item_ids)
        }
        grouped: dict[str, list[dict[str, Any]]] = {}
        for artifact in normalized:
            for field_name, value in artifact.payload.items():
                if value != "Other":
                    continue
                support = evidence.get(
                    artifact.item_id, Artifact("", "", "", {}, {}, "", None, None)
                ).payload.get(f"{field_name}_evidence")
                if not support:
                    continue
                grouped.setdefault(field_name, []).append(
                    {
                        "item_id": artifact.item_id,
                        "source_file": store.load_item(artifact.item_id)["item_label"],
                        "evidence_statement": support,
                    }
                )
        candidates: dict[str, list[dict[str, Any]]] = {}
        for field_name, evidence_list in grouped.items():
            vocabulary, description = self._vocabulary(field_name, normalization_schema)
            prompt = template.format(
                field_name=field_name,
                controlled_vocabulary=(
                    f"Description: {description}\nAllowed Terms: {json.dumps(vocabulary)}"
                ),
                evidence_list=json.dumps(evidence_list, ensure_ascii=False, indent=2),
            )
            response = _payload(
                self._llm_call(
                    prompt, candidate_schema, configuration["model_settings"]
                )
            )
            entries: list[dict[str, Any]] = []
            for candidate in response.get("candidates", []):
                if not isinstance(candidate, Mapping):
                    continue
                entries.append(
                    {
                        "original_term": str(
                            candidate.get("proposed_new_term", "Other")
                        ),
                        "term": candidate.get("proposed_new_term"),
                        "rationale": candidate.get("rationale"),
                        "supporting_evidence": [
                            {**entry, "item_id": entry["item_id"]}
                            for entry in evidence_list
                        ],
                    }
                )
            candidates[field_name] = entries
        store.record_discovery(
            execution_id,
            [artifact.artifact_id for artifact in normalized],
            candidates,
        )

    def _run_step4(
        self, store: CurationRunStore, execution_id: str, item_ids: Sequence[str]
    ) -> None:
        source = self._source_artifacts(store, "normalized", item_ids)
        changes = self._backfill_changes(store, source)
        changes_by_item: dict[str, list[BackfillChange]] = {}
        for change in changes:
            changes_by_item.setdefault(change.item_id, []).append(change)
        revision_id = store.create_backfill_revision(
            execution_id, [artifact.artifact_id for artifact in source]
        )
        for artifact in source:
            updated = dict(artifact.payload)
            for change in changes_by_item.get(artifact.item_id, []):
                updated[change.field_name] = change.new_value
            store.put_artifact(
                artifact.item_id,
                "backfilled",
                updated,
                {
                    **artifact.metadata,
                    "review_revision": store.summary().revision,
                    "source_kind": "normalized",
                },
                execution_id,
                parent_artifact_id=artifact.artifact_id,
                revision_id=revision_id,
            )

    def _run_step5(
        self, store: CurationRunStore, execution_id: str, item_ids: Sequence[str]
    ) -> None:
        source = self._final_source_artifacts(store, item_ids)
        for artifact in source:
            item = store.load_item(artifact.item_id)
            payload = merge_final_metadata_from_entries(
                artifact.payload,
                [entry["payload"] for entry in item["mavedb_snapshots"]],
            )
            store.put_artifact(
                artifact.item_id,
                "final",
                payload,
                {
                    **artifact.metadata,
                    "source_kind": artifact.kind,
                    "source_artifact_id": artifact.artifact_id,
                },
                execution_id,
                parent_artifact_id=artifact.artifact_id,
            )

    def _source_artifacts(
        self,
        store: CurationRunStore,
        kind: str,
        item_ids: Sequence[str] | None = None,
    ) -> list[Artifact]:
        return store.latest_artifacts(kind, item_ids)

    @staticmethod
    def _parallel_items(
        configuration: Mapping[str, Any],
        item_ids: Sequence[str],
        work: Callable[[str], None],
    ) -> None:
        """Run independent LLM item work with the run's persisted worker limit."""
        workers = int(configuration["model_settings"].get("max_workers", 1))
        if workers <= 1:
            for item_id in item_ids:
                work(item_id)
            return
        failures: list[str] = []
        with ThreadPoolExecutor(max_workers=min(workers, len(item_ids))) as executor:
            futures = {executor.submit(work, item_id): item_id for item_id in item_ids}
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    failures.append(f"{futures[future]}: {exc}")
        if failures:
            raise CurationRunError("; ".join(failures))

    def _final_source_artifacts(
        self, store: CurationRunStore, item_ids: Sequence[str]
    ) -> list[Artifact]:
        backfilled = {
            a.item_id: a for a in store.latest_artifacts("backfilled", item_ids)
        }
        normalized = {
            a.item_id: a for a in store.latest_artifacts("normalized", item_ids)
        }
        source: list[Artifact] = []
        for item_id in item_ids:
            artifact = backfilled.get(item_id) or normalized.get(item_id)
            if artifact is None:
                raise CurationRunError(
                    f"Step 5 requires Step 2 data for item {item_id}"
                )
            source.append(artifact)
        return source

    def _backfill_changes(
        self, store: CurationRunStore, source: Sequence[Artifact]
    ) -> list[BackfillChange]:
        approved = store.resolve_approved_terms()
        changes: list[BackfillChange] = []
        for artifact in source:
            for field_name, new_value in approved.get(artifact.item_id, {}).items():
                old_value = artifact.payload.get(field_name)
                if old_value == "Other" and old_value != new_value:
                    changes.append(
                        BackfillChange(
                            artifact.item_id,
                            artifact.artifact_id,
                            field_name,
                            old_value,
                            new_value,
                        )
                    )
        return changes

    @staticmethod
    def _provenance(item: Mapping[str, Any]) -> dict[str, Any]:
        urns = [entry["urn"] for entry in item["mavedb_snapshots"]]
        source_files = [entry["source_label"] for entry in item["mavedb_snapshots"]]
        return {
            "__source_urns": urns,
            "__source_files": source_files,
            "dataset_id": item.get("primary_dataset_id") or (urns[0] if urns else None),
        }

    @staticmethod
    def _prompt_context(item: Mapping[str, Any]) -> dict[str, Any]:
        return {
            **dict(item["prompt_context"]),
            "source_urns": [entry["urn"] for entry in item["mavedb_snapshots"]],
            "curated_mavedb_metadata": [
                entry["prompt_metadata"] for entry in item["mavedb_snapshots"]
            ],
        }

    @staticmethod
    def _artifact_metadata(
        store: CurationRunStore,
        item: Mapping[str, Any],
        prompt: str,
        schema_name: str,
    ) -> dict[str, Any]:
        configuration = store.configuration()
        return {
            "prompt_text": prompt,
            "prompt_hash": _hash(prompt),
            "model_settings": configuration["model_settings"],
            "schema_fingerprint": configuration["pydantic_schema_fingerprints"][
                schema_name
            ],
            "source_hashes": {
                "documents": [
                    document["content_hash"] for document in item["documents"]
                ],
                "mavedb": [entry["payload_hash"] for entry in item["mavedb_snapshots"]],
            },
        }

    @staticmethod
    def _controlled_vocabulary_hints(
        normalization_schema: Type[BaseModel],
    ) -> str:
        lines: list[str] = []
        for name in normalization_schema.model_fields:
            values, _ = CurationWorkflow._vocabulary(name, normalization_schema)
            values = [value for value in values if value != "Other"]
            if values:
                lines.append(
                    f"- `{name}_evidence`: {', '.join(f'`{value}`' for value in values)}"
                )
        return "\n".join(lines)

    @staticmethod
    def _vocabulary(
        field_name: str, normalization_schema: Type[BaseModel]
    ) -> tuple[list[str], str]:
        import types
        import typing

        field = normalization_schema.model_fields.get(field_name)
        if field is None:
            return [], ""

        def literals(annotation: Any) -> list[str]:
            origin = typing.get_origin(annotation)
            if origin is typing.Literal:
                return [str(value) for value in typing.get_args(annotation)]
            if origin in {typing.Union, types.UnionType}:
                return [
                    value
                    for arg in typing.get_args(annotation)
                    for value in literals(arg)
                ]
            return []

        return literals(field.annotation), field.description or ""

    @staticmethod
    def _schema_model(store: CurationRunStore, name: str) -> Type[BaseModel]:
        """Load a Pydantic model from the run's immutable schema source snapshot."""
        namespace: dict[str, Any] = {"__name__": "_curation_run_schema"}
        exec(store.configuration()["schema_source_text"], namespace)
        for value in namespace.values():
            if (
                isinstance(value, type)
                and value.__module__ == "_curation_run_schema"
                and issubclass(value, BaseModel)
            ):
                value.model_rebuild(_types_namespace=namespace)
        model = namespace.get(name)
        if not isinstance(model, type) or not issubclass(model, BaseModel):
            raise CurationRunError(
                f"Run schema snapshot does not define Pydantic model: {name}"
            )
        return model

    @staticmethod
    def _safe_filename(value: str) -> str:
        cleaned = "".join(
            char if char.isalnum() or char in "._-" else "_" for char in value
        )
        return cleaned.strip("._") or "study"

    @staticmethod
    def _call_instructor(
        prompt: str, response_model: Type[BaseModel], model_settings: Mapping[str, Any]
    ) -> BaseModel:
        import instructor

        model = str(model_settings.get("model", "google/gemini-3.7-flash"))
        client = instructor.from_provider(model, location="global", vertexai=True)
        return client.create(
            response_model=response_model,
            messages=[{"role": "user", "content": prompt}],
            thinking_config=dict(
                model_settings.get("thinking_config", {"thinking_level": "high"})
            ),
            generation_config=dict(
                model_settings.get("generation_config", {"temperature": 0.2})
            ),
        )
