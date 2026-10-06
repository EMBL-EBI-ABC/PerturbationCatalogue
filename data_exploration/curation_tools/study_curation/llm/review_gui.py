"""Operational dashboard for SQLite curation runs.

The dashboard retains the practical inspection and selection tools of the
original control centre, but all status tables and payload viewers query the
named run's SQLite state instead of intermediate workflow files.
"""

from __future__ import annotations

import uuid
from pathlib import Path
from typing import Callable

import pandas as pd
import streamlit as st

from curation_tools.study_curation.llm.curation_run_store import (
    DATABASE_FILENAME,
    CurationRunConflictError,
)
from curation_tools.study_curation.llm.source_ingestion import build_create_run_request
from curation_tools.study_curation.llm.workflow import CurationWorkflow
from curation_tools.study_curation.paths import (
    CURATION_RUNS_DIR,
    FULL_TEXT_MD_DIR,
    MAVEDB_METADATA_OUTPUT_DIR,
    MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
)

STEP_LABELS = {
    "step1": "Extract publication evidence",
    "step2": "Normalize controlled terms",
    "step3": "Discover ontology candidates",
    "step4": "Apply approved term mappings",
    "step5": "Assemble final metadata",
}
LOGO_PATH = Path(__file__).with_name("assets") / "CurateLab.svg"


def _run_names(runs_root: Path) -> list[str]:
    if not runs_root.is_dir():
        return []
    return sorted(
        path.name
        for path in runs_root.iterdir()
        if path.is_dir() and (path / DATABASE_FILENAME).is_file()
    )


def _session_id() -> str:
    if "curation_session_id" not in st.session_state:
        st.session_state.curation_session_id = str(uuid.uuid4())
    return st.session_state.curation_session_id


def _parse_urns(value: str) -> list[str]:
    return [urn.strip() for urn in value.replace(",", "\n").splitlines() if urn.strip()]


def _create_run_panel(workflow: CurationWorkflow) -> None:
    with st.expander("Create a new run", expanded=False):
        st.caption(
            "Source caches are read once here. Prompts, model settings, schema source, "
            "Markdown, and selected MaveDB records are then snapshotted into SQLite."
        )
        with st.form("create-run"):
            run_name = st.text_input("Run name")
            publication_dir = st.text_input(
                "Publication Markdown directory", value=str(FULL_TEXT_MD_DIR)
            )
            mavedb_dir = st.text_input(
                "MaveDB metadata directory", value=str(MAVEDB_METADATA_OUTPUT_DIR)
            )
            mapping_file = st.text_input(
                "URN-to-DOI mapping JSON", value=str(MAVEDB_URN_TO_DOIS_OUTPUT_FILE)
            )
            left, right = st.columns(2)
            model = left.text_input("Model", value="google/gemini-3.7-flash")
            max_workers = right.number_input(
                "Parallel LLM workers", min_value=1, max_value=64, value=8, step=1
            )
            selected_urns = st.text_area(
                "MaveDB URIs to include (optional; one per line or comma-separated)",
                help="Leave blank to create items for every linked Markdown publication.",
            )
            submitted = st.form_submit_button("Snapshot inputs and create run")
        if not submitted:
            return
        try:
            request = build_create_run_request(
                run_name=run_name,
                publication_dir=publication_dir,
                mavedb_metadata_dir=mavedb_dir,
                urn_to_dois_file=mapping_file,
                model=model,
                max_workers=int(max_workers),
                selected_urns=_parse_urns(selected_urns),
            )
            summary = workflow.create_run(request)
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            st.error(str(exc))
            return
        st.success(
            f"Created {summary.run_name}; all selected inputs are now persisted."
        )
        st.session_state.selected_native_run = summary.run_name
        st.rerun()


def _item_table_rows(
    workflow: CurationWorkflow, run_name: str
) -> list[dict[str, object]]:
    return [
        {
            **row,
            "mavedb_uris": "\n".join(row["source_urns"]),
            "publication_sources": "\n".join(row["source_documents"]),
        }
        for row in workflow.inspect_items(run_name)
    ]


def _select_items(
    rows: list[dict[str, object]],
    key: str,
    eligible: Callable[[dict[str, object]], bool],
) -> list[str]:
    """Render the curation queue using stable item IDs."""
    all_urns = sorted({urn for row in rows for urn in row["source_urns"]})
    uri_filter = st.multiselect("Filter MaveDB URIs", all_urns, key=f"{key}-uri-filter")
    displayed = [
        row
        for row in rows
        if not uri_filter or any(urn in uri_filter for urn in row["source_urns"])
    ]
    editable = []
    for row in displayed:
        can_run = eligible(row)
        editable.append(
            {
                "run": can_run,
                "eligible": can_run,
                "dataset_id": row["dataset_id"],
                "MaveDB URI(s)": row["mavedb_uris"],
                "publication source(s)": row["publication_sources"],
                "Step 1": row["evidence"],
                "Step 2": row["normalized"],
                "Step 4": row["backfilled"],
                "Step 5": row["final"],
                "item_id": row["item_id"],
            }
        )
    if not editable:
        st.info("No curation items match this filter.")
        return []
    edited = st.data_editor(
        pd.DataFrame(editable),
        hide_index=True,
        disabled=[column for column in editable[0] if column not in {"run"}],
        column_config={
            "run": st.column_config.CheckboxColumn("Run", help="Queue this item"),
            "MaveDB URI(s)": st.column_config.TextColumn(width="large"),
        },
        key=f"{key}-selection-table",
        use_container_width=True,
    )
    return edited.loc[(edited["run"]) & (edited["eligible"]), "item_id"].tolist()


def _artifact_inspector(workflow: CurationWorkflow, run_name: str, kind: str) -> None:
    artifacts = workflow.inspect_artifacts(run_name, kind)
    if not artifacts:
        st.info(f"No persisted {kind} artifacts yet.")
        return
    table = pd.DataFrame(
        [
            {
                "artifact_id": artifact.artifact_id,
                "dataset_id": artifact.dataset_id,
                "MaveDB URI(s)": "\n".join(artifact.source_urns),
                "publication source(s)": "\n".join(artifact.source_documents),
                "payload fields": len(artifact.payload),
            }
            for artifact in artifacts
        ]
    )
    st.dataframe(table, hide_index=True, use_container_width=True)
    by_id = {artifact.artifact_id: artifact for artifact in artifacts}
    selected_id = st.selectbox(
        f"Inspect {kind} payload",
        options=list(by_id),
        format_func=lambda value: (
            f"{by_id[value].dataset_id or by_id[value].item_label} — {value[:8]}"
        ),
        key=f"inspect-{kind}",
    )
    selected = by_id[selected_id]
    left, right = st.columns(2)
    with left:
        st.caption("Persisted payload")
        st.json(dict(selected.payload))
    with right:
        st.caption("Prompt/model/source provenance")
        st.json(dict(selected.metadata))


def _execute_step(
    workflow: CurationWorkflow,
    run_name: str,
    step: str,
    selected_item_ids: list[str],
    active: bool,
) -> None:
    result_key = f"{run_name}-{step}-last-result"
    previous = st.session_state.pop(result_key, None)
    if previous:
        message, is_failure = previous
        (st.error if is_failure else st.success)(message)
    if st.button(
        f"Execute {STEP_LABELS[step]}",
        type="primary",
        disabled=not active or not selected_item_ids,
        key=f"execute-{step}",
    ):
        result = workflow.run_step(run_name, step, selected_item_ids)
        if result.state == "failed":
            message = "; ".join(result.failures)
            st.session_state[result_key] = (message, True)
        elif step == "step3" and result.state == "skipped":
            st.session_state[result_key] = (
                "No new 'Other' terms with supporting evidence require candidate discovery.",
                False,
            )
        else:
            message = (
                f"{STEP_LABELS[step]}: {result.state}; "
                f"processed {len(result.processed_item_ids)}, skipped {len(result.skipped_item_ids)}."
            )
            st.session_state[result_key] = (message, False)
        st.rerun()


def _step_tab(
    workflow: CurationWorkflow,
    run_name: str,
    rows: list[dict[str, object]],
    step: str,
    kind: str,
    eligible: Callable[[dict[str, object]], bool],
    description: str,
    active: bool,
) -> None:
    st.header(STEP_LABELS[step])
    st.write(description)
    selected = _select_items(rows, step, eligible)
    st.caption(f"{len(selected)} eligible item(s) selected.")
    _execute_step(workflow, run_name, step, selected, active)
    with st.expander(f"Inspect persisted {kind} results", expanded=True):
        _artifact_inspector(workflow, run_name, kind)


def _review_panel(
    workflow: CurationWorkflow, run_name: str, item_labels: dict[str, str], active: bool
) -> None:
    review = workflow.load_review(run_name)
    st.header("Review candidates and update schema")
    if not review.candidates:
        st.info("Run Step 3 to create persisted ontology candidates.")
        return
    summary = pd.DataFrame(
        [
            {
                "field": row["field_name"],
                "proposed term": row["term"],
                "status": row["status"],
                "supporting datasets": len(row["supporting_evidence"]),
            }
            for row in review.candidates
        ]
    )
    st.dataframe(summary, hide_index=True, use_container_width=True)
    changes = []
    with st.form("candidate-review"):
        for row in review.candidates:
            st.markdown(f"#### `{row['field_name']}` — `{row['term']}`")
            left, middle, right = st.columns((2, 3, 2))
            decision = left.selectbox(
                "Global decision",
                ("Pending", "Approved", "Rejected", "Accepted"),
                index=("Pending", "Approved", "Rejected", "Accepted").index(
                    row["status"]
                ),
                key=f"status-{row['candidate_id']}",
                disabled=not active,
            )
            term = middle.text_input(
                "Effective term",
                value=row["term"] or "",
                key=f"term-{row['candidate_id']}",
                disabled=not active,
            )
            is_null = right.checkbox(
                "Explicit null",
                value=row["is_null"],
                key=f"null-{row['candidate_id']}",
                disabled=not active,
            )
            with st.expander("Supporting evidence and dataset overrides"):
                st.dataframe(
                    pd.DataFrame(row["supporting_evidence"]), use_container_width=True
                )
                for item_id in sorted(
                    {
                        str(evidence["item_id"])
                        for evidence in row["supporting_evidence"]
                        if evidence.get("item_id")
                    }
                ):
                    current = row["item_decisions"].get(
                        item_id,
                        {
                            "event_id": None,
                            "status": "Pending",
                            "term": row["term"],
                            "is_null": False,
                            "rationale": row["rationale"],
                        },
                    )
                    st.caption(f"Dataset override: {item_labels.get(item_id, item_id)}")
                    item_left, item_middle, item_right = st.columns((2, 3, 2))
                    override_status = item_left.selectbox(
                        "Decision",
                        ("Pending", "Approved", "Rejected", "Accepted"),
                        index=("Pending", "Approved", "Rejected", "Accepted").index(
                            current["status"]
                        ),
                        key=f"item-status-{row['candidate_id']}-{item_id}",
                        disabled=not active,
                    )
                    override_term = item_middle.text_input(
                        "Term",
                        value=current["term"] or "",
                        key=f"item-term-{row['candidate_id']}-{item_id}",
                        disabled=not active,
                    )
                    override_null = item_right.checkbox(
                        "Explicit null",
                        value=current["is_null"],
                        key=f"item-null-{row['candidate_id']}-{item_id}",
                        disabled=not active,
                    )
                    if (
                        override_status != current["status"]
                        or (override_term or None) != current["term"]
                        or override_null != current["is_null"]
                    ):
                        changes.append(
                            {
                                "candidate_id": row["candidate_id"],
                                "scope_type": "item",
                                "scope_key": item_id,
                                "status": override_status,
                                "term": override_term or None,
                                "is_null": override_null,
                                "rationale": current["rationale"],
                                "expected_event_id": current["event_id"],
                            }
                        )
            if (
                decision != row["status"]
                or (term or None) != row["term"]
                or is_null != row["is_null"]
            ):
                changes.append(
                    {
                        "candidate_id": row["candidate_id"],
                        "scope_type": "global",
                        "status": decision,
                        "term": term or None,
                        "is_null": is_null,
                        "rationale": row["rationale"],
                        "expected_event_id": row["event_id"],
                    }
                )
        submitted = st.form_submit_button("Commit review changes", disabled=not active)
    if submitted:
        try:
            workflow.commit_review_changes(run_name, changes, None, _session_id())
        except CurationRunConflictError as exc:
            st.error(f"A concurrent review changed these rows: {exc.conflicts}")
        else:
            st.success("Review events committed.")
            st.rerun()

    preview = workflow.preview_schema_update(run_name)
    st.subheader("Schema diff")
    st.json(preview["approved_terms"])
    st.code(preview["diff"] or "No schema changes are pending.", language="diff")
    if st.button("Apply approved terms to schema", disabled=not active):
        try:
            operation = workflow.apply_schema_terms(run_name, preview["schema_hash"])
        except RuntimeError as exc:
            st.error(str(exc))
        else:
            st.success(f"Schema operation {operation.operation_id} applied.")
            st.rerun()


def _step3_tab(
    workflow: CurationWorkflow,
    run_name: str,
    rows: list[dict[str, object]],
    active: bool,
) -> None:
    st.header(STEP_LABELS["step3"])
    st.write(
        "Find reusable ontology candidates from persisted Step 1 evidence and Step 2 ‘Other’ values."
    )
    selected = _select_items(rows, "step3", lambda row: row["normalized"] == "ready")
    st.caption(f"{len(selected)} normalized item(s) selected.")
    _execute_step(workflow, run_name, "step3", selected, active)
    review = workflow.load_review(run_name)
    if review.candidates:
        st.subheader("Persisted candidate table")
        st.dataframe(
            pd.DataFrame(
                [
                    {
                        "field": row["field_name"],
                        "term": row["term"],
                        "status": row["status"],
                        "evidence": len(row["supporting_evidence"]),
                    }
                    for row in review.candidates
                ]
            ),
            hide_index=True,
            use_container_width=True,
        )
    else:
        st.info("No persisted candidates yet.")


def _step4_tab(
    workflow: CurationWorkflow,
    run_name: str,
    rows: list[dict[str, object]],
    active: bool,
) -> None:
    st.header(STEP_LABELS["step4"])
    preview = workflow.preview_backfill(run_name)
    st.caption("This preview uses the exact SQLite resolver used when Step 4 commits.")
    if preview:
        st.dataframe(
            pd.DataFrame([change.__dict__ for change in preview]),
            hide_index=True,
            use_container_width=True,
        )
    else:
        st.info("No approved replacements are currently pending.")
    selected = _select_items(rows, "step4", lambda row: row["normalized"] == "ready")
    _execute_step(workflow, run_name, "step4", selected, active)
    with st.expander("Inspect persisted backfilled results", expanded=True):
        _artifact_inspector(workflow, run_name, "backfilled")


def _step5_tab(
    workflow: CurationWorkflow,
    run_name: str,
    rows: list[dict[str, object]],
    active: bool,
) -> None:
    st.header(STEP_LABELS["step5"])
    st.write(
        "Project the effective Step 4 or Step 2 payload and stored MaveDB snapshots into ObsSchema."
    )
    selected = _select_items(
        rows,
        "step5",
        lambda row: row["normalized"] == "ready" or row["backfilled"] == "ready",
    )
    _execute_step(workflow, run_name, "step5", selected, active)
    final_rows = workflow.load_final_metadata(run_name)
    if final_rows:
        st.subheader("Editable effective final metadata")
        st.dataframe(
            pd.DataFrame(final_rows), hide_index=True, use_container_width=True
        )
        with st.form("final-edit"):
            artifact_id = st.selectbox(
                "Final artifact", [row["artifact_id"] for row in final_rows]
            )
            field_name = st.text_input("Field name")
            new_value = st.text_input("New value")
            submitted = st.form_submit_button("Append final edit", disabled=not active)
        if submitted:
            if not field_name:
                st.error("Field name is required.")
            else:
                workflow.save_final_edits(
                    run_name,
                    [
                        {
                            "artifact_id": artifact_id,
                            "field_name": field_name,
                            "new_value": new_value,
                        }
                    ],
                    _session_id(),
                )
                st.rerun()
    with st.expander("Inspect final artifact provenance"):
        _artifact_inspector(workflow, run_name, "final")


def _overview_tab(workflow: CurationWorkflow, run_name: str, status, rows) -> None:
    st.header("Run overview")
    configuration = workflow.run_configuration(run_name)
    model_settings = configuration["model_settings"]
    first, second, third = st.columns(3)
    first.metric("Run state", status.summary.state)
    second.metric("Model", model_settings.get("model", ""))
    third.metric("Parallel LLM workers", model_settings.get("max_workers", 1))
    st.subheader("All persisted curation items")
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.subheader("Step execution history")
    st.dataframe(
        pd.DataFrame(status.executions), hide_index=True, use_container_width=True
    )
    st.subheader("Structured run events")
    st.dataframe(pd.DataFrame(status.events), hide_index=True, use_container_width=True)
    with st.expander("Immutable run configuration"):
        st.json(configuration)


def _completion_controls(workflow: CurationWorkflow, run_name: str, status) -> None:
    if status.summary.state == "active":
        failed = [entry for entry in status.executions if entry["state"] == "failed"]
        left, right = st.columns(2)
        if left.button("Resume latest failed execution", disabled=not bool(failed)):
            result = workflow.run_step(run_name, failed[-1]["step"], retry_failed=True)
            st.write(result)
            st.rerun()
        if right.button(
            "Seal run", disabled=not bool(workflow.load_final_metadata(run_name))
        ):
            try:
                workflow.seal_run(run_name)
            except RuntimeError as exc:
                st.error(str(exc))
            else:
                st.success(
                    "Run sealed. Only inspection and explicit export remain available."
                )
                st.rerun()
    else:
        export_dir = st.text_input(
            "Final deliverable directory",
            value=str(workflow.runs_root / run_name / "exports"),
        )
        overwrite = st.checkbox("Overwrite non-empty export directory")
        if st.button("Export final JSON and CSV"):
            try:
                result = workflow.export_final(run_name, export_dir, overwrite)
            except (FileExistsError, RuntimeError, ValueError) as exc:
                st.error(str(exc))
            else:
                st.success(
                    f"Exported {len(result.json_paths)} studies to {result.output_dir}"
                )


def main() -> None:
    st.set_page_config(page_title="CurateLab", layout="wide")
    st.image(str(LOGO_PATH), width=360)
    st.subheader("Experimental metadata curation for the Perturbation Catalogue.")
    runs_root = Path(
        st.sidebar.text_input("Runs directory", value=str(CURATION_RUNS_DIR))
    ).resolve()
    workflow = CurationWorkflow(runs_root)
    _create_run_panel(workflow)
    run_names = _run_names(runs_root)
    if not run_names:
        st.info(
            "Create a run to begin. JSON-era folders are intentionally unsupported."
        )
        return
    current = st.sidebar.selectbox(
        "Run", run_names, key="selected_native_run", index=None
    )
    if not current:
        return
    status = workflow.get_status(current)
    rows = _item_table_rows(workflow, current)
    active = status.summary.state == "active"
    st.caption(
        f"{status.summary.state} · revision {status.summary.revision} · {len(rows)} curation items"
    )
    tabs = st.tabs(
        (
            "📋 Overview",
            "📚 Extract evidence",
            "🏷️ Normalize terms",
            "💡 Discover candidates",
            "🔍 Review & update schema",
            "🔄 Apply approved terms",
            "🧬 Assemble metadata",
            "✅ Seal & export",
        )
    )
    with tabs[0]:
        _overview_tab(workflow, current, status, rows)
    with tabs[1]:
        _step_tab(
            workflow,
            current,
            rows,
            "step1",
            "evidence",
            lambda row: row["evidence"] == "pending",
            "Extract verbatim evidence from the run's snapshotted publication text and MaveDB context.",
            active,
        )
    with tabs[2]:
        _step_tab(
            workflow,
            current,
            rows,
            "step2",
            "normalized",
            lambda row: row["evidence"] == "ready" and row["normalized"] == "pending",
            "Map extracted evidence to the run's snapshotted controlled vocabulary schema.",
            active,
        )
    with tabs[3]:
        _step3_tab(workflow, current, rows, active)
    with tabs[4]:
        _review_panel(
            workflow,
            current,
            {row["item_id"]: str(row["item_label"]) for row in rows},
            active,
        )
    with tabs[5]:
        _step4_tab(workflow, current, rows, active)
    with tabs[6]:
        _step5_tab(workflow, current, rows, active)
    with tabs[7]:
        _completion_controls(workflow, current, status)


if __name__ == "__main__":
    main()
