"""Streamlit GUI for reviewing Step 3 ontology candidates, applying approved terms to SpecificTermExtractionSchema, and executing Step 4 backfill."""

import json
from pathlib import Path
import pandas as pd
import streamlit as st

from curation_tools.llm_curation.backfill_terms import (
    backfill_approved_terms,
    preview_backfill_changes,
)
from curation_tools.llm_curation.schema_updater import (
    DEFAULT_SCHEMA_PATH,
    apply_schema_update,
    get_existing_literals,
    get_schema_diff,
)

DEFAULT_CANDIDATES_PATH = (
    Path.cwd()
    / "test_output"
    / "step3_ontology_candidates"
    / "step3_ontology_candidates.json"
)

st.set_page_config(
    page_title="Step 3 & 4 Ontology Candidate Curator",
    page_icon="🧬",
    layout="wide",
)


def load_candidates_file(file_path: Path) -> dict:
    """Load JSON file containing Step 3 ontology candidates."""
    if not file_path.exists():
        return {}
    try:
        return json.loads(file_path.read_text(encoding="utf-8"))
    except Exception as e:
        st.error(f"Error reading {file_path}: {e}")
        return {}


def init_session_state(
    candidates_data: dict, candidates_path: Path | None = None
) -> None:
    """Initialize decision state for candidates, syncing with existing schema allowed values and prior audit logs."""
    if "candidates_raw" not in st.session_state:
        st.session_state["candidates_raw"] = candidates_data

    if "decisions" not in st.session_state:
        existing_literals = get_existing_literals(DEFAULT_SCHEMA_PATH)

        # Check for prior decision audit file
        saved_audit = {}
        if candidates_path:
            audit_path = candidates_path.parent / "approved_ontology_terms.json"
            if audit_path.exists():
                try:
                    saved_audit = json.loads(audit_path.read_text(encoding="utf-8"))
                except Exception:
                    saved_audit = {}

        decisions = {}
        for field, cand_list in candidates_data.items():
            decisions[field] = []
            allowed_vocab = set(existing_literals.get(field, []))

            # Build map of saved field decisions
            saved_field_map = {}
            if field in saved_audit and isinstance(saved_audit[field], list):
                for saved_entry in saved_audit[field]:
                    orig_k = saved_entry.get("original_term") or saved_entry.get("term")
                    if orig_k:
                        saved_field_map[orig_k] = saved_entry

            for item in cand_list:
                term = item.get("proposed_new_term") or item.get("proposed_label") or ""
                saved_entry = saved_field_map.get(term)

                # Determine effective term label and whether it is already present in the schema
                display_term = saved_entry.get("term", term) if saved_entry else term
                is_in_schema = term in allowed_vocab or display_term in allowed_vocab

                if is_in_schema:
                    status = "Approved"
                elif saved_entry:
                    status = saved_entry.get("status", "Pending")
                else:
                    status = "Pending"

                decisions[field].append(
                    {
                        "term": display_term,
                        "original_term": term,
                        "status": status,
                        "rationale": item.get("rationale", ""),
                        "supporting_evidence": item.get("supporting_evidence", []),
                    }
                )
        st.session_state["decisions"] = decisions


def main():
    st.title("🧬 Step 3 & 4 Ontology Candidate Curator & Backfill")
    st.caption(
        "Review proposed terms, update SpecificTermExtractionSchema, and backfill 'Other' values into Step 4 copies."
    )

    # Sidebar setup
    st.sidebar.header("📁 Data Source")
    file_path_str = st.sidebar.text_input(
        "Candidates JSON Path",
        value=str(DEFAULT_CANDIDATES_PATH),
        help="Path to step3_ontology_candidates.json",
    )
    candidates_path = Path(file_path_str).resolve()

    if st.sidebar.button("Reload Candidates & Refresh Schema"):
        st.session_state.clear()
        st.rerun()

    candidates_data = load_candidates_file(candidates_path)
    if not candidates_data:
        st.warning(
            f"No candidate data found at `{candidates_path}`. Please verify the file path."
        )
        st.stop()

    init_session_state(candidates_data, candidates_path=candidates_path)
    decisions = st.session_state["decisions"]

    # Compute overall statistics
    total_fields = len(decisions)
    total_candidates = sum(len(c_list) for c_list in decisions.values())
    approved_count = sum(
        1 for c_list in decisions.values() for c in c_list if c["status"] == "Approved"
    )
    rejected_count = sum(
        1 for c_list in decisions.values() for c in c_list if c["status"] == "Rejected"
    )
    pending_count = sum(
        1 for c_list in decisions.values() for c in c_list if c["status"] == "Pending"
    )

    st.sidebar.subheader("📊 Metrics")
    col_m1, col_m2 = st.sidebar.columns(2)
    col_m1.metric("Total Fields", total_fields)
    col_m2.metric("Total Candidates", total_candidates)

    col_m3, col_m4, col_m5 = st.sidebar.columns(3)
    col_m3.metric("Approved", approved_count)
    col_m4.metric("Rejected", rejected_count)
    col_m5.metric("Pending", pending_count)

    st.sidebar.divider()

    # Field selection navigation
    field_options = list(decisions.keys())
    field_labels = [f"{field} ({len(decisions[field])})" for field in field_options]
    selected_field_idx = st.sidebar.selectbox(
        "Select Field to Review",
        options=range(len(field_options)),
        format_func=lambda i: field_labels[i],
    )
    selected_field = field_options[selected_field_idx]

    status_filter = st.sidebar.radio(
        "Filter Candidates Status",
        options=["All", "Pending", "Approved", "Rejected"],
        index=0,
    )

    # Tabs
    tab_review, st_tab_diff, tab_backfill = st.tabs(
        [
            "🔍 Candidate Review",
            "📝 Approved Terms & Schema Diff",
            "🔄 Step 4: Backfill Approved Terms",
        ]
    )

    with tab_review:
        st.header(f"Field: `{selected_field}`")

        # Display existing literals for context
        existing_literals = get_existing_literals(DEFAULT_SCHEMA_PATH).get(
            selected_field, []
        )
        if existing_literals:
            with st.expander(
                f"📋 Current Allowed Vocabulary ({len(existing_literals)} terms)",
                expanded=False,
            ):
                st.write(", ".join([f"`{term}`" for term in existing_literals]))

        # Field-level bulk action buttons
        col_b1, col_b2, col_b3, _ = st.columns([1, 1, 1, 3])
        if col_b1.button("✅ Approve All in Field"):
            for item in decisions[selected_field]:
                item["status"] = "Approved"
            st.rerun()
        if col_b2.button("❌ Reject All in Field"):
            for item in decisions[selected_field]:
                item["status"] = "Rejected"
            st.rerun()
        if col_b3.button("🔄 Reset Field to Pending"):
            for item in decisions[selected_field]:
                item["status"] = "Pending"
            st.rerun()

        st.divider()

        # Render candidates for selected field
        candidates_to_render = []
        for idx, item in enumerate(decisions[selected_field]):
            if status_filter == "All" or item["status"] == status_filter:
                candidates_to_render.append((idx, item))

        if not candidates_to_render:
            st.info(
                f"No candidates matching filter status '{status_filter}' for `{selected_field}`."
            )

        for idx, item in candidates_to_render:
            with st.container(border=True):
                col_term, col_actions = st.columns([3, 2])

                with col_term:
                    new_term_val = st.text_input(
                        "Proposed Term Label",
                        value=item["term"],
                        key=f"term_input_{selected_field}_{idx}",
                        help="Edit the term label if needed before approving.",
                    )
                    item["term"] = new_term_val.strip()

                with col_actions:
                    st.write("**Decision:**")
                    btn_a, btn_r, btn_p = st.columns(3)

                    status = item["status"]

                    if btn_a.button(
                        "✅ Approve",
                        key=f"app_{selected_field}_{idx}",
                        type="primary" if status == "Approved" else "secondary",
                    ):
                        item["status"] = "Approved"
                        st.rerun()
                    if btn_r.button(
                        "❌ Reject",
                        key=f"rej_{selected_field}_{idx}",
                        type="primary" if status == "Rejected" else "secondary",
                    ):
                        item["status"] = "Rejected"
                        st.rerun()
                    if btn_p.button("↩️ Reset", key=f"rst_{selected_field}_{idx}"):
                        item["status"] = "Pending"
                        st.rerun()

                # Status tag
                if item["status"] == "Approved":
                    st.success("Status: Approved (in schema or marked for update)")
                elif item["status"] == "Rejected":
                    st.error("Status: Rejected")
                else:
                    st.info("Status: Pending Review")

                if item["rationale"]:
                    st.markdown(f"**Rationale:** {item['rationale']}")

                evidence = item["supporting_evidence"]
                if evidence:
                    with st.expander(
                        f"💬 Supporting Evidence Snippets ({len(evidence)})"
                    ):
                        for ev_idx, ev_item in enumerate(evidence):
                            if isinstance(ev_item, dict):
                                stmt = (
                                    ev_item.get("evidence_statement")
                                    or ev_item.get("evidence")
                                    or ""
                                )
                                src = ev_item.get("source_file", "unknown")
                                st.markdown(f"**{ev_idx + 1}. Source:** `{src}`")
                                st.caption(f'> "{stmt}"')
                            else:
                                st.caption(f'> "{ev_item}"')

    with st_tab_diff:
        st.header("📝 Approved Terms & Schema Diff")

        # Build approved terms mapping
        approved_map: dict[str, list[str]] = {}
        for field, c_list in decisions.items():
            approved_terms = [
                c["term"] for c in c_list if c["status"] == "Approved" and c["term"]
            ]
            if approved_terms:
                approved_map[field] = approved_terms

        if not approved_map:
            st.info(
                "No candidates have been approved yet. Switch to the 'Candidate Review' tab and approve terms to see proposed schema changes."
            )
        else:
            st.subheader("Summary of Approved Terms")
            for field, terms in approved_map.items():
                st.write(f"- **`{field}`**: " + ", ".join([f"`{t}`" for t in terms]))

            st.divider()
            st.subheader("Unified Code Diff (`llm_curation_schema.py`)")

            try:
                diff_str = get_schema_diff(approved_map, DEFAULT_SCHEMA_PATH)
                if diff_str.strip():
                    st.code(diff_str, language="diff")
                else:
                    st.info(
                        "All currently approved terms are already merged into the schema."
                    )
            except Exception as e:
                st.error(f"Error computing diff: {e}")

            st.divider()

            col_apply, col_save_json = st.columns(2)

            create_backup = col_apply.checkbox(
                "Create `.py.bak` backup before applying", value=True
            )
            if col_apply.button(
                "🚀 Apply Approved Terms to SpecificTermExtractionSchema",
                type="primary",
            ):
                try:
                    applied_diff = apply_schema_update(
                        approved_map,
                        DEFAULT_SCHEMA_PATH,
                        create_backup=create_backup,
                    )
                    # Automatically export audit log
                    output_path = (
                        candidates_path.parent / "approved_ontology_terms.json"
                    )
                    output_path.write_text(
                        json.dumps(decisions, indent=2), encoding="utf-8"
                    )

                    st.success(
                        "Successfully updated `SpecificTermExtractionSchema` in `llm_curation_schema.py` and saved audit log!"
                    )
                    st.balloons()
                except Exception as ex:
                    st.error(f"Failed to update schema: {ex}")

            if col_save_json.button("💾 Export Decision Audit Trail (JSON)"):
                output_path = candidates_path.parent / "approved_ontology_terms.json"
                output_path.write_text(
                    json.dumps(decisions, indent=2), encoding="utf-8"
                )
                st.success(f"Saved decisions to `{output_path}`")

    with tab_backfill:
        st.header("🔄 Step 4: Backfill Approved Terms")
        st.caption(
            "Visual preview of 'Other' replacements that will be applied to Step 2 copies when executing Step 4."
        )

        default_step2_dir = (
            candidates_path.parent.parent / "step2_normalized"
            if candidates_path.parent.parent.joinpath("step2_normalized").exists()
            else Path.cwd() / "test_output" / "step2_normalized"
        )
        default_step4_dir = (
            candidates_path.parent.parent / "step4_backfilled"
            if candidates_path.parent.parent.exists()
            else Path.cwd() / "test_output" / "step4_backfilled"
        )

        col_b_in, col_b_out = st.columns(2)
        step2_dir_input = col_b_in.text_input(
            "Step 2 Source Directory",
            value=str(default_step2_dir),
            help="Directory containing Step 2 normalized JSON files.",
        )
        step4_dir_input = col_b_out.text_input(
            "Step 4 Output Directory",
            value=str(default_step4_dir),
            help="Directory where backfilled Step 4 copies will be written.",
        )

        decisions_file_path = candidates_path.parent / "approved_ontology_terms.json"

        # Compute preview of planned changes BEFORE writing
        preview_records = preview_backfill_changes(
            step2_dir=Path(step2_dir_input),
            decisions_data=decisions,
        )

        st.divider()
        st.subheader(
            f"👀 Planned Backfill Changes Preview ({len(preview_records)} replacements planned)"
        )

        if preview_records:
            df_preview = pd.DataFrame(preview_records)

            # Preview Metrics
            unique_datasets = df_preview["dataset_id"].nunique()
            unique_files = df_preview["source_file"].nunique()
            unique_fields = df_preview["field_name"].nunique()

            col_p1, col_p2, col_p3, col_p4 = st.columns(4)
            col_p1.metric("Planned Replacements", len(preview_records))
            col_p2.metric("Datasets Affected", unique_datasets)
            col_p3.metric("Files Affected", unique_files)
            col_p4.metric("Fields Modified", unique_fields)

            st.write("**Filter Preview Instances:**")
            col_pf1, col_pf2 = st.columns(2)

            field_filter_opts = ["All Fields"] + sorted(
                df_preview["field_name"].unique().tolist()
            )
            selected_audit_field = col_pf1.selectbox(
                "Filter Preview by Field Name", options=field_filter_opts
            )

            dataset_filter_opts = ["All Datasets"] + sorted(
                df_preview["dataset_id"].unique().tolist()
            )
            selected_audit_dataset = col_pf2.selectbox(
                "Filter Preview by Dataset ID", options=dataset_filter_opts
            )

            df_filtered = df_preview.copy()
            if selected_audit_field != "All Fields":
                df_filtered = df_filtered[
                    df_filtered["field_name"] == selected_audit_field
                ]
            if selected_audit_dataset != "All Datasets":
                df_filtered = df_filtered[
                    df_filtered["dataset_id"] == selected_audit_dataset
                ]

            st.dataframe(
                df_filtered,
                column_config={
                    "source_file": st.column_config.TextColumn("Source File"),
                    "dataset_id": st.column_config.TextColumn("Dataset ID"),
                    "field_name": st.column_config.TextColumn("Field Name"),
                    "previous_value": st.column_config.TextColumn("Previous Value"),
                    "new_value": st.column_config.TextColumn("New Backfilled Value"),
                },
                hide_index=True,
                use_container_width=True,
            )
        else:
            st.info(
                "No pending 'Other' replacements found in Step 2 for the currently approved terms."
            )

        st.divider()

        if st.button(
            "🚀 Execute Step 4 Backfill & Write Files",
            type="primary",
            disabled=len(preview_records) == 0,
        ):
            # Save current decisions
            decisions_file_path.parent.mkdir(parents=True, exist_ok=True)
            decisions_file_path.write_text(
                json.dumps(decisions, indent=2), encoding="utf-8"
            )

            log_file = Path(step4_dir_input).parent / "step4_backfill.log"

            try:
                copied_count, fields_updated, audit_records = backfill_approved_terms(
                    step2_dir=Path(step2_dir_input),
                    output_dir=Path(step4_dir_input),
                    decisions_file=decisions_file_path,
                    log_file=log_file,
                    create_csv=True,
                )
                st.success(
                    f"Step 4 Backfill Complete! Copied {copied_count} files and replaced {fields_updated} 'Other' values in `{step4_dir_input}`."
                )
                st.balloons()
            except Exception as e:
                st.error(f"Backfill failed: {e}")


if __name__ == "__main__":
    main()
