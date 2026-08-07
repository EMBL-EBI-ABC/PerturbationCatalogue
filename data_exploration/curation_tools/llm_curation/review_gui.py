"""Consolidated Streamlit Curation Control Center for MaveDB LLM Metadata Pipeline (Steps 1 - 4)."""

import json
from pathlib import Path
import pandas as pd
import streamlit as st

from curation_tools.llm_curation.backfill_terms import (
    backfill_approved_terms,
    preview_backfill_changes,
)
from curation_tools.llm_curation.candidate_discovery import discover_candidates
from curation_tools.llm_curation.final_metadata import (
    finalize_metadata,
    get_final_csv_path,
    load_final_csv,
    save_final_csv_edits,
)
from curation_tools.llm_curation.gui_utils import (
    calculate_file_signature,
    get_mavedb_urn_status,
    get_step2_file_status,
    get_step3_other_corpus_summary,
    get_publication_dois_for_source_file,
    is_candidate_discovery_current,
    read_last_log_lines,
    record_pipeline_step,
    resolve_effective_normalized_dir,
)
from curation_tools.llm_curation.mavedb.mavedb_metadata_extraction_runner import (
    parse_target_urns,
)
from curation_tools.llm_curation.mavedb.processing import (
    FULL_TEXT_MD_DIR,
    MAVEDB_METADATA_OUTPUT_DIR,
    MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
    bulk_extract_evidence_for_mavedb_urns,
    format_urn_for_filename,
    get_dois_from_mavedb_entry,
    load_mavedb_urn_to_dois,
)
from curation_tools.llm_curation.schema_loading import load_extraction_schema
from curation_tools.llm_curation.schema_updater import (
    DEFAULT_SCHEMA_PATH,
    apply_schema_update,
    get_existing_literals,
    get_schema_diff,
)
from curation_tools.llm_curation.specific_term_extraction import (
    normalize_evidence_artifacts,
)

# Default Test Directory Configuration
ROOT_DIR = Path.cwd()
DEFAULT_MD_DIR = ROOT_DIR / "test_md_dir"
DEFAULT_TEST_OUTPUT = ROOT_DIR / "test_output"

DEFAULT_STEP1_OUT = DEFAULT_TEST_OUTPUT / "step1_evidence"
DEFAULT_STEP2_OUT = DEFAULT_TEST_OUTPUT / "step2_normalized"
DEFAULT_STEP3_OUT = DEFAULT_TEST_OUTPUT / "step3_ontology_candidates"
DEFAULT_STEP4_OUT = DEFAULT_TEST_OUTPUT / "step4_backfilled"
DEFAULT_STEP5_OUT = DEFAULT_TEST_OUTPUT / "step5_final"

DEFAULT_CANDIDATES_JSON = DEFAULT_STEP3_OUT / "step3_ontology_candidates.json"
DEFAULT_PROMPT_DIR = ROOT_DIR / "data_exploration" / "curation_tools" / "llm_curation"

_DEFAULT_TOOLBAR_VALUES = {
    "base_out_input": str(DEFAULT_TEST_OUTPUT),
    "s1_sub_input": "step1_evidence",
    "s2_sub_input": "step2_normalized",
    "s3_sub_input": "step3_ontology_candidates",
    "s4_sub_input": "step4_backfilled",
    "s5_sub_input": "step5_final",
    "pub_md_input": str(FULL_TEXT_MD_DIR),
    "mavedb_meta_input": str(MAVEDB_METADATA_OUTPUT_DIR),
    "prompt_dir_input": str(DEFAULT_PROMPT_DIR),
}

_DEFAULT_ACTIVE_PATH_VALUES = {
    "s1_input_dir_val": str(Path(FULL_TEXT_MD_DIR).resolve()),
    "s1_output_dir_val": str(DEFAULT_STEP1_OUT.resolve()),
    "s1_prompt_file_val": str(
        (DEFAULT_PROMPT_DIR / "step1_evidence_extraction_prompt.md").resolve()
    ),
    "s1_log_file_val": str(
        (DEFAULT_TEST_OUTPUT / "step1_evidence_extraction.log").resolve()
    ),
    "s2_input_dir_val": str(DEFAULT_STEP1_OUT.resolve()),
    "s2_output_dir_val": str(DEFAULT_STEP2_OUT.resolve()),
    "s2_prompt_file_val": str(
        (DEFAULT_PROMPT_DIR / "step2_specific_term_extraction.md").resolve()
    ),
    "s2_log_file_val": str(
        (DEFAULT_TEST_OUTPUT / "step2_specific_term_extraction.log").resolve()
    ),
    "s2_mavedb_dir_val": str(Path(MAVEDB_METADATA_OUTPUT_DIR).resolve()),
    "s3_step1_dir_val": str(DEFAULT_STEP1_OUT.resolve()),
    "s3_step2_dir_val": str(DEFAULT_STEP2_OUT.resolve()),
    "s3_output_dir_val": str(DEFAULT_STEP3_OUT.resolve()),
    "s3_log_file_val": str(
        (DEFAULT_TEST_OUTPUT / "step3_candidate_discovery.log").resolve()
    ),
    "s3_prompt_file_val": str(
        (DEFAULT_PROMPT_DIR / "step3_candidate_discovery_prompt.md").resolve()
    ),
    "s4_step2_dir_val": str(DEFAULT_STEP2_OUT.resolve()),
    "s4_step4_dir_val": str(DEFAULT_STEP4_OUT.resolve()),
    "s5_source_dir_val": str(DEFAULT_STEP4_OUT.resolve()),
    "s5_output_dir_val": str(DEFAULT_STEP5_OUT.resolve()),
    "s5_mavedb_dir_val": str(Path(MAVEDB_METADATA_OUTPUT_DIR).resolve()),
}

_TOOLBAR_PATH_KEYS = {
    "base_out_input": {
        "s1_output_dir_val": "step1_out",
        "s2_output_dir_val": "step2_out",
        "s3_output_dir_val": "step3_out",
        "s4_step4_dir_val": "step4_out",
        "s5_source_dir_val": "step4_out",
        "s5_output_dir_val": "step5_out",
        "s1_log_file_val": "step1_log",
        "s2_log_file_val": "step2_log",
        "s3_log_file_val": "step3_log",
    },
    "s1_sub_input": {
        "s1_output_dir_val": "step1_out",
        "s2_input_dir_val": "step1_out",
        "s3_step1_dir_val": "step1_out",
    },
    "s2_sub_input": {
        "s2_output_dir_val": "step2_out",
        "s3_step2_dir_val": "step2_out",
        "s4_step2_dir_val": "step2_out",
        "s5_source_dir_val": "step2_out",
    },
    "s3_sub_input": {
        "s3_output_dir_val": "step3_out",
    },
    "s4_sub_input": {
        "s4_step4_dir_val": "step4_out",
        "s5_source_dir_val": "step4_out",
    },
    "s5_sub_input": {
        "s5_output_dir_val": "step5_out",
    },
    "pub_md_input": {
        "s1_input_dir_val": "pub_md_dir",
    },
    "mavedb_meta_input": {
        "s2_mavedb_dir_val": "mavedb_meta_dir",
        "s5_mavedb_dir_val": "mavedb_meta_dir",
    },
    "prompt_dir_input": {
        "s1_prompt_file_val": "step1_prompt",
        "s2_prompt_file_val": "step2_prompt",
        "s3_prompt_file_val": "step3_prompt",
    },
}

st.set_page_config(
    page_title="MaveDB LLM Curation Control Center",
    page_icon="🧬",
    layout="wide",
)


def load_candidates_file(file_path: Path) -> dict:
    """Load JSON file containing Step 3 ontology candidates."""
    if not file_path.exists():
        return {}
    try:
        return json.loads(file_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
        st.error(f"Error reading {file_path}: {e}")
        return {}


def save_decision_audit_trail(file_path: Path, decisions_data: dict) -> None:
    """Save candidates decision dictionary to JSON audit file."""
    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_text(json.dumps(decisions_data, indent=2), encoding="utf-8")


def init_session_state(
    candidates_data: dict, candidates_path: Path | None = None
) -> None:
    """Initialize decision state for candidates, syncing with existing schema allowed values and prior audit logs."""
    candidate_identity = None
    if candidates_path:
        resolved_candidates_path = Path(candidates_path).resolve()
        candidate_identity = (
            str(resolved_candidates_path),
            calculate_file_signature(resolved_candidates_path),
        )

    if st.session_state.get("_candidate_identity") != candidate_identity:
        st.session_state.pop("decisions", None)
        st.session_state["candidates_raw"] = candidates_data
        st.session_state["_candidate_identity"] = candidate_identity

    if "decisions" not in st.session_state:
        existing_literals = get_existing_literals(DEFAULT_SCHEMA_PATH)

        saved_audit = {}
        if candidates_path:
            audit_path = candidates_path.parent / "approved_ontology_terms.json"
            if audit_path.exists():
                try:
                    saved_audit = json.loads(audit_path.read_text(encoding="utf-8"))
                except (FileNotFoundError, json.JSONDecodeError, OSError):
                    saved_audit = {}

        decisions = {}
        for field, cand_list in candidates_data.items():
            decisions[field] = []
            allowed_vocab = set(existing_literals.get(field, []))

            saved_field_map = {}
            if field in saved_audit and isinstance(saved_audit[field], list):
                for saved_entry in saved_audit[field]:
                    orig_k = saved_entry.get("original_term") or saved_entry.get("term")
                    if orig_k:
                        saved_field_map[orig_k] = saved_entry

            for item in cand_list:
                term = item.get("proposed_new_term") or item.get("proposed_label") or ""
                saved_entry = saved_field_map.get(term)

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


def reset_default_paths() -> None:
    """Restore toolbar and active-step path widgets to their default values."""
    st.session_state.update(_DEFAULT_TOOLBAR_VALUES)
    st.session_state.update(_DEFAULT_ACTIVE_PATH_VALUES)
    st.session_state["_toolbar_path_values"] = _DEFAULT_TOOLBAR_VALUES.copy()


def sync_active_pipeline_paths(
    toolbar_values: dict[str, str], paths: dict[str, Path]
) -> None:
    """Update step path widgets when their corresponding toolbar setting changes."""
    previous_toolbar_values = st.session_state.get("_toolbar_path_values")
    if previous_toolbar_values is not None:
        changed_toolbar_keys = {
            key
            for key, value in toolbar_values.items()
            if previous_toolbar_values.get(key) != value
        }

        for toolbar_key in changed_toolbar_keys:
            for widget_key, path_key in _TOOLBAR_PATH_KEYS.get(toolbar_key, {}).items():
                if path_key in paths:
                    st.session_state[widget_key] = str(paths[path_key])

    st.session_state["_toolbar_path_values"] = toolbar_values.copy()


def main():
    st.title("🧬 MaveDB LLM Curation Control Center")
    st.caption(
        "End-to-end pipeline execution and curation: Evidence Extraction (Step 1) ➔ Term Normalization (Step 2) ➔ Candidate Discovery (Step 3a) ➔ Candidate Review (Step 3b) ➔ Step 4 Backfill ➔ Step 5 Final Metadata"
    )

    if step4_message := st.session_state.pop("_step4_backfill_message", None):
        st.success(step4_message)
    if step3a_message := st.session_state.pop("_step3a_message", None):
        st.success(step3a_message)
    if step5_message := st.session_state.pop("_step5_message", None):
        st.success(step5_message)

    # Sidebar setup
    st.sidebar.header("⚙️ Global Execution Settings")

    selected_model = st.sidebar.selectbox(
        "LLM Model ID",
        options=[
            "google/gemini-3.6-flash",
            "google/gemini-3.5-flash-lite",
            "google/gemini-3.1-pro-preview",
        ],
        index=0,
    )

    max_workers_slider = st.sidebar.slider(
        "Concurrency (Max Workers)",
        min_value=1,
        max_value=16,
        value=8,
    )

    st.sidebar.divider()
    st.sidebar.header("📁 Workspace & Output Folders")

    base_output_str = st.sidebar.text_input(
        "Base Output Directory",
        value=_DEFAULT_TOOLBAR_VALUES["base_out_input"],
        help="Root output directory containing all step output folders and log files.",
        key="base_out_input",
    )
    base_output_path = Path(base_output_str).resolve()

    st.sidebar.subheader("📂 Step Output Folders")
    s1_sub_str = st.sidebar.text_input(
        "Step 1 - Evidence Extraction",
        value=_DEFAULT_TOOLBAR_VALUES["s1_sub_input"],
        help="Output folder name inside Base Output Directory for Step 1 evidence extraction.",
        key="s1_sub_input",
    )
    s2_sub_str = st.sidebar.text_input(
        "Step 2 - Specific Term Normalization",
        value=_DEFAULT_TOOLBAR_VALUES["s2_sub_input"],
        help="Output folder name inside Base Output Directory for Step 2 normalized terms.",
        key="s2_sub_input",
    )
    s3_sub_str = st.sidebar.text_input(
        "Step 3a - Candidate Discovery",
        value=_DEFAULT_TOOLBAR_VALUES["s3_sub_input"],
        help="Output folder name inside Base Output Directory for Step 3a candidate discovery.",
        key="s3_sub_input",
    )
    s4_sub_str = st.sidebar.text_input(
        "Step 4 - Approved Terms Backfill",
        value=_DEFAULT_TOOLBAR_VALUES["s4_sub_input"],
        help="Output folder name inside Base Output Directory for Step 4 backfilled results.",
        key="s4_sub_input",
    )
    s5_sub_str = st.sidebar.text_input(
        "Step 5 - Final Metadata",
        value=_DEFAULT_TOOLBAR_VALUES["s5_sub_input"],
        help="Output folder name inside Base Output Directory for final schema-projected metadata and CSV.",
        key="s5_sub_input",
    )

    st.sidebar.divider()
    st.sidebar.subheader("📖 Source Data & Templates")

    pub_md_str = st.sidebar.text_input(
        "Publication Text Dir (.md)",
        value=_DEFAULT_TOOLBAR_VALUES["pub_md_input"],
        help="Directory containing publication markdown files.",
        key="pub_md_input",
    )
    pub_md_path = Path(pub_md_str).resolve()

    mavedb_meta_str = st.sidebar.text_input(
        "MaveDB Metadata Dir",
        value=_DEFAULT_TOOLBAR_VALUES["mavedb_meta_input"],
        help="Directory containing cached MaveDB entry JSONs.",
        key="mavedb_meta_input",
    )
    mavedb_meta_path = Path(mavedb_meta_str).resolve()

    prompt_dir_str = st.sidebar.text_input(
        "Prompt Templates Dir",
        value=_DEFAULT_TOOLBAR_VALUES["prompt_dir_input"],
        help="Directory containing prompt template Markdown files.",
        key="prompt_dir_input",
    )
    prompt_dir_path = Path(prompt_dir_str).resolve()

    st.sidebar.button(
        "↩️ Reset Default Paths",
        help="Restore all toolbar and Active Step path fields to their defaults.",
        on_click=reset_default_paths,
    )

    # Compute full output paths by joining Base Output Directory + Subdirectory
    step1_out_path = (base_output_path / s1_sub_str.strip()).resolve()
    step2_out_path = (base_output_path / s2_sub_str.strip()).resolve()
    step3_out_path = (base_output_path / s3_sub_str.strip()).resolve()
    step4_out_path = (base_output_path / s4_sub_str.strip()).resolve()
    step5_out_path = (base_output_path / s5_sub_str.strip()).resolve()
    candidates_json_path = step3_out_path / "step3_ontology_candidates.json"

    # Consolidated Pipeline Paths Map
    paths = {
        "base_out": base_output_path,
        "pub_md_dir": pub_md_path,
        "mavedb_meta_dir": mavedb_meta_path,
        "step1_out": step1_out_path,
        "step2_out": step2_out_path,
        "step3_out": step3_out_path,
        "step4_out": step4_out_path,
        "step5_out": step5_out_path,
        "candidates_json": candidates_json_path,
        "prompt_dir": prompt_dir_path,
        "step1_prompt": prompt_dir_path / "step1_evidence_extraction_prompt.md",
        "step2_prompt": prompt_dir_path / "step2_specific_term_extraction.md",
        "step3_prompt": prompt_dir_path / "step3_candidate_discovery_prompt.md",
        "step1_log": base_output_path / "step1_evidence_extraction.log",
        "step2_log": base_output_path / "step2_specific_term_extraction.log",
        "step3_log": base_output_path / "step3_candidate_discovery.log",
        "step4_log": base_output_path / "step4_backfill.log",
        "step5_log": base_output_path / "step5_final_metadata.log",
        "pipeline_manifest": base_output_path / "pipeline_manifest.json",
    }

    sync_active_pipeline_paths(
        {
            "base_out_input": base_output_str,
            "s1_sub_input": s1_sub_str,
            "s2_sub_input": s2_sub_str,
            "s3_sub_input": s3_sub_str,
            "s4_sub_input": s4_sub_str,
            "s5_sub_input": s5_sub_str,
            "pub_md_input": pub_md_str,
            "mavedb_meta_input": mavedb_meta_str,
            "prompt_dir_input": prompt_dir_str,
        },
        paths,
    )

    candidates_path = paths["candidates_json"]

    candidate_step2_dir = Path(
        st.session_state.get("s3_step2_dir_val", paths["step2_out"])
    )
    candidate_step4_dir = Path(
        st.session_state.get("s4_step4_dir_val", paths["step4_out"])
    )
    effective_candidate_input_dir = resolve_effective_normalized_dir(
        candidate_step2_dir,
        candidate_step4_dir,
        manifest_path=paths["pipeline_manifest"],
        decisions_file=candidates_path.parent / "approved_ontology_terms.json",
    )
    st.session_state.setdefault("s5_source_dir_val", str(effective_candidate_input_dir))

    if st.sidebar.button("Reload Session Cache"):
        st.session_state.clear()
        st.rerun()

    # Load candidates if file exists
    candidates_data = load_candidates_file(candidates_path)
    candidates_are_current = is_candidate_discovery_current(
        candidates_path,
        effective_candidate_input_dir,
        paths["pipeline_manifest"],
    )
    if candidates_data and candidates_are_current:
        init_session_state(candidates_data, candidates_path=candidates_path)
        st.session_state["_candidate_data_stale"] = False
    else:
        st.session_state["_candidate_data_stale"] = bool(
            candidates_data and not candidates_are_current
        )
        st.session_state.pop("candidates_raw", None)
        st.session_state.pop("decisions", None)
        st.session_state.pop("_candidate_identity", None)

    # Keep the selected tab in Streamlit session state so button-triggered reruns
    # (for example, Step 2 normalization or Step 5 assembly) return to the tab
    # where the user initiated the action instead of resetting to Step 1.
    tab_step1, tab_step2, tab_step3a, tab_step3b, tab_step4, tab_step5 = st.tabs(
        [
            "⚡ Step 1: Evidence Extraction",
            "🏷️ Step 2: Term Normalization",
            "💡 Step 3a: Candidate Discovery",
            "🔍 Step 3b: Candidate Review & Schema Diff",
            "🔄 Step 4: Backfill Approved Terms",
            "🧬 Step 5: Final Metadata",
        ],
        key="pipeline_tabs",
        on_change="rerun",
    )

    # -----------------------------------------------------------------------------
    # TAB 1: Step 1 Evidence Extraction
    # -----------------------------------------------------------------------------
    with tab_step1:
        st.header("⚡ Step 1: Evidence Extraction")
        st.caption(
            "Locates experiment targets and extracts verbatim quotes from publication text without normalization."
        )

        with st.expander("📍 Active Step 1 Pipeline Paths", expanded=False):
            col1_s1, col2_s1 = st.columns(2)
            s1_input_dir = col1_s1.text_input(
                "Publication Full Text Directory (.md)",
                value=str(paths["pub_md_dir"]),
                key="s1_input_dir_val",
            )
            s1_output_dir = col2_s1.text_input(
                "Step 1 Output Directory",
                value=str(paths["step1_out"]),
                key="s1_output_dir_val",
            )

            col3_s1, col4_s1 = st.columns(2)
            s1_prompt_file = col3_s1.text_input(
                "Step 1 Prompt Template",
                value=str(paths["step1_prompt"]),
                key="s1_prompt_file_val",
            )
            s1_log_file = col4_s1.text_input(
                "Step 1 Log File",
                value=str(paths["step1_log"]),
                key="s1_log_file_val",
            )

        s1_schema_str = st.text_input(
            "Extraction Schema",
            value="curation_tools.llm_curation.llm_curation_schema:EvidenceExtractionSchema",
        )

        col_opt1, col_opt2, col_opt3 = st.columns(3)
        s1_overwrite = col_opt1.checkbox(
            "Overwrite existing outputs", value=True, key="s1_ov"
        )
        s1_create_csv = col_opt2.checkbox("Create merged CSV", value=True, key="s1_csv")
        s1_verbose = col_opt3.checkbox("Verbose logging", value=True, key="s1_verb")

        st.divider()

        # Pre-run File / URN Status Table
        s1_status_records = get_mavedb_urn_status(
            mapping_file=MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
            output_dir=Path(s1_output_dir),
            metadata_dir=MAVEDB_METADATA_OUTPUT_DIR,
            target_urns=None,
        )
        st.subheader(
            f"📂 MaveDB Datasets Queued ({len(s1_status_records)} URN datasets matched)"
        )

        s1_urn_filter_mode = st.radio(
            "URN Filter Mode",
            options=["All MaveDB URNs", "Selected URNs (checkboxes)"],
            index=0,
            key="s1_urn_filter_mode",
            horizontal=True,
        )
        s1_selected_urns = []

        if s1_status_records:
            df_s1 = pd.DataFrame(s1_status_records)
            all_s1_urns_selected = s1_urn_filter_mode == "All MaveDB URNs"
            df_s1.insert(0, "run", all_s1_urns_selected)
            disabled_s1_columns = [
                "urn",
                "title",
                "primary_dois",
                "status",
                "output_json",
            ]
            if all_s1_urns_selected:
                disabled_s1_columns.append("run")

            edited_s1_df = st.data_editor(
                df_s1,
                column_config={
                    "run": st.column_config.CheckboxColumn(
                        "Extract",
                        help="Select this MaveDB dataset for Step 1 evidence extraction.",
                        default=False,
                    ),
                    "urn": st.column_config.TextColumn("MaveDB URN"),
                    "title": st.column_config.TextColumn("Dataset Title"),
                    "primary_dois": st.column_config.TextColumn("Primary DOIs"),
                    "status": st.column_config.TextColumn("Status"),
                    "output_json": st.column_config.TextColumn("Generated Output"),
                },
                disabled=disabled_s1_columns,
                hide_index=True,
                width="stretch",
                key=(
                    "s1_urn_selection_table_all"
                    if all_s1_urns_selected
                    else "s1_urn_selection_table_selected"
                ),
            )
            s1_selected_urns = edited_s1_df.loc[
                edited_s1_df["run"].fillna(False), "urn"
            ].tolist()

        parsed_target_urns = (
            None
            if s1_urn_filter_mode == "All MaveDB URNs"
            else parse_target_urns(s1_selected_urns)
        )

        st.divider()

        if st.button(
            "🚀 Execute Step 1 Evidence Extraction",
            type="primary",
            key="btn_run_s1",
            disabled=(
                s1_urn_filter_mode == "Selected URNs (checkboxes)"
                and not s1_selected_urns
            ),
        ):
            with st.spinner("Extracting verbatim evidence quotes via LLM..."):
                try:
                    schema_cls = load_extraction_schema(s1_schema_str)
                    urn_to_dois = load_mavedb_urn_to_dois(
                        MAVEDB_URN_TO_DOIS_OUTPUT_FILE
                    )

                    if parsed_target_urns:
                        filtered_urn_to_dois: dict[str, list[str]] = {}
                        for target_urn in sorted(parsed_target_urns):
                            if target_urn in urn_to_dois:
                                filtered_urn_to_dois[target_urn] = urn_to_dois[
                                    target_urn
                                ]
                            else:
                                urn_stem = format_urn_for_filename(target_urn)
                                meta_file = (
                                    MAVEDB_METADATA_OUTPUT_DIR / f"{urn_stem}.json"
                                )
                                if meta_file.is_file():
                                    try:
                                        entry = json.loads(
                                            meta_file.read_text(encoding="utf-8")
                                        )
                                        dois = (
                                            get_dois_from_mavedb_entry(entry, log=False)
                                            or []
                                        )
                                        filtered_urn_to_dois[target_urn] = dois
                                    except (json.JSONDecodeError, OSError):
                                        filtered_urn_to_dois[target_urn] = []
                                else:
                                    filtered_urn_to_dois[target_urn] = []
                        urn_to_dois = filtered_urn_to_dois

                    bulk_extract_evidence_for_mavedb_urns(
                        urn_to_dois=urn_to_dois,
                        extraction_schema=schema_cls,
                        output_dir=Path(s1_output_dir),
                        log_file=Path(s1_log_file),
                        prompt_template_file=Path(s1_prompt_file),
                        publication_full_text_dir=Path(s1_input_dir),
                        max_workers=max_workers_slider,
                        overwrite=s1_overwrite,
                        model_name=selected_model,
                        create_csv=s1_create_csv,
                        verbose=s1_verbose,
                    )
                    record_pipeline_step(
                        paths["pipeline_manifest"],
                        "step1",
                        input_dir=Path(s1_input_dir),
                        output_dir=Path(s1_output_dir),
                        selected_items=s1_selected_urns,
                    )
                    st.success("Step 1 Evidence Extraction Complete!")
                    st.toast(
                        "Step 1 Evidence Extraction completed successfully!", icon="✅"
                    )
                except Exception as ex:
                    st.error(f"Step 1 Extraction Failed: {ex}")

        with st.expander("📜 Live Execution Log Stream"):
            st.code(read_last_log_lines(Path(s1_log_file), num_lines=30))

        # Output JSON Inspector
        s1_out_files = (
            sorted(list(Path(s1_output_dir).glob("*.json")))
            if Path(s1_output_dir).is_dir()
            else []
        )
        if s1_out_files:
            with st.expander("🔍 Inspect Step 1 Evidence Outputs"):
                selected_s1_json = st.selectbox(
                    "Select Output JSON", options=[p.name for p in s1_out_files]
                )
                if selected_s1_json:
                    p = Path(s1_output_dir) / selected_s1_json
                    st.json(json.loads(p.read_text(encoding="utf-8")))

    # -----------------------------------------------------------------------------
    # TAB 2: Step 2 Term Normalization
    # -----------------------------------------------------------------------------
    with tab_step2:
        st.header("🏷️ Step 2: Specific Term Normalization")
        st.caption(
            "Maps verbatim Step 1 evidence quotes to controlled vocabularies without needing full publication text."
        )

        with st.expander("📍 Active Step 2 Pipeline Paths", expanded=False):
            col1_s2, col2_s2 = st.columns(2)
            s2_input_dir = col1_s2.text_input(
                "Step 1 Evidence Directory (Input)",
                value=str(paths["step1_out"]),
                key="s2_input_dir_val",
            )
            s2_output_dir = col2_s2.text_input(
                "Step 2 Output Directory",
                value=str(paths["step2_out"]),
                key="s2_output_dir_val",
            )

            col3_s2, col4_s2 = st.columns(2)
            s2_prompt_file = col3_s2.text_input(
                "Step 2 Prompt Template",
                value=str(paths["step2_prompt"]),
                key="s2_prompt_file_val",
            )
            s2_log_file = col4_s2.text_input(
                "Step 2 Log File",
                value=str(paths["step2_log"]),
                key="s2_log_file_val",
            )

            s2_mavedb_dir = st.text_input(
                "MaveDB Metadata Directory (Optional)",
                value=str(paths["mavedb_meta_dir"]),
                key="s2_mavedb_dir_val",
            )

        col_s2_o1, col_s2_o2, col_s2_o3 = st.columns(3)
        s2_overwrite = col_s2_o1.checkbox(
            "Overwrite existing outputs", value=True, key="s2_ov"
        )
        s2_create_csv = col_s2_o2.checkbox(
            "Create merged CSV", value=True, key="s2_csv"
        )
        s2_verbose = col_s2_o3.checkbox("Verbose logging", value=True, key="s2_verb")

        st.divider()

        # Pre-run File Status Table
        s2_status_records = get_step2_file_status(
            Path(s2_input_dir), Path(s2_output_dir)
        )
        st.subheader(
            f"📂 Step 1 Evidence Files queued ({len(s2_status_records)} files found)"
        )

        s2_file_filter_mode = st.radio(
            "File Filter Mode",
            options=["All Step 1 Evidence Files", "Selected Files (checkboxes)"],
            index=0,
            key="s2_file_filter_mode",
            horizontal=True,
        )
        s2_selected_files = []

        if s2_status_records:
            df_s2 = pd.DataFrame(s2_status_records)
            all_s2_files_selected = s2_file_filter_mode == "All Step 1 Evidence Files"
            df_s2.insert(0, "run", all_s2_files_selected)
            disabled_s2_columns = [
                "file_name",
                "status",
                "other_fields_count",
                "output_path",
            ]
            if s2_file_filter_mode == "All Step 1 Evidence Files":
                disabled_s2_columns.append("run")

            edited_s2_df = st.data_editor(
                df_s2,
                column_config={
                    "run": st.column_config.CheckboxColumn(
                        "Normalize",
                        help="Select this file for Step 2 term normalization.",
                        default=False,
                    ),
                    "file_name": st.column_config.TextColumn("Evidence JSON File"),
                    "status": st.column_config.TextColumn("Status"),
                    "other_fields_count": st.column_config.TextColumn(
                        "'Other' Fields Count"
                    ),
                    "output_path": st.column_config.TextColumn(
                        "Normalized Output Path"
                    ),
                },
                disabled=disabled_s2_columns,
                hide_index=True,
                width="stretch",
                key=(
                    "s2_file_selection_table_all"
                    if all_s2_files_selected
                    else "s2_file_selection_table_selected"
                ),
            )
            s2_selected_files = edited_s2_df.loc[
                edited_s2_df["run"].fillna(False), "file_name"
            ].tolist()

        s2_target_files = (
            None
            if s2_file_filter_mode == "All Step 1 Evidence Files"
            else s2_selected_files
        )

        st.divider()

        if st.button(
            "🚀 Execute Step 2 Term Normalization",
            type="primary",
            key="btn_run_s2",
            disabled=(
                s2_file_filter_mode == "Selected Files (checkboxes)"
                and not s2_selected_files
            ),
        ):
            with st.spinner(
                "Normalizing evidence quotes to controlled vocabularies..."
            ):
                try:
                    mavedb_dir_p = (
                        Path(s2_mavedb_dir)
                        if s2_mavedb_dir and Path(s2_mavedb_dir).exists()
                        else None
                    )
                    normalize_evidence_artifacts(
                        step1_dir=Path(s2_input_dir),
                        output_dir=Path(s2_output_dir),
                        log_file=Path(s2_log_file),
                        prompt_template_file=Path(s2_prompt_file),
                        mavedb_metadata_dir=mavedb_dir_p,
                        selected_files=s2_target_files,
                        max_workers=max_workers_slider,
                        overwrite=s2_overwrite,
                        model_name=selected_model,
                        create_csv=s2_create_csv,
                        verbose=s2_verbose,
                    )
                    record_pipeline_step(
                        paths["pipeline_manifest"],
                        "step2",
                        input_dir=Path(s2_input_dir),
                        output_dir=Path(s2_output_dir),
                        selected_items=s2_selected_files,
                    )
                    st.success("Step 2 Term Normalization Complete!")
                    st.toast(
                        "Step 2 Term Normalization completed successfully!", icon="✅"
                    )
                except Exception as ex:
                    st.error(f"Step 2 Normalization Failed: {ex}")

        with st.expander("📜 Live Execution Log Stream"):
            st.code(read_last_log_lines(Path(s2_log_file), num_lines=30))

        # Output JSON Inspector
        s2_out_files = (
            sorted(list(Path(s2_output_dir).glob("*.json")))
            if Path(s2_output_dir).is_dir()
            else []
        )
        if s2_out_files:
            with st.expander("🔍 Inspect Step 2 Normalized Outputs"):
                selected_s2_json = st.selectbox(
                    "Select Normalized Output JSON",
                    options=[p.name for p in s2_out_files],
                )
                if selected_s2_json:
                    p = Path(s2_output_dir) / selected_s2_json
                    data_s2 = json.loads(p.read_text(encoding="utf-8"))

                    # Highlight 'Other' fields
                    other_fields = [k for k, v in data_s2.items() if v == "Other"]
                    if other_fields:
                        st.warning(
                            f"Fields classified as 'Other': `{', '.join(other_fields)}`"
                        )
                    st.json(data_s2)

    # -----------------------------------------------------------------------------
    # TAB 3a: Step 3a Candidate Discovery
    # -----------------------------------------------------------------------------
    with tab_step3a:
        st.header("💡 Step 3a: Candidate Discovery")
        st.caption(
            "Aggregates recurring 'Other' evidence snippets across the corpus and uses LLM synthesis to propose reusable ontology candidate terms."
        )

        with st.expander("📍 Active Step 3a Pipeline Paths", expanded=False):
            col1_s3, col2_s3 = st.columns(2)
            s3_step1_dir = col1_s3.text_input(
                "Step 1 Evidence Directory (Source)",
                value=str(paths["step1_out"]),
                key="s3_step1_dir_val",
            )
            s3_step2_dir = col2_s3.text_input(
                "Step 2 Normalized Directory (Filter)",
                value=str(paths["step2_out"]),
                key="s3_step2_dir_val",
            )

            col3_s3, col4_s3 = st.columns(2)
            s3_output_dir = col3_s3.text_input(
                "Step 3a Output Directory",
                value=str(paths["step3_out"]),
                key="s3_output_dir_val",
            )
            s3_log_file = col4_s3.text_input(
                "Step 3a Log File",
                value=str(paths["step3_log"]),
                key="s3_log_file_val",
            )

            s3_prompt_file = st.text_input(
                "Step 3a Prompt Template",
                value=str(paths["step3_prompt"]),
                key="s3_prompt_file_val",
            )

        s3_verbose = st.checkbox("Verbose prompt logging", value=True, key="s3_verb")

        st.divider()

        # Pre-run Corpus "Other" Analysis Summary
        step4_output_dir = Path(
            st.session_state.get("s4_step4_dir_val", paths["step4_out"])
        )
        decisions_file_path = candidates_path.parent / "approved_ontology_terms.json"
        s3_effective_step2_dir = resolve_effective_normalized_dir(
            Path(s3_step2_dir),
            step4_output_dir,
            manifest_path=paths["pipeline_manifest"],
            decisions_file=decisions_file_path,
        )
        st.caption(f"Analysis source: `{s3_effective_step2_dir}`")

        # Pre-run file selector. Only normalized outputs can contribute to discovery.
        s3_status_records = [
            record
            for record in get_step2_file_status(
                Path(s3_step1_dir), s3_effective_step2_dir
            )
            if record["status"] == "Completed"
        ]
        st.subheader(
            f"📂 Normalized Evidence Files queued ({len(s3_status_records)} files found)"
        )

        s3_file_filter_mode = st.radio(
            "File Filter Mode",
            options=["All Normalized Evidence Files", "Selected Files (checkboxes)"],
            index=0,
            key="s3_file_filter_mode",
            horizontal=True,
        )
        s3_selected_files = []

        if s3_status_records:
            df_s3_files = pd.DataFrame(s3_status_records)
            all_s3_files_selected = (
                s3_file_filter_mode == "All Normalized Evidence Files"
            )
            df_s3_files.insert(0, "run", all_s3_files_selected)
            disabled_s3_columns = [
                "file_name",
                "status",
                "other_fields_count",
                "output_path",
            ]
            if all_s3_files_selected:
                disabled_s3_columns.append("run")

            edited_s3_files_df = st.data_editor(
                df_s3_files,
                column_config={
                    "run": st.column_config.CheckboxColumn(
                        "Discover",
                        help="Select this file for Step 3a candidate discovery.",
                        default=False,
                    ),
                    "file_name": st.column_config.TextColumn("Evidence JSON File"),
                    "status": st.column_config.TextColumn("Status"),
                    "other_fields_count": st.column_config.TextColumn(
                        "'Other' Fields Count"
                    ),
                    "output_path": st.column_config.TextColumn(
                        "Normalized Output Path"
                    ),
                },
                disabled=disabled_s3_columns,
                hide_index=True,
                width="stretch",
                key=(
                    "s3_file_selection_table_all"
                    if all_s3_files_selected
                    else "s3_file_selection_table_selected"
                ),
            )
            s3_selected_files = edited_s3_files_df.loc[
                edited_s3_files_df["run"].fillna(False), "file_name"
            ].tolist()
        else:
            st.info(
                "No normalized evidence files are available for candidate discovery."
            )

        s3_target_files = (
            None
            if s3_file_filter_mode == "All Normalized Evidence Files"
            else s3_selected_files
        )

        s3_summary_records = get_step3_other_corpus_summary(
            Path(s3_step1_dir),
            s3_effective_step2_dir,
            selected_files=s3_selected_files,
        )
        st.subheader(
            f"📊 Corpus 'Other' Evidence Summary ({len(s3_summary_records)} fields with unmapped 'Other' evidence)"
        )

        if s3_summary_records:
            df_s3 = pd.DataFrame(s3_summary_records)
            st.dataframe(
                df_s3,
                column_config={
                    "field_name": st.column_config.TextColumn("Metadata Field"),
                    "other_instances_count": st.column_config.NumberColumn(
                        "'Other' Instances Count"
                    ),
                    "sample_evidence": st.column_config.TextColumn(
                        "Sample Verbatim Evidence"
                    ),
                },
                hide_index=True,
                width="stretch",
            )
        else:
            st.info(
                "No 'Other' fields currently detected in the selected normalized outputs, or the selected directories are missing."
            )

        st.divider()

        if st.button(
            "🚀 Execute Step 3a Candidate Discovery",
            type="primary",
            key="btn_run_s3",
            disabled=(
                s3_file_filter_mode == "Selected Files (checkboxes)"
                and not s3_selected_files
            ),
        ):
            with st.spinner(
                "Synthesizing ontology candidate terms from 'Other' evidence..."
            ):
                try:
                    out_p = discover_candidates(
                        step1_dir=Path(s3_step1_dir),
                        step2_dir=s3_effective_step2_dir,
                        output_dir=Path(s3_output_dir),
                        log_file=Path(s3_log_file),
                        prompt_template_file=Path(s3_prompt_file),
                        model_name=selected_model,
                        verbose=s3_verbose,
                        selected_files=s3_target_files,
                    )
                    record_pipeline_step(
                        paths["pipeline_manifest"],
                        "step3a",
                        input_dir=s3_effective_step2_dir,
                        output_dir=Path(s3_output_dir),
                        selected_items=s3_selected_files,
                    )
                    st.session_state["_step3a_message"] = (
                        f"Step 3a Candidate Discovery complete: report saved to `{out_p}`."
                    )
                    st.toast(
                        "Step 3a Candidate Discovery completed successfully!", icon="✅"
                    )
                    st.rerun()
                except Exception as ex:
                    st.error(f"Candidate Discovery Failed: {ex}")

        with st.expander("📜 Live Execution Log Stream"):
            st.code(read_last_log_lines(Path(s3_log_file), num_lines=30))

    # -----------------------------------------------------------------------------
    # TAB 3b: Step 3b Candidate Review & Schema Diff
    # -----------------------------------------------------------------------------
    with tab_step3b:
        if st.session_state.pop("_candidate_data_stale", False):
            st.warning(
                "The Step 3a candidate report is stale relative to the current effective normalized outputs. Run Step 3a again before reviewing candidates."
            )
        decisions = st.session_state.get("decisions", {})
        if not decisions:
            st.info(
                f"No candidate decisions loaded yet. Ensure candidate discovery report exists at `{candidates_path}` and click 'Reload Session Cache'."
            )
        else:
            field_options = list(decisions.keys())
            field_labels = [
                f"{field} ({len(decisions[field])})" for field in field_options
            ]
            selected_field_idx = st.selectbox(
                "Select Field to Review",
                options=range(len(field_options)),
                format_func=lambda i: field_labels[i],
            )
            selected_field = field_options[selected_field_idx]

            st.header(f"Field: `{selected_field}`")

            existing_literals = get_existing_literals(DEFAULT_SCHEMA_PATH).get(
                selected_field, []
            )
            if existing_literals:
                with st.expander(
                    f"📋 Current Allowed Vocabulary ({len(existing_literals)} terms)",
                    expanded=False,
                ):
                    st.write(", ".join([f"`{term}`" for term in existing_literals]))

            col_b1, col_b2, col_b3, col_filter = st.columns([1, 1, 1, 3])
            if col_b1.button("✅ Approve All in Field", key="rev_app_all"):
                for item in decisions[selected_field]:
                    item["status"] = "Approved"
                st.rerun()
            if col_b2.button("❌ Reject All in Field", key="rev_rej_all"):
                for item in decisions[selected_field]:
                    item["status"] = "Rejected"
                st.rerun()
            if col_b3.button("🔄 Reset Field to Pending", key="rev_rst_all"):
                for item in decisions[selected_field]:
                    item["status"] = "Pending"
                st.rerun()

            status_filter = col_filter.radio(
                "Filter Candidate Status",
                options=["All", "Pending", "Approved", "Rejected"],
                horizontal=True,
                key="rev_status_filter",
            )

            st.divider()

            candidates_to_render = []
            for idx, item in enumerate(decisions[selected_field]):
                if status_filter == "All" or item["status"] == status_filter:
                    candidates_to_render.append((idx, item))

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
                                    source_dois = get_publication_dois_for_source_file(
                                        src,
                                        s3_effective_step2_dir,
                                        MAVEDB_URN_TO_DOIS_OUTPUT_FILE,
                                    )
                                    if source_dois:
                                        doi_links = ", ".join(
                                            f"[`{doi}`](https://doi.org/{doi})"
                                            for doi in source_dois
                                        )
                                        st.markdown(f"**Publication DOI:** {doi_links}")
                                    else:
                                        st.caption("Publication DOI: Not found")
                                    st.caption(f'> "{stmt}"')
                                else:
                                    st.caption(f'> "{ev_item}"')

            st.divider()

            # Build approved terms mapping for Schema Diff
            approved_map: dict[str, list[str]] = {}
            for field, c_list in decisions.items():
                approved_terms = [
                    c["term"] for c in c_list if c["status"] == "Approved" and c["term"]
                ]
                if approved_terms:
                    approved_map[field] = approved_terms

            st.subheader("📝 Live Unified Code Diff (`llm_curation_schema.py`)")
            if approved_map:
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

                col_apply, col_save_json = st.columns(2)
                create_backup = col_apply.checkbox(
                    "Create `.py.bak` backup before applying", value=True
                )
                if col_apply.button(
                    "🚀 Apply Approved Terms to SpecificTermExtractionSchema",
                    type="primary",
                ):
                    try:
                        apply_schema_update(
                            approved_map,
                            DEFAULT_SCHEMA_PATH,
                            create_backup=create_backup,
                        )
                        output_path = (
                            candidates_path.parent / "approved_ontology_terms.json"
                        )
                        save_decision_audit_trail(output_path, decisions)
                        st.success(
                            "Successfully updated `SpecificTermExtractionSchema` in `llm_curation_schema.py` and saved audit log!"
                        )
                        st.toast("Schema updated and audit log saved!", icon="✅")
                    except Exception as ex:
                        st.error(f"Failed to update schema: {ex}")

                if col_save_json.button("💾 Export Decision Audit Trail (JSON)"):
                    output_path = (
                        candidates_path.parent / "approved_ontology_terms.json"
                    )
                    save_decision_audit_trail(output_path, decisions)
                    st.success(f"Saved decisions to `{output_path}`")
                    st.toast(
                        f"Decision audit trail saved to `{output_path.name}`", icon="💾"
                    )

    # -----------------------------------------------------------------------------
    # TAB 4: Step 4 Backfill Approved Terms
    # -----------------------------------------------------------------------------
    with tab_step4:
        st.header("🔄 Step 4: Backfill Approved Terms")
        st.caption(
            "Visual preview of 'Other' replacements that will be applied to Step 2 copies when executing Step 4."
        )

        with st.expander("📍 Active Step 4 Pipeline Paths", expanded=False):
            col_b_in, col_b_out = st.columns(2)
            step2_dir_input = col_b_in.text_input(
                "Step 2 Source Directory",
                value=str(paths["step2_out"]),
                help="Directory containing Step 2 normalized JSON files.",
                key="s4_step2_dir_val",
            )
            step4_dir_input = col_b_out.text_input(
                "Step 4 Output Directory",
                value=str(paths["step4_out"]),
                help="Directory where backfilled Step 4 copies will be written.",
                key="s4_step4_dir_val",
            )

        decisions = st.session_state.get("decisions", {})
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
                width="stretch",
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
            save_decision_audit_trail(decisions_file_path, decisions)

            log_file = paths["step4_log"]

            try:
                copied_count, fields_updated, audit_records = backfill_approved_terms(
                    step2_dir=Path(step2_dir_input),
                    output_dir=Path(step4_dir_input),
                    decisions_file=decisions_file_path,
                    log_file=log_file,
                    create_csv=True,
                )
                record_pipeline_step(
                    paths["pipeline_manifest"],
                    "step4",
                    input_dir=Path(step2_dir_input),
                    output_dir=Path(step4_dir_input),
                    decisions_file=decisions_file_path,
                    source_step="step2",
                )
                st.success(
                    f"Step 4 Backfill Complete! Copied {copied_count} files and replaced {fields_updated} 'Other' values in `{step4_dir_input}`."
                )
                st.toast("Step 4 Backfill completed successfully!", icon="✅")
                st.session_state["_step4_backfill_message"] = (
                    f"Step 4 Backfill complete: copied {copied_count} files and replaced "
                    f"{fields_updated} 'Other' values."
                )
                st.rerun()
            except Exception as e:
                st.error(f"Backfill failed: {e}")

    # -----------------------------------------------------------------------------
    # TAB 5: Step 5 Final Metadata and CSV Review
    # -----------------------------------------------------------------------------
    with tab_step5:
        st.header("🧬 Step 5: Final Metadata")
        st.caption(
            "Combines the curated LLM metadata with MaveDB supplementary metadata, projects each study onto ObsSchema, and creates the final editable CSV."
        )

        with st.expander("📍 Active Step 5 Pipeline Paths", expanded=False):
            col_f_in, col_f_out = st.columns(2)
            s5_source_dir = col_f_in.text_input(
                "Step 4/2 Source Directory",
                value=str(effective_candidate_input_dir),
                help="Directory containing the latest normalized or Step 4 backfilled JSON files.",
                key="s5_source_dir_val",
            )
            s5_output_dir = col_f_out.text_input(
                "Step 5 Output Directory",
                value=str(paths["step5_out"]),
                help="Directory where final JSON files and step5_final_metadata.csv are written.",
                key="s5_output_dir_val",
            )
            s5_mavedb_dir = st.text_input(
                "MaveDB Metadata Directory",
                value=str(paths["mavedb_meta_dir"]),
                help="Cached MaveDB entry JSONs used to fill fields that are not LLM-curated.",
                key="s5_mavedb_dir_val",
            )

        s5_overwrite = st.checkbox(
            "Overwrite existing Step 5 JSON files",
            value=False,
            key="s5_overwrite",
            help="Leave disabled to preserve existing final outputs and manual CSV edits; enable only when you intentionally want to regenerate them.",
        )
        if st.button("🚀 Assemble Step 5 Final Metadata", type="primary"):
            try:
                output_paths = finalize_metadata(
                    normalized_metadata_dir=Path(s5_source_dir),
                    output_dir=Path(s5_output_dir),
                    log_file=paths["step5_log"],
                    mavedb_metadata_dir=Path(s5_mavedb_dir),
                    overwrite=s5_overwrite,
                    create_csv=True,
                )
                record_pipeline_step(
                    paths["pipeline_manifest"],
                    "step5",
                    input_dir=Path(s5_source_dir),
                    output_dir=Path(s5_output_dir),
                    source_step=(
                        "step4"
                        if Path(s5_source_dir).resolve()
                        == Path(paths["step4_out"]).resolve()
                        else "step2"
                    ),
                )
                st.session_state["_step5_message"] = (
                    f"Step 5 complete: assembled {len(output_paths)} study records. "
                    f"Edit the CSV below when needed."
                )
                st.rerun()
            except Exception as exc:
                st.error(f"Step 5 finalization failed: {exc}")

        csv_path = get_final_csv_path(Path(s5_output_dir))
        st.divider()
        st.subheader("✏️ Edit Final CSV Cells")
        st.caption(
            "Click any editable metadata cell, change its value, and save. Dataset identifiers and provenance columns are protected. Changes are written to the CSV and recorded in step5_csv_edit_audit.json; the per-study JSON files are left unchanged."
        )
        if not csv_path.is_file():
            st.info(
                f"No Step 5 CSV found at `{csv_path}`. Run Step 5 finalization first."
            )
        else:
            try:
                csv_df = load_final_csv(csv_path)
                csv_signature = calculate_file_signature(csv_path) or "missing"
                editor_key = f"s5_csv_editor_{csv_signature}"
                protected_columns = [
                    column
                    for column in (
                        "dataset_id",
                        "source_json_file",
                        "__source_urns",
                        "__source_files",
                    )
                    if column in csv_df.columns
                ]
                edited_csv_df = st.data_editor(
                    csv_df,
                    key=editor_key,
                    hide_index=True,
                    num_rows="fixed",
                    disabled=protected_columns,
                    width="stretch",
                )
                if st.button("💾 Save CSV Cell Edits", type="secondary"):
                    changes = save_final_csv_edits(
                        csv_path,
                        csv_df,
                        edited_csv_df,
                    )
                    st.success(f"Saved {len(changes)} changed cell(s) to `{csv_path}`.")
                    st.toast("CSV edits saved and audited.", icon="💾")
                    st.rerun()
            except (OSError, ValueError, pd.errors.ParserError) as exc:
                st.error(f"Could not load or edit the final CSV: {exc}")


if __name__ == "__main__":
    main()
