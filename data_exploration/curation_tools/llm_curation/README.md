# MaveDB LLM Metadata Extraction & Curation Control Center

This folder contains an analysis pipeline for deriving structured MaveDB experiment metadata from publications linked to MaveDB entries. The workflow consists of four main stages and a consolidated Streamlit Curation Control Center:

1. **Step 1: Evidence Extraction** — Extracts verbatim quotes from publication text using LLM queries.
2. **Step 2: Specific Term Normalization** — Maps verbatim evidence to controlled vocabularies (`SpecificTermExtractionSchema`).
3. **Step 3a: Candidate Discovery** — Aggregates recurring `"Other"` evidence and synthesizes proposed ontology candidate terms.
4. **Step 3b: Candidate Review & Schema Update** — Interactive review interface to approve/edit terms and inject them into `llm_curation_schema.py`.
5. **Step 4: Approved Terms Backfill** — Copies normalized Step 2 JSON files to Step 4 and replaces `"Other"` values with approved terms.

---

## Interactive Streamlit Curation Dashboard

To launch the consolidated 5-tab Curation Control Center:

```bash
PYTHONPATH=data_exploration ./.venv/bin/streamlit run data_exploration/curation_tools/llm_curation/review_gui.py
```

### Dashboard Tabs
- **⚡ Step 1: Evidence Extraction:** Configure models/concurrency, inspect target `.md` file status, trigger bulk evidence extraction, and view live log streams.
- **🏷️ Step 2: Term Normalization:** Inspect Step 1 evidence files and `"Other"` field counts, configure parameters, trigger normalization, and inspect mapped outputs.
- **💡 Step 3a: Candidate Discovery:** Preview corpus-level `"Other"` evidence frequencies per field and run LLM candidate discovery.
- **🔍 Step 3b: Candidate Review & Schema Diff:** Review candidate terms, edit labels, approve/reject candidates, preview live AST code diffs of `llm_curation_schema.py`, and apply approved terms to the schema.
- **🔄 Step 4: Backfill Approved Terms:** View a live preview table of all planned `"Other"` replacements before writing, then execute backfill to generate Step 4 copies, audit JSON, and compiled CSV outputs.

---

## Command Line Execution

### Step 1: Evidence Extraction
```bash
PYTHONPATH=data_exploration ./.venv/bin/python -m curation_tools.llm_curation.mavedb.mavedb_metadata_extraction_runner \
  --extraction-schema curation_tools.llm_curation.llm_curation_schema:EvidenceExtractionSchema \
  --publication-full-text-dir test_md_dir \
  --output-dir test_output/step1_evidence \
  --log-file test_output/step1_evidence_extraction.log \
  --prompt-template-file data_exploration/curation_tools/llm_curation/step1_evidence_extraction_prompt.md \
  --llm-model google/gemini-3.6-flash \
  --max-workers 8 \
  --overwrite \
  --create-csv \
  --verbose
```

### Step 2: Specific Term Extraction
```bash
PYTHONPATH=data_exploration ./.venv/bin/python -m curation_tools.llm_curation.specific_term_extraction \
  --step1-dir test_output/step1_evidence \
  --output-dir test_output/step2_normalized \
  --log-file test_output/step2_specific_term_extraction.log \
  --prompt-template-file data_exploration/curation_tools/llm_curation/step2_specific_term_extraction.md \
  --mavedb-metadata-dir data_exploration/MaveDB/llm_metadata_extraction/mavedb_metadata \
  --llm-model google/gemini-3.6-flash \
  --max-workers 8 \
  --overwrite \
  --no-csv \
  --verbose
```

### Step 3a: Candidate Discovery
```bash
PYTHONPATH=data_exploration ./.venv/bin/python -m curation_tools.llm_curation.candidate_discovery \
  --step1-dir test_output/step1_evidence \
  --step2-dir test_output/step2_normalized \
  --output-dir test_output/step3_ontology_candidates \
  --log-file test_output/step3_candidate_discovery.log \
  --prompt-template-file data_exploration/curation_tools/llm_curation/step3_candidate_discovery_prompt.md \
  --llm-model google/gemini-3.6-flash \
  --verbose
```

### Step 4: Approved Terms Backfill
```bash
PYTHONPATH=data_exploration ./.venv/bin/python -m curation_tools.llm_curation.backfill_terms \
  --step2-dir test_output/step2_normalized \
  --output-dir test_output/step4_backfilled \
  --decisions-file test_output/step3_ontology_candidates/approved_ontology_terms.json \
  --log-file test_output/step4_backfill.log \
  --create-csv
```
