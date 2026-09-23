# Native SQLite Curation Runs

Each curation batch is a new self-contained run:

```text
curation_runs/<run-name>/curation_run.sqlite3
```

The database snapshots publication Markdown, MaveDB records, prompts, model
settings, schema fingerprints, LLM artifacts, review decisions, schema-update
journals, final edits, events, and export provenance. Step 1--5 never read or
write an intermediate JSON, CSV, manifest, audit, or log file.

Source cache locations are used only when creating a run. Once creation
finishes, changing or deleting those source files cannot change the run.

## CLI

Run from the repository root with `PYTHONPATH=data_exploration`:

```bash
python -m curation_tools.llm_curation.run_cli create \
  --run-name 20260918-mavedb \
  --publication-dir data_exploration/curation_tools/llm_curation/mavedb/pub_full_text_md \
  --mavedb-metadata-dir data_exploration/MaveDB/llm_metadata_extraction/mavedb_metadata \
  --urn-to-dois-file data_exploration/MaveDB/llm_metadata_extraction/mavedb_urn_to_dois.json

python -m curation_tools.llm_curation.run_cli step --run 20260918-mavedb --step step1
python -m curation_tools.llm_curation.run_cli step --run 20260918-mavedb --step step2
python -m curation_tools.llm_curation.run_cli step --run 20260918-mavedb --step step3
# Review candidates in the dashboard, then run Step 4 and Step 5.
python -m curation_tools.llm_curation.run_cli step --run 20260918-mavedb --step step4
python -m curation_tools.llm_curation.run_cli step --run 20260918-mavedb --step step5
python -m curation_tools.llm_curation.run_cli seal --run 20260918-mavedb
python -m curation_tools.llm_curation.run_cli export \
  --run 20260918-mavedb --output-dir deliverables/20260918-mavedb
```

`status` reports execution and event history. `step --retry-failed` retries only
the latest failed execution's incomplete items and never overwrites an existing
artifact. A sealed run is immutable; it permits status inspection and explicit
final export only.

`export` is the only operation that produces files: one JSON object per study
and `curation_final_metadata.csv`. It refuses a non-empty destination unless
`--overwrite` is supplied.

## Dashboard

```bash
PYTHONPATH=data_exploration .venv/bin/streamlit run \
  data_exploration/curation_tools/llm_curation/review_gui.py
```

The dashboard creates named runs, shows their immutable configuration, invokes
stored steps by stable item ID, reviews candidates from SQLite, previews the
same resolver used for backfill, records append-only final edits, seals runs,
and exports deliverables.

## Legacy runs

Folders containing `review_state.sqlite3`, step directories, or legacy markers
are intentionally unsupported. They are neither migrated nor changed. Create a
fresh native run for every new batch.
