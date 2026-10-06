# Study metadata curation

Start here for publication-based metadata extraction, human review, and MaveDB
dataset postprocessing. The workflow has three homes:

| Home | Contents |
| --- | --- |
| `data_exploration/curation_tools/study_curation/` | Reusable code, prompts, extraction schema, dashboard, and maintained mapping rules |
| `data_exploration/curation_cache/` | Downloaded publications, converted Markdown, MaveDB records, and source lookup maps |
| `curation_runs/<run-name>/` | Run database, exported metadata, and named postprocessing outputs |

Shared defaults are defined in [paths.py](paths.py) and resolve from the
repository location. Explicit input, output, and run-root overrides remain
available. The MaveDB dump, manual workbook, dataset curation notebooks, shared
AnnData schema, and ontologies retain their existing locations.

## Code layout

```text
study_curation/
├── paths.py
├── logging_utils.py
├── sources/
│   ├── mavedb.py
│   ├── publication_text.py
│   └── xml_parser.py
├── llm/
│   ├── review_gui.py
│   ├── run_cli.py
│   ├── workflow.py
│   ├── curation_run_store.py
│   ├── llm_curation_schema.py
│   ├── prompts/
│   └── assets/
└── mavedb/
    ├── cli.py
    ├── workflow.py
    ├── standardization.py
    └── resources/
```

`sources` retrieves publications and prepares MaveDB context. `llm` extracts
evidence, normalizes metadata, supports human review, and seals runs. `mavedb`
merges exported metadata with manual curation, processes score data, standardizes
terms, and validates outputs. These stages can be run independently.

## Run the workflow

Run these commands from the repository root using the project environment:

```bash
source .venv/bin/activate
export PYTHONPATH=data_exploration
```

1. **Prepare sources**, if the caches need to be built or refreshed:

   ```bash
   python -m curation_tools.study_curation.sources.mavedb
   ```

2. **Extract and review metadata** using the dashboard:

   ```bash
   streamlit run data_exploration/curation_tools/study_curation/llm/review_gui.py
   ```

   Create a run, execute the extraction steps, review candidates and final
   fields, and seal the run. The [LLM workflow guide](llm/README.md) documents
   each step and the CLI alternative. Vertex AI credentials are needed for
   model execution.

3. **Export a sealed run** to its own export directory:

   ```bash
   python -m curation_tools.study_curation.llm.run_cli export \
     --run my-run \
     --output-dir curation_runs/my-run/exports
   ```

4. **Postprocess MaveDB datasets**, choosing a fresh processing name:

   ```bash
   python -m curation_tools.study_curation.mavedb.cli run \
     --llm-metadata curation_runs/my-run/exports/curation_final_metadata.csv \
     --manual-metadata data_exploration/MaveDB/mavedb_studies.xlsx \
     --output-dir curation_runs/my-run/postprocessing/default
   ```

   See the [MaveDB postprocessing guide](mavedb/README.md) for mappings, QC,
   optional joint artifacts, and the separate publishing command. The
   [MaveDB workflow notebook](../../MaveDB/curation_notebooks/mavedb_curation_workflow.ipynb)
   calls the same functions and uses the same defaults.

## Run history and source caches

A native run uses this layout:

```text
curation_runs/<run-name>/
├── curation_run.sqlite3
├── exports/
│   ├── curation_final_metadata.csv
│   └── <dataset-id>.json
└── postprocessing/
    └── <processing-name>/
        ├── sources/
        ├── logs/
        ├── curated/
        ├── combined/
        ├── normalized/
        ├── qc/
        ├── final/
        └── run_manifest.json
```

The database snapshots source text, records, prompts, model settings, extraction
schema, intermediate results, and review history. Moving or refreshing live
caches does not change those snapshots. Export and postprocessing directories
are created explicitly; postprocessing requires a new or empty destination.

Existing file-based runs retain their historical layouts. Native runs created
before this reorganization continue to use their schema snapshots; shared-schema
preview and update resolve the former default schema path to its new location.
Historical extraction material remains in
`data_exploration/MaveDB/llm_metadata_extraction/`, with a legacy notice.

Downloaded cache contents and run artifacts are ignored by Git. The two source
lookup JSON files and postprocessing mapping rules remain versioned. The
[cache guide](../../curation_cache/README.md) lists the source locations.
