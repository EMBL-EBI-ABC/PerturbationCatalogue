# MaveDB LLM Metadata Extraction

This folder contains an analysis pipeline for deriving structured MaveDB experiment metadata from publications linked to MaveDB entries. The workflow has three main stages:

1. Fetch MaveDB entry JSON and resolve linked publication DOIs.
2. Download publication full text and convert it to Markdown.
3. Use an LLM plus a user-provided Pydantic schema to extract normalized metadata into JSON and CSV outputs.

For now, this workflow has been designed primarily for MaveDB dataset curation. Work to further develop the workflow for CRISPR and Perturb-seq curation is ongoing.

For an example MaveDB curation workflow of an individual entry, see `data_exploration/MaveDB/mavedb_example.ipynb`.

## Pipeline Overview

### 1. Collect MaveDB entries and publication text

`mavedb/processing.py`:

- reads unique MaveDB URNs from the dump at `data_exploration/MaveDB/Dump/.../csv`
- fetches the corresponding MaveDB API records
- extracts publication DOIs from each entry
- downloads publication full text with `paperscraper`
- converts PDF/XML full text into Markdown

By default it writes:

- MaveDB entry JSON to `mavedb_metadata/`
- URN-to-DOI mapping to `mavedb_urn_to_dois.json`
- downloaded raw full text to `pub_full_text_raw/`
- converted Markdown to `pub_full_text_md/`

This script currently has no CLI arguments. Running it launches the full data collection pipeline using the hard-coded/default paths in the script.

Note: publication retrieval in this step uses `paperscraper`, which reads publisher API credentials from the repository root `.env` file. Add these variables to `.env` before running the collection pipeline:

```bash
WILEY_TDM_API_TOKEN=<WILEY_TOKEN>
ELSEVIER_TDM_API_KEY=<ELSEVIER_TOKEN>
```

Run it from the repository root:

```bash
PYTHONPATH=data_exploration python -m curation_tools.llm_curation.mavedb.processing
```

### 2. Extract metadata with an LLM

`mavedb/metadata_extraction_runner.py` reads every Markdown file in `pub_full_text_md/`, augments the prompt with matching MaveDB entry metadata when available, and requests structured output using the Pydantic schema provided through `--extraction-schema`.

Outputs are written to:

- `extracted_metadata/clean/`: JSON with only the normalized fields
- `extracted_metadata/with_evidence/`: JSON including evidence quotes for each normalized field (mainly for debugging and sanity checks)
- `extracted_metadata/clean_metadata.csv`: optional flattened CSV built from the clean JSON files

Run:

```bash
PYTHONPATH=data_exploration python -m curation_tools.llm_curation.mavedb.metadata_extraction_runner \
  --extraction-schema path/to/schema.py:MyMetadataExtractionSchema
```

Useful flags:

```bash
PYTHONPATH=data_exploration python -m curation_tools.llm_curation.mavedb.metadata_extraction_runner \
  --extraction-schema path/to/schema.py:MyMetadataExtractionSchema \
  --llm-model google/gemini-2.5-flash \
  --max-workers 8 \
  --overwrite \
  --create-csv
```

Arguments:

- `--extraction-schema`: required Pydantic `BaseModel` subclass to use for extraction, formatted as `package.module:SchemaClass` or `path/to/schema.py:SchemaClass`
- `--llm-model`: model ID used for extraction; defaults to `LLM_MODEL_NAME` from the environment, or `google/gemini-3.5-flash` if unset
- `--max-workers`: number of worker threads used for bulk extraction; must be at least `1`
- `--overwrite`: overwrite existing outputs in `extracted_metadata/` instead of skipping files that already have clean JSON outputs
- `--create-csv`: create `extracted_metadata/clean_metadata.csv` from the curated JSON files in `extracted_metadata/clean/`

## Inputs and Outputs

### Inputs

- MaveDB dump CSVs in `data_exploration/MaveDB/Dump/.../csv`
- publication full text resolved via DOI (`.xml` or `.pdf` files)
- prompt template in `mavedb_metadata_extraction_prompt_template.md`
- user-provided Pydantic extraction schema

### Outputs

- cached MaveDB records in `mavedb_metadata/*.json`
- DOI lookup table in `mavedb_urn_to_dois.json`
- raw publication files in `pub_full_text_raw/`
- Markdown full text in `pub_full_text_md/*.md`
- extraction JSON in `extracted_metadata/with_evidence/*.json`
- curated extraction JSON in `extracted_metadata/clean/*.json`
- combined CSV in `extracted_metadata/clean_metadata.csv`
- progress logs in `pub_full_text_download.log` and `mavedb_metadata_extraction.log`

If multiple distinct MaveDB contexts map to the same publication, `mavedb/metadata_extraction_runner.py` writes separate output files with URN suffixes.

## Environment and Dependencies

This folder uses several external Python packages, including:

- `instructor`
- `pydantic`
- `tqdm`
- `requests`
- `paperscraper`
- `pymupdf4llm`
- `lxml` (optional but used when available for XML parsing)

The extraction script defaults to:

- model: `google/gemini-3.5-flash` (alternatively, e.g. `google/gemini-2.5-pro`)
- environment variable override: `LLM_MODEL_NAME`

The LLM client is initialized through `instructor.from_provider(..., vertexai=True, location='global')`, so you need working Vertex AI / Google credentials in the environment before running extraction.



Notes

- The Markdown conversion step tries to remove trailing reference sections conservatively.
- Publication downloads can yield either PDF or XML; both are supported.
- Existing cached files are reused unless overwrite behavior is explicitly enabled.
- One publication file, `10_1101_2024_04_26_591310.md`, is excluded from bulk LLM extraction by default in `mavedb/metadata_extraction_runner.py`, since it will yield more than 500 redundant calls to the LLM.
- Evidence quotes are intended to come from publication text, even when supplementary MaveDB metadata is injected into the prompt for disambiguation.

## Suggested Workflow

From the repository root:

```bash
PYTHONPATH=data_exploration python -m curation_tools.llm_curation.mavedb.processing
PYTHONPATH=data_exploration python -m curation_tools.llm_curation.mavedb.metadata_extraction_runner \
  --extraction-schema path/to/schema.py:MyMetadataExtractionSchema \
  --llm-model google/gemini-3.5-flash \
  --create-csv
```

Then inspect:

- `data_exploration/MaveDB/llm_metadata_extraction/extracted_metadata/clean_metadata.csv`
