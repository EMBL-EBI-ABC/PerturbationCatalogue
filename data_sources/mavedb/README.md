# MaveDB

## Background

MaveDB data can be downloaded from Zenodo and is licensed under CC0.

The ZIP archive contains:

- `main.json` — metadata
- `csv/` — data packaged per study
  - `*.scores.csv` — always present
  - `*.counts.csv` — sometimes present

## Inspect the dump schema

Run these commands from `data_sources/mavedb` with the project virtual
environment active. Install the `genson` dependency first if needed:

```bash
pip install -r requirements.txt
```

Infer a JSON Schema from the dump and save it to `tmp/inferred_schema.json`:

```bash
python utils/infer_schema.py --input tmp/mavedb_dump/main.json > tmp/inferred_schema.json
```

Write a readable key and value-type summary to `tmp/json_schema_summary.txt`:

```bash
python utils/describe_json.py --input tmp/mavedb_dump/main.json > tmp/json_schema_summary.txt
```

The inferred schema captures the structure observed in the dump. The summary
shows the nested keys and types, inspecting the first item of each list.

## Consolidated metadata and scores

`process_mavedb.py` is the single processing entry point. It incorporates the
metadata and variant extraction and the corrected MAVE/MPRA/splicing analysis
from `ot_perturbation_catalog_mave_experiments/bin/process_mavedb_data.R`.

Run from the repository root with the project environment:

```bash
.venv/bin/python data_sources/mavedb/scripts/process_mavedb.py \
  --input data_sources/mavedb/tmp/mavedb_dump \
  --output data_sources/mavedb/output/processed
```

With no arguments these are also the default paths, independently of the working
directory. `--input` accepts the dump directory, its `csv/` directory, the parent
containing `mavedb_dump/`, or `main.json` directly. `--scores` can override the CSV
directory.

The default output contains only these two tables, joined by **`dataset_id`**:

| File | Grain | Contents |
| --- | --- | --- |
| `scores.tsv` | One original variant record | Dataset ID, original variant number, source row number, nucleotide/protein/splice HGVS, simple-protein-variant parsing and score |
| `metadata.tsv` | One dataset | Source and unpacked target metadata, analysis types, variant counts and derived statistics including coverage/completeness |

### Dataset selection and complete scores

The same selection is applied to both tables. The two localization score sets
`urn:mavedb:00001237-a-1` and `urn:mavedb:00001237-a-2` are excluded. Other datasets
must have target metadata and at least one scored row eligible for the existing
analysis: protein-coding with protein HGVS (MAVE), regulatory (MPRA), or
protein-coding with splice HGVS / `other_noncoding` (splicing).

After a dataset qualifies, **all its original variant records are retained**,
including missing scores and rows outside an individual analysis subset. Rows
are never duplicated because a dataset has multiple targets or analysis types.
`source_row` is the one-based input record number and, with `dataset_id`, uniquely
identifies a row. `variant_num` retains the original accession suffix when
available. Malformed accessions are retained in `source_accession`.

`dataset_id` comes from the score filename. If a file's raw accession prefix
differs (including temporary accessions), it is recorded once as
`source_accession_prefix` in metadata. In the exceptional case of several
prefixes within one file, differing full accessions remain on their score rows.

`is_simple_p`, `ref_aa_abbr`, `position` and `alt_aa_abbr` follow the former
score extractor's simple protein-HGVS parsing. These are **source coordinates**
and three-letter amino-acid values; `=` and `*` are preserved. Scores are numeric;
missing values are `NA`. Invalid numeric values are retained in `score_raw` and
excluded from analysis. Custom auxiliary score/count columns are not exported.

`ref_aa_sym` and `alt_aa_sym` add single-letter amino-acid codes using the
reference R project's `AA_3_TO_1` mapping. Synonymous `=` uses the reference
amino-acid letter; `Ter` and `*` use `stop`. Both columns are `NA` for protein calls outside
the simple-HGVS parser, including multi-substitution calls.

### Non-redundant source metadata

Source fields are limited to those extracted by
`ot_perturbation_catalog_mave_experiments/bin/process_mavedb_data.R` or
`protvar-mave-import/scripts/extract_meta.py`: score-set title, description,
abstract, reported variant count, and target names, categories, sequences,
sequence labels/types, taxonomy and external identifiers/offsets. Experiment and
experiment-set IDs remain as parent relationships. Administrative fields,
methods, publications, dump provenance and arbitrary extra metadata are omitted.

Columns have an explicit order: dataset and parent IDs, title, description,
abstract, reported variant count, unpacked target fields, then all derived
columns. Columns are not sorted alphabetically.

Field names use snake case. Target fields are individual columns:
`target_name`, `target_category`, `sequence_type`, `sequence`, `label`, `taxid`,
`organism`, `identifier_uniprot`/`offset_uniprot`, and corresponding Ensembl and
RefSeq columns. `targets` is not exported. Single-target datasets have scalar
values; multi-target datasets have aligned JSON arrays, with `null` placeholders
for missing values. Multiple identifiers from the same database also retain
their order in arrays, including matching offsets.
NCBI taxonomy is normalized to `taxid` and enriched with `kingdom`.
Mapped UniProt metadata fills a missing external UniProt identifier, with its
provenance recorded in `identifier_uniprot_sources` and no offset invented. An
identical mapped/external ID is stored once with both sources. A differing mapped
ID remains in `uniprot_id_from_mapped_metadata`. In-frame protein-coding DNA gets a
`derived_protein_sequence`; original sequences and sequence types remain intact.
Other list/object cells are also JSON, with missing values represented by `null`
inside JSON and `NA` in TSV cells.

Each metadata record occupies exactly one physical TSV line. Tabs, carriage
returns and line breaks in scalar text are written as `\t`, `\r` and `\n`;
literal backslashes are escaped as `\\`. JSON cells use standard JSON escaping
and can still be read directly with `json.loads`. TSV quotation marks are handled
by a standard tab-delimited CSV reader. Metadata cells are checked against the
32,767-character spreadsheet limit before publishing; oversized values cause an
error rather than silent truncation or a broken spreadsheet import.

### Derived metadata

`n_variants_reported` is the original MaveDB count. `n_variants` counts exported
records; `n_scored_variants`, `n_missing_scores` and `n_invalid_scores` distinguish
their score status. `n_analyzed_variants` counts the primary analysis subset.

`analysis_types` lists all applicable analyses. `primary_analysis` prefers MAVE,
then MPRA, then splicing. For a single primary target, its corrected R analysis
metrics are flat columns: score summaries, mutation summaries, positional
coverage and completeness. Statistics use valid scores from that analysis
subset, with sample standard deviations and no score normalization.

Secondary analyses and multi-target positional metrics are omitted from
`metadata.tsv`; no `target_<analysis>_<metric>` columns are exported.
Multi-target datasets have no pooled scalar coverage: their variant
counts and general score summary count each original variant once.
`position_reference` states whether primary positional metrics use a known
UniProt offset or source coordinates. All derived target fields follow the
source columns, alongside the existing derived statistics.

### Optional development compatibility tables

Add `--compatibility-tables` to also generate:

```text
maves/singles.tsv
maves/doubles.tsv
maves/wide.tsv
```

These retain the R layout, including its legacy `score_urn` column. The canonical
tables use only `dataset_id`. Singles are averaged by mutation after applying
known UniProt offsets. Doubles include calls with at most two mutations from
score sets containing multiple mutants, including single mutants and WT
controls. Wide columns list the 20 amino acids, then `del` and `stop`, within
each score/count/standard-deviation block. Synonymous scores are selected by
reference amino-acid column name and do not depend on column order.

The new tables are staged before replacing existing outputs. A successful run
removes the processor's obsolete category-specific `full.tsv`, `meta.tsv` and
`meta_columns.tsv` files, plus compatibility files when the flag is absent.
Unrelated files in the output directory are preserved. The old `extract_*.py`
commands are deprecated aliases to this processor: they now write both tables
and accept an output **directory**, not an individual CSV filename.

The script caches NCBI's `categories.dmp` as `tmp/taxcat.dmp`, downloading it from
`https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz` when absent. Use
`--taxonomy /path/to/taxcat.dmp` for an existing pinned/offline copy. No FASTA,
sequence-mapping or count files are produced.

For the ingestion workflow below, set the following environment secrets:

- LAKE_BUCKET
- WAREHOUSE_BUCKET
- ELASTIC_ENDPOINT

## Ingest metadata

```bash
# Mirror.
wget -q -O mavedb-data.zip "https://zenodo.org/records/14172004/files/mavedb-dump.20241114101443.zip?download=1"
unzip -q -d mavedb mavedb-data.zip
gsutil -q -m rm -r "gs://${LAKE_BUCKET}/mavedb"
gsutil -q -m cp -r mavedb "gs://${LAKE_BUCKET}/"

# Process.
gsutil -q cp "gs://${LAKE_BUCKET}/mavedb/main.json" /tmp/main.json
python3 process.py \
  --input-filename /tmp/main.json \
  --output-filename /tmp/metadata.jsonl
gsutil -q -m rm -r "gs://${WAREHOUSE_BUCKET}/mavedb"
gsutil -q cp /tmp/metadata.jsonl "gs://${WAREHOUSE_BUCKET}/mavedb/metadata.jsonl"

# Ingest into Elastic.
python3 ../elastic_load.py \
  --elastic-endpoint "${ELASTIC_ENDPOINT}" \
  --elastic-index "mavedb" \
  --jsonl-data "gs://${WAREHOUSE_BUCKET}/mavedb/metadata.jsonl" \
  --id-field "urn" \
  --field-properties '{
    "geneCategory": {
      "type": "keyword"
    },
    "sequenceType": {
      "type": "keyword"
    },
    "normalisedGeneName": {
      "type": "keyword"
    }
  }'
```
