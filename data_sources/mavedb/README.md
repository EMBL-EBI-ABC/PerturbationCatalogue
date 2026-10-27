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

## Extract target gene identifiers

After downloading and unpacking the dump, extract one row per target gene from
`main.json`. The CSV includes the flattened fields from the
`sourceScoreSetMetadata` contract, along with the existing UniProt, RefSeq and
Ensembl identifier columns and their offsets. List and object values are
JSON-encoded in their CSV cells. The extractor keeps non-human target rows and
records their source taxonomy ID.

```bash
python data_sources/mavedb/scripts/download_zenodo.py \
  --url "https://zenodo.org/records/18511521/files/mavedb-dump.20260206153444.zip"
  --output data_sources/mavedb/tmp
python data_sources/mavedb/scripts/extract_metadata.py \
  --input data_sources/mavedb/tmp \
  --output data_sources/mavedb/tmp/mave_identifier.csv
```

The extractor also accepts the path to `main.json` directly. Use `--dry-run` to
print summary counts without writing the CSV.

## Extract variant scores

Run from `data_sources/mavedb` with the project virtual environment active. The
script reads `*.scores.csv` files from the dump's `csv/` directory and streams
one output row per variant. It splits variant numbers from accessions ending in
`#<number>`, detects the HGVS scheme in `hgvs_pro`, and parses simple protein
substitutions into reference amino acid, position, and alternate amino acid.
Missing score tokens are written as `NA`; `*.counts.csv` files are not read.

```bash
python scripts/extract_score.py \
  --input tmp/mavedb_dump/csv \
  --output tmp/mave_score.csv
```

The input can also be the dump directory `tmp/mavedb_dump`; the script selects
its `csv/` subdirectory automatically. With no arguments, it uses the paths
shown above.

## Process MAVE, MPRA and splicing data

`scripts/process_mavedb.py` reproduces the TSV analysis in
`ot_perturbation_catalog_mave_experiments/bin/process_mavedb_data.R`. Run from
the repository root after installing `data_sources/mavedb/requirements.txt`:

```bash
.venv/bin/python data_sources/mavedb/scripts/process_mavedb.py \
  --input data_sources/mavedb/tmp/mavedb_dump \
  --output data_sources/mavedb/output/processed
```

`--input` also accepts `main.json` directly. `--scores` overrides the sibling
`csv/` directory. The script caches NCBI's `categories.dmp` as
`data_sources/mavedb/tmp/taxcat.dmp`, downloading it from the NCBI taxonomy FTP
archive when it is missing. `--taxonomy` accepts an existing `taxcat.dmp` to
use instead, for example when running offline with a pinned copy. The cache is
reused on later runs. With no arguments, input and output use the paths shown
above, independently of the working directory.

The output directory contains:

| Directory | Files |
| --- | --- |
| `maves/` | `meta.tsv`, `meta_columns.tsv`, `full.tsv`, `singles.tsv`, `doubles.tsv`, `wide.tsv` |
| `mpras/` | `meta.tsv`, `meta_columns.tsv`, `full.tsv` |
| `splicing/` | `meta.tsv`, `meta_columns.tsv`, `full.tsv` |

Existing output files are replaced. The script processes one score file at a
time, ignores count files, drops missing scores, and excludes the same two
multi-gene localization score sets as the R workflow. Protein-coding rows with
protein HGVS calls enter MAVEs; regulatory targets enter MPRAs; protein-coding
rows with splice HGVS calls and `other_noncoding` targets enter splicing. A row
can enter both MAVEs and splicing.

Single mutants are aggregated after applying the UniProt offset, retaining the
mean score, replicate count and sample standard deviation. As in R,
`doubles.tsv` contains calls with **at most two mutations** from score sets
that contain multiple mutants, including their single mutants and WT controls.
`wide.tsv` includes all 20 standard amino acids, deletion and stop columns,
even for sparse inputs. It also accepts `*` stop notation. Missing values are
`NA`, logical values are `TRUE`/`FALSE`, and metadata text is quoted. DNA is
translated only when its length is divisible by three. Both legacy `taxId` and
newer `code` taxonomy fields are supported, and targets without external
identifiers or sequences are retained. No FASTA or sequence-mapping files are
written.

The wide-table column order follows the corrected R workflow (commit
`58eed32`): each score, count and standard-deviation block lists the 20 standard
amino acids, followed by deletion and stop. Synonymous scores are selected by
the reference amino acid's column name, so synonymous flags and derived
missense-completeness metrics do not depend on column order. The Python
wild-type code list also accepts `n.=` and `n.(=)` individually, alongside the
original R literal `n.=, n.(=)`. Sparse-input cases that R cannot process
(such as missing stop/deletion columns or nucleotide calls with no positions)
produce complete headers and `NA` metrics where needed.

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
