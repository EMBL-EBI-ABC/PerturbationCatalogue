# MaveDB

## Background

MaveDB data can be downloaded from Zenodo and is licensed under CC0.

The ZIP archive contains:

- `main.json` — metadata
- `csv/` — data packaged per study
  - `*.scores.csv` — always present
  - `*.counts.csv` — sometimes present

## GMM fitting

### Run the bounded smoke workflow

After the one-time environment setup below, these commands run real CmdStan
fits against the bundled 500-row fixture. They do not read the full processed
catalogue.

```bash
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile smoke --model flexible \
  --out data_sources/mavedb/output/gmm-smoke

.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile smoke --model modes \
  --out data_sources/mavedb/output/gmm-smoke

# Print the smoke selection and batch plan without fitting or creating outputs.
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile smoke --list-batches \
  --out data_sources/mavedb/output/gmm-smoke
```

The smoke profile selects at most one study and 500 rows for fitting and
post-processing. It uses two chains sequentially, 200 warmup and 200 retained
draws per chain, seed `123`, and a 2,048 MiB memory budget. Outputs are isolated
under `data_sources/mavedb/output/gmm-smoke/smoke/maves/<model>/`; QC records
effective settings, source identities, exclusions and sampler/PSIS warnings.
The fixture is a fixed, source-row-preserving reservoir sample of 500 from the
3,180 eligible rows of `urn:mavedb:00000001-a-1`. Its attribution and selection
details are in [the fixture notes](tests/fixtures/gmm/README.md).

The CLI resolves its bundled smoke fixtures and default canonical production
paths from the repository root, even when invoked from another working
directory. To invoke it elsewhere, use absolute paths for both the virtual
environment and script, for example:

```bash
repo=/path/to/PerturbationCatalogue
"$repo/.venv/bin/python" "$repo/data_sources/mavedb/scripts/fit_gmm_batch.py" \
  --profile smoke --model flexible --out "$repo/data_sources/mavedb/output/gmm-smoke"
```

Relative paths supplied explicitly with `--scores`, `--meta` or `--out` are
resolved from the current working directory; use absolute paths when launching
from elsewhere.

### One-time local setup

Use the repository's `.venv`. The validated environment is Python 3.12.8,
NumPy 2.2.6, pandas 2.2.3, SciPy 1.14.1, xarray 2025.3.0, ArviZ 0.22.0,
CmdStanPy 1.2.5 and stanio 0.5.1, with CmdStan 2.37.0. On the validated Mac,
Apple clang++ 21.0.0 was available for `arm64-apple-darwin25.6.0`. If `clang++`
is missing, install the Apple command-line developer tools with
`xcode-select --install`, then run from the repository root:

```bash
xcode-select -p
clang++ --version
.venv/bin/python -m pip install -r data_sources/mavedb/requirements-gmm.txt
mkdir -p data_sources/mavedb/tmp/cmdstan
.venv/bin/python -c 'from cmdstanpy import install_cmdstan; install_cmdstan(version="2.37.0", dir="data_sources/mavedb/tmp/cmdstan", cores=1, progress=True, verbose=True)'
export CMDSTAN="$PWD/data_sources/mavedb/tmp/cmdstan/cmdstan-2.37.0"
```

Provisioning CmdStan is an explicit developer setup step. The fitting runtime
does not download or install toolchains or dependencies. Stan executables are
cached below the ignored `data_sources/mavedb/tmp/` directory and are reused
when their sources/settings match.

### Outputs, limits and failure behavior

Each completed batch writes `score{index}.tsv`, `density{index}.tsv`,
`fitted_params{index}.tsv`, `posterior{index}.tsv`, `fit_summary{index}.tsv`
and `loo{index}.tsv`. For catalogue inputs the score table also carries
`dataset_id` and `source_row`; repeated observations remain separate. Study QC
is written to `qc/batch{index}/`, including `study_diagnostics.tsv`,
`run_metadata.json`, `eligibility.tsv` and `exclusions.tsv` (and the selected
identities for smoke runs).

Expected study-level numerical failures yield one raw-score placeholder per
eligible input row and identifier-preserving sentinel rows in the other five
tables; all required candidates must succeed for a study to be successful.
A metadata-eligible study with no valid score/protein-HGVS rows is retained as
a failed study with zero score rows, QC exclusions and compact-table
placeholders. Global input, dependency, compile and output errors exit
nonzero. The six TSVs are staged and validated before replacement, but separate
filesystem replacements are not a multi-file transaction. Downstream jobs
should wait for a successful process exit before reading the files. Convergence
and Pareto-k warnings remain visible in QC and do not change model selection.

### Batch planning and production inputs

Production defaults use the processed canonical `scores.tsv` and `metadata.tsv`
tables. Start by generating the current plan; it counts analysis-eligible rows,
orders studies by decreasing count then full ID, and assigns whole studies to
deterministic batches. Positive `--index` values are 1-based; `--index 0` means
all selected studies in one run. Oversized studies occupy their own batch.
`--list-batches` reports the plan without fitting or creating outputs.

```bash
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile production --model modes --list-batches --index 0 \
  --scores data_sources/mavedb/output/processed/scores.tsv \
  --meta data_sources/mavedb/output/processed/metadata.tsv

# Run one planned batch.
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile production --model modes --index 1 --variants_per_batch 50000 \
  --scores data_sources/mavedb/output/processed/scores.tsv \
  --meta data_sources/mavedb/output/processed/metadata.tsv \
  --out data_sources/mavedb/output/gmm-production

# Or explicitly request every planned study in one invocation.
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile production --model modes --index 0 \
  --scores data_sources/mavedb/output/processed/scores.tsv \
  --meta data_sources/mavedb/output/processed/metadata.tsv \
  --out data_sources/mavedb/output/gmm-production
```

Production defaults retain the R sampler settings: four chains and four
parallel chains, 1,000 warmup and 2,000 retained draws, `adapt_delta=0.9`,
`max_treedepth=10`, and no seed unless one is supplied. Fitting is capped at
256,000 rows per study; `--max-fit-rows` can lower that cap. Use
`--max-study-rows` only when intentionally limiting both fitting and
post-processing. `--variants_per_batch` controls batch size, not memory.

The default input format is `catalogue`. For the legacy score/meta tables,
provide both files explicitly and select `--input-format legacy`:

```bash
.venv/bin/python data_sources/mavedb/scripts/fit_gmm_batch.py \
  --profile production --input-format legacy --model modes --index 1 \
  --scores /absolute/path/to/legacy-scores.tsv \
  --meta /absolute/path/to/legacy-metadata.tsv \
  --out /absolute/path/to/gmm-production
```

The Python port preserves the R transforms, output-density constants,
classification order, candidate comparison and `p_loo < 12` behavior. It keeps
Stan's distinct likelihood constants and the positive candidate's exclusion
from winner selection. See the detailed [GMM acceptance record](analysis/VALIDATION.md)
and [implementation plan](../../development_docs/mavedb-gmm-python-implementation-plan.md)
for parity evidence, measured local resources and the remaining production
gate.

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
| `scores.tsv` | One original variant record | Dataset ID, original variant number, source row number, nucleotide/protein/splice HGVS, simple-protein-variant parsing, amino-acid codes, `wt_mut` and score |
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

`wt_mut` is true when `hgvs_pro` matches one of the `WT_HGVS_CODES` values:
`p.=`, `p.(=)`, `_wt`, `_sy`, `n.=`, `n.(=)` or `n.=, n.(=)`.

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

### Per-study GMM summary reports

The Python replacement for
`ot_perturbation_catalog_mave_experiments/bin/summary_plot_batch.R` lives here:

- `scripts/summary_plot_batch.py`: batch report CLI.
- `analysis/mave_summary/`: input adapters, report layout, and GMM diagnostics.
- `tests/summary*_test.py`: plotting and input regression tests.

It uses the plotting dependencies already listed in `requirements.txt`; R,
CmdStan, and neural-network packages are not needed to generate reports from
completed fits. From the repository root:

```bash
.venv/bin/python -m pip install -r data_sources/mavedb/requirements.txt
.venv/bin/python data_sources/mavedb/scripts/summary_plot_batch.py \
  --model modes \
  --gmm-dir data_sources/mavedb/output/gmm-production/maves/modes \
  --meta data_sources/mavedb/output/processed/metadata.tsv \
  --index 1
```

One of `--gmm-dir` or `--data-root` is required. For native outputs, pass the
directory containing the fitter's completed `score<N>.tsv`,
`density<N>.tsv`, `fitted_params<N>.tsv`, `posterior<N>.tsv`,
`fit_summary<N>.tsv`, and `loo<N>.tsv` files. The loader selects the best
candidate for each study while retaining all candidates for the LOO panel.
It accepts canonical metadata keyed by `dataset_id` and derives mutation
heatmaps from fitted score rows, averaging repeated observations and applying
unambiguous UniProt offsets. No separate singles export is required. Missing or
ambiguous target identifiers are not used to group correlation comparisons.

For the bounded smoke fits described above:

```bash
.venv/bin/python data_sources/mavedb/scripts/summary_plot_batch.py \
  --model flexible \
  --gmm-dir data_sources/mavedb/output/gmm-smoke/smoke/maves/flexible \
  --meta data_sources/mavedb/tests/fixtures/gmm/smoke/metadata.tsv
```

`--index 0` (the default) plots all available score sets. Positive indices select
consecutive, non-overlapping groups of 100 studies in input order;
`--batch-size` changes that group size. This is a report batch index, independent
of the fitting batch numbers in filenames. Use one coherent set of fitting
batches: overlapping studies across batch files are rejected rather than
counted twice.

The default output root is `data_sources/mavedb/output/figures`, with reports
under `maves/studies_<model>/`. Each score-set URN becomes a PDF filename with
colons and hyphens replaced by underscores. Override the destination using
`--out`. Metadata and output defaults resolve from the repository, while
explicitly supplied relative paths resolve from the calling directory.

Reports contain original and standardized mutation heatmaps, the fitted score
transformation and component densities, posterior parameter histograms colored
by Rhat, LOO model comparisons, and correlations with related score sets.
Failed fits retain their raw-score plots; studies with no score rows receive
empty-data panels. The layout uses Matplotlib and can differ visually from R.
MPRA and splicing reports remain unsupported.

The original consolidated R exports remain supported through an explicit
`--data-root PATH`: `PATH/maves/meta.tsv` plus the seven tables in
`PATH/maves/gmm_<model>/` (`scores`, `singles`, `density`, `posterior`,
`fit_summary`, `params`, and `loo`, each with a `.tsv` suffix).

Run the focused regression tests with:

```bash
.venv/bin/python -m pytest data_sources/mavedb/tests/summary*_test.py
```

### MAVE summary and QC figures

The Python port of `ot_perturbation_catalog_mave_experiments/analysis/mave_qc.R`
reads the canonical `scores.tsv` and `metadata.tsv`, plus the optional
`maves/singles.tsv` compatibility table. Generate the latter with
`process_mavedb.py --compatibility-tables` when processing the dump.

```bash
.venv/bin/python -m pip install -r data_sources/mavedb/requirements.txt
.venv/bin/python data_sources/mavedb/analysis/mave_qc.py
```

Defaults resolve relative to this module, so the command also works from another
working directory. Override inputs or the figure destination explicitly:

```bash
.venv/bin/python data_sources/mavedb/analysis/mave_qc.py \
  --scores data_sources/mavedb/output/processed/scores.tsv \
  --metadata data_sources/mavedb/output/processed/metadata.tsv \
  --singles data_sources/mavedb/output/processed/maves/singles.tsv \
  --output data_sources/mavedb/output/figures/maves/qc
```

The default run writes 13 figures in both PNG and PDF (26 files), replacing
existing figures with the same names. Use `--formats png` or `--formats pdf`
for a single format, `--dpi` for raster resolution, and `--chunk-size` for
input read batches. It does not modify the processed data or export analysis
tables. Neither `doubles.tsv` nor `wide.tsv` is required.

Figures cover dataset/organism counts, study sizes, data per gene and variant,
completeness and positional coverage, mutation backgrounds, correlations between
score sets, and correlations with BLOSUM100. The exact Biostrings BLOSUM100
matrix is bundled in `analysis/data/BLOSUM100.txt`; R is not needed at runtime.
The four per-variant panels use their respective counts, correcting the original
R script's repeated experiment-set axis.

Scores are selected for the MAVE analysis using finite numeric scores and
nonmissing protein HGVS, rather than counting every exported canonical record.
Meta-analysis score sets are omitted from coverage, per-variant counts, and
between-score-set correlations, as in the R script. Multi-target datasets remain
in global dataset counts; ambiguous target assignments are logged and omitted
from target-specific comparisons rather than duplicating their score rows.

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
