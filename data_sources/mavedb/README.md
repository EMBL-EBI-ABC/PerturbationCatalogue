# MaveDB

## Background

MaveDB data can be downloaded from Zenodo and is licensed under CC0.

The ZIP archive contains:

- `main.json` — metadata
- `csv/` — data packaged per study
  - `*.scores.csv` — always present
  - `*.counts.csv` — sometimes present

## Nextflow workflow

The [Nextflow pipeline](pipeline/README.md) implements raw-data ingest,
reference preparation and canonical metadata/score processing. It accepts a
pinned Zenodo download or local ZIP, validates and unpacks it, and prepares a
single reference bundle containing taxonomy plus ClinVar and GWAS variants.
Default raw-source runs also generate singles, doubles and wide tables under
`<outdir>/tables/maves/`. Outputs include raw, reference and processed
provenance manifests and a dataset selection audit. Use `--ingest_only` to stop after ingest.
Local setup, offline fixtures, reuse behavior and validation status are documented there. GMM
fitting is available with `--run_gmm`, including frozen batch plans, per-batch
tasks and verified collection. Completed results are reused by default from
ordinary stage directories under `<outdir>/`, independently of Nextflow work
directories and `-resume`. Use `--rerun merge` to repackage fits,
`--rerun gmm,merge` to refit and merge, or `--run_gmm --rerun all` to regenerate
all enabled stages. Changed inputs/settings or damaged saved outputs stop and
require explicit reruns rather than silently refitting. Preserve the entire
output directory, including each stage's `pipeline_manifest.json` and the small
root-level `pipeline_state.json`. Fitted collections are in `gmm-fit/`; packaged
batch and merged results are in `gmm/`.

The default `tables` stage has its own cache, so generating missing tables or
using `--rerun tables` does not repeat compatible completed fits. Optional
`--run_qc` and `--run_reports` add basic QC figures and per-study GMM reports.
QC uses the normal singles output; `--run_reports` enables GMM products if they
do not yet exist. The downstream-only `--processed_dir` mode imports canonical
scores and metadata for GMM/report work and does not generate source tables.
An interrupted GMM attempt reuses completed batches when continued without a
new `--rerun gmm`.
Use `--run_validation` to validate GMM results against the prepared ClinVar and
GWAS tables. Homology remains a separate command.

### Pipeline reference inputs

Normal source-mode runs prepare ClinVar and GWAS variants from Open Targets
coding variants release 25.06, evidence and credible-set releases 25.12, and
the checked-in `analysis/data/mavedb_uniprot_to_ensembl.tsv` map. These source
downloads are approximately 3.4 GB. The pipeline saves the complete, verified
bundle under `<outdir>/references/`; the standalone downloader and converter
below keep their existing all-source and all-table defaults.

For offline preparation, pass `--reference_sources_manifest` with a verified
version-3 raw source manifest. To import already prepared variant tables, pass
both `--validation_clinvar` and `--validation_gwas`; a partial pair is rejected
before any stages start. The pair skips raw source acquisition and preparation.
`--run_validation` consumes the prepared tables automatically. Processed-input
GMM-only runs do not prepare references; processed-input validation does.

Reference products are reused without checking for newer provider data. Use
`--rerun references` to regenerate from retained source snapshots, or add
`--refresh_references` to fetch a new source snapshot explicitly. Refreshing
requires the references rerun. The reference stage preserves taxonomy source
provenance and its exact bytes across reruns unless a different snapshot is
selected deliberately.

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
draws per chain, a dataset-derived seed (`10101` for `urn:mavedb:00000001-a-1`),
and a 2,048 MiB memory budget. Outputs are isolated
under `data_sources/mavedb/output/gmm-smoke/smoke/maves/<model>/`; QC records
effective settings, source identities, exclusions and sampler/PSIS warnings.
The fixture is a fixed, source-row-preserving reservoir sample of 500 from the
3,180 eligible rows of `urn:mavedb:00000001-a-1`. Its attribution and selection
details are in [the fixture notes](pipeline/assets/gmm-smoke/README.md).

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
.venv/bin/python -c 'from cmdstanpy import install_cmdstan; install_cmdstan(version="2.37.0", cores=1, progress=True, verbose=True)'
export CMDSTAN="$HOME/.cmdstan/cmdstan-2.37.0"
```

CmdStan is installed persistently under `$HOME/.cmdstan/`, which CmdStanPy also
discovers automatically when `CMDSTAN` is unset. The export above selects the
validated version explicitly.

Provisioning CmdStan is an explicit developer setup step. The fitting runtime
does not download or install toolchains or dependencies. Compiled model
executables are cached and reused when their sources and settings match.

### Outputs, limits and failure behavior

Each completed batch writes `score{index}.tsv`, `density{index}.tsv`,
`fitted_params{index}.tsv`, `posterior{index}.tsv`, `fit_summary{index}.tsv`
and `loo{index}.tsv`. For catalogue inputs the score table also carries
`dataset_id` and `source_row`; repeated observations remain separate. Study QC
is written to `qc/batch{index}/`, including `study_diagnostics.tsv`,
`run_metadata.json`, `eligibility.tsv` and `exclusions.tsv` (and the selected
identities for smoke runs).

The Nextflow packaging step places these under
`maves/<model>/batches/batch<N>/` using the unnumbered table names, and writes
the concatenated tables under `maves/<model>/merged/`. Batch-level
`score.tsv` files retain both `dataset_id` and `source_row`; the merged
`score.tsv` omits the redundant `dataset_id` column. Standalone batch runs
retain the indexed filenames and columns described above.

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

Production uses four parallel chains, 1,000 warmup draws and 2,000 retained
draws per chain, `adapt_delta=0.9`,
`max_treedepth=10`, and a deterministic dataset-derived seed unless `--seed` is
supplied. Fitting is capped at 256,000 rows per study; `--max-fit-rows` can
lower that cap. Use `--max-study-rows` only when intentionally limiting both
fitting and post-processing. `--variants_per_batch` controls batch size, not
memory.

When `--seed` is omitted, each study derives its Stan seed from its MaveDB URN
(`urn:mavedb:<exp_set>-<exp_letter>-<score_set>` ->
`exp_set * 10000 + letter_idx * 100 + score_set`, e.g. `urn:mavedb:00000984-a-1`
-> `9840101`, `urn:mavedb:00000223-c-2` -> `2230302`). Passing `--seed 28`
overrides this and passes the exact integer `28` to Stan for every study and
candidate. Row subsampling uses separate deterministic study-specific seeds.

The modes workflow compares the `ko`, `del` and `neut` candidates, then applies
the `p_loo < 12` cutoff for winner selection. The `pos` fit is retained for
reporting. See the [GMM validation record](analysis/VALIDATION.md) for measured
resource use and the scope of the smoke checks.

## External reference data

To download reference data for standalone analysis, run the following command.
Generate the canonical metadata first if you need Ensembl homology.

```bash
.venv/bin/python data_sources/mavedb/scripts/download_reference_data.py --allow-missing-homology
```

It downloads Open Targets coding variants, ClinVar evidence, credible sets,
pharmacogenomics and drug molecules; downloads and unpacks the ProtVar
`2024.05.28_pockets.tsv.gz`; and retrieves aligned Ensembl protein homologies.
Open Targets is pinned to `25.12/output/`, except for the coding-variant
source at `25.06/view/coding_variant_table/`. Ensembl uses the live REST endpoint.

All inputs needed for downstream processing live together under `references/`
(paths in this section are relative to `data_sources/mavedb/`):

```text
references/
├── manifest.json
├── opentargets/
│   ├── 25.06/coding_variants/
│   └── 25.12/{evidence_eva,credible_sets,pharmacogenomics,drugs}/
├── protvar/2024.05.28/
├── ensembl/{homology.tsv,metadata.tsv,retrieval_report.json}
└── mappings/mavedb_uniprot_to_ensembl.tsv
```

The version-3 manifest records source URLs, exact partition lists, file sizes and
SHA-256 identities. Every file path is relative to the manifest, so the entire
directory can be moved or copied and passed to either processor using `--sources`.
The metadata and normalized ID mapping are snapshots of the inputs used to build
the homology table. The retrieval report records skipped and failed Ensembl queries.
Open Targets and ProtVar directories retain `download.json` transfer records for
verified reuse and offline operation. Custom source URLs without a release/date
use a `custom-<URL digest>` directory beneath their provider.

Only disposable Ensembl response caches live under `tmp/reference_data/`.
Downstream processing does not read that directory. The downloader reuses verified
local files by default; `--refresh` fetches them again, and `--offline` prohibits
network access. Offline homology regeneration requires the Ensembl response cache;
processing an already published bundle does not. Processing verifies the exact
files in the manifest, so obsolete partitions cannot enter an analysis. Source
bytes may be updated in place during a refresh; a failed refresh retains the old
manifest but changed inputs are rejected until a complete bundle is published.

`--out` selects the complete reference directory. `--cache-dir` selects only the
disposable Ensembl cache. `--metadata` defaults to `output/processed/metadata.tsv`;
`--id-map` defaults to the checked-in mapping under `analysis/data/`. These
determine the unique Ensembl gene/species requests.
Each source URL can be changed with its `--*-url` flag. Ensembl requests default
to four workers, a global 0.2-second request interval, a 60-second timeout and
three retries. The corresponding flags are `--workers`, `--request-interval`,
`--timeout` and `--retries`.

The Nextflow references stage selects only ClinVar and GWAS inputs with
`download_reference_data.py --only clinvar,gwas`; the downloader's standalone
default remains all source families. `process_variant_data.py --only
clinvar,gwas` writes both prepared tables in one pass with a shared coding
index; its standalone default remains all five output tables.

Unavailable Ensembl HTTP 400/404 responses stop the run unless
`--allow-missing-homology` is explicitly supplied. Skips are recorded in the
bundle's homology report. `--allow-incomplete-homology` separately permits
exhausted transient failures. `--skip-homology` downloads only Open Targets and
ProtVar, without requiring processed metadata or querying Ensembl. A failed run
keeps the previously published manifest and homology tables; completed transfers
remain available for the next attempt.

Acquisition is implemented through three library modules:
`_reference_sources.py` defines sources and reads the bundle manifest;
`_reference_downloads.py` owns HTTP, retry, cache and gzip/Parquet validation;
`_ensembl_homology.py` plans, retrieves and flattens Ensembl homologies. The
downloader validates the homology plan before starting the other network transfers.
`_reference_io.py` shares streaming provenance hashes and transactional output
publication. The downloader calls these libraries directly.

## Variant, drug and pocket reference tables

After downloading, generate all five reference tables with:

```bash
.venv/bin/python data_sources/mavedb/scripts/process_variant_data.py
```

The processor reads `references/manifest.json` by default and performs
only local processing. It has no download, refresh or network fallback modes.
Use `--sources /path/to/references/manifest.json` for a different bundle, `--out` for a
different output directory, and `--batch-size` to control source scanning.
`--only` selects `clinvar`, `gwas`, `pharmacogenomics`, `drugs`, or `pockets`;
the default is `all`.

The default destination `output/variant_data/` contains:

- `mave_clinvar_variants.tsv`
- `mave_gwas_variants.tsv`
- `mave_pharmacogenomics.tsv`
- `mave_drugs.tsv`
- `mave_pockets.tsv`
- `variant_data_manifest.json`

ClinVar is rebuilt from raw evidence and coding-variant Parquet files on every
run. The processor constructs one temporary SQLite coding index for all selected
coding-dependent stages, rather than rescanning the source separately for
ClinVar, GWAS and pharmacogenomics. Parquet is scanned in bounded batches.
GWAS retains global coding-variant matches; pharmacogenomics matches target and
variant before filtering mapped UniProt accessions. Duplicate joins and nested
list expansions retain the original biological row semantics. Pocket confidence
uses the strict `score > 800` threshold. Missing values use `NA`, booleans use
`TRUE/FALSE`, and remaining nested cells use compact JSON.

Selected outputs are staged before publication, with rollback on I/O failure.
Unrelated files are preserved. Default paths resolve from the repository;
explicit relative paths resolve from the calling directory.

## GMM validation against controls, ClinVar and GWAS

`analysis/mave_gmm_validation.py` reads completed GMM fits, canonical
metadata/scores, and the two variant TSVs generated by
`scripts/process_variant_data.py`. It performs local analysis without fitting
models or downloading references.

```bash
.venv/bin/python data_sources/mavedb/analysis/mave_gmm_validation.py \
  --gmm-dir /absolute/path/to/completed/maves/modes \
  --out data_sources/mavedb/output/gmm_validation
```

The defaults are `output/processed/metadata.tsv`,
`output/processed/scores.tsv`, `output/variant_data/mave_clinvar_variants.tsv`
and `output/variant_data/mave_gwas_variants.tsv`, all below this MaveDB directory.
Override them with `--metadata`, `--source-scores`, `--clinvar` and `--gwas`.
Default paths resolve from the repository. Explicit relative paths resolve from
the calling directory. Pass the packaged model directory (`maves/<model>/`) or
a standalone batch directory with matching `score<N>.tsv` and
`fitted_params<N>.tsv` files. Validation reads batch score tables with both
`dataset_id` and `source_row` to recover WT flags from canonical source scores.
It selects top models and retains each original observation. It derives
synonymous and nonsense flags from protein HGVS and applies the metadata offset
to get UniProt positions.

The three benchmarks compare nonsense controls against WT/synonymous controls,
pathogenic/likely pathogenic ClinVar variants against benign/likely benign
variants, and GWAS matches against unmatched variants. ClinVar labels are
combined per protein/substitution.
Uncertain or conflicting labels remain in matched outputs and distributions but
are excluded from binary evaluation. GWAS is restricted to human score sets,
and only score sets with a hit enter its per-study and pooled ROC calculations.
Repeated observations and multiple GWAS associations retain their join
multiplicity. Matchable variants without a GWAS hit form the benchmark negative
class. Direct matching uses UniProt, amino-acid position, reference residue and
alternate residue.

| Predictor | Formula | Benchmarks |
| --- | --- | --- |
| `s` | Raw fitted input score | All |
| `ns` | `(s - median) / mad` | All |
| `w` | GMM-standardised score | All |
| `aw` | `abs(w - 1)` | ClinVar, GWAS |
| `p` | `w / sigma_neut` | All |
| `ap` | `abs(w - 1) / sigma_neut` | ClinVar, GWAS |

Lower signed scores and larger absolute deviations predict the positive class.
The benchmarks thin thresholds and integrate ROC and precision–recall
curves with left steps. Precision–recall figures plot recall on the x-axis.
Matplotlib uses a symmetric-log raw-score axis, and density overlays retain all
observations.

The default `output/gmm_validation/` contains eight TSVs:
`controls.tsv`, `clinvar_matches.tsv`, `gwas_matches.tsv`, `ctrl_roc.tsv`,
`clinvar_roc.tsv`, `gwas_roc.tsv`, `auc_summary.tsv` and `counts.tsv`.
The manifest records input and implementation SHA-256 hashes, package versions,
coverage, settings, warnings and output files. `auc_summary.tsv` includes the
per-predictor class counts; missing metrics remain `NA`. Ten figures cover raw
and standardised distributions plus each benchmark's labelled distributions,
ROC/PR curves and AUC distributions. They are written under `figures/` in PNG
and PDF by default; select `--formats png`, `--dpi 150` or `--no-plots` as needed.
Managed files are staged before publication with rollback on I/O failure;
unrelated files are preserved and obsolete managed figures are removed on rerun.
Wait for successful completion before reading a report.

Empty clinical or GWAS matches produce explicit coverage warnings and labelled
empty figures. Select fits with known overlap with the chosen references to
evaluate clinical or GWAS performance.

### Nextflow validation stage

The Nextflow pipeline can run the same validation as an optional `validation`
stage. Enable it with `--run_validation`; this also enables GMM fitting and
merge, reusing compatible saved products. The stage consumes the ClinVar and
GWAS TSVs from the pipeline's `references/variants/` bundle. To use prepared
tables, supply both `--validation_clinvar` and `--validation_gwas`. The pipeline
records the content identities of both tables for reuse.

The stage publishes its eight TSV tables and `manifest.json` under
`<outdir>/validation/`. Plots are written under `validation/figures/` by
default. `--validation_formats` accepts comma-separated `png`, `pdf`, or both
(default `png,pdf`); `--validation_dpi` sets raster resolution (default `150`),
and `--validation_no_plots` skips figures. Unchanged inputs, options and
implementation reuse the completed stage. Changes to references, options or
implementation require `--rerun validation` to recompute it.

## Homology comparisons and ClinVar evaluation

`analysis/mave_homology.py` consumes native fitted GMM batches, canonical
processed metadata and the downloaded source bundle. It compares standardized
variant effects between score sets for the same protein and aligned homologs,
and evaluates ClinVar labels transferred across those alignments.

```bash
.venv/bin/python data_sources/mavedb/analysis/mave_homology.py \
  --gmm-dir data_sources/mavedb/output/gmm-production/maves/modes
```

The analysis reads `references/manifest.json` and the generated
`output/variant_data/mave_clinvar_variants.tsv` by default. It requires completed
native score and fitted-parameter tables and canonical metadata. Pass the
collected model directory (`maves/<model>/`), its `merged/` directory, or a
standalone batch directory containing matching `score<N>.tsv` and
`fitted_params<N>.tsv` pairs. `--sources` selects another bundle and
`--clinvar` selects another generated ClinVar table. `--skip-clinvar`
explicitly runs only the between-study analysis.

Metadata coverage and fitted-score coverage are distinct. One fitted study
cannot produce between-study correlations. Two studies require at least three
overlapping substitutions with variation on both sides. The two UBE2I studies
`urn:mavedb:00000001-a-1` and `urn:mavedb:00000001-a-3` provide a small
same-protein comparison example. Select completed fits with ClinVar overlap for
clinical evaluation.

Outputs under `output/homology/` include `correlations.tsv`,
`combined_clinvar.tsv`, `clinvar_homologs.tsv`, counts, ROC curves, summary tables,
a manifest and figures. Comparisons report Pearson confidence intervals,
asymptotic Spearman tests, fixed-category Cohen's kappa and Benjamini–Hochberg
corrections. ROC and precision–recall areas use left-step integration. Undefined
metrics are missing.

The four ClinVar predictors are `w`, `abs(w - 1)`, `w / sigma_neut` and
`abs(w - 1) / sigma_neut`. Lower signed scores or larger absolute deviations
predict pathogenicity. Counts refer to matched observations and retain
uncertain/conflicting labels until binary classification is calculated.
Plot layout uses Matplotlib.

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

`scripts/process_mavedb.py` extracts score-set metadata and variant records
and computes MAVE, MPRA and splicing summaries.

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

Use `--taxonomy /path/to/taxcat.dmp` to supply a pinned NCBI taxonomy snapshot.
Otherwise, standalone processing uses `tmp/taxcat.dmp` and downloads it when
missing.

The default output contains these two canonical tables, joined by **`dataset_id`**,
and the singles, doubles and wide tables described below:

| File | Grain | Contents |
| --- | --- | --- |
| `scores.tsv` | One original variant record | Dataset ID, original variant number, source row number, nucleotide/protein/splice HGVS, simple-protein-variant parsing, amino-acid codes, `wt_mut` and score |
| `metadata.tsv` | One dataset | Source and unpacked target metadata, analysis types, variant counts and derived statistics including coverage/completeness |

Pass `--audit` to also write `dataset_processing.tsv`, containing one record per
input score file with inclusion/exclusion status and reason, and validate the
canonical outputs before publication. Nextflow calls this processor directly
with `--ingest-manifest PATH`, enabling the audit, output checks and provenance
manifest together. That mode requires an explicit `--taxonomy` input. Its
canonical stage passes `--no-tables` because the separate default `tables` stage
owns those products; standalone processing writes them unless `--no-tables` is
requested.

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

`is_simple_p`, `ref_aa_abbr`, `position` and `alt_aa_abbr` describe simple
protein HGVS calls in source coordinates, with three-letter amino-acid values
and literal `=` and `*` symbols. Scores are numeric;
missing values are `NA`. Invalid numeric values are retained in `score_raw` and
excluded from analysis. Custom auxiliary score/count columns are not exported.

`ref_aa_sym` and `alt_aa_sym` add single-letter amino-acid codes. Synonymous
`=` uses the reference amino-acid letter. Both `Ter` and `*` use `stop`.
Both columns are `NA` for protein calls outside the simple-HGVS parser,
including multi-substitution calls.

`wt_mut` is true when `hgvs_pro` matches one of the `WT_HGVS_CODES` values:
`p.=`, `p.(=)`, `_wt`, `_sy`, `n.=`, `n.(=)` or `n.=, n.(=)`.

### Non-redundant source metadata

Source fields include the score-set title, description, abstract, reported
variant count and target fields. The target fields record names, categories,
sequences, labels, types, taxonomy, external identifiers and offsets. Experiment and
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
then MPRA, then splicing. For a single primary target, the flat columns contain
score summaries, mutation summaries, positional coverage and completeness.
Statistics use valid scores from that analysis
subset, with sample standard deviations and no score normalization.

Secondary analyses and multi-target positional metrics are omitted from
`metadata.tsv`; no `target_<analysis>_<metric>` columns are exported.
Multi-target datasets have no pooled scalar coverage: their variant
counts and general score summary count each original variant once.
`position_reference` states whether primary positional metrics use a known
UniProt offset or source coordinates. All derived target fields follow the
source columns, alongside the existing derived statistics.

### Singles, doubles and wide tables

Standalone processing generates these standard outputs by default:

```text
maves/singles.tsv
maves/doubles.tsv
maves/wide.tsv
```

These tables identify datasets with `score_urn`; the canonical tables use
`dataset_id`. Singles are averaged by mutation after applying
known UniProt offsets. Doubles include calls with at most two mutations from
score sets containing multiple mutants, including single mutants and WT
controls. Wide columns list the 20 amino acids, then `del` and `stop`, within
each score/count/standard-deviation block. Synonymous scores are selected by
reference amino-acid column name and do not depend on column order.

Pass `--no-tables` to generate only the canonical tables in standalone processing.

In Nextflow, the separate default `tables` stage generates these outputs under
`<outdir>/tables/maves/`, using the matching raw source, taxonomy and canonical
selection through `scripts/generate_tables.py`. It reuses the same aggregation
helpers and does not rewrite canonical tables or invalidate existing GMM fits.
Use `--rerun tables` to regenerate this stage explicitly. Raw-source runs adopt
or generate missing tables automatically. External `--processed_dir` runs
bypass this stage because they do not supply its raw-source inputs.

### Per-study GMM summary reports

`scripts/summary_plot_batch.py` uses the plotting dependencies in
`requirements.txt` to generate reports from completed fits. Run it from the
repository root.

```bash
.venv/bin/python -m pip install -r data_sources/mavedb/requirements.txt
.venv/bin/python data_sources/mavedb/scripts/summary_plot_batch.py \
  --model modes \
  --gmm-dir data_sources/mavedb/output/gmm-production/maves/modes \
  --meta data_sources/mavedb/output/processed/metadata.tsv \
  --index 1
```

Pass the collected model directory (`maves/<model>/`), its `merged/` directory,
or a standalone batch directory containing the fitter's indexed `score<N>.tsv`,
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
  --meta data_sources/mavedb/pipeline/assets/gmm-smoke/metadata.tsv
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
empty-data panels. Per-study reports cover protein-coding MAVE score sets.

### MAVE summary and QC figures

`analysis/mave_qc.py` reads the canonical `scores.tsv` and `metadata.tsv`, plus
`maves/singles.tsv` table. Standalone processing writes this table by default;
Nextflow's raw-source runs publish it at `<outdir>/tables/maves/singles.tsv`.
Use `--run_qc` to add QC independently of GMM fitting. The Nextflow figures are
published under `<outdir>/qc/`.

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
score sets, and correlations with the BLOSUM100 matrix bundled in
`analysis/data/BLOSUM100.txt`.

Scores are selected for the MAVE analysis using finite numeric scores and
nonmissing protein HGVS, rather than counting every exported canonical record.
Meta-analysis score sets are omitted from coverage, per-variant counts, and
between-score-set correlations. Multi-target datasets remain
in global dataset counts; ambiguous target assignments are logged and omitted
from target-specific comparisons rather than duplicating their score rows.
