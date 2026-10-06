# MaveDB metadata post-processing

See the [study curation overview](../README.md) for the complete workflow and directory layout.

This package takes the sealed LLM metadata export and optional human curation
spreadsheet through the local MaveDB processing, term standardization,
validation, and parquet output stages. It keeps each run in a new output
directory. Manual values take precedence when both sources provide a value;
conflicts are written to a QC report.

The manual workbook may contain multiple rows for one `dataset_id` when a
score set contains multiple target genes. The merge keeps every manual row and
broadcasts dataset-level LLM metadata onto them. During score processing,
target-specific rows are assigned to variants by matching the gene prefix in
their HGVS identifiers to `perturbed_target_symbol`. LLM exports still need
one row per `dataset_id`. Conflict reports include the source workbook row
number so repeated dataset IDs can be distinguished.

Run commands from the repository root with the project environment active and
`data_exploration` on `PYTHONPATH`:

```bash
source .venv/bin/activate
export PYTHONPATH=data_exploration

python -m curation_tools.study_curation.mavedb.cli run \
  --llm-metadata curation_runs/my-run/exports/curation_final_metadata.csv \
  --manual-metadata data_exploration/MaveDB/mavedb_studies.xlsx \
  --output-dir curation_runs/my-run/postprocessing/default
```

Omit `--manual-metadata` for an LLM-only input. The selected output directory
must be new or empty. Per-dataset logs, merged input metadata, standardized
combined metadata, QC reports, parquet files, and the run manifest are saved
there. Term standardization runs after per-dataset processing and parquet
concatenation to preserve the existing stage order.

Pass `--save-joint-artifact` to also write
`final/mavedb_all_curated_data_with_metadata_postfilter.parquet`, joining each
score-data row to its validated metadata by `dataset_id` and, when a dataset
has target-specific metadata rows, by `perturbed_target_symbol` as well. This
output is omitted by default.

Term QC runs on the final standardized metadata. `qc/unmapped_terms.csv` reports
unrecognized labels and missing mapping-source columns, including fields that
have only ID lookups. `qc/missing_expected_ids.csv` separately reports recognized
labels whose expected IDs are absent. Both reports retain affected row counts
and an `issue_type` column. Checks respect dataset scopes, recognize normalized
mapping outputs, and treat explicit ID clearing as intentional. Treatment labels
and IDs are checked at corresponding pipe-delimited positions. The QC table
returned by `standardize_metadata` contains both report categories.

Review the final parquet files before publishing. BigQuery upload is a separate
command; it checks both destination schemas before updating the existing
`mavedb.metadata` and `mavedb.data` tables:

```bash
python -m curation_tools.study_curation.mavedb.cli publish \
  --metadata-parquet curation_runs/my-run/postprocessing/default/final/mavedb_all_curated_metadata_postvalidation.parquet \
  --data-parquet curation_runs/my-run/postprocessing/default/final/mavedb_all_curated_data_postfilter.parquet \
  --project-id prj-ext-dev-pertcat-437314
```

`resources/metadata_mappings.csv` contains controlled term lookups and direct
field/value set or clear rules. Equality-scoped overrides that fit this form
are stored here. Their `rule_id` and `reason` remain part of the QC audit.
`resources/metadata_overrides.csv` keeps rules that need other scope predicates
or conditional actions, such as prefix, not-equals, or fill-missing. The
notebook shows the intermediate QC tables and calls the same functions as the
CLI.

Mappings run once in ascending `mapping_order`; their physical row order is
irrelevant. Choose each new value to place the rule after the inputs it needs
and before any rule it should precede. `source_case_sensitive` preserves exact
matching for identifiers; other source values are case-insensitive. A
self-mapping (`source_field` and `target_field` are the same) must follow rules
that use the same source field and source value when their `dataset_id_prefix`
scopes overlap. Different source values and disjoint dataset scopes can be
ordered independently.
