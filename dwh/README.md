# Data Warehouse Pipeline

Automated pipeline for transforming and loading data from BigQuery to Postgres and Elasticsearch.

> **Before running any commands in this document**, set up the environment by running `pc_secrets dev`.

## Pipeline stages

The pipeline runs five stages sequentially:

| Stage | Directory | Description | Duration |
|-------|-----------|-------------|----------|
| 1. **Open Targets reference** | (Native) | Loads Open Targets Platform targets into the configured BQ reference table | ~minutes |
| 2. **dbt** | `bq_dbt/` | Transforms source BQ tables into final data mart tables | ~minutes |
| 3. **BQ → Postgres** | `bq_to_postgres/` | Loads final BQ data tables into Cloud SQL (Postgres) | ~hours |
| 4. **BQ → Elastic** | `bq_to_elastic/` | Loads summary tables into Elasticsearch | ~minutes |
| 5. **Release artifacts** | `release/` | Clusters source data, then fans out one Cloud Run task per dataset | variable |

Each stage depends on the previous one. If any stage fails, the pipeline stops.
For the ENSG dev stack, the Open Targets reference stage writes to
`BQ_DATASET.opentargets_targets`.

## Prerequisites

### 1. Google Cloud SDK

Install and authenticate the [gcloud CLI](https://cloud.google.com/sdk/docs/install):
```bash
gcloud auth login
gcloud auth application-default login
```

### 2. Enable required APIs

```bash
gcloud services enable cloudbuild.googleapis.com --project=$GCLOUD_PROJECT
gcloud services enable compute.googleapis.com --project=$GCLOUD_PROJECT
gcloud services enable servicenetworking.googleapis.com --project=$GCLOUD_PROJECT
```

### 3. Environment variables

The trigger script requires `GCLOUD_PROJECT`, `GCLOUD_REGION`, `BQ_DATASET`,
`BQ_LOCATION`, `CLOUD_TMP_BUCKET` (or legacy `GCLOUD_TMP_BUCKET`),
`PG_CONN_INTERNAL`, `ES_URL`, `ES_USERNAME` and `ES_PASSWORD`. Source and verify
the environment for the intended development deployment before invoking the
trigger; it prints the project and region and immediately submits the build. The
trigger refuses project IDs containing `prod`.

After `pc_secrets dev`, normalize the bucket variable for the manual commands:

```bash
export CLOUD_TMP_BUCKET="${CLOUD_TMP_BUCKET:-${GCLOUD_TMP_BUCKET:-}}"
test -n "$CLOUD_TMP_BUCKET"
```

Each build writes artifacts under its own temporary prefix:
`gs://$CLOUD_TMP_BUCKET/release/$BUILD_ID/`. Its manifest is under
`release-staging/$BUILD_ID/`. Existing objects under `release/` do not block a
build: Cloud Storage prefixes are virtual, and the preflight checks only that
this unique `release/$BUILD_ID/` path has no objects from an earlier run. A
failed build removes only its own staging and `release/$BUILD_ID/` objects.
Release staging tables and the Cloud Run Job are removed after the build. A
successful run's artifacts remain in the temporary prefix for manual review,
along with its small `release-staging/$BUILD_ID/manifest.json`. No global
`release/` prefix needs to be empty.

The release stage clusters source data by `dataset_id`, then runs one Cloud Run
Job task per selected dataset. It writes metadata JSON, CSV.GZ and Parquet files
for CRISPR and MAVE. Perturb-seq also gets `*.gsea.csv.gz` and `*.gsea.parquet`
alongside its DEA files. To review a specific build, list only its Perturb-seq
objects:

```bash
gcloud storage ls --long --recursive \
  "gs://$CLOUD_TMP_BUCKET/release/$BUILD_ID/perturb-seq/"
```

Compare the selected IDs with the successful run manifest and verify five
nonempty objects for each Perturb-seq dataset: metadata JSON, DEA CSV.GZ and
Parquet, and GSEA CSV.GZ and Parquet. A scoped run can be checked against its
input manifest with:

```bash
diff -u \
  <(jq -r '.datasets[].dataset_id' "$MANIFEST" | sort -u) \
  <(gcloud storage cat \
      "gs://$CLOUD_TMP_BUCKET/release-staging/$BUILD_ID/manifest.json" \
      | jq -r '.items[] | select(.modality == "perturb-seq") | .dataset_id' | sort -u)
```

After manual review, copy only the IDs in the input manifest to the shared
serving bucket under the release version prefix. For this release, the
development backend must use `RELEASE_BUCKET=perturbation-catalogue-release`
and `RELEASE_VERSION_PREFIX=2026.10`. Production keeps using the same bucket's
existing unprefixed paths by leaving `RELEASE_VERSION_PREFIX` unset. Copying
the versioned objects therefore leaves the current production artifacts in
place.

```bash
MANIFEST=data_sources/perturb-seq/pipeline/manifests/single-condition-20.json
: "${RELEASE_VERSION_PREFIX:=2026.10}"
files=()
while IFS= read -r dataset_id; do
  for suffix in metadata.json csv.gz parquet gsea.csv.gz gsea.parquet; do
    files+=("gs://$CLOUD_TMP_BUCKET/release/$BUILD_ID/perturb-seq/$dataset_id.$suffix")
  done
done < <(jq -r '.datasets[].dataset_id' "$MANIFEST")
gcloud storage cp "${files[@]}" \
  "gs://perturbation-catalogue-release/$RELEASE_VERSION_PREFIX/perturb-seq/"
```

For Perturb-seq, finish publication checks before setting reprocessed markers.
Run from the repository root. For this run, set `MANIFEST` to
`data_sources/perturb-seq/pipeline/manifests/single-condition-20.json`. The
Postgres postflight uses the external dev connection from `pc_secrets dev`
(`PG_CONN`); the Cloud Build trigger separately uses `PG_CONN_INTERNAL` to
reach Cloud SQL over its VPC connection.

```bash
MANIFEST=data_sources/perturb-seq/pipeline/manifests/single-condition-20.json
python3 - "$MANIFEST" <<'PY'
import json
import os
import sys

import psycopg2

with open(sys.argv[1]) as source:
    dataset_ids = [item["dataset_id"] for item in json.load(source)["datasets"]]
if len(dataset_ids) != 20 or len(set(dataset_ids)) != 20:
    raise SystemExit("Expected 20 unique dataset IDs")
with psycopg2.connect(os.environ["PG_CONN"]) as connection:
    with connection.cursor() as cursor:
        for table in ("perturb_seq_dea", "perturb_seq_gsea"):
            cursor.execute(
                f"SELECT dataset_id, COUNT(*) FROM {table} "
                "WHERE dataset_id = ANY(%s) GROUP BY dataset_id",
                (dataset_ids,),
            )
            counts = dict(cursor.fetchall())
            if set(counts) != set(dataset_ids) or any(count < 1 for count in counts.values()):
                raise SystemExit(f"Missing rows in {table}: {set(dataset_ids) - set(counts)}")
print("All selected IDs have DEA and GSEA rows in Postgres")
PY
```

Then verify the API and release artifacts. These checks confirm results and
per-dataset download links exist; the signed URLs can be opened manually to
confirm the downloaded files.

```bash
MANIFEST=data_sources/perturb-seq/pipeline/manifests/single-condition-20.json
: "${RELEASE_VERSION_PREFIX:=2026.10}"
API_BASE=${DEV_API_URL:-http://127.0.0.1:8000}
while IFS= read -r dataset_id; do
  curl -fsS "$API_BASE/v1/perturb-seq/$dataset_id/search?limit=1" \
    | jq -e '.total_rows_count > 0 and (.results | length) == 1' >/dev/null
  curl -fsS "$API_BASE/v1/perturb-seq/$dataset_id/gsea?limit=1" \
    | jq -e '.total_rows_count > 0 and (.results | length) == 1' >/dev/null
  for url in \
    "$API_BASE/v1/perturb-seq/$dataset_id/download?format=parquet" \
    "$API_BASE/v1/perturb-seq/$dataset_id/download?format=csv.gz" \
    "$API_BASE/v1/perturb-seq/$dataset_id/gsea/download?format=parquet" \
    "$API_BASE/v1/perturb-seq/$dataset_id/gsea/download?format=csv.gz"; do
    status=$(curl -sS -o /dev/null -w '%{http_code}' "$url")
    test "$status" -eq 307
  done
  for suffix in metadata.json csv.gz parquet gsea.csv.gz gsea.parquet; do
    size=$(gcloud storage objects describe \
      "gs://perturbation-catalogue-release/$RELEASE_VERSION_PREFIX/perturb-seq/$dataset_id.$suffix" \
      --format='value(size)')
    test "$size" -gt 0
  done
done < <(jq -r '.datasets[].dataset_id' "$MANIFEST")
```

Open representative dataset pages in the GSEA-enabled frontend and confirm the
GSEA table appears below DEA with Parquet and CSV.GZ download buttons. Check
that both DEA and GSEA download endpoints follow to nonempty files in the
development serving bucket. The API and frontend revisions must include the
dataset-level GSEA routes and display before this check can pass.

To regenerate only metadata JSONs, run `python3 release/metadata.py` with the
same project, dataset, location, and bucket options.

`OPENTARGETS_RELEASE` is optional and defaults to `26.03`. `ES_INDEX_SET` is
optional and defaults to empty. The projector creates a unique timestamped
index set on each run; the optional suffix is appended to index names and
aliases.

### 4. Grant IAM permissions to Cloud Build service account

The Cloud Build service account (`PROJECT_NUMBER@cloudbuild.gserviceaccount.com`) needs the following roles:

```bash
export CB_SA=$(gcloud projects describe $GCLOUD_PROJECT --format='value(projectNumber)')@cloudbuild.gserviceaccount.com

gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$CB_SA" \
    --role="roles/bigquery.dataEditor"

gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$CB_SA" \
    --role="roles/bigquery.jobUser"

gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$CB_SA" \
    --role="roles/storage.objectAdmin"

gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$CB_SA" \
    --role="roles/cloudsql.client"

gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$CB_SA" \
    --role="roles/run.admin"

gcloud iam service-accounts add-iam-policy-binding \
    "$(gcloud projects describe $GCLOUD_PROJECT --format='value(projectNumber)')-compute@developer.gserviceaccount.com" \
    --member="serviceAccount:$CB_SA" \
    --role="roles/iam.serviceAccountUser"

# The Cloud Run Job's default compute service account also needs BigQuery job
# and read access plus write access to CLOUD_TMP_BUCKET.
export RELEASE_SA="$(gcloud projects describe $GCLOUD_PROJECT --format='value(projectNumber)')-compute@developer.gserviceaccount.com"
gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$RELEASE_SA" \
    --role="roles/bigquery.jobUser"
gcloud projects add-iam-policy-binding $GCLOUD_PROJECT \
    --member="serviceAccount:$RELEASE_SA" \
    --role="roles/bigquery.dataViewer"
gcloud storage buckets add-iam-policy-binding "gs://$CLOUD_TMP_BUCKET" \
    --member="serviceAccount:$RELEASE_SA" \
    --role="roles/storage.objectAdmin"
```

### 5. Create Cloud Build private worker pool

The pipeline connects to Cloud SQL via its internal (VPC) IP. This requires a Cloud Build private worker pool connected to your VPC.

**One-time setup:**

```bash
gcloud builds worker-pools create dwh-pipeline-pool \
    --project=$GCLOUD_PROJECT \
    --region=$GCLOUD_REGION \
    --peered-network=projects/$GCLOUD_PROJECT/global/networks/default \
    --worker-machine-type=e2-highmem-4
```

## Running the pipeline

```bash
./dwh/trigger_pipeline.sh
```

### Suppress specific datasets

To exclude datasets from metadata tables (while keeping them in data tables):

```bash
./dwh/trigger_pipeline.sh --suppress-datasets "dataset_id_1,dataset_id_2"
```

### Force-refresh replaced Perturb-seq datasets

When selected DEA/GSEA rows were replaced in BigQuery, pass their exact IDs so
Postgres reloads those partitions even if `max_ingested_at` did not change:

```bash
./dwh/trigger_pipeline.sh --force-dataset-ids "dataset_id_1,dataset_id_2"
```

The preflight requires every ID to exist in `perturb_seq.metadata`; an ID may
have zero result rows in either result table. It replaces only those Perturb-seq
partitions, clearing stale rows when a result table is empty; other modalities
and unselected Perturb-seq data use the normal incremental sync.

To restrict release-artifact generation to a selected dataset set, also pass
`--release-dataset-ids`. For this publication run, use the same manifest IDs
for both options:

```bash
./dwh/trigger_pipeline.sh \
  --force-dataset-ids "$DATASET_IDS" \
  --release-dataset-ids "$DATASET_IDS"
```

The release filter applies to the staged metadata, DEA, and GSEA tables, so the
Cloud Run Job creates artifacts only for the selected IDs. Without this option,
release artifacts are generated for all available datasets.

### What happens

1. The script submits a Cloud Build job and starts streaming logs
2. You can **safely close your laptop** — the build continues in Google Cloud
3. To re-attach to logs later:
   ```bash
   gcloud builds log --stream BUILD_ID --region=$GCLOUD_REGION --project=$GCLOUD_PROJECT
   ```
4. To list recent builds:
   ```bash
   gcloud builds list --region=$GCLOUD_REGION --project=$GCLOUD_PROJECT --limit=5
   ```

## Running stages individually

For debugging or partial re-runs, you can run each stage manually.

### Open Targets reference

To reload the Open Targets reference table manually, you can execute a native serverless load command directly via the BigQuery CLI:

```bash
bq load --source_format=PARQUET \
    --location=$BQ_LOCATION \
    --replace \
    $GCLOUD_PROJECT:$BQ_DATASET.opentargets_targets \
    "gs://open-targets-data-releases/$OPENTARGETS_RELEASE/output/target/*.parquet"
```

### dbt

```bash
cd dwh
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cd bq_dbt
dbt run --profiles-dir .
```

To suppress datasets: `dbt run --profiles-dir . --vars '{"suppress_datasets": "id1,id2"}'`

Additional dbt commands:
- Specific model + dependencies: `dbt run --profiles-dir . --select +dataset_summary`
- Full refresh (non-incremental): `dbt run --profiles-dir . --full-refresh`

### BQ → Postgres

Cloud Build uses `PG_CONN_INTERNAL` to reach Cloud SQL over the private VPC.
For a manual run from an approved network path, use `PG_CONN` instead.

```bash
cd dwh
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 bq_to_postgres/bq_to_postgres.py \
    --bq-dataset $BQ_DATASET \
    --bq-location $BQ_LOCATION \
    --pg-conn "$PG_CONN" \
    --gcs-bucket "$GCLOUD_TMP_BUCKET" \
    --drop-and-recreate-indexes
```

For datasets whose BigQuery rows were explicitly replaced, pass their exact
comma-separated IDs with `--force-dataset-ids`. This drops and reloads only
those dataset partitions in `perturb_seq_dea` and `perturb_seq_gsea`, even when
their `max_ingested_at` timestamps are unchanged. Without this option the
incremental sync selects replacements by timestamp.

### BQ → Elasticsearch

```bash
cd dwh
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 bq_to_elastic/bq_to_es_projector.py --dataset-metadata ../be/dataset_metadata.json
```

The projector loads all three summary tables into new timestamped indexes,
then moves the
`dataset-summary`, `target-summary` and `landing-page-summary` aliases only
after all three loads succeed. It keeps the three newest indexes per family.
With a custom suffix such as `-ensg-dev`, the dated names and aliases include
that suffix.

## Mark datasets as reprocessed

Do not check or set `perturb_seq.reprocessed_datasets` before publication. Once
BigQuery replacement, Postgres/API synchronization, Elasticsearch metadata and
serving-bucket artifact checks have all succeeded, add the exact IDs using an
idempotent BigQuery merge:

```sql
MERGE `PROJECT_ID.perturb_seq.reprocessed_datasets` AS target
USING (
  SELECT dataset_id
  FROM UNNEST(["dataset_id_1", "dataset_id_2"]) AS dataset_id
) AS source
ON target.dataset_id = source.dataset_id
WHEN NOT MATCHED THEN
  INSERT (dataset_id) VALUES (source.dataset_id);
```

Replace `PROJECT_ID` and the example IDs with the development project and the
exact IDs from the successful run manifest. Then refresh the dataset mart and
search index so the reprocessed flag is visible in dataset discovery:

```bash
cd dwh/bq_dbt
dbt run --profiles-dir . --select dataset_summary
cd ..
python3 bq_to_elastic/bq_to_es_projector.py \
  --dataset-metadata ../be/dataset_metadata.json
cd ..
```

The published metadata JSONs also contain the dataset summary. Regenerate
metadata into a unique temporary prefix, then replace only the selected
Perturb-seq metadata objects in the versioned serving prefix:

```bash
set -euo pipefail
CLOUD_TMP_BUCKET="${CLOUD_TMP_BUCKET:-${GCLOUD_TMP_BUCKET:-}}"
MANIFEST=data_sources/perturb-seq/pipeline/manifests/single-condition-20.json
: "${RELEASE_VERSION_PREFIX:=2026.10}"
METADATA_PREFIX="release/marker-refresh-$(date +%s)"
python3 release/metadata.py --prefix "$METADATA_PREFIX"
while IFS= read -r dataset_id; do
  gcloud storage cp \
    "gs://$CLOUD_TMP_BUCKET/$METADATA_PREFIX/perturb-seq/$dataset_id.metadata.json" \
    "gs://perturbation-catalogue-release/$RELEASE_VERSION_PREFIX/perturb-seq/$dataset_id.metadata.json"
done < <(jq -r '.datasets[].dataset_id' "$MANIFEST")
gcloud storage rm --recursive "gs://$CLOUD_TMP_BUCKET/$METADATA_PREFIX/**"
```

Finally verify the dataset API responses, result-table row counts and the
`perturb_seq_reprocessed` value in the development dataset-summary index.

## Creating the PostgreSQL instance

Create the database at: https://console.cloud.google.com/sql/instances/create;engine=PostgreSQL;template=POSTGRES_ENTERPRISE_PLUS_DATA_CACHE_ENABLED_DEV_TEMPLATE

Settings:
* Cloud SQL edition: Enterprise
* Edition preset: Development
* Instance ID: `$PG_INSTANCE_ID`
* Password: `$PG_PASSWORD`
* Zonal availability: Single zone
* Customise your instance:
  - Machine configuration: Dedicated core, 4 vCPU, 16 GB
  - Data protection: disable "Automated daily backups" and "Enable point-in-time recovery"
  - Flags and parameters:
    - `temp_file_limit` = 104857600
    - `maintenance_work_mem` = 4194304
    - `max_parallel_maintenance_workers` = 4

### Materialized views

The script expects these materialized views to exist for `perturb_seq_dea`. Create them once after the initial data load:

```sql
CREATE MATERIALIZED VIEW perturb_seq_summary_perturbation AS
SELECT
    dataset_id,
    perturbed_target_ensg,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up
FROM perturb_seq_dea
WHERE padj <= 0.05 AND perturbed_target_ensg IS NOT NULL
GROUP BY dataset_id, perturbed_target_ensg;

CREATE UNIQUE INDEX idx_perturb_seq_summary_perturbation_pk ON perturb_seq_summary_perturbation (dataset_id, perturbed_target_ensg);

CREATE MATERIALIZED VIEW perturb_seq_summary_effect AS
SELECT
    dataset_id,
    effect_gene_ensg,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up,
    AVG(score_value) AS avg_score
FROM perturb_seq_dea
WHERE padj <= 0.05 AND effect_gene_ensg IS NOT NULL
GROUP BY dataset_id, effect_gene_ensg;

CREATE UNIQUE INDEX idx_perturb_seq_summary_effect_pk ON perturb_seq_summary_effect (dataset_id, effect_gene_ensg);

CREATE MATERIALIZED VIEW perturb_seq_summary_dataset AS
SELECT
    dataset_id,
    COUNT(*) AS n_total
FROM perturb_seq_dea
WHERE effect_gene_ensg IS NOT NULL
GROUP BY dataset_id;

CREATE UNIQUE INDEX idx_perturb_seq_summary_dataset_pk ON perturb_seq_summary_dataset (dataset_id);
```

## Dev → Prod migration

See [dev-to-prod.md](bq_to_postgres/dev-to-prod.md) for promoting a dev database to production.
