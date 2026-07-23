# Data Warehouse Pipeline

Automated pipeline for transforming and loading data from BigQuery to Postgres and Elasticsearch.

> **Before running any commands in this document**, set up the environment by running `dev_secrets`.

## Pipeline stages

The pipeline runs four stages sequentially:

| Stage | Directory | Description | Duration |
|-------|-----------|-------------|----------|
| 1. **Open Targets reference** | (Native) | Loads Open Targets Platform targets into the configured BQ reference table | ~minutes |
| 2. **dbt** | `bq_dbt/` | Transforms source BQ tables into final data mart tables | ~minutes |
| 3. **BQ → Postgres** | `bq_to_postgres/` | Loads final BQ data tables into Cloud SQL (Postgres) | ~hours |
| 4. **BQ → Elastic** | `bq_to_elastic/` | Loads summary tables into Elasticsearch | ~minutes |

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

The trigger script requires the following variables (all provided by `dev_secrets`): `GCLOUD_PROJECT`, `GCLOUD_REGION`, `BQ_DATASET`, `BQ_LOCATION`, `GCLOUD_TMP_BUCKET`, `PG_CONN_INTERNAL`, `ES_URL`, `ES_USERNAME`, `ES_PASSWORD`.

`OPENTARGETS_RELEASE` is optional and defaults to `26.03`. `ES_INDEX_SET` is
optional and defaults to empty. Its value is appended directly to all three
Elasticsearch index names.

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

### BQ → Elasticsearch

```bash
cd dwh
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 bq_to_elastic/bq_to_es_projector.py --dataset-metadata ../be/dataset_metadata.json
```

With an empty `ES_INDEX_SET`, projections use dated indexes and update the
`dataset-summary`, `target-summary`, and `landing-page-summary` aliases, keeping
the three newest dated indexes in each family. With a custom suffix, for example
`-ensg-dev`, projections overwrite the standalone `dataset-summary-ensg-dev`,
`target-summary-ensg-dev`, and `landing-page-summary-ensg-dev` indexes without
updating aliases or pruning old indexes.

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
