# Data Warehouse Pipeline

Automated pipeline for transforming and loading data from BigQuery to Postgres and Elasticsearch.

## Pipeline stages

The pipeline runs three stages sequentially:

| Stage | Directory | Description | Duration |
|-------|-----------|-------------|----------|
| 1. **dbt** | `bq_dbt/` | Transforms source BQ tables into final data mart tables | ~minutes |
| 2. **BQ → Postgres** | `bq_to_postgres/` | Loads final BQ data tables into Cloud SQL (Postgres) | ~hours |
| 3. **BQ → Elastic** | `bq_to_elastic/` | Loads summary tables into Elasticsearch | ~minutes |

Each stage depends on the previous one. If any stage fails, the pipeline stops.

## Prerequisites

1. **Google Cloud SDK** (`gcloud`) installed and authenticated
2. **Cloud Build API** enabled in your GCP project
3. **Environment variables** — source your secrets before running:
   ```bash
   dev_secrets
   ```
   Required variables: `GCLOUD_PROJECT`, `BQ_DATASET`, `BQ_LOCATION`, `GCLOUD_TMP_BUCKET`, `PG_CONN_INTERNAL`, `ES_URL`, `ES_USERNAME`, `ES_PASSWORD`

4. **Cloud Build permissions** — the Cloud Build service account needs:
   - BigQuery Data Editor & Job User
   - Cloud Storage Object Admin (for temp GCS files)
   - Cloud SQL Client (for Postgres access)

5. **Network access** — if your Postgres instance uses a private IP, you need either:
   - A [Cloud Build private pool](https://cloud.google.com/build/docs/private-pools/create-manage-private-pools) with VPC access, or
   - The Cloud SQL instance's public IP enabled

## Running the pipeline

```bash
dev_secrets
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
   gcloud builds log --stream BUILD_ID --project=$GCLOUD_PROJECT
   ```
4. To list recent builds:
   ```bash
   gcloud builds list --project=$GCLOUD_PROJECT --limit=5
   ```

## Running stages individually

For debugging or partial re-runs, you can run each stage manually.

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
python3 bq_to_elastic/bq_to_es_projector.py
```

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
    perturbed_target_symbol,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up
FROM perturb_seq_dea
WHERE padj <= 0.05
GROUP BY dataset_id, perturbed_target_symbol;

CREATE UNIQUE INDEX idx_perturb_seq_summary_perturbation_pk ON perturb_seq_summary_perturbation (dataset_id, perturbed_target_symbol);

CREATE MATERIALIZED VIEW perturb_seq_summary_effect AS
SELECT
    dataset_id,
    gene,
    COUNT(*) AS n_total,
    COUNT(*) FILTER (WHERE log2foldchange < 0) AS n_down,
    COUNT(*) FILTER (WHERE log2foldchange > 0) AS n_up,
    AVG(score_value) AS avg_score
FROM perturb_seq_dea
WHERE padj <= 0.05
GROUP BY dataset_id, gene;

CREATE UNIQUE INDEX idx_perturb_seq_summary_effect_pk ON perturb_seq_summary_effect (dataset_id, gene);

CREATE MATERIALIZED VIEW perturb_seq_summary_dataset AS
SELECT
    dataset_id,
    COUNT(*) AS n_total
FROM perturb_seq_dea
WHERE gene IS NOT NULL
GROUP BY dataset_id;

CREATE UNIQUE INDEX idx_perturb_seq_summary_dataset_pk ON perturb_seq_summary_dataset (dataset_id);
```

## Dev → Prod migration

See [dev-to-prod.md](bq_to_postgres/dev-to-prod.md) for promoting a dev database to production.
