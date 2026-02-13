# Create PostgreSQL instance

Create the database using the link: https://console.cloud.google.com/sql/instances/create;engine=PostgreSQL;template=POSTGRES_ENTERPRISE_PLUS_DATA_CACHE_ENABLED_DEV_TEMPLATE

Change the settings:
* Cloud SQL edition: Enterprise
* Edition preset: Development
* Instance ID: `$PG_INSTANCE_ID`
* Password: `$PG_PASSWORD`
* Zonal availability: Single zone
* Customise your instance
  + Machine configuration: Dedicated core, 4 vCPU, 16 GB
  + Data protection: disable "Automated daily backups" and "Enable point-in-time recovery"
  + Flags and parameters:
    - `temp_file_limit` = 104857600
    - `maintenance_work_mem` = 4194304
    - `max_parallel_maintenance_workers` = 4

# Migrate data from BigQuery

## 1. Create a Google Cloud VM
```bash
dev_secrets
gcloud compute instances create bq-to-pg-projector \
    --project=${GCLOUD_PROJECT} \
    --zone=${GCLOUD_ZONE} \
    --machine-type=e2-highmem-4 \
    --network=default \
    --scopes=https://www.googleapis.com/auth/cloud-platform
```

## 2. Copy files and SSH into the VM
```bash
gcloud compute scp --project=${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE} dwh/postgres/requirements.txt dwh/postgres/bq_to_postgres.py bq-to-pg-projector:~
gcloud compute ssh bq-to-pg-projector --project=${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE}
```

## 3. Create permanent session
Run `screen` so that the commands above can run for a long time and will not be interrupted by a connection failure. (If you are getting a terminfo issue, which can happen with newer terminals such as `foot`, run `export TERM=xterm-256color` before running `screen`.)

## 4. Install dependencies
```bash
sudo apt update
sudo apt install -y python3-pip python3-venv postgresql-client
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

## 5. Run the script
Note: you should set `$PG_CONN` to `$PG_CONN_INTERNAL` from the list of secrets, as the VM is connected to the VPC and should connect to the SQL instance via its private IP.

Standard mode (keeps indexes, updates data):
```bash
python3 bq_to_postgres.py \
    --bq-dataset ${BQ_DATASET} \
    --bq-location ${BQ_LOCATION} \
    --pg-conn "${PG_CONN}" \
    --gcs-bucket "${GCLOUD_TMP_BUCKET}"
```

Mode with index dropping (drops indexes -> updates data -> recreates indexes):
Use this for large updates where updating indexes row-by-row is too slow.
```bash
python3 bq_to_postgres.py \
    --bq-dataset ${BQ_DATASET} \
    --bq-location ${BQ_LOCATION} \
    --pg-conn "${PG_CONN}" \
    --gcs-bucket "${GCLOUD_TMP_BUCKET}" \
    --drop-and-recreate-indexes
```

## 6. Remove the VM
Once the ingestion is complete (including any index creation as described above), exit the session and remove the instance:
```bash
gcloud compute instances delete bq-to-pg-projector --project ${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE}
```

# Materialized Views

The script expects the following materialized views to exist for `perturb_seq_dea`. It will refresh them concurrently after data sync.

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

# Migrate tables and indexes from development to production

Export and import is done via a shared temporary bucket. 

## Export in development environment
```bash
dev_secrets
export PG_TABLE=...
# Start export
gcloud sql export sql --async ${PG_INSTANCE_ID} \
  gs://${GCLOUD_TMP_BUCKET}/postgres_dev_to_prod/${PG_TABLE}.sql \
  --project=${GCLOUD_PROJECT} \
  --database=${PG_DB} \
  --table=${PG_TABLE}
# Wait for the operation to complete
gcloud beta sql operations wait --timeout 86400 $(gcloud beta sql operations list --instance=${PG_INSTANCE_ID} --filter STATUS=RUNNING --format "value(name)")
```

## Import in production environment
```bash
prod_secrets
export PG_TABLE=...
# Start import
gcloud sql import sql --async ${PG_INSTANCE_ID} \
  gs://${GCLOUD_TMP_BUCKET}/postgres_dev_to_prod/${PG_TABLE}.sql \
  --database=${PG_DB}
# Wait for the operation to complete
gcloud beta sql operations wait --timeout 86400 $(gcloud beta sql operations list --instance=${PG_INSTANCE_ID} --filter STATUS=RUNNING --format "value(name)")
```

## Clean up
```bash
gsutil rm gs://${GCLOUD_TMP_BUCKET}/postgres_dev_to_prod/${PG_TABLE}.sql
```

## Note on service accounts
Cloud SQL service accounts of *both* projects need to have access to the bucket. Run the below command in both projects. It only needs to be run once per project:

```bash
gsutil iam ch serviceAccount:$(gcloud sql instances describe ${PG_INSTANCE_ID} --project=${GCLOUD_PROJECT} --format="value(serviceAccountEmailAddress)"):roles/storage.objectAdmin gs://${GCLOUD_TMP_BUCKET}
```
