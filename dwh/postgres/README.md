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
gcloud compute instances create bq-to-pg-projector \
    --project=${GCLOUD_PROJECT} \
    --zone=${GCLOUD_ZONE} \
    --machine-type=e2-medium \
    --network=default \
    --scopes=https://www.googleapis.com/auth/cloud-platform
```

## 2. Copy files and SSH into the VM
```bash
gcloud compute scp --project=${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE} postgres/requirements.txt postgres/bq_to_postgres.py bq-to-pg-projector:~
gcloud compute ssh bq-to-pg-projector --project=${GCLOUD_PROJECT} --zone=${GCLOUD_ZONE}
```

## 3. Create permanent session
Run `screen` so that the commands above can run for a long time and will not be interrupted by a connection failure.

## 4. Install dependencies
```bash
sudo apt install -y python3-pip python3-venv
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

## 5. Run the script
```bash
export BQ_TABLE=...
export PG_TABLE=...
python3 bq_to_postgres.py \
    --bq-dataset ${BQ_DATASET} \
    --bq-table ${BQ_TABLE} \
    --bq-location ${BQ_LOCATION} \
    --pg-conn "${PG_CONN}" \
    --pg-table ${PG_TABLE} \
    --gcs-bucket ${GCLOUD_TMP_BUCKET}
```

# Create indexes
After running the command above, run `psql $PG_CONN` and create the indexes.

## Perturb-Seq
```sql
\timing on
CREATE INDEX CONCURRENTLY idx_perturbation
  ON perturb_seq (perturbed_target_symbol, dataset_id, padj, basemean, log2foldchange);

CREATE INDEX CONCURRENTLY idx_phenotype
  ON perturb_seq (gene, dataset_id, padj, basemean, log2foldchange);
```

## CRISPR
```sql
\timing on
CREATE INDEX idx_crispr_data_dataset ON public.crispr_data USING btree (dataset_id);
CREATE INDEX idx_crispr_data_target ON public.crispr_data USING btree (perturbed_target_symbol);
```

## Monitoring
You can use this query in a separate psql session to monitor the progress of index creation:
```sql
SELECT
    p.pid,
    p.datname,
    c.relname AS table_name,
    i.relname AS index_name,
    p.phase,
    p.blocks_total,
    p.blocks_done,
    ROUND(100.0 * p.blocks_done / NULLIF(p.blocks_total, 0), 1) AS pct_done
FROM pg_stat_progress_create_index p
JOIN pg_class c ON p.relid = c.oid
LEFT JOIN pg_class i ON p.index_relid = i.oid;
\watch 10
```
