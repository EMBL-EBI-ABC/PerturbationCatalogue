# dbt project for data marts

## Environment configuration
These commands should be run from the parent directory, `dwh`.

```python
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cd bq_dbt
```

## dbt commands

These commands should be run from this directory, `bq_dbt`.

### Initialize dbt
- `dbt init`
  1. Select `bigquery`
  2. Select `oauth`
  3. For project, enter the value of dev_secrets → $GCLOUD_PROJECT
  4. For dataset, enter the value of dev_secrets → $BQ_DATASET
  5. For threads, enter the number of CPUs available on your machine
  6. For job_execution_timeout_seconds, enter 3600
  7. For desired location option, choose EU
  8. Then, in ~/.dbt/profiles.yml, change `location: EU` to `location: europe-west2`

### Run dbt
- All models except listed: `dbt run --exclude target_summary gene_summary`
- Only listed models: `dbt run --select target_summary gene_summary`
- A specific model and its upstream dependencies: `dbt run --select +dataset_summary`
- A specific model and its upstream dependencies with full refresh: `dbt run --full-refresh --select +dataset_summary`. For example, this is needed if the schema has changed and all data needs to be reingested, not only the new rows.
