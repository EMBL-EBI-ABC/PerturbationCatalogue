# BigQuery to Elasticsearch projector

## Running the projector

These commands should be run from the parent directory, `dwh`.

The script is automated and will:
1. Load the three tables (`dataset_summary`, `target_summary`, and `landing_page_summary`) into Elastic under the format `YYYY-MM-DD-index-name`.
2. If the sync is successful, move aliases such as `dataset-summary` to point to the latest index version.
3. If the sync is successful, prune old index versions to keep only the live one + up to two earlier versions. 

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cd bq_to_es_projector
dev_secrets
python3 bq_to_es_projector.py
```

## Additional configuration

Optional parameters that can be set using environmental variables:
```bash
BULK_CHUNK_SIZE=2000
BULK_MAX_RETRIES=5
BULK_TIMEOUT=120
```
