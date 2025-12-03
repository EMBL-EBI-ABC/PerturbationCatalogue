# BigQuery to Elasticsearch projector

## Environment configuration

These commands should be run from the parent directory, `dwh`.

```python
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cd bq_to_es_projector
```

## Running the projector

These commands should be run from the current directory, `bq_to_es_projector`.

First, load environment variables using `dev_secrets`.

Then, additionally set the following variables:
* BQ_TABLE: which BQ table to ingest from, for example `dataset_summary`.
* ES_INDEX: which Elastic index to ingest to, for example `dataset-summary-v3`.

Optional parameters that can also be set using environmental variables include:
```bash
BULK_CHUNK_SIZE=2000
BULK_MAX_RETRIES=5
BULK_TIMEOUT=120
```

Finally, run the projector: `python bq_to_es_projector.py`.
