# MaveDB

## Background

MaveDB data can be downloaded from Zenodo and is licensed under CC0.

The ZIP archive contains:

- `main.json` — metadata
- `csv/` — data packaged per study
  - `*.scores.csv` — always present
  - `*.counts.csv` — sometimes present

Before processing, set the following environment secrets:
* LAKE_BUCKET
* WAREHOUSE_BUCKET
* ELASTIC_ENDPOINT

## Ingest metadata

```bash
# Mirror.
wget -q -O mavedb-data.zip "https://zenodo.org/records/14172004/files/mavedb-dump.20241114101443.zip?download=1"
unzip -q -d mavedb mavedb-data.zip
gsutil -q -m rm -r "gs://${LAKE_BUCKET}/mavedb"
gsutil -q -m cp -r mavedb "gs://${LAKE_BUCKET}/"

# Process.
gsutil -q cp "gs://${LAKE_BUCKET}/mavedb/main.json" /tmp/main.json
python3 process.py \
  --input-filename /tmp/main.json \
  --output-filename /tmp/metadata.jsonl
gsutil -q -m rm -r "gs://${WAREHOUSE_BUCKET}/mavedb"
gsutil -q cp /tmp/metadata.jsonl "gs://${WAREHOUSE_BUCKET}/mavedb/metadata.jsonl"

# Ingest into Elastic.
## 1. Remove the index.
curl -X DELETE "${ELASTIC_ENDPOINT}/mavedb" -H "Content-Type: application/json"
## 2. Create the index.
curl -X PUT "${ELASTIC_ENDPOINT}/mavedb" -H "Content-Type: application/json" -d'{
  "settings": {
    "index": {
      "number_of_shards": 1,
      "number_of_replicas": 1
    }
  }
}'
## 3. Prepare data for loading.
gsutil -q cp "gs://${WAREHOUSE_BUCKET}/mavedb/metadata.jsonl" /tmp/metadata.json
awk '{print "{ \"index\" : { \"_index\" : \"mavedb\", \"_id\" : \"" NR "\" } }"; print}' /tmp/metadata.jsonl > /tmp/formatted_metadata.jsonl
## 4. Load data into index.
curl -X POST "${ELASTIC_ENDPOINT}/mavedb/_bulk" -H "Content-Type: application/x-ndjson" --data-binary @/tmp/formatted_metadata.jsonl
```
