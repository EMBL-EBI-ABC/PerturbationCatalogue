# MaveDB

## Background

MaveDB data can be downloaded from Zenodo and is licensed under CC0.

The ZIP archive contains:

- `main.json` — metadata
- `csv/` — data packaged per study
  - `*.scores.csv` — always present
  - `*.counts.csv` — sometimes present

## Mirror

```bash
export LAKE_BUCKET=...
wget -O mavedb-data.zip "https://zenodo.org/records/11201737/files/mavedb-data.20240520.zip?download=1"
unzip -d mavedb mavedb-data.zip
gsutil -m cp -r mavedb gs://${LAKE_BUCKET}/
```

## Process

```bash
export WAREHOUSE_BUCKET=...
gsutil cp "gs://${LAKE_BUCKET}/mavedb/main.json" /tmp/main.json
python3 process.py \
  --input-filename /tmp/main.json \
  --output-filename /tmp/metadata.jsonl
gsutil cp /tmp/metadata.jsonl "gs://${WAREHOUSE_BUCKET}/mavedb/metadata.jsonl"
```

## Ingest into Elastic

```bash
gsutil cp "gs://${WAREHOUSE_BUCKET}/mavedb/metadata.jsonl" /tmp/metadata.json
# Create the index.
curl -X PUT "${ELASTIC_ENDPOINT}/mavedb" -H "Content-Type: application/json" -d'{
  "settings": {
    "index": {
      "number_of_shards": 1,
      "number_of_replicas": 1
    }
  }
}'
# Verify the index has been created – should contain "mavedb".
curl -X GET "${ELASTIC_ENDPOINT}/_cat/indices?v"
# Prepare data for loading.
awk '{print "{ \"index\" : { \"_index\" : \"mavedb\", \"_id\" : \"" NR "\" } }"; print}' /tmp/metadata.jsonl > /tmp/formatted_metadata.jsonl
# Load data into index.
curl -X POST "${ELASTIC_ENDPOINT}/mavedb/_bulk" -H "Content-Type: application/x-ndjson" --data-binary @/tmp/formatted_metadata.jsonl
# Verify data load.
curl -X GET "${ELASTIC_ENDPOINT}/mavedb/_search?filter_path=hits.hits._source" -H "Content-Type: application/json" -d'{
  "size": 10,
  "query": {
    "match_all": {}
  }
}'
```
