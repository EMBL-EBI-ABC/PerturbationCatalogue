#!/bin/bash
set -euo pipefail

export ELASTIC_ENDPOINT="$1"
export ELASTIC_INDEX="$2"
export JSONL_DATA="$3"

echo "Removing the index..."
curl -X DELETE "${ELASTIC_ENDPOINT}/${ELASTIC_INDEX}" -H "Content-Type: application/json"

echo "Creating the index..."
curl -X PUT "${ELASTIC_ENDPOINT}/${ELASTIC_INDEX}" -H "Content-Type: application/json" -d'{
  "settings": {
    "index": {
      "number_of_shards": 1,
      "number_of_replicas": 1
    }
  }
}'

echo "Preparing the data for loading..."
gsutil -q cp "${JSONL_DATA}" /tmp/metadata.jsonl
awk \
  -v elastic_index=${ELASTIC_INDEX} \
  '{print "{ \"index\" : { \"_index\" : \"" elastic_index "\", \"_id\" : \"" NR "\" } }"; print}' \
  /tmp/metadata.jsonl \
  > /tmp/formatted_metadata.jsonl

echo "Loading data into index..."
curl -X POST "${ELASTIC_ENDPOINT}/${ELASTIC_INDEX}/_bulk" -H "Content-Type: application/x-ndjson" --data-binary @/tmp/formatted_metadata.jsonl

echo "All done."
