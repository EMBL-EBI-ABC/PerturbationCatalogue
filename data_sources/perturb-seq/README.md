# Perturb-Seq

## Background

Perturb-Seq data is processed from a CSV file containing differential expression results.

Before processing, set the following environment secrets:

- WAREHOUSE_BUCKET
- ELASTIC_ENDPOINT

## Ingest metadata

```bash
# Process.
python3 process.py \
  --input-filename /tmp/perturb-seq.csv \
  --output-filename /tmp/perturb-seq.jsonl
gsutil -q -m rm -r "gs://${WAREHOUSE_BUCKET}/perturb-seq"
gsutil -q cp /tmp/perturb-seq.jsonl "gs://${WAREHOUSE_BUCKET}/perturb-seq/metadata.jsonl"

# Ingest into Elastic.
python3 ../elastic_load.py \
  --elastic-endpoint "${ELASTIC_ENDPOINT}" \
  --elastic-index "perturb-seq" \
  --jsonl-data "gs://${WAREHOUSE_BUCKET}/perturb-seq/metadata.jsonl" \
  --id-field "record_id" \
  --field-properties '{
    "study_id": {
      "type": "keyword"
    },
    "perturbation": {
      "type": "keyword"
    },
    "gene": {
      "type": "keyword"
    }
  }'
```

## Stats for the four files

These are the four studies curated and currently used:

```
Processing complete:
  Records written: 29476
  Records skipped: 2803366
  Filter criteria:
    padj <= 0.05
    log2fc < -1.0 or > 1.0
```
