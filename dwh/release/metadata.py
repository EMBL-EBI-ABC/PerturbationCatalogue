"""Write the dataset-level metadata artifacts from the DWH summary table."""

import argparse
import json
import os

from google.cloud import bigquery, storage

from release import DATASET_METADATA, table

MODALITIES = {
    "CRISPR screen": "crispr",
    "MAVE": "mave",
    "Perturb-seq": "perturb-seq",
}


def generate(project, dataset, location, bucket_name, prefix="release"):
    rows = (
        bigquery.Client(project=project)
        .query(
            f"SELECT * FROM {table(project, dataset, DATASET_METADATA)}",
            location=location,
        )
        .result()
    )
    bucket = storage.Client(project=project).bucket(bucket_name)
    for row in rows:
        if row.dataset_id:
            modality = MODALITIES[row.data_modalities[0]]
            bucket.blob(
                f"{prefix}/{modality}/" f"{row.dataset_id}.metadata.json"
            ).upload_from_string(
                json.dumps(dict(row), default=str, separators=(",", ":")),
                content_type="application/json",
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default=os.getenv("GCLOUD_PROJECT"))
    parser.add_argument("--dataset", default=os.getenv("BQ_DATASET"))
    parser.add_argument("--location", default=os.getenv("BQ_LOCATION"))
    parser.add_argument(
        "--bucket",
        default=os.getenv("CLOUD_TMP_BUCKET") or os.getenv("GCLOUD_TMP_BUCKET"),
    )
    parser.add_argument("--prefix", default="release")
    args = parser.parse_args()
    missing = [
        name
        for name in ("project", "dataset", "location", "bucket")
        if not getattr(args, name)
    ]
    if missing:
        parser.error("missing: " + ", ".join(missing))
    generate(args.project, args.dataset, args.location, args.bucket, args.prefix)


if __name__ == "__main__":
    main()
