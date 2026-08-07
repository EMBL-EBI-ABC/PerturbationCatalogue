"""Create clustered release staging tables and a dataset task manifest."""

import argparse
import json
import logging
import re
from datetime import datetime, timedelta, timezone

from google.cloud import bigquery, storage

from release import DATASETS, data_query, table


def _table_id(run_id, modality, kind):
    safe_id = re.sub(r"[^a-z0-9_]", "_", run_id.lower())
    return f"release_{safe_id}_{modality.replace('-', '_')}_{kind}"


def _stage(client, project, dataset, location, name, query):
    relation = table(project, dataset, name)
    client.query(
        f"CREATE OR REPLACE TABLE {relation} CLUSTER BY dataset_id AS {query}",
        location=location,
    ).result()
    staged = client.get_table(f"{project}.{dataset}.{name}")
    staged.expires = datetime.now(timezone.utc) + timedelta(days=2)
    client.update_table(staged, ["expires"])
    return f"{project}.{dataset}.{name}"


def create_staging(project, dataset, location, bucket_name, run_id):
    client = bigquery.Client(project=project)
    items = []
    for modality, config in DATASETS.items():
        data_name = _table_id(run_id, modality, "data")
        metadata_name = _table_id(run_id, modality, "metadata")
        data_table = _stage(
            client,
            project,
            dataset,
            location,
            data_name,
            data_query(project, dataset, modality),
        )
        metadata_table = _stage(
            client,
            project,
            dataset,
            location,
            metadata_name,
            f"SELECT * FROM {table(project, dataset, config['metadata'])}",
        )
        rows = client.query(
            "SELECT dataset_id FROM "
            f"{table(project, dataset, data_name)} WHERE dataset_id IS NOT NULL "
            "UNION DISTINCT SELECT dataset_id FROM "
            f"{table(project, dataset, metadata_name)} WHERE dataset_id IS NOT NULL",
            location=location,
        ).result()
        items.extend(
            {
                "modality": modality,
                "dataset_id": row.dataset_id,
                "data_table": data_table,
                "metadata_table": metadata_table,
            }
            for row in rows
        )

    manifest = {"run_id": run_id, "items": items}
    path = f"release-staging/{run_id}/manifest.json"
    storage.Client(project=project).bucket(bucket_name).blob(path).upload_from_string(
        json.dumps(manifest), content_type="application/json"
    )
    return len(items)


def cleanup(project, dataset, location, bucket_name, run_id):
    client = bigquery.Client(project=project)
    for modality in DATASETS:
        for kind in ("data", "metadata"):
            client.delete_table(
                f"{project}.{dataset}.{_table_id(run_id, modality, kind)}",
                not_found_ok=True,
            )
    bucket = storage.Client(project=project).bucket(bucket_name)
    for blob in bucket.list_blobs(prefix=f"release-staging/{run_id}/"):
        blob.delete()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("create", "cleanup"))
    parser.add_argument("--project", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--location", required=True)
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if args.mode == "create":
        print(
            create_staging(
                args.project,
                args.dataset,
                args.location,
                args.bucket,
                args.run_id,
            )
        )
    else:
        cleanup(args.project, args.dataset, args.location, args.bucket, args.run_id)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )
    main()
