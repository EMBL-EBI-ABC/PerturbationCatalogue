"""Create clustered release staging tables and a dataset task manifest."""

import argparse
import json
import logging
import re
from datetime import datetime, timedelta, timezone
from typing import Optional

from google.cloud import bigquery, storage

from release import DATASET_METADATA, DATASETS, data_query, gsea_query, table


def _table_id(run_id, modality, kind):
    safe_id = re.sub(r"[^a-z0-9_]", "_", run_id.lower())
    return f"release_{safe_id}_{modality.replace('-', '_')}_{kind}"


def _validate_dataset_ids(dataset_ids):
    if dataset_ids is None:
        return None
    dataset_ids = [dataset_id.strip() for dataset_id in dataset_ids]
    if not dataset_ids or any(not dataset_id for dataset_id in dataset_ids):
        raise ValueError("dataset IDs must be a nonempty comma-separated list")
    if len(set(dataset_ids)) != len(dataset_ids):
        raise ValueError("dataset IDs must be unique")
    return dataset_ids


def _parse_dataset_ids(value):
    return _validate_dataset_ids(value.split(","))


def _stage(client, project, dataset, location, name, query, dataset_ids=None):
    relation = table(project, dataset, name)
    query_kwargs = {"location": location}
    if dataset_ids is not None:
        query = (
            f"SELECT * FROM ({query}) AS release_rows "
            "WHERE dataset_id IN UNNEST(@dataset_ids)"
        )
        query_kwargs["job_config"] = bigquery.QueryJobConfig(
            query_parameters=[
                bigquery.ArrayQueryParameter("dataset_ids", "STRING", dataset_ids)
            ]
        )
    client.query(
        f"CREATE OR REPLACE TABLE {relation} CLUSTER BY dataset_id AS {query}",
        **query_kwargs,
    ).result()
    staged = client.get_table(f"{project}.{dataset}.{name}")
    staged.expires = datetime.now(timezone.utc) + timedelta(days=2)
    client.update_table(staged, ["expires"])
    return f"{project}.{dataset}.{name}"


def create_staging(
    project, dataset, location, bucket_name, run_id, dataset_ids: Optional[list] = None
):
    dataset_ids = _validate_dataset_ids(dataset_ids)
    client = bigquery.Client(project=project)
    items = []
    metadata_name = _table_id(run_id, "dataset", "metadata")
    metadata_table = _stage(
        client,
        project,
        dataset,
        location,
        metadata_name,
        f"SELECT * FROM {table(project, dataset, DATASET_METADATA)}",
        dataset_ids,
    )
    # ponytail: cloudbuild's cleanup list only covers per-modality data tables;
    # this extra stage expires after two days until that list is generalized.
    gsea_table = _stage(
        client,
        project,
        dataset,
        location,
        _table_id(run_id, "perturb-seq", "gsea"),
        gsea_query(project, dataset),
        dataset_ids,
    )
    for modality in DATASETS:
        data_name = _table_id(run_id, modality, "data")
        data_table = _stage(
            client,
            project,
            dataset,
            location,
            data_name,
            data_query(project, dataset, modality),
            dataset_ids,
        )
        id_sources = [
            f"SELECT DISTINCT dataset_id FROM {table(project, dataset, data_name)} "
            "WHERE dataset_id IS NOT NULL",
            f"SELECT dataset_id FROM {table(project, dataset, metadata_name)} "
            "WHERE dataset_id IS NOT NULL",
        ]
        if modality == "perturb-seq":
            id_sources.append(
                f"SELECT dataset_id FROM {gsea_table} WHERE dataset_id IS NOT NULL"
            )
        dataset_id_query = " UNION DISTINCT ".join(id_sources)
        rows = client.query(
            dataset_id_query,
            location=location,
        ).result()
        for row in rows:
            item = {
                "modality": modality,
                "dataset_id": row.dataset_id,
                "data_table": data_table,
                "metadata_table": metadata_table,
            }
            if modality == "perturb-seq":
                item["gsea_table"] = gsea_table
            items.append(item)

    manifest = {"run_id": run_id, "items": items}
    path = f"release-staging/{run_id}/manifest.json"
    storage.Client(project=project).bucket(bucket_name).blob(path).upload_from_string(
        json.dumps(manifest), content_type="application/json"
    )
    return len(items)


def cleanup(project, dataset, location, bucket_name, run_id):
    client = bigquery.Client(project=project)
    for modality in DATASETS:
        client.delete_table(
            f"{project}.{dataset}.{_table_id(run_id, modality, 'data')}",
            not_found_ok=True,
        )
    client.delete_table(
        f"{project}.{dataset}.{_table_id(run_id, 'perturb-seq', 'gsea')}",
        not_found_ok=True,
    )
    client.delete_table(
        f"{project}.{dataset}.{_table_id(run_id, 'dataset', 'metadata')}",
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
    parser.add_argument("--dataset-ids", type=_parse_dataset_ids)
    args = parser.parse_args()
    if args.mode == "create":
        print(
            create_staging(
                args.project,
                args.dataset,
                args.location,
                args.bucket,
                args.run_id,
                args.dataset_ids,
            )
        )
    else:
        cleanup(args.project, args.dataset, args.location, args.bucket, args.run_id)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )
    main()
