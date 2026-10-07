#!/usr/bin/env python3
"""Replace selected Perturb-seq DEA/GSEA datasets in BigQuery.

Manifest format:
{
  "datasets": [
    {
      "dataset_id": "example_2026",
      "dea_parquet": "example_2026.dea.parquet",
      "gsea_parquet": "example_2026.gsea.parquet"
    }
  ]
}

The default run validates local Parquet files only. Pass --apply and --run-id to
stage and replace rows. For apply, --project (or GCLOUD_PROJECT), BQ_LOCATION, and
CLOUD_TMP_BUCKET (or GCLOUD_TMP_BUCKET) must point at the dev environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from io_schemas import DEA_SCHEMA, GSEA_SCHEMA


TARGETS = {
    "dea": ("pertpy_dea", DEA_SCHEMA),
    "gsea": ("pertpy_gsea", GSEA_SCHEMA),
}
DATASET_ID_PATTERN = re.compile(r"[A-Za-z0-9_-]+\Z")
RUN_ID_PATTERN = re.compile(r"[a-z0-9_]{1,40}\Z")
PROJECT_PATTERN = re.compile(r"[A-Za-z0-9-]+\Z")
STAGING_TTL = timedelta(days=2)


@dataclass(frozen=True)
class DatasetFiles:
    dataset_id: str
    dea_path: Path
    gsea_path: Path
    dea_rows: int
    gsea_rows: int


def _parquet_rows(path: Path, dataset_id: str, expected_schema: pa.Schema) -> int:
    parquet = pq.ParquetFile(path)
    schema = parquet.schema_arrow
    if not schema.equals(expected_schema, check_metadata=False):
        raise ValueError(
            f"{path}: schema does not match the pipeline schema for "
            f"{path.name.split('.')[-2].upper()}"
        )

    row_count = parquet.metadata.num_rows
    for batch in parquet.iter_batches(columns=["dataset_id"], batch_size=1_000_000):
        values = batch.column("dataset_id").to_pylist()
        if any(value is None or value != dataset_id for value in values):
            raise ValueError(f"{path}: dataset_id values must all equal {dataset_id!r}")
    return row_count


def validate_manifest(manifest_path: str | Path) -> list[DatasetFiles]:
    manifest_path = Path(manifest_path).resolve()
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read manifest {manifest_path}: {error}") from error

    datasets = manifest.get("datasets") if isinstance(manifest, dict) else None
    if not isinstance(datasets, list) or not datasets:
        raise ValueError("Manifest must contain a nonempty 'datasets' list")

    seen: set[str] = set()
    output: list[DatasetFiles] = []
    for item in datasets:
        if not isinstance(item, dict):
            raise ValueError("Each manifest entry must be an object")
        dataset_id = item.get("dataset_id")
        if not isinstance(dataset_id, str) or not DATASET_ID_PATTERN.fullmatch(
            dataset_id
        ):
            raise ValueError(f"Invalid dataset_id in manifest: {dataset_id!r}")
        if dataset_id in seen:
            raise ValueError(f"Duplicate dataset_id in manifest: {dataset_id!r}")
        seen.add(dataset_id)

        paths: dict[str, Path] = {}
        for key, suffix in (("dea_parquet", "dea"), ("gsea_parquet", "gsea")):
            value = item.get(key)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{dataset_id}: manifest requires {key}")
            path = Path(value)
            if not path.is_absolute():
                path = manifest_path.parent / path
            path = path.resolve()
            expected_name = f"{dataset_id}.{suffix}.parquet"
            if path.name != expected_name:
                raise ValueError(
                    f"{dataset_id}: expected file name {expected_name!r}, "
                    f"got {path.name!r}"
                )
            if not path.is_file():
                raise ValueError(f"{dataset_id}: Parquet file does not exist: {path}")
            paths[suffix] = path

        dea_rows = _parquet_rows(paths["dea"], dataset_id, DEA_SCHEMA)
        gsea_rows = _parquet_rows(paths["gsea"], dataset_id, GSEA_SCHEMA)
        output.append(
            DatasetFiles(
                dataset_id=dataset_id,
                dea_path=paths["dea"],
                gsea_path=paths["gsea"],
                dea_rows=dea_rows,
                gsea_rows=gsea_rows,
            )
        )
    return output


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _run_key(
    run_id: str, files: list[DatasetFiles]
) -> tuple[str, dict[tuple[str, str], str]]:
    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise ValueError(
            "run-id must contain 1-40 lowercase letters, digits, or underscores"
        )
    hashes = {
        (item.dataset_id, kind): _sha256(path)
        for item in files
        for kind, path in (("dea", item.dea_path), ("gsea", item.gsea_path))
    }
    identity = {
        "run_id": run_id,
        "inputs": [
            {
                "dataset_id": item.dataset_id,
                "dea_sha256": hashes[(item.dataset_id, "dea")],
                "gsea_sha256": hashes[(item.dataset_id, "gsea")],
            }
            for item in files
        ],
    }
    key = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:32]
    return key, hashes


def _schema_signature(fields: list[Any]) -> tuple[tuple[str, str, str], ...]:
    aliases = {"FLOAT": "FLOAT64", "INTEGER": "INT64", "BOOLEAN": "BOOL"}
    return tuple(
        sorted(
            (
                field.name,
                aliases.get(field.field_type.upper(), field.field_type.upper()),
                (field.mode or "NULLABLE").upper(),
            )
            for field in fields
        )
    )


def _expected_bq_schema(schema: pa.Schema, bigquery: Any) -> list[Any]:
    fields = []
    for field in schema:
        value_type = field.type
        mode = "REQUIRED" if not field.nullable else "NULLABLE"
        if pa.types.is_list(value_type) or pa.types.is_large_list(value_type):
            value_type = value_type.value_type
            mode = "REPEATED"
        if pa.types.is_string(value_type):
            bq_type = "STRING"
        elif pa.types.is_float64(value_type):
            bq_type = "FLOAT64"
        elif pa.types.is_int64(value_type):
            bq_type = "INT64"
        elif pa.types.is_timestamp(value_type):
            bq_type = "TIMESTAMP"
        else:
            raise ValueError(
                f"Unsupported Parquet field type: {field.name} {field.type}"
            )
        fields.append(bigquery.SchemaField(field.name, bq_type, mode=mode))
    return fields


def _check_table(client: Any, project: str, location: str, bigquery: Any) -> None:
    dataset = client.get_dataset(f"{project}.perturb_seq")
    if (dataset.location or "").casefold() != location.casefold():
        raise ValueError(
            f"BQ dataset location {dataset.location!r} does not match BQ_LOCATION {location!r}"
        )
    for kind, (table_name, schema) in TARGETS.items():
        table = client.get_table(f"{project}.perturb_seq.{table_name}")
        if table.table_type != "TABLE":
            raise ValueError(f"Target must be a native table: {table.full_table_id}")
        if table.require_partition_filter:
            raise ValueError(
                "Target requires a partition filter; refusing an unbounded dataset "
                f"replacement: {table.full_table_id}"
            )
        if table.streaming_buffer:
            raise ValueError(
                "Target has a streaming buffer; wait for it to clear before "
                f"replacement: {table.full_table_id}"
            )
        if _schema_signature(table.schema) != _schema_signature(
            _expected_bq_schema(schema, bigquery)
        ):
            raise ValueError(
                f"Target schema does not match pipeline {kind.upper()} schema: "
                f"{table.full_table_id}"
            )


def _table_id(project: str, table_name: str) -> str:
    return f"{project}.perturb_seq.{table_name}"


def _stage_id(project: str, run_key: str, kind: str) -> str:
    return _table_id(project, f"ingest_{run_key}_{kind}")


def transaction_sql(project: str, dea_stage: str, gsea_stage: str) -> str:
    operations = []
    for kind, (target_name, schema) in TARGETS.items():
        target = f"`{_table_id(project, target_name)}`"
        stage = f"`{_table_id(project, dea_stage if kind == 'dea' else gsea_stage)}`"
        columns = ", ".join(f"`{field.name}`" for field in schema)
        operations.extend(
            [
                f"DELETE FROM {target} WHERE dataset_id IN UNNEST(@dataset_ids);",
                f"INSERT INTO {target} ({columns}) SELECT {columns} FROM {stage} "
                "WHERE dataset_id IN UNNEST(@dataset_ids);",
            ]
        )
    # ponytail: BigQuery transactions do not serialize external writers; serialize
    # replacement runs operationally while this transaction protects both tables.
    return "BEGIN TRANSACTION;\n" + "\n".join(operations) + "\nCOMMIT TRANSACTION;"


def _counts(
    client: Any,
    table_id: str,
    location: str,
    bigquery: Any,
    dataset_ids: list[str] | None = None,
) -> dict[str, int]:
    where = ""
    job_config = None
    if dataset_ids is not None:
        where = " WHERE dataset_id IN UNNEST(@dataset_ids)"
        job_config = bigquery.QueryJobConfig(
            query_parameters=[
                bigquery.ArrayQueryParameter("dataset_ids", "STRING", dataset_ids)
            ]
        )
    result = client.query(
        f"SELECT dataset_id, COUNT(*) AS row_count FROM `{table_id}`{where} GROUP BY dataset_id",
        job_config=job_config,
        location=location,
    ).result()
    return {row["dataset_id"]: int(row["row_count"]) for row in result}


def _verify_counts(
    expected: dict[str, int], actual: dict[str, int], label: str
) -> None:
    unexpected = set(actual) - set(expected)
    if unexpected:
        raise RuntimeError(
            f"{label}: unexpected dataset IDs in staged/target data: "
            f"{sorted(unexpected)}"
        )
    mismatches = {
        dataset_id: (expected_count, actual.get(dataset_id, 0))
        for dataset_id, expected_count in expected.items()
        if actual.get(dataset_id, 0) != expected_count
    }
    if mismatches:
        raise RuntimeError(
            f"{label}: row-count mismatch (expected, actual): {mismatches}"
        )


def _upload_blob(
    storage_client: Any,
    bucket_name: str,
    name: str,
    path: Path,
    metadata: dict[str, str],
) -> str:
    from google.api_core.exceptions import PreconditionFailed

    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(name)
    blob.metadata = metadata
    try:
        blob.upload_from_filename(str(path), if_generation_match=0, checksum="auto")
    except PreconditionFailed:
        existing = bucket.get_blob(name)
        if (
            existing is None
            or existing.size != path.stat().st_size
            or not existing.metadata
            or existing.metadata.get("sha256") != metadata["sha256"]
        ):
            raise ValueError(
                "GCS object already exists with different content: "
                f"gs://{bucket_name}/{name}"
            )
    return f"gs://{bucket_name}/{name}"


def _expire_stage(client: Any, table_id: str) -> None:
    table = client.get_table(table_id)
    table.expires = datetime.now(timezone.utc) + STAGING_TTL
    client.update_table(table, ["expires"])


def _cleanup_uploaded_objects(
    storage_client: Any, bucket_name: str, names: list[str]
) -> None:
    for name in names:
        try:
            storage_client.bucket(bucket_name).blob(name).delete()
        except Exception as error:
            print(
                f"Warning: could not remove gs://{bucket_name}/{name}: {error}",
                file=sys.stderr,
            )


def apply_replacement(
    files: list[DatasetFiles], run_id: str, project: str | None = None
) -> None:
    project = project or os.getenv("GCLOUD_PROJECT", "")
    location = os.getenv("BQ_LOCATION", "")
    bucket_name = os.getenv("CLOUD_TMP_BUCKET") or os.getenv("GCLOUD_TMP_BUCKET", "")
    if not project or not PROJECT_PATTERN.fullmatch(project):
        raise ValueError("Set --project or GCLOUD_PROJECT to the dev project ID")
    if "prod" in project.casefold():
        raise ValueError(f"Refusing a production-like project: {project}")
    if not location:
        raise ValueError("Set BQ_LOCATION from the dev environment")
    if not bucket_name or "/" in bucket_name or bucket_name.startswith("gs:"):
        raise ValueError("Set CLOUD_TMP_BUCKET (or GCLOUD_TMP_BUCKET) to a bucket name")

    from google.api_core.exceptions import Conflict
    from google.cloud import bigquery, storage

    client = bigquery.Client(project=project)
    storage_client = storage.Client(project=project)
    _check_table(client, project, location, bigquery)
    bucket = storage_client.get_bucket(bucket_name)
    if (bucket.location or "").casefold() != location.casefold():
        raise ValueError(
            "GCS bucket location does not exactly match BQ_LOCATION: "
            f"{bucket.location!r} != {location!r}"
        )

    run_key, hashes = _run_key(run_id, files)
    dataset_ids = [item.dataset_id for item in files]
    stage_ids = {kind: _stage_id(project, run_key, kind) for kind in TARGETS}
    uris: dict[str, list[str]] = {kind: [] for kind in TARGETS}
    object_names: list[str] = []
    try:
        for item in files:
            for kind, path in (("dea", item.dea_path), ("gsea", item.gsea_path)):
                object_name = (
                    f"perturb-seq-ingest/{run_key}/{item.dataset_id}/{path.name}"
                )
                metadata = {
                    "sha256": hashes[(item.dataset_id, kind)],
                    "dataset_id": item.dataset_id,
                    "kind": kind,
                }
                uri = _upload_blob(
                    storage_client, bucket_name, object_name, path, metadata
                )
                object_names.append(object_name)
                uris[kind].append(uri)

        job_config = {
            "dea": bigquery.LoadJobConfig(
                source_format=bigquery.SourceFormat.PARQUET,
                write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
            ),
            "gsea": bigquery.LoadJobConfig(
                source_format=bigquery.SourceFormat.PARQUET,
                write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
            ),
        }
        local_counts = {
            "dea": {item.dataset_id: item.dea_rows for item in files},
            "gsea": {item.dataset_id: item.gsea_rows for item in files},
        }
        for kind, (target_name, schema) in TARGETS.items():
            load = client.load_table_from_uri(
                uris[kind],
                stage_ids[kind],
                job_config=job_config[kind],
                location=location,
            )
            load.result()
            if load.errors:
                raise RuntimeError(
                    f"BigQuery {kind.upper()} load had row errors: {load.errors}"
                )
            _expire_stage(client, stage_ids[kind])
            stage = client.get_table(stage_ids[kind])
            expected_schema = _expected_bq_schema(schema, bigquery)
            if _schema_signature(stage.schema) != _schema_signature(expected_schema):
                raise RuntimeError(
                    f"Staging schema mismatch for {kind.upper()}: {stage_ids[kind]}"
                )
            if stage.num_rows != sum(local_counts[kind].values()):
                raise RuntimeError(
                    f"Staging row count mismatch for {kind.upper()}: "
                    f"expected {sum(local_counts[kind].values())}, got {stage.num_rows}"
                )
            _verify_counts(
                local_counts[kind],
                _counts(client, stage_ids[kind], location, bigquery),
                f"{kind.upper()} staging",
            )

    finally:
        _cleanup_uploaded_objects(storage_client, bucket_name, object_names)

    sql = transaction_sql(
        project,
        stage_ids["dea"].split(".")[-1],
        stage_ids["gsea"].split(".")[-1],
    )
    query_config = bigquery.QueryJobConfig(
        query_parameters=[
            bigquery.ArrayQueryParameter("dataset_ids", "STRING", dataset_ids)
        ]
    )
    job_id = f"replace_{run_key}"
    try:
        transaction = client.query(
            sql, job_config=query_config, location=location, job_id=job_id
        )
    except Conflict:
        # The same content-derived job ID means an earlier attempt may have committed
        # before its caller lost the response; inspect that job first. If it failed,
        # replay the idempotent scoped replacement with a fresh attempt ID.
        transaction = client.get_job(job_id, location=location)
        if transaction.query != sql:
            raise RuntimeError(f"BigQuery job ID collision: {job_id}")
        if transaction.error_result:
            job_id = f"replace_{run_key}_{uuid.uuid4().hex[:8]}"
            transaction = client.query(
                sql, job_config=query_config, location=location, job_id=job_id
            )
    transaction.result()

    for kind, (target_name, _schema) in TARGETS.items():
        target = _table_id(project, target_name)
        _verify_counts(
            local_counts[kind],
            _counts(client, target, location, bigquery, dataset_ids),
            f"{kind.upper()} target",
        )

    print(
        f"Replaced {len(files)} dataset IDs in {project}.perturb_seq; run key {run_key}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        required=True,
        help="JSON manifest of dataset IDs and paired Parquet paths",
    )
    parser.add_argument(
        "--project",
        default=os.getenv("GCLOUD_PROJECT", ""),
        help="dev project ID (defaults to GCLOUD_PROJECT)",
    )
    parser.add_argument(
        "--run-id", help="Stable lowercase run label; required with --apply"
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Upload and replace rows; default is local dry run",
    )
    args = parser.parse_args(argv)

    files = validate_manifest(args.manifest)
    print(f"Validated {len(files)} datasets; local row counts:")
    for item in files:
        print(f"  {item.dataset_id}: DEA {item.dea_rows:,}, GSEA {item.gsea_rows:,}")
    if not args.apply:
        print("Dry run only. Pass --apply and --run-id to load these datasets.")
        return 0
    if not args.run_id:
        parser.error("--run-id is required with --apply")
    apply_replacement(files, args.run_id, project=args.project)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
