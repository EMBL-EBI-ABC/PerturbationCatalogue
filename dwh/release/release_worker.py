"""Generate one dataset's release artifacts in a Cloud Run Job task."""

import argparse
import csv
import gzip
import io
import json
import logging
import os
from contextlib import ExitStack, contextmanager
from itertools import chain

from google.api_core.exceptions import NotFound
from google.cloud import bigquery, storage

try:
    from google.cloud import bigquery_storage
except ImportError:  # The REST fallback keeps local checks dependency-light.
    bigquery_storage = None
import pyarrow as pa
import pyarrow.parquet as pq

from release import DATASETS, table


BATCH_SIZE = 100_000
CHUNK_SIZE = 8 * 1024 * 1024


def _split_uri(uri):
    bucket, _, path = uri.removeprefix("gs://").partition("/")
    if not bucket or not path:
        raise ValueError(f"Invalid GCS URI: {uri}")
    return bucket, path


def _arrow_type(field):
    types = {
        "BOOL": pa.bool_(),
        "BYTES": pa.binary(),
        "DATE": pa.date32(),
        "FLOAT": pa.float64(),
        "FLOAT64": pa.float64(),
        "INT64": pa.int64(),
        "INTEGER": pa.int64(),
        "NUMERIC": pa.decimal128(38, 9),
        "STRING": pa.string(),
        "TIME": pa.time64("us"),
        "TIMESTAMP": pa.timestamp("us", tz="UTC"),
    }
    return types.get(field.field_type, pa.string())


def _empty_schema(client, relation):
    fields = client.get_table(relation).schema
    return pa.schema(
        [
            (field.name, _arrow_type(field))
            for field in fields
            if field.name != "dataset_id"
        ]
    )


def _batches(result, bqstorage_client):
    try:
        yield from result.to_arrow_iterable(
            bqstorage_client=bqstorage_client, max_stream_count=1
        )
        return
    except (ImportError, ValueError):
        rows = []
        for row in result:
            rows.append(dict(row))
            if len(rows) == BATCH_SIZE:
                yield from pa.Table.from_pylist(rows).to_batches()
                rows = []
        if rows:
            yield from pa.Table.from_pylist(rows).to_batches()


@contextmanager
def _query(
    client, relation, dataset_id, location, include_dataset_id=False, kind="data"
):
    config = bigquery.QueryJobConfig(
        query_parameters=[
            bigquery.ScalarQueryParameter("dataset_id", "STRING", dataset_id)
        ]
    )
    project, dataset, _ = relation.split(".", 2)
    safe_id = "".join(
        char if char.isalnum() or char == "_" else "_" for char in dataset_id
    )
    destination = f"{project}.{dataset}.release_query_{safe_id}_{kind}"
    config.destination = destination
    config.write_disposition = bigquery.WriteDisposition.WRITE_TRUNCATE
    job = client.query(
        f"SELECT {'*' if include_dataset_id else '* EXCEPT(dataset_id)'} "
        f"FROM `{relation}` "
        "WHERE dataset_id = @dataset_id",
        job_config=config,
        location=location,
    )
    job.result(max_results=1)
    try:
        yield client.list_rows(destination)
    finally:
        client.delete_table(destination, not_found_ok=True)


def _write_data(bucket, prefix, client, item, location, bqstorage_client):
    config = DATASETS[item["modality"]]
    with _query(client, item["data_table"], item["dataset_id"], location) as result:
        batches = iter(_batches(result, bqstorage_client))
        first = next(batches, None)
        schema = (
            first.schema
            if first is not None
            else _empty_schema(client, item["data_table"])
        )
        parquet_blob = bucket.blob(f"{prefix}.parquet")
        csv_blob = bucket.blob(f"{prefix}.csv.gz")
        with ExitStack() as stack:
            parquet_stream = stack.enter_context(
                parquet_blob.open(
                    "wb",
                    chunk_size=CHUNK_SIZE,
                    ignore_flush=True,
                    content_type="application/vnd.apache.parquet",
                )
            )
            parquet_writer = stack.enter_context(
                pq.ParquetWriter(parquet_stream, schema, compression="zstd")
            )
            csv_stream = stack.enter_context(
                csv_blob.open(
                    "wb",
                    chunk_size=CHUNK_SIZE,
                    ignore_flush=True,
                    content_type="application/gzip",
                )
            )
            gzip_stream = stack.enter_context(
                gzip.GzipFile(fileobj=csv_stream, mode="wb")
            )
            csv_text = stack.enter_context(io.TextIOWrapper(gzip_stream, newline=""))
            writer = csv.writer(csv_text)
            writer.writerow([label for _, label in config["fields"]])
            for batch in chain((first,), batches):
                if batch is None:
                    continue
                parquet_writer.write_batch(batch)
                writer.writerows(
                    [row.get(field) for field, _ in config["fields"]]
                    for row in batch.to_pylist()
                )


def _write_metadata(bucket, prefix, client, item, location, bqstorage_client):
    with _query(
        client,
        item["metadata_table"],
        item["dataset_id"],
        location,
        include_dataset_id=True,
        kind="metadata",
    ) as result:
        with bucket.blob(f"{prefix}.metadata.json").open(
            "wb",
            chunk_size=CHUNK_SIZE,
            ignore_flush=True,
            content_type="application/json",
        ) as output:
            row = next(
                (
                    row
                    for batch in _batches(result, bqstorage_client)
                    for row in batch.to_pylist()
                ),
                {},
            )
            output.write(json.dumps(row, default=str, separators=(",", ":")).encode())


def generate(item, bucket, prefix, client, location, bqstorage_client):
    names = [f"{prefix}.{suffix}" for suffix in ("metadata.json", "csv.gz", "parquet")]
    try:
        _write_data(bucket, prefix, client, item, location, bqstorage_client)
        _write_metadata(bucket, prefix, client, item, location, bqstorage_client)
    except Exception:
        for name in names:
            try:
                bucket.blob(name).delete()
            except NotFound:
                pass
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest-uri", default=os.getenv("RELEASE_MANIFEST_URI"), required=False
    )
    parser.add_argument(
        "--project", default=os.getenv("GCLOUD_PROJECT"), required=False
    )
    parser.add_argument("--location", default=os.getenv("BQ_LOCATION"), required=False)
    parser.add_argument("--bucket", default=os.getenv("RELEASE_BUCKET"), required=False)
    parser.add_argument("--prefix", default=os.getenv("RELEASE_PREFIX", "release"))
    args = parser.parse_args()
    if not all((args.manifest_uri, args.project, args.location, args.bucket)):
        parser.error("manifest URI, project, location and bucket are required")

    task_index = int(os.getenv("CLOUD_RUN_TASK_INDEX", "0"))
    manifest_bucket, manifest_path = _split_uri(args.manifest_uri)
    storage_client = storage.Client(project=args.project)
    manifest = json.loads(
        storage_client.bucket(manifest_bucket).blob(manifest_path).download_as_bytes()
    )
    if task_index >= len(manifest["items"]):
        logging.info("No dataset assigned to task %d", task_index)
        return
    item = manifest["items"][task_index]
    prefix = f"{args.prefix}/{item['modality']}/{item['dataset_id']}"
    client = bigquery.Client(project=args.project)
    bqstorage_client = (
        bigquery_storage.BigQueryReadClient() if bigquery_storage else None
    )
    try:
        logging.info("Building %s", prefix)
        generate(
            item,
            storage_client.bucket(args.bucket),
            prefix,
            client,
            args.location,
            bqstorage_client,
        )
    finally:
        if bqstorage_client:
            getattr(bqstorage_client, "close", lambda: None)()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )
    main()
