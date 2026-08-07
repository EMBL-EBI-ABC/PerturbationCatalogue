"""Build one metadata, CSV.GZ and Parquet artifact per dataset."""

import argparse
import csv
import gzip
import json
import logging
import os
import shutil
import tempfile
import uuid
from collections import defaultdict
from pathlib import Path

from google.cloud import bigquery, storage
import pyarrow as pa
import pyarrow.parquet as pq


DATASETS = {
    "crispr": {
        "data": "crispr.data",
        "metadata": "crispr.metadata",
        "fields": [
            ("perturbed_target_ensg", "Perturbed Target ENSG"),
            ("perturbed_target_name", "Perturbed Target Name"),
            ("score_name", "Score Name"),
            ("score_value", "Score Value"),
            ("significant", "Significant"),
            ("significance_criteria", "Significance Criteria"),
        ],
    },
    "perturb-seq": {
        "data": "perturb_seq.pertpy_dea",
        "metadata": "perturb_seq.metadata",
        "fields": [
            ("perturbed_target_ensg", "Perturbed Target ENSG"),
            ("perturbed_target_name", "Perturbed Target Name"),
            ("effect_gene_ensg", "Effect Gene ENSG"),
            ("effect_gene_name", "Effect Gene Name"),
            ("log2foldchange", "Log2FC"),
            ("padj", "Padj"),
            ("score_name", "Score Name"),
            ("score_value", "Score Value"),
            ("cell_type", "Cell Type"),
        ],
    },
    "mave": {
        "data": "mavedb.data",
        "metadata": "mavedb.metadata",
        "fields": [
            ("perturbed_target_ensg", "Perturbed Target ENSG"),
            ("perturbed_target_name", "Perturbed Target Name"),
            ("perturbation_name", "Perturbation Name"),
            ("perturbation_position", "Position"),
            ("perturbation_aa_wt", "AA WT"),
            ("perturbation_aa_change", "AA Change"),
            ("score_name", "Score Name"),
            ("score_value", "Score Value"),
        ],
    },
}


def _table(project, dataset, name):
    return f"`{project}.{name}`" if "." in name else f"`{project}.{dataset}.{name}`"


def _data_query(project, dataset, modality):
    config = DATASETS[modality]
    source = _table(project, dataset, config["data"])
    fields = [name for name, _ in config["fields"]]
    select = []
    for field in fields:
        if field == "perturbed_target_name":
            select.append("d.perturbed_target_symbol AS perturbed_target_name")
        elif field == "effect_gene_name":
            select.append("d.effect_gene_symbol AS effect_gene_name")
        elif field == "perturbation_position":
            select.append(
                "SAFE_CAST(REGEXP_EXTRACT(d.perturbation_name, r'p\\.[A-Za-z]+(\\d+)') AS INT64) AS perturbation_position"
            )
        elif field == "perturbation_aa_wt":
            select.append(
                "REGEXP_EXTRACT(d.perturbation_name, r'p\\.([A-Za-z]+)\\d+') AS perturbation_aa_wt"
            )
        elif field == "perturbation_aa_change":
            select.append(
                "REGEXP_EXTRACT(d.perturbation_name, r'p\\.[A-Za-z]+\\d+([A-Za-z=]+)') AS perturbation_aa_change"
            )
        else:
            select.append(f"d.{field}")
    return "SELECT " + ", ".join(["d.dataset_id"] + select) + f" FROM {source} d"


def _metadata_query(project, dataset, modality):
    return f"SELECT * FROM {_table(project, dataset, DATASETS[modality]['metadata'])}"


def _export_data_shards(client, query, bq_dataset, location, bucket_name):
    table = client.dataset(bq_dataset).table(f"release_{uuid.uuid4().hex}")
    prefix = f"release-tmp/{uuid.uuid4().hex}"
    try:
        job_config = bigquery.QueryJobConfig(
            destination=table,
            write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
        )
        client.query(query, job_config=job_config, location=location).result()
        extract = client.extract_table(
            table,
            f"gs://{bucket_name}/{prefix}/part-*.parquet",
            location=location,
            job_config=bigquery.ExtractJobConfig(destination_format="PARQUET"),
        )
        extract.result()
    finally:
        client.delete_table(table, not_found_ok=True)
    return prefix


def _rows(client, query, location):
    result = client.query(query, location=location).result()
    return [dict(row) for row in result]


def _write_data(shards, fields, dataset_ids, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_files = {}
    csv_writers = {}
    parquet_writers = {}
    paths = {}
    data_schema = None
    try:
        for dataset_id in dataset_ids:
            dataset_dir = output_dir / dataset_id
            dataset_dir.mkdir()
            paths[dataset_id] = {
                "parquet": dataset_dir / f"{dataset_id}.parquet",
                "csv.gz": dataset_dir / f"{dataset_id}.csv.gz",
            }
            csv_files[dataset_id] = gzip.open(
                paths[dataset_id]["csv.gz"], "wt", newline=""
            )
            csv_writers[dataset_id] = csv.writer(csv_files[dataset_id])
            csv_writers[dataset_id].writerow([label for _, label in fields])

        for blob in shards:
            shard_path = output_dir / blob.name.rsplit("/", 1)[-1]
            blob.download_to_filename(shard_path)
            parquet = pq.ParquetFile(shard_path)
            data_schema = parquet.schema_arrow
            data_schema = data_schema.remove(data_schema.get_field_index("dataset_id"))
            for batch in parquet.iter_batches(batch_size=100_000):
                grouped = defaultdict(list)
                for row in batch.to_pylist():
                    grouped[row.pop("dataset_id")].append(row)
                for dataset_id, rows in grouped.items():
                    if dataset_id not in paths:
                        continue
                    table = pa.Table.from_pylist(rows, schema=data_schema)
                    writer = parquet_writers.get(dataset_id)
                    if writer is None:
                        writer = pq.ParquetWriter(
                            paths[dataset_id]["parquet"],
                            data_schema,
                            compression="zstd",
                        )
                        parquet_writers[dataset_id] = writer
                    writer.write_table(table)
                    csv_writers[dataset_id].writerows(
                        [row.get(field) for field, _ in fields] for row in rows
                    )
    finally:
        for writer in parquet_writers.values():
            writer.close()
        for csv_file in csv_files.values():
            csv_file.close()
        if data_schema is not None:
            for dataset_id in dataset_ids:
                if dataset_id not in parquet_writers:
                    pq.ParquetWriter(
                        paths[dataset_id]["parquet"], data_schema, compression="zstd"
                    ).close()
    return paths


def _write_metadata(rows, output_dir):
    path = output_dir / "metadata.json"
    with path.open("w") as handle:
        json.dump(rows, handle, indent=2, default=str)
    return path


def build_release(bq_project, bq_dataset, bq_location, bucket_name, release_prefix):
    client = bigquery.Client(project=bq_project)
    bucket = storage.Client(project=bq_project).bucket(bucket_name)
    for modality, config in DATASETS.items():
        data_table = _table(bq_project, bq_dataset, config["data"])
        metadata_table = _table(bq_project, bq_dataset, config["metadata"])
        ids = {
            row.dataset_id
            for row in client.query(
                f"SELECT DISTINCT dataset_id FROM {data_table} UNION DISTINCT SELECT DISTINCT dataset_id FROM {metadata_table}",
                location=bq_location,
            ).result()
            if row.dataset_id
        }
        if not ids:
            continue
        logging.info("Building %s artifacts for %d datasets", modality, len(ids))
        prefix = _export_data_shards(
            client,
            _data_query(bq_project, bq_dataset, modality),
            bq_dataset,
            bq_location,
            bucket_name,
        )
        output_dir = Path(tempfile.mkdtemp(prefix="release-"))
        try:
            paths = _write_data(
                list(bucket.list_blobs(prefix=prefix)),
                config["fields"],
                sorted(ids),
                output_dir,
            )
            metadata_by_dataset = defaultdict(list)
            for row in _rows(
                client, _metadata_query(bq_project, bq_dataset, modality), bq_location
            ):
                if row.get("dataset_id") in ids:
                    metadata_by_dataset[row["dataset_id"]].append(row)
            for dataset_id in sorted(ids):
                metadata_path = _write_metadata(
                    metadata_by_dataset[dataset_id], output_dir / dataset_id
                )
                object_prefix = f"{release_prefix}/{modality}/{dataset_id}"
                for path, suffix, content_type in (
                    (metadata_path, "metadata.json", "application/json"),
                    (paths[dataset_id]["csv.gz"], "csv.gz", "application/gzip"),
                    (
                        paths[dataset_id]["parquet"],
                        "parquet",
                        "application/vnd.apache.parquet",
                    ),
                ):
                    bucket.blob(f"{object_prefix}.{suffix}").upload_from_filename(
                        path, content_type=content_type
                    )
        finally:
            for blob in bucket.list_blobs(prefix=prefix):
                blob.delete()
            shutil.rmtree(output_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bq-project", default=os.getenv("GCLOUD_PROJECT"), required=False
    )
    parser.add_argument("--bq-dataset", default=os.getenv("BQ_DATASET"), required=False)
    parser.add_argument(
        "--bq-location", default=os.getenv("BQ_LOCATION"), required=False
    )
    parser.add_argument(
        "--gcs-bucket",
        default=os.getenv("CLOUD_TMP_BUCKET") or os.getenv("GCLOUD_TMP_BUCKET"),
        required=False,
    )
    parser.add_argument("--release-prefix", default="release")
    args = parser.parse_args()
    missing = [
        name
        for name in ("bq_project", "bq_dataset", "bq_location", "gcs_bucket")
        if not getattr(args, name)
    ]
    if missing:
        parser.error("missing: " + ", ".join(missing))
    build_release(
        args.bq_project,
        args.bq_dataset,
        args.bq_location,
        args.gcs_bucket,
        args.release_prefix,
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )
    main()
