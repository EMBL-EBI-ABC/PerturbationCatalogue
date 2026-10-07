#!/usr/bin/env python3
"""Upload the selected cluster Parquets to GCS and replace their BQ rows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import uuid

RESULTS_ROOT = Path("/hps/nobackup/mfreeberg/perturb_seq_fastq/results")
EXPECTED_DATASET_COUNT = 20
DATASET_ID = re.compile(r"[A-Za-z0-9_-]+\Z")
TABLES = {"dea": "pertpy_dea", "gsea": "pertpy_gsea"}


def read_ids(manifest: Path) -> list[str]:
    try:
        values = [
            item["dataset_id"] for item in json.loads(manifest.read_text())["datasets"]
        ]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
        raise ValueError(f"Cannot read dataset manifest: {error}") from error
    if (
        len(values) != EXPECTED_DATASET_COUNT
        or any(
            not isinstance(value, str) or not DATASET_ID.fullmatch(value)
            for value in values
        )
        or len(set(values)) != EXPECTED_DATASET_COUNT
    ):
        raise ValueError(
            f"Manifest must contain {EXPECTED_DATASET_COUNT} unique dataset IDs"
        )
    return values


def result_files(dataset_ids: list[str]) -> list[Path]:
    files = []
    for dataset_id in dataset_ids:
        for kind in TABLES:
            name = f"{dataset_id}.{kind}.parquet"
            matches = [
                *RESULTS_ROOT.glob(f"*/dea_gsea/{name}"),
                *RESULTS_ROOT.glob(f"*/*/dea_gsea/{name}"),
            ]
            if len(matches) != 1 or not matches[0].is_file():
                raise FileNotFoundError(
                    f"Expected one {kind.upper()} Parquet for {dataset_id}; found {len(matches)}"
                )
            files.append(matches[0])
    return files


def publish(dataset_ids: list[str], project: str, location: str, bucket: str) -> str:
    if not project or "prod" in project.casefold():
        raise ValueError(f"Refusing non-development project: {project!r}")
    if not location or not bucket:
        raise ValueError("BigQuery location and temporary bucket are required")
    if (
        len(dataset_ids) != EXPECTED_DATASET_COUNT
        or len(set(dataset_ids)) != EXPECTED_DATASET_COUNT
    ):
        raise ValueError(f"Expected {EXPECTED_DATASET_COUNT} unique dataset IDs")

    run_prefix = f"perturb-seq-ingest/{uuid.uuid4().hex}"
    bucket_prefix = f"gs://{bucket}/{run_prefix}"
    files = result_files(dataset_ids)
    print(f"Uploading {len(files)} Parquets to {bucket_prefix}/", flush=True)
    subprocess.run(
        [
            "gcloud",
            "storage",
            "cp",
            *(str(path) for path in files),
            bucket_prefix + "/",
        ],
        check=True,
    )

    ids_sql = ", ".join(f"'{dataset_id}'" for dataset_id in dataset_ids)
    delete_sql = f"""BEGIN TRANSACTION;
DELETE FROM `{project}.perturb_seq.pertpy_dea` WHERE dataset_id IN ({ids_sql});
DELETE FROM `{project}.perturb_seq.pertpy_gsea` WHERE dataset_id IN ({ids_sql});
COMMIT TRANSACTION;"""
    subprocess.run(
        [
            "bq",
            f"--project_id={project}",
            f"--location={location}",
            "query",
            "--use_legacy_sql=false",
            delete_sql,
        ],
        check=True,
    )

    for kind, table in TABLES.items():
        subprocess.run(
            [
                "bq",
                f"--project_id={project}",
                f"--location={location}",
                "load",
                "--source_format=PARQUET",
                "--parquet_enable_list_inference=true",
                f"{project}:perturb_seq.{table}",
                f"{bucket_prefix}/*.{kind}.parquet",
            ],
            check=True,
        )

    subprocess.run(
        ["gcloud", "storage", "rm", "--recursive", bucket_prefix + "/"], check=True
    )
    print("BigQuery loads succeeded; temporary Parquets removed.")
    return bucket_prefix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--location", required=True)
    parser.add_argument("--bucket", required=True)
    args = parser.parse_args()
    publish(read_ids(args.manifest), args.project, args.location, args.bucket)


if __name__ == "__main__":
    main()
