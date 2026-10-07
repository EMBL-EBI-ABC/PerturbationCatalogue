#!/usr/bin/env python3
"""Validate a cluster upload receipt and replace its selected BigQuery rows."""

from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlsplit

DEA_SCHEMA = [
    {"name": "dataset_id", "type": "STRING", "mode": "NULLABLE"},
    {"name": "perturbed_target_symbol", "type": "STRING", "mode": "NULLABLE"},
    {"name": "perturbed_target_ensg", "type": "STRING", "mode": "NULLABLE"},
    {"name": "effect_gene_symbol", "type": "STRING", "mode": "NULLABLE"},
    {"name": "effect_gene_ensg", "type": "STRING", "mode": "NULLABLE"},
    {"name": "padj", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "log2foldchange", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "score_name", "type": "STRING", "mode": "NULLABLE"},
    {"name": "score_value", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "cell_type", "type": "STRING", "mode": "NULLABLE"},
    {"name": "max_ingested_at", "type": "TIMESTAMP", "mode": "NULLABLE"},
]
GSEA_SCHEMA = [
    {"name": "dataset_id", "type": "STRING", "mode": "NULLABLE"},
    {"name": "term", "type": "STRING", "mode": "NULLABLE"},
    {"name": "perturbed_target_symbol", "type": "STRING", "mode": "NULLABLE"},
    {"name": "perturbed_target_ensg", "type": "STRING", "mode": "NULLABLE"},
    {"name": "es", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "nes", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "pval", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "sidak", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "fdr", "type": "FLOAT64", "mode": "NULLABLE"},
    {"name": "geneset_size", "type": "INT64", "mode": "NULLABLE"},
    {"name": "leading_edge", "type": "STRING", "mode": "REPEATED"},
    {"name": "cell_type", "type": "STRING", "mode": "NULLABLE"},
    {"name": "max_ingested_at", "type": "TIMESTAMP", "mode": "NULLABLE"},
]
TARGETS = {"dea": ("pertpy_dea", DEA_SCHEMA), "gsea": ("pertpy_gsea", GSEA_SCHEMA)}
DATASET_ID_PATTERN = re.compile(r"[A-Za-z0-9_-]+\Z")
PROJECT_PATTERN = re.compile(r"[A-Za-z0-9-]+\Z")
RUN_KEY_PATTERN = re.compile(r"[a-f0-9]{32}\Z")
EXPECTED_DATASET_COUNT = 20
STAGING_TTL_SECONDS = 2 * 24 * 60 * 60


@dataclass(frozen=True)
class GcsParquet:
    uri: str
    object_name: str
    size_bytes: int
    rows: int
    sha256: str
    md5_hash: str
    crc32c: str | None
    generation: int


@dataclass(frozen=True)
class GcsDatasetFiles:
    dataset_id: str
    dea: GcsParquet
    gsea: GcsParquet


@dataclass(frozen=True)
class GcsManifest:
    run_key: str
    bucket: str
    datasets: list[GcsDatasetFiles]


def _run_cli(args: list[str]) -> str:
    try:
        result = subprocess.run(args, check=True, capture_output=True, text=True)
    except FileNotFoundError as error:
        raise RuntimeError(f"Required cluster CLI is not on PATH: {args[0]}") from error
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "").strip()
        raise RuntimeError(f"{args[0]} failed: {detail or error.returncode}") from error
    return result.stdout


def _json_cli(args: list[str]) -> dict:
    try:
        value = json.loads(_run_cli(args))
    except json.JSONDecodeError as error:
        raise RuntimeError(f"{args[0]} returned invalid JSON") from error
    if not isinstance(value, dict):
        raise RuntimeError(f"{args[0]} returned an unexpected response")
    return value


def _bq_args(project: str, location: str, command: str) -> list[str]:
    return [
        "bq",
        f"--project_id={project}",
        f"--location={location}",
        command,
    ]


def validate_gcs_manifest(manifest_path: str | Path) -> GcsManifest:
    manifest_path = Path(manifest_path).resolve()
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read GCS receipt {manifest_path}: {error}") from error
    if not isinstance(manifest, dict):
        raise ValueError("GCS receipt must be a JSON object")
    run_key = manifest.get("run_key")
    bucket = manifest.get("bucket")
    entries = manifest.get("datasets")
    if not isinstance(run_key, str) or not RUN_KEY_PATTERN.fullmatch(run_key):
        raise ValueError("GCS receipt has an invalid run_key")
    if not isinstance(bucket, str) or not re.fullmatch(r"[A-Za-z0-9._-]+", bucket):
        raise ValueError("GCS receipt has an invalid bucket name")
    if not isinstance(entries, list) or not entries:
        raise ValueError("GCS receipt must contain a nonempty datasets list")

    seen: set[str] = set()
    datasets = []
    for item in entries:
        dataset_id = item.get("dataset_id") if isinstance(item, dict) else None
        if not isinstance(dataset_id, str) or not DATASET_ID_PATTERN.fullmatch(
            dataset_id
        ):
            raise ValueError(f"Invalid dataset_id in GCS receipt: {dataset_id!r}")
        if dataset_id in seen:
            raise ValueError(f"Duplicate dataset_id in GCS receipt: {dataset_id!r}")
        seen.add(dataset_id)
        parquet_files = {}
        for kind in ("dea", "gsea"):
            value = item.get(kind)
            if not isinstance(value, dict):
                raise ValueError(f"{dataset_id}: missing {kind} upload receipt")
            name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
            )
            uri = value.get("uri")
            parsed = urlsplit(uri) if isinstance(uri, str) else None
            if (
                parsed is None
                or parsed.scheme != "gs"
                or parsed.netloc != bucket
                or parsed.path != "/" + name
                or parsed.query
                or parsed.fragment
                or uri != f"gs://{bucket}/{name}"
            ):
                raise ValueError(f"{dataset_id}: unexpected GCS URI for {kind}")
            size = value.get("size_bytes")
            rows = value.get("row_count")
            generation = value.get("generation")
            if not isinstance(size, int) or isinstance(size, bool) or size <= 0:
                raise ValueError(f"{dataset_id}: invalid {kind} object size")
            if not isinstance(rows, int) or isinstance(rows, bool) or rows < 0:
                raise ValueError(f"{dataset_id}: invalid {kind} row count")
            try:
                generation = int(generation)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{dataset_id}: invalid {kind} object generation"
                ) from error
            if generation <= 0:
                raise ValueError(f"{dataset_id}: invalid {kind} object generation")
            sha256 = value.get("sha256")
            if not isinstance(sha256, str) or not re.fullmatch(r"[a-f0-9]{64}", sha256):
                raise ValueError(f"{dataset_id}: invalid {kind} SHA-256")
            md5_hash = value.get("md5_hash")
            if not isinstance(md5_hash, str):
                raise ValueError(f"{dataset_id}: invalid {kind} md5_hash")
            try:
                decoded_md5 = base64.b64decode(md5_hash, validate=True)
            except (ValueError, TypeError) as error:
                raise ValueError(f"{dataset_id}: invalid {kind} md5_hash") from error
            if len(decoded_md5) != 16:
                raise ValueError(f"{dataset_id}: invalid {kind} md5_hash")
            crc32c = value.get("crc32c")
            if crc32c is not None:
                if not isinstance(crc32c, str):
                    raise ValueError(f"{dataset_id}: invalid {kind} crc32c")
                try:
                    decoded_crc32c = base64.b64decode(crc32c, validate=True)
                except (ValueError, TypeError) as error:
                    raise ValueError(f"{dataset_id}: invalid {kind} crc32c") from error
                if len(decoded_crc32c) != 4:
                    raise ValueError(f"{dataset_id}: invalid {kind} crc32c")
            parquet_files[kind] = GcsParquet(
                uri=uri,
                object_name=name,
                size_bytes=size,
                rows=rows,
                sha256=sha256,
                md5_hash=md5_hash,
                crc32c=crc32c,
                generation=generation,
            )
        datasets.append(
            GcsDatasetFiles(dataset_id, parquet_files["dea"], parquet_files["gsea"])
        )
    return GcsManifest(run_key, bucket, datasets)


def validate_expected_datasets(id_manifest_path: str | Path) -> list[str]:
    try:
        entries = json.loads(Path(id_manifest_path).read_text())["datasets"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
        raise ValueError(f"Cannot read selected dataset manifest: {error}") from error
    if not isinstance(entries, list) or len(entries) != EXPECTED_DATASET_COUNT:
        raise ValueError(
            f"Selected manifest must contain exactly {EXPECTED_DATASET_COUNT} datasets"
        )
    dataset_ids = [
        item.get("dataset_id") if isinstance(item, dict) else None for item in entries
    ]
    if any(
        not isinstance(value, str) or not DATASET_ID_PATTERN.fullmatch(value)
        for value in dataset_ids
    ):
        raise ValueError("Selected manifest contains an invalid dataset_id")
    if len(set(dataset_ids)) != EXPECTED_DATASET_COUNT:
        raise ValueError("Selected manifest dataset IDs must be unique")
    return dataset_ids


def _schema_signature(fields: list[dict]) -> tuple[tuple[str, str, str], ...]:
    aliases = {"FLOAT": "FLOAT64", "INTEGER": "INT64", "BOOLEAN": "BOOL"}
    return tuple(
        sorted(
            (
                field["name"],
                aliases.get(str(field["type"]).upper(), str(field["type"]).upper()),
                str(field.get("mode") or "NULLABLE").upper(),
            )
            for field in fields
        )
    )


def _check_table(project: str, location: str) -> None:
    dataset = _json_cli(
        [
            "bq",
            f"--project_id={project}",
            "show",
            "--dataset=true",
            "--format=prettyjson",
            f"{project}:perturb_seq",
        ]
    )
    if (dataset.get("location") or "").casefold() != location.casefold():
        raise ValueError(
            f"BQ dataset location {dataset.get('location')!r} does not match BQ_LOCATION {location!r}"
        )
    for kind, (table_name, schema) in TARGETS.items():
        table = _json_cli(
            [
                "bq",
                f"--project_id={project}",
                "show",
                "--format=prettyjson",
                f"{project}:perturb_seq.{table_name}",
            ]
        )
        if table.get("type") != "TABLE":
            raise ValueError(
                f"Target must be a native table: {project}.perturb_seq.{table_name}"
            )
        if table.get("requirePartitionFilter"):
            raise ValueError(
                f"Target requires a partition filter; refusing an unbounded dataset replacement: {project}.perturb_seq.{table_name}"
            )
        if table.get("streamingBuffer"):
            raise ValueError(
                f"Target has a streaming buffer; wait for it to clear before replacement: {project}.perturb_seq.{table_name}"
            )
        if _schema_signature(
            table.get("schema", {}).get("fields", [])
        ) != _schema_signature(schema):
            raise ValueError(
                f"Target schema does not match pipeline {kind.upper()} schema: {project}.perturb_seq.{table_name}"
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
        columns = ", ".join(f"`{field['name']}`" for field in schema)
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


def _query(
    project: str,
    location: str,
    sql: str,
    dataset_ids: list[str] | None = None,
    output_format: str = "prettyjson",
) -> str:
    command = _bq_args(project, location, "query")
    command.append("--use_legacy_sql=false")
    command.append(f"--format={output_format}")
    if dataset_ids is not None:
        command.append(
            "--parameter=dataset_ids:ARRAY<STRING>:" + json.dumps(dataset_ids)
        )
    command.append(sql)
    return _run_cli(command)


def _counts(
    project: str,
    table_id: str,
    location: str,
    dataset_ids: list[str] | None = None,
) -> dict[str, int]:
    where = ""
    if dataset_ids is not None:
        where = " WHERE dataset_id IN UNNEST(@dataset_ids)"
    rows = json.loads(
        _query(
            project,
            location,
            f"SELECT dataset_id, COUNT(*) AS row_count FROM `{table_id}`{where} GROUP BY dataset_id",
            dataset_ids,
        )
    )
    if not isinstance(rows, list):
        raise RuntimeError(f"Unexpected BigQuery count response for {table_id}")
    return {row["dataset_id"]: int(row["row_count"]) for row in rows}


def _verify_counts(
    expected: dict[str, int], actual: dict[str, int], label: str
) -> None:
    unexpected = set(actual) - set(expected)
    if unexpected:
        raise RuntimeError(
            f"{label}: unexpected dataset IDs in staged/target data: {sorted(unexpected)}"
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


def _verify_gcs_objects(manifest: GcsManifest, project: str) -> dict[str, int]:
    generations = {}
    for dataset in manifest.datasets:
        for kind, parquet in (("dea", dataset.dea), ("gsea", dataset.gsea)):
            blob = _json_cli(
                [
                    "gcloud",
                    "storage",
                    "objects",
                    "describe",
                    parquet.uri,
                    f"--project={project}",
                    "--format=json",
                ]
            )
            metadata = blob.get("metadata") or {}
            if (
                int(blob.get("generation", 0)) != parquet.generation
                or int(blob.get("size", -1)) != parquet.size_bytes
                or (
                    blob.get("md5Hash") is not None
                    and blob["md5Hash"] != parquet.md5_hash
                )
                or blob.get("crc32c") != parquet.crc32c
                or metadata.get("run_key") != manifest.run_key
                or metadata.get("dataset_id") != dataset.dataset_id
                or metadata.get("kind") != kind
                or metadata.get("sha256") != parquet.sha256
                or metadata.get("row_count") != str(parquet.rows)
            ):
                raise ValueError(
                    f"GCS object verification failed for gs://{manifest.bucket}/{parquet.object_name}"
                )
            generations[parquet.object_name] = parquet.generation
    return generations


def _load_stage(
    project: str,
    location: str,
    run_key: str,
    kind: str,
    schema: list[dict[str, str]],
    files: list[GcsParquet],
) -> str:
    stage = _stage_id(project, run_key, kind)
    _run_cli(
        _bq_args(project, location, "load")
        + [
            "--source_format=PARQUET",
            "--replace=true",
            "--parquet_enable_list_inference=true",
            stage.replace(".", ":", 1),
            ",".join(item.uri for item in files),
        ]
    )
    _run_cli(
        _bq_args(project, location, "update")
        + [f"--expiration={STAGING_TTL_SECONDS}", stage.replace(".", ":", 1)]
    )
    table = _json_cli(
        [
            "bq",
            f"--project_id={project}",
            "show",
            "--format=prettyjson",
            stage.replace(".", ":", 1),
        ]
    )
    if _schema_signature(
        table.get("schema", {}).get("fields", [])
    ) != _schema_signature(schema):
        raise RuntimeError(f"Staging schema mismatch for {kind.upper()}: {stage}")
    expected_rows = sum(item.rows for item in files)
    if int(table.get("numRows", -1)) != expected_rows:
        raise RuntimeError(
            f"Staging row count mismatch for {kind.upper()}: "
            f"expected {expected_rows}, got {table.get('numRows')}"
        )
    return stage


def _replace_rows(
    project: str,
    location: str,
    run_key: str,
    dataset_ids: list[str],
    row_counts: dict[str, dict[str, int]],
    dea_stage: str,
    gsea_stage: str,
) -> None:
    for kind, stage in (("dea", dea_stage), ("gsea", gsea_stage)):
        _verify_counts(
            row_counts[kind],
            _counts(project, stage, location),
            f"{kind.upper()} staging",
        )
    _query(
        project,
        location,
        transaction_sql(
            project, dea_stage.rsplit(".", 1)[-1], gsea_stage.rsplit(".", 1)[-1]
        ),
        dataset_ids,
    )
    for kind, (target_name, _) in TARGETS.items():
        _verify_counts(
            row_counts[kind],
            _counts(project, _table_id(project, target_name), location, dataset_ids),
            f"{kind.upper()} target",
        )


def apply_gcs_replacement(manifest: GcsManifest, project: str | None = None) -> None:
    project = project or os.getenv("GCLOUD_PROJECT", "")
    location = os.getenv("BQ_LOCATION", "")
    bucket_name = os.getenv("CLOUD_TMP_BUCKET") or os.getenv("GCLOUD_TMP_BUCKET", "")
    if not PROJECT_PATTERN.fullmatch(project) or "prod" in project.casefold():
        raise ValueError("Set --project to the dev project ID; production is refused")
    if project != os.getenv("GCLOUD_PROJECT"):
        raise ValueError("--project must match GCLOUD_PROJECT from the dev environment")
    if not location:
        raise ValueError("Set BQ_LOCATION from the dev environment")
    if bucket_name != manifest.bucket:
        raise ValueError("GCS receipt bucket does not match CLOUD_TMP_BUCKET")
    dataset_ids = [item.dataset_id for item in manifest.datasets]
    if (
        len(dataset_ids) != EXPECTED_DATASET_COUNT
        or len(set(dataset_ids)) != EXPECTED_DATASET_COUNT
    ):
        raise ValueError(
            f"Replacement must contain exactly {EXPECTED_DATASET_COUNT} unique datasets"
        )

    _check_table(project, location)
    generations = _verify_gcs_objects(manifest, project)
    files = {
        "dea": [item.dea for item in manifest.datasets],
        "gsea": [item.gsea for item in manifest.datasets],
    }
    row_counts = {
        kind: {item.dataset_id: getattr(item, kind).rows for item in manifest.datasets}
        for kind in TARGETS
    }
    stages = {
        kind: _load_stage(
            project, location, manifest.run_key, kind, schema, files[kind]
        )
        for kind, (_, schema) in TARGETS.items()
    }
    current = _verify_gcs_objects(manifest, project)
    if current != generations:
        raise ValueError("GCS source generations changed during BigQuery staging")

    _replace_rows(
        project,
        location,
        manifest.run_key,
        dataset_ids,
        row_counts,
        stages["dea"],
        stages["gsea"],
    )
    for name, generation in generations.items():
        try:
            _run_cli(
                [
                    "gcloud",
                    "storage",
                    "rm",
                    f"gs://{bucket_name}/{name}",
                    f"--project={project}",
                    f"--if-generation-match={generation}",
                ]
            )
        except RuntimeError as error:
            print(
                f"Warning: could not remove gs://{bucket_name}/{name}: {error}",
                file=sys.stderr,
            )
    print(
        f"Replaced {len(dataset_ids)} dataset IDs in {project}.perturb_seq; "
        f"run key {manifest.run_key}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gcs-manifest", required=True, help="URL-free cluster receipt"
    )
    parser.add_argument(
        "--id-manifest", required=True, help="the selected 20-dataset manifest"
    )
    parser.add_argument(
        "--project", default=os.getenv("GCLOUD_PROJECT", ""), help="dev project ID"
    )
    parser.add_argument(
        "--apply", action="store_true", help="load and replace these datasets"
    )
    args = parser.parse_args(argv)

    manifest = validate_gcs_manifest(args.gcs_manifest)
    expected_ids = validate_expected_datasets(args.id_manifest)
    receipt_ids = [item.dataset_id for item in manifest.datasets]
    if set(receipt_ids) != set(expected_ids):
        raise ValueError(
            "GCS receipt dataset IDs do not match the selected 20-dataset manifest"
        )
    print(f"Validated {len(manifest.datasets)} datasets; cluster-verified row counts:")
    for item in manifest.datasets:
        print(f"  {item.dataset_id}: DEA {item.dea.rows:,}, GSEA {item.gsea.rows:,}")
    if not args.apply:
        print("Dry run only. Pass --apply to load these datasets.")
        return 0
    apply_gcs_replacement(manifest, project=args.project)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, OSError) as error:
        print(f"BigQuery replacement failed: {error}", file=sys.stderr)
        raise SystemExit(1)
