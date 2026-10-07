#!/usr/bin/env python3
"""Validate and upload selected DEA/GSEA Parquets from the cluster."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import base64
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import uuid

RESULTS_ROOT = "/hps/nobackup/mfreeberg/perturb_seq_fastq/results"
CACHE_ROOT = "/hps/nobackup/mfreeberg/cache"
DATASET_ID = re.compile(r"[A-Za-z0-9_-]+\Z")
PROJECT = re.compile(r"[A-Za-z0-9-]+\Z")
RUN_KEY = re.compile(r"[a-f0-9]{32}\Z")
EXPECTED_DATASET_COUNT = 20
CLI_WORKERS = 4


class TransferError(RuntimeError):
    pass


def _schema_signatures() -> dict[str, list[list[object]]]:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from io_schemas import DEA_SCHEMA, GSEA_SCHEMA

    return {
        "dea": _schema_signature(DEA_SCHEMA),
        "gsea": _schema_signature(GSEA_SCHEMA),
    }


def _schema_signature(schema) -> list[list[object]]:
    import pyarrow as pa

    def type_signature(value) -> str:
        if pa.types.is_list(value):
            return f"list<{value.value_type}>"
        if pa.types.is_large_list(value):
            return f"large_list<{value.value_type}>"
        return str(value)

    return [
        [field.name, type_signature(field.type), field.nullable] for field in schema
    ]


def _read_source_manifest_from_stream(stream) -> list[dict[str, object]]:
    try:
        datasets = json.load(stream)["datasets"]
    except (json.JSONDecodeError, KeyError, TypeError) as error:
        raise ValueError(
            f"Cannot read source manifest: {type(error).__name__}"
        ) from error
    if not isinstance(datasets, list) or not datasets:
        raise ValueError("Source manifest must contain a nonempty datasets list")

    seen: set[str] = set()
    output = []
    for item in datasets:
        if not isinstance(item, dict):
            raise ValueError("Each source manifest entry must be an object")
        dataset_id = item.get("dataset_id")
        if not isinstance(dataset_id, str) or not DATASET_ID.fullmatch(dataset_id):
            raise ValueError(f"Invalid dataset_id in source manifest: {dataset_id!r}")
        if dataset_id in seen:
            raise ValueError(f"Duplicate dataset_id in source manifest: {dataset_id}")
        seen.add(dataset_id)
        entry: dict[str, object] = {"dataset_id": dataset_id}
        for kind in ("dea", "gsea"):
            key = f"{kind}_parquet"
            size_key = f"{kind}_size_bytes"
            value = item.get(key)
            size = item.get(size_key)
            rows = item.get(f"{kind}_rows")
            md5_hash = item.get(f"{kind}_md5_hash")
            sha256 = item.get(f"{kind}_sha256")
            source_stat = item.get(f"{kind}_stat")
            path_value = PurePosixPath(value) if isinstance(value, str) else None
            if (
                path_value is None
                or not path_value.is_absolute()
                or ".." in path_value.parts
                or not str(path_value).startswith(RESULTS_ROOT + "/")
                or path_value.name != f"{dataset_id}.{kind}.parquet"
                or not isinstance(size, int)
                or isinstance(size, bool)
                or size <= 0
                or not isinstance(rows, int)
                or isinstance(rows, bool)
                or rows < 0
                or not isinstance(md5_hash, str)
                or not re.fullmatch(r"[A-Za-z0-9+/]{22}==", md5_hash)
                or not isinstance(sha256, str)
                or not re.fullmatch(r"[a-f0-9]{64}", sha256)
                or not isinstance(source_stat, dict)
                or any(
                    not isinstance(source_stat.get(key), int)
                    for key in ("dev", "ino", "mtime_ns")
                )
            ):
                raise ValueError(
                    f"{dataset_id}: invalid cluster source path or size for {kind}"
                )
            entry[key] = str(path_value)
            entry[size_key] = size
            entry[f"{kind}_rows"] = rows
            entry[f"{kind}_md5_hash"] = md5_hash
            entry[f"{kind}_sha256"] = sha256
            entry[f"{kind}_stat"] = source_stat
        output.append(entry)
    return output


def _private_cache_path(value: str) -> Path:
    path = Path(value)
    if (
        not path.is_absolute()
        or ".." in path.parts
        or not path.is_relative_to(CACHE_ROOT)
        or not path.parent.resolve(strict=True).is_relative_to(CACHE_ROOT)
    ):
        raise TransferError(f"Expected a path under {CACHE_ROOT}")
    return path


def inventory(dataset_ids: list[str], output: Path) -> None:
    if len(dataset_ids) != EXPECTED_DATASET_COUNT:
        raise ValueError(f"Expected exactly {EXPECTED_DATASET_COUNT} dataset IDs")
    if not dataset_ids or any(not DATASET_ID.fullmatch(value) for value in dataset_ids):
        raise ValueError("Provide one or more comma-separated dataset IDs")
    if len(dataset_ids) != len(set(dataset_ids)):
        raise ValueError("Dataset IDs must be unique")
    output = _private_cache_path(str(output))
    if output.exists() or output.is_symlink():
        raise ValueError(f"Refusing to overwrite source manifest: {output}")

    root = Path(RESULTS_ROOT)
    schemas = _schema_signatures()
    datasets = []
    for dataset_id in dataset_ids:
        entry: dict[str, object] = {"dataset_id": dataset_id}
        paths = {}
        for kind in ("dea", "gsea"):
            name = f"{dataset_id}.{kind}.parquet"
            matches = [
                *root.glob(f"*/dea_gsea/{name}"),
                *root.glob(f"*/*/dea_gsea/{name}"),
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected exactly one {kind.upper()} Parquet for {dataset_id}; "
                    f"found {len(matches)}"
                )
            path = matches[0]
            info = path.lstat()
            if path.is_symlink() or not path.is_file() or info.st_size <= 0:
                raise ValueError(f"Invalid result file for {dataset_id}: {path}")
            resolved = path.resolve(strict=True)
            if not resolved.is_relative_to(RESULTS_ROOT):
                raise ValueError(f"Result path escapes the results tree: {path}")
            stream, rows, source_stat = _verify_source(
                str(resolved), dataset_id, kind, schemas[kind]
            )
            try:
                if (
                    info.st_dev,
                    info.st_ino,
                    info.st_size,
                    info.st_mtime_ns,
                ) != (
                    source_stat.st_dev,
                    source_stat.st_ino,
                    source_stat.st_size,
                    source_stat.st_mtime_ns,
                ):
                    raise TransferError(
                        f"Source changed during inventory: {dataset_id} {kind}"
                    )
                md5 = hashlib.md5(usedforsecurity=False)
                sha256 = hashlib.sha256()
                for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                    md5.update(block)
                    sha256.update(block)
                current = os.fstat(stream.fileno())
                if (
                    source_stat.st_dev,
                    source_stat.st_ino,
                    source_stat.st_size,
                    source_stat.st_mtime_ns,
                ) != (
                    current.st_dev,
                    current.st_ino,
                    current.st_size,
                    current.st_mtime_ns,
                ):
                    raise TransferError(
                        f"Source changed during inventory: {dataset_id} {kind}"
                    )
            finally:
                stream.close()
            entry[f"{kind}_parquet"] = str(resolved)
            entry[f"{kind}_size_bytes"] = info.st_size
            entry[f"{kind}_rows"] = rows
            entry[f"{kind}_md5_hash"] = base64.b64encode(md5.digest()).decode("ascii")
            entry[f"{kind}_sha256"] = sha256.hexdigest()
            entry[f"{kind}_stat"] = {
                "dev": source_stat.st_dev,
                "ino": source_stat.st_ino,
                "mtime_ns": source_stat.st_mtime_ns,
            }
            paths[kind] = resolved.parent
        if paths["dea"] != paths["gsea"]:
            raise ValueError(
                f"DEA and GSEA products are in different directories for {dataset_id}"
            )
        datasets.append(entry)

    raw = json.dumps({"datasets": datasets}, separators=(",", ":")).encode()
    fd = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "wb") as stream:
        stream.write(raw)
    print(
        f"Inventoried {len(datasets)} datasets ({len(raw)} bytes of path/size metadata)"
    )


def _run_cli(args: list[str]) -> str:
    try:
        result = subprocess.run(args, check=True, capture_output=True, text=True)
    except FileNotFoundError as error:
        raise TransferError(
            f"Required cluster CLI is not on PATH: {args[0]}"
        ) from error
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "").strip()
        raise TransferError(
            f"{args[0]} failed: {detail or error.returncode}"
        ) from error
    return result.stdout


def _read_json_cli(args: list[str]) -> dict:
    try:
        value = json.loads(_run_cli(args))
    except json.JSONDecodeError as error:
        raise TransferError(f"{args[0]} returned invalid JSON") from error
    if not isinstance(value, dict):
        raise TransferError(f"{args[0]} returned an unexpected response")
    return value


def _verify_source(path_value: str, dataset_id: str, kind: str, expected_schema):
    import pyarrow.parquet as pq

    path = Path(path_value)
    if not path.is_absolute() or ".." in path.parts:
        raise TransferError(f"Invalid source path for {dataset_id} {kind}")
    original = path.lstat()
    if (
        not path.is_relative_to(RESULTS_ROOT)
        or path.name != f"{dataset_id}.{kind}.parquet"
    ):
        raise TransferError(
            f"Source path is outside the permitted results tree: {path}"
        )
    if not path.is_file() or path.is_symlink():
        raise TransferError(f"Source is not a regular file: {path}")
    if not path.resolve(strict=True).is_relative_to(RESULTS_ROOT):
        raise TransferError(
            f"Source resolves outside the permitted results tree: {path}"
        )

    parquet = pq.ParquetFile(path)
    signature = _schema_signature(parquet.schema_arrow)
    if signature != expected_schema:
        raise TransferError(f"Parquet schema mismatch for {dataset_id} {kind}")
    rows = parquet.metadata.num_rows
    for batch in parquet.iter_batches(columns=["dataset_id"], batch_size=1_000_000):
        if any(value != dataset_id for value in batch.column(0).to_pylist()):
            raise TransferError(f"Wrong dataset_id in {dataset_id} {kind} Parquet")

    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    stream = os.fdopen(fd, "rb")
    current = os.fstat(stream.fileno())
    if (
        original.st_dev,
        original.st_ino,
        original.st_size,
        original.st_mtime_ns,
    ) != (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns):
        stream.close()
        raise TransferError(f"Source changed during validation: {path}")
    return stream, rows, current


def _upload_one(
    run_key: str,
    bucket: str,
    project: str,
    dataset_id: str,
    kind: str,
    source: dict,
) -> dict[str, object]:
    path_value = source.get(f"{kind}_parquet")
    if not isinstance(path_value, str):
        raise TransferError(f"Invalid source path for {dataset_id} {kind}")
    path = Path(path_value)
    if (
        not path.is_absolute()
        or ".." in path.parts
        or not path.is_relative_to(RESULTS_ROOT)
        or path.name != f"{dataset_id}.{kind}.parquet"
        or path.is_symlink()
        or not path.is_file()
        or not path.resolve(strict=True).is_relative_to(RESULTS_ROOT)
    ):
        raise TransferError(f"Invalid source path for {dataset_id} {kind}")
    original = path.stat()
    expected_stat = source.get(f"{kind}_stat")
    if (
        original.st_size != source.get(f"{kind}_size_bytes")
        or not isinstance(expected_stat, dict)
        or (original.st_dev, original.st_ino, original.st_mtime_ns)
        != (
            expected_stat.get("dev"),
            expected_stat.get("ino"),
            expected_stat.get("mtime_ns"),
        )
    ):
        raise TransferError(f"Source changed after validation: {dataset_id} {kind}")
    rows = source[f"{kind}_rows"]
    md5_hash = source[f"{kind}_md5_hash"]
    sha256_hash = source[f"{kind}_sha256"]
    object_name = (
        f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
    )
    uri = f"gs://{bucket}/{object_name}"
    metadata = {
        "run_key": run_key,
        "dataset_id": dataset_id,
        "kind": kind,
        "sha256": sha256_hash,
        "row_count": str(rows),
    }
    _run_cli(
        [
            "gcloud",
            "storage",
            "cp",
            f"--project={project}",
            "--if-generation-match=0",
            f"--content-md5={md5_hash}",
            "--content-type=application/vnd.apache.parquet",
            "--custom-metadata="
            + ",".join(f"{key}={value}" for key, value in metadata.items()),
            str(path),
            uri,
        ]
    )
    uploaded = _read_json_cli(
        [
            "gcloud",
            "storage",
            "objects",
            "describe",
            uri,
            f"--project={project}",
            "--format=json",
        ]
    )
    generation = int(uploaded.get("generation", 0))
    if (
        uploaded.get("name") != object_name
        or uploaded.get("bucket") != bucket
        or int(uploaded.get("size", -1)) != original.st_size
        or generation <= 0
        or (uploaded.get("md5Hash") is not None and uploaded["md5Hash"] != md5_hash)
        or any(
            (uploaded.get("metadata") or {}).get(key) != value
            for key, value in metadata.items()
        )
    ):
        raise TransferError(f"GCS verification failed for {dataset_id} {kind}")
    return {
        "uri": uri,
        "size_bytes": original.st_size,
        "row_count": rows,
        "sha256": sha256_hash,
        "md5_hash": md5_hash,
        "crc32c": uploaded.get("crc32c"),
        "generation": str(generation),
    }


def _delete_uploaded(
    project: str, completed: list[tuple[str, dict[str, object]]]
) -> None:
    for uri, receipt in completed:
        try:
            _run_cli(
                [
                    "gcloud",
                    "storage",
                    "rm",
                    uri,
                    f"--project={project}",
                    f"--if-generation-match={receipt['generation']}",
                ]
            )
        except TransferError as error:
            print(f"Warning: could not remove {uri}: {error}", file=sys.stderr)


def upload(
    manifest_path: Path,
    receipt_path: Path,
    project: str,
    bucket: str,
    delete_worker: bool = False,
) -> None:
    manifest_path = _private_cache_path(str(manifest_path))
    receipt_path = _private_cache_path(str(receipt_path))
    if not PROJECT.fullmatch(project) or "prod" in project.casefold():
        raise ValueError("Refusing an empty or production-like Google Cloud project")
    if project != os.getenv("GCLOUD_PROJECT"):
        raise ValueError("--project must match GCLOUD_PROJECT from the dev environment")
    configured_bucket = os.getenv("CLOUD_TMP_BUCKET") or os.getenv(
        "GCLOUD_TMP_BUCKET", ""
    )
    if (
        not bucket
        or not re.fullmatch(r"[A-Za-z0-9._-]+", bucket)
        or bucket != configured_bucket
    ):
        raise ValueError(
            "--bucket must match CLOUD_TMP_BUCKET from the dev environment"
        )
    location = os.getenv("BQ_LOCATION", "")
    if not location:
        raise ValueError("Set BQ_LOCATION from the development environment")
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise ValueError("Expected an inventory manifest")
    if receipt_path.exists() or receipt_path.is_symlink():
        raise ValueError(f"Refusing to overwrite receipt: {receipt_path}")
    try:
        fd = os.open(manifest_path, os.O_RDONLY | os.O_NOFOLLOW)
        with os.fdopen(fd, "r") as stream:
            info = os.fstat(stream.fileno())
            if not stat.S_ISREG(info.st_mode) or info.st_mode & 0o077:
                raise ValueError("Inventory manifest must be a private regular file")
            datasets = _read_source_manifest_from_stream(stream)
    except (OSError, ValueError) as error:
        raise TransferError(f"Cannot read inventory manifest: {error}") from error
    if len(datasets) != EXPECTED_DATASET_COUNT:
        raise TransferError(
            f"Expected exactly {EXPECTED_DATASET_COUNT} datasets in inventory"
        )

    bucket_info = _read_json_cli(
        [
            "gcloud",
            "storage",
            "buckets",
            "describe",
            f"gs://{bucket}",
            f"--project={project}",
            "--format=json",
        ]
    )
    if (bucket_info.get("location") or "").casefold() != location.casefold():
        raise ValueError(
            f"GCS bucket location {bucket_info.get('location')!r} does not match "
            f"BQ_LOCATION {location!r}"
        )

    run_key = uuid.uuid4().hex
    work = [
        (item["dataset_id"], kind, item)
        for item in datasets
        for kind in ("dea", "gsea")
    ]
    completed: list[tuple[str, dict[str, object]]] = []
    errors = []
    receipts = {item["dataset_id"]: {} for item in datasets}
    with ThreadPoolExecutor(max_workers=CLI_WORKERS) as pool:
        futures = {
            pool.submit(
                _upload_one,
                run_key,
                bucket,
                project,
                dataset_id,
                kind,
                source,
            ): (dataset_id, kind)
            for dataset_id, kind, source in work
        }
        for future in as_completed(futures):
            dataset_id, kind = futures[future]
            try:
                value = future.result()
                receipts[dataset_id][kind] = value
                completed.append((value["uri"], value))
                print(f"Uploaded {dataset_id} {kind}", flush=True)
            except Exception as error:
                errors.append(error)
    if errors:
        _delete_uploaded(project, completed)
        raise errors[0]

    receipt = {
        "run_key": run_key,
        "bucket": bucket,
        "datasets": [
            {"dataset_id": item["dataset_id"], **receipts[item["dataset_id"]]}
            for item in datasets
        ],
    }
    raw = json.dumps(receipt, separators=(",", ":")).encode()
    fd = os.open(
        receipt_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600
    )
    with os.fdopen(fd, "wb") as stream:
        stream.write(raw)

    if delete_worker:
        worker = _private_cache_path(str(Path(__file__).resolve()))
        if worker.is_symlink() or not worker.is_file():
            raise TransferError("Temporary cluster worker is not a regular cache file")
        worker.unlink()
    print(f"Verified and uploaded {len(work)} Parquet files; receipt {receipt_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    inventory_parser = subparsers.add_parser(
        "inventory", help="locate final products and record sizes on cluster"
    )
    inventory_parser.add_argument("--dataset-ids", required=True)
    inventory_parser.add_argument("--output", type=Path, required=True)
    upload_parser = subparsers.add_parser(
        "upload", help="run in a gcloud-auth cluster job"
    )
    upload_parser.add_argument("--manifest", type=Path, required=True)
    upload_parser.add_argument("--receipt", type=Path, required=True)
    upload_parser.add_argument("--project", default=os.getenv("GCLOUD_PROJECT", ""))
    upload_parser.add_argument("--bucket", default=os.getenv("CLOUD_TMP_BUCKET", ""))
    upload_parser.add_argument("--delete-worker", action="store_true")
    args = parser.parse_args()
    if args.mode == "inventory":
        inventory(args.dataset_ids.split(","), args.output)
    else:
        upload(
            args.manifest, args.receipt, args.project, args.bucket, args.delete_worker
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (TransferError, ValueError, OSError) as error:
        print(f"Cluster Parquet transfer failed: {error}", file=sys.stderr)
        raise SystemExit(1)
