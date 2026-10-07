#!/usr/bin/env python3
"""Prepare resumable GCS sessions locally and upload Parquet files from Slurm."""

from __future__ import annotations

import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener
import uuid

import pyarrow.parquet as pq


RESULTS_ROOT = "/hps/nobackup/mfreeberg/perturb_seq_fastq/results"
CACHE_ROOT = "/hps/nobackup/mfreeberg/cache"
DATASET_ID = re.compile(r"[A-Za-z0-9_-]+\Z")
RUN_KEY = re.compile(r"[a-f0-9]{32}\Z")
SHA256 = re.compile(r"[a-f0-9]{64}\Z")
CHUNK_SIZE = 32 * 1024 * 1024  # GCS requires multiples of 256 KiB.
MAX_MANIFEST_BYTES = 32 * 1024  # Restricted cluster connector --upload limit.
MAX_RETRIES = 6
HTTP_TIMEOUT = 300


class TransferError(RuntimeError):
    pass


class RetryableTransferError(TransferError):
    pass


def _schema_signatures() -> dict[str, list[list[object]]]:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from io_schemas import DEA_SCHEMA, GSEA_SCHEMA

    return {
        "dea": [[field.name, str(field.type), field.nullable] for field in DEA_SCHEMA],
        "gsea": [
            [field.name, str(field.type), field.nullable] for field in GSEA_SCHEMA
        ],
    }


def _read_source_manifest(path: Path) -> list[dict[str, object]]:
    try:
        datasets = json.loads(path.read_text())["datasets"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
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
            path_value = PurePosixPath(value) if isinstance(value, str) else None
            if (
                path_value is None
                or not path_value.is_absolute()
                or ".." in path_value.parts
                or not str(path_value).startswith(RESULTS_ROOT + "/")
                or path_value.name != f"{dataset_id}.{kind}.parquet"
                or not isinstance(size, int)
                or size <= 0
            ):
                raise ValueError(
                    f"{dataset_id}: invalid cluster source path or size for {kind}"
                )
            entry[key] = str(path_value)
            entry[size_key] = size
        output.append(entry)
    return output


def inventory(dataset_ids: list[str], output: Path) -> None:
    if not dataset_ids or any(not DATASET_ID.fullmatch(value) for value in dataset_ids):
        raise ValueError("Provide one or more comma-separated dataset IDs")
    if len(dataset_ids) != len(set(dataset_ids)):
        raise ValueError("Dataset IDs must be unique")
    output = _private_cache_path(str(output))
    if output.exists() or output.is_symlink():
        raise ValueError(f"Refusing to overwrite source manifest: {output}")

    root = Path(RESULTS_ROOT)
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
            entry[f"{kind}_parquet"] = str(resolved)
            entry[f"{kind}_size_bytes"] = info.st_size
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


def prepare(
    source_manifest: Path, project: str, bucket_name: str, output: Path
) -> None:
    if not project or "prod" in project.casefold():
        raise ValueError("Refusing an empty or production-like Google Cloud project")
    if project != os.getenv("GCLOUD_PROJECT"):
        raise ValueError("--project must match GCLOUD_PROJECT from the dev environment")
    if not bucket_name or "/" in bucket_name or bucket_name.startswith("gs:"):
        raise ValueError("Provide a bucket name, not a gs:// URI")
    configured_bucket = os.getenv("CLOUD_TMP_BUCKET") or os.getenv(
        "GCLOUD_TMP_BUCKET", ""
    )
    if not configured_bucket or bucket_name != configured_bucket:
        raise ValueError(
            "--bucket must match CLOUD_TMP_BUCKET from the dev environment"
        )
    location = os.getenv("BQ_LOCATION", "")
    if not location:
        raise ValueError("Set BQ_LOCATION from the development environment")
    if output.exists() or output.is_symlink():
        raise ValueError(f"Refusing to overwrite session manifest: {output}")

    datasets = _read_source_manifest(source_manifest)
    from google.cloud import storage

    storage_client = storage.Client(project=project)
    bucket = storage_client.get_bucket(bucket_name)
    if (bucket.location or "").casefold() != location.casefold():
        raise ValueError(
            f"GCS bucket location {bucket.location!r} does not match BQ_LOCATION {location!r}"
        )

    run_key = uuid.uuid4().hex
    sessions = []
    for item in datasets:
        dataset_id = item["dataset_id"]
        entry: dict[str, object] = {"dataset_id": dataset_id}
        for kind in ("dea", "gsea"):
            source_path = item[f"{kind}_parquet"]
            size = item[f"{kind}_size_bytes"]
            object_name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{Path(source_path).name}"
            )
            blob = bucket.blob(object_name)
            blob.metadata = {"run_key": run_key, "dataset_id": dataset_id, "kind": kind}
            session_uri = blob.create_resumable_upload_session(
                content_type="application/vnd.apache.parquet",
                size=size,
                if_generation_match=0,
            )
            entry[kind] = {
                "source_path": source_path,
                "size_bytes": size,
                "object_name": object_name,
                "session_uri": session_uri,
            }
        sessions.append(entry)

    content = json.dumps(
        {
            "run_key": run_key,
            "bucket": bucket_name,
            "schemas": _schema_signatures(),
            "datasets": sessions,
        },
        separators=(",", ":"),
    ).encode()
    if len(content) > MAX_MANIFEST_BYTES:
        raise ValueError(
            f"Private transfer manifest is {len(content)} bytes; connector limit is "
            f"{MAX_MANIFEST_BYTES} bytes"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "wb") as stream:
        stream.write(content)
    print(
        f"Prepared {len(sessions) * 2} resumable sessions for {len(datasets)} datasets; "
        f"run key {run_key}, private manifest {output} ({len(content)} bytes)."
    )


def cleanup(manifest_path: Path, project: str, bucket_name: str) -> None:
    if not project or "prod" in project.casefold():
        raise ValueError("Refusing cleanup for an empty or production-like project")
    if project != os.getenv("GCLOUD_PROJECT"):
        raise ValueError("--project must match GCLOUD_PROJECT from the dev environment")
    configured_bucket = os.getenv("CLOUD_TMP_BUCKET") or os.getenv(
        "GCLOUD_TMP_BUCKET", ""
    )
    if not configured_bucket or bucket_name != configured_bucket:
        raise ValueError(
            "--bucket must match CLOUD_TMP_BUCKET from the dev environment"
        )
    if (
        manifest_path.is_symlink()
        or not manifest_path.is_file()
        or manifest_path.stat().st_mode & 0o077
    ):
        raise ValueError("Expected the original private transfer manifest")
    manifest = json.loads(manifest_path.read_text())
    run_key = manifest.get("run_key")
    if (
        not isinstance(run_key, str)
        or not re.fullmatch(r"[a-f0-9]{32}", run_key)
        or manifest.get("bucket") != bucket_name
        or not isinstance(manifest.get("datasets"), list)
    ):
        raise ValueError("Transfer manifest does not match the requested dev bucket")
    from google.cloud import storage

    bucket = storage.Client(project=project).bucket(bucket_name)
    removed = 0
    for item in manifest["datasets"]:
        dataset_id = item.get("dataset_id") if isinstance(item, dict) else None
        if not isinstance(dataset_id, str) or not DATASET_ID.fullmatch(dataset_id):
            raise ValueError("Invalid dataset ID in transfer manifest")
        for kind in ("dea", "gsea"):
            source = item.get(kind)
            name = (
                f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
            )
            if not isinstance(source, dict) or source.get("object_name") != name:
                raise ValueError("Unexpected object name in transfer manifest")
            blob = bucket.get_blob(name)
            if blob is None:
                continue
            metadata = blob.metadata or {}
            if (
                metadata.get("run_key") != run_key
                or metadata.get("dataset_id") != dataset_id
                or metadata.get("kind") != kind
            ):
                raise ValueError(
                    f"Refusing to delete an unverified object: gs://{bucket_name}/{name}"
                )
            blob.delete(if_generation_match=blob.generation)
            removed += 1
    print(f"Removed {removed} completed objects from this transfer run")


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, msg, headers, new_url):
        return None


_OPENER = build_opener(_NoRedirect())


def _http_put(url: str, body: bytes, headers: dict[str, str]):
    try:
        with _OPENER.open(
            Request(url, data=body, headers=headers, method="PUT"), timeout=HTTP_TIMEOUT
        ) as response:
            return response.status, response.headers, response.read()
    except HTTPError as error:
        return error.code, error.headers, error.read()
    except (URLError, TimeoutError, OSError) as error:
        raise RetryableTransferError(type(error).__name__) from None


def _offset(headers, total: int) -> int:
    value = headers.get("Range")
    if value is None:
        return 0
    match = re.fullmatch(r"bytes=0-(\d+)", value.strip())
    if not match:
        raise TransferError("GCS returned an invalid resumable-upload Range")
    offset = int(match.group(1)) + 1
    if offset > total:
        raise TransferError("GCS resumable-upload offset exceeds source size")
    return offset


def _response(status, headers, body: bytes, total: int):
    if status == 308:
        return _offset(headers, total), None
    if status in (200, 201):
        try:
            value = json.loads(body)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            raise TransferError(
                "GCS final upload response was not valid JSON"
            ) from error
        return total, value
    if status in (408, 429, 500, 502, 503, 504):
        raise RetryableTransferError(f"GCS HTTP {status}")
    raise TransferError(f"GCS HTTP {status}")


def _status(session_uri: str, total: int):
    for attempt in range(MAX_RETRIES):
        try:
            status, headers, body = _http_put(
                session_uri,
                b"",
                {"Content-Length": "0", "Content-Range": f"bytes */{total}"},
            )
            return _response(status, headers, body, total)
        except RetryableTransferError:
            if attempt + 1 == MAX_RETRIES:
                raise
            time.sleep(min(2**attempt, 16))
    raise AssertionError("unreachable")


def _hash_prefix(stream, length: int):
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    stream.seek(0)
    remaining = length
    while remaining:
        chunk = stream.read(min(8 * 1024 * 1024, remaining))
        if not chunk:
            raise TransferError("Source became shorter while verifying upload state")
        md5.update(chunk)
        sha256.update(chunk)
        remaining -= len(chunk)
    return md5, sha256


def _verify_source(path_value: str, dataset_id: str, kind: str, expected_schema):
    path = Path(path_value)
    if not path.is_absolute() or ".." in path.parts:
        raise TransferError(f"Invalid source path for {dataset_id} {kind}")
    original = path.lstat()
    if (
        not path.is_relative_to(RESULTS_ROOT)
        or not path.name == f"{dataset_id}.{kind}.parquet"
    ):
        raise TransferError(
            f"Source path is outside the permitted results tree: {path}"
        )
    if not path.is_file() or path.is_symlink():
        raise TransferError(f"Source is not a regular file: {path}")
    resolved = path.resolve(strict=True)
    if not resolved.is_relative_to(RESULTS_ROOT):
        raise TransferError(
            f"Source resolves outside the permitted results tree: {path}"
        )

    parquet = pq.ParquetFile(path)
    signature = [
        [field.name, str(field.type), field.nullable] for field in parquet.schema_arrow
    ]
    if signature != expected_schema:
        raise TransferError(f"Parquet schema mismatch for {dataset_id} {kind}")
    rows = parquet.metadata.num_rows
    for batch in parquet.iter_batches(columns=["dataset_id"], batch_size=1_000_000):
        if any(value != dataset_id for value in batch.column(0).to_pylist()):
            raise TransferError(f"Wrong dataset_id in {dataset_id} {kind} Parquet")

    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    stream = os.fdopen(fd, "rb")
    current = os.fstat(stream.fileno())
    if (original.st_dev, original.st_ino) != (current.st_dev, current.st_ino):
        stream.close()
        raise TransferError(f"Source changed during validation: {path}")
    return stream, rows, current


def _upload_one(
    run_key: str, bucket: str, dataset_id: str, kind: str, source: dict, schema
):
    session_uri = source.get("session_uri")
    parsed = urlsplit(session_uri) if isinstance(session_uri, str) else None
    if (
        parsed is None
        or parsed.scheme != "https"
        or parsed.hostname not in {"storage.googleapis.com", "www.googleapis.com"}
        or parsed.port not in (None, 443)
        or parsed.username
        or parsed.password
        or parsed.fragment
    ):
        raise TransferError(f"Invalid GCS upload session for {dataset_id} {kind}")
    expected_name = (
        f"perturb-seq-ingest/{run_key}/{dataset_id}/{dataset_id}.{kind}.parquet"
    )
    session_query = parse_qs(parsed.query, keep_blank_values=True)
    if (
        source.get("object_name") != expected_name
        or parsed.path != f"/upload/storage/v1/b/{bucket}/o"
        or session_query.get("uploadType") != ["resumable"]
        or session_query.get("name") != [expected_name]
        or len(session_query.get("upload_id", [])) != 1
        or not session_query["upload_id"][0]
    ):
        raise TransferError(f"Unexpected GCS object name for {dataset_id} {kind}")

    stream, rows, original = _verify_source(
        source.get("source_path", ""), dataset_id, kind, schema
    )
    try:
        total = original.st_size
        if total <= 0 or source.get("size_bytes") != total:
            raise TransferError(f"Source size changed for {dataset_id} {kind}")
        offset, final_object = _status(session_uri, total)
        md5, sha256 = _hash_prefix(stream, offset)
        stream.seek(offset)
        retries_without_progress = 0
        while final_object is None:
            start = offset
            chunk = stream.read(min(CHUNK_SIZE, total - start))
            if not chunk:
                raise TransferError(f"Unexpected end of {dataset_id} {kind} Parquet")
            headers = {
                "Content-Length": str(len(chunk)),
                "Content-Range": f"bytes {start}-{start + len(chunk) - 1}/{total}",
            }
            if start + len(chunk) == total:
                full_md5 = md5.copy()
                full_md5.update(chunk)
                headers["X-Goog-Hash"] = "md5=" + base64.b64encode(
                    full_md5.digest()
                ).decode("ascii")
            try:
                status, response_headers, body = _http_put(session_uri, chunk, headers)
                offset, final_object = _response(status, response_headers, body, total)
            except RetryableTransferError:
                time.sleep(1)
                offset, final_object = _status(session_uri, total)

            if offset < start or offset > start + len(chunk):
                raise TransferError(f"Invalid upload progress for {dataset_id} {kind}")
            accepted = offset - start
            md5.update(chunk[:accepted])
            sha256.update(chunk[:accepted])
            if accepted == 0 and final_object is None:
                retries_without_progress += 1
                if retries_without_progress >= MAX_RETRIES:
                    raise TransferError(
                        f"GCS made no upload progress for {dataset_id} {kind}"
                    )
                time.sleep(min(2 ** (retries_without_progress - 1), 16))
            else:
                retries_without_progress = 0
            stream.seek(offset)

        current = os.fstat(stream.fileno())
        if (
            original.st_dev,
            original.st_ino,
            original.st_size,
            original.st_mtime_ns,
        ) != (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns):
            raise TransferError(f"Source changed during upload: {dataset_id} {kind}")
        md5_hash = base64.b64encode(md5.digest()).decode("ascii")
        response_md5 = final_object.get("md5Hash")
        if (
            final_object.get("name") != expected_name
            or final_object.get("bucket") != bucket
            or int(final_object.get("size", -1)) != total
            or (response_md5 is not None and response_md5 != md5_hash)
            or int(final_object.get("generation", 0)) <= 0
        ):
            raise TransferError(
                f"GCS size or MD5 verification failed for {dataset_id} {kind}"
            )
        return {
            "uri": f"gs://{bucket}/{expected_name}",
            "size_bytes": total,
            "row_count": rows,
            "sha256": sha256.hexdigest(),
            "md5_hash": md5_hash,
            "crc32c": final_object.get("crc32c"),
            "generation": str(final_object.get("generation", "")),
        }
    finally:
        stream.close()


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


def upload(
    manifest_path: Path, receipt_path: Path, delete_worker: bool = False
) -> None:
    manifest_path = _private_cache_path(str(manifest_path))
    receipt_path = _private_cache_path(str(receipt_path))
    try:
        fd = os.open(manifest_path, os.O_RDONLY | os.O_NOFOLLOW)
        with os.fdopen(fd, "r") as stream:
            info = os.fstat(stream.fileno())
            if not stat.S_ISREG(info.st_mode) or info.st_mode & 0o077:
                raise TransferError(
                    "Private transfer manifest must be a mode-0600 regular file"
                )
            manifest = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise TransferError(
            f"Cannot read private transfer manifest: {type(error).__name__}"
        ) from error
    finally:
        try:
            manifest_path.unlink()
        except FileNotFoundError:
            pass
        except OSError as error:
            raise TransferError(
                "Could not remove private transfer manifest after reading"
            ) from error
    if not isinstance(manifest, dict):
        raise TransferError("Invalid transfer manifest")
    run_key = manifest.get("run_key")
    bucket = manifest.get("bucket")
    datasets = manifest.get("datasets")
    schemas = manifest.get("schemas")
    if not isinstance(run_key, str) or not RUN_KEY.fullmatch(run_key):
        raise TransferError("Invalid run key in transfer manifest")
    if (
        not isinstance(bucket, str)
        or not bucket
        or not isinstance(datasets, list)
        or not datasets
    ):
        raise TransferError("Invalid bucket or datasets in transfer manifest")
    expected_schemas = {"dea", "gsea"}
    if not isinstance(schemas, dict) or set(schemas) != expected_schemas:
        raise TransferError("Invalid Parquet schema signatures in transfer manifest")

    work = []
    seen = set()
    for item in datasets:
        dataset_id = item.get("dataset_id") if isinstance(item, dict) else None
        if (
            not isinstance(dataset_id, str)
            or not DATASET_ID.fullmatch(dataset_id)
            or dataset_id in seen
        ):
            raise TransferError("Invalid or duplicate dataset_id in transfer manifest")
        seen.add(dataset_id)
        for kind in ("dea", "gsea"):
            source = item.get(kind)
            if not isinstance(source, dict):
                raise TransferError(f"Missing {kind} upload session for {dataset_id}")
            work.append((dataset_id, kind, source))

    receipts = {item["dataset_id"]: {} for item in datasets}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(
                _upload_one, run_key, bucket, dataset_id, kind, source, schemas[kind]
            ): (
                dataset_id,
                kind,
            )
            for dataset_id, kind, source in work
        }
        for future in as_completed(futures):
            dataset_id, kind = futures[future]
            receipts[dataset_id][kind] = future.result()
            print(f"Uploaded {dataset_id} {kind}", flush=True)

    receipt = {
        "run_key": run_key,
        "bucket": bucket,
        "datasets": [
            {"dataset_id": dataset_id, **receipts[dataset_id]}
            for dataset_id in (item["dataset_id"] for item in datasets)
        ],
    }
    raw = json.dumps(receipt, separators=(",", ":")).encode()
    fd = os.open(
        receipt_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600
    )
    with os.fdopen(fd, "wb") as stream:
        stream.write(raw)

    remove = []
    if delete_worker:
        worker = _private_cache_path(str(Path(__file__).resolve()))
        if worker.is_symlink() or not worker.is_file():
            raise TransferError("Temporary cluster worker is not a regular cache file")
        remove.append(worker)
    failures = []
    for path in remove:
        try:
            path.unlink()
        except OSError:
            failures.append(path.name)
    if failures:
        raise TransferError(
            "Upload succeeded but temporary secret/code cleanup failed: "
            + ", ".join(failures)
        )
    print(f"Verified and uploaded {len(work)} Parquet files; receipt {receipt_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    prepare_parser = subparsers.add_parser(
        "prepare", help="create private GCS sessions on the loader machine"
    )
    prepare_parser.add_argument("--source-manifest", type=Path, required=True)
    prepare_parser.add_argument("--project", default=os.getenv("GCLOUD_PROJECT", ""))
    prepare_parser.add_argument("--bucket", default=os.getenv("CLOUD_TMP_BUCKET", ""))
    prepare_parser.add_argument("--output", type=Path, required=True)
    inventory_parser = subparsers.add_parser(
        "inventory", help="locate final products and record sizes on cluster"
    )
    inventory_parser.add_argument("--dataset-ids", required=True)
    inventory_parser.add_argument("--output", type=Path, required=True)
    cleanup_parser = subparsers.add_parser(
        "cleanup", help="remove only completed objects from a failed transfer"
    )
    cleanup_parser.add_argument("--manifest", type=Path, required=True)
    cleanup_parser.add_argument("--project", default=os.getenv("GCLOUD_PROJECT", ""))
    cleanup_parser.add_argument("--bucket", default=os.getenv("CLOUD_TMP_BUCKET", ""))
    upload_parser = subparsers.add_parser("upload", help="run inside the cluster SIF")
    upload_parser.add_argument("--manifest", type=Path, required=True)
    upload_parser.add_argument("--receipt", type=Path, required=True)
    upload_parser.add_argument("--delete-worker", action="store_true")
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.source_manifest, args.project, args.bucket, args.output)
    elif args.mode == "inventory":
        inventory(args.dataset_ids.split(","), args.output)
    elif args.mode == "cleanup":
        cleanup(args.manifest, args.project, args.bucket)
    else:
        upload(args.manifest, args.receipt, args.delete_worker)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (TransferError, ValueError, OSError) as error:
        print(f"Cluster Parquet transfer failed: {error}", file=sys.stderr)
        raise SystemExit(1)
