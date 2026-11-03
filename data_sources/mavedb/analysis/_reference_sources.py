"""Shared source locations and the cached ProtVar pocket download."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
from http.client import HTTPException
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

from data_sources.mavedb.analysis._homology_downloads import _pause, write_json
from data_sources.mavedb.analysis._homology_references import (
    DEFAULT_CODING,
    DEFAULT_EVIDENCE,
    MAVEDB_ROOT,
)

DEFAULT_SOURCE_ROOT = MAVEDB_ROOT / "output" / "reference_data"
SOURCE_MANIFEST = "downloads.json"
DEFAULT_SOURCE_MANIFEST = DEFAULT_SOURCE_ROOT / SOURCE_MANIFEST
OT_ROOT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output"
SOURCE_URLS = {
    "coding_variants": DEFAULT_CODING,
    "evidence": DEFAULT_EVIDENCE,
    "credible_sets": f"{OT_ROOT}/credible_set/",
    "pharmacogenomics": f"{OT_ROOT}/pharmacogenomics/",
    "drugs": f"{OT_ROOT}/drug_molecule/",
}
POCKET_FILENAME = "2024.05.28_pockets.tsv"
DEFAULT_POCKETS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/pockets/"
    f"{POCKET_FILENAME}.gz"
)
POCKET_COLUMNS = {
    "struct_id",
    "pocket_id",
    "pocket_resid",
    "pocket_score_combined_scaled",
    "pocket_rad_gyration",
    "pocket_plddt_mean",
}


def file_record(path: Path) -> dict[str, Any]:
    """Hash a file without collecting it in memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _matches(path: Path, record: dict[str, Any]) -> bool:
    if (
        not isinstance(record, dict)
        or not path.is_file()
        or path.stat().st_size != record.get("bytes")
    ):
        return False
    return file_record(path)["sha256"] == record.get("sha256")


def _unpack_pockets(compressed: Path, destination: Path) -> None:
    """Validate the complete gzip stream and TSV header before replacing a file."""
    temporary: Path | None = None
    try:
        with gzip.open(compressed, "rb") as source, tempfile.NamedTemporaryFile(
            dir=destination.parent, suffix=".part", delete=False
        ) as target:
            temporary = Path(target.name)
            while block := source.read(1024 * 1024):
                target.write(block)
        with temporary.open(newline="", encoding="utf-8") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
        if missing := POCKET_COLUMNS - set(header):
            raise ValueError(
                f"ProtVar pocket TSV is missing columns: {sorted(missing)}"
            )
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def download_pockets(
    url: str,
    cache: Path,
    *,
    offline: bool = False,
    refresh: bool = False,
    timeout: float = 60,
    retries: int = 3,
) -> tuple[Path, dict[str, Any]]:
    """Download, verify and unpack the dated ProtVar gzip, with offline reuse."""
    directory = cache / "pockets" / hashlib.sha256(url.encode()).hexdigest()[:16]
    directory.mkdir(parents=True, exist_ok=True)
    compressed = directory / f"{POCKET_FILENAME}.gz"
    destination = directory / POCKET_FILENAME
    manifest_path = directory / "download.json"
    manifest = (
        json.loads(manifest_path.read_text())
        if manifest_path.is_file() and not refresh
        else {}
    )
    if not isinstance(manifest, dict):
        raise ValueError(f"Invalid ProtVar pocket cache manifest: {manifest_path}")
    same_source = manifest.get("url") == url
    if (
        not refresh
        and same_source
        and _matches(compressed, manifest.get("compressed", {}))
    ):
        if not _matches(destination, manifest.get("uncompressed", {})):
            _unpack_pockets(compressed, destination)
            manifest["uncompressed"] = file_record(destination)
            write_json(manifest_path, manifest)
        return destination, manifest
    if offline:
        raise FileNotFoundError(
            f"No verified ProtVar pocket cache for {url}; rerun online"
        )
    for attempt in range(retries + 1):
        temporary: Path | None = None
        try:
            request = Request(
                url, headers={"User-Agent": "PerturbationCatalogue-MaveDB"}
            )
            with urlopen(
                request, timeout=timeout
            ) as response, tempfile.NamedTemporaryFile(
                dir=directory, suffix=".gz.part", delete=False
            ) as target:
                temporary = Path(target.name)
                expected = response.headers.get("Content-Length")
                size = 0
                while block := response.read(1024 * 1024):
                    target.write(block)
                    size += len(block)
            if expected is not None and size != int(expected):
                raise OSError(f"Incomplete pocket download: {size}/{expected} bytes")
            _unpack_pockets(temporary, destination)
            os.replace(temporary, compressed)
            manifest = {
                "url": url,
                "compressed": file_record(compressed),
                "uncompressed": file_record(destination),
            }
            write_json(manifest_path, manifest)
            return destination, manifest
        except (URLError, OSError, TimeoutError, EOFError, HTTPException) as exc:
            if attempt == retries:
                raise
            _pause(exc, attempt, url=url)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    raise RuntimeError("Unreachable pocket download state")


def load_sources(path: Path) -> dict[str, Any]:
    """Read a completed download manifest; local paths are stored absolutely."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise ValueError(f"Unsupported source manifest: {path}")
    if not isinstance(value.get("sources"), dict):
        raise ValueError(f"Source manifest is missing its sources: {path}")
    for name, record in value["sources"].items():
        if (
            not isinstance(record, dict)
            or not isinstance(record.get("path"), str)
            or not Path(record["path"]).is_absolute()
        ):
            raise ValueError(
                f"Source manifest entry {name} requires an absolute path: {path}"
            )
    return value
