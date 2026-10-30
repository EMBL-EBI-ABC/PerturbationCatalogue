"""Cached HTTP inputs for MAVE reference generation."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
import time
from html.parser import HTMLParser
from http.client import HTTPException
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlsplit
from urllib.request import Request, urlopen

LOGGER = logging.getLogger(__name__)
_USER_AGENT = "PerturbationCatalogue MAVE reference generator"


def write_json(path: Path, value: Any) -> None:
    """Publish a JSON cache entry only after the complete write succeeds."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.close()
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)


def _pause(error: OSError | HTTPException, attempt: int) -> None:
    if isinstance(error, HTTPError) and error.code not in (429, 500, 502, 503, 504):
        raise error
    retry_after = (
        (error.headers or {}).get("Retry-After", "")
        if isinstance(error, HTTPError)
        else ""
    )
    try:
        delay = float(retry_after)
    except ValueError:
        delay = 2**attempt
    LOGGER.warning(
        "HTTP request failed; retrying in %.1f seconds: %s", min(delay, 30), error
    )
    time.sleep(max(0, min(delay, 30)))


def request_bytes(
    url: str, timeout: float, retries: int, *, accept: str = "application/json"
) -> bytes:
    for attempt in range(retries + 1):
        try:
            headers = {"User-Agent": _USER_AGENT, "Accept": accept}
            if accept == "application/json":
                headers["Content-Type"] = "application/json"
            request = Request(url, headers=headers)
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except (URLError, OSError, TimeoutError, HTTPException) as exc:
            if attempt == retries:
                raise
            _pause(exc, attempt)
    raise RuntimeError("Unreachable retry state")


def json_cache_path(cache: Path, url: str) -> Path:
    return cache / "ensembl" / f"{hashlib.sha256(url.encode()).hexdigest()}.json"


def cached_json(
    url: str,
    cache: Path,
    *,
    offline: bool,
    refresh: bool,
    timeout: float,
    retries: int,
    interval: float,
) -> dict[str, Any]:
    path = json_cache_path(cache, url)
    if path.is_file() and not refresh:
        entry = json.loads(path.read_text())
        if entry.get("url") == url and entry.get("http_status") in (400, 404):
            raise HTTPError(
                url,
                entry["http_status"],
                entry.get("reason", "Unavailable Ensembl gene"),
                {},
                None,
            )
        if entry.get("url") != url or not isinstance(entry.get("payload"), dict):
            raise ValueError(f"Invalid Ensembl cache entry: {path}")
        return entry["payload"]
    if offline:
        raise FileNotFoundError(
            f"No cached Ensembl response for {url}; rerun without --offline"
        )
    time.sleep(interval)
    try:
        payload = json.loads(request_bytes(url, timeout, retries))
    except HTTPError as exc:
        if exc.code in (400, 404):
            write_json(
                path, {"url": url, "http_status": exc.code, "reason": str(exc.reason)}
            )
        raise
    if not isinstance(payload, dict) or not isinstance(payload.get("data"), list):
        raise ValueError(f"Invalid Ensembl homology response from {url}")
    write_json(path, {"url": url, "payload": payload})
    return payload


class _Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            self.links.extend(value for key, value in attrs if key == "href" and value)


def parquet_links(base: str, listing: bytes) -> dict[str, str]:
    """Select direct Parquet files from an EBI directory index."""
    parser = _Links()
    parser.feed(listing.decode("utf-8"))
    root = urlsplit(base)
    files = {}
    for href in parser.links:
        url = urljoin(base, href)
        parsed = urlsplit(url)
        name = unquote(parsed.path.rsplit("/", 1)[-1])
        if (parsed.scheme, parsed.netloc) != (root.scheme, root.netloc):
            continue
        if parsed.path.rsplit("/", 1)[0] != root.path.rstrip("/"):
            continue
        if name.endswith(".parquet") and Path(name).name == name:
            if name in files and files[name] != url:
                raise ValueError(f"Duplicate Parquet filename in {base}: {name}")
            files[name] = url
    if not files:
        raise ValueError(f"No Parquet files found in directory listing: {base}")
    return dict(sorted(files.items()))


def _parquet_valid(path: Path) -> bool:
    import pyarrow as pa
    import pyarrow.parquet as pq

    try:
        pq.ParquetFile(path)
        return True
    except (OSError, pa.ArrowInvalid):
        return False


def _download(
    url: str, destination: Path, timeout: float, retries: int
) -> dict[str, Any]:
    for attempt in range(retries + 1):
        temporary: Path | None = None
        try:
            digest = hashlib.sha256()
            count = 0
            request = Request(url, headers={"User-Agent": _USER_AGENT})
            with urlopen(request, timeout=timeout) as response:
                length = response.headers.get("Content-Length")
                with tempfile.NamedTemporaryFile(
                    dir=destination.parent, suffix=".part", delete=False
                ) as handle:
                    temporary = Path(handle.name)
                    while block := response.read(1024 * 1024):
                        handle.write(block)
                        digest.update(block)
                        count += len(block)
            if length is not None and count != int(length):
                raise OSError(
                    f"Incomplete download of {url}: received {count}/{length} bytes"
                )
            if not _parquet_valid(temporary):
                raise ValueError(f"Downloaded file is not valid Parquet: {url}")
            os.replace(temporary, destination)
            return {"bytes": count, "sha256": digest.hexdigest()}
        except (URLError, OSError, TimeoutError, HTTPException) as exc:
            if attempt == retries:
                raise
            _pause(exc, attempt)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    raise RuntimeError("Unreachable retry state")


def download_dataset(
    url: str,
    cache: Path,
    *,
    offline: bool = False,
    refresh: bool = False,
    timeout: float = 60,
    retries: int = 3,
) -> list[Path]:
    """Download a complete partitioned dataset, resuming verified local files.

    Cache locations are isolated by source URL. The complete expected file list
    is saved before downloading, so an offline run cannot use an incomplete cache.
    """
    url = url.rstrip("/") + "/"
    directory = cache / "parquet" / hashlib.sha256(url.encode()).hexdigest()[:16]
    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory / "download.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    if offline:
        if manifest.get("url") != url or not manifest.get("expected"):
            raise FileNotFoundError(
                f"No complete download cache for {url}; supply local Parquet inputs or rerun online"
            )
        links = {name: urljoin(url, name) for name in manifest["expected"]}
    else:
        links = parquet_links(
            url, request_bytes(url, timeout, retries, accept="text/html")
        )
    if any(Path(name).name != name or not name.endswith(".parquet") for name in links):
        raise ValueError(f"Invalid filename in download cache: {manifest_path}")
    old_files = manifest.get("files", {}) if manifest.get("url") == url else {}
    manifest = {"url": url, "expected": list(links), "files": old_files}
    if not offline:
        write_json(manifest_path, manifest)
    paths = []
    for index, (name, file_url) in enumerate(links.items(), 1):
        path = directory / name
        info = old_files.get(name, {})
        valid = (
            path.is_file()
            and info.get("bytes") == path.stat().st_size
            and _parquet_valid(path)
        )
        if refresh or not valid:
            if offline:
                raise FileNotFoundError(
                    f"Missing or invalid cached Parquet file: {path}"
                )
            LOGGER.info("Downloading Parquet %d/%d: %s", index, len(links), file_url)
            manifest["files"][name] = _download(file_url, path, timeout, retries)
            write_json(manifest_path, manifest)
        paths.append(path)
    return paths
