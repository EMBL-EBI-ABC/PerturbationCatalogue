"""Cached HTTP inputs for MAVE reference generation."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import tempfile
import time
from concurrent.futures import CancelledError
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from http.client import HTTPException
from pathlib import Path
from threading import Event, Lock
from typing import Any, Callable
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlsplit
from urllib.request import Request, urlopen

LOGGER = logging.getLogger(__name__)
_USER_AGENT = "PerturbationCatalogue MAVE reference generator"
RETRYABLE_HTTP_CODES = frozenset({408, 429, 500, 502, 503, 504})


class RequestRateLimiter:
    """Space HTTP attempts globally and share any server-requested cooldown."""

    def __init__(self, interval: float) -> None:
        if not math.isfinite(interval) or interval < 0:
            raise ValueError("Request interval must be finite and nonnegative")
        self.interval = interval
        self._lock = Lock()
        self._next_request_at = 0.0
        self._cooldown_until = 0.0

    def acquire(self, cancellation_event: Event | None = None) -> None:
        """Wait for a fresh slot, rechecking cooldown after every wait."""
        while True:
            if cancellation_event is not None and cancellation_event.is_set():
                raise CancelledError("Request cancelled after another query failed")
            with self._lock:
                now = time.monotonic()
                ready_at = max(self._next_request_at, self._cooldown_until)
                if ready_at <= now:
                    self._next_request_at = now + self.interval
                    return
                delay = ready_at - now
            _interruptible_sleep(min(delay, 0.25), cancellation_event)

    def defer(self, seconds: float) -> None:
        """Extend the shared no-request window without reserving future slots."""
        if seconds <= 0:
            return
        with self._lock:
            self._cooldown_until = max(self._cooldown_until, time.monotonic() + seconds)


def _interruptible_sleep(delay: float, cancellation_event: Event | None) -> None:
    deadline = time.monotonic() + max(0, delay)
    while True:
        if cancellation_event is not None and cancellation_event.is_set():
            raise CancelledError("Request cancelled after another query failed")
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return
        time.sleep(min(remaining, 0.25))


def _retry_delay(error: OSError | HTTPException, attempt: int) -> float:
    if isinstance(error, HTTPError):
        retry_after = (error.headers or {}).get("Retry-After", "")
        if retry_after:
            try:
                delay = float(retry_after)
                if math.isfinite(delay):
                    return max(0.0, delay)
            except ValueError:
                pass
            try:
                retry_at = parsedate_to_datetime(retry_after)
                if retry_at.tzinfo is None:
                    retry_at = retry_at.replace(tzinfo=timezone.utc)
                return max(0.0, (retry_at - datetime.now(timezone.utc)).total_seconds())
            except (TypeError, ValueError, OverflowError):
                pass
    return float(min(2**attempt, 30))


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


def _pause(
    error: OSError | HTTPException,
    attempt: int,
    *,
    rate_limiter: RequestRateLimiter | None = None,
    cancellation_event: Event | None = None,
    url: str | None = None,
) -> None:
    if isinstance(error, HTTPError) and error.code not in RETRYABLE_HTTP_CODES:
        raise error
    delay = _retry_delay(error, attempt)
    LOGGER.warning(
        "HTTP request failed for %s; retrying in %.1f seconds: %s",
        url or getattr(error, "url", "unknown URL"),
        delay,
        error,
    )
    if rate_limiter is not None and _shared_cooldown(error):
        rate_limiter.defer(delay)
    else:
        _interruptible_sleep(delay, cancellation_event)


def _shared_cooldown(error: OSError | HTTPException) -> bool:
    """Only rate limits or explicit server instructions pause other requests."""
    return isinstance(error, HTTPError) and (
        error.code == 429 or bool((error.headers or {}).get("Retry-After"))
    )


def request_bytes(
    url: str,
    timeout: float,
    retries: int,
    *,
    accept: str = "application/json",
    rate_limiter: RequestRateLimiter | None = None,
    cancellation_event: Event | None = None,
) -> bytes:
    limiter = rate_limiter or RequestRateLimiter(0)
    for attempt in range(retries + 1):
        try:
            limiter.acquire(cancellation_event)
            headers = {"User-Agent": _USER_AGENT, "Accept": accept}
            if accept == "application/json":
                headers["Content-Type"] = "application/json"
            request = Request(url, headers=headers)
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except (URLError, OSError, TimeoutError, HTTPException) as exc:
            if attempt == retries:
                if _shared_cooldown(exc):
                    limiter.defer(_retry_delay(exc, attempt))
                raise
            _pause(
                exc,
                attempt,
                rate_limiter=limiter,
                cancellation_event=cancellation_event,
                url=url,
            )
    raise RuntimeError("Unreachable retry state")


def json_cache_path(cache: Path, url: str) -> Path:
    return cache / "ensembl" / f"{hashlib.sha256(url.encode()).hexdigest()}.json"


def _error_body(error: HTTPError) -> str:
    """Keep bounded provider diagnostics and release the failed connection."""
    try:
        return error.read(8192).decode("utf-8", errors="replace")
    except (OSError, HTTPException):
        return ""
    finally:
        error.close()


def _contextual_http_error(
    url: str,
    status: int,
    reason: str,
    body: str = "",
    headers: Any = None,
) -> HTTPError:
    detail = body
    if body:
        try:
            parsed = json.loads(body)
            if isinstance(parsed, dict) and parsed.get("error") is not None:
                detail = parsed["error"]
        except json.JSONDecodeError:
            pass
    if not isinstance(detail, str):
        detail = json.dumps(detail)
    detail = " ".join(detail.split())[:2000]
    message = f"{reason}: {detail}" if detail else reason
    return HTTPError(url, status, f"{message} [URL: {url}]", headers or {}, None)


def cached_json(
    url: str,
    cache: Path,
    *,
    offline: bool,
    refresh: bool,
    timeout: float,
    retries: int,
    interval: float,
    rate_limiter: RequestRateLimiter | None = None,
    cancellation_event: Event | None = None,
    on_cache_hit: Callable[[], None] | None = None,
) -> dict[str, Any]:
    path = json_cache_path(cache, url)
    if path.is_file() and not refresh:
        entry = json.loads(path.read_text())
        if entry.get("url") == url and entry.get("http_status") in (400, 404):
            if on_cache_hit is not None:
                on_cache_hit()
            raise _contextual_http_error(
                url,
                entry["http_status"],
                entry.get("reason", "Unavailable Ensembl gene"),
                entry.get("error_body", ""),
            )
        if entry.get("url") != url or not isinstance(entry.get("payload"), dict):
            raise ValueError(f"Invalid Ensembl cache entry: {path}")
        if on_cache_hit is not None:
            on_cache_hit()
        return entry["payload"]
    if offline:
        raise FileNotFoundError(
            f"No cached Ensembl response for {url}; rerun without --offline"
        )
    try:
        payload = json.loads(
            request_bytes(
                url,
                timeout,
                retries,
                rate_limiter=rate_limiter or RequestRateLimiter(interval),
                cancellation_event=cancellation_event,
            )
        )
    except HTTPError as exc:
        body = _error_body(exc)
        if exc.code in (400, 404):
            write_json(
                path,
                {
                    "url": url,
                    "http_status": exc.code,
                    "reason": str(exc.reason),
                    "error_body": body,
                },
            )
        raise _contextual_http_error(
            url, exc.code, str(exc.reason), body, exc.headers
        ) from exc
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
            _pause(exc, attempt, url=url)
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
