"""Build the Ensembl homology reference table used by MaveDB analyses."""

from __future__ import annotations

import csv
import json
import logging
import math
import shutil
import tempfile
import threading
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from http.client import HTTPException
from pathlib import Path
from typing import Any, TextIO
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode

import pandas as pd

from data_sources.mavedb.analysis._homology_inputs import read_metadata
from data_sources.mavedb.analysis._reference_downloads import (
    RETRYABLE_HTTP_CODES,
    RequestRateLimiter,
    cached_json,
)
from data_sources.mavedb.analysis._reference_io import write_json
from data_sources.mavedb.analysis._reference_sources import DEFAULT_ENSEMBL

LOGGER = logging.getLogger(__name__)
PROGRESS_INTERVAL_SECONDS = 5.0


class HomologyFetchError(ValueError):
    """Some gene requests exhausted their retries; the cache can be resumed."""


HOMOLOGY_COLUMNS = (
    "taxonomy_level",
    "target_species",
    "source_species",
    "type",
    "target_id",
    "source_id",
    "target_align_seq",
    "source_align_seq",
    "target_perc_pos",
    "target_perc_id",
)


def _ensembl_url(base: str, species: str, gene: str) -> str:
    return (
        f"{base.rstrip('/')}/homology/id/{quote(species, safe='')}/{quote(gene, safe='')}?"
        + urlencode({"aligned": 1, "sequence": "protein", "format": "full"})
    )


def _homology_queries(metadata: Path, mapping: pd.DataFrame) -> pd.DataFrame:
    meta = read_metadata(metadata)
    incomplete = meta["identifier_uniprot"].isna() | meta["organism"].isna()
    if incomplete.any():
        LOGGER.warning(
            "%d metadata score sets lack an unambiguous UniProt ID or organism and are excluded from homology requests",
            int(incomplete.sum()),
        )
    pairs = (
        mapping.merge(
            meta[["identifier_uniprot", "organism"]].drop_duplicates(),
            left_on="From",
            right_on="identifier_uniprot",
            how="inner",
        )[["To", "organism"]]
        .dropna()
        .drop_duplicates()
    )
    pairs["organism"] = (
        pairs["organism"]
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"\s+", "_", regex=True)
    )
    pairs = pairs.loc[pairs["organism"].ne("")]
    if pairs.empty:
        raise ValueError(
            "No mapped Ensembl genes have an unambiguous organism in metadata"
        )
    missing = set(meta["identifier_uniprot"].dropna()) - set(mapping["From"])
    if missing:
        LOGGER.warning(
            "%d metadata UniProt IDs are absent from the ID map; supply --id-map for wider coverage",
            len(missing),
        )
    return pairs.rename(columns={"To": "gene", "organism": "species"}).reset_index(
        drop=True
    )


def _flatten_homology(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten Ensembl's nested source/target objects without swapping endpoints."""
    if not isinstance(payload.get("data"), list):
        raise ValueError("Ensembl homology payload must contain a data array")
    rows = []
    for gene in payload["data"]:
        if not isinstance(gene, dict) or not isinstance(gene.get("homologies"), list):
            raise ValueError("Ensembl data entry must contain a homologies array")
        for homology in gene["homologies"]:
            if not isinstance(homology, dict):
                raise ValueError("Invalid Ensembl homology entry")
            row = {
                key: value
                for key, value in homology.items()
                if key not in ("source", "target")
            }
            for side in ("source", "target"):
                endpoint = homology.get(side)
                if not isinstance(endpoint, dict):
                    raise ValueError(f"Ensembl homology is missing its {side} endpoint")
                row.update({f"{side}_{key}": value for key, value in endpoint.items()})
            # Ensembl omits this annotation for valid alternate-allele links.
            # Keep a nullable column instead of inventing a taxonomy label.
            row.setdefault("taxonomy_level", None)
            missing = [
                column
                for column in HOMOLOGY_COLUMNS
                if column != "taxonomy_level" and row.get(column) is None
            ]
            context = (
                f"{row.get('source_id') or '?'} -> {row.get('target_id') or '?'} "
                f"({row.get('type') or 'unknown type'})"
            )
            if missing:
                raise ValueError(
                    f"Ensembl homology {context} is missing required alignment fields: "
                    + ", ".join(missing)
                )
            target, source = row["target_align_seq"], row["source_align_seq"]
            if (
                not isinstance(target, str)
                or not isinstance(source, str)
                or len(target) != len(source)
            ):
                raise ValueError(
                    f"Ensembl alignment sequences for {context} must be equal-length strings"
                )
            rows.append(row)
    return rows


def generate_homology(
    metadata: Path,
    mapping: pd.DataFrame,
    cache: Path,
    destination: Path,
    *,
    base_url: str = DEFAULT_ENSEMBL,
    offline: bool = False,
    refresh: bool = False,
    timeout: float = 60,
    retries: int = 3,
    interval: float = 0.2,
    workers: int = 4,
    allow_missing: bool = False,
    allow_incomplete: bool = False,
) -> dict[str, Any]:
    if workers <= 0:
        raise ValueError("Homology workers must be a positive integer")
    queries = _homology_queries(metadata, mapping)
    # Spool completed queries to disk so ordering never blocks new requests.
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HOMOLOGY_COLUMNS)
        with tempfile.TemporaryDirectory(
            prefix=".homology-results-", dir=destination.parent
        ) as directory:
            count, missing, failed = _write_homology_queries(
                queries,
                handle,
                base_url,
                cache,
                offline,
                refresh,
                timeout,
                retries,
                interval,
                workers,
                RequestRateLimiter(interval),
                allow_missing,
                Path(directory),
            )
    report_path = (cache / "homology_failures.json").resolve()
    write_json(
        report_path,
        {
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "api": base_url,
            "queries": len(queries),
            "successful_queries": len(queries) - len(failed) - len(missing),
            "complete": not (failed or missing),
            "failed": failed,
            "skipped": missing,
        },
    )
    if failed and not allow_incomplete:
        raise HomologyFetchError(
            f"{len(failed)}/{len(queries)} Ensembl gene requests failed after retries. "
            f"Successful responses are cached; rerun to retry the remaining genes. "
            f"Failure details: {report_path}"
        )
    if len(missing) + len(failed) == len(queries):
        raise ValueError(
            "Every Ensembl query was unavailable; check the gene IDs and species"
        )
    if failed:
        LOGGER.warning(
            "Publishing explicitly permitted incomplete homology: %d failed genes; see %s",
            len(failed),
            report_path,
        )
    return {
        "queries": len(queries),
        "rows": count,
        "skipped": missing,
        "failed": failed,
        "complete": not (failed or missing),
        "failure_report": str(report_path),
        "api": base_url,
    }


def _write_homology_queries(
    queries: pd.DataFrame,
    handle: TextIO,
    base_url: str,
    cache: Path,
    offline: bool,
    refresh: bool,
    timeout: float,
    retries: int,
    interval: float,
    workers: int,
    rate_limiter: RequestRateLimiter,
    allow_missing: bool,
    spool: Path,
) -> tuple[int, list[dict[str, Any]], list[dict[str, Any]]]:
    count = 0
    total = len(queries)
    started = time.monotonic()
    progress_lock = threading.Lock()
    progress = {"completed": 0, "cache_hits": 0, "in_flight": 0, "failed": 0}
    log_every = max(1, total // 20)
    last_progress_log = started
    cancellation_event = threading.Event()
    query_rows = iter(queries.itertuples(index=False))
    next_index = 1
    pending: dict[Future[tuple[Path, int]], tuple[int, Any]] = {}
    results: dict[int, tuple[Path, int]] = {}
    skipped_by_index: dict[int, dict[str, Any]] = {}
    failed_by_index: dict[int, dict[str, Any]] = {}

    LOGGER.info(
        "Homology query plan: %d genes, %d workers, %.3g s global request interval",
        total,
        workers,
        interval,
    )

    def cache_hit() -> None:
        with progress_lock:
            progress["cache_hits"] += 1

    def report_progress(*, force: bool = False) -> None:
        nonlocal last_progress_log
        now = time.monotonic()
        with progress_lock:
            if not force and now - last_progress_log < PROGRESS_INTERVAL_SECONDS:
                return
            last_progress_log = now
            completed_count = progress["completed"]
            cache_hits = progress["cache_hits"]
            in_flight = progress["in_flight"]
            failed_count = progress["failed"]
        LOGGER.info(
            "Homology progress: %d/%d complete (%d cache hits, %d in flight; %.1f s elapsed; %d failed)",
            completed_count,
            total,
            cache_hits,
            in_flight,
            now - started,
            failed_count,
        )

    def completed(_future: Future[tuple[Path, int]]) -> None:
        with progress_lock:
            progress["completed"] += 1
            progress["in_flight"] -= 1
            completed_count = progress["completed"]
        report_progress(
            force=completed_count % log_every == 0 or completed_count == total
        )

    def fetch_and_spool(index: int, row: Any) -> tuple[Path, int]:
        payload = cached_json(
            _ensembl_url(base_url, row.species, row.gene),
            cache,
            offline=offline,
            refresh=refresh,
            timeout=timeout,
            retries=retries,
            interval=interval,
            rate_limiter=rate_limiter,
            cancellation_event=cancellation_event,
            on_cache_hit=cache_hit,
        )
        path = spool / f"{index}.tsv"
        rows = 0
        with path.open("w", newline="") as part:
            writer = csv.writer(part, delimiter="\t", lineterminator="\n")
            for homology in _flatten_homology(payload):
                writer.writerow(
                    [_tsv_value(homology[column]) for column in HOMOLOGY_COLUMNS]
                )
                rows += 1
        return path, rows

    executor = ThreadPoolExecutor(max_workers=workers)

    def submit_available() -> None:
        nonlocal next_index
        while len(pending) < workers:
            try:
                row = next(query_rows)
            except StopIteration:
                return
            with progress_lock:
                progress["in_flight"] += 1
            try:
                future = executor.submit(
                    fetch_and_spool,
                    next_index,
                    row,
                )
            except BaseException:
                with progress_lock:
                    progress["in_flight"] -= 1
                raise
            future.add_done_callback(completed)
            pending[future] = (next_index, row)
            next_index += 1

    try:
        submit_available()
        while pending:
            done, _ = wait(
                pending,
                timeout=PROGRESS_INTERVAL_SECONDS,
                return_when=FIRST_COMPLETED,
            )
            if not done:
                report_progress()
                continue

            for future in done:
                index, row = pending.pop(future)
                try:
                    results[index] = future.result()
                except Exception as exc:
                    if (
                        isinstance(exc, HTTPError)
                        and allow_missing
                        and exc.code in (400, 404)
                    ):
                        LOGGER.warning(
                            "Skipping unavailable Ensembl gene %s (%s): HTTP %d",
                            row.gene,
                            row.species,
                            exc.code,
                        )
                        skipped_by_index[index] = {
                            "gene": row.gene,
                            "species": row.species,
                            "http_status": exc.code,
                            "url": _ensembl_url(base_url, row.species, row.gene),
                            "error": str(exc),
                        }
                    elif _transient_homology_error(exc):
                        failed_by_index[index] = {
                            "gene": row.gene,
                            "species": row.species,
                            "url": _ensembl_url(base_url, row.species, row.gene),
                            "http_status": (
                                exc.code if isinstance(exc, HTTPError) else None
                            ),
                            "error": str(exc),
                        }
                        with progress_lock:
                            progress["failed"] += 1
                        LOGGER.warning(
                            "Ensembl gene request exhausted retries: %s (%s): %s; continuing other genes",
                            row.gene,
                            row.species,
                            exc,
                        )
                    else:
                        if isinstance(exc, HTTPError) and exc.code in (400, 404):
                            raise HTTPError(
                                exc.url,
                                exc.code,
                                f"{exc.reason}; gene={row.gene}, species={row.species}. "
                                "Use --allow-missing-homology only to explicitly permit skipping unavailable entries.",
                                exc.headers or {},
                                None,
                            ) from exc
                        raise
            submit_available()
        for index in sorted(results):
            path, rows = results[index]
            with path.open(newline="") as part:
                shutil.copyfileobj(part, handle)
            path.unlink()
            count += rows
        report_progress(force=True)
    except BaseException:
        cancellation_event.set()
        for future in pending:
            future.cancel()
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)

    return (
        count,
        [skipped_by_index[i] for i in sorted(skipped_by_index)],
        [failed_by_index[i] for i in sorted(failed_by_index)],
    )


def _transient_homology_error(error: Exception) -> bool:
    if isinstance(error, HTTPError):
        return error.code in RETRYABLE_HTTP_CODES
    return isinstance(error, (URLError, TimeoutError, ConnectionError, HTTPException))


def _tsv_value(value: Any) -> Any:
    if (
        value is None
        or value is pd.NA
        or isinstance(value, float)
        and math.isnan(value)
    ):
        return "NA"
    if isinstance(value, (dict, list)):
        return json.dumps(value, separators=(",", ":"), sort_keys=True)
    return value
