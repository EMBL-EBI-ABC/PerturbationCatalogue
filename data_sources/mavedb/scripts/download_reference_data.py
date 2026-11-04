#!/usr/bin/env python3
"""Download the external source bundle consumed by MaveDB analyses."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from http.client import HTTPException
from pathlib import Path
from typing import Any, Sequence
from urllib.error import URLError

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis._ensembl_homology import generate_homology
from data_sources.mavedb.analysis._reference_downloads import (
    download_dataset,
    download_pockets,
)
from data_sources.mavedb.analysis._reference_io import (
    file_record,
    nonnegative_float,
    nonnegative_int,
    positive_float,
    positive_int,
    publish_tables,
    write_json,
)
from data_sources.mavedb.analysis._reference_sources import (
    DEFAULT_CACHE,
    DEFAULT_ENSEMBL,
    DEFAULT_ID_MAP,
    DEFAULT_METADATA,
    DEFAULT_POCKETS_URL,
    DEFAULT_SOURCE_ROOT,
    SOURCE_MANIFEST,
    SOURCE_SCHEMA_VERSION,
    SOURCE_URLS,
    load_id_map,
    reference_directory,
)

LOGGER = logging.getLogger("download_reference_data")
ID_MAP_FILENAME = "mappings/mavedb_uniprot_to_ensembl.tsv"
HOMOLOGY_FILENAME = "ensembl/homology.tsv"
METADATA_FILENAME = "ensembl/metadata.tsv"
REPORT_FILENAME = "ensembl/retrieval_report.json"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_SOURCE_ROOT,
        help="Root containing all references and their portable manifest",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE,
        help="Disposable Ensembl API response cache (not required for processing)",
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--id-map", type=Path, default=DEFAULT_ID_MAP)
    parser.add_argument("--ensembl-url", default=DEFAULT_ENSEMBL)
    parser.add_argument("--pockets-url", default=DEFAULT_POCKETS_URL)
    for name, url in SOURCE_URLS.items():
        parser.add_argument("--" + name.replace("_", "-") + "-url", default=url)
    parser.add_argument("--workers", type=positive_int, default=4)
    parser.add_argument("--request-interval", type=nonnegative_float, default=0.2)
    parser.add_argument("--timeout", type=positive_float, default=60)
    parser.add_argument("--retries", type=nonnegative_int, default=3)
    parser.add_argument(
        "--skip-homology",
        action="store_true",
        help="Download Open Targets and ProtVar without Ensembl requests",
    )
    parser.add_argument("--allow-missing-homology", action="store_true")
    parser.add_argument("--allow-incomplete-homology", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--offline",
        action="store_true",
        help="Reuse verified local references and Ensembl responses without network access",
    )
    mode.add_argument("--refresh", action="store_true", help="Redownload every source")
    return parser


def _published_record(staged: Path, relative: str) -> dict[str, Any]:
    return file_record(staged) | {"path": relative}


def _relative_source(record: dict[str, Any], root: Path) -> dict[str, Any]:
    """Keep stored source paths relative to the complete reference directory."""
    result = record | {
        "path": Path(record["path"]).resolve().relative_to(root.resolve()).as_posix()
    }
    if "files" in record:
        result["files"] = [_relative_source(item, root) for item in record["files"]]
    return result


def run(args: argparse.Namespace) -> dict[str, Any]:
    mapping = load_id_map(args.id_map)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".reference-data-", dir=args.out.parent
    ) as folder:
        staged = Path(folder)
        sources: dict[str, Any] = {}
        homology_result = None
        (staged / "mappings").mkdir()
        if not args.skip_homology:
            (staged / "ensembl").mkdir()
            shutil.copyfile(args.metadata, staged / METADATA_FILENAME)
            LOGGER.info("Fetching Ensembl homologies")
            homology_result = generate_homology(
                staged / METADATA_FILENAME,
                mapping,
                args.cache_dir,
                staged / HOMOLOGY_FILENAME,
                base_url=args.ensembl_url,
                offline=args.offline,
                refresh=args.refresh,
                timeout=args.timeout,
                retries=args.retries,
                interval=args.request_interval,
                workers=args.workers,
                allow_missing=args.allow_missing_homology,
                allow_incomplete=args.allow_incomplete_homology,
            )
            # Publish the retrieval record with the reference, independent of cache lifetime.
            homology_result.pop("failure_report", None)
            write_json(staged / REPORT_FILENAME, homology_result)
        for name in SOURCE_URLS:
            url = getattr(args, name + "_url")
            LOGGER.info("Fetching Open Targets %s", name)
            sources[name] = _relative_source(
                download_dataset(
                    url,
                    reference_directory(args.out, name, url),
                    offline=args.offline,
                    refresh=args.refresh,
                    timeout=args.timeout,
                    retries=args.retries,
                ),
                args.out,
            )
        LOGGER.info("Fetching ProtVar pockets")
        sources["pockets"] = _relative_source(
            download_pockets(
                args.pockets_url,
                reference_directory(args.out, "pockets", args.pockets_url),
                offline=args.offline,
                refresh=args.refresh,
                timeout=args.timeout,
                retries=args.retries,
            ),
            args.out,
        )
        mapping.to_csv(staged / ID_MAP_FILENAME, sep="\t", index=False)
        sources["id_map"] = _published_record(staged / ID_MAP_FILENAME, ID_MAP_FILENAME)
        if not args.skip_homology:
            sources["homology"] = _published_record(
                staged / HOMOLOGY_FILENAME, HOMOLOGY_FILENAME
            )
            sources["metadata"] = _published_record(
                staged / METADATA_FILENAME, METADATA_FILENAME
            )
            sources["homology_report"] = _published_record(
                staged / REPORT_FILENAME, REPORT_FILENAME
            )
        manifest = {
            "schema_version": SOURCE_SCHEMA_VERSION,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "sources": sources,
            "homology": homology_result,
        }
        write_json(staged / SOURCE_MANIFEST, manifest)
        publish_tables(
            staged,
            args.out,
            manifest=SOURCE_MANIFEST,
            managed={
                Path(SOURCE_MANIFEST),
                Path(ID_MAP_FILENAME),
                Path(HOMOLOGY_FILENAME),
                Path(METADATA_FILENAME),
                Path(REPORT_FILENAME),
            },
        )
    LOGGER.info("Source bundle ready: %s", args.out / SOURCE_MANIFEST)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    try:
        run(args)
    except (OSError, URLError, ValueError, ImportError, HTTPException, EOFError) as exc:
        LOGGER.error("Reference download failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
