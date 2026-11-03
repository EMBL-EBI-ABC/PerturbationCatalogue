#!/usr/bin/env python3
"""Download all Open Targets, Ensembl and ProtVar inputs for MaveDB analyses."""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime, timezone
from http.client import HTTPException
from pathlib import Path
from typing import Any, Sequence
from urllib.error import URLError

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis._homology_downloads import (
    download_dataset,
    write_json,
)
from data_sources.mavedb.analysis._homology_references import (
    DEFAULT_CACHE,
    DEFAULT_ENSEMBL,
    DEFAULT_ID_MAP,
    DEFAULT_METADATA,
    DEFAULT_REFERENCES,
    homology_queries,
    load_id_map,
)
from data_sources.mavedb.analysis._reference_sources import (
    DEFAULT_POCKETS_URL,
    DEFAULT_SOURCE_ROOT,
    SOURCE_MANIFEST,
    SOURCE_URLS,
    download_pockets,
    file_record,
)
from data_sources.mavedb.scripts import prepare_homology_references as homology

LOGGER = logging.getLogger("download_reference_data")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_SOURCE_ROOT,
        help="Destination for downloads.json (raw files stay in the shared cache)",
    )
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--metadata", "--meta", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--id-map", type=Path, default=DEFAULT_ID_MAP)
    parser.add_argument("--homology-out", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--ensembl-url", default=DEFAULT_ENSEMBL)
    parser.add_argument("--pockets-url", default=DEFAULT_POCKETS_URL)
    for name, url in SOURCE_URLS.items():
        parser.add_argument("--" + name.replace("_", "-") + "-url", default=url)
    parser.add_argument("--workers", type=homology._positive_int, default=4)
    parser.add_argument("--request-interval", type=homology._nonnegative, default=0.2)
    parser.add_argument("--timeout", type=homology._positive, default=60)
    parser.add_argument("--retries", type=homology._nonnegative_int, default=3)
    parser.add_argument(
        "--skip-homology",
        action="store_true",
        help="Download Open Targets and ProtVar only; no metadata/Ensembl requests",
    )
    parser.add_argument("--allow-missing-homology", action="store_true")
    parser.add_argument("--allow-incomplete-homology", action="store_true")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--offline",
        action="store_true",
        help="Require complete local download/API caches; never access the network",
    )
    modes.add_argument(
        "--refresh",
        action="store_true",
        help="Redownload every selected source, including Ensembl responses",
    )
    return parser


def run(args: argparse.Namespace) -> dict[str, Any]:
    # Detect absent/ambiguous inputs before starting any large source downloads.
    mapping = load_id_map(args.id_map, deduplicate=False)
    if not args.skip_homology:
        homology_queries(args.metadata, mapping)
    manifest_path = args.out / SOURCE_MANIFEST
    sources = {}
    for name in SOURCE_URLS:
        url = getattr(args, name + "_url")
        LOGGER.info("Fetching Open Targets %s", name)
        paths = download_dataset(
            url,
            args.cache_dir,
            offline=args.offline,
            refresh=args.refresh,
            timeout=args.timeout,
            retries=args.retries,
        )
        if not paths or len({path.parent.resolve() for path in paths}) != 1:
            raise ValueError(f"Expected one cache directory for {name}")
        sources[name] = {
            "url": url,
            "path": str(paths[0].parent.resolve()),
            "files": [str(path.resolve()) for path in paths],
            "download_manifest": str((paths[0].parent / "download.json").resolve()),
        }
    LOGGER.info("Fetching and unpacking ProtVar pockets")
    pockets, pocket_info = download_pockets(
        args.pockets_url,
        args.cache_dir,
        offline=args.offline,
        refresh=args.refresh,
        timeout=args.timeout,
        retries=args.retries,
    )
    sources["pockets"] = {"url": args.pockets_url, **pocket_info["uncompressed"]}
    sources["id_map"] = file_record(args.id_map)
    homology_result = None
    if not args.skip_homology:
        LOGGER.info("Fetching Ensembl homologies")
        options = argparse.Namespace(
            **(vars(args) | {"only": "homology", "out": args.homology_out})
        )
        homology_result = homology.run(options)
        sources["homology"] = {
            "path": str((args.homology_out / "homology.tsv").resolve()),
            "references": str((args.homology_out / "references.json").resolve()),
        }
        sources["metadata"] = file_record(args.metadata)
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "sources": sources,
        "homology": homology_result,
    }
    write_json(manifest_path, manifest)
    LOGGER.info("All requested source downloads are ready: %s", manifest_path)
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
