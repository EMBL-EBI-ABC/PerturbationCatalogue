#!/usr/bin/env python3
"""Generate local Ensembl homology and Open Targets ClinVar reference TSVs."""

from __future__ import annotations

import argparse
import logging
import math
import sys
import tempfile
from http.client import HTTPException
from pathlib import Path
from typing import Sequence
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
    DEFAULT_CODING,
    DEFAULT_ENSEMBL,
    DEFAULT_EVIDENCE,
    DEFAULT_ID_MAP,
    DEFAULT_METADATA,
    DEFAULT_REFERENCES,
    generate_clinvar,
    generate_homology,
    load_id_map,
    publish_references,
    reference_manifest,
)

LOGGER = logging.getLogger("prepare_homology_references")


def _positive(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("Must be positive")
    return number


def _positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return number


def _nonnegative(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("Must be nonnegative")
    return number


def _nonnegative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("Must be a nonnegative integer")
    return number


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=("all", "homology", "clinvar"), default="all")
    parser.add_argument("--metadata", "--meta", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--id-map", type=Path, default=DEFAULT_ID_MAP)
    parser.add_argument("--out", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--evidence",
        type=Path,
        help="Existing evidence_eva Parquet file/directory; otherwise download it",
    )
    parser.add_argument(
        "--coding-variants",
        type=Path,
        help="Existing coding_variant_table Parquet file/directory; otherwise download it",
    )
    parser.add_argument("--evidence-url", default=DEFAULT_EVIDENCE)
    parser.add_argument("--coding-url", default=DEFAULT_CODING)
    parser.add_argument("--ensembl-url", default=DEFAULT_ENSEMBL)
    parser.add_argument("--batch-size", type=_positive_int, default=65_536)
    parser.add_argument(
        "--workers",
        type=_positive_int,
        default=4,
        help="Maximum concurrent Ensembl gene queries (default: 4)",
    )
    parser.add_argument("--timeout", type=_positive, default=60)
    parser.add_argument(
        "--retries",
        type=_nonnegative_int,
        default=3,
        help="Retries after transient HTTP failures",
    )
    parser.add_argument(
        "--request-interval",
        type=_nonnegative,
        default=0.2,
        help="Minimum global seconds between uncached Ensembl HTTP attempts (default: 0.2)",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--offline",
        action="store_true",
        help="Require local Parquet files and cached Ensembl responses; never use the network",
    )
    mode.add_argument(
        "--refresh",
        action="store_true",
        help="Replace cached API responses and downloads",
    )
    parser.add_argument(
        "--allow-missing-homology",
        action="store_true",
        help="Skip HTTP 400/404 Ensembl genes and record every skip in references.json",
    )
    parser.add_argument(
        "--allow-incomplete-homology",
        action="store_true",
        help="Explicitly publish partial homology after exhausted transient failures; record all failed genes",
    )
    return parser


def run(args: argparse.Namespace) -> dict:
    mapping = load_id_map(args.id_map)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results = {}
    with tempfile.TemporaryDirectory(
        prefix=".homology-references-", dir=args.out.parent
    ) as directory:
        staged = Path(directory)
        mapping.to_csv(staged / "mavedb_uniprot_to_ensembl.tsv", sep="\t", index=False)
        if args.only in ("all", "homology"):
            results["homology"] = generate_homology(
                args.metadata,
                mapping,
                args.cache_dir,
                staged / "homology.tsv",
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
        if args.only in ("all", "clinvar"):
            sources = {}
            for name, local, url in (
                ("evidence", args.evidence, args.evidence_url),
                ("coding", args.coding_variants, args.coding_url),
            ):
                if local is not None:
                    if not local.exists():
                        raise FileNotFoundError(
                            f"Local {name} Parquet input not found: {local}"
                        )
                    sources[name] = local
                else:
                    sources[name] = download_dataset(
                        url,
                        args.cache_dir,
                        offline=args.offline,
                        refresh=args.refresh,
                        timeout=args.timeout,
                        retries=args.retries,
                    )
            results["clinvar"] = generate_clinvar(
                sources["evidence"],
                sources["coding"],
                mapping,
                staged / "mave_clinvar_variants.tsv",
                batch_size=args.batch_size,
            )
            results["clinvar"]["inputs"] = {
                name: (
                    {"path": str(path.resolve())}
                    if isinstance(path, Path)
                    else {
                        "url": (
                            args.evidence_url if name == "evidence" else args.coding_url
                        ),
                        "files": [str(file.resolve()) for file in path],
                    }
                )
                for name, path in sources.items()
            }
        manifest = reference_manifest(
            args.id_map, args.metadata if args.only != "clinvar" else None, results
        )
        write_json(staged / "references.json", manifest)
        publish_references(staged, args.out)
    LOGGER.info("Saved references to %s", args.out)
    for name, result in results.items():
        LOGGER.info("%s: %d rows", name, result["rows"])
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, URLError, ValueError, ImportError, HTTPException) as exc:
        LOGGER.error("Reference generation failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
