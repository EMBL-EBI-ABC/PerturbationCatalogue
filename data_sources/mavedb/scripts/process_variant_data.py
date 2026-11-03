#!/usr/bin/env python3
"""Build MAVE variant, drug and pocket reference tables from local source data."""

from __future__ import annotations

import argparse
import hashlib
import logging
import math
import os
import sys
import tempfile
from datetime import datetime, timezone
from http.client import HTTPException
from pathlib import Path
from typing import Any, Sequence
from urllib.error import URLError

import pyarrow.dataset as ds

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis._homology_downloads import (
    download_dataset,
    write_json,
)
from data_sources.mavedb.analysis._homology_references import (
    DEFAULT_CACHE,
    DEFAULT_ID_MAP,
    generate_clinvar,
    load_id_map,
)
from data_sources.mavedb.analysis._reference_sources import (
    DEFAULT_SOURCE_MANIFEST,
    SOURCE_URLS,
    load_sources,
)
from data_sources.mavedb.analysis._variant_annotations import (
    generate_drugs,
    generate_pharmacogenomics,
    generate_pockets,
)
from data_sources.mavedb.analysis._variant_data import generate_gwas

LOGGER = logging.getLogger("process_variant_data")
DEFAULT_OUT = REPOSITORY_ROOT / "data_sources/mavedb/output/variant_data"
MANIFEST = "variant_data_manifest.json"
STAGES = ("clinvar", "gwas", "pharmacogenomics", "drugs", "pockets")
SOURCE_NAMES = {
    "coding_variants": "coding_variant_table",
    "evidence": "evidence_eva",
    "credible_sets": "credible_set",
    "pharmacogenomics": "pharmacogenomics",
    "drugs": "drug_molecule",
}
REQUIRED_SOURCES = {
    "clinvar": ("coding_variants", "evidence"),
    "gwas": ("coding_variants", "credible_sets"),
    "pharmacogenomics": ("coding_variants", "pharmacogenomics"),
    "drugs": ("drugs",),
    "pockets": (),
}
OUTPUT_NAMES = {
    "clinvar": "mave_clinvar_variants.tsv",
    "gwas": "mave_gwas_variants.tsv",
    "pharmacogenomics": "mave_pharmacogenomics.tsv",
    "drugs": "mave_drugs.tsv",
    "pockets": "mave_pockets.tsv",
}


def _positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return number


def _positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("Must be finite and positive")
    return number


def _nonnegative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("Must be nonnegative")
    return number


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=("all", *STAGES), default="all")
    parser.add_argument("--id-map", type=Path)
    parser.add_argument(
        "--sources",
        type=Path,
        help="downloads.json from download_reference_data.py; detected automatically at its default location",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--input-root",
        "--ot-root",
        type=Path,
        help="Directory containing the original Open Targets dataset subdirectories",
    )
    for name, directory in SOURCE_NAMES.items():
        flags = ["--" + name.replace("_", "-")]
        if name == "credible_sets":
            flags.append("--credible-set")
        parser.add_argument(*flags, type=Path, help=f"Local {directory} Parquet input")
        parser.add_argument(
            "--" + name.replace("_", "-") + "-url",
            default=SOURCE_URLS[name],
            help=f"Optional FTP/HTTP directory listing for {directory} Parquet files",
        )
    parser.add_argument("--pockets", type=Path, help="Local ProtVar pocket TSV input")
    parser.add_argument("--batch-size", type=_positive_int, default=65_536)
    parser.add_argument("--timeout", type=_positive_float, default=60)
    parser.add_argument("--retries", type=_nonnegative_int, default=3)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--offline",
        action="store_true",
        help="Use local inputs or complete caches only",
    )
    modes.add_argument(
        "--refresh", action="store_true", help="Refresh downloaded source datasets"
    )
    return parser


def _source_plan(
    args: argparse.Namespace, stages: Sequence[str]
) -> dict[str, Path | str]:
    """Validate every selected stage before starting any downloads or processing."""
    names = dict.fromkeys(name for stage in stages for name in REQUIRED_SOURCES[stage])
    plan: dict[str, Path | str] = {}
    for name in names:
        local = getattr(args, name)
        if local is None and args.input_root is not None:
            local = args.input_root / SOURCE_NAMES[name]
        if local is not None:
            if not local.exists():
                raise FileNotFoundError(
                    f"Local {SOURCE_NAMES[name]} input not found: {local}"
                )
            dataset = ds.dataset(
                str(local), format="parquet", exclude_invalid_files=True
            )
            if next(dataset.get_fragments(), None) is None:
                raise ValueError(
                    f"Local {SOURCE_NAMES[name]} input contains no Parquet files: {local}"
                )
            plan[name] = local
        else:
            url = getattr(args, name + "_url")
            if not url:
                flag = name.replace("_", "-")
                raise ValueError(f"Provide --{flag}, --input-root, or --{flag}-url")
            plan[name] = url
    if "pockets" in stages and (args.pockets is None or not args.pockets.is_file()):
        raise FileNotFoundError(
            "The pockets stage requires --pockets with an existing ProtVar TSV"
        )
    return plan


def _source_record(paths: Path | Sequence[Path], source: Path | str) -> dict[str, Any]:
    files = (
        sorted(paths.rglob("*.parquet"))
        if isinstance(paths, Path) and paths.is_dir()
        else ([paths] if isinstance(paths, Path) else list(paths))
    )
    return {
        "path" if isinstance(source, Path) else "url": (
            str(source.resolve()) if isinstance(source, Path) else source
        ),
        "files": [
            {"path": str(path.resolve()), "bytes": path.stat().st_size}
            for path in files
        ],
    }


def _publish(staged: Path, output: Path) -> None:
    """Publish all selected tables together, restoring previous files on error."""
    output.mkdir(parents=True, exist_ok=True)
    names = sorted(path.name for path in staged.iterdir() if path.name != MANIFEST)
    names.append(MANIFEST)
    with tempfile.TemporaryDirectory(
        prefix=".variant-data-backup-", dir=output.parent
    ) as directory:
        backup = Path(directory)
        saved, installed = [], []
        try:
            for name in names:
                destination = output / name
                if destination.exists():
                    if not destination.is_file():
                        raise OSError(f"Output is not a regular file: {destination}")
                    os.replace(destination, backup / name)
                    saved.append(name)
                os.replace(staged / name, destination)
                installed.append(name)
        except (OSError, KeyboardInterrupt):
            for name in reversed(installed):
                (output / name).unlink(missing_ok=True)
            for name in reversed(saved):
                os.replace(backup / name, output / name)
            raise


def _resolve_sources(args: argparse.Namespace) -> None:
    path = args.sources
    if path is None and args.input_root is None and DEFAULT_SOURCE_MANIFEST.is_file():
        path = DEFAULT_SOURCE_MANIFEST
    sources = load_sources(path)["sources"] if path is not None else {}
    if args.input_root is None:
        for name in SOURCE_NAMES:
            if (
                getattr(args, name) is None
                and name in sources
                and getattr(args, name + "_url") == SOURCE_URLS[name]
            ):
                setattr(args, name, Path(sources[name]["path"]))
    if args.pockets is None and "pockets" in sources:
        args.pockets = Path(sources["pockets"]["path"])
    if args.id_map is None:
        args.id_map = (
            Path(sources["id_map"]["path"]) if "id_map" in sources else DEFAULT_ID_MAP
        )


def run(args: argparse.Namespace) -> dict[str, Any]:
    _resolve_sources(args)
    stages = list(STAGES) if args.only == "all" else [args.only]
    plan = _source_plan(args, stages)
    mapping = (
        load_id_map(args.id_map, deduplicate=False)
        if any(stage != "gwas" for stage in stages)
        else None
    )
    sources = {
        name: (
            source
            if isinstance(source, Path)
            else download_dataset(
                source,
                args.cache_dir,
                offline=args.offline,
                refresh=args.refresh,
                timeout=args.timeout,
                retries=args.retries,
            )
        )
        for name, source in plan.items()
    }
    provenance = {
        name: _source_record(paths, plan[name]) for name, paths in sources.items()
    }
    if mapping is not None:
        provenance["id_map"] = {
            "path": str(args.id_map.resolve()),
            "sha256": hashlib.sha256(args.id_map.read_bytes()).hexdigest(),
        }
    if "pockets" in stages:
        provenance["pockets"] = {
            "path": str(args.pockets.resolve()),
            "sha256": hashlib.sha256(args.pockets.read_bytes()).hexdigest(),
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results = {}
    with tempfile.TemporaryDirectory(
        prefix=".variant-data-", dir=args.out.parent
    ) as directory:
        staged = Path(directory)
        for stage in stages:
            destination = staged / OUTPUT_NAMES[stage]
            options = {"batch_size": args.batch_size}
            if stage == "clinvar":
                result = generate_clinvar(
                    sources["evidence"],
                    sources["coding_variants"],
                    mapping,
                    destination,
                    **options,
                )
            elif stage == "gwas":
                result = generate_gwas(
                    sources["credible_sets"],
                    sources["coding_variants"],
                    destination,
                    **options,
                )
            elif stage == "pharmacogenomics":
                result = generate_pharmacogenomics(
                    sources["pharmacogenomics"],
                    sources["coding_variants"],
                    mapping,
                    destination,
                    **options,
                )
            elif stage == "drugs":
                result = generate_drugs(
                    sources["drugs"], mapping, destination, **options
                )
            else:
                result = generate_pockets(args.pockets, mapping, destination, **options)
            results[stage] = {**result, "file": OUTPUT_NAMES[stage]}
            LOGGER.info("%s: %d rows", stage, result["rows"])
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "stages": stages,
            "sources": provenance,
            "results": results,
        }
        write_json(staged / MANIFEST, manifest)
        _publish(staged, args.out)
    LOGGER.info("Saved variant references to %s", args.out)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    try:
        run(args)
    except (OSError, URLError, ValueError, ImportError, HTTPException) as exc:
        LOGGER.error("Variant reference generation failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
