#!/usr/bin/env python3
"""Build MAVE variant, drug and pocket tables from a downloaded source bundle."""

from __future__ import annotations

import argparse
import logging
import sqlite3
import sys
import tempfile
from pathlib import Path
from typing import Any, Sequence

import pandas as pd
import pyarrow.dataset as ds

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]

if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis._reference_io import (
    positive_int,
    publish_tables,
    verify_file,
    write_json,
)
from data_sources.mavedb.analysis._reference_sources import (
    DEFAULT_SOURCE_MANIFEST,
    load_id_map,
    load_sources,
)
from data_sources.mavedb.analysis._variant_annotations import (
    generate_drugs,
    generate_pharmacogenomics,
    generate_pockets,
)
from data_sources.mavedb.analysis._variant_data import (
    coding_index,
    generate_clinvar,
    generate_gwas,
)

LOGGER = logging.getLogger("process_variant_data")
DEFAULT_OUT = REPOSITORY_ROOT / "data_sources/mavedb/output/variant_data"
MANIFEST = "variant_data_manifest.json"
STAGES = ("clinvar", "gwas", "pharmacogenomics", "drugs", "pockets")
PARQUET_SOURCES = {
    "coding_variants",
    "evidence",
    "credible_sets",
    "pharmacogenomics",
    "drugs",
}
REQUIRED_SOURCES = {
    "clinvar": ("coding_variants", "evidence", "id_map"),
    "gwas": ("coding_variants", "credible_sets"),
    "pharmacogenomics": ("coding_variants", "pharmacogenomics", "id_map"),
    "drugs": ("drugs", "id_map"),
    "pockets": ("pockets", "id_map"),
}
OUTPUT_NAMES = {
    "clinvar": "mave_clinvar_variants.tsv",
    "gwas": "mave_gwas_variants.tsv",
    "pharmacogenomics": "mave_pharmacogenomics.tsv",
    "drugs": "mave_drugs.tsv",
    "pockets": "mave_pockets.tsv",
}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCE_MANIFEST)
    parser.add_argument("--only", choices=("all", *STAGES), default="all")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--batch-size", type=positive_int, default=65_536)
    return parser


def _source_path(sources: dict[str, Any], name: str) -> Path | tuple[Path, ...]:
    """Resolve and validate one source from the completed download manifest."""
    record = sources.get(name)
    if not isinstance(record, dict):
        raise ValueError(f"Downloads manifest is missing required source {name!r}")
    if name in PARQUET_SOURCES:
        files = record.get("files")
        if not isinstance(files, list) or not files:
            raise ValueError(
                f"Downloads manifest source {name!r} requires a nonempty files list"
            )
        if any(not isinstance(descriptor, dict) for descriptor in files):
            raise ValueError(
                f"Downloads manifest source {name!r} requires file descriptors"
            )
        paths = tuple(verify_file(descriptor) for descriptor in files)
        dataset = ds.dataset(
            [str(file) for file in paths], format="parquet", exclude_invalid_files=False
        )
        if next(dataset.get_fragments(), None) is None:
            raise ValueError(f"Downloaded {name} source contains no Parquet files")
        return paths
    return verify_file(record)


def _resolve_sources(
    manifest_path: Path, stages: Sequence[str]
) -> tuple[dict[str, Path | tuple[Path, ...]], dict[str, Any]]:
    manifest = load_sources(manifest_path)
    records = manifest["sources"]
    names = tuple(
        dict.fromkeys(name for stage in stages for name in REQUIRED_SOURCES[stage])
    )
    paths = {name: _source_path(records, name) for name in names}
    return paths, {name: records[name] for name in names}


def _generate_stages(
    stages: Sequence[str],
    sources: dict[str, Path | tuple[Path, ...]],
    mapping: pd.DataFrame | None,
    connection: sqlite3.Connection | None,
    staged: Path,
    batch_size: int,
) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    for stage in stages:
        destination = staged / OUTPUT_NAMES[stage]
        options = {"batch_size": batch_size}
        if stage == "clinvar":
            if mapping is None or connection is None:
                raise ValueError(
                    "ClinVar generation requires an ID map and coding index"
                )
            result = generate_clinvar(
                sources["evidence"], connection, mapping, destination, **options
            )
        elif stage == "gwas":
            if connection is None:
                raise ValueError("GWAS generation requires a coding index")
            result = generate_gwas(
                sources["credible_sets"], connection, destination, **options
            )
        elif stage == "pharmacogenomics":
            if mapping is None or connection is None:
                raise ValueError(
                    "Pharmacogenomics generation requires an ID map and coding index"
                )
            result = generate_pharmacogenomics(
                sources["pharmacogenomics"], connection, mapping, destination, **options
            )
        elif stage == "drugs":
            if mapping is None:
                raise ValueError("Drug generation requires an ID map")
            result = generate_drugs(sources["drugs"], mapping, destination, **options)
        else:
            if mapping is None:
                raise ValueError("Pocket generation requires an ID map")
            result = generate_pockets(
                sources["pockets"], mapping, destination, **options
            )
        results[stage] = {**result, "file": OUTPUT_NAMES[stage]}
        LOGGER.info("%s: %d rows", stage, result["rows"])
    return results


def run(args: argparse.Namespace) -> dict[str, Any]:
    stages = list(STAGES) if args.only == "all" else [args.only]
    sources, provenance = _resolve_sources(args.sources, stages)
    mapping = (
        load_id_map(sources["id_map"])
        if any(stage != "gwas" for stage in stages)
        else None
    )
    coding_stages = {"clinvar", "gwas", "pharmacogenomics"}.intersection(stages)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".variant-data-", dir=args.out.parent) as d:
        staged = Path(d)
        if coding_stages:
            genes = None
            if not {"gwas", "pharmacogenomics"}.intersection(stages):
                if mapping is None:
                    raise ValueError("ClinVar generation requires an ID map")
                genes = sorted(set(mapping["To"].dropna().astype(str)))
            with coding_index(
                sources["coding_variants"],
                staged,
                batch_size=args.batch_size,
                genes=genes,
            ) as connection:
                results = _generate_stages(
                    stages,
                    sources,
                    mapping,
                    connection,
                    staged,
                    args.batch_size,
                )
        else:
            results = _generate_stages(
                stages, sources, mapping, None, staged, args.batch_size
            )
        manifest = {
            "source_manifest": str(args.sources.resolve()),
            "stages": stages,
            "sources": provenance,
            "results": results,
        }
        write_json(staged / MANIFEST, manifest)
        publish_tables(
            staged,
            args.out,
            manifest=MANIFEST,
            managed={Path(OUTPUT_NAMES[stage]) for stage in stages},
        )
    LOGGER.info("Saved variant references to %s", args.out)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError, ImportError) as exc:
        LOGGER.error("Variant reference processing failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
