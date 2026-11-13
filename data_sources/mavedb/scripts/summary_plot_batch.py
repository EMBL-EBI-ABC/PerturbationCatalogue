#!/usr/bin/env python3
"""Generate summary PDFs for MAVE GMM score sets."""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
MAVEDB_ROOT = REPOSITORY_ROOT / "data_sources" / "mavedb"
DEFAULT_META = MAVEDB_ROOT / "output" / "processed" / "metadata.tsv"
DEFAULT_OUT = MAVEDB_ROOT / "output" / "figures"
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis.mave_summary.data import (
    TableDataError,
    load_native_tables,
    load_tables,
)

LOGGER = logging.getLogger("summary_plot_batch")


def _nonnegative_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return parsed


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot a batch of MAVE studies using a fitted GMM model."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--gmm-dir",
        type=Path,
        help=(
            "Collected GMM model directory, its merged/ subdirectory, or a "
            "legacy directory of indexed native TSVs."
        ),
    )
    source.add_argument(
        "--data-root",
        type=Path,
        help="Explicit legacy export root containing maves/meta.tsv and gmm_<model> tables.",
    )
    parser.add_argument("--meta", type=Path, default=DEFAULT_META)
    parser.add_argument(
        "--study_type",
        default="maves",
        choices=("maves", "mpras", "splicing"),
        help="Study type to process (default: maves).",
    )
    parser.add_argument(
        "--model",
        default="modes",
        choices=("modes", "flexible"),
        help="GMM model to process (default: modes).",
    )
    parser.add_argument(
        "--index",
        type=_nonnegative_int,
        default=0,
        help="Batch index: 0 processes all score sets; positive indices are 1-based.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help="Output root directory (default: this repository's MaveDB output/figures).",
    )
    parser.add_argument(
        "--batch-size",
        type=_positive_int,
        default=100,
        help="Number of score sets per positive batch index (default: 100).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.study_type != "maves":
        parser.error(
            f"study type {args.study_type!r} is not implemented; only 'maves' is supported"
        )

    try:
        if args.gmm_dir is not None:
            tables, urns = load_native_tables(
                gmm_dir=args.gmm_dir,
                meta_path=args.meta,
                model=args.model,
                index=args.index,
                batch_size=args.batch_size,
            )
        else:
            tables, urns = load_tables(
                data_root=args.data_root,
                model=args.model,
                index=args.index,
                batch_size=args.batch_size,
            )
    except (TableDataError, FileNotFoundError) as exc:
        LOGGER.error("Input data error: %s", exc)
        return 2

    output_root = args.out / args.study_type / f"studies_{args.model}"
    output_root.mkdir(parents=True, exist_ok=True)
    if not urns:
        LOGGER.info(
            "No score sets selected for batch index %d; nothing to plot.", args.index
        )
        return 0

    from data_sources.mavedb.analysis.mave_summary.plots import make_summary_plot
    from matplotlib import pyplot as plt

    total = len(urns)
    for position, urn in enumerate(urns, start=1):
        LOGGER.info("Generating figure %d/%d: %s", position, total, urn)
        figure = None
        try:
            figure = make_summary_plot(urn, tables)
            filename = re.sub(r"[:-]", "_", urn) + ".pdf"
            figure.savefig(output_root / filename, format="pdf")
        finally:
            if figure is not None:
                plt.close(figure)

    LOGGER.info("Saved %d summary figure(s) to %s", total, output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
