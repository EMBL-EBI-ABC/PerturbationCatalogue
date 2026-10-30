#!/usr/bin/env python3
"""Compare MAVE studies and assess ClinVar discrimination through homologs.

Python implementation of ``analysis/mave_homology.R`` from
``ot_perturbation_catalog_mave_experiments``. Run with ``--help`` for input modes.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import tempfile
from pathlib import Path
from typing import Sequence

import pandas as pd

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis._homology_clinvar import evaluate_clinvar
from data_sources.mavedb.analysis._homology_inputs import (
    input_fingerprint,
    load_inputs,
    read_tsv,
)
from data_sources.mavedb.analysis._mave_homology import (
    build_homology_key,
    calculate_correlations,
)

LOGGER = logging.getLogger("mave_homology")
MAVEDB_ROOT = REPOSITORY_ROOT / "data_sources" / "mavedb"
_MANAGED_TABLES = (
    "correlations",
    "combined_clinvar",
    "clinvar_homologs",
    "clinvar_counts",
    "clinvar_roc",
    "clinvar_summary",
)
_MANAGED_FIGURES = (
    "overlapping_variants",
    "pearson_distribution",
    "spearman_distribution",
    "pearson_bars",
    "correlation_vs_percent_identity",
    "correlation_vs_percent_identity_filtered",
    "correlation_vs_percent_identity_boxplot",
    "clinvar_roc_curves",
    "clinvar_roc_reduced",
)


def _publish(staged: Path, output: Path) -> None:
    """Replace managed files, remove obsolete ones, and roll back on I/O errors.

    Individual renames are atomic; concurrent readers must still wait for a
    successful process exit. Unrelated files in the destination are preserved.
    """
    managed = {Path(f"{name}.tsv") for name in _MANAGED_TABLES}
    managed.update(
        Path("figures") / f"{name}.{extension}"
        for name in _MANAGED_FIGURES
        for extension in ("png", "pdf")
    )
    managed.add(Path("manifest.json"))
    ordered = sorted(managed - {Path("manifest.json")}) + [Path("manifest.json")]
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".homology-previous-", dir=output.parent
    ) as previous:
        backup = Path(previous)
        saved: list[Path] = []
        installed: list[Path] = []
        try:
            for relative in ordered:
                destination = output / relative
                if destination.exists():
                    if not destination.is_file():
                        raise OSError(
                            f"Managed output is not a regular file: {destination}"
                        )
                    saved_path = backup / relative
                    saved_path.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(destination, saved_path)
                    saved.append(relative)
                source = staged / relative
                if source.is_file():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(source, destination)
                    installed.append(relative)
        except (OSError, KeyboardInterrupt):
            for relative in reversed(installed):
                (output / relative).unlink(missing_ok=True)
            for relative in reversed(saved):
                os.replace(backup / relative, output / relative)
            raise


def _positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return number


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument(
        "--gmm-dir", type=Path, help="Native score<N>/fitted_params<N> batch directory"
    )
    inputs.add_argument(
        "--data-root", type=Path, help="Legacy R data directory containing maves/"
    )
    parser.add_argument(
        "--reference-root",
        type=Path,
        help="Original R project root; supplies legacy data and reference defaults",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        help="Generated Python reference directory (default: output/homology_references)",
    )
    parser.add_argument(
        "--metadata",
        "--meta",
        type=Path,
        help="Canonical metadata.tsv or legacy meta.tsv",
    )
    parser.add_argument("--model", choices=("modes", "flexible"), default="modes")
    parser.add_argument(
        "--homology",
        type=Path,
        help="Ensembl homology.tsv including aligned protein sequences",
    )
    parser.add_argument(
        "--id-map",
        type=Path,
        help="UniProt-to-Ensembl TSV (From/To or uniprot_id/ensembl_id)",
    )
    parser.add_argument(
        "--clinvar", type=Path, help="mave_clinvar_variants.tsv from the R pipeline"
    )
    parser.add_argument(
        "--skip-clinvar",
        action="store_true",
        help="Explicitly run only between-study comparisons",
    )
    parser.add_argument(
        "--out", "--output", type=Path, default=MAVEDB_ROOT / "output" / "homology"
    )
    parser.add_argument(
        "--formats", choices=("png", "pdf"), nargs="+", default=["png", "pdf"]
    )
    parser.add_argument("--dpi", type=_positive_int, default=150)
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Write analysis tables without rendering figures",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Recompute correlations even when input content matches",
    )
    return parser


def _resolve_paths(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    reference = args.reference_root
    if reference is not None and args.reference_dir is not None:
        parser.error("choose --reference-root or --reference-dir")
    if args.gmm_dir is None and args.data_root is None:
        if reference is None:
            parser.error("provide --gmm-dir, --data-root or --reference-root")
        args.data_root = reference / "data"
    if args.metadata is None:
        args.metadata = (
            args.data_root / "maves" / "meta.tsv"
            if args.data_root is not None
            else MAVEDB_ROOT / "output" / "processed" / "metadata.tsv"
        )
    for name, relative in (
        ("homology", "intermediates/ensembl/homology.tsv"),
        ("id_map", "meta/mavedb_uniprot_to_ensembl.tsv"),
        ("clinvar", "intermediates/ot_variant_data/mave_clinvar_variants.tsv"),
    ):
        if getattr(args, name) is None and reference is not None:
            setattr(args, name, reference / relative)
    local_references = (
        args.reference_dir or MAVEDB_ROOT / "output" / "homology_references"
    )
    for name, filename in (
        ("homology", "homology.tsv"),
        ("id_map", "mavedb_uniprot_to_ensembl.tsv"),
        ("clinvar", "mave_clinvar_variants.tsv"),
    ):
        if getattr(args, name) is None and not (
            name == "clinvar" and args.skip_clinvar
        ):
            setattr(args, name, local_references / filename)
    for name in ("homology", "id_map"):
        if getattr(args, name) is None:
            parser.error(f"provide --{name.replace('_', '-')} or --reference-root")
    if args.clinvar is None and not args.skip_clinvar:
        parser.error("provide --clinvar or explicitly use --skip-clinvar")
    if args.skip_clinvar and args.clinvar is not None and reference is None:
        parser.error("--clinvar and --skip-clinvar cannot be combined")


def run(args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    """Calculate all results before replacing any previously completed outputs."""
    inputs = load_inputs(
        metadata=args.metadata,
        gmm_dir=args.gmm_dir,
        data_root=args.data_root,
        model=args.model,
    )
    homology = read_tsv(args.homology)
    id_map = read_tsv(args.id_map)
    homology_key = build_homology_key(homology, id_map)
    paths = [*inputs.paths, args.homology, args.id_map]
    # Include the implementation itself, so algorithm changes invalidate the cache.
    paths.extend(
        Path(__file__).parent / name
        for name in (
            "_homology_inputs.py",
            "_mave_homology.py",
            "gmm.py",
            "mave_homology.py",
            "mave_summary/data.py",
            "_mave_qc.py",
            "../scripts/process_mavedb.py",
        )
    )
    fingerprints = input_fingerprint(paths)
    manifest_path = args.out / "manifest.json"
    correlation_path = args.out / "correlations.tsv"
    cached = False
    if (
        not args.refresh_cache
        and manifest_path.is_file()
        and correlation_path.is_file()
    ):
        try:
            old_manifest = json.loads(manifest_path.read_text())
            cached = old_manifest.get("correlation_inputs") == fingerprints
        except (json.JSONDecodeError, AttributeError):
            LOGGER.warning("Ignoring invalid homology cache manifest")
    if cached:
        LOGGER.info("Loading content-matched correlation cache")
        correlations = read_tsv(
            correlation_path, {"pearson", "spearman", "relationship"}
        )
    else:
        LOGGER.info("Calculating correlations")
        correlations = calculate_correlations(
            inputs.meta, inputs.scores, inputs.singles, inputs.params, homology_key
        )
    tables = {"correlations": correlations}
    if not args.skip_clinvar:
        clinvar = read_tsv(args.clinvar)
        tables.update(
            evaluate_clinvar(clinvar, inputs.scores, inputs.params, homology_key)
        )
    clinvar_roc = tables.get("clinvar_roc", pd.DataFrame())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".mave-homology-", dir=args.out.parent
    ) as staging:
        staged = Path(staging)
        for name, frame in tables.items():
            frame.to_csv(staged / f"{name}.tsv", sep="\t", index=False, na_rep="NA")
        if not args.no_plots:
            import matplotlib.pyplot as plt

            from data_sources.mavedb.analysis._homology_plots import (
                build_plots,
                save_plots,
            )

            plots = build_plots(correlations, clinvar_roc)
            try:
                save_plots(
                    plots, staged / "figures", formats=args.formats, dpi=args.dpi
                )
            finally:
                for figure in plots.values():
                    plt.close(figure)
        manifest = {
            "source": "ot_perturbation_catalog_mave_experiments/analysis/mave_homology.R",
            "correlation_inputs": fingerprints,
            "clinvar_inputs": (
                input_fingerprint([args.clinvar]) if not args.skip_clinvar else {}
            ),
            "correlations_from_cache": cached,
            "tables": {name: len(frame) for name, frame in tables.items()},
            "directed_pairs": True,
            "files": sorted(
                ["manifest.json"]
                + [
                    str(path.relative_to(staged))
                    for path in staged.rglob("*")
                    if path.is_file()
                ]
            ),
        }
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        _publish(staged, args.out)
    LOGGER.info("Saved %d directed comparisons to %s", len(correlations), args.out)
    if "clinvar_summary" in tables:
        LOGGER.info(
            "Pooled ClinVar results:\n%s",
            tables["clinvar_summary"].to_string(index=False),
        )
    return tables


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = _parser()
    args = parser.parse_args(argv)
    _resolve_paths(args, parser)
    try:
        run(args)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        LOGGER.error("Homology analysis failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
