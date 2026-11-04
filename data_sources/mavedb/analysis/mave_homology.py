#!/usr/bin/env python3
"""Compare MAVE studies and assess ClinVar discrimination through homologs.

Python implementation of ``analysis/mave_homology.R`` from
``ot_perturbation_catalog_mave_experiments``. Run with ``--help`` for input modes.
"""

from __future__ import annotations

import argparse
import json
import logging
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
)
from data_sources.mavedb.analysis._mave_homology import (
    build_homology_key,
    calculate_correlations,
)
from data_sources.mavedb.analysis._reference_io import (
    publish_tables,
    read_tsv,
    verify_file,
)
from data_sources.mavedb.analysis._reference_sources import (
    DEFAULT_SOURCE_MANIFEST,
    load_sources,
)

LOGGER = logging.getLogger("mave_homology")
MAVEDB_ROOT = REPOSITORY_ROOT / "data_sources" / "mavedb"
DEFAULT_METADATA = MAVEDB_ROOT / "output" / "processed" / "metadata.tsv"
DEFAULT_CLINVAR = MAVEDB_ROOT / "output" / "variant_data" / "mave_clinvar_variants.tsv"
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


def _positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return number


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gmm-dir",
        type=Path,
        required=True,
        help="Native score<N>/fitted_params<N> batch directory",
    )
    parser.add_argument(
        "--sources",
        type=Path,
        default=DEFAULT_SOURCE_MANIFEST,
        help="references/manifest.json from download_reference_data.py",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Canonical processed metadata (default: manifest metadata, then output/processed/metadata.tsv)",
    )
    parser.add_argument("--clinvar", type=Path, help="Generated ClinVar variant TSV")
    parser.add_argument(
        "--skip-clinvar",
        action="store_true",
        help="Run only between-study comparisons",
    )
    parser.add_argument("--out", type=Path, default=MAVEDB_ROOT / "output" / "homology")
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


def _resolve_paths(args: argparse.Namespace) -> None:
    """Resolve canonical references from the unified local download manifest."""
    sources = load_sources(args.sources)["sources"]
    for name in ("homology", "id_map"):
        if name not in sources:
            raise ValueError(f"Downloads manifest is missing required source {name!r}")
        setattr(args, name, verify_file(sources[name]))
    if args.metadata is None:
        args.metadata = (
            verify_file(sources["metadata"])
            if "metadata" in sources
            else DEFAULT_METADATA
        )
    if args.clinvar is None and not args.skip_clinvar:
        args.clinvar = DEFAULT_CLINVAR
    if args.skip_clinvar and args.clinvar is not None:
        raise ValueError("--clinvar and --skip-clinvar cannot be combined")


def run(args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    """Calculate all results before replacing any previously completed outputs."""
    inputs = load_inputs(
        metadata=args.metadata,
        gmm_dir=args.gmm_dir,
    )
    fitted_urns = sorted(inputs.scores["score_urn"].dropna().unique().tolist())
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
    if correlations.empty:
        if len(fitted_urns) < 2:
            LOGGER.warning(
                "No correlations: at least two fitted score sets are required; "
                "loaded %d. Fit additional comparable studies before rerunning.",
                len(fitted_urns),
            )
        else:
            LOGGER.warning(
                "No correlations: no eligible pairs share at least three "
                "aligned variants with nonconstant scores in both studies."
            )
    tables = {"correlations": correlations}
    if not args.skip_clinvar:
        clinvar = read_tsv(args.clinvar)
        tables.update(
            evaluate_clinvar(clinvar, inputs.scores, inputs.params, homology_key)
        )
        if tables["clinvar_summary"].empty:
            fitted_proteins = inputs.scores["identifier_uniprot"].dropna()
            targets = homology_key.loc[
                homology_key["source_uniprot"].isin(fitted_proteins),
                "target_uniprot",
            ].dropna()
            if targets.empty:
                reason = "no mapped homolog proteins for the fitted score sets"
            elif (
                not tables["combined_clinvar"]["identifier_uniprot"].isin(targets).any()
            ):
                reason = "no ClinVar variants for the mapped homolog proteins"
            elif tables["clinvar_homologs"].empty:
                reason = "no ClinVar substitutions match the aligned fitted variants"
            else:
                reason = "matched variants do not support a binary ClinVar evaluation"
            LOGGER.warning("No pooled ClinVar results: %s.", reason)
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
            "source_manifest": str(args.sources.resolve()),
            "clinvar_inputs": (
                input_fingerprint([args.clinvar]) if not args.skip_clinvar else {}
            ),
            "correlations_from_cache": cached,
            "input_coverage": {
                "eligible_metadata_score_sets": len(inputs.meta),
                "fitted_score_sets": len(fitted_urns),
                "fitted_score_urns": fitted_urns,
            },
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
        managed = {Path(f"{name}.tsv") for name in _MANAGED_TABLES}
        managed.update(
            Path("figures") / f"{name}.{extension}"
            for name in _MANAGED_FIGURES
            for extension in ("png", "pdf")
        )
        publish_tables(
            staged,
            args.out,
            manifest="manifest.json",
            managed=managed,
        )
    LOGGER.info("Saved %d directed comparisons to %s", len(correlations), args.out)
    if "clinvar_summary" in tables and not tables["clinvar_summary"].empty:
        LOGGER.info(
            "Pooled ClinVar results:\n%s",
            tables["clinvar_summary"].to_string(index=False),
        )
    return tables


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        _resolve_paths(args)
        run(args)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        LOGGER.error("Homology analysis failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
