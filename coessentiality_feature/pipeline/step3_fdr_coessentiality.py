#!/usr/bin/env python3
"""
GLS FDR filtering and co-essentiality network pipeline.

Reads the three outputs of the GLS co-essentiality pipeline (p-value matrix,
sign matrix, and gene list), applies Benjamini-Hochberg FDR correction, filters
by a user-defined FDR threshold, and writes the significant gene-pair network
to a CSV file.

Usage:
    python step3_fdr_coessentiality.py \\
        --gls-p  <prefix_GLS_p.npy> \\
        --gls-sign <prefix_GLS_sign.npy> \\
        --genes  <prefix_genes.txt> \\
        --fdr    0.05 \\
        --output <network_output.csv>

Outputs:
    <output>  — CSV with columns: source, target, pvalue, pvalue_adj, corr_genes
"""

import argparse
import importlib.util
import os
import subprocess
import sys


def _ensure_deps():
    packages = {"numpy": "numpy", "pandas": "pandas", "statsmodels": "statsmodels"}
    missing = [pkg for mod, pkg in packages.items() if importlib.util.find_spec(mod) is None]
    if missing:
        print(f"Installing missing packages: {', '.join(missing)} ...")
        subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing)
        print("Done.")

_ensure_deps()

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def load_inputs(gls_p_path, gls_sign_path, genes_path):
    print(f"[1/4] Loading inputs ...")
    genes = pd.read_csv(genes_path, header=None).squeeze()
    GLS_p = pd.DataFrame(np.load(gls_p_path), columns=genes, index=genes)
    GLS_sign = pd.DataFrame(np.load(gls_sign_path), columns=genes, index=genes)
    print(f"      Genes: {len(genes)}")
    print(f"      GLS_p shape:    {GLS_p.shape}")
    print(f"      GLS_sign shape: {GLS_sign.shape}")
    return GLS_p, GLS_sign, genes


def apply_fdr(GLS_p, GLS_sign, fdr_threshold):
    print(f"[2/4] Stacking matrices and applying BH FDR correction ...")

    # Keep only upper triangle to avoid duplicate pairs and self-interactions
    stacked_p = GLS_p.stack()
    stacked_p = stacked_p[
        stacked_p.index.get_level_values(0) < stacked_p.index.get_level_values(1)
    ]
    print(f"      Unique gene pairs: {len(stacked_p):,}")
    print(f"      Median raw p-value: {stacked_p.median():.4e}")

    stacked_sign = GLS_sign.stack()
    stacked_sign = stacked_sign[
        stacked_sign.index.get_level_values(0) < stacked_sign.index.get_level_values(1)
    ]

    fdr = pd.Series(
        multipletests(stacked_p, method="fdr_bh")[1], index=stacked_p.index
    )

    combined = pd.concat([stacked_p, fdr, stacked_sign], axis=1)
    combined.columns = ["pvalue", "pvalue_adj", "corr_genes"]

    print(f"[3/4] Filtering at FDR <= {fdr_threshold} ...")
    significant = combined[combined["pvalue_adj"] <= fdr_threshold]
    print(f"      Significant pairs: {len(significant):,} of {len(combined):,}")

    return significant


def build_network_df(significant):
    return pd.DataFrame(
        {
            "source": significant.index.get_level_values(0),
            "target": significant.index.get_level_values(1),
            "pvalue": significant["pvalue"].values,
            "pvalue_adj": significant["pvalue_adj"].values,
            "corr_genes": significant["corr_genes"].values,
        }
    )


def run_pipeline(gls_p_path, gls_sign_path, genes_path, fdr_threshold, output_path):
    GLS_p, GLS_sign, genes = load_inputs(gls_p_path, gls_sign_path, genes_path)
    significant = apply_fdr(GLS_p, GLS_sign, fdr_threshold)

    print(f"[4/4] Writing network to {output_path} ...")
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    network_df = build_network_df(significant)
    network_df.to_csv(output_path, index=False)
    print(f"      {len(network_df):,} edges written.")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description="GLS FDR filtering and co-essentiality network pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--gls-p",
        required=True,
        metavar="FILE",
        help=".npy file with GLS p-value matrix (n_genes x n_genes)",
    )
    parser.add_argument(
        "--gls-sign",
        required=True,
        metavar="FILE",
        help=".npy file with GLS sign matrix (n_genes x n_genes)",
    )
    parser.add_argument(
        "--genes",
        required=True,
        metavar="FILE",
        help=".txt file with gene names, one per line (no header)",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        metavar="THRESHOLD",
        help="FDR threshold for Benjamini-Hochberg correction (default: 0.05)",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Output CSV file path for the significant co-essentiality network",
    )

    args = parser.parse_args()

    if not (0 < args.fdr <= 1):
        parser.error(f"--fdr must be between 0 and 1, got {args.fdr}")

    for path, flag in [
        (args.gls_p, "--gls-p"),
        (args.gls_sign, "--gls-sign"),
        (args.genes, "--genes"),
    ]:
        if not os.path.isfile(path):
            parser.error(f"{flag}: file not found: {path}")

    run_pipeline(args.gls_p, args.gls_sign, args.genes, args.fdr, args.output)


if __name__ == "__main__":
    main()
