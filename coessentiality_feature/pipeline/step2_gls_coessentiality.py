#!/usr/bin/env python3
"""
GLS co-essentiality pipeline.

Reads an essentiality matrix CSV, cleans gene symbol column names,
drops genes with any NA values, then computes GLS p-values and sign matrix
for all gene pairs.

Usage:
    python step2_gls_coessentiality.py <input_csv> [--output-dir DIR] [--prefix NAME]

Outputs (written to output_dir):
    <prefix>_genes.txt       — gene list after NA filtering
    <prefix>_GLS_p.npy      — n x n matrix of two-sided p-values
    <prefix>_GLS_sign.npy   — n x n matrix of correlation signs (+1 / -1)
    <prefix>_metadata.json  — n_cell_lines, n_genes_profiled
"""

import argparse
import json
import os
import re

import numpy as np
import pandas as pd
from scipy.special import stdtr


def clean_column_names(df):
    """Strip Entrez gene IDs from DepMap-style column names.

    Converts 'GENE_SYMBOL (ENTREZ_ID)' -> 'GENE_SYMBOL' for every column.
    Columns that don't match the pattern are left unchanged.

    Raises ValueError if stripping the Entrez ID produces duplicate gene
    symbols (e.g. aliased/repeated symbols across different Entrez IDs) —
    silently merging them would corrupt the covariance/GLS computation.
    """
    cleaned = df.rename(columns=lambda x: re.sub(r"\s*\(\d+\)$", "", x))
    counts = cleaned.columns.value_counts()
    duplicates = counts[counts > 1]
    if not duplicates.empty:
        raise ValueError(
            f"Stripping Entrez IDs produced {len(duplicates)} duplicate gene symbol(s): "
            f"{duplicates.index.tolist()}. Resolve the underlying alias collision before proceeding."
        )
    return cleaned


def linear_regression(warped_screens, warped_intercept):
    n_genes = len(warped_screens)
    GLS_coef = np.empty((n_genes, n_genes))
    GLS_se = np.empty((n_genes, n_genes))
    ys = warped_screens.T

    for gene_index in range(n_genes):
        X = np.stack((warped_intercept, warped_screens[gene_index]), axis=1)
        coef, residues = np.linalg.lstsq(X, ys, rcond=None)[:2]
        df = warped_screens.shape[1] - 2
        GLS_coef[gene_index] = coef[1]
        GLS_se[gene_index] = np.sqrt(np.linalg.pinv(X.T @ X)[1, 1] * residues / df)

    return GLS_coef, GLS_se


def run_pipeline(input_file, output_dir, prefix):
    os.makedirs(output_dir, exist_ok=True)

    # ── Load ────────────────────────────────────────────────────────────────
    print(f"[1/6] Loading data from {input_file} ...")
    data = pd.read_csv(input_file, index_col=0)
    n_cell_lines = data.shape[0]
    print(f"      {data.shape[0]} screens x {data.shape[1]} genes")

    # ── Clean column names ───────────────────────────────────────────────────
    print("[2/6] Cleaning column names ...")
    before = list(data.columns[:3])
    data = clean_column_names(data)
    after = list(data.columns[:3])
    if before != after:
        print(f"      Example: {before[0]!r} -> {after[0]!r}")
    else:
        print("      Columns already clean, no changes made.")

    # ── NA handling ─────────────────────────────────────────────────────────
    print("[3/6] Handling NA values ...")
    columns_with_na = data.columns[data.isna().any()]
    if len(columns_with_na) == 0:
        print("      No NA values found.")
    else:
        print(f"      {len(columns_with_na)} genes with NAs — dropping all.")
        data = data.drop(columns=columns_with_na)

    assert not data.isnull().any().any(), "NA values remain after cleaning."

    # ── Export genes ────────────────────────────────────────────────────────
    print("[4/6] Exporting gene list ...")
    genes_path = os.path.join(output_dir, f"{prefix}_genes.txt")
    data.T.index.to_series().to_csv(genes_path, index=False, header=False)
    print(f"      {len(data.columns)} genes -> {genes_path}")

    metadata_path = os.path.join(output_dir, f"{prefix}_metadata.json")
    with open(metadata_path, "w") as fh:
        json.dump({"n_cell_lines": n_cell_lines, "n_genes_profiled": len(data.columns)}, fh, indent=2)
    print(f"      metadata -> {metadata_path}")

    # ── GLS decomposition ───────────────────────────────────────────────────
    print("[5/6] Computing GLS (Cholesky of pseudoinverse of covariance) ...")
    cov_pinv = np.linalg.pinv(np.cov(data))
    # MKL/Accelerate BLAS can push the one theoretically-zero eigenvalue of the
    # pseudoinverse slightly negative, breaking Cholesky. Clipping via eigh fixes this.
    eigvals, eigvecs = np.linalg.eigh(cov_pinv)
    cholsigmainv = np.linalg.cholesky(eigvecs @ np.diag(np.clip(eigvals, 0, None)) @ eigvecs.T)
    warped_screens = data.T.values @ cholsigmainv      # (n_genes, n_screens)
    warped_intercept = cholsigmainv.sum(axis=0)        # (n_screens,)

    print(f"[6/6] Running GLS regression ({warped_screens.shape[0]} genes x {warped_screens.shape[1]} screens) ...")
    GLS_coef, GLS_se = linear_regression(warped_screens, warped_intercept)

    df = warped_screens.shape[1] - 2
    GLS_p = 2 * stdtr(df, -np.abs(GLS_coef / GLS_se))
    np.fill_diagonal(GLS_p, 1)

    # ── Save outputs ────────────────────────────────────────────────────────
    gls_p_path = os.path.join(output_dir, f"{prefix}_GLS_p.npy")
    gls_sign_path = os.path.join(output_dir, f"{prefix}_GLS_sign.npy")
    np.save(gls_p_path, GLS_p)
    np.save(gls_sign_path, np.sign(GLS_coef))
    print(f"      GLS_p    -> {gls_p_path}")
    print(f"      GLS_sign -> {gls_sign_path}")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description="GLS co-essentiality pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input",
        help="Essentiality matrix CSV (rows = screens, columns = genes)",
    )
    parser.add_argument(
        "--output-dir",
        default="./output",
        help="Directory for output files (default: ./output)",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Prefix for output filenames (default: input filename without extension)",
    )
    args = parser.parse_args()

    if args.prefix is None:
        args.prefix = os.path.splitext(os.path.basename(args.input))[0]

    run_pipeline(args.input, args.output_dir, args.prefix)


if __name__ == "__main__":
    main()
