#!/usr/bin/env python3
"""
Validate the GLS co-essentiality pipeline against CORUM protein complexes.

For PR #365 Comment 4 ("we need *some* sort of benchmark which takes pipeline
results as input and produces a numerical metric"). Replicates the spirit of
Wainberg et al.'s Figure 2 / Extended Data Figure 4 ("GLS improves recall of
known functional interactions in co-essential gene pairs"): for each gene with
at least one known CORUM complex-mate, rank all other genes by GLS p-value and
check whether its top-N partners (N=1..10) are enriched for true CORUM pairs
relative to chance.

Gold standard: Enrichr's "CORUM" gene-set library (1,658 human complexes,
HGNC gene symbols — fetched once via gseapy, no ID-mapping needed since these
already match DepMap's gene symbol convention).

Metric ("enrichment", matching the paper's definition):
    enrichment(N) = (% of top-N pairs that are true CORUM pairs)
                  / (% of all possible pairs that are true CORUM pairs)
A value of 1.0 means GLS ranking is no better than random; the paper reports
several-fold enrichment for GLS at low N.

Usage:
    python validate_corum_enrichment.py \\
        --gls-p   ../required_data/depmap_26Q1_GLS_p.npy \\
        --genes   ../required_data/depmap_26Q1_genes.txt \\
        --corum   CORUM.gmt
"""

import argparse
import os

import gseapy as gp
import numpy as np

MAX_N = 10


def load_corum_partners(corum_path, organism, panel_genes):
    """Return {gene: set(complex-mate genes)}, restricted to the gene panel."""
    if os.path.isfile(corum_path):
        library = gp.get_library(name=corum_path)
    else:
        library = gp.get_library(name="CORUM", organism=organism, save=corum_path)

    panel_set = set(panel_genes)
    partners = {g: set() for g in panel_genes}
    for complex_genes in library.values():
        members = [g for g in complex_genes if g in panel_set]
        for i, gi in enumerate(members):
            for gj in members:
                if gi != gj:
                    partners[gi].add(gj)
    return partners


def compute_enrichment(gls_p_path, genes_path, corum_path, organism="Human"):
    print("[1/4] Loading GLS p-value matrix and gene list ...")
    genes = [g.strip() for g in open(genes_path) if g.strip()]
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)
    GLS_p = np.load(gls_p_path, mmap_mode="r")
    print(f"      {n_genes} genes, matrix shape {GLS_p.shape}")

    print("[2/4] Loading CORUM gold-standard complexes ...")
    corum_partners = load_corum_partners(corum_path, organism, genes)
    n_with_partners = sum(1 for v in corum_partners.values() if v)
    total_corum_pairs = sum(len(v) for v in corum_partners.values()) // 2
    total_possible_pairs = n_genes * (n_genes - 1) // 2
    background_rate = total_corum_pairs / total_possible_pairs
    print(f"      {n_with_partners} genes have >=1 CORUM partner in the panel")
    print(f"      {total_corum_pairs:,} true CORUM pairs / {total_possible_pairs:,} "
          f"possible pairs -> background rate = {background_rate:.3e}")

    print(f"[3/4] Ranking each gene's top-{MAX_N} GLS partners ...")
    hits_at_n = np.zeros(MAX_N, dtype=np.int64)
    n_evaluated_genes = 0
    for gene, true_partners in corum_partners.items():
        if not true_partners:
            continue
        idx = gene_to_idx[gene]
        row = np.array(GLS_p[idx, :], dtype=np.float64, copy=True)
        row[idx] = np.inf  # exclude self
        top_idx = np.argpartition(row, MAX_N)[:MAX_N]
        top_idx = top_idx[np.argsort(row[top_idx])]
        top_genes = [genes[i] for i in top_idx]
        for n in range(1, MAX_N + 1):
            if top_genes[n - 1] in true_partners:
                hits_at_n[n - 1:] += 1
        n_evaluated_genes += 1

    print(f"      Evaluated {n_evaluated_genes} genes")

    print("[4/4] Computing enrichment ...")
    print()
    print(f"{'N':>3}  {'observed rate':>15}  {'enrichment':>11}")
    cumulative_slots = 0
    results = []
    for n in range(1, MAX_N + 1):
        cumulative_slots = n_evaluated_genes * n
        observed_rate = hits_at_n[n - 1] / cumulative_slots
        enrichment = observed_rate / background_rate
        results.append((n, observed_rate, enrichment))
        print(f"{n:>3}  {observed_rate:>15.3e}  {enrichment:>11.2f}")

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gls-p", required=True, help="Path to GLS_p.npy")
    parser.add_argument("--genes", required=True, help="Path to genes.txt")
    parser.add_argument("--corum", default="CORUM.gmt",
                         help="Path to cache the CORUM .gmt (downloaded once if absent)")
    parser.add_argument("--organism", default="Human")
    args = parser.parse_args()

    compute_enrichment(args.gls_p, args.genes, args.corum, args.organism)


if __name__ == "__main__":
    main()
