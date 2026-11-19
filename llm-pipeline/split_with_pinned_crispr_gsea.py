"""
Gene-level split with CRISPR/GSEA assignments pinned to a reference split.

Problem this solves: pipeline.py --split computes a fresh random gene-level
split across the WHOLE combined corpus every time it runs. Since CRISPR, DEA,
and GSEA share one gene universe in the split, changing DEA's data (or
instructions) reshuffles which CRISPR/GSEA genes land in train vs. test --
even when CRISPR/GSEA's own training data never changed. This makes
CRISPR/GSEA metric comparisons across experiments unreliable.

Fix: read a reference split manifest (e.g. exp_002's), reuse its exact
gene -> {train, val, test} assignment for any gene that's a CRISPR or GSEA
gene in the new corpus. Only genes that are DEA-only, or genes not present in
the reference manifest at all, get freshly assigned via the same 80/10/10
random split logic as before.

Usage:

    python3 split_with_pinned_crispr_gsea.py \
        --corpus data/full_corpus_exp00X.jsonl \
        --reference_manifest runs/exp_002/split_manifest.json \
        --output_dir splits_exp00X_pinned/
"""
import json
import argparse
import random
from pathlib import Path
from collections import defaultdict


def load_records(path):
    return [json.loads(l) for l in open(path) if l.strip()]


def build_gene_modality_map(records):
    """gene -> set of modalities it appears under in this corpus."""
    gene_modalities = defaultdict(set)
    for r in records:
        gene = r["metadata"]["gene"]
        modality = r["metadata"].get("modality", "")
        gene_modalities[gene].add(modality)
    return gene_modalities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--reference_manifest", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--train_frac", type=float, default=0.8)
    parser.add_argument("--val_frac", type=float, default=0.1)
    args = parser.parse_args()

    records = load_records(args.corpus)
    print(f"Loaded {len(records)} records from {args.corpus}")

    with open(args.reference_manifest) as f:
        ref_manifest = json.load(f)
    # split_manifest.json format confirmed by direct inspection (Aug 19):
    # a FLAT dict of {gene_name: split_name}, e.g. {"ARPIN": "train", ...}
    # NOT nested {"train": [...], "val": [...], "test": [...]} as first assumed.
    ref_gene_split = dict(ref_manifest)
    print(f"Reference manifest: {len(ref_gene_split)} genes with a known split assignment")

    gene_modalities = build_gene_modality_map(records)
    crispr_gsea_modalities = {"CRISPR_screen", "scPerturb-seq_GSEA"}

    pinned = {}
    to_assign = []
    for gene, modalities in gene_modalities.items():
        is_crispr_or_gsea = bool(modalities & crispr_gsea_modalities)
        if is_crispr_or_gsea and gene in ref_gene_split:
            pinned[gene] = ref_gene_split[gene]
        else:
            to_assign.append(gene)

    print(f"Pinned (CRISPR/GSEA genes matched to reference): {len(pinned)}")
    print(f"To freshly assign (new/DEA-only genes): {len(to_assign)}")

    rng = random.Random(args.seed)
    rng.shuffle(to_assign)
    n_train = int(len(to_assign) * args.train_frac)
    n_val = int(len(to_assign) * args.val_frac)

    fresh_split = {}
    for gene in to_assign[:n_train]:
        fresh_split[gene] = "train"
    for gene in to_assign[n_train:n_train + n_val]:
        fresh_split[gene] = "val"
    for gene in to_assign[n_train + n_val:]:
        fresh_split[gene] = "test"

    gene_split = {**pinned, **fresh_split}

    splits = {"train": [], "val": [], "test": []}
    for r in records:
        gene = r["metadata"]["gene"]
        split_name = gene_split.get(gene)
        if split_name is None:
            continue
        splits[split_name].append(r)

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    for split_name, split_records in splits.items():
        out_path = Path(args.output_dir) / f"{split_name}.jsonl"
        with open(out_path, "w") as f:
            for r in split_records:
                f.write(json.dumps(r) + "\n")
        print(f"{split_name}: {len(split_records)} records -> {out_path}")

    # Write in the SAME flat {gene: split_name} format as the reference
    # manifest, confirmed by direct inspection of exp_002's real file --
    # so this output can itself be used as a --reference_manifest later.
    manifest = dict(gene_split)
    stats_path = Path(args.output_dir) / "split_stats.json"
    with open(stats_path, "w") as f:
        json.dump({
            "n_pinned_from_reference": len(pinned),
            "n_freshly_assigned": len(to_assign),
            "n_train": sum(1 for v in manifest.values() if v == "train"),
            "n_val": sum(1 for v in manifest.values() if v == "val"),
            "n_test": sum(1 for v in manifest.values() if v == "test"),
        }, f, indent=2)
    print(f"Split stats saved to {stats_path}")
    manifest_path = Path(args.output_dir) / "split_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest saved to {manifest_path}")

    # Sanity check: verify CRISPR/GSEA genes actually match the reference
    mismatches = 0
    for gene, modalities in gene_modalities.items():
        if (modalities & crispr_gsea_modalities) and gene in ref_gene_split:
            if gene_split.get(gene) != ref_gene_split[gene]:
                mismatches += 1
    print(f"\nSanity check -- CRISPR/GSEA genes not matching reference split: {mismatches} (should be 0)")


if __name__ == "__main__":
    main()