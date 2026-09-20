"""
Generate DEA training records with a configurable gene-selection strategy.

Standalone from prepare_training_data.py by design: the selection strategy is
the only variable under test here, so the fetch/grouping/filtering/output-
template logic is duplicated verbatim rather than parameterizing the main
pipeline. This keeps the baseline pipeline untouched and this script fully
self-contained for comparison runs.

Selection strategies:

  flat_top_k    Always take the top K genes by |log2fc|, regardless of how
                many genes actually passed the significance filter. This is
                the strategy used in prepare_training_data.py.

  adaptive_k    Take min(K, n_significant) genes. A perturbation with only
                3 significant DEGs yields a 3-gene record instead of being
                padded to 10; a perturbation with 300 significant DEGs is
                still capped at K, not truncated further than necessary.

Usage:

    python3 generate_dea_variants_new.py \
        --dataset_id nadig_2025_hepg2 \
        --strategy adaptive \
        --output data/dea_adaptive/nadig_2025_hepg2.jsonl

    python3 generate_dea_variants_new.py \
        --dataset_id nadig_2025_hepg2 \
        --strategy flat \
        --output data/dea_flat/nadig_2025_hepg2.jsonl

Run once per DEA dataset, concatenate the outputs, then run
pipeline.py --split on the combined corpus before fine-tuning.
"""
import json
import logging
import argparse
from pathlib import Path

import pandas as pd

from catalogue_api_new import query_perturb_seq, get_dataset_metadata

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


def flat_top_k(gene_list, k=10):
    """Always return the top k genes by |log2fc|."""
    return sorted(gene_list, key=lambda x: abs(x[1]), reverse=True)[:k]


def adaptive_k(gene_list, k_max=10):
    """Return min(k_max, len(gene_list)) genes by |log2fc| — no padding."""
    n = min(k_max, len(gene_list))
    return sorted(gene_list, key=lambda x: abs(x[1]), reverse=True)[:n]


SELECTION_STRATEGIES = {
    "flat": flat_top_k,
    "adaptive": adaptive_k,
}


def build_dea_records(dataset_id, strategy_fn, max_records=5000, include_condition=True):
    """
    Fetch DEA data for a dataset and build training records, selecting
    up/down genes per perturbation using strategy_fn.
    """
    log.info(f"Fetching DEA data for {dataset_id}...")
    raw = query_perturb_seq(dataset_id=dataset_id, max_records=max_records)
    if not raw:
        log.warning(f"No data returned for {dataset_id}")
        return [], pd.DataFrame()

    metadata = get_dataset_metadata(dataset_id)
    log.info(f"Dataset metadata: {metadata}")

    perturbation_effects = {}
    for r in raw:
        # FIX (Aug 18): API renamed "gene_name" to "gene_symbol" — see
        # catalogue_api_new.py for the same fix and verification.
        raw_gene_name = r.get("perturbation", {}).get("gene_symbol", "unknown")
        if "|" in str(raw_gene_name):
            parts = raw_gene_name.split("|")
            if "control_nontargeting" in parts[1].lower():
                perturbed_gene = parts[0]
            else:
                continue  # true combo perturbation, skip
        else:
            perturbed_gene = raw_gene_name
        effect = r.get("effect", {})

        record_cell_type = effect.get("cell_type") or "unknown"
        group_key = (perturbed_gene, record_cell_type)
        if group_key not in perturbation_effects:
            perturbation_effects[group_key] = {
                "up_genes": [],
                "down_genes": [],
                "cell_type": record_cell_type,
                "n_total": r.get("perturbation", {}).get("n_total", 0),
                "n_up": r.get("perturbation", {}).get("n_up", 0),
                "n_down": r.get("perturbation", {}).get("n_down", 0),
            }

        direction = effect.get("direction", "")
        log2fc = effect.get("log2fc", 0.0)
        # FIX (Aug 18): same rename applies to the affected-gene field.
        affected_gene = effect.get("gene_symbol", "unknown")
        padj = effect.get("padj", 1.0)

        if padj < 0.05 and not affected_gene.startswith("ENSG"):
            if direction == "increased":
                perturbation_effects[group_key]["up_genes"].append(
                    (affected_gene, round(log2fc, 3))
                )
            elif direction == "decreased":
                perturbation_effects[group_key]["down_genes"].append(
                    (affected_gene, round(log2fc, 3))
                )

    log.info(f"Grouped into {len(perturbation_effects)} unique perturbations")

    records = []
    rows = []
    cell_line = metadata.get("cell_line", "unknown")
    disease = metadata.get("disease", "unknown")
    condition = disease if disease != "unknown" else "standard growth"
    perturbation_type = metadata.get("library_perturbation_type", "knockout") or "knockout"

    for (perturbed_gene, record_cell_type), effects in perturbation_effects.items():
        cell_line_display = cell_line if cell_line != "unknown" else record_cell_type

        up_genes = strategy_fn(effects["up_genes"])
        down_genes = strategy_fn(effects["down_genes"])

        up_str = (
            ", ".join([f"{g} ({fc:+.2f})" for g, fc in up_genes])
            if up_genes else "none detected"
        )
        down_str = (
            ", ".join([f"{g} ({fc:+.2f})" for g, fc in down_genes])
            if down_genes else "none detected"
        )

        n_shown_up = len(up_genes)
        n_shown_down = len(down_genes)
        n_significant_up = len(effects["up_genes"])
        n_significant_down = len(effects["down_genes"])

        output_text = (
            f"CRISPR-mediated {perturbation_type} of {perturbed_gene} in {cell_line_display} causes "
            f"upregulation of: {up_str}; "
            f"and downregulation of: {down_str}. "
            f"Top differentially expressed genes shown: "
            f"{n_shown_up} upregulated, {n_shown_down} downregulated "
            f"(filtered by padj < 0.05, ranked by absolute log2fc)."
        )

        record = {
            "instruction": (
                f"Predict the transcriptional response to CRISPR-mediated {perturbation_type} "
                f"of gene {perturbed_gene} in {cell_line_display} under {condition}. "
                f"Describe the key upregulated and downregulated genes."
            ),
            "input": (
                f"Gene: {perturbed_gene}. "
                f"Cell line: {cell_line_display}."
                + (f" Condition: {condition}." if include_condition and condition != "standard growth" else "")
            ),
            "output": output_text,
            "metadata": {
                "gene": perturbed_gene,
                "dataset_id": dataset_id,
                "cell_line": cell_line_display,
                "disease": disease,
                "top_up_genes": [g for g, _ in up_genes],
                "top_down_genes": [g for g, _ in down_genes],
                "n_shown_up": n_shown_up,
                "n_shown_down": n_shown_down,
                "n_significant_up": n_significant_up,
                "n_significant_down": n_significant_down,
                "n_total_api": effects["n_total"],
                "n_up_api": effects["n_up"],
                "n_down_api": effects["n_down"],
                "modality": "scPerturb-seq",
                "source": "perturbation_catalogue_api",
            },
        }
        records.append(record)
        rows.append({
            "gene": perturbed_gene,
            "cell_line": cell_line_display,
            "n_significant_up": n_significant_up,
            "n_significant_down": n_significant_down,
            "n_shown_up": n_shown_up,
            "n_shown_down": n_shown_down,
        })

    log.info(f"Built {len(records)} training records from {dataset_id}")
    return records, pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Generate DEA training records with a selectable gene-selection strategy")
    parser.add_argument("--dataset_id", required=True)
    parser.add_argument("--strategy", choices=list(SELECTION_STRATEGIES.keys()), required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max_records", type=int, default=5000)
    parser.add_argument("--no_condition", action="store_true")
    args = parser.parse_args()

    strategy_fn = SELECTION_STRATEGIES[args.strategy]
    records, df = build_dea_records(
        dataset_id=args.dataset_id,
        strategy_fn=strategy_fn,
        max_records=args.max_records,
        include_condition=not args.no_condition,
    )

    if records:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            for record in records:
                f.write(json.dumps(record) + "\n")
        log.info(f"Saved {len(records)} records to {args.output}")

        if not df.empty:
            avg_shown_up = df["n_shown_up"].mean()
            avg_sig_up = df["n_significant_up"].mean()
            print(f"\n{args.strategy} strategy — {args.dataset_id}")
            print(f"Records: {len(records)}")
            print(f"Avg significant up-genes per record: {avg_sig_up:.1f}")
            print(f"Avg shown up-genes per record: {avg_shown_up:.1f}")
            if args.strategy == "adaptive":
                capped = (df["n_significant_up"] > 10).sum()
                unpadded = (df["n_significant_up"] < 10).sum()
                print(f"Records capped at 10 (more than 10 significant): {capped}")
                print(f"Records with fewer than 10 significant (no longer padded): {unpadded}")

    print(f"\nSaved {len(records)} DEA training records to {args.output}")


if __name__ == "__main__":
    main()