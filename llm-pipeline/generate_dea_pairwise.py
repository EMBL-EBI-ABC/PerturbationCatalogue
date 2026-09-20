"""
Generate one training record per (perturbed gene, affected gene) pair,
using every significant gene per perturbation (padj < 0.05), no top-k cap.

Instruction format: "After CRISPR-mediated inhibition of gene A in cell line
DEF, what happens to gene B?" -> "gene B is downregulated (log2fc = -1.58)"

Usage:

    python3 generate_dea_pairwise.py \
        --dataset_id nadig_2025_hepg2 \
        --output data/dea_pairwise/nadig_2025_hepg2.jsonl
"""
import json
import logging
import argparse
from pathlib import Path

import pandas as pd

import requests
from catalogue_api_new import query_perturb_seq, get_dataset_metadata

BASE_URL = "https://perturbation-catalogue-be-328296435987.europe-west2.run.app"


def resolve_ensg_to_symbol(ensg_id, cache):
    """Resolve an ENSG ID to a gene symbol via the /v1/target endpoint, caching results."""
    if ensg_id in cache:
        return cache[ensg_id]
    try:
        resp = requests.get(f"{BASE_URL}/v1/target/{ensg_id}", timeout=15)
        resp.raise_for_status()
        symbol = resp.json().get("approved_symbol")
        cache[ensg_id] = symbol
        return symbol
    except Exception:
        cache[ensg_id] = None
        return None

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


def build_pairwise_dea_records(dataset_id, max_records=5000, include_condition=True):
    log.info(f"Fetching DEA data for {dataset_id}...")
    raw = query_perturb_seq(dataset_id=dataset_id, max_records=max_records)
    if not raw:
        log.warning(f"No data returned for {dataset_id}")
        return [], pd.DataFrame()

    metadata = get_dataset_metadata(dataset_id)
    log.info(f"Dataset metadata: {metadata}")

    cell_line = metadata.get("cell_line", "unknown")
    disease = metadata.get("disease", "unknown")
    condition = disease if disease != "unknown" else "standard growth"
    perturbation_type = metadata.get("library_perturbation_type", "knockout") or "knockout"

    records = []
    rows = []
    n_skipped_combo = 0
    n_below_significance = 0
    n_skipped_invalid_perturbed_gene = 0
    n_resolved_via_ensg = 0
    ensg_cache = {}

    for r in raw:
        raw_gene_name = r.get("perturbation", {}).get("gene_symbol")
        if raw_gene_name is None:
            ensg_field = r.get("perturbation", {}).get("perturbed_target_ensg")
            if ensg_field:
                parts = [p for p in str(ensg_field).split("|") if p.startswith("ENSG")]
                if parts:
                    symbol = resolve_ensg_to_symbol(parts[0], ensg_cache)
                    if symbol:
                        raw_gene_name = symbol
                        n_resolved_via_ensg += 1
            if raw_gene_name is None:
                n_skipped_invalid_perturbed_gene += 1
                continue
        elif str(raw_gene_name).startswith("ENSG"):
            symbol = resolve_ensg_to_symbol(raw_gene_name, ensg_cache)
            if symbol:
                raw_gene_name = symbol
                n_resolved_via_ensg += 1
            else:
                n_skipped_invalid_perturbed_gene += 1
                continue
        if "|" in str(raw_gene_name):
            parts = raw_gene_name.split("|")
            if "control_nontargeting" in parts[1].lower():
                perturbed_gene = parts[0]
            else:
                n_skipped_combo += 1
                continue
        else:
            perturbed_gene = raw_gene_name

        effect = r.get("effect", {})
        record_cell_type = effect.get("cell_type") or "unknown"
        cell_line_display = cell_line if cell_line != "unknown" else record_cell_type

        direction = effect.get("direction", "")
        log2fc = effect.get("log2fc", 0.0)
        affected_gene = effect.get("gene_symbol", "unknown")
        padj = effect.get("padj", 1.0)

        if padj >= 0.05 or affected_gene.startswith("ENSG") or direction not in ("increased", "decreased"):
            n_below_significance += 1
            continue

        direction_word = "upregulated" if direction == "increased" else "downregulated"

        output_text = (
            f"After CRISPR-mediated {perturbation_type} of {perturbed_gene} in "
            f"{cell_line_display}, gene {affected_gene} is {direction_word} "
            f"(log2fc = {log2fc:+.2f})."
        )

        record = {
            "instruction": (
                f"After CRISPR-mediated {perturbation_type} of gene {perturbed_gene} "
                f"in {cell_line_display} under {condition}, what happens to gene "
                f"{affected_gene}?"
            ),
            "input": (
                f"Perturbed gene: {perturbed_gene}. Affected gene: {affected_gene}. "
                f"Cell line: {cell_line_display}."
                + (f" Condition: {condition}." if include_condition and condition != "standard growth" else "")
            ),
            "output": output_text,
            "metadata": {
                "perturbed_gene": perturbed_gene,
                "affected_gene": affected_gene,
                "gene": perturbed_gene,
                "dataset_id": dataset_id,
                "cell_line": cell_line_display,
                "disease": disease,
                "direction": direction_word,
                "log2fc": round(log2fc, 3),
                "padj": padj,
                "modality": "scPerturb-seq_pairwise",
                "source": "perturbation_catalogue_api",
            },
        }
        records.append(record)
        rows.append({
            "perturbed_gene": perturbed_gene,
            "affected_gene": affected_gene,
            "direction": direction_word,
            "log2fc": log2fc,
        })

    log.info(f"Built {len(records)} pairwise records from {dataset_id} "
              f"(skipped {n_skipped_combo} combos, {n_below_significance} non-significant, "
              f"{n_skipped_invalid_perturbed_gene} unresolvable perturbed_gene, "
              f"{n_resolved_via_ensg} resolved via ENSG lookup)")
    return records, pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Generate pairwise DEA training records")
    parser.add_argument("--dataset_id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max_records", type=int, default=5000)
    parser.add_argument("--no_condition", action="store_true")
    args = parser.parse_args()

    records, df = build_pairwise_dea_records(
        dataset_id=args.dataset_id,
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
            n_perturbations = df["perturbed_gene"].nunique()
            print(f"\n{args.dataset_id}: {len(records)} pairwise records "
                  f"across {n_perturbations} unique perturbed genes "
                  f"({len(records)/n_perturbations:.1f} affected genes per perturbation, avg)")

    print(f"\nSaved {len(records)} pairwise DEA training records to {args.output}")


if __name__ == "__main__":
    main()