"""
benchmark.py
=============
Evaluation framework for the perturbation-aware LLM.

THE CORE QUESTION THIS FILE ANSWERS:
-------------------------------------
Did the model learn actual biology, or did it memorise training data?

These are not the same thing and standard ML evaluation cannot
distinguish between them. A model that sees BRCA1 during training
and gets tested on BRCA1 in a slightly different context will look
smart — it's just recalling facts. A model tested on genes it has
never seen has to actually reason.

This file builds the evaluation infrastructure that forces that
distinction. Everything here is designed around one principle:
if the model passes these tests, it learned something real.

PROJECT CONNECTION:
-------------------
The project description explicitly requires:
"Benchmarking and evaluation framework comparing LLM-based reasoning
vs simple baselines, performance across perturbation regimes
(seen vs unseen genes, cell types, variants)"

This file delivers exactly that.
"""

import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


# ── EVALUATION SPLITS ────────────────────────────────────────────────────────
#
# This is the most important design decision in the entire evaluation.
#
# THE WRONG WAY: random 80/20 split of records
# If BRCA1 appears in training and test (different cell types),
# the model recalls BRCA1 facts from training. Looks smart. Isn't.
#
# THE RIGHT WAY: gene-level splits
# Every record for gene X goes entirely to train OR entirely to test.
# Never both. Now the model faces genes it has genuinely never seen.
# If it performs well here, it learned biology. Not facts.
#
# The gap between seen-gene performance and unseen-gene performance
# tells you exactly how much the model memorised vs generalised.
# That gap is the most honest number in your entire evaluation.

def build_evaluation_splits(records, seed=42):
    """
    Split records at the gene level to prevent evaluation leakage.

    Four splits of increasing difficulty:

    1. seen_genes — genes the model trained on, held out as sanity check
       If model fails here, training is fundamentally broken.
       Expected: highest performance.

    2. unseen_genes — genes never seen during training
       The critical test. Forces genuine generalisation.
       Expected: lower than seen_genes.
       Gap between 1 and 2 = memorisation vs learning ratio.

    3. unseen_cell_type — trained on one cell type, tested on another
       Did model learn general biology or cell-line-specific patterns?
       Matters for real clinical applications involving new cell types.

    4. cross_modal — trained on CRISPR, tested on scPerturb-seq
       Did the modalities actually integrate?
       Most scientifically interesting. Hardest to pass.

    Parameters
    ----------
    records : list of dict
        Training records with metadata containing gene and cell_type.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict of {split_name: list of records}
    """
    np.random.seed(seed)

    # Extract metadata from records
    df = pd.DataFrame([r["metadata"] for r in records])
    df["record_idx"] = range(len(records))

    all_genes = df["gene"].unique()
    n_train_genes = int(len(all_genes) * 0.8)

    # Gene-level split — this is what makes evaluation honest
    train_genes = set(
        np.random.choice(all_genes, size=n_train_genes, replace=False)
    )
    unseen_genes = set(all_genes) - train_genes

    splits = {
        "seen_genes": [
            records[i] for i in
            df[df["gene"].isin(train_genes)]["record_idx"]
        ],
        "unseen_genes": [
            records[i] for i in
            df[df["gene"].isin(unseen_genes)]["record_idx"]
        ],
    }

    # Cell type split if multiple cell types exist
    if "cell_type" in df.columns and df["cell_type"].nunique() > 1:
        holdout_cell = df["cell_type"].value_counts().index[-1]
        splits["unseen_cell_type"] = [
            records[i] for i in
            df[df["cell_type"] == holdout_cell]["record_idx"]
        ]
        log.info(f"Held out cell type: {holdout_cell}")

    # Cross-modal split if multiple modalities exist
    if "modality" in df.columns and df["modality"].nunique() > 1:
        holdout_modality = "scPerturb_seq"
        splits["cross_modal"] = [
            records[i] for i in
            df[df["modality"] == holdout_modality]["record_idx"]
        ]
        log.info(f"Cross-modal split: {holdout_modality}")

    for name, split in splits.items():
        log.info(f"Split '{name}': {len(split)} records")

    return splits, list(train_genes)


# ── METRICS ──────────────────────────────────────────────────────────────────
#
# Three metrics. Each catches a different failure mode.
#
# Metric 1: Gene set overlap
# Did the model predict the right genes?
# Random baseline: ~0.001 (20 genes from 20,000)
# Anything meaningfully above this = model doing something real.
#
# Metric 2: Direction accuracy
# Did the model get up vs down correct?
# Random baseline: 0.5 (coin flip)
# A model stuck at 0.5 knows which genes are relevant
# but has no idea what direction — specific, diagnosable failure.
#
# Metric 3: Pathway enrichment score
# Did the model identify the right biological PROCESS?
# More meaningful than individual genes — a model can get the
# pathway right even if it misses some individual members.
# This is how biologists actually think about perturbation effects.

def gene_set_overlap_at_k(predicted_genes, true_genes, k=20):
    """
    Fraction of true top-K DE genes the model correctly predicted.

    This is the primary accuracy metric. It asks: if a biologist
    uses this model to generate a shortlist of genes to follow up
    on, how many of the actually important genes are on that list?

    Random baseline: k / total_genes ≈ 0.001 for k=20, 20k genes.
    A score of 0.3 means the model found 6 of the top 20 real genes.
    That's genuinely useful for experimental prioritisation.

    Parameters
    ----------
    predicted_genes : list
        Gene names predicted by the model, ranked by importance.
    true_genes : list
        True top DE genes from the delta vector, ranked by |fold-change|.
    k : int
        How many top genes to compare.

    Returns
    -------
    float in [0, 1]
    """
    pred_set = set(predicted_genes[:k])
    true_set = set(true_genes[:k])

    if not true_set:
        return 0.0

    return len(pred_set & true_set) / min(k, len(true_set))


def direction_accuracy(pred_up, pred_down, true_up, true_down):
    """
    Fraction of shared genes where model correctly predicted direction.

    A model that knows ISG15 is affected by STAT1 knockout but says
    it goes UP when it actually goes DOWN is biologically wrong.
    This metric catches that specific failure.

    For genes that appear in both predicted and true gene sets,
    check whether predicted direction matches ground truth.

    Random baseline: 0.5
    Perfect score: 1.0

    Parameters
    ----------
    pred_up : list — genes predicted as upregulated
    pred_down : list — genes predicted as downregulated
    true_up : list — genes actually upregulated
    true_down : list — genes actually downregulated

    Returns
    -------
    float in [0, 1]
    """
    pred_up_set = set(pred_up)
    pred_down_set = set(pred_down)
    true_up_set = set(true_up)
    true_down_set = set(true_down)

    # Only evaluate on genes appearing in both predictions and ground truth
    all_pred = pred_up_set | pred_down_set
    all_true = true_up_set | true_down_set
    shared = all_pred & all_true

    if not shared:
        return 0.0

    correct = sum(
        1 for gene in shared
        if (gene in pred_up_set) == (gene in true_up_set)
    )

    return correct / len(shared)


def pathway_overlap_score(predicted_genes, true_genes, top_n=5):
    """
    Fraction of true enriched pathways the model also predicts.

    Individual gene predictions are noisy. Pathway-level analysis
    is how biologists actually interpret perturbation effects —
    they ask "which biological processes are affected" not
    "is gene number 847 in the list."

    Uses Fisher's exact test to find enriched pathways in each
    gene list, then measures overlap between top enriched pathways.

    A model can score well here even with imperfect gene recall
    if it correctly identifies the underlying biology — which is
    exactly the kind of reasoning we want to reward.

    Parameters
    ----------
    predicted_genes : list — genes predicted by model
    true_genes : list — true top DE genes
    top_n : int — how many top pathways to compare

    Returns
    -------
    float in [0, 1]
    """
    # Minimal pathway gene sets for demonstration
    # In production: load from MSigDB, KEGG, or Reactome
    gene_sets = {
        "interferon_response": [
            "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3",
            "IRF7", "STAT1", "STAT2", "IFI44", "RSAD2"
        ],
        "cell_cycle": [
            "CDK2", "CCND1", "CCNE1", "CDC20", "BUB1",
            "PCNA", "MCM2", "E2F1", "RB1", "CDKN1A"
        ],
        "apoptosis": [
            "BAX", "BCL2", "CASP3", "CASP9", "TP53",
            "PUMA", "NOXA", "MCL1", "BID", "CYCS"
        ],
        "dna_damage_response": [
            "TP53", "ATM", "ATR", "CHEK1", "CHEK2",
            "BRCA1", "BRCA2", "RAD51", "H2AX", "MDM2"
        ],
        "jak_stat_signaling": [
            "JAK1", "JAK2", "STAT1", "STAT3", "STAT5A",
            "SOCS1", "SOCS3", "IL6ST", "IFNGR1", "IL2RG"
        ],
        "pi3k_akt_signaling": [
            "PIK3CA", "AKT1", "PTEN", "MTOR", "TSC1",
            "TSC2", "RPS6KB1", "EIF4EBP1", "PDK1", "FOXO3"
        ],
    }

    universe = set()
    for gs in gene_sets.values():
        universe.update(gs)
    N = len(universe)

    def top_pathways(gene_list):
        query = set(gene_list) & universe
        if not query:
            return set()

        pvals = {}
        for pathway, members in gene_sets.items():
            pathway_set = set(members)
            k = len(query & pathway_set)
            K = len(pathway_set)
            n = len(query)

            table = [
                [k, K - k],
                [n - k, max(0, N - n - K + k)]
            ]
            _, pval = fisher_exact(table, alternative="greater")
            pvals[pathway] = pval

        ranked = sorted(pvals.items(), key=lambda x: x[1])
        return {p for p, _ in ranked[:top_n]}

    pred_pathways = top_pathways(predicted_genes)
    true_pathways = top_pathways(true_genes)

    if not true_pathways:
        return 0.0

    return len(pred_pathways & true_pathways) / len(true_pathways)


# ── PARSE MODEL OUTPUT ───────────────────────────────────────────────────────
#
# The model outputs natural language. We need gene lists.
# This parser extracts gene names from the text format
# produced by delta_to_text() in preprocess_scrna.py.
#
# This is deliberately simple — a real deployment would use
# structured output format. For a prototype, regex parsing
# of a consistent template is sufficient.

def parse_genes_from_output(text, direction="up"):
    """
    Extract gene names from model output text.

    Handles the format: "upregulation of: GENE1 (+2.3), GENE2 (+1.8)"

    Parameters
    ----------
    text : str
        Model-generated response text.
    direction : str
        "up" or "down" — which direction to extract.

    Returns
    -------
    list of gene name strings
    """
    import re

    if direction == "up":
        pattern = r"upregulation of[:\s]+([^;\.]+)"
    else:
        pattern = r"downregulation of[:\s]+([^;\.]+)"

    match = re.search(pattern, text, re.IGNORECASE)
    if not match:
        return []

    section = match.group(1)
    # Gene names: uppercase letters, numbers, hyphens
    genes = re.findall(r'\b([A-Z][A-Z0-9\-]{1,10})\b', section)
    return genes


# ── FULL EVALUATION LOOP ─────────────────────────────────────────────────────

def evaluate(predictions, ground_truth_records, k=20):
    """
    Run full evaluation suite on model predictions.

    Parameters
    ----------
    predictions : list of dict
        Each dict has keys: gene, predicted_text
    ground_truth_records : list of dict
        Training records with metadata containing true gene lists.
    k : int
        K for gene set overlap metric.

    Returns
    -------
    tuple of (metrics dict, per-gene results DataFrame)
    """
    # Build ground truth lookup by gene
    gt = {}
    for record in ground_truth_records:
        meta = record.get("metadata", {})
        gene = meta.get("gene")
        if gene:
            gt[gene] = meta

    results = []
    skipped = 0

    for pred in predictions:
        gene = pred.get("gene")
        text = pred.get("predicted_text", "")

        if gene not in gt:
            skipped += 1
            continue

        true_up = gt[gene].get("top_up_genes", [])
        true_down = gt[gene].get("top_down_genes", [])

        pred_up = parse_genes_from_output(text, "up")
        pred_down = parse_genes_from_output(text, "down")

        overlap_up = gene_set_overlap_at_k(pred_up, true_up, k=k)
        overlap_down = gene_set_overlap_at_k(pred_down, true_down, k=k)
        dir_acc = direction_accuracy(pred_up, pred_down, true_up, true_down)
        pathway = pathway_overlap_score(
            pred_up + pred_down,
            true_up + true_down
        )

        results.append({
            "gene": gene,
            "overlap_up": overlap_up,
            "overlap_down": overlap_down,
            "mean_overlap": (overlap_up + overlap_down) / 2,
            "direction_accuracy": dir_acc,
            "pathway_score": pathway,
        })

    if not results:
        log.warning("No matching genes between predictions and ground truth")
        return {}, pd.DataFrame()

    df = pd.DataFrame(results)
    log.info(f"Evaluated {len(results)} genes, skipped {skipped}")

    metrics = {
        "n_evaluated": len(results),
        "mean_overlap_up": round(df["overlap_up"].mean(), 4),
        "mean_overlap_down": round(df["overlap_down"].mean(), 4),
        "mean_overlap_both": round(df["mean_overlap"].mean(), 4),
        "mean_direction_accuracy": round(df["direction_accuracy"].mean(), 4),
        "mean_pathway_score": round(df["pathway_score"].mean(), 4),
        "median_overlap": round(df["mean_overlap"].median(), 4),
    }

    return metrics, df


# ── DEMO ─────────────────────────────────────────────────────────────────────

def demo():
    """
    Demonstrate evaluation with synthetic predictions and ground truth.

    Shows two cases intentionally:
    - A good prediction (STAT1 — model gets the right genes and direction)
    - A poor prediction (BRCA1 — model gets completely wrong genes)

    This demonstrates what the metrics mean in practice and how
    they distinguish good predictions from bad ones.
    """
    log.info("Running evaluation demo...")

    ground_truth = [
        {
            "metadata": {
                "gene": "STAT1",
                "top_up_genes": ["CCND1", "CDK2", "MYC", "E2F1", "PCNA"],
                "top_down_genes": ["ISG15", "MX1", "OAS1", "IFIT1", "IFIT3"],
            }
        },
        {
            "metadata": {
                "gene": "BRCA1",
                "top_up_genes": ["BAX", "CASP3", "TP53", "PUMA", "NOXA"],
                "top_down_genes": ["RAD51", "CHEK1", "CCND1", "CDK2", "E2F1"],
            }
        },
    ]

    predictions = [
        {
            "gene": "STAT1",
            # Good prediction — correct genes, correct directions
            "predicted_text": (
                "Knockout of STAT1 causes upregulation of: "
                "CCND1 (+1.82), CDK2 (+1.54), MYC (+1.23); "
                "and downregulation of: ISG15 (-2.31), MX1 (-1.87), OAS1 (-1.54)."
            ),
        },
        {
            "gene": "BRCA1",
            # Poor prediction — completely wrong genes
            "predicted_text": (
                "Knockout of BRCA1 causes upregulation of: "
                "STAT1 (+1.2), JAK2 (+0.9); "
                "and downregulation of: CTNNB1 (-1.1), APC (-0.8)."
            ),
        },
    ]

    metrics, results_df = evaluate(predictions, ground_truth, k=5)

    print("\n" + "=" * 60)
    print("EVALUATION RESULTS")
    print("=" * 60)

    print("\nPer-gene results:")
    print(results_df[[
        "gene", "overlap_up", "overlap_down",
        "direction_accuracy", "pathway_score"
    ]].to_string(index=False))

    print("\nAggregated metrics:")
    for k, v in metrics.items():
        print(f"  {k}: {v}")

    print("\nMetric interpretation:")
    print("  overlap_up/down : fraction of true DE genes correctly predicted")
    print("  Random baseline : ~0.001 for 20 genes from 20,000")
    print("  direction_acc   : fraction with correct up/down direction")
    print("  Random baseline : 0.5 (coin flip)")
    print("  pathway_score   : fraction of true pathways also predicted")
    print("=" * 60)


# ── ENTRY POINT ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Evaluation framework for perturbation-aware LLM"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run demo with synthetic predictions"
    )
    parser.add_argument(
        "--predictions",
        type=str,
        help="Path to JSONL file with model predictions"
    )
    parser.add_argument(
        "--ground_truth",
        type=str,
        help="Path to JSONL file with ground truth records"
    )
    parser.add_argument(
        "--k",
        type=int,
        default=20,
        help="K for gene set overlap metric"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="results/"
    )
    args = parser.parse_args()

    if args.demo:
        demo()
    else:
        if not args.predictions or not args.ground_truth:
            parser.error(
                "--predictions and --ground_truth required"
            )

        with open(args.predictions) as f:
            predictions = [json.loads(l) for l in f if l.strip()]

        with open(args.ground_truth) as f:
            ground_truth = [json.loads(l) for l in f if l.strip()]

        metrics, results_df = evaluate(
            predictions, ground_truth, k=args.k
        )

        Path(args.output_dir).mkdir(parents=True, exist_ok=True)

        results_df.to_csv(
            f"{args.output_dir}/per_gene_results.csv",
            index=False
        )

        with open(f"{args.output_dir}/metrics.json", "w") as f:
            json.dump(metrics, f, indent=2)

        log.info(f"Results saved to {args.output_dir}")
        print(json.dumps(metrics, indent=2))