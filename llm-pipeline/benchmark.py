import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


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

    df = pd.DataFrame([r["metadata"] for r in records])
    df["record_idx"] = range(len(records))

    all_genes = df["gene"].unique()
    n_train_genes = int(len(all_genes) * 0.8)

    # Gene-level split
    train_genes = set(np.random.choice(all_genes, size=n_train_genes, replace=False))
    unseen_genes = set(all_genes) - train_genes

    splits = {
        "seen_genes": [
            records[i] for i in df[df["gene"].isin(train_genes)]["record_idx"]
        ],
        "unseen_genes": [
            records[i] for i in df[df["gene"].isin(unseen_genes)]["record_idx"]
        ],
    }

    if "cell_type" in df.columns and df["cell_type"].nunique() > 1:
        holdout_cell = df["cell_type"].value_counts().index[-1]
        splits["unseen_cell_type"] = [
            records[i] for i in df[df["cell_type"] == holdout_cell]["record_idx"]
        ]
        log.info(f"Held out cell type: {holdout_cell}")

    if "modality" in df.columns and df["modality"].nunique() > 1:
        holdout_modality = "scPerturb-seq"
        splits["cross_modal"] = [
            records[i] for i in df[df["modality"] == holdout_modality]["record_idx"]
        ]
        log.info(f"Cross-modal split: {holdout_modality}")

    for name, split in splits.items():
        log.info(f"Split '{name}': {len(split)} records")

    return splits, list(train_genes)


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

    all_pred = pred_up_set | pred_down_set
    all_true = true_up_set | true_down_set
    shared = all_pred & all_true

    if not shared:
        return 0.0

    correct = sum(
        1 for gene in shared if (gene in pred_up_set) == (gene in true_up_set)
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

    gene_sets = {
        "interferon_response": [
            "ISG15",
            "MX1",
            "OAS1",
            "IFIT1",
            "IFIT3",
            "IRF7",
            "STAT1",
            "STAT2",
            "IFI44",
            "RSAD2",
        ],
        "cell_cycle": [
            "CDK2",
            "CCND1",
            "CCNE1",
            "CDC20",
            "BUB1",
            "PCNA",
            "MCM2",
            "E2F1",
            "RB1",
            "CDKN1A",
        ],
        "apoptosis": [
            "BAX",
            "BCL2",
            "CASP3",
            "CASP9",
            "TP53",
            "PUMA",
            "NOXA",
            "MCL1",
            "BID",
            "CYCS",
        ],
        "dna_damage_response": [
            "TP53",
            "ATM",
            "ATR",
            "CHEK1",
            "CHEK2",
            "BRCA1",
            "BRCA2",
            "RAD51",
            "H2AX",
            "MDM2",
        ],
        "jak_stat_signaling": [
            "JAK1",
            "JAK2",
            "STAT1",
            "STAT3",
            "STAT5A",
            "SOCS1",
            "SOCS3",
            "IL6ST",
            "IFNGR1",
            "IL2RG",
        ],
        "pi3k_akt_signaling": [
            "PIK3CA",
            "AKT1",
            "PTEN",
            "MTOR",
            "TSC1",
            "TSC2",
            "RPS6KB1",
            "EIF4EBP1",
            "PDK1",
            "FOXO3",
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

            table = [[k, K - k], [n - k, max(0, N - n - K + k)]]
            _, pval = fisher_exact(table, alternative="greater")
            pvals[pathway] = pval

        ranked = sorted(pvals.items(), key=lambda x: x[1])
        return {p for p, _ in ranked[:top_n]}

    pred_pathways = top_pathways(predicted_genes)
    true_pathways = top_pathways(true_genes)

    if not true_pathways:
        return 0.0

    return len(pred_pathways & true_pathways) / len(true_pathways)


def parse_genes_from_output(text, direction="up"):
    """
    Extract gene names from model output text.

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
        pattern = r"upregulation of[:\s]+([^;]+?)(?:;|$)"
    else:
        pattern = r"downregulation of[:\s]+([^;]+?)(?:;|\.$|$)"

    match = re.search(pattern, text, re.IGNORECASE)
    if not match:
        return []

    section = match.group(1)
    genes = re.findall(r"\b([A-Z][A-Z0-9\-]{1,10})\b", section)
    return genes


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
        pathway = pathway_overlap_score(pred_up + pred_down, true_up + true_down)

        results.append(
            {
                "gene": gene,
                "overlap_up": overlap_up,
                "overlap_down": overlap_down,
                "mean_overlap": (overlap_up + overlap_down) / 2,
                "direction_accuracy": dir_acc,
                "pathway_score": pathway,
            }
        )

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


if __name__ == "__main__":
    print("benchmark.py can be imported to use evaluation functions.")
