"""
preprocess_crispr.py
=====================
Preprocessing pipeline for CRISPR pooled fitness screen data.

BIOLOGICAL CONTEXT:
-------------------
A CRISPR pooled fitness screen works like this:
You take millions of cells. You use CRISPR to knock out a different
gene in each cell. Then you apply selective pressure — a drug, a
disease condition, or just normal growth. After several days you
sequence what's left and count which knockouts survived.

Genes whose knockout caused cells to die = essential genes
Genes whose knockout helped cells survive = anti-essential (tumour suppressors)
Genes whose knockout did nothing = neutral

The output is a log fold change (LFC) per gene — how much did cells
with that knockout increase or decrease in the population over time.
Negative LFC = cells dropped out = gene was essential for survival.
Positive LFC = cells enriched = gene was suppressing growth.

This is processed by MAGeCK — the standard tool for CRISPR screen
analysis. MAGeCK takes raw guide RNA counts and outputs gene-level
scores with statistical significance.

PROJECT CONNECTION:
-------------------
CRISPR screen data is one of the three modalities in the training
corpus. It answers a different question than scPerturb-seq:
not "what transcriptional changes happen" but "does this gene
matter for survival at all." Together they give the model both
the fitness consequence AND the molecular mechanism.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import logging
import json

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


# ── STEP 1: LOAD DATA ────────────────────────────────────────────────────────
#
# MAGeCK outputs a tab-separated file with one row per gene.
# Key columns we care about:
#   neg|lfc  — log fold change in negative selection (dropout)
#   neg|fdr  — false discovery rate for negative selection
#   pos|lfc  — log fold change in positive selection (enrichment)
#   pos|fdr  — false discovery rate for positive selection
#
# FDR (false discovery rate) is the statistical significance measure.
# FDR < 0.05 means less than 5% chance this is a false positive.
# This is the standard threshold in genomics.


def load_mageck_output(filepath):
    """
    Load MAGeCK gene_summary.txt output.

    MAGeCK is the standard tool for analysing pooled CRISPR screens.
    It takes raw guide RNA counts and outputs gene-level statistics.

    Parameters
    ----------
    filepath : str
        Path to MAGeCK gene_summary.txt file.

    Returns
    -------
    pd.DataFrame with standardised column names.
    """
    log.info(f"Loading MAGeCK output from {filepath}")
    df = pd.read_csv(filepath, sep="\t")

    # Standardise column names — MAGeCK uses pipe notation (neg|lfc)
    # which is awkward to work with. Rename to clean snake_case.
    rename_map = {
        "Gene": "gene",
        "neg|lfc": "neg_lfc",
        "neg|fdr": "neg_fdr",
        "pos|lfc": "pos_lfc",
        "pos|fdr": "pos_fdr",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    log.info(f"Loaded {len(df)} genes")
    return df


# ── STEP 2: NORMALISE ────────────────────────────────────────────────────────
#
# Raw LFC values are not comparable across screens.
# One screen might have a mean LFC of -0.5, another -2.0,
# depending on sequencing depth and library complexity.
#
# Z-score normalisation fixes this:
# subtract the mean, divide by standard deviation.
# Now every screen has mean=0 and std=1.
# A z-score of -3 means "3 standard deviations below average"
# regardless of which screen it came from.


def normalise_lfc(df, lfc_col="neg_lfc"):
    """
    Z-score normalise LFC values within a screen.

    Essential for cross-screen comparisons — raw LFCs are not
    comparable between screens with different sequencing depths.

    Parameters
    ----------
    df : pd.DataFrame
        MAGeCK output with LFC column.
    lfc_col : str
        Column to normalise.

    Returns
    -------
    pd.DataFrame with added z-score column.
    """
    lfc = df[lfc_col].copy()
    df[f"{lfc_col}_zscore"] = stats.zscore(lfc, nan_policy="omit")

    log.info(f"Normalised {lfc_col}: mean={lfc.mean():.3f}, std={lfc.std():.3f}")
    return df


# ── STEP 3: CLASSIFY FITNESS EFFECT ─────────────────────────────────────────
#
# Three categories:
#
# essential — strong negative LFC, significant FDR
#   Knocking out this gene kills or severely impairs the cell.
#   These are genes the cell cannot live without.
#   Example: core DNA replication machinery, ribosomal genes.
#
# anti_essential — strong positive LFC, significant FDR
#   Knocking out this gene HELPS the cell survive.
#   These are often tumour suppressors — genes that normally
#   put the brakes on cell growth. Cancer cells benefit from
#   losing them. Examples: TP53, PTEN, RB1.
#
# neutral — no significant effect
#   The cell doesn't care whether this gene is present
#   under these specific conditions.


def classify_fitness_effect(
    df, lfc_col="neg_lfc", fdr_col="neg_fdr", fdr_threshold=0.05, lfc_threshold=0.5
):
    """
    Classify each gene's fitness effect into three categories.

    This classification becomes the label in CRISPR training records.
    The model learns to predict which category a gene falls into
    and why — connecting fitness phenotype to biological function.

    Parameters
    ----------
    df : pd.DataFrame
        Normalised MAGeCK output.
    lfc_col : str
        LFC column to use for classification.
    fdr_col : str
        FDR column for significance filtering.
    fdr_threshold : float
        Significance cutoff. 0.05 is standard in genomics.
    lfc_threshold : float
        Minimum absolute LFC to call an effect real.
        Filters out statistically significant but tiny effects.

    Returns
    -------
    pd.DataFrame with added fitness_class column.
    """
    conditions = [
        (df[fdr_col] < fdr_threshold) & (df[lfc_col] < -lfc_threshold),
        (df[fdr_col] < fdr_threshold) & (df[lfc_col] > lfc_threshold),
    ]
    choices = ["essential", "anti_essential"]

    df["fitness_class"] = np.select(conditions, choices, default="neutral")

    counts = df["fitness_class"].value_counts()
    log.info(f"Fitness classification: {counts.to_dict()}")
    return df


# ── STEP 4: BUILD TRAINING RECORDS ───────────────────────────────────────────
#
# Same structure as scPerturb-seq records — instruction/input/output.
# This consistency matters: the model sees the same format regardless
# of which modality the data came from. That's what makes cross-modal
# reasoning possible — unified representation.


def fitness_class_to_text(gene, fitness_class, lfc, cell_line, condition):
    """
    Convert CRISPR fitness classification to natural language.

    The text explains not just the result but the biological
    interpretation — connecting the fitness phenotype to what
    it means mechanistically. This is what teaches the model
    to reason rather than just classify.
    """
    descriptions = {
        "essential": (
            f"{gene} is essential for survival of {cell_line} "
            f"under {condition} conditions (LFC: {lfc:.2f}). "
            f"Knockout causes significant cell depletion, indicating "
            f"this gene is required for cell viability or proliferation. "
            f"This is consistent with a role in core cellular processes "
            f"such as DNA replication, transcription, or metabolism."
        ),
        "anti_essential": (
            f"{gene} acts as a fitness suppressor in {cell_line} "
            f"under {condition} conditions (LFC: {lfc:.2f}). "
            f"Knockout causes cell enrichment, meaning cells without "
            f"this gene grow faster. This pattern is consistent with "
            f"tumour suppressor activity — the gene normally restrains "
            f"cell proliferation."
        ),
        "neutral": (
            f"{gene} shows no significant fitness effect in {cell_line} "
            f"cells under {condition} conditions (LFC: {lfc:.2f}). "
            f"Knockout does not substantially alter cell survival or "
            f"proliferation under these specific conditions, though "
            f"effects may emerge under different stresses."
        ),
    }
    return descriptions.get(fitness_class, "Fitness effect unknown.")


def build_crispr_training_record(
    row, screen_id, cell_line, condition, lfc_col="neg_lfc", fdr_col="neg_fdr"
):
    """
    Build one instruction-tuning record from a CRISPR screen gene.

    Parameters
    ----------
    row : pd.Series
        One row from the classified MAGeCK DataFrame.
    screen_id : str
        Identifier for this screen — used for traceability.
    cell_line : str
        Cell line used in the screen. e.g. "K562", "HeLa"
    condition : str
        Experimental condition. e.g. "standard_growth", "drug_treatment"

    Returns
    -------
    dict with instruction, input, output, metadata
    """
    gene = row["gene"]
    lfc = row.get(lfc_col, 0.0)
    fdr = row.get(fdr_col, 1.0)
    fitness_class = row.get("fitness_class", "neutral")

    return {
        "instruction": (
            f"What is the fitness effect of knocking out gene {gene} "
            f"in {cell_line} cells under {condition} conditions? "
            f"Describe the phenotype and its biological interpretation."
        ),
        "input": (
            f"Gene: {gene}. "
            f"Cell line: {cell_line}. "
            f"Screen condition: {condition}. "
            f"Screen ID: {screen_id}. "
            f"Modality: CRISPR pooled fitness screen."
        ),
        "output": fitness_class_to_text(gene, fitness_class, lfc, cell_line, condition),
        "metadata": {
            "gene": gene,
            "screen_id": screen_id,
            "cell_line": cell_line,
            "condition": condition,
            "lfc": float(lfc),
            "fdr": float(fdr),
            "fitness_class": fitness_class,
            "modality": "CRISPR_screen",
        },
    }


# ── FULL PIPELINE ────────────────────────────────────────────────────────────


def run_pipeline(
    input_path,
    screen_id,
    cell_line,
    condition,
    output_path,
    fdr_threshold=0.05,
    lfc_threshold=0.5,
):
    """
    Full pipeline: MAGeCK output → JSONL training records.
    """
    df = load_mageck_output(input_path)
    df = normalise_lfc(df)
    df = classify_fitness_effect(
        df, fdr_threshold=fdr_threshold, lfc_threshold=lfc_threshold
    )

    records = []
    for _, row in df.iterrows():
        record = build_crispr_training_record(row, screen_id, cell_line, condition)
        records.append(record)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for record in records:
            f.write(json.dumps(record) + "\n")

    log.info(f"Saved {len(records)} records to {output_path}")
    return records


# ── DEMO ─────────────────────────────────────────────────────────────────────


def demo():
    """
    Run on synthetic data — no real datasets needed.

    Simulates a realistic CRISPR screen with known essential genes,
    anti-essential genes, and neutral genes. Verifies the full
    pipeline works end to end.
    """
    log.info("Running CRISPR demo with synthetic data...")
    np.random.seed(0)

    n_genes = 300
    genes = [f"GENE_{i:04d}" for i in range(n_genes)]

    # Simulate realistic LFC distribution
    # Most genes are neutral — LFC around zero
    neg_lfc = np.random.normal(0, 0.5, n_genes)

    # ~10% essential — strong negative LFC
    essential_idx = np.random.choice(n_genes, size=30, replace=False)
    neg_lfc[essential_idx] = np.random.normal(-2.5, 0.4, 30)

    # ~5% anti-essential — strong positive LFC
    remaining = [i for i in range(n_genes) if i not in essential_idx]
    anti_idx = np.random.choice(remaining, size=15, replace=False)
    neg_lfc[anti_idx] = np.random.normal(2.0, 0.4, 15)

    # FDR: essential and anti-essential genes get low FDR
    neg_fdr = np.ones(n_genes) * 0.5
    neg_fdr[essential_idx] = np.random.uniform(0.001, 0.04, 30)
    neg_fdr[anti_idx] = np.random.uniform(0.001, 0.04, 15)

    df = pd.DataFrame({"gene": genes, "neg_lfc": neg_lfc, "neg_fdr": neg_fdr})

    df = normalise_lfc(df)
    df = classify_fitness_effect(df)

    # Build example records
    records = []
    for _, row in df.head(5).iterrows():
        record = build_crispr_training_record(
            row, screen_id="demo_screen", cell_line="K562", condition="standard_growth"
        )
        records.append(record)

    # Show one example
    r = records[0]
    print("\n" + "=" * 60)
    print("EXAMPLE CRISPR TRAINING RECORD")
    print("=" * 60)
    print(f"\nINSTRUCTION:\n{r['instruction']}")
    print(f"\nINPUT:\n{r['input']}")
    print(f"\nOUTPUT:\n{r['output']}")
    print(f"\nMETADATA:\n{r['metadata']}")
    print("=" * 60)

    counts = df["fitness_class"].value_counts()
    print(f"\nDataset summary: {counts.to_dict()}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="CRISPR screen preprocessing pipeline")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run on synthetic data — no real data needed",
    )
    parser.add_argument("--input", type=str)
    parser.add_argument("--screen_id", type=str, default="screen_001")
    parser.add_argument("--cell_line", type=str, default="K562")
    parser.add_argument("--condition", type=str, default="standard_growth")
    parser.add_argument("--output", type=str, default="output/crispr_records.jsonl")
    parser.add_argument("--fdr_threshold", type=float, default=0.05)
    parser.add_argument("--lfc_threshold", type=float, default=0.5)
    args = parser.parse_args()

    if args.demo:
        demo()
    else:
        if not args.input:
            parser.error("--input required unless --demo is specified")
        run_pipeline(
            input_path=args.input,
            screen_id=args.screen_id,
            cell_line=args.cell_line,
            condition=args.condition,
            output_path=args.output,
            fdr_threshold=args.fdr_threshold,
            lfc_threshold=args.lfc_threshold,
        )
