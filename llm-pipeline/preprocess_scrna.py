"""
preprocess_scrna.py
====================
Preprocessing pipeline for scPerturb-seq data.

WHAT THIS FILE IS FOR:
----------------------
The project requires a "curated multimodal training corpus derived from 
the Perturbation Catalogue, including scPerturb-seq perturbation-response 
profiles." This file builds exactly that — it takes raw single-cell 
perturbation data and converts it into training examples an LLM can learn from.

BIOLOGICAL CONTEXT:
-------------------
scPerturb-seq (single-cell Perturbation sequencing) is an experimental 
technology that combines two things:

  1. CRISPR knockout — breaking one specific gene in each cell
  2. scRNA-seq — reading the activity level of every gene in that cell

The result: for each cell, you know WHICH gene was knocked out AND 
HOW the rest of the genome responded. This is the richest perturbation
data type available — it captures the full transcriptional consequence
of each perturbation at single-cell resolution.

The key output of this file is the PERTURBATION DELTA — the difference
in gene expression between perturbed cells and normal cells. This delta
is the biological ground truth that the LLM learns to predict.

PROJECT CONNECTION:
-------------------
- Supports: "Encoding perturbation experiments into LLM-friendly representations"
- Produces: scPerturb-seq perturbation-response profiles for the training corpus
- Enables: "What happens if gene X is knocked out in cell type Y?" queries
"""

import numpy as np
import pandas as pd
from scipy import sparse, stats
from pathlib import Path
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)



#
# Two filters catch most of the garbage:
#   1. Minimum genes detected — real cells express hundreds of genes.
#      If a "cell" only has 50 genes detected, it's probably an empty droplet.
#   2. Mitochondrial fraction — when a cell is dying, it loses cytoplasmic
#      RNA but mitochondrial RNA stays behind (mitochondria have their own
#      membrane). So a cell with >20% mitochondrial reads is probably dying.


def run_qc(adata, min_genes=200, max_mito_frac=0.2):
    """
    Filter out low-quality cells before any analysis.
    
    Parameters
    ----------
    adata : AnnData
        Raw single-cell data. Rows are cells, columns are genes.
        This is the standard format for single-cell data in Python.
    min_genes : int
        Minimum number of genes a cell must express to be kept.
        Cells below this are likely empty droplets.
    max_mito_frac : float
        Maximum allowed fraction of mitochondrial reads.
        Cells above this are likely dying.
    
    Returns
    -------
    AnnData with bad cells removed.
    """
    import anndata as ad
    import scanpy as sc

    log.info(f"Before QC: {adata.shape[0]} cells, {adata.shape[1]} genes")

    # Find mitochondrial genes — in human data they start with "MT-"
    # These are genes encoded in the mitochondrial genome, not the nucleus
    adata.var["is_mito"] = adata.var_names.str.startswith("MT-")

    # Calculate QC metrics — this adds columns to adata.obs like
    # n_genes_by_counts (how many genes detected) and 
    # pct_counts_is_mito (what % of reads are mitochondrial)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["is_mito"],
        inplace=True
    )

    # Apply the filters
    # We do these separately so we can see how many cells each filter removes
    before = adata.shape[0]
    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes].copy()
    after_gene_filter = adata.shape[0]

    adata = adata[adata.obs["pct_counts_is_mito"] <= (max_mito_frac * 100)].copy()
    after_mito_filter = adata.shape[0]

    log.info(f"Gene filter removed {before - after_gene_filter} cells")
    log.info(f"Mito filter removed {after_gene_filter - after_mito_filter} cells")
    log.info(f"After QC: {adata.shape[0]} cells remaining")

    return adata

def normalise(adata):
    """
    Library-size normalisation followed by log1p transform.

    Two steps:
      1. Scale every cell to 10,000 total counts — makes cells
         comparable regardless of how deeply they were sequenced
      2. Log1p transform — compresses large values, expands small ones
         so biological differences at low expression aren't drowned out

    We use log1p (log of count + 1) not log because some genes have
    zero counts and log(0) is undefined.
    """
    import scanpy as sc

    # Store raw counts before we touch them
    adata.layers["raw_counts"] = adata.X.copy()

    # Scale to 10,000 counts per cell
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform
    sc.pp.log1p(adata)

    # Save normalised values as a named layer
    adata.layers["log1p_norm"] = adata.X.copy()

    log.info("Normalisation complete: log1p normalised counts stored")
    return adata


def compute_pseudobulk(adata, perturbation_col, control_label):
    """
    Average expression across all cells sharing the same perturbation.
    
    Parameters
    ----------
    adata : AnnData
        Normalised single-cell data.
    perturbation_col : str
        Column in adata.obs that contains perturbation labels.
        e.g. "perturbation" where values are gene names like "STAT1", "BRCA1"
        or "non-targeting" for control cells.
    control_label : str
        The label used for control cells — cells where no real gene
        was knocked out. Typically "non-targeting".
    
    Returns
    -------
    pd.DataFrame of shape (n_perturbations, n_genes)
        Each row is the mean expression profile for one perturbation.
        Rows are perturbation labels, columns are gene names.
    """
    log.info("Computing pseudobulk profiles...")

    groups = adata.obs[perturbation_col].unique()
    gene_names = adata.var_names.tolist()
    pseudobulk = {}

    for group in groups:
        # Get all cells with this perturbation
        mask = adata.obs[perturbation_col] == group
        cells = adata[mask]

        # Get the expression matrix for these cells
        # AnnData can store matrices as sparse (mostly zeros) to save memory
        # We convert to dense (regular numpy array) for averaging
        X = cells.X
        if sparse.issparse(X):
            X = X.toarray()

        # Average across all cells — rows are cells, so axis=0 averages
        # across cells giving us one value per gene
        pseudobulk[group] = X.mean(axis=0)

        log.info(f"  {group}: averaged {mask.sum()} cells")

    # Build a clean DataFrame — perturbations as rows, genes as columns
    df = pd.DataFrame(pseudobulk, index=gene_names).T

    log.info(f"Pseudobulk matrix: {df.shape[0]} perturbations x {df.shape[1]} genes")
    return df


def compute_deltas(pseudobulk, control_label):
    """
    Compute perturbation effect as delta = perturbed - control.
    
    This is the ground truth signal the LLM learns to predict.
    Every training record, every evaluation metric, every biological
    interpretation in this project is built on these delta vectors.
    
    Parameters
    ----------
    pseudobulk : pd.DataFrame
        Output of compute_pseudobulk(). Rows are perturbations, 
        columns are genes.
    control_label : str
        The row label for control cells — e.g. "non-targeting".
        This row becomes the baseline everything is compared against.
    
    Returns
    -------
    pd.DataFrame of shape (n_perturbations - 1, n_genes)
        Same as pseudobulk but with control row removed and every
        value now representing CHANGE relative to control.
        
    Example
    -------
    If deltas.loc["STAT1", "ISG15"] = -2.3, that means:
    In cells where STAT1 was knocked out, ISG15 was expressed 2.3 units
    LOWER than in control cells. STAT1 normally activates ISG15 —
    remove STAT1, ISG15 goes quiet.
    """
    if control_label not in pseudobulk.index:
        available = pseudobulk.index.tolist()[:10]
        raise ValueError(
            f"Control label '{control_label}' not found. "
            f"Available labels: {available}"
        )

    # Extract the control row — this is our baseline
    # Every other perturbation gets compared against this
    control_expr = pseudobulk.loc[control_label]

    # Remove control row — we don't want a "control vs control = 0" 
    # training example, that teaches the model nothing
    perturbed = pseudobulk.drop(index=control_label)

    # Subtract control from every row simultaneously
    # axis=1 means subtract column-by-column (gene by gene)
    # so for each perturbation: delta[gene] = perturbed[gene] - control[gene]
    deltas = perturbed.subtract(control_expr, axis=1)

    log.info(f"Delta matrix: {deltas.shape[0]} perturbations x {deltas.shape[1]} genes")
    log.info(f"Delta range: [{deltas.values.min():.2f}, {deltas.values.max():.2f}]")

    return deltas



def get_top_de_genes(delta_row, n_top=50):
    """
    Extract the most strongly affected genes from a delta vector.
    
    We sort by magnitude of change and take the top n_top genes
    in each direction — the most upregulated and most downregulated.
    
    Parameters
    ----------
    delta_row : pd.Series
        One row from the delta matrix. Index is gene names,
        values are fold-changes relative to control.
    n_top : int
        How many genes to keep in each direction.
        50 is a reasonable default — captures the main signal
        without overwhelming the text representation.
    
    Returns
    -------
    dict with keys:
        top_up   : list of (gene, fold_change) — most upregulated
        top_down : list of (gene, fold_change) — most downregulated
    """
    sorted_delta = delta_row.sort_values(ascending=False)

    # Positive delta = upregulated (gene got more active after knockout)
    top_up = [
        (gene, round(fc, 3))
        for gene, fc in sorted_delta.head(n_top).items()
        if fc > 0
    ]

    # Negative delta = downregulated (gene got less active after knockout)
    # tail() gives smallest values (most negative) — reverse so strongest first
    top_down = [
        (gene, round(fc, 3))
        for gene, fc in sorted_delta.tail(n_top).items()
        if fc < 0
    ]
    top_down = list(reversed(top_down))

    return {"top_up": top_up, "top_down": top_down}


def delta_to_text(gene, cell_type, top_de, n_display=10):
    """
    Convert top DE genes into natural language for LLM training.
    
    This is the actual text the model learns to generate.
    Format is consistent across all training examples so the model
    learns a reliable pattern — gene names, direction, magnitude.
    
    The fold-change values are kept as numbers in the text because
    they carry directional and magnitude information. A model trained
    on enough examples should learn that +2.3 means strongly up
    and -0.1 means barely changed.
    
    Parameters
    ----------
    gene : str
        The knocked-out gene. e.g. "STAT1"
    cell_type : str
        Cell type the experiment was run in. e.g. "K562"
    top_de : dict
        Output of get_top_de_genes()
    n_display : int
        How many genes to include in the text summary.
        10 balances informativeness with sequence length.
    
    Returns
    -------
    str — natural language description of the perturbation response
    """
    up_genes = top_de["top_up"][:n_display]
    down_genes = top_de["top_down"][:n_display]

    up_str = ", ".join(
        [f"{g} ({fc:+.2f})" for g, fc in up_genes]
    ) if up_genes else "none detected"

    down_str = ", ".join(
        [f"{g} ({fc:+.2f})" for g, fc in down_genes]
    ) if down_genes else "none detected"

    return (
        f"Knockout of {gene} in {cell_type} cells causes "
        f"upregulation of: {up_str}; "
        f"and downregulation of: {down_str}."
    )


def build_training_record(gene, cell_type, perturbation_type, delta_row, n_top=50):
    """
    Build one complete instruction-tuning training example.
    
    This is the format the LLM sees during fine-tuning.
    Three parts:
      - instruction: the question being asked
      - input: the context/metadata
      - output: the answer the model should learn to generate
    
    During training: model sees all three parts.
    During inference: model gets instruction + input, generates output.
    
    The metadata field is NOT fed to the model — it's kept for
    evaluation so we can compare model predictions against the
    actual gene lists from the delta vector.
    
    Parameters
    ----------
    gene : str
        Knocked-out gene name.
    cell_type : str
        Cell type of the experiment.
    perturbation_type : str
        e.g. "CRISPR knockout", "CRISPR interference"
    delta_row : pd.Series
        One row from the delta matrix for this gene.
    n_top : int
        Number of top DE genes to include.
    
    Returns
    -------
    dict with keys: instruction, input, output, metadata
    """
    top_de = get_top_de_genes(delta_row, n_top=n_top)
    text_summary = delta_to_text(gene, cell_type, top_de)

    # Count significantly changed genes for the input context
    # Using 0.5 as threshold — smaller changes are often noise
    n_up = len([g for g, fc in top_de["top_up"] if fc > 0.5])
    n_down = len([g for g, fc in top_de["top_down"] if fc < -0.5])

    return {
        "instruction": (
            f"Predict the transcriptional response to {perturbation_type} "
            f"of gene {gene} in {cell_type} cells. Describe the key "
            f"upregulated and downregulated genes and their biological "
            f"significance."
        ),
        "input": (
            f"Gene: {gene}. "
            f"Perturbation type: {perturbation_type}. "
            f"Cell type: {cell_type}. "
            f"Significantly upregulated genes (|log2FC| > 0.5): {n_up}. "
            f"Significantly downregulated genes (|log2FC| > 0.5): {n_down}."
        ),
        "output": text_summary,

        # Kept for evaluation — not seen by the model during training
        "metadata": {
            "gene": gene,
            "cell_type": cell_type,
            "perturbation_type": perturbation_type,
            "top_up_genes": [g for g, _ in top_de["top_up"][:20]],
            "top_down_genes": [g for g, _ in top_de["top_down"][:20]],
            "n_sig_up": n_up,
            "n_sig_down": n_down,
        }
    }

    

def run_pipeline(
    input_path,
    perturbation_col,
    control_label,
    cell_type,
    output_path,
    min_genes=200,
    max_mito_frac=0.2,
    n_top_de=50
):
    """
    Full pipeline: raw h5ad → JSONL training records.
    
    Parameters
    ----------
    input_path : str
        Path to .h5ad file from scPerturb or Replogle et al.
    perturbation_col : str
        Column name in adata.obs containing perturbation labels.
    control_label : str
        Label for control cells in perturbation_col.
    cell_type : str
        Human-readable cell type name for training records.
    output_path : str
        Where to save the JSONL training records.
    """
    import anndata as ad

    log.info(f"Loading {input_path}")
    adata = ad.read_h5ad(input_path)

    adata = run_qc(adata, min_genes=min_genes, max_mito_frac=max_mito_frac)
    adata = normalise(adata)

    pseudobulk = compute_pseudobulk(adata, perturbation_col, control_label)
    deltas = compute_deltas(pseudobulk, control_label)

    log.info("Building training records...")
    records = []
    for gene in deltas.index:
        record = build_training_record(
            gene=gene,
            cell_type=cell_type,
            perturbation_type="CRISPR knockout",
            delta_row=deltas.loc[gene],
            n_top=n_top_de
        )
        records.append(record)

    # Save as JSONL — one training record per line
    # JSONL is standard for LLM training datasets
    import json
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for record in records:
            f.write(json.dumps(record) + "\n")

    log.info(f"Saved {len(records)} training records to {output_path}")
    return records


def demo():
    """
    Run the full pipeline on synthetic data — no real datasets needed.
    
    This exists so anyone can clone the repo and immediately verify
    the pipeline works without downloading gigabyte-scale datasets.
    Good software practice: always provide a runnable demo.
    """
    import anndata as ad

    log.info("Running demo with synthetic data...")
    np.random.seed(42)

    n_cells = 600
    n_genes = 500
    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]

    # Simulate 5 perturbations + 1 control
    perturbations = [
        "STAT1", "BRCA1", "TP53", "EGFR", "MYC", "non-targeting"
    ]
    cell_labels = np.random.choice(perturbations, size=n_cells)

    # Simulate count data — negative binomial is realistic for RNA counts
    counts = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)

    # Add a real biological signal for STAT1 knockout:
    # STAT1 is a transcription factor that activates interferon genes
    # Knocking it out should suppress interferon-stimulated genes
    stat1_mask = cell_labels == "STAT1"
    counts[stat1_mask, 10:20] *= 0.2   # suppress downstream targets
    counts[stat1_mask, 50:60] *= 3.0   # activate compensatory genes

    # Add MT- genes for QC demo
    mt_names = ["MT-ND1", "MT-ND2", "MT-CO1"] + gene_names[3:]

    import pandas as pd
    adata = ad.AnnData(
        X=sparse.csr_matrix(counts),
        obs=pd.DataFrame(
            {"perturbation": cell_labels},
            index=[f"cell_{i}" for i in range(n_cells)]
        ),
        var=pd.DataFrame(index=mt_names)
    )

    # Run pipeline with relaxed QC thresholds for synthetic data
    adata = run_qc(adata, min_genes=5, max_mito_frac=0.9)
    adata = normalise(adata)
    pseudobulk = compute_pseudobulk(adata, "perturbation", "non-targeting")
    deltas = compute_deltas(pseudobulk, "non-targeting")

    # Show one example training record
    example = build_training_record(
        gene="STAT1",
        cell_type="K562",
        perturbation_type="CRISPR knockout",
        delta_row=deltas.loc["STAT1"]
    )

    print("\n" + "=" * 60)
    print("EXAMPLE TRAINING RECORD")
    print("=" * 60)
    print(f"\nINSTRUCTION:\n{example['instruction']}")
    print(f"\nINPUT:\n{example['input']}")
    print(f"\nOUTPUT:\n{example['output']}")
    print(f"\nMETADATA (for evaluation, not fed to model):")
    print(f"  Top upregulated: {example['metadata']['top_up_genes'][:5]}")
    print(f"  Top downregulated: {example['metadata']['top_down_genes'][:5]}")
    print("=" * 60)

    return deltas, example



if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="scPerturb-seq preprocessing pipeline"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run on synthetic data — no real data needed"
    )
    parser.add_argument("--input", type=str, help="Path to .h5ad file")
    parser.add_argument("--perturbation_col", type=str, default="perturbation")
    parser.add_argument("--control_label", type=str, default="non-targeting")
    parser.add_argument("--cell_type", type=str, default="unknown")
    parser.add_argument("--output", type=str, default="output/scrna_records.jsonl")
    parser.add_argument("--n_top_de", type=int, default=50)
    args = parser.parse_args()

    if args.demo:
        demo()
    else:
        if not args.input:
            parser.error("--input required unless --demo is specified")
        run_pipeline(
            input_path=args.input,
            perturbation_col=args.perturbation_col,
            control_label=args.control_label,
            cell_type=args.cell_type,
            output_path=args.output,
            n_top_de=args.n_top_de
        )
