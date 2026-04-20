import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import sys
import scipy.sparse as sp
import gc
import json
import textwrap

# For nice Jupyter rendering
try:
    from IPython.display import display
except ImportError:
    display = print


# ==============================================================================
# 1. FUNCTIONS AND UTILITIES
# ==============================================================================


def clean_barcodes(adata):
    """Removes common suffixes from barcodes (e.g., '-1') for alignment."""
    adata.obs_names = adata.obs_names.str.split("-").str[0]
    if adata.obs_names.duplicated().any():
        print(f"  Warning: Found {adata.obs_names.duplicated().sum()} duplicate barcodes. Suffixing...")
        adata.obs_names_make_unique()
    return adata


def get_raw_counts(adata):
    """Safely extracts raw counts."""
    X = adata.X
    if hasattr(X, "data"):
        sample = X.data[:2000]
    else:
        sample = X.flatten()[:2000]
    return np.all(np.equal(np.mod(sample, 1), 0))


def preprocess_adata(adata, name, target_sum=1e4, n_top_genes=2000):
    """Standardized preprocessing: Raw -> Norm -> Log -> HVG -> Scale -> PCA."""
    print(f"Starting standardized preprocessing for {name}...")
    is_raw = get_raw_counts(adata)
    if not is_raw:
        print(f"  [{name}] WARNING: Data does not appear to be raw counts.")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat", subset=False)
    return adata


def plot_scatter_comparison(
    df, x_col, y_col, title, xlabel, ylabel, filename, log_scale=False, description=""
):
    """
    Standardized scatter plot with professional layout:
    - Bold Title
    - Subtitle with stats
    - Wrapped Description at the bottom
    """
    fig, ax = plt.subplots(figsize=(10, 11))

    plot_df = df.copy().dropna(subset=[x_col, y_col])
    
    # Calculate deviation stats (within 10%)
    diff = np.abs(plot_df[x_col] - plot_df[y_col])
    # Use max to avoid division by zero
    denom = plot_df[[x_col, y_col]].mean(axis=1).replace(0, 1)
    rel_diff = diff / denom
    pct_deviant = (rel_diff > 0.1).mean() * 100

    if log_scale:
        plot_df[x_col] = plot_df[x_col] + 1
        plot_df[y_col] = plot_df[y_col] + 1

    # Main plot
    if len(plot_df) > 5000:
        ax.scatter(plot_df[x_col], plot_df[y_col], alpha=0.05, s=1, color='teal', rasterized=True)
    else:
        sns.scatterplot(data=plot_df, x=x_col, y=y_col, alpha=0.3, s=10, ax=ax)

    # Identity line
    min_val = min(plot_df[x_col].min(), plot_df[y_col].min())
    max_val = max(plot_df[x_col].max(), plot_df[y_col].max())
    ax.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="Identity")

    # Stats
    pearson, _ = stats.pearsonr(plot_df[x_col], plot_df[y_col])
    spearman, _ = stats.spearmanr(plot_df[x_col], plot_df[y_col])
    
    # Titling
    ax.set_title(title, fontsize=18, fontweight='bold', pad=35)
    ax.text(0.5, 1.02, f"Pearson r = {pearson:.4f} | Spearman rho = {spearman:.4f} | Deviants (>10%): {pct_deviant:.1f}%", 
            transform=ax.transAxes, ha='center', va='bottom', fontsize=11, style='italic')

    ax.set_xlabel(f"{xlabel} {'(+1 for log)' if log_scale else ''}", fontsize=12)
    ax.set_ylabel(f"{ylabel} {'(+1 for log)' if log_scale else ''}", fontsize=12)
    
    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.legend(loc='upper left')

    # Description at the bottom
    wrapped_desc = "\n".join(textwrap.wrap(description, width=100))
    fig.text(0.5, 0.02, wrapped_desc, ha='center', va='bottom', fontsize=10, linespacing=1.4)

    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    plt.savefig(filename, dpi=300)
    plt.show()
    plt.close()

    return {"pearson": pearson, "spearman": spearman, "pct_deviant": pct_deviant}


def call_guides(adata, count_threshold=5):
    """Dual-guide caller."""
    if "guides" not in adata.obsm:
        return None
    guide_matrix = adata.obsm["guides"]
    if sp.issparse(guide_matrix):
        guide_matrix = guide_matrix.tocsr()
    
    guide_names = np.array(adata.uns.get("guide_names", [f"guide_{i}" for i in range(guide_matrix.shape[1])]))
    calls = []
    for i in range(guide_matrix.shape[0]):
        row = np.array(guide_matrix[i].toarray()).flatten()
        top_idx = np.where(row >= count_threshold)[0]
        if len(top_idx) == 0:
            calls.append("None")
        else:
            top_idx = top_idx[np.argsort(row[top_idx])[::-1][:2]]
            names = sorted(guide_names[top_idx])
            calls.append("|".join(names))
    return pd.Series(calls, index=adata.obs_names)


def compare_perturbations(adata_cur, adata_rep, common_cells):
    """Compares perturbation assignments with deep mismatch diagnostics."""
    print("Comparing perturbation assignments...")

    cur_pert_col = next((c for c in ["sgID_AB", "perturbation"] if c in adata_cur.obs.columns), None)
    rep_labels = call_guides(adata_rep) if "guides" in adata_rep.obsm else None

    if cur_pert_col and rep_labels is not None:
        p_cur = adata_cur.obs.loc[common_cells, cur_pert_col].astype(str)
        p_rep = rep_labels.loc[common_cells].astype(str)

        normalize = lambda s: "|".join(sorted([x.strip() for x in s.split("|")])) if "|" in s else s.strip()
        p_cur, p_rep = p_cur.apply(normalize), p_rep.apply(normalize)
        
        # Deep Diagnostics
        has_guide_idx = np.where(p_rep != "None")[0]
        if len(has_guide_idx) > 0:
            print(f"  Mismatched cells diagnostics (Cells where Reprocessed has guides):")
            for idx in has_guide_idx[:10]:
                print(f"    {common_cells[idx]}: Curated='{p_cur.iloc[idx]}' vs Rep='{p_rep.iloc[idx]}'")

        accuracy = (p_cur == p_rep).mean()
        overlap_df = pd.DataFrame({"Curated": p_cur, "Reprocessed": p_rep})
        top_perts = p_cur.value_counts().head(20).index
        sub_df = overlap_df[overlap_df["Curated"].isin(top_perts)]
        ct = pd.crosstab(sub_df["Curated"], sub_df["Reprocessed"])

        plt.figure(figsize=(15, 13))
        sns.heatmap(ct, annot=False, cmap="YlGnBu")
        plt.suptitle("Perturbation Confusion Matrix (Top 20)", fontsize=18, fontweight='bold', y=0.98)
        
        desc = ("This heatmap compares guide assignments between the original study (Y) and our reprocessed pipeline (X). "
                "A strong diagonal indicates consistent guide recovery and labeling across pipelines. "
                "Mismatches may suggest barcode collisions or thresholding differences.")
        wrapped_desc = "\n".join(textwrap.wrap(desc, width=110))
        plt.title(f"Overall Match: {accuracy:.4%}\n{wrapped_desc}", 
                  fontsize=10, pad=15, style='italic', loc='center')
        
        plt.xticks(rotation=45, ha='right', fontsize=7)
        plt.yticks(fontsize=7)
        plt.tight_layout(rect=[0, 0.05, 1, 0.94])
        plt.savefig("comparison_results/perturbation_confusion_matrix.png")
        plt.show()
        plt.close()

        return {"accuracy": accuracy, "n_rep_with_guides": len(has_guide_idx)}
    return None


# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================
os.makedirs("comparison_results", exist_ok=True)
print("Loading and aligning...")
adata_cur = clean_barcodes(sc.read_h5ad("GSE264667_jurkat_raw_singlecell_01.h5ad"))
adata_rep = clean_barcodes(sc.read_h5ad("experiment_final.h5ad"))

if adata_rep.var_names.str.contains(r"\.").any():
    adata_rep.var_names = adata_rep.var_names.str.split(".").str[0]
    adata_rep.var_names_make_unique()

common_cells = np.intersect1d(adata_cur.obs_names, adata_rep.obs_names)
common_genes = np.intersect1d(adata_cur.var_names, adata_rep.var_names)

print(f"  Common cells: {len(common_cells)}")
print(f"  Common genes: {len(common_genes)}")

cur_sub = preprocess_adata(adata_cur[common_cells, common_genes].copy(), "Curated")
rep_sub = preprocess_adata(adata_rep[common_cells, common_genes].copy(), "Reprocessed")

results_summary = {"cell_metrics": {}, "gene_metrics": {}}

# --- Overlap Distributions (3 Histograms) ---
print("Generating distribution histograms...")
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

# 1. Cell UMI Counts distribution
rep_total_counts = np.array(adata_rep.X.sum(axis=1)).flatten()
common_cell_mask = adata_rep.obs_names.isin(common_cells)
sns.histplot(rep_total_counts[common_cell_mask], label="Common Cells", color="blue", alpha=0.5, log_scale=True, ax=ax1)
sns.histplot(rep_total_counts[~common_cell_mask], label="Reprocessed-only", color="orange", alpha=0.5, log_scale=True, ax=ax1)
ax1.set_title("Total UMIs per Cell"); ax1.legend()

# 2. Gene n_cells distribution (how many cells have each gene)
rep_gene_sums = np.array((adata_rep.X > 0).sum(axis=0)).flatten()
common_gene_mask = adata_rep.var_names.isin(common_genes)
sns.histplot(rep_gene_sums[common_gene_mask], label="Common Genes", color="blue", alpha=0.5, log_scale=True, ax=ax2)
sns.histplot(rep_gene_sums[~common_gene_mask], label="Reprocessed-only", color="orange", alpha=0.5, log_scale=True, ax=ax2)
ax2.set_title("Detection Freq per Gene"); ax2.legend()

# 3. Distinct Genes per Cell distribution
rep_n_genes = np.array((adata_rep.X > 0).sum(axis=1)).flatten()
sns.histplot(rep_n_genes[common_cell_mask], label="Common Cells", color="blue", alpha=0.5, log_scale=True, ax=ax3)
sns.histplot(rep_n_genes[~common_cell_mask], label="Reprocessed-only", color="orange", alpha=0.5, log_scale=True, ax=ax3)
ax3.set_title("Distinct Genes per Cell"); ax3.legend()

plt.suptitle("Comparison of Global Distribution Overlaps", fontsize=16, fontweight='bold')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("comparison_results/overlap_distributions.png")
plt.show()
plt.close()

# --- Cell-wise Scatter Plots ---
results_summary["cell_metrics"]["total_counts"] = plot_scatter_comparison(
    pd.DataFrame({"cur": cur_sub.obs["total_counts"], "rep": rep_sub.obs["total_counts"]}),
    "cur", "rep", "Total UMI Counts", "Curated", "Reprocessed", "comparison_results/counts_comparison.png", 
    log_scale=True, 
    description="Correlation of total UMI counts per cell. Kallisto+Bustools (Reprocessed) generally identifies more UMIs per cell, especially in the high-sensitivity regime, likely due to transcript-level mapping vs gene-level quantification."
)

results_summary["cell_metrics"]["n_genes"] = plot_scatter_comparison(
    pd.DataFrame({"cur": cur_sub.obs["n_genes_by_counts"], "rep": rep_sub.obs["n_genes_by_counts"]}),
    "cur", "rep", "Number of Detected Genes", "Curated", "Reprocessed", "comparison_results/genes_comparison.png",
    description="Correlation of unique genes detected per cell. The reprocessed pipeline recovers significantly more unique transcripts in a subset of cells, indicating higher detection sensitivity."
)

# --- Gene-wise Scatter Plots (Restored) ---
gene_metrics_df = pd.DataFrame({
    "mean_cur": cur_sub.var["mean_counts"],
    "mean_rep": rep_sub.var["mean_counts"],
    "dropout_cur": 100 - (np.array((cur_sub.layers["counts"] > 0).sum(axis=0)).flatten() / cur_sub.n_obs * 100),
    "dropout_rep": 100 - (np.array((rep_sub.layers["counts"] > 0).sum(axis=0)).flatten() / rep_sub.n_obs * 100)
})

results_summary["gene_metrics"]["mean_expression"] = plot_scatter_comparison(
    gene_metrics_df, "mean_cur", "mean_rep", "Mean Gene Expression", "Curated", "Reprocessed", "comparison_results/gene_expression_mean.png",
    log_scale=True,
    description="Average expression level for each shared gene. High correlation confirms that the biological signal is preserved across quantification methods."
)

results_summary["gene_metrics"]["dropout_rate"] = plot_scatter_comparison(
    gene_metrics_df, "dropout_cur", "dropout_rep", "Gene Dropout Rate (%)", "Curated", "Reprocessed", "comparison_results/sparsity_comparison.png",
    description="Percentage of cells where a gene has zero counts. Lower values in Reprocessed indicate higher sensitivity for those genes."
)

# --- Cell-wise Correlation (Restored Summary) ---
print("Calculating cell-wise correlations...")
def get_cell_corrs(adata1, adata2):
    m1 = adata1.X.toarray() if sp.issparse(adata1.X) else adata1.X
    m2 = adata2.X.toarray() if sp.issparse(adata2.X) else adata2.X
    corrs = []
    # Sample if too many cells
    n_sample = min(10000, m1.shape[0])
    idx = np.random.choice(m1.shape[0], n_sample, replace=False)
    for i in idx:
        r, _ = stats.pearsonr(m1[i], m2[i])
        corrs.append(r)
    return np.array(corrs)

cell_corrs = get_cell_corrs(cur_sub, rep_sub)
results_summary["cell_wise_corr"] = {
    "median": float(np.nanmedian(cell_corrs)),
    "mean": float(np.nanmean(cell_corrs))
}

# --- Perturbations ---
results_summary["perturbation"] = compare_perturbations(adata_cur, adata_rep, common_cells)

# --- Final Report ---
print("\n" + "=" * 50)
print("COMPARISON COMPLETE")
print("=" * 50)

with open("comparison_results/summary_report.txt", "w") as f:
    f.write(json.dumps(results_summary, indent=4))

print("\n[LLM_DATA_START]")
print(json.dumps(results_summary))
print("[LLM_DATA_END]\n")

display(pd.DataFrame([
    {"Metric": "Cell Overlap", "Value": f"{len(common_cells)} ({len(common_cells)/adata_cur.n_obs:.1%})"},
    {"Metric": "Gene Overlap", "Value": f"{len(common_genes)} ({len(common_genes)/adata_cur.n_vars:.1%})"},
    {"Metric": "Counts Correlation", "Value": f"{results_summary['cell_metrics']['total_counts']['pearson']:.4f}"},
    {"Metric": "Expression Correlation (Mean)", "Value": f"{results_summary['gene_metrics']['mean_expression']['pearson']:.4f}"},
    {"Metric": "Median Cell Correlation", "Value": f"{results_summary['cell_wise_corr']['median']:.4f}"},
    {"Metric": "Perturbation Accuracy", "Value": f"{results_summary['perturbation']['accuracy']:.2%}" if results_summary["perturbation"] else "N/A"},
]))
