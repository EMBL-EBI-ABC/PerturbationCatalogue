import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import sys
from concurrent.futures import ThreadPoolExecutor
import scipy.sparse as sp
import gc

# Ensure LAKE_BUCKET is available
lake_bucket = os.environ.get("LAKE_BUCKET")
if not lake_bucket:
    print("Error: LAKE_BUCKET environment variable is not set.")
    sys.exit(1)

# ==============================================================================
# 1. FUNCTIONS AND UTILITIES
# ==============================================================================


def download_only(gs_path, local_path):
    """Downloads an H5AD file from Google Cloud Storage if it doesn't exist."""
    if os.path.exists(local_path):
        print(f"File {local_path} already exists. Skipping download.")
    else:
        print(f"Downloading {gs_path} to {local_path}...")
        os.system(f"gsutil -m cp {gs_path} {local_path}")


def preprocess_adata(adata, is_raw, name, run_hvg=False):
    """Normalization, QC metrics, and PCA in one go for parallel execution."""
    print(f"Starting preprocessing for {name}...")
    if is_raw:
        print(f"  [{name}] Normalizing and log-transforming...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    print(f"  [{name}] Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    if run_hvg:
        print(f"  [{name}] Identifying highly variable genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")

    return adata


def is_raw_counts(adata):
    """Heuristic to check if data contains raw counts (integers)."""
    X = adata.X
    if hasattr(X, "data"):
        # Sparse matrix
        sample = X.data[:1000]
    else:
        # Dense matrix
        sample = X.flatten()[:1000]
    return np.all(np.equal(np.mod(sample, 1), 0))


def plot_scatter_comparison(df, x_col, y_col, title, xlabel, ylabel, filename):
    """Creates a scatter plot with correlation info."""
    plt.figure(figsize=(8, 8))

    # Use hexbin or scatter with low alpha if many points
    if len(df) > 10000:
        plt.hexbin(df[x_col], df[y_col], gridsize=50, cmap="Blues", bins="log")
        plt.colorbar(label="log10(count)")
    else:
        sns.scatterplot(data=df, x=x_col, y=y_col, alpha=0.5, s=10)

    # Add identity line
    min_val = min(df[x_col].min(), df[y_col].min())
    max_val = max(df[x_col].max(), df[y_col].max())
    plt.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="Identity")

    # Calculate correlation
    mask = ~(df[x_col].isna() | df[y_col].isna())
    corr, p_val = stats.pearsonr(df.loc[mask, x_col], df.loc[mask, y_col])
    spearman, _ = stats.spearmanr(df.loc[mask, x_col], df.loc[mask, y_col])

    plt.title(f"{title}\nPearson r = {corr:.4f}, Spearman rho = {spearman:.4f}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()
    plt.close()


def compare_perturbations(adata_cur, adata_rep, common_cells):
    """Compares perturbation assignments between curated and reprocessed datasets."""
    print("Comparing perturbation assignments...")

    # Try to find perturbation columns
    cur_pert_col = None
    for col in ["perturbation", "condition", "gene_target", "guide"]:
        if col in adata_cur.obs.columns:
            cur_pert_col = col
            break

    rep_pert_col = None
    for col in ["perturbation", "guide", "target", "gene"]:
        if col in adata_rep.obs.columns:
            rep_pert_col = col
            break

    if cur_pert_col and rep_pert_col:
        p_cur = adata_cur.obs.loc[common_cells, cur_pert_col].astype(str)
        p_rep = adata_rep.obs.loc[common_cells, rep_pert_col].astype(str)

        # Create a confusion matrix
        overlap_df = pd.DataFrame({"Curated": p_cur, "Reprocessed": p_rep})

        # Top 20 perturbations by frequency in curated
        top_perts = p_cur.value_counts().head(20).index
        mask = overlap_df["Curated"].isin(top_perts) & overlap_df["Reprocessed"].isin(
            top_perts
        )

        if mask.any():
            ct = pd.crosstab(
                overlap_df.loc[mask, "Curated"], overlap_df.loc[mask, "Reprocessed"]
            )
            plt.figure(figsize=(12, 10))
            sns.heatmap(ct, annot=False, cmap="YlGnBu")
            plt.title("Perturbation Assignment Overlap (Top 20)")
            plt.tight_layout()
            plt.savefig("comparison_results/perturbation_confusion_matrix.png")
            plt.show()
            plt.close()

            # Accuracy (simple matching)
            accuracy = (p_cur == p_rep).mean()
            print(f"Perturbation assignment match: {accuracy:.2%}")
            return accuracy
    else:
        print("Could not find matching perturbation columns in both datasets.")
        return None


# ==============================================================================
# 2. CONFIGURATION AND DOWNLOAD
# ==============================================================================
# File paths
curated_gs = f"gs://{lake_bucket}/perturbseq/curated/nadig_2025_jurkat_curated.h5ad"
curated_local = "nadig_2025_jurkat_curated.h5ad"
reprocessed_summed_local = "nadig_2025_jurkat_reprocessed_summed.h5ad"

os.makedirs("comparison_results", exist_ok=True)

# We expect reprocessed_summed_local to be generated by process_reprocessed.py
if not os.path.exists(reprocessed_summed_local):
    print(f"Error: {reprocessed_summed_local} not found.")
    print(f"Please run 'python process_reprocessed.py' first.")
    sys.exit(1)

print("Downloading curated file...")
download_only(curated_gs, curated_local)

# ==============================================================================
# 3. PARALLEL DATA LOADING
# ==============================================================================
print("Loading datasets into memory in parallel...")
with ThreadPoolExecutor(max_workers=2) as executor:
    f_cur = executor.submit(sc.read_h5ad, curated_local)
    f_rep = executor.submit(sc.read_h5ad, reprocessed_summed_local)
    adata_cur = f_cur.result()
    adata_rep_sum = f_rep.result()

# ==============================================================================
# 4. HEURISTIC CHECK FOR NORMALIZATION
# ==============================================================================
raw_cur = is_raw_counts(adata_cur)
raw_rep = is_raw_counts(adata_rep_sum)

print(f"Curated is raw: {raw_cur}")
print(f"Reprocessed is raw: {raw_rep}")

# ==============================================================================
# 6. RAW DISTRIBUTION HISTOGRAMS (CELLS AND GENES)
# ==============================================================================
print("Computing raw distribution statistics in parallel...")


def get_sums(adata):
    """Computes total counts per cell and per gene."""
    cell_sums = np.array(adata.X.sum(axis=1)).flatten()
    gene_sums = np.array(adata.X.sum(axis=0)).flatten()
    return cell_sums, gene_sums


with ThreadPoolExecutor(max_workers=2) as executor:
    f_cur_sums = executor.submit(get_sums, adata_cur)
    f_rep_sums = executor.submit(get_sums, adata_rep_sum)

    cur_cell_sums, cur_gene_sums = f_cur_sums.result()
    rep_cell_sums, rep_gene_sums = f_rep_sums.result()

# Plotting histograms
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Cell counts distribution
sns.histplot(
    cur_cell_sums,
    bins=100,
    label="Curated",
    ax=ax1,
    alpha=0.5,
    log_scale=True,
    color="blue",
)
sns.histplot(
    rep_cell_sums,
    bins=100,
    label="Reprocessed (Summed)",
    ax=ax1,
    alpha=0.5,
    log_scale=True,
    color="orange",
)
ax1.set_title("Total Counts per Cell (Log Scale)")
ax1.set_xlabel("Total UMI Counts")
ax1.legend()

# Gene counts distribution
sns.histplot(
    cur_gene_sums,
    bins=100,
    label="Curated",
    ax=ax2,
    alpha=0.5,
    log_scale=True,
    color="blue",
)
sns.histplot(
    rep_gene_sums,
    bins=100,
    label="Reprocessed (Summed)",
    ax=ax2,
    alpha=0.5,
    log_scale=True,
    color="orange",
)
ax2.set_title("Total Counts per Gene (Log Scale)")
ax2.set_xlabel("Total UMI Counts")
ax2.legend()

plt.tight_layout()
plt.savefig("comparison_results/raw_distributions_histogram.png", dpi=300)
plt.show()

# ==============================================================================
# 7. BASIC SUMMARY STATISTICS
# ==============================================================================
summary_df = pd.DataFrame(
    {
        "Metric": [
            "Total Cells (n_obs)",
            "Total Genes (n_vars)",
            "Obs Columns",
            "Var Columns",
            "Layers",
            "Unstructured (uns) Keys",
        ],
        "Curated": [
            adata_cur.n_obs,
            adata_cur.n_vars,
            ", ".join(adata_cur.obs.columns),
            ", ".join(adata_cur.var.columns),
            ", ".join(adata_cur.layers.keys()),
            ", ".join(adata_cur.uns.keys()),
        ],
        "Reprocessed (summed)": [
            adata_rep_sum.n_obs,
            adata_rep_sum.n_vars,
            ", ".join(adata_rep_sum.obs.columns),
            ", ".join(adata_rep_sum.var.columns),
            ", ".join(adata_rep_sum.layers.keys()),
            ", ".join(adata_rep_sum.uns.keys()),
        ],
    }
)
print("\n### Basic Dataset Comparison ###")
with pd.option_context("display.max_colwidth", None, "display.max_rows", None):
    print(summary_df)

# ==============================================================================
# 8. DATAFRAME EXPLORATION
# ==============================================================================
print("\n### Curated - Obs (first 5 rows) ###")
with pd.option_context("display.max_colwidth", None, "display.max_rows", None):
    print(adata_cur.obs.head())

print("\n### Curated - Var (first 5 rows) ###")
with pd.option_context("display.max_colwidth", None, "display.max_rows", None):
    print(adata_cur.var.head())

print("\n" + "=" * 40)

print("\n### Reprocessed (summed) - Obs (first 5 rows) ###")
with pd.option_context("display.max_colwidth", None, "display.max_rows", None):
    print(adata_rep_sum.obs.head())

print("\n### Reprocessed (summed) - Var (first 5 rows) ###")
with pd.option_context("display.max_colwidth", None, "display.max_rows", None):
    print(adata_rep_sum.var.head())

# ==============================================================================
# 9. GENE AND CELL ALIGNMENT
# ==============================================================================
# Gene Alignment
common_genes = np.intersect1d(adata_cur.var_names, adata_rep_sum.var_names)
if len(common_genes) == 0:
    print("No direct gene overlap. Attempting to align via var columns...")
    for col in ["gene_symbols", "symbols", "gene_name"]:
        if col in adata_cur.var.columns:
            adata_cur.var_names = adata_cur.var[col].astype(str)
            break
    for col in ["gene_symbols", "symbols", "gene_name"]:
        if col in adata_rep_sum.var.columns:
            adata_rep_sum.var_names = adata_rep_sum.var[col].astype(str)
            break
    common_genes = np.intersect1d(adata_cur.var_names, adata_rep_sum.var_names)

# Cell Alignment
common_cells = np.intersect1d(adata_cur.obs_names, adata_rep_sum.obs_names)

print(f"Common cells: {len(common_cells)}")
print(f"Common genes: {len(common_genes)}")

if len(common_cells) == 0 or len(common_genes) == 0:
    print("ERROR: No overlap found.")
else:
    # Subset to common elements
    cur_sub = adata_cur[common_cells, common_genes].copy()
    rep_sub = adata_rep_sum[common_cells, common_genes].copy()

# ==============================================================================
# 10. PREPROCESSING
# ==============================================================================
print("Preprocessing subsets in parallel...")
with ThreadPoolExecutor(max_workers=2) as executor:
    f_cur = executor.submit(preprocess_adata, cur_sub, raw_cur, "Curated", run_hvg=True)
    f_rep = executor.submit(
        preprocess_adata, rep_sub, raw_rep, "Reprocessed", run_hvg=False
    )
    cur_sub = f_cur.result()
    rep_sub = f_rep.result()

# ==============================================================================
# 11. QC AND GENE EXPRESSION COMPARISON
# ==============================================================================
metrics_df = pd.DataFrame(
    {
        "total_counts_cur": cur_sub.obs["total_counts"],
        "total_counts_rep": rep_sub.obs["total_counts"],
        "n_genes_cur": cur_sub.obs["n_genes_by_counts"],
        "n_genes_rep": rep_sub.obs["n_genes_by_counts"],
    }
)

plot_scatter_comparison(
    metrics_df,
    "total_counts_cur",
    "total_counts_rep",
    "Total Counts (Normalized/Log) Correlation",
    "Curated",
    "Reprocessed",
    "comparison_results/counts_comparison.png",
)

plot_scatter_comparison(
    metrics_df,
    "n_genes_cur",
    "n_genes_rep",
    "Detected Genes Correlation",
    "Curated",
    "Reprocessed",
    "comparison_results/genes_comparison.png",
)

# Gene expression correlation
gene_metrics = pd.DataFrame(
    {"mean_cur": cur_sub.var["mean_counts"], "mean_rep": rep_sub.var["mean_counts"]}
)

plot_scatter_comparison(
    gene_metrics,
    "mean_cur",
    "mean_rep",
    "Mean Gene Expression Correlation",
    "Curated",
    "Reprocessed",
    "comparison_results/gene_expression_mean.png",
)

# ==============================================================================
# 12. PERTURBATION COMPARISON
# ==============================================================================
pert_acc = compare_perturbations(adata_cur, adata_rep_sum, common_cells)

# ==============================================================================
# 13. STRUCTURAL COMPARISON (PCA)
# ==============================================================================
print("Performing structural comparison (PCA)...")
# Use the same highly variable genes for both to ensure comparability
rep_sub.var["highly_variable"] = cur_sub.var["highly_variable"]

print("  Calculating PCA in parallel...")
with ThreadPoolExecutor(max_workers=2) as executor:
    f_cur = executor.submit(
        sc.tl.pca,
        cur_sub,
        n_comps=30,
        use_highly_variable=True,
        svd_solver="arpack",
    )
    f_rep = executor.submit(
        sc.tl.pca,
        rep_sub,
        n_comps=30,
        use_highly_variable=True,
        svd_solver="arpack",
    )
    f_cur.result()
    f_rep.result()

# Plot side-by-side PCA
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.pca(cur_sub, ax=ax1, show=False, title="PCA Curated (on Curated HVGs)")
sc.pl.pca(rep_sub, ax=ax2, show=False, title="PCA Reprocessed (on Curated HVGs)")
plt.tight_layout()
plt.savefig("comparison_results/pca_comparison.png")
plt.show()
plt.close()

# ==============================================================================
# 14. SUMMARY REPORT
# ==============================================================================
with open("comparison_results/summary_report.txt", "w") as f:
    f.write("Nadig 2025 Jurkat Comparison Report\n")
    f.write("===================================\n\n")
    f.write(f"Curated dataset: {adata_cur.n_obs} cells, {adata_cur.n_vars} genes\n")
    f.write(
        f"Reprocessed dataset: {adata_rep_sum.n_obs} cells, {adata_rep_sum.n_vars} genes\n"
    )
    f.write(
        f"Cell overlap: {len(common_cells)} ({len(common_cells)/adata_cur.n_obs:.1%} of curated)\n"
    )
    f.write(
        f"Gene overlap: {len(common_genes)} ({len(common_genes)/adata_cur.n_vars:.1%} of curated)\n"
    )

    counts_corr, _ = stats.pearsonr(
        metrics_df["total_counts_cur"], metrics_df["total_counts_rep"]
    )
    genes_corr, _ = stats.pearsonr(metrics_df["n_genes_cur"], metrics_df["n_genes_rep"])
    expr_corr, _ = stats.pearsonr(gene_metrics["mean_cur"], gene_metrics["mean_rep"])

    f.write(f"\nPearson Correlations for common elements:\n")
    f.write(f"- Total counts per cell: {counts_corr:.4f}\n")
    f.write(f"- Number of genes per cell: {genes_corr:.4f}\n")
    f.write(f"- Average gene expression: {expr_corr:.4f}\n")

    if pert_acc is not None:
        f.write(f"- Perturbation assignment match: {pert_acc:.2%}\n")

print("\nComparison finished. Check the 'comparison_results' folder for details.")
