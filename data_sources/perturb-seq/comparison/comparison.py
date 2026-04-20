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
import json

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
    # Handle possible duplicates after cleaning
    if adata.obs_names.duplicated().any():
        print(
            f"  Warning: Found {adata.obs_names.duplicated().sum()} duplicate barcodes after cleaning. Suffixing..."
        )
        adata.obs_names_make_unique()
    return adata


def get_raw_counts(adata):
    """Safely extracts raw counts or identifies if data is already normalized."""
    X = adata.X
    if hasattr(X, "data"):
        sample = X.data[:2000]
    else:
        sample = X.flatten()[:2000]

    is_integers = np.all(np.equal(np.mod(sample, 1), 0))
    return is_integers


def preprocess_adata(adata, name, target_sum=1e4, n_top_genes=2000):
    """Standardized preprocessing: Raw -> Norm -> Log -> HVG -> Scale -> PCA."""
    print(f"Starting standardized preprocessing for {name}...")

    # Check if raw
    is_raw = get_raw_counts(adata)
    if not is_raw:
        print(
            f"  [{name}] WARNING: Data does not appear to be raw counts. Preprocessing might be biased."
        )

    # Store raw in layer if possible
    adata.layers["counts"] = adata.X.copy()

    # Calculate QC on raw
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Normalize and Log
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # HVG
    print(f"  [{name}] Identifying highly variable genes...")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat", subset=False
    )

    return adata


def plot_scatter_comparison(
    df, x_col, y_col, title, xlabel, ylabel, filename, log_scale=False
):
    """Creates a scatter plot with correlation info and identity line."""
    plt.figure(figsize=(8, 8))

    plot_df = df.copy().dropna(subset=[x_col, y_col])

    if len(plot_df) > 10000:
        plt.hexbin(
            plot_df[x_col], plot_df[y_col], gridsize=50, cmap="viridis", bins="log"
        )
        plt.colorbar(label="log10(count)")
    else:
        sns.scatterplot(data=plot_df, x=x_col, y=y_col, alpha=0.3, s=10)

    # Add identity line
    min_val = min(plot_df[x_col].min(), plot_df[y_col].min())
    max_val = max(plot_df[x_col].max(), plot_df[y_col].max())
    plt.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="Identity")

    # Calculate correlations
    pearson, _ = stats.pearsonr(plot_df[x_col], plot_df[y_col])
    spearman, _ = stats.spearmanr(plot_df[x_col], plot_df[y_col])

    plt.title(f"{title}\nPearson r = {pearson:.4f}, Spearman rho = {spearman:.4f}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()
    plt.close()

    return {"pearson": pearson, "spearman": spearman}


def call_guides(adata, count_threshold=5):
    """Dual-guide caller: identifies top 2 guides and joins them with '|'."""
    if "guides" not in adata.obsm:
        return None

    guide_matrix = adata.obsm["guides"]
    guide_names = np.array(
        adata.uns.get("guide_names", [f"guide_{i}" for i in range(guide_matrix.shape[1])])
    )

    # Convert to CSR if not already for efficient row slicing
    if not sp.iscsr(guide_matrix):
        guide_matrix = sp.csr_matrix(guide_matrix)

    print(f"  Calling dual guides (threshold={count_threshold})...")
    calls = []
    for i in range(guide_matrix.shape[0]):
        row = guide_matrix.getrow(i).toarray().flatten()
        # Find indices of guides above threshold
        top_idx = np.where(row >= count_threshold)[0]
        
        if len(top_idx) == 0:
            calls.append("None")
        else:
            # Sort by counts descending and take top 2
            top_idx = top_idx[np.argsort(row[top_idx])[::-1][:2]]
            # Sort names alphabetically for consistent pipe joining
            names = sorted(guide_names[top_idx])
            calls.append("|".join(names))

    return pd.Series(calls, index=adata.obs_names)


def compare_perturbations(adata_cur, adata_rep, common_cells):
    """Compares perturbation assignments, handling both .obs and .obsm['guides']."""
    print("Comparing perturbation assignments...")

    # 1. Get Curated Labels
    cur_pert_col = None
    for col in ["perturbation", "condition", "gene_target", "guide", "target"]:
        if col in adata_cur.obs.columns:
            cur_pert_col = col
            break

    # 2. Get Reprocessed Labels
    rep_labels = None
    if "perturbation" in adata_rep.obs.columns:
        rep_labels = adata_rep.obs["perturbation"]
    elif "guides" in adata_rep.obsm:
        print("  Calling guides from .obsm['guides'] for reprocessed data...")
        rep_labels = call_guides(adata_rep)

    if cur_pert_col and rep_labels is not None:
        p_cur = adata_cur.obs.loc[common_cells, cur_pert_col].astype(str)
        p_rep = rep_labels.loc[common_cells].astype(str)

        overlap_df = pd.DataFrame({"Curated": p_cur, "Reprocessed": p_rep})

        # Calculate overall accuracy
        match_mask = p_cur == p_rep
        accuracy = match_mask.mean()

        # Top 20 perturbations in Curated
        top_perts = p_cur.value_counts().head(20).index
        sub_df = overlap_df[overlap_df["Curated"].isin(top_perts)]

        ct = pd.crosstab(sub_df["Curated"], sub_df["Reprocessed"])

        plt.figure(figsize=(12, 10))
        sns.heatmap(ct, annot=False, cmap="YlGnBu")
        plt.title(
            f"Perturbation Confusion Matrix (Top 20)\nOverall Match: {accuracy:.2%}"
        )
        plt.tight_layout()
        plt.savefig("comparison_results/perturbation_confusion_matrix.png")
        plt.show()
        plt.close()

        print(f"  Perturbation assignment match: {accuracy:.2%}")
        return {
            "accuracy": accuracy,
            "cur_col": cur_pert_col,
            "n_matched": match_mask.sum(),
            "n_total": len(common_cells),
        }
    else:
        print("  Could not find comparable perturbation data.")
        return None


# ==============================================================================
# 2. CONFIGURATION AND PATHS
# ==============================================================================
curated_local = "GSE264667_jurkat_raw_singlecell_01.h5ad"
reprocessed_local = "experiment_final.h5ad"
os.makedirs("comparison_results", exist_ok=True)

# ==============================================================================
# 3. DATA LOADING AND ALIGNMENT
# ==============================================================================
print("Loading datasets...")
adata_cur = sc.read_h5ad(curated_local)
adata_rep = sc.read_h5ad(reprocessed_local)

print("Aligning barcodes and genes...")
adata_cur = clean_barcodes(adata_cur)
adata_rep = clean_barcodes(adata_rep)

# Strip Ensembl versions from reprocessed data if present
if adata_rep.var_names.str.contains(r"\.").any():
    print("  Stripping Ensembl versions from reprocessed var_names...")
    adata_rep.var_names = adata_rep.var_names.str.split(".").str[0]
    adata_rep.var_names_make_unique()

# Gene Alignment
common_genes = np.intersect1d(adata_cur.var_names, adata_rep.var_names)
if len(common_genes) < 100:
    print("  Low gene overlap. Attempting alignment via curated 'gene_name' vs reprocessed index...")
    # Curated has gene_name, Reprocessed has index (Ensembl)
    # This might happen if Curated var_names are Symbols. 
    # Let's check if Curated var names look like Ensembl
    is_ensembl_cur = adata_cur.var_names.str.startswith("ENSG").any()
    
    if not is_ensembl_cur and "gene_name" in adata_cur.var.columns:
        print("    Curated var_names are likely symbols. Switching to Ensembl IDs from index...")
        # (Actually, diagnostics showed Curated var_names ARE Ensembl IDs)
        pass

    # Fallback search
    cur_cols = [c for c in ["gene_symbols", "symbols", "gene_name", "gene_ids"] if c in adata_cur.var.columns] + ["index"]
    rep_cols = [c for c in ["gene_symbols", "symbols", "gene_name", "gene_ids"] if c in adata_rep.var.columns] + ["index"]
    
    best_overlap = len(common_genes)
    best_pair = ("index", "index")

    for c1 in cur_cols:
        v1 = adata_cur.var_names if c1 == "index" else adata_cur.var[c1].astype(str)
        for c2 in rep_cols:
            v2 = adata_rep.var_names if c2 == "index" else adata_rep.var[c2].astype(str)
            overlap = np.intersect1d(v1, v2)
            if len(overlap) > best_overlap:
                best_overlap = len(overlap)
                best_pair = (c1, c2)
    
    if best_overlap > len(common_genes):
        print(f"    Best alignment: Curated['{best_pair[0]}'] vs Reprocessed['{best_pair[1]}'] ({best_overlap} genes)")
        if best_pair[0] != "index":
            adata_cur.var_names = adata_cur.var[best_pair[0]].astype(str)
        if best_pair[1] != "index":
            adata_rep.var_names = adata_rep.var[best_pair[1]].astype(str)
        
        adata_cur.var_names_make_unique()
        adata_rep.var_names_make_unique()
        common_genes = np.intersect1d(adata_cur.var_names, adata_rep.var_names)

common_cells = np.intersect1d(adata_cur.obs_names, adata_rep.obs_names)

print(f"Alignment Summary:")
print(f"  Common cells: {len(common_cells)} ({len(common_cells)/adata_cur.n_obs:.1%})")
print(f"  Common genes: {len(common_genes)} ({len(common_genes)/adata_cur.n_vars:.1%})")

if len(common_cells) == 0 or len(common_genes) == 0:
    print("CRITICAL ERROR: Zero overlap. Stopping.")
    sys.exit(1)

# Subset to common elements
cur_sub = adata_cur[common_cells, common_genes].copy()
rep_sub = adata_rep[common_cells, common_genes].copy()
gc.collect()

# ==============================================================================
# 4. PREPROCESSING
# ==============================================================================
cur_sub = preprocess_adata(cur_sub, "Curated")
rep_sub = preprocess_adata(rep_sub, "Reprocessed")

# ==============================================================================
# 5. METRICS CALCULATION
# ==============================================================================
results_summary = {
    "dimensions": {
        "cur": [adata_cur.n_obs, adata_cur.n_vars],
        "rep": [adata_rep.n_obs, adata_rep.n_vars],
        "overlap": [len(common_cells), len(common_genes)],
    }
}

# 1. Total Counts & Gene Count Correlation
metrics_df = pd.DataFrame(
    {
        "counts_cur": cur_sub.obs["total_counts"],
        "counts_rep": rep_sub.obs["total_counts"],
        "genes_cur": cur_sub.obs["n_genes_by_counts"],
        "genes_rep": rep_sub.obs["n_genes_by_counts"],
    }
)

results_summary["cell_metrics"] = {
    "total_counts": plot_scatter_comparison(
        metrics_df,
        "counts_cur",
        "counts_rep",
        "Total UMI Counts",
        "Curated",
        "Reprocessed",
        "comparison_results/counts_comparison.png",
        log_scale=True,
    ),
    "n_genes": plot_scatter_comparison(
        metrics_df,
        "genes_cur",
        "genes_rep",
        "Number of Detected Genes",
        "Curated",
        "Reprocessed",
        "comparison_results/genes_comparison.png",
    ),
}

# 2. Gene-wise Statistics (Mean & Sparsity)
gene_metrics = pd.DataFrame(
    {
        "mean_cur": cur_sub.var["mean_counts"],
        "mean_rep": rep_sub.var["mean_counts"],
        "pct_dropout_cur": 100
        - (
            np.array((cur_sub.layers["counts"] > 0).sum(axis=0)).flatten()
            / cur_sub.n_obs
            * 100
        ),
        "pct_dropout_rep": 100
        - (
            np.array((rep_sub.layers["counts"] > 0).sum(axis=0)).flatten()
            / rep_sub.n_obs
            * 100
        ),
    }
)

results_summary["gene_metrics"] = {
    "mean_expression": plot_scatter_comparison(
        gene_metrics,
        "mean_cur",
        "mean_rep",
        "Mean Gene Expression",
        "Curated",
        "Reprocessed",
        "comparison_results/gene_expression_mean.png",
        log_scale=True,
    ),
    "sparsity": plot_scatter_comparison(
        gene_metrics,
        "pct_dropout_cur",
        "pct_dropout_rep",
        "Gene Dropout Rate (%)",
        "Curated",
        "Reprocessed",
        "comparison_results/sparsity_comparison.png",
    ),
}

# 3. HVG Overlap (Jaccard)
hvg_cur = set(cur_sub.var_names[cur_sub.var["highly_variable"]])
hvg_rep = set(rep_sub.var_names[rep_sub.var["highly_variable"]])
jaccard = len(hvg_cur & hvg_rep) / len(hvg_cur | hvg_rep) if hvg_cur or hvg_rep else 0
results_summary["hvg"] = {
    "jaccard": jaccard,
    "n_cur": len(hvg_cur),
    "n_rep": len(hvg_rep),
    "n_shared": len(hvg_cur & hvg_rep),
}
print(f"HVG Jaccard Similarity: {jaccard:.4f}")

# 4. Cell-wise Correlation
print("Calculating cell-wise correlations...")


def get_cell_corrs(adata1, adata2):
    m1 = adata1.X
    m2 = adata2.X
    if sp.issparse(m1):
        m1 = m1.toarray()
    if sp.issparse(m2):
        m2 = m2.toarray()

    # Pearson corr per row
    corrs = []
    for i in range(m1.shape[0]):
        r, _ = stats.pearsonr(m1[i], m2[i])
        corrs.append(r)
    return np.array(corrs)


# Sample if too large
if cur_sub.n_obs > 10000:
    idx = np.random.choice(cur_sub.n_obs, 10000, replace=False)
    cell_corrs = get_cell_corrs(cur_sub[idx], rep_sub[idx])
else:
    cell_corrs = get_cell_corrs(cur_sub, rep_sub)

results_summary["cell_wise_corr"] = {
    "median": float(np.nanmedian(cell_corrs)),
    "mean": float(np.nanmean(cell_corrs)),
}

plt.figure(figsize=(8, 5))
sns.histplot(cell_corrs, bins=50, color="teal")
plt.axvline(
    np.nanmedian(cell_corrs),
    color="red",
    linestyle="--",
    label=f"Median: {np.nanmedian(cell_corrs):.4f}",
)
plt.title("Distribution of Cell-wise Expression Correlations")
plt.xlabel("Pearson r")
plt.tight_layout()
plt.savefig("comparison_results/cell_wise_correlation_hist.png")
plt.show()

# 5. Perturbations
results_summary["perturbation"] = compare_perturbations(
    adata_cur, adata_rep, common_cells
)

# 6. Structural Comparison (PCA)
print("Performing PCA comparison...")
sc.pp.scale(cur_sub, max_value=10)
sc.pp.scale(rep_sub, max_value=10)

# Calculate PCA on common HVGs for direct comparison
common_hvgs = list(hvg_cur & hvg_rep)
if len(common_hvgs) < 50:
    common_hvgs = list(hvg_cur | hvg_rep)[:500]

sc.tl.pca(cur_sub, n_comps=30, use_highly_variable=False)
cur_sub.uns["pca_on_shared"] = common_hvgs  # tag it
sc.tl.pca(rep_sub, n_comps=30, use_highly_variable=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.pca(cur_sub, ax=ax1, show=False, title="PCA Curated (Common HVGs)")
sc.pl.pca(rep_sub, ax=ax2, show=False, title="PCA Reprocessed (Common HVGs)")
plt.tight_layout()
plt.savefig("comparison_results/pca_comparison.png")
plt.show()

# ==============================================================================
# 6. FINAL REPORT
# ==============================================================================
print("\n" + "=" * 50)
print("COMPARISON COMPLETE")
print("=" * 50)

# Write human-readable summary
with open("comparison_results/summary_report.txt", "w") as f:
    f.write(json.dumps(results_summary, indent=4))

# LLM-Readable Output
print("\n[LLM_DATA_START]")
print(json.dumps(results_summary))
print("[LLM_DATA_END]\n")

display(
    pd.DataFrame(
        [
            {
                "Metric": "Cell Overlap",
                "Value": f"{len(common_cells)} ({len(common_cells)/adata_cur.n_obs:.1%})",
            },
            {
                "Metric": "Gene Overlap",
                "Value": f"{len(common_genes)} ({len(common_genes)/adata_cur.n_vars:.1%})",
            },
            {
                "Metric": "Counts Correlation",
                "Value": f"{results_summary['cell_metrics']['total_counts']['pearson']:.4f}",
            },
            {
                "Metric": "Expression Correlation (Mean)",
                "Value": f"{results_summary['gene_metrics']['mean_expression']['pearson']:.4f}",
            },
            {
                "Metric": "Median Cell Correlation",
                "Value": f"{results_summary['cell_wise_corr']['median']:.4f}",
            },
            {
                "Metric": "HVG Jaccard",
                "Value": f"{results_summary['hvg']['jaccard']:.4f}",
            },
            {
                "Metric": "Perturbation Accuracy",
                "Value": (
                    f"{results_summary['perturbation']['accuracy']:.2%}"
                    if results_summary["perturbation"]
                    else "N/A"
                ),
            },
        ]
    )
)
