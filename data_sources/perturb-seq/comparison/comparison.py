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


CURATED_H5AD_PATH = (
    "/hps/nobackup/mfreeberg/perturb_seq_fastq/source_h5ad/nadig_2025_jurkat.h5ad"
)
REPROCESSED_H5AD_PATH = (
    "/hps/nobackup/mfreeberg/perturb_seq_fastq/results/"
    "nadig_2025_jurkat/experiment_final.h5ad"
)
GUIDE_COUNT_THRESHOLDS = [1, 2, 3, 4, 5]
GENE_CALL_OUTCOMES = ["0_genes", "1_gene_1_probe", "1_gene_2_probes", ">1_gene"]


# ==============================================================================
# 1. FUNCTIONS AND UTILITIES
# ==============================================================================


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, set):
        return sorted(value)
    return str(value)


def log_record(event, **fields):
    print(json.dumps({"event": event, **fields}, sort_keys=True, default=json_default))


def filter_unique_barcodes(adata, name):
    """Strips suffixes and keeps only globally unique barcodes."""
    original_count = adata.n_obs

    # Strip the suffix (e.g., "-1" or "-8")
    base_barcodes = adata.obs_names.str.split("-").str[0]

    # Identify barcodes that are truly unique across the dataset
    barcode_counts = base_barcodes.value_counts()
    unique_barcodes = set(barcode_counts[barcode_counts == 1].index)

    # Filter the matrix
    mask = base_barcodes.isin(unique_barcodes)
    adata = adata[mask].copy()

    # Re-index with the base barcode for direct comparison
    adata.obs_names = base_barcodes[mask]

    filtered_count = original_count - adata.n_obs
    log_record(
        "barcode_filter",
        dataset=name,
        original_cells=original_count,
        removed_collision_cells=filtered_count,
        removed_collision_pct=(
            (filtered_count / original_count) * 100 if original_count else 0.0
        ),
        remaining_cells=adata.n_obs,
    )
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
    is_raw = get_raw_counts(adata)
    log_record(
        "preprocess_input",
        dataset=name,
        cells=adata.n_obs,
        genes=adata.n_vars,
        raw_counts_detected=bool(is_raw),
        normalization_target_sum=target_sum,
        hvg_n_top_genes=n_top_genes,
    )

    adata.layers["counts"] = adata.X.copy()
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat", subset=False
    )
    return adata


def plot_scatter_comparison(
    df,
    x_col,
    y_col,
    title,
    xlabel,
    ylabel,
    filename,
    log_scale=False,
    description="",
    deviation_on_log=False,
):
    """
    Standardized scatter plot with professional layout:
    - Bold Title
    - Subtitle with stats
    - Wrapped Description at the bottom
    """
    fig, ax = plt.subplots(figsize=(10, 11))

    plot_df = df.copy().dropna(subset=[x_col, y_col])

    # Calculate deviation stats (within 10%). For UMI counts, use the log-value so
    # proportional differences at very high depth do not dominate this summary.
    deviation_df = plot_df[[x_col, y_col]].astype(float)
    if deviation_on_log:
        deviation_df = np.log1p(deviation_df)

    diff = np.abs(deviation_df[x_col] - deviation_df[y_col])
    # Avoid division by zero for cells/genes with zero values in both datasets.
    denom = deviation_df[[x_col, y_col]].mean(axis=1).replace(0, 1)
    rel_diff = diff / denom
    pct_deviant = (rel_diff > 0.1).mean() * 100
    deviant_label = (
        "Deviants (>10% on log1p values)" if deviation_on_log else "Deviants (>10%)"
    )

    if log_scale:
        plot_df[x_col] = plot_df[x_col] + 1
        plot_df[y_col] = plot_df[y_col] + 1

    # Main plot
    if len(plot_df) > 5000:
        ax.scatter(
            plot_df[x_col],
            plot_df[y_col],
            alpha=0.05,
            s=1,
            color="teal",
            rasterized=True,
        )
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
    ax.set_title(title, fontsize=18, fontweight="bold", pad=35)
    ax.text(
        0.5,
        1.02,
        f"Pearson r = {pearson:.4f} | Spearman rho = {spearman:.4f} | {deviant_label}: {pct_deviant:.1f}%",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=11,
        style="italic",
    )

    ax.set_xlabel(f"{xlabel} {'(+1 for log)' if log_scale else ''}", fontsize=12)
    ax.set_ylabel(f"{ylabel} {'(+1 for log)' if log_scale else ''}", fontsize=12)

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.legend(loc="upper left")

    # Description at the bottom
    wrapped_desc = "\n".join(textwrap.wrap(description, width=100))
    fig.text(
        0.5, 0.02, wrapped_desc, ha="center", va="bottom", fontsize=10, linespacing=1.4
    )

    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    plt.savefig(filename, dpi=300)
    plt.show()
    plt.close()

    return {"pearson": pearson, "spearman": spearman, "pct_deviant": pct_deviant}


def canonical_guide_id(guide_name):
    """Normalize guide IDs to match the curated Nadig label formatting."""
    return str(guide_name).strip().replace(",", "-")


def guide_aliases(guide_name):
    """Return all guide IDs represented by a potentially semicolon-joined KITE feature."""
    aliases = [
        canonical_guide_id(part)
        for part in str(guide_name).split(";")
        if str(part).strip()
    ]
    return aliases or [canonical_guide_id(guide_name)]


def guide_target_name(guide_name):
    """Extract the target gene symbol from a probe ID."""
    guide_name = canonical_guide_id(guide_name)
    if guide_name.startswith("non-targeting"):
        return "non-targeting"
    return guide_name.split("_", 1)[0]


def guide_target_names(guide_name):
    return {guide_target_name(alias) for alias in guide_aliases(guide_name)}


def is_empty_guide_label(label):
    label = str(label).strip()
    return label.lower() in {"", "none", "nan"}


def canonical_label_for_display(label):
    label = str(label).strip()
    if is_empty_guide_label(label):
        return "None"
    parts = [canonical_guide_id(part) for part in label.split("|") if part.strip()]
    return "|".join(sorted(parts)) if parts else "None"


def label_to_probe_list(label):
    label = canonical_label_for_display(label)
    if label == "None":
        return []
    return [part for part in label.split("|") if part]


def genes_from_probe_list(probes):
    genes = set()
    for probe in probes:
        for gene in guide_target_names(probe):
            if gene and gene != "non-targeting":
                genes.add(gene)
    return sorted(genes)


def gene_call_outcome(probes, genes):
    if len(genes) == 0:
        return "0_genes"
    if len(genes) > 1:
        return ">1_gene"

    gene = genes[0]
    n_gene_probes = sum(gene in guide_target_names(probe) for probe in probes)
    if n_gene_probes <= 1:
        return "1_gene_1_probe"
    return "1_gene_2_probes"


def build_gene_call_table(labels):
    records = []
    for cell_id, label in labels.items():
        probe_label = canonical_label_for_display(label)
        probes = label_to_probe_list(probe_label)
        genes = genes_from_probe_list(probes)
        records.append(
            {
                "cell_id": cell_id,
                "probe_label": probe_label,
                "n_probes": len(probes),
                "genes": genes,
                "gene_label": "|".join(genes) if genes else "None",
                "n_genes": len(genes),
                "outcome": gene_call_outcome(probes, genes),
            }
        )
    return pd.DataFrame.from_records(records, index=labels.index)


def outcome_count_dict(outcomes):
    counts = outcomes.value_counts().reindex(GENE_CALL_OUTCOMES, fill_value=0)
    return {outcome: int(counts.loc[outcome]) for outcome in GENE_CALL_OUTCOMES}


def outcome_matrix_dict(cur_outcomes, rep_outcomes):
    matrix = pd.crosstab(cur_outcomes, rep_outcomes)
    matrix = matrix.reindex(
        index=GENE_CALL_OUTCOMES, columns=GENE_CALL_OUTCOMES, fill_value=0
    )
    return {
        cur_outcome: {
            rep_outcome: int(matrix.loc[cur_outcome, rep_outcome])
            for rep_outcome in GENE_CALL_OUTCOMES
        }
        for cur_outcome in GENE_CALL_OUTCOMES
    }


def single_gene_match_summary(cur_calls, rep_calls):
    cur_single = cur_calls["n_genes"] == 1
    rep_single = rep_calls["n_genes"] == 1
    both_single = cur_single & rep_single
    same_gene = (
        cur_calls.loc[both_single, "gene_label"]
        == rep_calls.loc[both_single, "gene_label"]
    )

    n_both_single = int(both_single.sum())
    n_same = int(same_gene.sum())
    n_different = n_both_single - n_same
    return {
        "n_cells_curated_single_gene": int(cur_single.sum()),
        "n_cells_reprocessed_single_gene": int(rep_single.sum()),
        "n_cells_both_single_gene": n_both_single,
        "n_same_gene": n_same,
        "n_different_gene": n_different,
        "pct_same_gene": (n_same / n_both_single * 100) if n_both_single else 0.0,
        "pct_different_gene": (
            n_different / n_both_single * 100 if n_both_single else 0.0
        ),
    }


def summarize_numeric(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"min": 0.0, "p50": 0.0, "p90": 0.0, "p99": 0.0, "max": 0.0}
    return {
        "min": float(np.min(values)),
        "p50": float(np.percentile(values, 50)),
        "p90": float(np.percentile(values, 90)),
        "p99": float(np.percentile(values, 99)),
        "max": float(np.max(values)),
    }


def fit_poisson_gaussian_mixture(counts, max_iter=100, tol=1e-4):
    """Fits a 2-component mixture (Poisson bg, Gaussian sig) on non-zero UMI counts."""
    vals, freqs = np.unique(counts, return_counts=True)

    threshold_init = np.percentile(counts, 90)
    sig_mask = vals >= threshold_init
    bg_mask = vals < threshold_init

    if sig_mask.sum() == 0 or bg_mask.sum() == 0:
        sig_mask = vals >= np.median(vals)
        bg_mask = ~sig_mask

    pi_bg = np.sum(freqs[bg_mask]) / np.sum(freqs)
    lambda_ = np.average(vals[bg_mask], weights=freqs[bg_mask])
    mu = np.average(vals[sig_mask], weights=freqs[sig_mask])
    var_sig = np.average((vals[sig_mask] - mu) ** 2, weights=freqs[sig_mask])
    sigma = np.sqrt(var_sig) if var_sig > 0 else 1.0

    lambda_ = max(lambda_, 0.1)
    sigma = max(sigma, 1.0)

    log_likelihood = -np.inf
    W_sig = np.zeros_like(vals, dtype=float)
    W_bg = np.ones_like(vals, dtype=float)

    for i in range(max_iter):
        L_bg = stats.poisson.pmf(vals, mu=lambda_)
        L_sig = stats.norm.pdf(vals, loc=mu, scale=sigma)

        Z_bg = pi_bg * L_bg
        Z_sig = (1 - pi_bg) * L_sig
        Z_total = Z_bg + Z_sig
        Z_total[Z_total == 0] = 1e-12

        W_bg = Z_bg / Z_total
        W_sig = Z_sig / Z_total

        N_bg = np.sum(freqs * W_bg)
        N_sig = np.sum(freqs * W_sig)

        if N_bg > 0:
            pi_bg = N_bg / np.sum(freqs)
            lambda_ = np.sum(freqs * W_bg * vals) / N_bg
        if N_sig > 0:
            mu = np.sum(freqs * W_sig * vals) / N_sig
            var_sig = np.sum(freqs * W_sig * (vals - mu) ** 2) / N_sig
            sigma = np.sqrt(var_sig) if var_sig > 0 else 0.1

        lambda_ = max(lambda_, 0.1)
        sigma = max(sigma, 0.1)

        new_log_likelihood = np.sum(freqs * np.log(Z_total))
        if np.abs(new_log_likelihood - log_likelihood) < tol:
            break
        log_likelihood = new_log_likelihood

    crossover_vals = vals[W_sig > W_bg]
    if len(crossover_vals) > 0 and mu > lambda_:
        decision_threshold = max(2, int(np.min(crossover_vals)))
    else:
        decision_threshold = 5

    return {
        "pi_bg": float(pi_bg),
        "lambda_bg": float(lambda_),
        "mu_sig": float(mu),
        "sigma_sig": float(sigma),
        "decision_threshold": decision_threshold,
    }


def call_probe_lists(
    adata, method="gaussian_poisson", count_threshold=None, return_diagnostics=False
):
    """Call all probe features above threshold for each cell."""
    if "guides" not in adata.obsm:
        return (None, None) if return_diagnostics else None

    guide_matrix = adata.obsm["guides"]
    if sp.issparse(guide_matrix):
        guide_matrix = guide_matrix.tocsr()

    mixture_res = None
    if method == "gaussian_poisson":
        if sp.issparse(guide_matrix):
            non_zeros = guide_matrix.data
        else:
            non_zeros = guide_matrix[guide_matrix > 0]

        if len(non_zeros) > 0:
            mixture_res = fit_poisson_gaussian_mixture(non_zeros)
            count_threshold = mixture_res["decision_threshold"]
        else:
            count_threshold = 5
    elif method == "threshold":
        if count_threshold is None:
            raise ValueError("count_threshold is required when method='threshold'")
    else:
        raise ValueError(f"Unsupported guide calling method: {method}")

    guide_names = np.array(
        adata.uns.get(
            "guide_names", [f"guide_{i}" for i in range(guide_matrix.shape[1])]
        )
    )
    total_guide_umis = np.asarray(guide_matrix.sum(axis=1)).ravel()
    positive_guide_counts = []
    calls = []

    for i in range(guide_matrix.shape[0]):
        if sp.issparse(guide_matrix):
            row = guide_matrix.getrow(i)
            keep = row.data >= count_threshold
            positive_idx = row.indices[keep]
        else:
            row = np.asarray(guide_matrix[i]).ravel()
            positive_idx = np.where(row >= count_threshold)[0]

        positive_guide_counts.append(len(positive_idx))
        if len(positive_idx) == 0:
            calls.append("None")
        else:
            names = sorted(canonical_guide_id(guide_names[idx]) for idx in positive_idx)
            calls.append("|".join(names))

    calls = pd.Series(calls, index=adata.obs_names)
    positive_guide_counts = np.asarray(positive_guide_counts)
    count_values, count_freqs = np.unique(positive_guide_counts, return_counts=True)
    positive_cell_mask = positive_guide_counts > 0
    gene_calls = build_gene_call_table(calls)

    diagnostics = {
        "method": method,
        "count_threshold": int(count_threshold),
        "n_cells": int(guide_matrix.shape[0]),
        "n_guides": int(guide_matrix.shape[1]),
        "n_cells_with_any_called_probe": int((positive_guide_counts >= 1).sum()),
        "n_cells_with_two_or_more_called_probes": int(
            (positive_guide_counts >= 2).sum()
        ),
        "n_cells_with_any_called_gene": int((gene_calls["n_genes"] >= 1).sum()),
        "called_probe_count_distribution": {
            str(int(count)): int(freq) for count, freq in zip(count_values, count_freqs)
        },
        "gene_call_outcome_counts": outcome_count_dict(gene_calls["outcome"]),
        "total_guide_umi_distribution_all_cells": summarize_numeric(total_guide_umis),
        "total_guide_umi_distribution_cells_with_any_called_probe": summarize_numeric(
            total_guide_umis[positive_cell_mask]
        ),
    }
    if mixture_res is not None:
        diagnostics["gaussian_poisson_parameters"] = mixture_res

    if return_diagnostics:
        return calls, diagnostics
    return calls


def compare_perturbations(adata_cur, adata_rep, common_cells):
    """Compare perturbation assignments at gene level, after probe calling."""
    cur_pert_col = next(
        (c for c in ["sgID_AB", "perturbation"] if c in adata_cur.obs.columns), None
    )

    if not cur_pert_col:
        log_record(
            "perturbation_comparison_skipped",
            reason="curated_perturbation_column_missing",
            candidate_columns=["sgID_AB", "perturbation"],
        )
        return None

    if "guides" not in adata_rep.obsm:
        log_record(
            "perturbation_comparison_skipped",
            reason="reprocessed_guide_matrix_missing",
            expected_obsm_key="guides",
        )
        return None

    p_cur = (
        adata_cur.obs.loc[common_cells, cur_pert_col]
        .astype(str)
        .apply(canonical_label_for_display)
    )
    cur_calls = build_gene_call_table(p_cur)
    curated_summary = {
        "perturbation_column": cur_pert_col,
        "n_common_cells": int(len(common_cells)),
        "n_cells_with_any_probe": int((cur_calls["n_probes"] > 0).sum()),
        "n_cells_with_any_gene": int((cur_calls["n_genes"] > 0).sum()),
        "outcome_counts": outcome_count_dict(cur_calls["outcome"]),
        "probe_count_distribution": {
            str(int(count)): int(freq)
            for count, freq in cur_calls["n_probes"].value_counts().sort_index().items()
        },
    }
    log_record("curated_gene_call_summary", **curated_summary)

    model_specs = [
        {
            "model": f"threshold_{threshold}",
            "method": "threshold",
            "count_threshold": threshold,
        }
        for threshold in GUIDE_COUNT_THRESHOLDS
    ]
    model_specs.append({"model": "gaussian_poisson", "method": "gaussian_poisson"})

    model_results = {}
    for spec in model_specs:
        rep_labels, rep_diagnostics = call_probe_lists(
            adata_rep,
            method=spec["method"],
            count_threshold=spec.get("count_threshold"),
            return_diagnostics=True,
        )
        p_rep = (
            rep_labels.loc[common_cells].astype(str).apply(canonical_label_for_display)
        )
        rep_calls = build_gene_call_table(p_rep)

        result = {
            "method": spec["method"],
            "count_threshold": rep_diagnostics["count_threshold"],
            "probe_call_diagnostics": rep_diagnostics,
            "reprocessed_outcome_counts": outcome_count_dict(rep_calls["outcome"]),
            "outcome_matrix_rows_curated_columns_reprocessed": outcome_matrix_dict(
                cur_calls["outcome"], rep_calls["outcome"]
            ),
            "single_gene_match": single_gene_match_summary(cur_calls, rep_calls),
        }
        model_results[spec["model"]] = result
        log_record("perturbation_gene_comparison", model=spec["model"], **result)

    return {
        "curated": curated_summary,
        "models": model_results,
        "outcome_definitions": {
            "0_genes": "no gene-targeting probes called",
            "1_gene_1_probe": "one gene called from one gene-targeting probe",
            "1_gene_2_probes": "one gene called from two or more gene-targeting probes",
            ">1_gene": "more than one gene called",
        },
    }


# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================
os.makedirs("comparison_results", exist_ok=True)
log_record(
    "input_paths",
    curated_h5ad_path=CURATED_H5AD_PATH,
    reprocessed_h5ad_path=REPROCESSED_H5AD_PATH,
)
adata_cur = filter_unique_barcodes(sc.read_h5ad(CURATED_H5AD_PATH), "Curated")
adata_rep = filter_unique_barcodes(sc.read_h5ad(REPROCESSED_H5AD_PATH), "Reprocessed")

if adata_rep.var_names.str.contains(r"\.").any():
    adata_rep.var_names = adata_rep.var_names.str.split(".").str[0]
    adata_rep.var_names_make_unique()

common_cells = np.intersect1d(adata_cur.obs_names, adata_rep.obs_names)
common_genes = np.intersect1d(adata_cur.var_names, adata_rep.var_names)

overlap_summary = {
    "common_cells": int(len(common_cells)),
    "common_genes": int(len(common_genes)),
    "curated_cells": int(adata_cur.n_obs),
    "curated_genes": int(adata_cur.n_vars),
    "reprocessed_cells": int(adata_rep.n_obs),
    "reprocessed_genes": int(adata_rep.n_vars),
    "common_cells_pct_of_curated": (
        len(common_cells) / adata_cur.n_obs * 100 if adata_cur.n_obs else 0.0
    ),
    "common_genes_pct_of_curated": (
        len(common_genes) / adata_cur.n_vars * 100 if adata_cur.n_vars else 0.0
    ),
}
log_record("overlap_summary", **overlap_summary)

cur_sub = preprocess_adata(adata_cur[common_cells, common_genes].copy(), "Curated")
rep_sub = preprocess_adata(adata_rep[common_cells, common_genes].copy(), "Reprocessed")

results_summary = {
    "input_paths": {
        "curated_h5ad_path": CURATED_H5AD_PATH,
        "reprocessed_h5ad_path": REPROCESSED_H5AD_PATH,
    },
    "overlap": overlap_summary,
    "cell_metrics": {},
    "gene_metrics": {},
}

# --- Overlap Distributions (3 Histograms) ---
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

# 1. Cell UMI Counts distribution
rep_total_counts = np.array(adata_rep.X.sum(axis=1)).flatten()
common_cell_mask = adata_rep.obs_names.isin(common_cells)
sns.histplot(
    rep_total_counts[common_cell_mask],
    label="Common Cells",
    color="blue",
    alpha=0.5,
    log_scale=True,
    ax=ax1,
)
sns.histplot(
    rep_total_counts[~common_cell_mask],
    label="Reprocessed-only",
    color="orange",
    alpha=0.5,
    log_scale=True,
    ax=ax1,
)
ax1.set_title("Total UMIs per Cell")
ax1.legend()

# 2. Gene n_cells distribution (how many cells have each gene)
rep_gene_sums = np.array((adata_rep.X > 0).sum(axis=0)).flatten()
common_gene_mask = adata_rep.var_names.isin(common_genes)
sns.histplot(
    rep_gene_sums[common_gene_mask],
    label="Common Genes",
    color="blue",
    alpha=0.5,
    log_scale=True,
    ax=ax2,
)
sns.histplot(
    rep_gene_sums[~common_gene_mask],
    label="Reprocessed-only",
    color="orange",
    alpha=0.5,
    log_scale=True,
    ax=ax2,
)
ax2.set_title("Detection Freq per Gene")
ax2.legend()

# 3. Distinct Genes per Cell distribution
rep_n_genes = np.array((adata_rep.X > 0).sum(axis=1)).flatten()
sns.histplot(
    rep_n_genes[common_cell_mask],
    label="Common Cells",
    color="blue",
    alpha=0.5,
    log_scale=True,
    ax=ax3,
)
sns.histplot(
    rep_n_genes[~common_cell_mask],
    label="Reprocessed-only",
    color="orange",
    alpha=0.5,
    log_scale=True,
    ax=ax3,
)
ax3.set_title("Distinct Genes per Cell")
ax3.legend()

plt.suptitle(
    "Comparison of Global Distribution Overlaps", fontsize=16, fontweight="bold"
)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("comparison_results/overlap_distributions.png")
plt.show()
plt.close()

# --- Cell-wise Scatter Plots ---
results_summary["cell_metrics"]["total_counts"] = plot_scatter_comparison(
    pd.DataFrame(
        {"cur": cur_sub.obs["total_counts"], "rep": rep_sub.obs["total_counts"]}
    ),
    "cur",
    "rep",
    "Total UMI Counts",
    "Original",
    "Reprocessed",
    "comparison_results/counts_comparison.png",
    log_scale=True,
    deviation_on_log=True,
    description="Each point is one shared cell. The x-axis shows the original study total UMI count; the y-axis shows the reprocessed total UMI count. Values are plotted on log-scaled axes after adding 1.",
)
log_record(
    "comparison_metric",
    metric_group="cell_metrics",
    metric="total_counts",
    **results_summary["cell_metrics"]["total_counts"],
)

results_summary["cell_metrics"]["n_genes"] = plot_scatter_comparison(
    pd.DataFrame(
        {
            "cur": cur_sub.obs["n_genes_by_counts"],
            "rep": rep_sub.obs["n_genes_by_counts"],
        }
    ),
    "cur",
    "rep",
    "Number of Detected Genes",
    "Original",
    "Reprocessed",
    "comparison_results/genes_comparison.png",
    description="Each point is one shared cell. The x-axis shows the number of genes detected in the original study; the y-axis shows the number of genes detected in the reprocessed data.",
)
log_record(
    "comparison_metric",
    metric_group="cell_metrics",
    metric="n_genes",
    **results_summary["cell_metrics"]["n_genes"],
)

# --- Gene-wise Scatter Plots (Restored) ---
gene_metrics_df = pd.DataFrame(
    {
        "mean_cur": cur_sub.var["mean_counts"],
        "mean_rep": rep_sub.var["mean_counts"],
        "dropout_cur": 100
        - (
            np.array((cur_sub.layers["counts"] > 0).sum(axis=0)).flatten()
            / cur_sub.n_obs
            * 100
        ),
        "dropout_rep": 100
        - (
            np.array((rep_sub.layers["counts"] > 0).sum(axis=0)).flatten()
            / rep_sub.n_obs
            * 100
        ),
    }
)

results_summary["gene_metrics"]["mean_expression"] = plot_scatter_comparison(
    gene_metrics_df,
    "mean_cur",
    "mean_rep",
    "Mean Gene Expression",
    "Original",
    "Reprocessed",
    "comparison_results/gene_expression_mean.png",
    log_scale=True,
    description="Each point is one shared gene. The x-axis shows mean expression across shared cells in the original study; the y-axis shows mean expression across shared cells in the reprocessed data. Values are plotted on log-scaled axes after adding 1.",
)
log_record(
    "comparison_metric",
    metric_group="gene_metrics",
    metric="mean_expression",
    **results_summary["gene_metrics"]["mean_expression"],
)

results_summary["gene_metrics"]["dropout_rate"] = plot_scatter_comparison(
    gene_metrics_df,
    "dropout_cur",
    "dropout_rep",
    "Gene Dropout Rate (%)",
    "Original",
    "Reprocessed",
    "comparison_results/sparsity_comparison.png",
    description="Each point is one shared gene. The x-axis shows the percentage of shared cells with zero counts in the original study; the y-axis shows the percentage of shared cells with zero counts in the reprocessed data.",
)
log_record(
    "comparison_metric",
    metric_group="gene_metrics",
    metric="dropout_rate",
    **results_summary["gene_metrics"]["dropout_rate"],
)

# --- Cell-wise Correlation (Restored Summary) ---


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
    "mean": float(np.nanmean(cell_corrs)),
}
log_record("cell_wise_correlation", **results_summary["cell_wise_corr"])

# --- Perturbations ---
results_summary["perturbation"] = compare_perturbations(
    adata_cur, adata_rep, common_cells
)

# --- Final Report ---
with open("comparison_results/summary_report.txt", "w") as f:
    f.write(json.dumps(results_summary, indent=4, default=json_default))
log_record("summary_report_written", path="comparison_results/summary_report.txt")

perturbation_summary = results_summary["perturbation"]
if perturbation_summary and "gaussian_poisson" in perturbation_summary["models"]:
    gaussian_poisson_match = perturbation_summary["models"]["gaussian_poisson"][
        "single_gene_match"
    ]
    gaussian_poisson_match_value = (
        f"{gaussian_poisson_match['pct_same_gene']:.2f}% "
        f"({gaussian_poisson_match['n_same_gene']}/"
        f"{gaussian_poisson_match['n_cells_both_single_gene']})"
    )
else:
    gaussian_poisson_match_value = "N/A"

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
                "Metric": "Gaussian-Poisson single-gene match",
                "Value": gaussian_poisson_match_value,
            },
        ]
    )
)

# ==============================================================================
# 4. SUPPLEMENTARY VISUALIZATIONS (For Jupyter Notebook / Panel Presentation)
# ==============================================================================
# The following code block is designed to be copy-pasted into a new Jupyter Cell.
# It assumes the main execution block above has already run and variables like
# `cur_sub`, `rep_sub`, and `cell_corrs` are in memory.

import os
import scipy.sparse as sp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import scanpy as sc

os.makedirs("comparison_results/supplementary", exist_ok=True)

# ------------------------------------------------------------------------------
# Graph 1: Distribution of Cell-wise Correlations
# ------------------------------------------------------------------------------
plt.figure(figsize=(9, 6))
sns.histplot(cell_corrs, bins=100, kde=True, color="purple", alpha=0.4)
median_val = np.nanmedian(cell_corrs)
plt.axvline(
    median_val,
    color="red",
    linestyle="dashed",
    linewidth=2,
    label=f"Median: {median_val:.3f}",
)
plt.title(
    "Preservation of Single-Cell Identity\n(Cell-wise Pearson Correlation)",
    fontsize=16,
    fontweight="bold",
    pad=15,
)
plt.xlabel("Pearson Correlation (Original vs Reprocessed cell)", fontsize=12)
plt.ylabel("Number of Cells", fontsize=12)
plt.legend()
desc1 = (
    "This distribution shows the correlation of the full expression profile for "
    "each individual cell against its exact counterpart in the reprocessed dataset. "
    "A strong peak near 1.0 confirms that single-cell identities are highly preserved."
)
plt.figtext(
    0.5,
    -0.05,
    desc1,
    wrap=True,
    horizontalalignment="center",
    fontsize=10,
    style="italic",
)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig(
    "comparison_results/supplementary/cellwise_correlation_dist.png",
    bbox_inches="tight",
    dpi=300,
)
plt.show()

# ------------------------------------------------------------------------------
# Graph 2: Principal Component Alignment Heatmap
# ------------------------------------------------------------------------------
# Calculate PCA independently for both datasets to ensure structure is inherent
sc.tl.pca(cur_sub, n_comps=10)
sc.tl.pca(rep_sub, n_comps=10)

pc_corr = np.zeros((10, 10))
for i in range(10):
    for j in range(10):
        # We use absolute correlation because the sign (direction) of a PC is arbitrary
        corr, _ = stats.pearsonr(
            cur_sub.obsm["X_pca"][:, i], rep_sub.obsm["X_pca"][:, j]
        )
        pc_corr[i, j] = np.abs(corr)

off_diagonal = pc_corr[~np.eye(pc_corr.shape[0], dtype=bool)]
log_record(
    "pca_alignment",
    diagonal_abs_correlations=[float(value) for value in np.diag(pc_corr)],
    max_off_diagonal_abs_correlation=float(np.max(off_diagonal)),
)

plt.figure(figsize=(9, 7))
sns.heatmap(
    pc_corr,
    annot=True,
    cmap="YlGnBu",
    fmt=".2f",
    vmin=0,
    vmax=1,
    xticklabels=[f"Rep PC{i+1}" for i in range(10)],
    yticklabels=[f"Cur PC{i+1}" for i in range(10)],
)
plt.title(
    "Latent Structural Integrity\n(Alignment of Top 10 Principal Components)",
    fontsize=16,
    fontweight="bold",
    pad=15,
)
desc2 = (
    "Absolute Pearson correlation between the top 10 independent Principal Components of "
    "both datasets. A strong diagonal demonstrates that the global biological covariance "
    "structure and major axes of variation remain intact."
)
plt.figtext(
    0.5,
    -0.05,
    desc2,
    wrap=True,
    horizontalalignment="center",
    fontsize=10,
    style="italic",
)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig(
    "comparison_results/supplementary/pca_alignment_heatmap.png",
    bbox_inches="tight",
    dpi=300,
)
plt.show()

# ------------------------------------------------------------------------------
# Graph 3: Gene Variance (Dispersion) Scatter Plot
# ------------------------------------------------------------------------------


def calc_variance(matrix):
    if sp.issparse(matrix):
        # E[X^2] - (E[X])^2 for sparse matrices to avoid dense memory explosion
        mean = matrix.mean(axis=0).A.squeeze()
        sq_mean = matrix.multiply(matrix).mean(axis=0).A.squeeze()
        return sq_mean - (mean**2)
    else:
        return np.var(matrix, axis=0)


cur_var = calc_variance(cur_sub.X)
rep_var = calc_variance(rep_sub.X)

plt.figure(figsize=(9, 9))
# Add 1e-4 pseudocount for log-scale plotting
plt.scatter(cur_var + 1e-4, rep_var + 1e-4, alpha=0.3, s=15, color="darkgreen")

# Identity line
min_val = min(np.min(cur_var), np.min(rep_var)) + 1e-4
max_val = max(np.max(cur_var), np.max(rep_var)) + 1e-4
plt.plot(
    [min_val, max_val], [min_val, max_val], "r--", linewidth=2, label="Identity (y=x)"
)

pearson_var, _ = stats.pearsonr(cur_var, rep_var)
spearman_var, _ = stats.spearmanr(cur_var, rep_var)
log_record(
    "gene_variance_metric",
    pearson=float(pearson_var),
    spearman=float(spearman_var),
)

plt.xscale("log")
plt.yscale("log")
plt.title(
    "Statistical Noise Preservation\n(Gene Variance Comparison)",
    fontsize=16,
    fontweight="bold",
    pad=35,
)
plt.text(
    0.5,
    1.02,
    f"Pearson r = {pearson_var:.4f} | Spearman rho = {spearman_var:.4f}",
    transform=plt.gca().transAxes,
    ha="center",
    va="bottom",
    fontsize=12,
    style="italic",
)
plt.xlabel("Gene Variance in Original Data (+ 1e-4)", fontsize=13)
plt.ylabel("Gene Variance in Reprocessed Data (+ 1e-4)", fontsize=13)
plt.legend(loc="upper left")

desc3 = (
    "Compares the variance of each gene across all cells. High correlation indicates "
    "that the biological overdispersion and noise characteristics required for rigorous "
    "differential expression modeling (like DESeq2/TRADE) are fully preserved."
)
plt.figtext(
    0.5,
    -0.05,
    desc3,
    wrap=True,
    horizontalalignment="center",
    fontsize=10,
    style="italic",
)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig(
    "comparison_results/supplementary/gene_variance_scatter.png",
    bbox_inches="tight",
    dpi=300,
)
plt.show()
