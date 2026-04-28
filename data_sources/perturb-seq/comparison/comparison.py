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
    print(
        f"[{name}] Filtered out {filtered_count} cells due to barcode collisions ({(filtered_count/original_count)*100:.2f}%). Remaining: {adata.n_obs}"
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
    print(f"Starting standardized preprocessing for {name}...")
    is_raw = get_raw_counts(adata)
    if not is_raw:
        print(f"  [{name}] WARNING: Data does not appear to be raw counts.")

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
    """Extract the target label used to validate same-gene dual-guide calls."""
    guide_name = canonical_guide_id(guide_name)
    if guide_name.startswith("non-targeting"):
        return "non-targeting"
    return guide_name.split("_", 1)[0]


def guide_target_names(guide_name):
    return {guide_target_name(alias) for alias in guide_aliases(guide_name)}


def label_guide_alias_sets(label):
    label = str(label).strip()
    if label.lower() in {"", "none", "nan"}:
        return []
    return [set(guide_aliases(part)) for part in label.split("|") if part.strip()]


def label_target_name(label):
    """Extract the target label from a pipe-delimited guide assignment."""
    label = str(label).strip()
    if label.lower() in {"", "none", "nan"}:
        return "None"

    per_guide_targets = [
        guide_target_names(part) for part in label.split("|") if part.strip()
    ]
    if not per_guide_targets:
        return "None"

    shared_targets = set.intersection(*per_guide_targets)
    targets = shared_targets if shared_targets else set.union(*per_guide_targets)
    if len(targets) == 1:
        return next(iter(targets))
    return "mixed:" + "|".join(sorted(targets))


def labels_exact_match(cur_label, rep_label):
    cur_parts = label_guide_alias_sets(cur_label)
    rep_parts = label_guide_alias_sets(rep_label)
    if len(cur_parts) != len(rep_parts):
        return False
    if len(cur_parts) == 0:
        return True
    if len(cur_parts) == 1:
        return bool(cur_parts[0] & rep_parts[0])
    if len(cur_parts) == 2:
        return (
            bool(cur_parts[0] & rep_parts[0]) and bool(cur_parts[1] & rep_parts[1])
        ) or (bool(cur_parts[0] & rep_parts[1]) and bool(cur_parts[1] & rep_parts[0]))

    unmatched = list(rep_parts)
    for cur_aliases in cur_parts:
        for i, rep_aliases in enumerate(unmatched):
            if cur_aliases & rep_aliases:
                unmatched.pop(i)
                break
        else:
            return False
    return True


def labels_target_match(cur_label, rep_label):
    cur_target = label_target_name(cur_label)
    rep_target = label_target_name(rep_label)
    if cur_target == "None" or rep_target == "None":
        return cur_target == rep_target
    cur_targets = set(cur_target.removeprefix("mixed:").split("|"))
    rep_targets = set(rep_target.removeprefix("mixed:").split("|"))
    return bool(cur_targets & rep_targets)


def canonical_label_for_display(label):
    label = str(label).strip()
    if label.lower() in {"", "none", "nan"}:
        return "None"
    parts = [canonical_guide_id(part) for part in label.split("|") if part.strip()]
    return "|".join(sorted(parts)) if parts else "None"


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


def build_guide_alias_index(guide_names):
    alias_to_indices = {}
    for idx, guide_name in enumerate(guide_names):
        for alias in guide_aliases(guide_name):
            alias_to_indices.setdefault(alias, set()).add(idx)
    return {alias: sorted(indices) for alias, indices in alias_to_indices.items()}


def label_feature_index_sets(label, alias_to_indices):
    index_sets = []
    for guide_part in str(label).split("|"):
        guide_part = guide_part.strip()
        if not guide_part or guide_part.lower() in {"none", "nan"}:
            continue

        indices = set()
        for alias in guide_aliases(guide_part):
            indices.update(alias_to_indices.get(alias, []))
        index_sets.append(indices)
    return index_sets


def max_count_for_indices(row, indices):
    if not indices:
        return np.nan

    if sp.issparse(row):
        data_by_col = {
            int(col): float(value) for col, value in zip(row.indices, row.data)
        }
        return max(data_by_col.get(int(idx), 0.0) for idx in indices)

    row = np.asarray(row).ravel()
    return float(np.max(row[list(indices)]))


def summarize_curated_pair_counts(count_a, count_b, mask, thresholds):
    mask = np.asarray(mask, dtype=bool)
    resolved = mask & np.isfinite(count_a) & np.isfinite(count_b)
    a = count_a[resolved]
    b = count_b[resolved]
    min_counts = np.minimum(a, b)
    max_counts = np.maximum(a, b)

    summary = {
        "n_cells": int(mask.sum()),
        "n_cells_with_both_curated_guides_in_reference": int(resolved.sum()),
        "min_curated_guide_count_distribution": summarize_numeric(min_counts),
        "max_curated_guide_count_distribution": summarize_numeric(max_counts),
        "total_curated_pair_count_distribution": summarize_numeric(a + b),
        "thresholds": {},
    }

    for threshold in thresholds:
        both = min_counts >= threshold
        one = (max_counts >= threshold) & ~both
        neither = max_counts < threshold
        summary["thresholds"][str(threshold)] = {
            "n_cells_with_both_curated_guides_ge_threshold": int(both.sum()),
            "n_cells_with_exactly_one_curated_guide_ge_threshold": int(one.sum()),
            "n_cells_with_neither_curated_guide_ge_threshold": int(neither.sum()),
        }

    return summary


def curated_guide_count_diagnostics(
    adata_rep, common_cells, curated_labels, valid_cur, valid_rep
):
    if "guides" not in adata_rep.obsm:
        return None

    guide_matrix = adata_rep.obsm["guides"]
    if sp.issparse(guide_matrix):
        guide_matrix = guide_matrix.tocsr()

    guide_names = np.array(
        adata_rep.uns.get(
            "guide_names", [f"guide_{i}" for i in range(guide_matrix.shape[1])]
        )
    )
    alias_to_indices = build_guide_alias_index(guide_names)
    obs_to_pos = pd.Series(np.arange(adata_rep.n_obs), index=adata_rep.obs_names)
    row_positions = obs_to_pos.loc[common_cells].to_numpy()

    count_a = np.full(len(common_cells), np.nan)
    count_b = np.full(len(common_cells), np.nan)
    resolved_parts = np.zeros(len(common_cells), dtype=int)

    for i, (row_pos, label) in enumerate(zip(row_positions, curated_labels)):
        index_sets = label_feature_index_sets(label, alias_to_indices)
        if len(index_sets) < 2:
            continue

        first, second = index_sets[:2]
        resolved_parts[i] = int(bool(first)) + int(bool(second))
        row = (
            guide_matrix.getrow(row_pos)
            if sp.issparse(guide_matrix)
            else guide_matrix[row_pos, :]
        )
        count_a[i] = max_count_for_indices(row, first)
        count_b[i] = max_count_for_indices(row, second)

    thresholds = [1, 2, 3, 4, 5, 10]
    valid_cur_np = np.asarray(valid_cur, dtype=bool)
    valid_rep_np = np.asarray(valid_rep, dtype=bool)
    missing_rep_np = valid_cur_np & ~valid_rep_np
    resolved_values, resolved_counts = np.unique(
        resolved_parts[valid_cur_np], return_counts=True
    )

    return {
        "curated_label_feature_resolved_distribution": {
            str(int(value)): int(count)
            for value, count in zip(resolved_values, resolved_counts)
        },
        "all_curated_guide_cells": summarize_curated_pair_counts(
            count_a, count_b, valid_cur_np, thresholds
        ),
        "curated_cells_missing_reprocessed_call": summarize_curated_pair_counts(
            count_a, count_b, missing_rep_np, thresholds
        ),
    }


def call_guides(
    adata,
    count_threshold=5,
    require_dual_same_target=True,
    max_guides=2,
    return_diagnostics=False,
):
    """Call guides, optionally requiring two above-threshold guides to the same target."""
    if "guides" not in adata.obsm:
        return (None, None) if return_diagnostics else None

    guide_matrix = adata.obsm["guides"]
    if sp.issparse(guide_matrix):
        guide_matrix = guide_matrix.tocsr()

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
            positive_counts = row.data[keep]
        else:
            row = np.asarray(guide_matrix[i]).ravel()
            positive_idx = np.where(row >= count_threshold)[0]
            positive_counts = row[positive_idx]

        positive_guide_counts.append(len(positive_idx))
        if len(positive_idx) == 0:
            calls.append("None")
        elif require_dual_same_target:
            candidate_pairs = []
            entries = list(zip(positive_idx, positive_counts))
            for left_pos, (left_idx, left_count) in enumerate(entries):
                left_targets = guide_target_names(guide_names[left_idx])
                for right_idx, right_count in entries[left_pos + 1 :]:
                    shared_targets = left_targets & guide_target_names(
                        guide_names[right_idx]
                    )
                    if not shared_targets:
                        continue
                    score = left_count + right_count
                    target = sorted(shared_targets)[0]
                    candidate_pairs.append((score, target, [left_idx, right_idx]))

            if not candidate_pairs:
                calls.append("None")
                continue

            _, _, top_idx = sorted(
                candidate_pairs, key=lambda item: (-item[0], item[1])
            )[0]
            names = sorted(guide_names[top_idx])
            calls.append("|".join(names))
        else:
            top_idx = positive_idx[np.argsort(positive_counts)[::-1][:max_guides]]
            names = sorted(guide_names[top_idx])
            calls.append("|".join(names))

    calls = pd.Series(calls, index=adata.obs_names)
    positive_guide_counts = np.asarray(positive_guide_counts)
    count_values, count_freqs = np.unique(positive_guide_counts, return_counts=True)
    positive_cell_mask = positive_guide_counts > 0

    diagnostics = {
        "count_threshold": int(count_threshold),
        "require_dual_same_target": bool(require_dual_same_target),
        "n_cells": int(guide_matrix.shape[0]),
        "n_guides": int(guide_matrix.shape[1]),
        "n_cells_with_any_positive_guide": int((positive_guide_counts >= 1).sum()),
        "n_cells_with_two_or_more_positive_guides": int(
            (positive_guide_counts >= 2).sum()
        ),
        "n_cells_with_dual_same_target_call": int((calls != "None").sum()),
        "positive_guide_count_distribution": {
            str(int(count)): int(freq) for count, freq in zip(count_values, count_freqs)
        },
        "total_guide_umi_distribution_all_cells": summarize_numeric(total_guide_umis),
        "total_guide_umi_distribution_cells_with_any_positive_guide": summarize_numeric(
            total_guide_umis[positive_cell_mask]
        ),
    }

    if return_diagnostics:
        return calls, diagnostics
    return calls


def compare_perturbations(adata_cur, adata_rep, common_cells):
    """Compares perturbation assignments with deep mismatch diagnostics."""
    print("Comparing perturbation assignments...")

    cur_pert_col = next(
        (c for c in ["sgID_AB", "perturbation"] if c in adata_cur.obs.columns), None
    )
    if "guides" in adata_rep.obsm:
        rep_labels, rep_diagnostics = call_guides(adata_rep, return_diagnostics=True)
    else:
        rep_labels, rep_diagnostics = None, None

    if cur_pert_col and rep_labels is not None:
        p_cur = adata_cur.obs.loc[common_cells, cur_pert_col].astype(str)
        p_rep = rep_labels.loc[common_cells].astype(str)

        p_cur = p_cur.apply(canonical_label_for_display)
        p_rep = p_rep.apply(canonical_label_for_display)

        # Deep Diagnostics
        valid_cur = ~p_cur.str.lower().isin(["", "none", "nan"])
        valid_rep = ~p_rep.str.lower().isin(["", "none", "nan"])
        print(
            f"  Curated cells with guides: {valid_cur.sum()} / {len(common_cells)} ({valid_cur.mean():.1%})"
        )
        print(
            "  Reprocessed cells with valid dual same-target guide calls: "
            f"{valid_rep.sum()} / {len(common_cells)} ({valid_rep.mean():.1%})"
        )
        print("  Reprocessed guide matrix diagnostics:")
        print(textwrap.indent(json.dumps(rep_diagnostics, indent=2), "    "))

        exact_matches = pd.Series(
            [
                labels_exact_match(cur_label, rep_label)
                for cur_label, rep_label in zip(p_cur, p_rep)
            ],
            index=p_cur.index,
        )
        target_matches = pd.Series(
            [
                labels_target_match(cur_label, rep_label)
                for cur_label, rep_label in zip(p_cur, p_rep)
            ],
            index=p_cur.index,
        )

        missing_rep = valid_cur & ~valid_rep
        n_valid_cur = int(valid_cur.sum())
        missing_rate = (int(missing_rep.sum()) / n_valid_cur) if n_valid_cur else 0.0
        missing_call_diagnostics = curated_guide_count_diagnostics(
            adata_rep, common_cells, p_cur, valid_cur, valid_rep
        )
        print(
            "  Curated-guide cells missing a valid Reprocessed dual same-target call: "
            f"{missing_rep.sum()} / {n_valid_cur} ({missing_rate:.1%})"
        )
        print(
            f"  Alias-aware exact pair matches: {exact_matches.sum()} / {len(common_cells)} ({exact_matches.mean():.1%})"
        )
        print(
            f"  Alias-aware target matches: {target_matches.sum()} / {len(common_cells)} ({target_matches.mean():.1%})"
        )
        print("  Curated-guide count diagnostics:")
        print(textwrap.indent(json.dumps(missing_call_diagnostics, indent=2), "    "))

        mismatches = np.where(valid_rep & ~exact_matches.to_numpy())[0]
        if len(mismatches) > 0:
            print(
                f"  Mismatched cells diagnostics (Sample of mismatches where Reprocessed has guides):"
            )
            print(f"  Total such mismatches: {len(mismatches)}")
            for idx in mismatches[:10]:
                print(
                    f"    {common_cells[idx]}: Curated='{p_cur.iloc[idx]}' vs Rep='{p_rep.iloc[idx]}'"
                )

        accuracy = exact_matches.mean()
        target_accuracy = target_matches.mean()
        overlap_df = pd.DataFrame({"Curated": p_cur, "Reprocessed": p_rep})
        top_perts = p_cur.value_counts().head(20).index
        sub_df = overlap_df[overlap_df["Curated"].isin(top_perts)]
        ct = pd.crosstab(sub_df["Reprocessed"], sub_df["Curated"])

        plt.figure(figsize=(15, 13))
        ax = sns.heatmap(
            ct,
            annot=False,
            cmap="YlGnBu",
            cbar_kws={"label": "Number of shared cells"},
        )
        plt.suptitle(
            "Perturbation Confusion Matrix (Top 20)",
            fontsize=18,
            fontweight="bold",
            y=0.98,
        )

        desc = (
            "Columns are original study / Curated guide-pair assignments; rows are Reprocessed dual same-target calls. "
            "Each tile is the number of shared cells with that assignment pair. "
            "The diagonal is exact guide-pair agreement; the Reprocessed 'None' row contains cells where no valid dual same-target call passed the count threshold."
        )
        wrapped_desc = "\n".join(textwrap.wrap(desc, width=110))
        plt.title(
            f"Exact Pair Match: {accuracy:.4%} | Target Match: {target_accuracy:.4%}\n{wrapped_desc}",
            fontsize=10,
            pad=15,
            style="italic",
            loc="center",
        )

        ax.set_xlabel("Original study / Curated assignment")
        ax.set_ylabel("Reprocessed assignment")
        plt.xticks(rotation=45, ha="right", fontsize=7)
        plt.yticks(fontsize=7)
        plt.tight_layout(rect=[0, 0.05, 1, 0.94])
        plt.savefig("comparison_results/perturbation_confusion_matrix.png")
        plt.show()
        plt.close()

        return {
            "accuracy": float(accuracy),
            "target_accuracy": float(target_accuracy),
            "n_rep_with_guides": int(valid_rep.sum()),
            "n_curated_cells_missing_reprocessed_call": int(missing_rep.sum()),
            "curated_guide_count_diagnostics": missing_call_diagnostics,
            "guide_call_diagnostics": rep_diagnostics,
        }
    return None


# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================
os.makedirs("comparison_results", exist_ok=True)
print("Loading and aligning...")
adata_cur = filter_unique_barcodes(
    sc.read_h5ad("GSE264667_jurkat_raw_singlecell_01.h5ad"), "Curated"
)
adata_rep = filter_unique_barcodes(sc.read_h5ad("experiment_final.h5ad"), "Reprocessed")

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
    "mean": float(np.nanmean(cell_corrs)),
}

# --- Perturbations ---
results_summary["perturbation"] = compare_perturbations(
    adata_cur, adata_rep, common_cells
)

# --- Final Report ---
print("\n" + "=" * 50)
print("COMPARISON COMPLETE")
print("=" * 50)

with open("comparison_results/summary_report.txt", "w") as f:
    f.write(json.dumps(results_summary, indent=4))

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
