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
import gzip

# For nice Jupyter rendering
try:
    from IPython.display import display
except ImportError:
    display = print

# Set only this parameter for a particular dataset
DATASET_ID = os.environ.get("PERTURBSEQ_DATASET_ID", "nadig_2025_jurkat")

# Parameters below should not need to be modified between datasets
CURATED_H5AD_PATH = (
    f"/hps/nobackup/mfreeberg/perturb_seq_fastq/source_h5ad/{DATASET_ID}.h5ad"
)
REPROCESSED_H5AD_PATH = (
    f"/hps/nobackup/mfreeberg/perturb_seq_fastq/results/"
    f"{DATASET_ID}/experiment_final.h5ad"
)
FILTERED_REPROCESSED_H5AD_PATH = (
    f"/hps/nobackup/mfreeberg/perturb_seq_fastq/results/"
    f"{DATASET_ID}/experiment_final.filtered.h5ad"
)
CELL_QC_LOWER_QUANTILE = 0.01
CELL_MIN_COUNTS_FLOOR = 1000
CELL_MIN_GENES_FLOOR = 200
GENE_MIN_CELLS_FLOOR = 10
GENE_MIN_CELLS_PCT = 0.01
GENE_CALL_OUTCOMES = ["0_genes", "1_gene_1_probe", "1_gene_2_probes", ">1_gene"]
CONTROL_TARGET_SYMBOL = "non-targeting"
REFERENCE_GTF_PATH = (
    "/hps/nobackup/mfreeberg/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz"
)
CELL_TOTAL_COUNT_COLUMNS = [
    "UMI_count",
    "umi_count",
    "qc_total_counts",
    "total_counts",
    "n_counts",
    "nCount_RNA",
    "n_umi",
    "core_adjusted_UMI_count",
]
CELL_DETECTED_GENE_COLUMNS = [
    "n_genes_by_counts",
    "qc_n_genes_by_counts",
    "n_genes",
    "nFeature_RNA",
    "genes_detected",
    "detected_genes",
]
GENE_MEAN_COLUMNS = ["mean_counts", "mean_expression", "mean"]
GENE_DETECTED_CELL_COLUMNS = [
    "n_cells_by_counts",
    "qc_n_cells_by_counts",
    "n_cells",
    "num_cells_expressed",
]
PERTURBATION_COLUMN_CANDIDATES = [
    "sgID_AB",
    "perturbation",
    "guide_ids",
    "guide_id",
    "sgRNA",
    "sgRNA_ID",
    "gRNA",
    "grna",
    "probe_label",
]


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


def format_pct(value):
    return f"{value:.1f}%"


def format_count_pct(count, total):
    pct = (count / total * 100) if total else 0.0
    return f"{count}/{total} ({pct:.1f}%)"


def format_outcome_counts(counts):
    return ", ".join(
        f"{outcome}: {counts.get(outcome, 0)}" for outcome in GENE_CALL_OUTCOMES
    )


def log_record(event, **fields):
    if event == "input_paths":
        print("Inputs")
        print(f"  Curated H5AD: {fields['curated_h5ad_path']}")
        print(f"  Reprocessed H5AD: {fields['reprocessed_h5ad_path']}")
        print(
            f"  Filtered reprocessed H5AD: {fields['filtered_reprocessed_h5ad_path']}"
        )
    elif event == "barcode_filter":
        print(
            f"Barcode filtering - {fields['dataset']}: "
            f"removed {format_count_pct(fields['removed_collision_cells'], fields['original_cells'])}; "
            f"remaining {fields['remaining_cells']}"
        )
    elif event == "preprocess_input":
        raw_status = "yes" if fields["raw_counts_detected"] else "no"
        print(
            f"Preprocessing input - {fields['dataset']}: "
            f"{fields['cells']} cells, {fields['genes']} genes, raw counts detected: {raw_status}"
        )
        if fields.get("expression_value_kind"):
            print(
                f"  Expression values: {fields['expression_value_kind']}; "
                f"normalization: {fields.get('normalization_action', 'n/a')}"
            )
    elif event == "comparison_metric_source":
        print(f"Metric sources - {fields['dataset']}")
        for metric, source in fields["sources"].items():
            print(f"  {metric}: {source}")
    elif event == "comparison_metric_skipped":
        print(
            f"{fields['metric_group']}.{fields['metric']}: skipped "
            f"({fields['reason']})"
        )
    elif event == "overlap_summary":
        print("Overlap")
        print(
            f"  Cells: {fields['common_cells']} common; "
            f"curated {fields['curated_cells']}, reprocessed {fields['reprocessed_cells']}; "
            f"{format_pct(fields['common_cells_pct_of_curated'])} of curated"
        )
        print(
            f"  Genes: {fields['common_genes']} common; "
            f"curated {fields['curated_genes']}, reprocessed {fields['reprocessed_genes']}; "
            f"{format_pct(fields['common_genes_pct_of_curated'])} of curated"
        )
    elif event == "comparison_metric":
        if fields.get("skipped"):
            print(
                f"{fields['metric_group']}.{fields['metric']}: skipped "
                f"({fields.get('reason', 'not enough comparable values')})"
            )
        else:
            print(
                f"{fields['metric_group']}.{fields['metric']}: "
                f"Pearson {fields['pearson']:.4f}, Spearman {fields['spearman']:.4f}, "
                f"deviant {fields['pct_deviant']:.1f}%"
            )
    elif event == "cell_wise_correlation":
        print(
            f"Cell-wise expression correlation: median {fields['median']:.4f}, "
            f"mean {fields['mean']:.4f}"
        )
    elif event == "curated_gene_call_summary":
        print("Perturbation gene calls")
        print(
            f"  Curated column: {fields['perturbation_column']}; "
            f"cells with probes {format_count_pct(fields['n_cells_with_any_probe'], fields['n_common_cells'])}; "
            f"cells with genes {format_count_pct(fields['n_cells_with_any_gene'], fields['n_common_cells'])}"
        )
        print(f"  Curated outcomes: {format_outcome_counts(fields['outcome_counts'])}")
    elif event == "perturbation_gene_comparison":
        match = fields["single_gene_match"]
        common_summary = fields["reprocessed_common_call_summary"]
        print(
            f"  {fields['model']}: threshold >= {fields['count_threshold']} UMI; "
            f"called probes in "
            f"{format_count_pct(common_summary['n_cells_with_any_probe'], common_summary['n_common_cells'])}; "
            f"called genes in "
            f"{format_count_pct(common_summary['n_cells_with_any_gene'], common_summary['n_common_cells'])}"
        )
        print(
            f"    Outcomes: {format_outcome_counts(fields['reprocessed_outcome_counts'])}"
        )
        print(
            f"    Same single-gene calls: "
            f"{format_count_pct(match['n_same_gene'], match['n_cells_both_single_gene'])}; "
            f"different {format_count_pct(match['n_different_gene'], match['n_cells_both_single_gene'])}"
        )
        if fields.get("plot_path"):
            print(f"    Matrix plot: {fields['plot_path']}")
    elif event == "perturbation_comparison_skipped":
        print(f"Perturbation comparison skipped: {fields['reason']}")
    elif event == "summary_report_written":
        print(f"Summary report: {fields['path']}")
    elif event == "qc_filter":
        print("QC filtering - Reprocessed")
        print(
            f"  Cells: {fields['cells_before']} -> {fields['cells_after']} "
            f"(removed {format_count_pct(fields['cells_removed'], fields['cells_before'])})"
        )
        print(
            f"  Cell thresholds: total counts >= {fields['min_total_counts']}, "
            f"detected genes >= {fields['min_genes_by_counts']}"
        )
        print(
            f"  Genes: {fields['genes_before']} -> {fields['genes_after']} "
            f"(removed {format_count_pct(fields['genes_removed'], fields['genes_before'])})"
        )
        print(
            f"  Gene threshold: detected in >= {fields['min_cells_by_counts']} cells "
            f"({fields['min_cells_pct']:.3f}% of QC-passing cells)"
        )
    elif event == "knockout_annotation":
        print(
            "Knockout annotation - Reprocessed: "
            f"Gaussian-Poisson threshold >= {fields['count_threshold']} UMI; "
            f"single-gene/no-control calls {format_count_pct(fields['n_single_gene_calls'], fields['n_cells'])}; "
            f"control-only calls {fields['n_control_calls']}; "
            f"multi-gene calls {fields['n_multi_gene_calls']}; "
            f"no-guide calls {fields['n_no_guide_calls']}"
        )
    elif event == "gene_symbol_annotation":
        print(
            f"Gene symbol annotation - {fields['dataset']}: "
            f"mapped {format_count_pct(fields['n_mapped'], fields['n_genes'])}; "
            f"source {fields['source']}"
        )
    elif event == "gene_symbol_filter":
        print(
            f"Gene symbol filtering - {fields['dataset']}: "
            f"removed {format_count_pct(fields['n_removed'], fields['n_before'])} "
            "features without gene_name in the reference GTF"
        )
    elif event == "filtered_h5ad_written":
        print(
            f"Filtered H5AD written: {fields['path']} "
            f"({fields['cells']} cells x {fields['genes']} genes; "
            f"compression={fields['compression']})"
        )
    else:
        print(f"{event}: {fields}")


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


def matrix_sum(matrix, axis):
    return np.asarray(matrix.sum(axis=axis)).ravel()


def matrix_nnz(matrix, axis):
    return np.asarray((matrix > 0).sum(axis=axis)).ravel()


def lower_quantile_threshold(values, floor, quantile):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return int(floor)
    return int(max(floor, np.floor(np.quantile(values, quantile))))


def normalize_column_name(value):
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def numeric_series_or_none(series):
    values = pd.to_numeric(series, errors="coerce")
    if values.notna().sum() == 0:
        return None
    return values


def find_numeric_column(
    df, preferred_names, require_nonnegative=True, required_name_terms=None
):
    """Find a numeric metadata column by exact/normalized name, then by heuristics."""
    normalized_to_column = {
        normalize_column_name(column): column for column in df.columns
    }

    for name in preferred_names:
        column = (
            name
            if name in df.columns
            else normalized_to_column.get(normalize_column_name(name))
        )
        if column is None:
            continue
        values = numeric_series_or_none(df[column])
        if values is None:
            continue
        finite = values[np.isfinite(values)]
        if finite.empty:
            continue
        if require_nonnegative and (finite < 0).any():
            continue
        return column, values.astype(float)

    best = None
    best_score = -np.inf
    for column in df.columns:
        values = numeric_series_or_none(df[column])
        if values is None:
            continue
        finite = values[np.isfinite(values)]
        if finite.empty:
            continue
        if require_nonnegative and (finite < 0).any():
            continue

        name = normalize_column_name(column)
        if required_name_terms and not any(
            term in name for term in required_name_terms
        ):
            continue
        score = 0
        if "total" in name:
            score += 3
        if "umi" in name:
            score += 3
        if "count" in name:
            score += 2
        if "gene" in name or "feature" in name or "detected" in name:
            score += 2
        if "percent" in name or "pct" in name or "mito" in name:
            score -= 5
        if name.startswith("z") or "zscore" in name or "scalefactor" in name:
            score -= 5

        if score > best_score:
            best = (column, values.astype(float))
            best_score = score

    if best is not None and best_score > 0:
        return best
    return None, None


def sample_matrix_values(matrix, max_values=200000):
    """Sample matrix values without assuming dense, sparse, or backed storage."""
    if sp.issparse(matrix):
        values = matrix.data
        if values.size == 0:
            return np.array([0.0])
        if values.size <= max_values:
            return np.asarray(values, dtype=float)
        sample_idx = np.linspace(0, values.size - 1, max_values, dtype=np.intp)
        return np.asarray(values[sample_idx], dtype=float)

    values = np.asarray(matrix)
    if values.size == 0:
        return np.array([], dtype=float)
    flat = values.reshape(-1)
    if flat.size <= max_values:
        return np.asarray(flat, dtype=float)
    sample_idx = np.linspace(0, flat.size - 1, max_values, dtype=np.intp)
    return np.asarray(flat[sample_idx], dtype=float)


def profile_expression_matrix(matrix):
    values = sample_matrix_values(matrix)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "value_kind": "unknown",
            "raw_counts_detected": False,
            "min": np.nan,
            "max": np.nan,
            "integer_fraction": 0.0,
            "negative_fraction": 0.0,
            "below_minus_one_fraction": 0.0,
        }

    integer_fraction = float(np.mean(np.isclose(values, np.round(values))))
    negative_fraction = float(np.mean(values < 0))
    below_minus_one_fraction = float(np.mean(values < -1))
    min_value = float(np.min(values))
    max_value = float(np.max(values))

    if negative_fraction > 0:
        value_kind = "signed_transformed"
    elif integer_fraction > 0.999:
        value_kind = "raw_counts"
    else:
        value_kind = "nonnegative_transformed"

    return {
        "value_kind": value_kind,
        "raw_counts_detected": value_kind == "raw_counts",
        "min": min_value,
        "max": max_value,
        "integer_fraction": integer_fraction,
        "negative_fraction": negative_fraction,
        "below_minus_one_fraction": below_minus_one_fraction,
    }


def matrix_mean_vector(matrix, axis):
    return np.asarray(matrix.mean(axis=axis)).ravel()


def matrix_variance_vector(matrix, axis=0):
    if axis != 0:
        raise NotImplementedError("Only per-gene variance is needed here.")
    if sp.issparse(matrix):
        means = matrix_mean_vector(matrix, axis=0)
        mean_squares = matrix_mean_vector(matrix.power(2), axis=0)
        variances = mean_squares - means**2
        variances[variances < 0] = 0
        return variances
    return np.nanvar(np.asarray(matrix), axis=0)


def set_highly_variable_by_variance(adata, n_top_genes):
    n_top = int(min(n_top_genes, adata.n_vars))
    variances = matrix_variance_vector(adata.X, axis=0)
    order = np.argsort(-np.nan_to_num(variances, nan=-np.inf))
    highly_variable = np.zeros(adata.n_vars, dtype=bool)
    highly_variable[order[:n_top]] = True

    adata.var["highly_variable"] = highly_variable
    adata.var["means"] = matrix_mean_vector(adata.X, axis=0)
    adata.var["dispersions"] = variances
    adata.var["dispersions_norm"] = variances
    adata.uns["hvg"] = {
        "flavor": "variance",
        "n_top_genes": n_top,
        "reason": "matrix values are already transformed or Seurat HVG failed",
    }


def assign_comparison_metric_columns(adata, dataset_name, expression_profile):
    sources = {}

    total_column, total_values = find_numeric_column(
        adata.obs, CELL_TOTAL_COUNT_COLUMNS, require_nonnegative=True
    )
    if total_values is not None:
        adata.obs["comparison_total_counts"] = total_values.to_numpy(dtype=float)
        sources["cell_total_counts"] = f"obs['{total_column}']"
    elif expression_profile["value_kind"] != "signed_transformed":
        adata.obs["comparison_total_counts"] = matrix_sum(adata.X, axis=1)
        sources["cell_total_counts"] = "X row sums"
    else:
        adata.obs["comparison_total_counts"] = np.nan
        sources["cell_total_counts"] = "unavailable"

    n_gene_column, n_gene_values = find_numeric_column(
        adata.obs,
        CELL_DETECTED_GENE_COLUMNS,
        require_nonnegative=True,
        required_name_terms=["gene", "feature", "detected"],
    )
    if n_gene_values is not None:
        adata.obs["comparison_n_genes_by_counts"] = n_gene_values.to_numpy(dtype=float)
        sources["cell_detected_genes"] = f"obs['{n_gene_column}']"
    elif expression_profile["value_kind"] != "signed_transformed":
        adata.obs["comparison_n_genes_by_counts"] = matrix_nnz(adata.X, axis=1)
        sources["cell_detected_genes"] = "X > 0 row counts"
    else:
        adata.obs["comparison_n_genes_by_counts"] = np.nan
        sources["cell_detected_genes"] = "unavailable for signed transformed X"

    mean_column, mean_values = find_numeric_column(
        adata.var, GENE_MEAN_COLUMNS, require_nonnegative=False
    )
    if mean_values is not None:
        adata.var["comparison_mean_expression"] = mean_values.to_numpy(dtype=float)
        sources["gene_mean_expression"] = f"var['{mean_column}']"
    else:
        adata.var["comparison_mean_expression"] = matrix_mean_vector(adata.X, axis=0)
        sources["gene_mean_expression"] = "X column means"

    n_cell_column, n_cell_values = find_numeric_column(
        adata.var,
        GENE_DETECTED_CELL_COLUMNS,
        require_nonnegative=True,
        required_name_terms=["cell", "detected", "expressed"],
    )
    if n_cell_values is not None:
        detected_cells = n_cell_values.to_numpy(dtype=float)
        adata.var["comparison_dropout_pct"] = 100 - (detected_cells / adata.n_obs * 100)
        sources["gene_dropout_pct"] = f"var['{n_cell_column}']"
    elif expression_profile["value_kind"] != "signed_transformed":
        detected_cells = matrix_nnz(adata.X, axis=0)
        adata.var["comparison_dropout_pct"] = 100 - (detected_cells / adata.n_obs * 100)
        sources["gene_dropout_pct"] = "X > 0 column counts"
    else:
        adata.var["comparison_dropout_pct"] = np.nan
        sources["gene_dropout_pct"] = "unavailable for signed transformed X"

    adata.uns["comparison_metric_sources"] = sources
    log_record("comparison_metric_source", dataset=dataset_name, sources=sources)
    return sources


def detect_perturbation_column(obs):
    normalized_to_column = {
        normalize_column_name(column): column for column in obs.columns
    }
    for name in PERTURBATION_COLUMN_CANDIDATES:
        column = (
            name
            if name in obs.columns
            else normalized_to_column.get(normalize_column_name(name))
        )
        if column is not None:
            return column

    best_column = None
    best_score = -np.inf
    for column in obs.columns:
        series = obs[column]
        if pd.api.types.is_numeric_dtype(series):
            continue

        values = series.astype(str)
        sample = values[values.str.len() > 0].head(1000)
        if sample.empty:
            continue

        name = normalize_column_name(column)
        score = 0
        if any(term in name for term in ["perturb", "guide", "grna", "sgrna", "sgid"]):
            score += 8
        if any(term in name for term in ["probe", "target"]):
            score += 4
        if "geneid" in name or name == "gene":
            score -= 5

        sample_lower = sample.str.lower()
        if sample.str.contains(r"\|", regex=True).mean() > 0.05:
            score += 4
        if sample.str.contains("_", regex=False).mean() > 0.05:
            score += 2
        if sample_lower.str.contains(
            "non-targeting|control|scramble", regex=True
        ).any():
            score += 2
        if sample.str.startswith("ENSG").mean() > 0.1:
            score -= 6

        if score > best_score:
            best_score = score
            best_column = column

    return best_column if best_score > 0 else None


def filter_low_signal_cells_and_genes(adata):
    """Apply liberal expression QC before comparison and downstream export."""
    cells_before = adata.n_obs
    genes_before = adata.n_vars

    total_counts = matrix_sum(adata.X, axis=1)
    n_genes_by_counts = matrix_nnz(adata.X, axis=1)
    min_total_counts = lower_quantile_threshold(
        total_counts, CELL_MIN_COUNTS_FLOOR, CELL_QC_LOWER_QUANTILE
    )
    min_genes_by_counts = lower_quantile_threshold(
        n_genes_by_counts, CELL_MIN_GENES_FLOOR, CELL_QC_LOWER_QUANTILE
    )

    adata.obs["qc_total_counts"] = total_counts
    adata.obs["qc_n_genes_by_counts"] = n_genes_by_counts

    cell_mask = (total_counts >= min_total_counts) & (
        n_genes_by_counts >= min_genes_by_counts
    )
    if not np.any(cell_mask):
        raise ValueError(
            "Expression QC removed all cells. "
            f"Thresholds were total_counts >= {min_total_counts}, "
            f"n_genes_by_counts >= {min_genes_by_counts}."
        )

    adata = adata[cell_mask].copy()

    n_cells_by_counts = matrix_nnz(adata.X, axis=0)
    min_cells_by_counts = int(
        max(GENE_MIN_CELLS_FLOOR, np.ceil(GENE_MIN_CELLS_PCT * adata.n_obs))
    )
    gene_mask = n_cells_by_counts >= min_cells_by_counts
    if not np.any(gene_mask):
        raise ValueError(
            "Expression QC removed all genes. "
            f"Threshold was detection in >= {min_cells_by_counts} cells."
        )

    adata.var["qc_n_cells_by_counts"] = n_cells_by_counts
    adata.var["qc_pct_cells_by_counts"] = (
        n_cells_by_counts / adata.n_obs * 100 if adata.n_obs else 0.0
    )
    adata = adata[:, gene_mask].copy()

    summary = {
        "cells_before": int(cells_before),
        "cells_after": int(adata.n_obs),
        "cells_removed": int(cells_before - adata.n_obs),
        "genes_before": int(genes_before),
        "genes_after": int(adata.n_vars),
        "genes_removed": int(genes_before - adata.n_vars),
        "min_total_counts": int(min_total_counts),
        "min_genes_by_counts": int(min_genes_by_counts),
        "min_cells_by_counts": int(min_cells_by_counts),
        "min_cells_pct": float(min_cells_by_counts / adata.n_obs * 100),
        "cell_qc_lower_quantile": float(CELL_QC_LOWER_QUANTILE),
    }
    adata.uns["expression_qc"] = summary
    log_record("qc_filter", **summary)
    return adata, summary


def parse_gtf_attributes(raw):
    fields = {}
    for item in raw.rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        fields[key] = value.strip().strip('"')
    return fields


def load_gtf_gene_symbols(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Reference GTF not found: {path}")

    opener = gzip.open if str(path).endswith(".gz") else open
    mapping = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            if gene_id and gene_name:
                mapping[gene_id] = gene_name
                mapping.setdefault(gene_id.split(".", 1)[0], gene_name)
    if not mapping:
        raise ValueError(
            f"No gene_id/gene_name records loaded from reference GTF: {path}"
        )
    return mapping


def make_unique_index(values):
    seen = {}
    unique = []
    for value in values:
        base = str(value).strip()
        if not base:
            raise ValueError(
                "Cannot create feature index: empty gene symbol encountered"
            )
        count = seen.get(base, 0)
        candidate = base if count == 0 else f"{base}-{count}"
        while candidate in seen:
            count += 1
            candidate = f"{base}-{count}"
        seen[base] = count + 1
        seen[candidate] = 1
        unique.append(candidate)
    return pd.Index(unique)


def annotate_expression_gene_symbols(adata, dataset_name):
    """Ensure expression features carry gene symbols and use them as var_names."""
    original_index = pd.Index(adata.var_names.astype(str))
    if "gene_id" not in adata.var.columns:
        adata.var["gene_id"] = original_index.to_numpy()

    if "gene_name" in adata.var.columns:
        symbols = adata.var["gene_name"].astype("string").str.strip()
        source = "var['gene_name']"
        missing = symbols.isna() | (symbols.str.len() == 0)
        if missing.any():
            examples = original_index[missing.to_numpy()][:10].tolist()
            raise ValueError(
                f"{dataset_name} has {int(missing.sum())} empty var['gene_name'] "
                f"values. Examples: {examples}"
            )
    else:
        mapping = load_gtf_gene_symbols(REFERENCE_GTF_PATH)
        symbols_by_row = []
        unresolved = []
        for feature_name, gene_id in zip(
            original_index, adata.var["gene_id"].astype(str)
        ):
            gene_id = str(gene_id).strip()
            feature_name = str(feature_name).strip()
            symbol = mapping.get(gene_id) or mapping.get(gene_id.split(".", 1)[0])
            if symbol is None:
                unresolved.append({"feature_name": feature_name, "gene_id": gene_id})
                symbols_by_row.append(None)
            else:
                symbols_by_row.append(symbol)

        if unresolved:
            keep_mask = np.asarray([symbol is not None for symbol in symbols_by_row])
            log_record(
                "gene_symbol_filter",
                dataset=dataset_name,
                n_removed=int((~keep_mask).sum()),
                n_before=int(len(keep_mask)),
            )
            adata._inplace_subset_var(keep_mask)
            original_index = pd.Index(adata.var_names.astype(str))
            symbols_by_row = [symbol for symbol in symbols_by_row if symbol is not None]
            if not symbols_by_row:
                raise ValueError(
                    f"No {dataset_name} expression features could be resolved to "
                    f"gene symbols using {REFERENCE_GTF_PATH}."
                )

        symbols = pd.Series(symbols_by_row, index=adata.var.index, dtype="string")
        source = REFERENCE_GTF_PATH

    adata.var["gene_symbol"] = symbols.to_numpy()
    adata.var["feature_id"] = original_index.to_numpy()
    adata.var_names = make_unique_index(adata.var["gene_symbol"])
    adata.var.index.name = "gene_symbol_unique"

    n_mapped = int(
        sum(
            bool(symbol)
            and not str(symbol).startswith("ENSG")
            and str(symbol) != str(gene_id)
            for symbol, gene_id in zip(adata.var["gene_symbol"], adata.var["gene_id"])
        )
    )
    log_record(
        "gene_symbol_annotation",
        dataset=dataset_name,
        n_mapped=n_mapped,
        n_genes=adata.n_vars,
        source=source,
    )
    return adata


def preprocess_adata(adata, name, target_sum=1e4, n_top_genes=2000):
    """Standardized preprocessing with automatic count/transformed matrix handling."""
    expression_profile = profile_expression_matrix(adata.X)
    is_raw = expression_profile["raw_counts_detected"]
    normalization_action = (
        "normalize_total + log1p" if is_raw else "using transformed X as provided"
    )
    log_record(
        "preprocess_input",
        dataset=name,
        cells=adata.n_obs,
        genes=adata.n_vars,
        raw_counts_detected=bool(is_raw),
        expression_value_kind=expression_profile["value_kind"],
        normalization_action=normalization_action,
        normalization_target_sum=target_sum,
        hvg_n_top_genes=n_top_genes,
    )

    adata.uns["comparison_expression_profile"] = expression_profile
    if is_raw:
        adata.layers["counts"] = adata.X.copy()
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        try:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=min(n_top_genes, adata.n_vars),
                flavor="seurat",
                subset=False,
            )
        except ValueError as exc:
            print(
                f"HVG selection - {name}: falling back to variance ranking after "
                f"Seurat HVG failed ({exc})"
            )
            set_highly_variable_by_variance(adata, n_top_genes)
    else:
        adata.layers["transformed_expression"] = adata.X.copy()
        set_highly_variable_by_variance(adata, n_top_genes)

    assign_comparison_metric_columns(adata, name, expression_profile)
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

    plot_df = df.copy().replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])

    if log_scale:
        plot_df = plot_df[(plot_df[x_col] > -1) & (plot_df[y_col] > -1)]

    if len(plot_df) < 2:
        plt.close(fig)
        return {
            "pearson": np.nan,
            "spearman": np.nan,
            "pct_deviant": np.nan,
            "skipped": True,
            "reason": "fewer than two finite comparable values",
        }

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
    if plot_df[x_col].nunique() < 2 or plot_df[y_col].nunique() < 2:
        pearson = np.nan
        spearman = np.nan
    else:
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
    if guide_name.startswith(CONTROL_TARGET_SYMBOL):
        return CONTROL_TARGET_SYMBOL
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
            if gene and gene != CONTROL_TARGET_SYMBOL:
                genes.add(gene)
    return sorted(genes)


def control_probes_from_probe_list(probes):
    control_probes = []
    for probe in probes:
        if CONTROL_TARGET_SYMBOL in guide_target_names(probe):
            control_probes.append(probe)
    return sorted(control_probes)


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


def perturbation_call_type(n_genes, n_control_probes, n_probes):
    if n_genes == 0 and n_control_probes > 0:
        return "control"
    if n_genes == 1 and n_control_probes == 0:
        return "single_gene"
    if n_genes > 1 and n_control_probes == 0:
        return "multi_gene"
    if n_genes > 0 and n_control_probes > 0:
        return "mixed_gene_control"
    if n_probes == 0:
        return "unassigned"
    return "unclassified_probe"


def build_gene_call_table(labels):
    records = []
    for cell_id, label in labels.items():
        probe_label = canonical_label_for_display(label)
        probes = label_to_probe_list(probe_label)
        genes = genes_from_probe_list(probes)
        control_probes = control_probes_from_probe_list(probes)
        control_label = "|".join(control_probes) if control_probes else "None"
        records.append(
            {
                "cell_id": cell_id,
                "probe_label": probe_label,
                "n_probes": len(probes),
                "control_probe_label": control_label,
                "n_control_probes": len(control_probes),
                "genes": genes,
                "gene_label": "|".join(genes) if genes else "None",
                "n_genes": len(genes),
                "outcome": gene_call_outcome(probes, genes),
                "perturbation_call_type": perturbation_call_type(
                    len(genes), len(control_probes), len(probes)
                ),
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


def outcome_matrix_dataframe(matrix_dict):
    return (
        pd.DataFrame.from_dict(matrix_dict, orient="index")
        .reindex(index=GENE_CALL_OUTCOMES, columns=GENE_CALL_OUTCOMES, fill_value=0)
        .astype(int)
    )


def plot_gene_outcome_matrix(
    model_name, matrix_dict, single_gene_match, count_threshold
):
    matrix_df = outcome_matrix_dataframe(matrix_dict)
    display_labels = ["0 genes", "1 gene\n1 probe", "1 gene\n2+ probes", ">1 gene"]
    plot_path = f"comparison_results/{DATASET_ID}/perturbation_gene_outcome_matrix_{model_name}.png"

    fig, ax = plt.subplots(figsize=(8.5, 7))
    sns.heatmap(
        matrix_df,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        cbar_kws={"label": "Cells"},
        xticklabels=display_labels,
        yticklabels=display_labels,
        ax=ax,
    )

    same = single_gene_match["n_same_gene"]
    total = single_gene_match["n_cells_both_single_gene"]
    same_pct = (same / total * 100) if total else 0.0
    ax.set_title(
        f"{model_name.replace('_', ' ')}: Gene-Call Outcome Matrix\n"
        f"threshold >= {count_threshold} UMI | same single-gene calls: "
        f"{same}/{total} ({same_pct:.1f}%)",
        fontsize=13,
        pad=12,
    )
    ax.set_xlabel("Reprocessed outcome")
    ax.set_ylabel("Curated outcome")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)
    return plot_path


def single_gene_match_summary(cur_calls, rep_calls):
    cur_single = (cur_calls["n_genes"] == 1) & (cur_calls["n_control_probes"] == 0)
    rep_single = (rep_calls["n_genes"] == 1) & (rep_calls["n_control_probes"] == 0)
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


def call_probe_lists(adata, return_diagnostics=False):
    """Call all probe features above the Gaussian-Poisson-derived threshold."""
    if "guides" not in adata.obsm:
        return (None, None) if return_diagnostics else None

    guide_matrix = adata.obsm["guides"]
    if sp.issparse(guide_matrix):
        guide_matrix = guide_matrix.tocsr()

    if sp.issparse(guide_matrix):
        non_zeros = guide_matrix.data
    else:
        non_zeros = guide_matrix[guide_matrix > 0]

    if len(non_zeros) > 0:
        mixture_res = fit_poisson_gaussian_mixture(non_zeros)
        count_threshold = mixture_res["decision_threshold"]
    else:
        mixture_res = None
        count_threshold = 5

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
        "method": "gaussian_poisson",
        "count_threshold": int(count_threshold),
        "n_cells": int(guide_matrix.shape[0]),
        "n_guides": int(guide_matrix.shape[1]),
        "n_cells_with_any_called_probe": int((positive_guide_counts >= 1).sum()),
        "n_cells_with_two_or_more_called_probes": int(
            (positive_guide_counts >= 2).sum()
        ),
        "n_cells_with_any_called_gene": int((gene_calls["n_genes"] >= 1).sum()),
        "n_cells_with_any_called_control_probe": int(
            (gene_calls["n_control_probes"] >= 1).sum()
        ),
        "n_control_only_calls": int(
            (gene_calls["perturbation_call_type"] == "control").sum()
        ),
        "n_single_gene_no_control_calls": int(
            (gene_calls["perturbation_call_type"] == "single_gene").sum()
        ),
        "called_probe_count_distribution": {
            str(int(count)): int(freq) for count, freq in zip(count_values, count_freqs)
        },
        "gene_call_outcome_counts": outcome_count_dict(gene_calls["outcome"]),
        "perturbation_call_type_counts": {
            str(call_type): int(count)
            for call_type, count in gene_calls["perturbation_call_type"]
            .value_counts()
            .sort_index()
            .items()
        },
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


def annotate_knockout_genes(adata):
    """Add per-cell Gaussian-Poisson probe and knockout-gene calls to adata.obs."""
    probe_labels, diagnostics = call_probe_lists(adata, return_diagnostics=True)
    if probe_labels is None:
        log_record(
            "perturbation_comparison_skipped",
            reason="reprocessed_guide_matrix_missing",
            expected_obsm_key="guides",
        )
        return None

    call_table = build_gene_call_table(probe_labels)
    gene_labels = call_table["gene_label"].replace("None", "")
    single_gene_no_control = call_table["perturbation_call_type"] == "single_gene"
    control_only = call_table["perturbation_call_type"] == "control"

    adata.obs["called_probe_label"] = call_table["probe_label"].astype(str).to_numpy()
    adata.obs["called_probe_count"] = call_table["n_probes"].astype(int).to_numpy()
    adata.obs["called_control_probe_label"] = (
        call_table["control_probe_label"].replace("None", "").astype(str).to_numpy()
    )
    adata.obs["called_control_probe_count"] = (
        call_table["n_control_probes"].astype(int).to_numpy()
    )
    adata.obs["has_called_control_probe"] = (
        call_table["n_control_probes"].astype(int).to_numpy() > 0
    )
    adata.obs["called_knockout_genes"] = gene_labels.astype(str).to_numpy()
    adata.obs["called_knockout_gene_count"] = (
        call_table["n_genes"].astype(int).to_numpy()
    )
    adata.obs["knockout_call_outcome"] = call_table["outcome"].astype(str).to_numpy()
    adata.obs["perturbation_call_type"] = (
        call_table["perturbation_call_type"].astype(str).to_numpy()
    )
    adata.obs["is_control"] = control_only.to_numpy()
    adata.obs["is_single_gene_perturbation"] = single_gene_no_control.to_numpy()
    adata.obs["perturbed_target_symbol"] = np.where(
        single_gene_no_control.to_numpy(), call_table["gene_label"].to_numpy(), ""
    )
    adata.obs["perturbation_call_method"] = "gaussian_poisson"
    adata.obs["dataset_id"] = DATASET_ID

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    adata.uns["dataset_id"] = DATASET_ID
    adata.uns["perturbation_calling"] = diagnostics
    adata.uns["downstream_columns"] = {
        "expression_counts": "layers['counts']",
        "perturbed_target_symbol": "obs['perturbed_target_symbol']",
        "control_indicator": "obs['is_control']",
        "perturbation_indicator": "obs['is_single_gene_perturbation']",
        "control_probe_count": "obs['called_control_probe_count']",
        "all_called_knockout_genes": "obs['called_knockout_genes']",
        "call_outcome": "obs['knockout_call_outcome']",
        "perturbation_call_type": "obs['perturbation_call_type']",
    }

    summary = {
        "n_cells": int(adata.n_obs),
        "count_threshold": int(diagnostics["count_threshold"]),
        "n_single_gene_calls": int(single_gene_no_control.sum()),
        "n_single_gene_with_control_calls": int(
            ((call_table["n_genes"] == 1) & (call_table["n_control_probes"] > 0)).sum()
        ),
        "n_control_calls": int(control_only.sum()),
        "n_multi_gene_calls": int((call_table["n_genes"] > 1).sum()),
        "n_no_gene_calls": int((call_table["n_genes"] == 0).sum()),
        "n_no_guide_calls": int((call_table["n_probes"] == 0).sum()),
        "n_gene_control_mixed_calls": int(
            (call_table["perturbation_call_type"] == "mixed_gene_control").sum()
        ),
        "outcome_counts": outcome_count_dict(call_table["outcome"]),
        "perturbation_call_type_counts": {
            str(call_type): int(count)
            for call_type, count in call_table["perturbation_call_type"]
            .value_counts()
            .sort_index()
            .items()
        },
    }
    log_record("knockout_annotation", **summary)
    return summary


def write_filtered_h5ad(adata, path):
    output_dir = os.path.dirname(path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")
    log_record(
        "filtered_h5ad_written",
        path=path,
        cells=adata.n_obs,
        genes=adata.n_vars,
        compression="gzip",
    )


def compare_perturbations(adata_cur, adata_rep, common_cells):
    """Compare perturbation assignments at gene level, after probe calling."""
    cur_pert_col = detect_perturbation_column(adata_cur.obs)

    if not cur_pert_col:
        log_record(
            "perturbation_comparison_skipped",
            reason="curated_perturbation_column_missing",
            candidate_columns=PERTURBATION_COLUMN_CANDIDATES,
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
        "n_cells_with_any_control_probe": int(
            (cur_calls["n_control_probes"] > 0).sum()
        ),
        "n_valid_single_gene_no_control_cells": int(
            ((cur_calls["n_genes"] == 1) & (cur_calls["n_control_probes"] == 0)).sum()
        ),
        "outcome_counts": outcome_count_dict(cur_calls["outcome"]),
        "perturbation_call_type_counts": {
            str(call_type): int(count)
            for call_type, count in cur_calls["perturbation_call_type"]
            .value_counts()
            .sort_index()
            .items()
        },
        "probe_count_distribution": {
            str(int(count)): int(freq)
            for count, freq in cur_calls["n_probes"].value_counts().sort_index().items()
        },
    }
    log_record("curated_gene_call_summary", **curated_summary)

    model_specs = [{"model": "gaussian_poisson"}]

    model_results = {}
    for spec in model_specs:
        if "called_probe_label" in adata_rep.obs:
            rep_diagnostics = adata_rep.uns.get("perturbation_calling", {})
            p_rep = (
                adata_rep.obs.loc[common_cells, "called_probe_label"]
                .astype(str)
                .apply(canonical_label_for_display)
            )
        else:
            rep_labels, rep_diagnostics = call_probe_lists(
                adata_rep,
                return_diagnostics=True,
            )
            p_rep = (
                rep_labels.loc[common_cells]
                .astype(str)
                .apply(canonical_label_for_display)
            )
        rep_calls = build_gene_call_table(p_rep)
        reprocessed_common_call_summary = {
            "n_common_cells": int(len(common_cells)),
            "n_cells_with_any_probe": int((rep_calls["n_probes"] > 0).sum()),
            "n_cells_with_any_gene": int((rep_calls["n_genes"] > 0).sum()),
            "n_cells_with_any_control_probe": int(
                (rep_calls["n_control_probes"] > 0).sum()
            ),
            "n_valid_single_gene_no_control_cells": int(
                (
                    (rep_calls["n_genes"] == 1) & (rep_calls["n_control_probes"] == 0)
                ).sum()
            ),
        }
        outcome_matrix = outcome_matrix_dict(cur_calls["outcome"], rep_calls["outcome"])
        single_gene_match = single_gene_match_summary(cur_calls, rep_calls)
        plot_path = plot_gene_outcome_matrix(
            spec["model"],
            outcome_matrix,
            single_gene_match,
            rep_diagnostics.get("count_threshold"),
        )

        result = {
            "method": "gaussian_poisson",
            "count_threshold": rep_diagnostics.get("count_threshold"),
            "probe_call_diagnostics": rep_diagnostics,
            "reprocessed_common_call_summary": reprocessed_common_call_summary,
            "reprocessed_outcome_counts": outcome_count_dict(rep_calls["outcome"]),
            "outcome_matrix_rows_curated_columns_reprocessed": outcome_matrix,
            "single_gene_match": single_gene_match,
            "plot_path": plot_path,
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
        "control_definition": (
            "control means one or more non-targeting probes called and zero "
            "gene-targeting probes called"
        ),
        "valid_single_gene_definition": (
            "exactly one gene-targeting gene called and zero non-targeting probes called"
        ),
    }


# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================
os.makedirs(f"comparison_results/{DATASET_ID}", exist_ok=True)
log_record(
    "input_paths",
    curated_h5ad_path=CURATED_H5AD_PATH,
    reprocessed_h5ad_path=REPROCESSED_H5AD_PATH,
    filtered_reprocessed_h5ad_path=FILTERED_REPROCESSED_H5AD_PATH,
)
adata_cur = filter_unique_barcodes(sc.read_h5ad(CURATED_H5AD_PATH), "Curated")
adata_rep = filter_unique_barcodes(sc.read_h5ad(REPROCESSED_H5AD_PATH), "Reprocessed")

adata_cur = annotate_expression_gene_symbols(adata_cur, "Curated")
adata_rep = annotate_expression_gene_symbols(adata_rep, "Reprocessed")

adata_rep, qc_summary = filter_low_signal_cells_and_genes(adata_rep)
knockout_annotation_summary = annotate_knockout_genes(adata_rep)
write_filtered_h5ad(adata_rep, FILTERED_REPROCESSED_H5AD_PATH)

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
        "filtered_reprocessed_h5ad_path": FILTERED_REPROCESSED_H5AD_PATH,
    },
    "expression_qc": qc_summary,
    "knockout_annotation": knockout_annotation_summary,
    "overlap": overlap_summary,
    "preprocessing": {
        "curated": cur_sub.uns.get("comparison_expression_profile", {}),
        "reprocessed": rep_sub.uns.get("comparison_expression_profile", {}),
    },
    "metric_sources": {
        "curated": cur_sub.uns.get("comparison_metric_sources", {}),
        "reprocessed": rep_sub.uns.get("comparison_metric_sources", {}),
    },
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
plt.savefig(f"comparison_results/{DATASET_ID}/overlap_distributions.png")
plt.show()
plt.close()

# --- Cell-wise Scatter Plots ---
results_summary["cell_metrics"]["total_counts"] = plot_scatter_comparison(
    pd.DataFrame(
        {
            "cur": cur_sub.obs["comparison_total_counts"],
            "rep": rep_sub.obs["comparison_total_counts"],
        }
    ),
    "cur",
    "rep",
    "Total UMI Counts",
    "Original",
    "Reprocessed",
    f"comparison_results/{DATASET_ID}/counts_comparison.png",
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
            "cur": cur_sub.obs["comparison_n_genes_by_counts"],
            "rep": rep_sub.obs["comparison_n_genes_by_counts"],
        }
    ),
    "cur",
    "rep",
    "Number of Detected Genes",
    "Original",
    "Reprocessed",
    f"comparison_results/{DATASET_ID}/genes_comparison.png",
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
        "mean_cur": cur_sub.var["comparison_mean_expression"],
        "mean_rep": rep_sub.var["comparison_mean_expression"],
        "dropout_cur": cur_sub.var["comparison_dropout_pct"],
        "dropout_rep": rep_sub.var["comparison_dropout_pct"],
    }
)

results_summary["gene_metrics"]["mean_expression"] = plot_scatter_comparison(
    gene_metrics_df,
    "mean_cur",
    "mean_rep",
    "Mean Gene Expression",
    "Original",
    "Reprocessed",
    f"comparison_results/{DATASET_ID}/gene_expression_mean.png",
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
    f"comparison_results/{DATASET_ID}/sparsity_comparison.png",
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
with open(f"comparison_results/{DATASET_ID}/summary_report.txt", "w") as f:
    f.write(json.dumps(results_summary, indent=4, default=json_default))
log_record(
    "summary_report_written", path=f"comparison_results/{DATASET_ID}/summary_report.txt"
)

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

os.makedirs(f"comparison_results/{DATASET_ID}/supplementary", exist_ok=True)

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
    f"comparison_results/{DATASET_ID}/supplementary/cellwise_correlation_dist.png",
    bbox_inches="tight",
    dpi=300,
)
plt.show()
