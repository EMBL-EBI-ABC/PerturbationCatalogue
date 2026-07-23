import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba")

import gzip
import json
import re
import textwrap
from collections import Counter, defaultdict
import sys

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from scipy import stats

try:
    from IPython.display import display
except ImportError:
    display = print


DATASET_ID = sys.argv[1]
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
REFERENCE_GTF_PATH = (
    "/hps/nobackup/mfreeberg/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz"
)

ROW_CHUNK_SIZE = int(os.environ.get("PERTURBSEQ_ROW_CHUNK_SIZE", "2048"))
SCATTER_MAX_POINTS = int(os.environ.get("PERTURBSEQ_SCATTER_MAX_POINTS", "200000"))
CELL_CORR_SAMPLE_SIZE = int(os.environ.get("PERTURBSEQ_CELL_CORR_SAMPLE_SIZE", "10000"))
RANDOM_SEED = int(os.environ.get("PERTURBSEQ_RANDOM_SEED", "1"))

CELL_QC_LOWER_QUANTILE = 0.01
CELL_MIN_COUNTS_FLOOR = 1000
CELL_MIN_GENES_FLOOR = 200
GENE_MIN_CELLS_FLOOR = 10
GENE_MIN_CELLS_PCT = 0.01
TARGET_SUM = 1e4
CONTROL_TARGET_SYMBOL = "non-targeting"
GENE_CALL_OUTCOMES = ["0_genes", "1_gene_1_probe", "1_gene_2_probes", ">1_gene"]

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


def format_count_pct(count, total):
    pct = (count / total * 100) if total else 0.0
    return f"{count}/{total} ({pct:.1f}%)"


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
    elif event == "gene_symbol_filter":
        print(
            f"Gene symbol filtering - {fields['dataset']}: "
            f"removed {format_count_pct(fields['n_removed'], fields['n_before'])} "
            "features without gene_name in the reference GTF"
        )
    elif event == "gene_symbol_annotation":
        print(
            f"Gene symbol annotation - {fields['dataset']}: "
            f"mapped {format_count_pct(fields['n_mapped'], fields['n_genes'])}; "
            f"source {fields['source']}"
        )
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
    elif event == "filtered_h5ad_written":
        print(
            f"Filtered H5AD written: {fields['path']} "
            f"({fields['cells']} cells x {fields['genes']} genes; compression={fields['compression']})"
        )
    elif event == "overlap_summary":
        print("Overlap")
        print(
            f"  Cells: {fields['common_cells']} common; curated {fields['curated_cells']}, "
            f"reprocessed {fields['reprocessed_cells']}; "
            f"{fields['common_cells_pct_of_curated']:.1f}% of curated"
        )
        print(
            f"  Genes: {fields['common_genes']} common; curated {fields['curated_genes']}, "
            f"reprocessed {fields['reprocessed_genes']}; "
            f"{fields['common_genes_pct_of_curated']:.1f}% of curated"
        )
    elif event == "preprocess_input":
        raw_status = "yes" if fields["raw_counts_detected"] else "no"
        print(
            f"Preprocessing input - {fields['dataset']}: {fields['cells']} cells, "
            f"{fields['genes']} genes, raw counts detected: {raw_status}"
        )
        print(
            f"  Expression values: {fields['expression_value_kind']}; "
            f"normalization: {fields['normalization_action']}"
        )
    elif event == "comparison_metric_source":
        print(f"Metric sources - {fields['dataset']}")
        for metric, source in fields["sources"].items():
            print(f"  {metric}: {source}")
    elif event == "comparison_metric":
        if fields.get("skipped"):
            print(
                f"{fields['metric_group']}.{fields['metric']}: skipped ({fields['reason']})"
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
    elif event == "perturbation_gene_comparison":
        match = fields["single_gene_match"]
        common = fields["reprocessed_common_call_summary"]
        print(
            f"  {fields['model']}: threshold >= {fields['count_threshold']} UMI; "
            f"called probes in {format_count_pct(common['n_cells_with_any_probe'], common['n_common_cells'])}; "
            f"called genes in {format_count_pct(common['n_cells_with_any_gene'], common['n_common_cells'])}"
        )
        print(
            f"    Same single-gene calls: "
            f"{format_count_pct(match['n_same_gene'], match['n_cells_both_single_gene'])}; "
            f"different {format_count_pct(match['n_different_gene'], match['n_cells_both_single_gene'])}"
        )
    elif event == "summary_report_written":
        print(f"Summary report: {fields['path']}")
    else:
        print(f"{event}: {fields}")


def normalize_column_name(value):
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def chunk_array(values, chunk_size=ROW_CHUNK_SIZE):
    values = np.asarray(values)
    for start in range(0, len(values), chunk_size):
        yield values[start : start + chunk_size]


def as_csr(block):
    return block.tocsr() if sp.issparse(block) else sp.csr_matrix(np.asarray(block))


def read_block(matrix, row_pos, col_pos=None):
    row_pos = np.asarray(row_pos, dtype=np.int64)
    if len(row_pos) == 0:
        n_cols = len(col_pos) if col_pos is not None else matrix.shape[1]
        return sp.csr_matrix((0, n_cols))

    order = np.argsort(row_pos)
    rows = row_pos[order]
    block = matrix[rows, :]
    if col_pos is not None:
        block = block[:, np.asarray(col_pos, dtype=np.int64)]
    if not np.all(order == np.arange(len(order))):
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        block = block[inv, :]
    return block


class H5CSRMatrix:
    """Minimal row-slicing wrapper for an H5AD CSR group."""

    def __init__(self, group):
        self.group = group
        self.data = group["data"]
        self.indices = group["indices"]
        self.indptr = group["indptr"]
        self.shape = tuple(int(x) for x in group.attrs["shape"])
        self.dtype = self.data.dtype

    def __getitem__(self, key):
        rows, cols = key
        if isinstance(rows, slice):
            row_idx = np.arange(*rows.indices(self.shape[0]), dtype=np.int64)
        else:
            row_idx = np.asarray(rows, dtype=np.int64)

        if isinstance(cols, slice) and cols == slice(None):
            col_idx = None
            n_cols = self.shape[1]
        else:
            col_idx = np.asarray(cols, dtype=np.int64)
            n_cols = len(col_idx)
            remap = {int(old): new for new, old in enumerate(col_idx)}

        data_parts = []
        index_parts = []
        indptr = np.zeros(len(row_idx) + 1, dtype=np.int64)
        nnz = 0
        for out_i, row in enumerate(row_idx):
            start = int(self.indptr[row])
            end = int(self.indptr[row + 1])
            row_data = self.data[start:end]
            row_indices = self.indices[start:end]
            if col_idx is not None:
                keep = np.isin(row_indices, col_idx)
                row_data = row_data[keep]
                row_indices = np.fromiter(
                    (remap[int(i)] for i in row_indices[keep]),
                    dtype=np.int64,
                    count=int(keep.sum()),
                )
            data_parts.append(np.asarray(row_data))
            index_parts.append(np.asarray(row_indices, dtype=np.int64))
            nnz += len(row_data)
            indptr[out_i + 1] = nnz

        data = (
            np.concatenate(data_parts).astype(self.dtype, copy=False)
            if nnz
            else np.array([], dtype=self.dtype)
        )
        indices = (
            np.concatenate(index_parts).astype(np.int64, copy=False)
            if nnz
            else np.array([], dtype=np.int64)
        )
        return sp.csr_matrix((data, indices, indptr), shape=(len(row_idx), n_cols))


def matrix_row_sum_nnz(matrix, row_pos, col_pos):
    sums = np.zeros(len(row_pos), dtype=np.float64)
    nnz = np.zeros(len(row_pos), dtype=np.int64)
    offset = 0
    for rows in chunk_array(row_pos):
        block = read_block(matrix, rows, col_pos)
        sums[offset : offset + len(rows)] = np.asarray(block.sum(axis=1)).ravel()
        if sp.issparse(block):
            nnz[offset : offset + len(rows)] = block.getnnz(axis=1)
        else:
            nnz[offset : offset + len(rows)] = np.asarray(block > 0).sum(axis=1)
        offset += len(rows)
    return sums, nnz


def matrix_col_sum_nnz(matrix, row_pos, col_pos):
    sums = np.zeros(len(col_pos), dtype=np.float64)
    nnz = np.zeros(len(col_pos), dtype=np.int64)
    for rows in chunk_array(row_pos):
        block = read_block(matrix, rows, col_pos)
        sums += np.asarray(block.sum(axis=0)).ravel()
        if sp.issparse(block):
            nnz += block.getnnz(axis=0)
        else:
            nnz += np.asarray(block > 0).sum(axis=0)
    return sums, nnz


def sample_matrix_profile(matrix, row_pos, col_pos, max_rows=64, max_cols=512):
    rows = np.asarray(
        row_pos[
            np.linspace(
                0, len(row_pos) - 1, min(len(row_pos), max_rows), dtype=np.int64
            )
        ],
        dtype=np.int64,
    )
    cols = np.asarray(
        col_pos[
            np.linspace(
                0, len(col_pos) - 1, min(len(col_pos), max_cols), dtype=np.int64
            )
        ],
        dtype=np.int64,
    )
    block = read_block(matrix, rows, cols)
    values = block.data if sp.issparse(block) else np.asarray(block).ravel()
    values = np.asarray(values, dtype=float)
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
    value_kind = (
        "signed_transformed"
        if negative_fraction > 0
        else "raw_counts" if integer_fraction > 0.999 else "nonnegative_transformed"
    )
    return {
        "value_kind": value_kind,
        "raw_counts_detected": value_kind == "raw_counts",
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "integer_fraction": integer_fraction,
        "negative_fraction": negative_fraction,
        "below_minus_one_fraction": float(np.mean(values < -1)),
    }


def lower_quantile_threshold(values, floor, quantile):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return int(floor)
    return int(max(floor, np.floor(np.quantile(values, quantile))))


def parse_gtf_attributes(raw):
    fields = {}
    for item in raw.rstrip(";").split(";"):
        item = item.strip()
        if item and " " in item:
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
        raise ValueError(f"No gene_id/gene_name records loaded from {path}")
    return mapping


def make_unique_index(values):
    seen = {}
    result = []
    for value in values:
        base = str(value).strip()
        if not base:
            raise ValueError("Cannot create feature index from an empty gene symbol")
        i = seen.get(base, 0)
        candidate = base if i == 0 else f"{base}-{i}"
        while candidate in seen:
            i += 1
            candidate = f"{base}-{i}"
        seen[base] = i + 1
        seen[candidate] = 1
        result.append(candidate)
    return pd.Index(result)


def numeric_series_or_none(series):
    values = pd.to_numeric(series, errors="coerce")
    return None if values.notna().sum() == 0 else values.astype(float)


def find_numeric_column(df, names, require_nonnegative=True, required_name_terms=None):
    normalized = {normalize_column_name(c): c for c in df.columns}
    candidates = []
    for name in names:
        column = (
            name if name in df.columns else normalized.get(normalize_column_name(name))
        )
        if column is not None:
            candidates.append(column)
    preferred = set(candidates)
    candidates.extend(c for c in df.columns if c not in candidates)

    best = None
    best_score = -np.inf
    for column in candidates:
        values = numeric_series_or_none(df[column])
        if values is None:
            continue
        finite = values[np.isfinite(values)]
        if finite.empty or (require_nonnegative and (finite < 0).any()):
            continue
        if column in preferred:
            return column, values
        name = normalize_column_name(column)
        if required_name_terms and not any(t in name for t in required_name_terms):
            continue
        score = (
            3 * ("total" in name)
            + 3 * ("umi" in name)
            + 2 * ("count" in name)
            + 2 * any(t in name for t in ["gene", "feature", "detected"])
            - 5 * any(t in name for t in ["percent", "pct", "mito", "zscore"])
            - 5 * name.startswith("z")
        )
        if score > best_score:
            best = (column, values)
            best_score = score
    return best if best is not None and best_score > 0 else (None, None)


class DatasetView:
    def __init__(self, path, label):
        self.path = path
        self.label = label
        self.adata = ad.read_h5ad(path, backed="r")
        self.h5 = h5py.File(path, "r")
        self.obs = self.adata.obs.copy()
        self.var = self.adata.var.copy()
        self.guide_matrix = None
        if "obsm" in self.h5 and "guides" in self.h5["obsm"]:
            guide_obj = self.h5["obsm/guides"]
            if (
                isinstance(guide_obj, h5py.Group)
                and guide_obj.attrs.get("encoding-type") == "csr_matrix"
            ):
                self.guide_matrix = H5CSRMatrix(guide_obj)
        self._prepare_barcodes()
        self._prepare_gene_symbols()
        self.expression_profile = sample_matrix_profile(
            self.adata.X, self.obs_pos, self.var_pos
        )
        log_record(
            "preprocess_input",
            dataset=label,
            cells=len(self.obs_pos),
            genes=len(self.var_pos),
            raw_counts_detected=self.expression_profile["raw_counts_detected"],
            expression_value_kind=self.expression_profile["value_kind"],
            normalization_action=(
                "normalize_total + log1p"
                if self.expression_profile["raw_counts_detected"]
                else "using transformed X as provided"
            ),
        )

    def close(self):
        self.adata.file.close()
        self.h5.close()

    def _prepare_barcodes(self):
        base = pd.Index(self.adata.obs_names.astype(str)).str.split("-").str[0]
        counts = base.value_counts(sort=False)
        mask = np.asarray(base.isin(counts[counts == 1].index), dtype=bool)
        self.obs_pos = np.flatnonzero(mask)
        self.obs_names = pd.Index(base[mask], name="barcode")
        self.obs_pos_by_name = pd.Series(self.obs_pos, index=self.obs_names)
        log_record(
            "barcode_filter",
            dataset=self.label,
            original_cells=self.adata.n_obs,
            removed_collision_cells=int((~mask).sum()),
            remaining_cells=len(self.obs_pos),
        )

    def _prepare_gene_symbols(self):
        original_index = pd.Index(self.adata.var_names.astype(str))
        var = self.var.copy()
        if "gene_id" not in var.columns:
            var["gene_id"] = original_index.to_numpy()

        if "gene_name" in var.columns:
            symbols = var["gene_name"].astype("string").str.strip()
            source = "var['gene_name']"
            keep = ~(symbols.isna() | (symbols.str.len() == 0))
            if not keep.all():
                raise ValueError(
                    f"{self.label} has {int((~keep).sum())} empty var['gene_name'] values"
                )
        else:
            mapping = load_gtf_gene_symbols(REFERENCE_GTF_PATH)
            resolved = []
            for feature_name, gene_id in zip(
                original_index, var["gene_id"].astype(str)
            ):
                gene_id = gene_id.strip()
                resolved.append(
                    mapping.get(gene_id) or mapping.get(gene_id.split(".", 1)[0])
                )
            symbols = pd.Series(resolved, index=var.index, dtype="string")
            keep = symbols.notna() & (symbols.str.len() > 0)
            if not keep.all():
                log_record(
                    "gene_symbol_filter",
                    dataset=self.label,
                    n_removed=int((~keep).sum()),
                    n_before=len(keep),
                )
            source = REFERENCE_GTF_PATH

        self.var_pos = np.flatnonzero(keep.to_numpy())
        symbols = symbols.iloc[self.var_pos].astype(str).to_numpy()
        self.var_names = make_unique_index(symbols)
        self.var_pos_by_name = pd.Series(self.var_pos, index=self.var_names)
        self.var_table = pd.DataFrame(
            {
                "gene_id": var["gene_id"].iloc[self.var_pos].astype(str).to_numpy(),
                "gene_symbol": symbols,
                "feature_id": original_index[self.var_pos].to_numpy(),
            },
            index=self.var_names,
        )
        self.var_table.index.name = "gene_symbol_unique"
        n_mapped = int(
            sum(
                bool(symbol)
                and not str(symbol).startswith("ENSG")
                and str(symbol) != str(gene_id)
                for symbol, gene_id in zip(
                    self.var_table["gene_symbol"], self.var_table["gene_id"]
                )
            )
        )
        log_record(
            "gene_symbol_annotation",
            dataset=self.label,
            n_mapped=n_mapped,
            n_genes=len(self.var_pos),
            source=source,
        )


def detect_perturbation_column(obs):
    normalized = {normalize_column_name(c): c for c in obs.columns}
    for name in PERTURBATION_COLUMN_CANDIDATES:
        column = (
            name if name in obs.columns else normalized.get(normalize_column_name(name))
        )
        if column is not None:
            return column
    best_column = None
    best_score = -np.inf
    for column in obs.columns:
        if pd.api.types.is_numeric_dtype(obs[column]):
            continue
        sample = obs[column].astype(str).head(1000)
        name = normalize_column_name(column)
        score = (
            8 * any(t in name for t in ["perturb", "guide", "grna", "sgrna", "sgid"])
            + 4 * any(t in name for t in ["probe", "target"])
            + 4 * (sample.str.contains(r"\|", regex=True).mean() > 0.05)
            + 2 * (sample.str.contains("_", regex=False).mean() > 0.05)
            - 6 * (sample.str.startswith("ENSG").mean() > 0.1)
        )
        if score > best_score:
            best_column = column
            best_score = score
    return best_column if best_score > 0 else None


def filter_low_signal_cells_and_genes(rep):
    totals, n_genes = matrix_row_sum_nnz(rep.adata.X, rep.obs_pos, rep.var_pos)
    min_total = lower_quantile_threshold(
        totals, CELL_MIN_COUNTS_FLOOR, CELL_QC_LOWER_QUANTILE
    )
    min_genes = lower_quantile_threshold(
        n_genes, CELL_MIN_GENES_FLOOR, CELL_QC_LOWER_QUANTILE
    )
    cell_keep = (totals >= min_total) & (n_genes >= min_genes)
    if not np.any(cell_keep):
        raise ValueError("Expression QC removed all cells")

    cell_pos = rep.obs_pos[cell_keep]
    gene_sums, gene_ncells = matrix_col_sum_nnz(rep.adata.X, cell_pos, rep.var_pos)
    min_cells = int(
        max(GENE_MIN_CELLS_FLOOR, np.ceil(GENE_MIN_CELLS_PCT * len(cell_pos)))
    )
    gene_keep = gene_ncells >= min_cells
    if not np.any(gene_keep):
        raise ValueError("Expression QC removed all genes")

    rep.final_obs_pos = cell_pos
    rep.final_obs_names = rep.obs_names[cell_keep]
    rep.final_obs_pos_by_name = pd.Series(rep.final_obs_pos, index=rep.final_obs_names)
    rep.qc_total_counts = totals[cell_keep]
    rep.qc_n_genes_by_counts = n_genes[cell_keep]
    rep.final_var_pos = rep.var_pos[gene_keep]
    rep.final_var_names = rep.var_names[gene_keep]
    rep.final_var_pos_by_name = pd.Series(rep.final_var_pos, index=rep.final_var_names)
    rep.qc_n_cells_by_counts = gene_ncells[gene_keep]
    rep.qc_gene_sums = gene_sums[gene_keep]

    summary = {
        "cells_before": int(len(rep.obs_pos)),
        "cells_after": int(len(rep.final_obs_pos)),
        "cells_removed": int(len(rep.obs_pos) - len(rep.final_obs_pos)),
        "genes_before": int(len(rep.var_pos)),
        "genes_after": int(len(rep.final_var_pos)),
        "genes_removed": int(len(rep.var_pos) - len(rep.final_var_pos)),
        "min_total_counts": int(min_total),
        "min_genes_by_counts": int(min_genes),
        "min_cells_by_counts": int(min_cells),
        "min_cells_pct": float(min_cells / len(rep.final_obs_pos) * 100),
        "cell_qc_lower_quantile": float(CELL_QC_LOWER_QUANTILE),
    }
    log_record("qc_filter", **summary)
    rep.expression_qc_summary = summary
    return summary


def guide_aliases(guide_name):
    aliases = [
        str(part).strip().replace(",", "-")
        for part in str(guide_name).split(";")
        if str(part).strip()
    ]
    return aliases or [str(guide_name).strip().replace(",", "-")]


def guide_target_name(guide_name):
    guide_name = str(guide_name).strip().replace(",", "-")
    return (
        CONTROL_TARGET_SYMBOL
        if guide_name.startswith(CONTROL_TARGET_SYMBOL)
        else guide_name.split("_", 1)[0]
    )


def strip_ensembl_version(value):
    value = str(value).strip()
    return value.split(".", 1)[0] if value.startswith("ENSG") else ""


def guide_target_ensg(guide_name):
    match = re.search(r"ENSG\d+(?:\.\d+)?", str(guide_name))
    return strip_ensembl_version(match.group(0)) if match else ""


def call_info(label):
    label = str(label).strip()
    if label.lower() in {"", "none", "nan"}:
        probes = []
    else:
        probes = sorted(
            p.strip().replace(",", "-") for p in label.split("|") if p.strip()
        )
    gene_names = set()
    gene_ensgs = defaultdict(set)
    for probe in probes:
        for alias in guide_aliases(probe):
            gene = guide_target_name(alias)
            if gene == CONTROL_TARGET_SYMBOL:
                continue
            gene_names.add(gene)
            gene_ensg = guide_target_ensg(alias)
            if gene_ensg:
                gene_ensgs[gene].add(gene_ensg)
    genes = sorted(gene_names)
    controls = sorted(
        {
            probe
            for probe in probes
            if any(
                guide_target_name(alias) == CONTROL_TARGET_SYMBOL
                for alias in guide_aliases(probe)
            )
        }
    )
    if len(genes) == 0:
        outcome = "0_genes"
    elif len(genes) > 1:
        outcome = ">1_gene"
    else:
        outcome = (
            "1_gene_2_probes"
            if sum(
                genes[0] == guide_target_name(alias)
                for p in probes
                for alias in guide_aliases(p)
            )
            > 1
            else "1_gene_1_probe"
        )
    if len(genes) == 0 and controls:
        call_type = "control"
    elif len(genes) == 1 and not controls:
        call_type = "single_gene"
    elif len(genes) > 1 and not controls:
        call_type = "multi_gene"
    elif genes and controls:
        call_type = "mixed_gene_control"
    elif not probes:
        call_type = "unassigned"
    else:
        call_type = "unclassified_probe"
    return {
        "probe_label": "|".join(probes) if probes else "None",
        "n_probes": len(probes),
        "control_probe_label": "|".join(controls) if controls else "None",
        "n_control_probes": len(controls),
        "gene_label": "|".join(genes) if genes else "None",
        "gene_ensg_label": "|".join(
            sorted({ensg for values in gene_ensgs.values() for ensg in values})
        )
        or "None",
        "n_genes": len(genes),
        "outcome": outcome,
        "perturbation_call_type": call_type,
    }


def fit_poisson_gaussian_mixture_from_hist(vals, freqs, max_iter=100, tol=1e-4):
    vals = np.asarray(vals, dtype=float)
    freqs = np.asarray(freqs, dtype=float)
    order = np.argsort(vals)
    vals = vals[order]
    freqs = freqs[order]
    threshold_init = vals[np.searchsorted(np.cumsum(freqs), 0.9 * np.sum(freqs))]
    sig = vals >= threshold_init
    bg = ~sig
    if sig.sum() == 0 or bg.sum() == 0:
        sig = vals >= np.median(vals)
        bg = ~sig
    pi_bg = np.sum(freqs[bg]) / np.sum(freqs)
    lambda_ = max(np.average(vals[bg], weights=freqs[bg]), 0.1)
    mu = np.average(vals[sig], weights=freqs[sig])
    sigma = max(np.sqrt(np.average((vals[sig] - mu) ** 2, weights=freqs[sig])), 1.0)
    log_likelihood = -np.inf
    w_sig = np.zeros_like(vals)
    w_bg = np.ones_like(vals)
    for _ in range(max_iter):
        l_bg = stats.poisson.pmf(vals, mu=lambda_)
        l_sig = stats.norm.pdf(vals, loc=mu, scale=sigma)
        z_bg = pi_bg * l_bg
        z_sig = (1 - pi_bg) * l_sig
        z = np.maximum(z_bg + z_sig, 1e-12)
        w_bg = z_bg / z
        w_sig = z_sig / z
        n_bg = np.sum(freqs * w_bg)
        n_sig = np.sum(freqs * w_sig)
        if n_bg > 0:
            pi_bg = n_bg / np.sum(freqs)
            lambda_ = max(np.sum(freqs * w_bg * vals) / n_bg, 0.1)
        if n_sig > 0:
            mu = np.sum(freqs * w_sig * vals) / n_sig
            sigma = max(np.sqrt(np.sum(freqs * w_sig * (vals - mu) ** 2) / n_sig), 0.1)
        new_ll = np.sum(freqs * np.log(z))
        if abs(new_ll - log_likelihood) < tol:
            break
        log_likelihood = new_ll
    crossover = vals[w_sig > w_bg]
    threshold = max(2, int(np.min(crossover))) if len(crossover) and mu > lambda_ else 5
    return {
        "pi_bg": float(pi_bg),
        "lambda_bg": float(lambda_),
        "mu_sig": float(mu),
        "sigma_sig": float(sigma),
        "decision_threshold": threshold,
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


def call_probe_lists(rep):
    guide_matrix = rep.guide_matrix
    if guide_matrix is None and "guides" in rep.adata.obsm:
        guide_matrix = rep.adata.obsm["guides"]
        rep.guide_matrix = guide_matrix
    if guide_matrix is None:
        return None, None, None
    guide_names = np.asarray(
        rep.adata.uns.get(
            "guide_names", [f"guide_{i}" for i in range(guide_matrix.shape[1])]
        ),
        dtype=str,
    )

    hist = Counter()
    for rows in chunk_array(rep.final_obs_pos):
        block = as_csr(read_block(guide_matrix, rows))
        vals, freqs = np.unique(block.data, return_counts=True)
        hist.update({float(v): int(f) for v, f in zip(vals, freqs) if v > 0})
    mixture = (
        fit_poisson_gaussian_mixture_from_hist(list(hist.keys()), list(hist.values()))
        if hist
        else None
    )
    threshold = mixture["decision_threshold"] if mixture else 5

    labels = []
    total_umis = np.zeros(len(rep.final_obs_pos), dtype=np.float64)
    positive_counts = np.zeros(len(rep.final_obs_pos), dtype=np.int16)
    offset = 0
    for rows in chunk_array(rep.final_obs_pos):
        block = as_csr(read_block(guide_matrix, rows))
        total_umis[offset : offset + len(rows)] = np.asarray(block.sum(axis=1)).ravel()
        for i in range(block.shape[0]):
            row = block.getrow(i)
            keep = row.data >= threshold
            idx = row.indices[keep]
            positive_counts[offset + i] = len(idx)
            labels.append(
                "|".join(sorted(str(guide_names[j]).replace(",", "-") for j in idx))
                if len(idx)
                else "None"
            )
        offset += len(rows)

    labels = pd.Series(labels, index=rep.final_obs_names, name="called_probe_label")
    infos = [call_info(label) for label in labels]
    call_table = pd.DataFrame(infos, index=rep.final_obs_names)
    count_values, count_freqs = np.unique(positive_counts, return_counts=True)
    diagnostics = {
        "method": "gaussian_poisson",
        "count_threshold": int(threshold),
        "n_cells": int(len(labels)),
        "n_guides": int(guide_matrix.shape[1]),
        "n_cells_with_any_called_probe": int((positive_counts >= 1).sum()),
        "n_cells_with_two_or_more_called_probes": int((positive_counts >= 2).sum()),
        "n_cells_with_any_called_gene": int((call_table["n_genes"] >= 1).sum()),
        "n_cells_with_any_called_control_probe": int(
            (call_table["n_control_probes"] >= 1).sum()
        ),
        "n_control_only_calls": int(
            (call_table["perturbation_call_type"] == "control").sum()
        ),
        "n_single_gene_no_control_calls": int(
            (call_table["perturbation_call_type"] == "single_gene").sum()
        ),
        "called_probe_count_distribution": {
            str(int(c)): int(f) for c, f in zip(count_values, count_freqs)
        },
        "gene_call_outcome_counts": {
            outcome: int((call_table["outcome"] == outcome).sum())
            for outcome in GENE_CALL_OUTCOMES
        },
        "perturbation_call_type_counts": {
            str(k): int(v)
            for k, v in call_table["perturbation_call_type"]
            .value_counts()
            .sort_index()
            .items()
        },
        "total_guide_umi_distribution_all_cells": summarize_numeric(total_umis),
        "total_guide_umi_distribution_cells_with_any_called_probe": summarize_numeric(
            total_umis[positive_counts > 0]
        ),
    }
    if mixture is not None:
        diagnostics["gaussian_poisson_parameters"] = mixture

    summary = {
        "n_cells": int(len(labels)),
        "count_threshold": int(threshold),
        "n_single_gene_calls": int(
            (call_table["perturbation_call_type"] == "single_gene").sum()
        ),
        "n_control_calls": int(
            (call_table["perturbation_call_type"] == "control").sum()
        ),
        "n_multi_gene_calls": int((call_table["n_genes"] > 1).sum()),
        "n_no_guide_calls": int((call_table["n_probes"] == 0).sum()),
        "outcome_counts": diagnostics["gene_call_outcome_counts"],
        "perturbation_call_type_counts": diagnostics["perturbation_call_type_counts"],
    }
    log_record("knockout_annotation", **summary)
    rep.called_probe_labels = labels
    rep.call_table = call_table
    rep.guide_call_diagnostics = diagnostics
    rep.knockout_annotation_summary = summary
    rep.guide_names = guide_names
    return labels, call_table, summary


def output_obs(rep):
    obs = rep.obs.iloc[rep.final_obs_pos].copy()
    obs.index = rep.final_obs_names
    obs.index.name = "barcode"
    call_table = rep.call_table
    obs["qc_total_counts"] = rep.qc_total_counts
    obs["qc_n_genes_by_counts"] = rep.qc_n_genes_by_counts
    obs["called_probe_label"] = call_table["probe_label"].astype(str).to_numpy()
    obs["called_probe_count"] = call_table["n_probes"].astype(np.int16).to_numpy()
    obs["called_control_probe_label"] = (
        call_table["control_probe_label"].replace("None", "").astype(str).to_numpy()
    )
    obs["called_control_probe_count"] = (
        call_table["n_control_probes"].astype(np.int16).to_numpy()
    )
    obs["has_called_control_probe"] = obs["called_control_probe_count"].to_numpy() > 0
    obs["called_knockout_genes"] = (
        call_table["gene_label"].replace("None", "").astype(str).to_numpy()
    )
    obs["called_knockout_gene_count"] = (
        call_table["n_genes"].astype(np.int16).to_numpy()
    )
    obs["knockout_call_outcome"] = call_table["outcome"].astype(str).to_numpy()
    obs["perturbation_call_type"] = (
        call_table["perturbation_call_type"].astype(str).to_numpy()
    )
    obs["is_control"] = obs["perturbation_call_type"].to_numpy() == "control"
    obs["is_single_gene_perturbation"] = (
        obs["perturbation_call_type"].to_numpy() == "single_gene"
    )
    obs["perturbed_target_symbol"] = np.where(
        obs["is_single_gene_perturbation"].to_numpy(),
        call_table["gene_label"].to_numpy(),
        "",
    )
    obs["perturbed_target_ensg"] = np.where(
        obs["is_single_gene_perturbation"].to_numpy(),
        call_table["gene_ensg_label"].replace("None", "").to_numpy(),
        "",
    )
    obs["perturbation_call_method"] = "gaussian_poisson"
    obs["dataset_id"] = DATASET_ID
    return obs


def output_var(rep):
    var = rep.var_table.loc[rep.final_var_names].copy()
    var["gene_ensg"] = var["gene_id"].map(strip_ensembl_version)
    var["qc_n_cells_by_counts"] = rep.qc_n_cells_by_counts
    var["qc_pct_cells_by_counts"] = (
        rep.qc_n_cells_by_counts / len(rep.final_obs_pos) * 100
    )
    return var


def create_shell_h5ad(path, obs, var, uns, n_guides, x_dtype):
    output_dir = os.path.dirname(path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    shell = ad.AnnData(
        X=sp.csr_matrix((obs.shape[0], var.shape[0]), dtype=x_dtype),
        obs=obs,
        var=var,
        uns=uns,
    )
    if n_guides:
        shell.obsm["guides"] = sp.csr_matrix((obs.shape[0], n_guides), dtype=x_dtype)
    shell.layers["counts"] = shell.X
    shell.write_h5ad(path, compression="gzip")


def replace_sparse_group(path, key, matrix, row_pos, col_pos=None, dtype=None):
    row_pos = np.asarray(row_pos, dtype=np.int64)
    n_cols = len(col_pos) if col_pos is not None else matrix.shape[1]
    dtype = np.dtype(dtype or getattr(matrix, "dtype", np.float32))
    with h5py.File(path, "r+") as handle:
        parts = key.split("/")
        parent = handle
        for part in parts[:-1]:
            parent = parent.require_group(part)
        if parts[-1] in parent:
            del parent[parts[-1]]
        group = parent.create_group(parts[-1])
        group.attrs["encoding-type"] = "csr_matrix"
        group.attrs["encoding-version"] = "0.1.0"
        group.attrs["shape"] = np.array([len(row_pos), n_cols], dtype=np.int64)
        data_ds = group.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=dtype,
            chunks=True,
            compression="gzip",
        )
        indices_ds = group.create_dataset(
            "indices",
            shape=(0,),
            maxshape=(None,),
            dtype=np.int64,
            chunks=True,
            compression="gzip",
        )
        indptr_ds = group.create_dataset(
            "indptr", shape=(len(row_pos) + 1,), dtype=np.int64, compression="gzip"
        )
        indptr_ds[0] = 0
        nnz = 0
        out_row = 0
        for rows in chunk_array(row_pos):
            block = as_csr(read_block(matrix, rows, col_pos))
            block.eliminate_zeros()
            new_nnz = nnz + block.nnz
            data_ds.resize((new_nnz,))
            indices_ds.resize((new_nnz,))
            data_ds[nnz:new_nnz] = block.data.astype(dtype, copy=False)
            indices_ds[nnz:new_nnz] = block.indices.astype(np.int64, copy=False)
            indptr_ds[out_row + 1 : out_row + len(rows) + 1] = (
                block.indptr[1:].astype(np.int64) + nnz
            )
            nnz = new_nnz
            out_row += len(rows)


def write_filtered_h5ad(rep, path):
    obs = output_obs(rep)
    var = output_var(rep)
    uns = {
        "dataset_id": DATASET_ID,
        "guide_names": rep.guide_names,
        "expression_qc": rep.expression_qc_summary,
        "perturbation_calling": rep.guide_call_diagnostics,
        "downstream_columns": {
            "expression_counts": "layers['counts']",
            "perturbed_target_symbol": "obs['perturbed_target_symbol']",
            "perturbed_target_ensg": "obs['perturbed_target_ensg']",
            "control_indicator": "obs['is_control']",
            "perturbation_indicator": "obs['is_single_gene_perturbation']",
            "control_probe_count": "obs['called_control_probe_count']",
            "all_called_knockout_genes": "obs['called_knockout_genes']",
            "call_outcome": "obs['knockout_call_outcome']",
            "perturbation_call_type": "obs['perturbation_call_type']",
        },
    }
    n_guides = len(rep.guide_names) if hasattr(rep, "guide_names") else 0
    create_shell_h5ad(
        path, obs, var, uns, n_guides, getattr(rep.adata.X, "dtype", np.float32)
    )
    replace_sparse_group(path, "X", rep.adata.X, rep.final_obs_pos, rep.final_var_pos)
    if n_guides:
        replace_sparse_group(path, "obsm/guides", rep.guide_matrix, rep.final_obs_pos)
    with h5py.File(path, "r+") as handle:
        if "counts" in handle["layers"]:
            del handle["layers"]["counts"]
        handle["layers"]["counts"] = handle["X"]
    log_record(
        "filtered_h5ad_written",
        path=path,
        cells=len(rep.final_obs_pos),
        genes=len(rep.final_var_pos),
        compression="gzip",
    )


def value_by_cells(ds, cells, names, transformed_ok=False):
    if names == "total":
        column, values = find_numeric_column(ds.obs, CELL_TOTAL_COUNT_COLUMNS)
    else:
        column, values = find_numeric_column(
            ds.obs,
            CELL_DETECTED_GENE_COLUMNS,
            required_name_terms=["gene", "feature", "detected"],
        )
    if values is not None:
        pos = ds.obs_pos_by_name.loc[cells].to_numpy()
        return values.iloc[pos].to_numpy(dtype=float), f"obs['{column}']"
    if (
        ds.expression_profile["value_kind"] == "signed_transformed"
        and not transformed_ok
    ):
        return np.full(len(cells), np.nan), "unavailable for signed transformed X"
    pos = ds.obs_pos_by_name.loc[cells].to_numpy()
    sums, nnz = matrix_row_sum_nnz(ds.adata.X, pos, ds.var_pos)
    return (sums if names == "total" else nnz.astype(float)), "X row metrics"


def rep_value_by_cells(rep, cells, names):
    rel = (
        pd.Series(np.arange(len(rep.final_obs_names)), index=rep.final_obs_names)
        .loc[cells]
        .to_numpy()
    )
    if names == "total":
        return rep.qc_total_counts[rel], "obs['qc_total_counts']"
    return rep.qc_n_genes_by_counts[rel].astype(float), "obs['qc_n_genes_by_counts']"


def gene_metric_values(ds, genes, metric):
    if metric == "mean":
        column, values = find_numeric_column(
            ds.var, GENE_MEAN_COLUMNS, require_nonnegative=False
        )
    else:
        column, values = find_numeric_column(
            ds.var,
            GENE_DETECTED_CELL_COLUMNS,
            required_name_terms=["cell", "detected", "expressed"],
        )
    if values is not None:
        pos = ds.var_pos_by_name.loc[genes].to_numpy()
        vals = values.iloc[pos].to_numpy(dtype=float)
        if metric == "dropout":
            vals = 100 - vals / len(ds.obs_pos) * 100
        return vals, f"var['{column}']"
    if (
        metric == "dropout"
        and ds.expression_profile["value_kind"] == "signed_transformed"
    ):
        return np.full(len(genes), np.nan), "unavailable for signed transformed X"
    cell_pos = ds.obs_pos_by_name.loc[ds.obs_names].to_numpy()
    gene_pos = ds.var_pos_by_name.loc[genes].to_numpy()
    sums, nnz = matrix_col_sum_nnz(ds.adata.X, cell_pos, gene_pos)
    return (
        sums / len(cell_pos) if metric == "mean" else 100 - nnz / len(cell_pos) * 100
    ), "X column metrics"


def rep_gene_metric_values(rep, genes, metric):
    rel = (
        pd.Series(np.arange(len(rep.final_var_names)), index=rep.final_var_names)
        .loc[genes]
        .to_numpy()
    )
    if metric == "mean":
        return rep.qc_gene_sums[rel] / len(rep.final_obs_pos), "X column means after QC"
    return (
        100 - rep.qc_n_cells_by_counts[rel] / len(rep.final_obs_pos) * 100,
        "X detection after QC",
    )


def subsample_df(df, max_points=SCATTER_MAX_POINTS):
    if len(df) <= max_points:
        return df
    return df.sample(max_points, random_state=RANDOM_SEED)


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
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])
    df[x_col] = pd.to_numeric(df[x_col], errors="coerce")
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])
    if log_scale:
        df = df[(df[x_col] > -1) & (df[y_col] > -1)]
    if len(df) < 2:
        return {
            "pearson": np.nan,
            "spearman": np.nan,
            "pct_deviant": np.nan,
            "skipped": True,
            "reason": "fewer than two finite comparable values",
            "n_points": int(len(df)),
        }
    total_points = len(df)
    df = subsample_df(df)
    deviation = df[[x_col, y_col]].astype(float)
    if deviation_on_log:
        deviation = np.log1p(deviation)
    rel_diff = np.abs(deviation[x_col] - deviation[y_col]) / deviation.mean(
        axis=1
    ).replace(0, 1)
    pct_deviant = float((rel_diff > 0.1).mean() * 100)
    pearson, _ = (
        stats.pearsonr(df[x_col], df[y_col])
        if df[x_col].nunique() > 1 and df[y_col].nunique() > 1
        else (np.nan, None)
    )
    spearman, _ = (
        stats.spearmanr(df[x_col], df[y_col])
        if df[x_col].nunique() > 1 and df[y_col].nunique() > 1
        else (np.nan, None)
    )

    plot_df = df.copy()
    if log_scale:
        plot_df[x_col] += 1
        plot_df[y_col] += 1
    fig, ax = plt.subplots(figsize=(10, 11))
    ax.scatter(
        plot_df[x_col], plot_df[y_col], alpha=0.05, s=1, color="teal", rasterized=True
    )
    min_val = min(plot_df[x_col].min(), plot_df[y_col].min())
    max_val = max(plot_df[x_col].max(), plot_df[y_col].max())
    ax.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="Identity")
    ax.set_title(title, fontsize=18, fontweight="bold", pad=35)
    label = "Deviants (>10% on log1p values)" if deviation_on_log else "Deviants (>10%)"
    ax.text(
        0.5,
        1.02,
        f"Pearson r = {pearson:.4f} | Spearman rho = {spearman:.4f} | {label}: {pct_deviant:.1f}% | n = {len(df):,}/{total_points:,}",
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
    fig.text(
        0.5,
        0.02,
        "\n".join(textwrap.wrap(description, width=100)),
        ha="center",
        va="bottom",
        fontsize=10,
        linespacing=1.4,
    )
    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    return {
        "pearson": float(pearson),
        "spearman": float(spearman),
        "pct_deviant": pct_deviant,
        "n_points": int(len(df)),
        "n_total_points": int(total_points),
    }


def transform_expression_block(block, profile):
    if not profile["raw_counts_detected"]:
        return block
    if sp.issparse(block):
        block = block.astype(np.float64).tocsr(copy=True)
        row_sums = np.asarray(block.sum(axis=1)).ravel()
        scale = np.divide(
            TARGET_SUM, row_sums, out=np.zeros_like(row_sums), where=row_sums > 0
        )
        block = block.multiply(scale[:, None]).tocsr()
        block.data = np.log1p(block.data)
        return block
    block = np.asarray(block, dtype=np.float64)
    row_sums = block.sum(axis=1)
    scale = np.divide(
        TARGET_SUM, row_sums, out=np.zeros_like(row_sums), where=row_sums > 0
    )
    return np.log1p(block * scale[:, None])


def cell_correlations(cur, rep, common_cells, common_genes):
    n = min(CELL_CORR_SAMPLE_SIZE, len(common_cells))
    if n == 0 or len(common_genes) < 2:
        return np.array([])
    rng = np.random.default_rng(RANDOM_SEED)
    sample_idx = np.sort(rng.choice(len(common_cells), n, replace=False))
    cells = common_cells[sample_idx]
    cur_rows = cur.obs_pos_by_name.loc[cells].to_numpy()
    rep_rows = rep.final_obs_pos_by_name.loc[cells].to_numpy()
    cur_cols = cur.var_pos_by_name.loc[common_genes].to_numpy()
    rep_cols = rep.final_var_pos_by_name.loc[common_genes].to_numpy()
    corrs = []
    for i in range(0, n, 256):
        c_block = transform_expression_block(
            read_block(cur.adata.X, cur_rows[i : i + 256], cur_cols),
            cur.expression_profile,
        )
        r_block = transform_expression_block(
            read_block(rep.adata.X, rep_rows[i : i + 256], rep_cols),
            rep.expression_profile,
        )
        c_arr = c_block.toarray() if sp.issparse(c_block) else np.asarray(c_block)
        r_arr = r_block.toarray() if sp.issparse(r_block) else np.asarray(r_block)
        for c_row, r_row in zip(c_arr, r_arr):
            mask = np.isfinite(c_row) & np.isfinite(r_row)
            if mask.sum() < 2:
                corrs.append(np.nan)
            else:
                c_f = c_row[mask]
                r_f = r_row[mask]
                if np.std(c_f) == 0 or np.std(r_f) == 0:
                    corrs.append(np.nan)
                else:
                    corrs.append(stats.pearsonr(c_f, r_f)[0])
    return np.asarray(corrs, dtype=float)


def outcome_count_dict(counter):
    return {outcome: int(counter.get(outcome, 0)) for outcome in GENE_CALL_OUTCOMES}


def plot_gene_outcome_matrix(matrix_dict, single_gene_match, count_threshold):
    matrix_df = (
        pd.DataFrame.from_dict(matrix_dict, orient="index")
        .reindex(index=GENE_CALL_OUTCOMES, columns=GENE_CALL_OUTCOMES, fill_value=0)
        .astype(int)
    )
    plot_path = f"comparison_results/{DATASET_ID}/perturbation_gene_outcome_matrix_gaussian_poisson.png"
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 7))
    sns.heatmap(
        matrix_df,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        cbar_kws={"label": "Cells"},
        xticklabels=["0 genes", "1 gene\n1 probe", "1 gene\n2+ probes", ">1 gene"],
        yticklabels=["0 genes", "1 gene\n1 probe", "1 gene\n2+ probes", ">1 gene"],
        ax=ax,
    )
    same = single_gene_match["n_same_gene"]
    total = single_gene_match["n_cells_both_single_gene"]
    ax.set_title(
        f"Gaussian Poisson: Gene-Call Outcome Matrix\nthreshold >= {count_threshold} UMI | "
        f"same single-gene calls: {same}/{total} ({same / total * 100 if total else 0:.1f}%)",
        fontsize=13,
        pad=12,
    )
    ax.set_xlabel("Reprocessed outcome")
    ax.set_ylabel("Curated outcome")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return plot_path


def compare_perturbations(cur, rep, common_cells):
    column = detect_perturbation_column(cur.obs)
    if not column or not hasattr(rep, "called_probe_labels"):
        return None
    cur_values = pd.Series(
        cur.obs.iloc[cur.obs_pos_by_name.loc[common_cells].to_numpy()][column]
        .astype(str)
        .to_numpy(),
        index=common_cells,
    )
    rep_values = rep.called_probe_labels.loc[common_cells]
    cur_outcomes = Counter()
    rep_outcomes = Counter()
    matrix = defaultdict(Counter)
    call_type_cur = Counter()
    call_type_rep = Counter()
    cur_probe = cur_gene = cur_control = cur_single = 0
    rep_probe = rep_gene = rep_control = rep_single = 0
    both_single = same = different = 0
    for cur_label, rep_label in zip(cur_values, rep_values):
        c = call_info(cur_label)
        r = call_info(rep_label)
        cur_outcomes[c["outcome"]] += 1
        rep_outcomes[r["outcome"]] += 1
        matrix[c["outcome"]][r["outcome"]] += 1
        call_type_cur[c["perturbation_call_type"]] += 1
        call_type_rep[r["perturbation_call_type"]] += 1
        cur_probe += c["n_probes"] > 0
        cur_gene += c["n_genes"] > 0
        cur_control += c["n_control_probes"] > 0
        cur_single += c["perturbation_call_type"] == "single_gene"
        rep_probe += r["n_probes"] > 0
        rep_gene += r["n_genes"] > 0
        rep_control += r["n_control_probes"] > 0
        rep_single += r["perturbation_call_type"] == "single_gene"
        if (
            c["perturbation_call_type"] == "single_gene"
            and r["perturbation_call_type"] == "single_gene"
        ):
            both_single += 1
            same += c["gene_label"] == r["gene_label"]
            different += c["gene_label"] != r["gene_label"]
    curated = {
        "perturbation_column": column,
        "n_common_cells": int(len(common_cells)),
        "n_cells_with_any_probe": int(cur_probe),
        "n_cells_with_any_gene": int(cur_gene),
        "n_cells_with_any_control_probe": int(cur_control),
        "n_valid_single_gene_no_control_cells": int(cur_single),
        "outcome_counts": outcome_count_dict(cur_outcomes),
        "perturbation_call_type_counts": {
            str(k): int(v) for k, v in call_type_cur.items()
        },
    }
    log_record("curated_gene_call_summary", **curated)
    matrix_dict = {
        c: {r: int(matrix[c][r]) for r in GENE_CALL_OUTCOMES}
        for c in GENE_CALL_OUTCOMES
    }
    single_gene_match = {
        "n_cells_curated_single_gene": int(cur_single),
        "n_cells_reprocessed_single_gene": int(rep_single),
        "n_cells_both_single_gene": int(both_single),
        "n_same_gene": int(same),
        "n_different_gene": int(different),
        "pct_same_gene": float(same / both_single * 100) if both_single else 0.0,
        "pct_different_gene": (
            float(different / both_single * 100) if both_single else 0.0
        ),
    }
    model = {
        "method": "gaussian_poisson",
        "count_threshold": rep.guide_call_diagnostics.get("count_threshold"),
        "probe_call_diagnostics": rep.guide_call_diagnostics,
        "reprocessed_common_call_summary": {
            "n_common_cells": int(len(common_cells)),
            "n_cells_with_any_probe": int(rep_probe),
            "n_cells_with_any_gene": int(rep_gene),
            "n_cells_with_any_control_probe": int(rep_control),
            "n_valid_single_gene_no_control_cells": int(rep_single),
        },
        "reprocessed_outcome_counts": outcome_count_dict(rep_outcomes),
        "outcome_matrix_rows_curated_columns_reprocessed": matrix_dict,
        "single_gene_match": single_gene_match,
        "plot_path": plot_gene_outcome_matrix(
            matrix_dict,
            single_gene_match,
            rep.guide_call_diagnostics.get("count_threshold"),
        ),
    }
    log_record("perturbation_gene_comparison", model="gaussian_poisson", **model)
    return {
        "curated": curated,
        "models": {"gaussian_poisson": model},
        "control_definition": "control means one or more non-targeting probes called and zero gene-targeting probes called",
        "valid_single_gene_definition": "exactly one gene-targeting gene called and zero non-targeting probes called",
    }


def plot_overlap_distributions(rep, common_cells, common_genes):
    common_cell_mask = rep.final_obs_names.isin(common_cells)
    common_gene_mask = rep.final_var_names.isin(common_genes)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
    sns.histplot(
        rep.qc_total_counts[common_cell_mask],
        label="Common Cells",
        color="blue",
        alpha=0.5,
        log_scale=True,
        ax=ax1,
    )
    sns.histplot(
        rep.qc_total_counts[~common_cell_mask],
        label="Reprocessed-only",
        color="orange",
        alpha=0.5,
        log_scale=True,
        ax=ax1,
    )
    ax1.set_title("Total UMIs per Cell")
    ax1.legend()
    sns.histplot(
        rep.qc_n_cells_by_counts[common_gene_mask],
        label="Common Genes",
        color="blue",
        alpha=0.5,
        log_scale=True,
        ax=ax2,
    )
    sns.histplot(
        rep.qc_n_cells_by_counts[~common_gene_mask],
        label="Reprocessed-only",
        color="orange",
        alpha=0.5,
        log_scale=True,
        ax=ax2,
    )
    ax2.set_title("Detection Freq per Gene")
    ax2.legend()
    sns.histplot(
        rep.qc_n_genes_by_counts[common_cell_mask],
        label="Common Cells",
        color="blue",
        alpha=0.5,
        log_scale=True,
        ax=ax3,
    )
    sns.histplot(
        rep.qc_n_genes_by_counts[~common_cell_mask],
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
    plt.close(fig)


def main():
    os.makedirs(f"comparison_results/{DATASET_ID}", exist_ok=True)
    log_record(
        "input_paths",
        curated_h5ad_path=CURATED_H5AD_PATH,
        reprocessed_h5ad_path=REPROCESSED_H5AD_PATH,
        filtered_reprocessed_h5ad_path=FILTERED_REPROCESSED_H5AD_PATH,
    )
    cur = DatasetView(CURATED_H5AD_PATH, "Curated")
    rep = DatasetView(REPROCESSED_H5AD_PATH, "Reprocessed")
    try:
        qc_summary = filter_low_signal_cells_and_genes(rep)
        _, _, knockout_summary = call_probe_lists(rep)
        write_filtered_h5ad(rep, FILTERED_REPROCESSED_H5AD_PATH)

        common_cells = np.intersect1d(cur.obs_names, rep.final_obs_names)
        common_genes = np.intersect1d(cur.var_names, rep.final_var_names)
        overlap_summary = {
            "common_cells": int(len(common_cells)),
            "common_genes": int(len(common_genes)),
            "curated_cells": int(len(cur.obs_names)),
            "curated_genes": int(len(cur.var_names)),
            "reprocessed_cells": int(len(rep.final_obs_names)),
            "reprocessed_genes": int(len(rep.final_var_names)),
            "common_cells_pct_of_curated": (
                float(len(common_cells) / len(cur.obs_names) * 100)
                if len(cur.obs_names)
                else 0.0
            ),
            "common_genes_pct_of_curated": (
                float(len(common_genes) / len(cur.var_names) * 100)
                if len(cur.var_names)
                else 0.0
            ),
        }
        log_record("overlap_summary", **overlap_summary)

        plot_overlap_distributions(rep, common_cells, common_genes)

        cur_total, cur_total_src = value_by_cells(cur, common_cells, "total")
        rep_total, rep_total_src = rep_value_by_cells(rep, common_cells, "total")
        cur_ngenes, cur_ngenes_src = value_by_cells(cur, common_cells, "n_genes")
        rep_ngenes, rep_ngenes_src = rep_value_by_cells(rep, common_cells, "n_genes")
        cur_mean, cur_mean_src = gene_metric_values(cur, common_genes, "mean")
        rep_mean, rep_mean_src = rep_gene_metric_values(rep, common_genes, "mean")
        cur_dropout, cur_dropout_src = gene_metric_values(cur, common_genes, "dropout")
        rep_dropout, rep_dropout_src = rep_gene_metric_values(
            rep, common_genes, "dropout"
        )

        log_record(
            "comparison_metric_source",
            dataset="Curated",
            sources={
                "cell_total_counts": cur_total_src,
                "cell_detected_genes": cur_ngenes_src,
                "gene_mean_expression": cur_mean_src,
                "gene_dropout_pct": cur_dropout_src,
            },
        )
        log_record(
            "comparison_metric_source",
            dataset="Reprocessed",
            sources={
                "cell_total_counts": rep_total_src,
                "cell_detected_genes": rep_ngenes_src,
                "gene_mean_expression": rep_mean_src,
                "gene_dropout_pct": rep_dropout_src,
            },
        )

        results_summary = {
            "input_paths": {
                "curated_h5ad_path": CURATED_H5AD_PATH,
                "reprocessed_h5ad_path": REPROCESSED_H5AD_PATH,
                "filtered_reprocessed_h5ad_path": FILTERED_REPROCESSED_H5AD_PATH,
            },
            "expression_qc": qc_summary,
            "knockout_annotation": knockout_summary,
            "overlap": overlap_summary,
            "preprocessing": {
                "curated": cur.expression_profile,
                "reprocessed": rep.expression_profile,
            },
            "metric_sources": {
                "curated": {
                    "cell_total_counts": cur_total_src,
                    "cell_detected_genes": cur_ngenes_src,
                    "gene_mean_expression": cur_mean_src,
                    "gene_dropout_pct": cur_dropout_src,
                },
                "reprocessed": {
                    "cell_total_counts": rep_total_src,
                    "cell_detected_genes": rep_ngenes_src,
                    "gene_mean_expression": rep_mean_src,
                    "gene_dropout_pct": rep_dropout_src,
                },
            },
            "cell_metrics": {},
            "gene_metrics": {},
        }

        results_summary["cell_metrics"]["total_counts"] = plot_scatter_comparison(
            pd.DataFrame({"cur": cur_total, "rep": rep_total}),
            "cur",
            "rep",
            "Total UMI Counts",
            "Original",
            "Reprocessed",
            f"comparison_results/{DATASET_ID}/counts_comparison.png",
            log_scale=True,
            deviation_on_log=True,
            description="Each point is one shared cell. For large datasets, points and correlations are computed on a deterministic sample.",
        )
        log_record(
            "comparison_metric",
            metric_group="cell_metrics",
            metric="total_counts",
            **results_summary["cell_metrics"]["total_counts"],
        )
        results_summary["cell_metrics"]["n_genes"] = plot_scatter_comparison(
            pd.DataFrame({"cur": cur_ngenes, "rep": rep_ngenes}),
            "cur",
            "rep",
            "Number of Detected Genes",
            "Original",
            "Reprocessed",
            f"comparison_results/{DATASET_ID}/genes_comparison.png",
            description="Each point is one shared cell. For signed transformed curated matrices this metric is skipped unless metadata provides it.",
        )
        log_record(
            "comparison_metric",
            metric_group="cell_metrics",
            metric="n_genes",
            **results_summary["cell_metrics"]["n_genes"],
        )
        gene_df = pd.DataFrame(
            {
                "mean_cur": cur_mean,
                "mean_rep": rep_mean,
                "dropout_cur": cur_dropout,
                "dropout_rep": rep_dropout,
            }
        )
        results_summary["gene_metrics"]["mean_expression"] = plot_scatter_comparison(
            gene_df,
            "mean_cur",
            "mean_rep",
            "Mean Gene Expression",
            "Original",
            "Reprocessed",
            f"comparison_results/{DATASET_ID}/gene_expression_mean.png",
            log_scale=True,
            description="Each point is one shared gene.",
        )
        log_record(
            "comparison_metric",
            metric_group="gene_metrics",
            metric="mean_expression",
            **results_summary["gene_metrics"]["mean_expression"],
        )
        results_summary["gene_metrics"]["dropout_rate"] = plot_scatter_comparison(
            gene_df,
            "dropout_cur",
            "dropout_rep",
            "Gene Dropout Rate (%)",
            "Original",
            "Reprocessed",
            f"comparison_results/{DATASET_ID}/sparsity_comparison.png",
            description="Each point is one shared gene. For signed transformed curated matrices this metric is skipped unless metadata provides it.",
        )
        log_record(
            "comparison_metric",
            metric_group="gene_metrics",
            metric="dropout_rate",
            **results_summary["gene_metrics"]["dropout_rate"],
        )

        cell_corrs = cell_correlations(cur, rep, common_cells, common_genes)
        results_summary["cell_wise_corr"] = {
            "median": float(np.nanmedian(cell_corrs)) if len(cell_corrs) else np.nan,
            "mean": float(np.nanmean(cell_corrs)) if len(cell_corrs) else np.nan,
            "n_cells_sampled": int(len(cell_corrs)),
        }
        log_record("cell_wise_correlation", **results_summary["cell_wise_corr"])
        if len(cell_corrs):
            os.makedirs(f"comparison_results/{DATASET_ID}/supplementary", exist_ok=True)
            fig, ax = plt.subplots(figsize=(9, 6))
            sns.histplot(
                cell_corrs[np.isfinite(cell_corrs)],
                bins=100,
                kde=True,
                color="purple",
                alpha=0.4,
                ax=ax,
            )
            ax.axvline(
                np.nanmedian(cell_corrs), color="red", linestyle="dashed", linewidth=2
            )
            ax.set_title(
                "Preservation of Single-Cell Identity\n(Cell-wise Pearson Correlation)",
                fontsize=16,
                fontweight="bold",
                pad=15,
            )
            ax.set_xlabel(
                "Pearson Correlation (Original vs Reprocessed cell)", fontsize=12
            )
            ax.set_ylabel("Number of Cells", fontsize=12)
            plt.tight_layout()
            plt.savefig(
                f"comparison_results/{DATASET_ID}/supplementary/cellwise_correlation_dist.png",
                bbox_inches="tight",
                dpi=300,
            )
            plt.close(fig)

        results_summary["perturbation"] = compare_perturbations(cur, rep, common_cells)
        report_path = f"comparison_results/{DATASET_ID}/summary_report.txt"
        with open(report_path, "w") as handle:
            handle.write(json.dumps(results_summary, indent=4, default=json_default))
        log_record("summary_report_written", path=report_path)

        pert = results_summary["perturbation"]
        gp_match = (
            pert["models"]["gaussian_poisson"]["single_gene_match"] if pert else None
        )
        display(
            pd.DataFrame(
                [
                    {
                        "Metric": "Cell Overlap",
                        "Value": f"{len(common_cells)} ({len(common_cells) / len(cur.obs_names):.1%})",
                    },
                    {
                        "Metric": "Gene Overlap",
                        "Value": f"{len(common_genes)} ({len(common_genes) / len(cur.var_names):.1%})",
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
                        "Value": (
                            f"{gp_match['pct_same_gene']:.2f}% ({gp_match['n_same_gene']}/{gp_match['n_cells_both_single_gene']})"
                            if gp_match
                            else "N/A"
                        ),
                    },
                ]
            )
        )
    finally:
        cur.close()
        rep.close()


if __name__ == "__main__":
    main()
