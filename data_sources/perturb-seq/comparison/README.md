# Perturb-seq Curated vs. Reprocessed Comparison

## Running the pipeline

### Install required dependencies:
```bash
pip install anndata h5py pandas numpy matplotlib seaborn scipy
```

### Run

`datasets.txt` should contain list of dataset IDs, one per line.

```bash
mkdir -p logs
cat datasets.txt | parallel \
  sbatch --mem=128G --time=12:00:00 \
    --output="logs/{}.comparison.log" \
    --error="logs/{}.comparison.err" \
    --wrap=\"python3 comparison.py {}\"
```

## Metrics and Visualizations

The script generates the following outputs in the `comparison_results/` folder:

### Overall Metrics
*   **counts_comparison.png**: Scatter plot of total UMI counts per shared cell. Very large cell-level plots use a deterministic sample.
*   **genes_comparison.png**: Scatter plot of the number of unique genes detected per cell.
*   **summary_report.txt**: Quantitative summary of overlap and correlations (Pearson and Spearman).

### Gene Expression and Structure
*   **gene_expression_mean.png**: Correlation of average expression values for all common genes.
*   **supplementary/cellwise_correlation_dist.png**: Distribution of sampled cell-wise expression correlations.

### Perturbation Assignment
*   **perturbation_gene_outcome_matrix_gaussian_poisson.png**: Heatmap showing curated-vs-reprocessed perturbation gene-call outcomes.
*   **summary_report.txt**: Includes guide matrix diagnostics for the reprocessed data, including guide UMI distributions, explicit non-targeting control calls, valid single-gene/no-control calls, mixed gene/control calls, and multi-gene calls.

## 4. How the Pipeline Works
*   **Backed/Chunked Loading**: H5AD files are opened in backed mode. Expression and guide matrices are processed in row chunks controlled by `PERTURBSEQ_ROW_CHUNK_SIZE` instead of loading full matrices into RAM.
*   **Filtered Output**: The filtered reprocessed H5AD is written as on-disk CSR chunks. `layers["counts"]` is a hard link to `X`, avoiding a second on-disk copy of the expression matrix.
*   **Auto-Normalization**: The script samples expression values to detect raw integer counts versus already transformed expression. Raw counts are total-normalized and log-transformed; signed transformed matrices, such as gemgroup Z-normalized Replogle expression, are used as supplied to avoid invalid `log1p` transformations.
*   **Aggressive Alignment**: Gene names are aligned even if they are stored in different `var` columns (e.g., `gene_symbols` vs index).
*   **Structural Validation**: Sampled cell-wise correlations are computed on shared cells and genes. Raw count matrices are normalized/log-transformed per sampled chunk before correlation.
*   **Control Annotation**: `non-targeting_*` guides are recorded separately from gene-targeting guides. A control cell is one with at least one non-targeting guide and zero gene-targeting guides. A valid perturbation cell is one with exactly one gene-targeting gene and zero non-targeting guides.
*   **Gene Symbols**: The filtered H5AD stores expression feature symbols in `var["gene_symbol"]` and uses symbol-based `var_names`. Author-supplied `var["gene_name"]` is used when present; otherwise symbols are resolved from `/hps/nobackup/mfreeberg/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz`. The script fails if that GTF is missing. Gene-ID-only reprocessed features whose GTF records do not contain `gene_name` are removed before QC.
*   **Metric Sources**: Cell and gene metric columns are detected from metadata when available (`obs["UMI_count"]`, `obs["qc_total_counts"]`, `var["mean"]`, etc.). Metrics that cannot be recovered from a transformed curated matrix, such as detected genes or dropout, are reported as skipped instead of inferred from signed values.
*   **Sampling Controls**: `PERTURBSEQ_SCATTER_MAX_POINTS`, `PERTURBSEQ_CELL_CORR_SAMPLE_SIZE`, and `PERTURBSEQ_RANDOM_SEED` control deterministic sampling for large cell-level visualizations and correlations.
