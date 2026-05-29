# Nadig 2025 Jurkat: Curated vs. Reprocessed Comparison

This directory contains tools to compare two versions of the Nadig 2025 Jurkat Perturb-seq dataset:
1.  **Curated:** Author-supplied H5AD (`gs://${LAKE_BUCKET}/perturbseq/curated/nadig_2025_jurkat_curated.h5ad`).
2.  **Reprocessed:** FASTQ-to-H5AD reprocessed version (`gs://${LAKE_BUCKET}/perturbseq/fastq-reprocess/nadig_2025_jurkat.h5ad`).

## 1. Environment Setup (Google Cloud)

### Create a Vertex AI Workbench Instance
1.  **Source your development environment variables:**
    ```bash
    dev_secrets
    ```

2.  **Create the instance via CLI:**
    Run the following command in your terminal to create a high-memory Workbench instance (at least 64GB RAM recommended for full loading):

    ```bash
    gcloud workbench instances create nadig-comparison-notebook \
        --project=$GCLOUD_PROJECT \
        --location=$GCLOUD_ZONE \
        --machine-type=n1-highmem-32 \
        --vm-image-project=cloud-notebooks-managed \
        --vm-image-family=workbench-instances \
        --data-disk-size=1500
    ```

### Access JupyterLab
1.  Once the instance is "Active", go to the [Vertex AI Workbench Console](https://console.cloud.google.com/vertex-ai/workbench).
2.  Click **OPEN JUPYTERLAB** next to your instance name.

## 2. Running the Comparison

1.  Clone the repository or copy the `comparison.py` script to the VM.
2.  Install required dependencies:
    ```bash
    pip install scanpy pandas numpy matplotlib seaborn scipy
    ```
3.  Execute the script:
    ```bash
    python3 comparison.py
    ```

## 3. Metrics and Visualizations

The script generates the following outputs in the `comparison_results/` folder:

### Overall Metrics
*   **counts_comparison.png**: Scatter plot of total UMI counts per cell (normalized and log-transformed if necessary to match scales).
*   **genes_comparison.png**: Scatter plot of the number of unique genes detected per cell.
*   **summary_report.txt**: Quantitative summary of overlap and correlations (Pearson and Spearman).

### Gene Expression and Structure
*   **gene_expression_mean.png**: Correlation of average expression values for all common genes.
*   **pca_comparison.png**: Side-by-side PCA plots to compare the global structure and variance of the two datasets.

### Perturbation Assignment
*   **perturbation_confusion_matrix.png**: Heatmap showing the overlap of dual same-target guide assignments between the two versions for the top 20 most frequent perturbations.
*   **summary_report.txt**: Includes guide matrix diagnostics for the reprocessed data, including guide UMI distributions, explicit non-targeting control calls, valid single-gene/no-control calls, mixed gene/control calls, and multi-gene calls.

## 4. How the Pipeline Works
*   **Full Loading**: Datasets are loaded fully into memory for faster processing and more complex analyses.
*   **Auto-Normalization**: The script detects if datasets are raw counts or log-normalized and applies necessary transformations to ensure they are on a comparable scale.
*   **Aggressive Alignment**: Gene names are aligned even if they are stored in different `var` columns (e.g., `gene_symbols` vs index).
*   **Structural Validation**: PCA is used to verify that the reprocessed data preserves the biological structure of the original curated dataset.
*   **Control Annotation**: `non-targeting_*` guides are recorded separately from gene-targeting guides. A control cell is one with at least one non-targeting guide and zero gene-targeting guides. A valid perturbation cell is one with exactly one gene-targeting gene and zero non-targeting guides.
*   **Gene Symbols**: The filtered H5AD stores expression feature symbols in `var["gene_symbol"]` and uses symbol-based `var_names`. Author-supplied `var["gene_name"]` is used when present; otherwise symbols are resolved from `/hps/nobackup/mfreeberg/cache/reference/Homo_sapiens.GRCh38.115.gtf.gz`. The script fails if that GTF is missing or if gene-ID-only features cannot be resolved to symbols.

# Raw data from source
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667%5Fjurkat%5Fraw%5Fsinglecell%5F01.h5ad
```
