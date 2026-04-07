import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import sys
import gc


def download_file(gs_path, local_path):
    """Downloads a file from Google Cloud Storage using gsutil."""
    if os.path.exists(local_path):
        print(f"File {local_path} already exists. Skipping download.")
    else:
        print(f"Downloading {gs_path} to {local_path}...")
        # -m for multi-threaded/multi-processing copy
        ret = os.system(f"gsutil -m cp {gs_path} {local_path}")
        if ret != 0:
            print(f"Error: gsutil failed with exit code {ret}")
            sys.exit(1)


def aggregate_reprocessed(input_path, output_path):
    """
    Loads the reprocessed H5AD, aggregates counts by cell barcode (summing SRR runs),
    and saves the result. Optimized for memory efficiency.
    """
    print(f"Loading {input_path}...")
    # Load with backed mode initially to check size/metadata if needed,
    # but for matrix multiplication we'll need it in memory or use a chunked approach.
    # Given the original script loaded it fully, we'll do the same but with GC care.
    adata = sc.read_h5ad(input_path)

    print("Aggregating counts by barcode (summing multiple runs per cell)...")

    # Extract barcode (part before '-') from index
    # This is a common pattern for datasets with multiple SRR runs per cell
    adata.obs["barcode"] = adata.obs.index.str.split("-").str[0]

    # Grouping barcodes
    unique_barcodes, group_indices = np.unique(
        adata.obs["barcode"], return_inverse=True
    )
    n_groups = len(unique_barcodes)

    print(f"Found {n_groups} unique barcodes from {adata.n_obs} total observations.")

    # Create an aggregation matrix: (n_groups x n_obs)
    # This is a sparse mapping matrix that sums observations belonging to the same barcode.
    aggregation_matrix = sp.csr_matrix(
        (
            np.ones(adata.n_obs, dtype=np.float32),
            (group_indices, np.arange(adata.n_obs)),
        ),
        shape=(n_groups, adata.n_obs),
    )

    # Sum counts across runs for each barcode
    # sparse @ sparse or sparse @ dense is generally efficient and multi-threaded in some environments
    print("Performing matrix multiplication for aggregation...")
    summed_X = aggregation_matrix @ adata.X

    # To be memory efficient, we extract metadata before creating the new object
    # groupby().first() is used to keep the first occurrence of metadata for each barcode
    print("Aggregating observation metadata...")
    obs_summed = adata.obs.groupby("barcode").first()

    # Preserve variable metadata
    var_copy = adata.var.copy()

    # Create the new aggregated AnnData object
    adata_sum = sc.AnnData(X=summed_X, obs=obs_summed, var=var_copy)
    adata_sum.obs_names = unique_barcodes

    # Free original data as soon as possible
    del adata
    gc.collect()

    print(f"Saving aggregated data to {output_path}...")
    adata_sum.write(output_path)
    print("Aggregation complete.")


if __name__ == "__main__":
    # Ensure LAKE_BUCKET is available
    lake_bucket = os.environ.get("LAKE_BUCKET")
    if not lake_bucket:
        print("Error: LAKE_BUCKET environment variable is not set.")
        print(
            "Please run 'export LAKE_BUCKET=your-bucket-name' or source your secrets."
        )
        sys.exit(1)

    reprocessed_gs = (
        f"gs://{lake_bucket}/perturbseq/fastq-reprocess/nadig_2025_jurkat.h5ad"
    )
    reprocessed_local = "nadig_2025_jurkat_reprocessed.h5ad"
    output_local = "nadig_2025_jurkat_reprocessed_summed.h5ad"

    # 1. Download
    download_file(reprocessed_gs, reprocessed_local)

    # 2. Aggregate and Save
    aggregate_reprocessed(reprocessed_local, output_local)

    # 3. Optional: Cleanup the large original file to save disk space
    # os.remove(reprocessed_local)
