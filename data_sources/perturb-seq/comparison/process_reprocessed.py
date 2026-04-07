import os
import sys

# Set thread limits BEFORE importing numpy/scipy to ensure they are respected
n_cpus = os.cpu_count() or 1
os.environ["OMP_NUM_THREADS"] = str(n_cpus)
os.environ["MKL_NUM_THREADS"] = str(n_cpus)
os.environ["OPENBLAS_NUM_THREADS"] = str(n_cpus)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_cpus)
os.environ["NUMEXPR_NUM_THREADS"] = str(n_cpus)

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gc
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor


def download_file(gs_path, local_path):
    """Downloads a file from Google Cloud Storage using gsutil."""
    if os.path.exists(local_path):
        print(f"File {local_path} already exists. Skipping download.")
        return

    print(f"Downloading {gs_path} to {local_path}...")
    # -m for multi-threaded/multi-processing copy
    ret = os.system(f"gsutil -m cp {gs_path} {local_path}")
    if ret != 0:
        print(f"Error: gsutil failed with exit code {ret}")
        sys.exit(1)


def aggregate_reprocessed(input_path, output_path):
    """
    Memory-efficient and multi-threaded aggregation of reprocessed counts.
    Uses backed mode and chunked processing to handle large datasets on limited RAM.
    """
    print(f"Opening {input_path} in backed mode...")
    # Backed mode 'r' avoids loading the matrix into RAM, only metadata is loaded.
    adata = sc.read_h5ad(input_path, backed="r")

    n_obs = adata.n_obs
    n_vars = adata.n_vars
    print(f"Dataset dimensions: {n_obs} cells x {n_vars} genes")

    # 1. Prepare aggregation map
    print("Reading metadata and preparing aggregation map...")
    obs = adata.obs.copy()
    # Extract barcode (part before '-')
    obs["barcode"] = obs.index.str.split("-").str[0]
    unique_barcodes, group_indices = np.unique(obs["barcode"], return_inverse=True)
    n_unique = len(unique_barcodes)
    print(f"Unique barcodes: {n_unique} (Reduction factor: {n_obs/n_unique:.2f}x)")

    # 2. Initialize dense accumulator
    # We use a dense float32 array for the fastest possible accumulation via BLAS.
    # 1M cells x 30k genes @ float32 = ~111 GB. 340k cells = ~38 GB.
    expected_gb = (n_unique * n_vars * 4) / (1024**3)
    print(f"Allocating {expected_gb:.2f} GB for the dense accumulator...")

    try:
        summed_X = np.zeros((n_unique, n_vars), dtype=np.float32)
    except MemoryError:
        print(f"CRITICAL ERROR: Failed to allocate {expected_gb:.2f} GB RAM.")
        sys.exit(1)

    # 3. Chunked Processing
    # Row-wise chunks are mandatory for efficient reading from CSR HDF5 files.
    # 20,000 cells * 30,000 genes * 4 bytes = ~2.4 GB per chunk.
    row_chunk_size = 20000
    print(
        f"Processing {n_obs} cells in chunks of {row_chunk_size} using {n_cpus} CPUs..."
    )

    # Helper for parallel addition across columns to utilize all cores
    def parallel_add(target_indices, data_chunk):
        col_step = (n_vars + n_cpus - 1) // n_cpus

        def worker(c_start):
            c_end = min(c_start + col_step, n_vars)
            # Advanced indexing with += is vectorized.
            # Since target_indices are unique (pre-processed), this is thread-safe.
            summed_X[target_indices, c_start:c_end] += data_chunk[:, c_start:c_end]

        with ThreadPoolExecutor(max_workers=n_cpus) as executor:
            list(executor.map(worker, range(0, n_vars, col_step)))

    for start in tqdm(range(0, n_obs, row_chunk_size)):
        end = min(start + row_chunk_size, n_obs)

        # Disk I/O: Sequential read of the row chunk
        X_chunk = adata.X[start:end, :]
        if sp.issparse(X_chunk):
            X_chunk = X_chunk.toarray().astype(np.float32)
        else:
            X_chunk = np.array(X_chunk, dtype=np.float32)

        chunk_targets = group_indices[start:end]

        # Intra-chunk aggregation: Reduce updates to the large matrix by summing
        # multiple runs for the same cell that happen to be in this chunk.
        u_in_chunk, inv_in_chunk = np.unique(chunk_targets, return_inverse=True)

        if len(u_in_chunk) < len(chunk_targets):
            # Efficiently sum rows within the chunk using a small sparse mapping matrix
            agg_small = sp.csr_matrix(
                (
                    np.ones(len(chunk_targets), dtype=np.float32),
                    (inv_in_chunk, np.arange(len(chunk_targets))),
                ),
                shape=(len(u_in_chunk), len(chunk_targets)),
            )
            # Sparse-Dense multiplication is highly optimized in SciPy/BLAS
            res_chunk = agg_small @ X_chunk
            parallel_add(u_in_chunk, res_chunk)
        else:
            # All cells in chunk belong to unique barcodes
            parallel_add(chunk_targets, X_chunk)

    print("Aggregation complete.")

    # 4. Conversion and Save
    print("Converting dense accumulator to sparse CSR format...")
    # This step is single-threaded in SciPy but memory efficient.
    summed_X_sparse = sp.csr_matrix(summed_X)
    del summed_X
    gc.collect()

    print("Constructing final AnnData object...")
    # Keep the first metadata entry for each barcode
    obs_summed = obs.groupby("barcode").first()
    adata_sum = sc.AnnData(X=summed_X_sparse, obs=obs_summed, var=adata.var.copy())
    adata_sum.obs_names = unique_barcodes

    print(f"Saving aggregated data to {output_path}...")
    adata_sum.write(output_path)
    print("Pipeline finished successfully!")


if __name__ == "__main__":
    lake_bucket = os.environ.get("LAKE_BUCKET")
    if not lake_bucket:
        print("Error: LAKE_BUCKET environment variable is not set.")
        sys.exit(1)

    gs_path = f"gs://{lake_bucket}/perturbseq/fastq-reprocess/nadig_2025_jurkat.h5ad"
    local_in = "nadig_2025_jurkat_reprocessed.h5ad"
    local_out = "nadig_2025_jurkat_reprocessed_summed.h5ad"

    download_file(gs_path, local_in)
    aggregate_reprocessed(local_in, local_out)
