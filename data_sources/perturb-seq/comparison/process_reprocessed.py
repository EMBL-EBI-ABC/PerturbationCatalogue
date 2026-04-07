import os
import sys
import tempfile
import shutil
import gc
import numpy as np
import scipy.sparse as sp
import pandas as pd
import scanpy as sc
from tqdm import tqdm

# Set thread limits BEFORE importing numpy/scipy to ensure they are respected
n_cpus = os.cpu_count() or 1
os.environ["OMP_NUM_THREADS"] = str(n_cpus)
os.environ["MKL_NUM_THREADS"] = str(n_cpus)
os.environ["OPENBLAS_NUM_THREADS"] = str(n_cpus)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_cpus)
os.environ["NUMEXPR_NUM_THREADS"] = str(n_cpus)


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
    Memory-efficient aggregation of reprocessed counts.
    Uses backed mode and out-of-core disk partitioning to handle large datasets on limited RAM.
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

    # 2. Partition target barcodes into bins
    # Target ~50,000 barcodes per bin to keep memory usage very low during aggregation
    n_bins = max(1, n_unique // 50000)
    barcodes_per_bin = int(np.ceil(n_unique / n_bins))

    # Use current directory to avoid /tmp tmpfs RAM limits
    temp_dir = tempfile.mkdtemp(prefix="agg_temp_dir_", dir=os.getcwd())
    print(f"Using {n_bins} temporary bins in {temp_dir} for out-of-core aggregation...")

    try:
        bin_files_t = [
            open(os.path.join(temp_dir, f"bin_{i}_t.dat"), "wb") for i in range(n_bins)
        ]
        bin_files_c = [
            open(os.path.join(temp_dir, f"bin_{i}_c.dat"), "wb") for i in range(n_bins)
        ]
        bin_files_v = [
            open(os.path.join(temp_dir, f"bin_{i}_v.dat"), "wb") for i in range(n_bins)
        ]

        # 3. Chunked Processing & Partitioning
        row_chunk_size = 20000
        print(f"Processing {n_obs} cells in chunks of {row_chunk_size}...")

        group_indices = group_indices.astype(np.int32)

        for start in tqdm(range(0, n_obs, row_chunk_size)):
            end = min(start + row_chunk_size, n_obs)

            # Sequential read of the row chunk
            X_chunk = adata.X[start:end, :]

            if sp.issparse(X_chunk):
                X_chunk = X_chunk.tocoo()
                r = X_chunk.row.astype(np.int32)
                c = X_chunk.col.astype(np.int32)
                v = X_chunk.data.astype(np.float32)
            else:
                X_chunk = np.array(X_chunk, dtype=np.float32)
                r, c = X_chunk.nonzero()
                r = r.astype(np.int32)
                c = c.astype(np.int32)
                v = X_chunk[r, c]

            if len(r) == 0:
                continue

            # Map row to target barcode
            t = group_indices[start + r]

            # Determine bin for each element
            bin_idx = t // barcodes_per_bin

            # Sort by bin_idx to group writes
            sort_idx = np.argsort(bin_idx)
            bin_idx_sorted = bin_idx[sort_idx]
            t_sorted = t[sort_idx]
            c_sorted = c[sort_idx]
            v_sorted = v[sort_idx]

            # Find boundaries
            unique_bins, bin_starts = np.unique(bin_idx_sorted, return_index=True)
            bin_ends = np.append(bin_starts[1:], len(bin_idx_sorted))

            for ub, b_start, b_end in zip(unique_bins, bin_starts, bin_ends):
                bin_files_t[ub].write(t_sorted[b_start:b_end].tobytes())
                bin_files_c[ub].write(c_sorted[b_start:b_end].tobytes())
                bin_files_v[ub].write(v_sorted[b_start:b_end].tobytes())

        # Close all temp files
        for f_list in (bin_files_t, bin_files_c, bin_files_v):
            for f in f_list:
                f.close()

        print("Finished writing to temporary bins. Aggregating bins...")

        # 4. Bin Aggregation
        all_data = []
        all_indices = []
        all_indptr = [np.array([0], dtype=np.int64)]
        current_nnz = 0

        for i in tqdm(range(n_bins), desc="Aggregating bins"):
            path_t = os.path.join(temp_dir, f"bin_{i}_t.dat")
            path_c = os.path.join(temp_dir, f"bin_{i}_c.dat")
            path_v = os.path.join(temp_dir, f"bin_{i}_v.dat")

            bin_start_barcode = i * barcodes_per_bin
            bin_end_barcode = min((i + 1) * barcodes_per_bin, n_unique)
            num_barcodes = bin_end_barcode - bin_start_barcode

            if os.path.getsize(path_t) == 0:
                # Empty bin
                all_indptr.append(np.full(num_barcodes, current_nnz, dtype=np.int64))
                continue

            # Localize target indices for this bin
            t_bin = np.fromfile(path_t, dtype=np.int32) - bin_start_barcode
            c_bin = np.fromfile(path_c, dtype=np.int32)
            v_bin = np.fromfile(path_v, dtype=np.float32)

            # Create COO and sum duplicates within this bin
            coo = sp.coo_matrix((v_bin, (t_bin, c_bin)), shape=(num_barcodes, n_vars))
            coo.sum_duplicates()
            csr = coo.tocsr()

            all_data.append(csr.data)
            all_indices.append(csr.indices)
            all_indptr.append(csr.indptr[1:] + current_nnz)

            current_nnz += len(csr.data)

            # Free memory early
            del t_bin, c_bin, v_bin, coo, csr
            gc.collect()

    finally:
        # 5. Clean up temp files
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("Constructing final sparse matrix...")
    if all_data:
        final_data = np.concatenate(all_data)
        final_indices = np.concatenate(all_indices)
    else:
        final_data = np.array([], dtype=np.float32)
        final_indices = np.array([], dtype=np.int32)

    final_indptr = np.concatenate(all_indptr)

    del all_data, all_indices, all_indptr
    gc.collect()

    summed_X_sparse = sp.csr_matrix(
        (final_data, final_indices, final_indptr), shape=(n_unique, n_vars)
    )

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
