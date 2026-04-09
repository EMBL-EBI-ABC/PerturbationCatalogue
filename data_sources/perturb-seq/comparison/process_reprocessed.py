import os
import sys
import tempfile
import shutil
import gc
import uuid
import numpy as np
import scipy.sparse as sp
import pandas as pd
import scanpy as sc
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set thread limits to 1 BEFORE importing numpy/scipy to prevent thread thrashing
# since we will be using process-level parallelism (multiprocessing)
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

n_cpus = os.cpu_count() or 1


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


def process_chunk_batch(args):
    """Worker function for Phase 1: Reads chunks of data and partitions them into bin files."""
    (
        input_path,
        start_indices,
        chunk_size,
        group_indices_path,
        barcodes_per_bin,
        n_bins,
        temp_dir,
    ) = args
    batch_id = uuid.uuid4().hex

    # Each worker opens the file locally
    adata = sc.read_h5ad(input_path, backed="r")

    # Load shared group_indices array using memory-mapping to save memory
    group_indices = np.load(group_indices_path, mmap_mode="r")

    bin_files_t = [
        open(os.path.join(temp_dir, f"bin_{i}_{batch_id}_t.dat"), "wb")
        for i in range(n_bins)
    ]
    bin_files_c = [
        open(os.path.join(temp_dir, f"bin_{i}_{batch_id}_c.dat"), "wb")
        for i in range(n_bins)
    ]
    bin_files_v = [
        open(os.path.join(temp_dir, f"bin_{i}_{batch_id}_v.dat"), "wb")
        for i in range(n_bins)
    ]

    n_obs = adata.n_obs
    for start in start_indices:
        end = min(start + chunk_size, n_obs)
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

    for f in bin_files_t:
        f.close()
    for f in bin_files_c:
        f.close()
    for f in bin_files_v:
        f.close()

    return batch_id


def aggregate_bin(args):
    """Worker function for Phase 2: Reads all batch files for a specific bin and aggregates them."""
    bin_idx, temp_dir, batch_ids, barcodes_per_bin, n_vars, n_unique = args

    bin_start_barcode = bin_idx * barcodes_per_bin
    bin_end_barcode = min((bin_idx + 1) * barcodes_per_bin, n_unique)
    num_barcodes = bin_end_barcode - bin_start_barcode

    t_list, c_list, v_list = [], [], []
    for batch_id in batch_ids:
        path_t = os.path.join(temp_dir, f"bin_{bin_idx}_{batch_id}_t.dat")
        path_c = os.path.join(temp_dir, f"bin_{bin_idx}_{batch_id}_c.dat")
        path_v = os.path.join(temp_dir, f"bin_{bin_idx}_{batch_id}_v.dat")

        if os.path.exists(path_t) and os.path.getsize(path_t) > 0:
            t_list.append(np.fromfile(path_t, dtype=np.int32) - bin_start_barcode)
            c_list.append(np.fromfile(path_c, dtype=np.int32))
            v_list.append(np.fromfile(path_v, dtype=np.float32))

    if not t_list:
        # Empty bin
        csr = sp.csr_matrix((num_barcodes, n_vars), dtype=np.float32)
    else:
        t_bin = np.concatenate(t_list)
        c_bin = np.concatenate(c_list)
        v_bin = np.concatenate(v_list)

        # Create COO and sum duplicates within this bin
        coo = sp.coo_matrix((v_bin, (t_bin, c_bin)), shape=(num_barcodes, n_vars))
        coo.sum_duplicates()
        csr = coo.tocsr()

    out_path = os.path.join(temp_dir, f"bin_{bin_idx}_csr.npz")
    sp.save_npz(out_path, csr)

    # Delete raw dat files to save space immediately
    for batch_id in batch_ids:
        for suffix in ["t", "c", "v"]:
            p = os.path.join(temp_dir, f"bin_{bin_idx}_{batch_id}_{suffix}.dat")
            if os.path.exists(p):
                os.remove(p)

    return bin_idx


def aggregate_reprocessed(input_path, output_path):
    """
    Highly-parallel, memory-efficient aggregation of reprocessed counts.
    Uses backed mode, ProcessPoolExecutor, and out-of-core disk partitioning.
    """
    print(f"Opening {input_path} in backed mode...")
    adata = sc.read_h5ad(input_path, backed="r")

    n_obs = adata.n_obs
    n_vars = adata.n_vars
    print(f"Dataset dimensions: {n_obs} cells x {n_vars} genes")

    # 1. Prepare aggregation map
    print("Reading metadata and preparing aggregation map...")
    obs = adata.obs.copy()
    obs["barcode"] = obs.index.str.split("-").str[0]
    unique_barcodes, group_indices = np.unique(obs["barcode"], return_inverse=True)
    n_unique = len(unique_barcodes)
    print(f"Unique barcodes: {n_unique} (Reduction factor: {n_obs/n_unique:.2f}x)")

    # Target ~50,000 barcodes per bin
    n_bins = max(1, n_unique // 50000)
    barcodes_per_bin = int(np.ceil(n_unique / n_bins))

    # Use current directory to avoid /tmp tmpfs RAM limits
    temp_dir = tempfile.mkdtemp(prefix="agg_temp_dir_", dir=os.getcwd())
    print(
        f"Using {n_bins} temporary bins in {temp_dir} for parallel out-of-core aggregation..."
    )

    # Save group_indices to mmap file for fast worker sharing (zero IPC overhead)
    group_indices = group_indices.astype(np.int32)
    group_indices_path = os.path.join(temp_dir, "group_indices.npy")
    np.save(group_indices_path, group_indices)

    try:
        # 2. Phase 1: Partitioning using multiple processes
        row_chunk_size = 20000
        chunks_per_batch = 50  # Each batch processes ~1M cells
        chunk_starts = list(range(0, n_obs, row_chunk_size))
        batches = [
            chunk_starts[i : i + chunks_per_batch]
            for i in range(0, len(chunk_starts), chunks_per_batch)
        ]

        args_list = [
            (
                input_path,
                batch,
                row_chunk_size,
                group_indices_path,
                barcodes_per_bin,
                n_bins,
                temp_dir,
            )
            for batch in batches
        ]

        print(
            f"Phase 1: Partitioning {n_obs} cells across {len(batches)} batches using {n_cpus} CPUs..."
        )
        batch_ids = []
        with ProcessPoolExecutor(max_workers=n_cpus) as executor:
            futures = {
                executor.submit(process_chunk_batch, arg): arg for arg in args_list
            }
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Partitioning"
            ):
                batch_ids.append(future.result())

        print("Phase 1 complete. Phase 2: Aggregating bins...")

        # 3. Phase 2: Bin Aggregation using multiple processes
        agg_args_list = [
            (i, temp_dir, batch_ids, barcodes_per_bin, n_vars, n_unique)
            for i in range(n_bins)
        ]

        with ProcessPoolExecutor(max_workers=n_cpus) as executor:
            futures = {
                executor.submit(aggregate_bin, arg): arg for arg in agg_args_list
            }
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Aggregating bins"
            ):
                _ = future.result()

        print("Phase 2 complete. Constructing final sparse matrix...")

        # 4. Final Matrix Assembly
        all_data = []
        all_indices = []
        all_indptr = [np.array([0], dtype=np.int64)]
        current_nnz = 0

        for i in tqdm(range(n_bins), desc="Merging final matrix"):
            out_path = os.path.join(temp_dir, f"bin_{i}_csr.npz")
            csr = sp.load_npz(out_path)
            all_data.append(csr.data)
            all_indices.append(csr.indices)
            all_indptr.append(csr.indptr[1:] + current_nnz)
            current_nnz += csr.nnz

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
        obs_summed = obs.groupby("barcode").first()
        adata_sum = sc.AnnData(X=summed_X_sparse, obs=obs_summed, var=adata.var.copy())
        adata_sum.obs_names = unique_barcodes

        print(f"Saving aggregated data to {output_path}...")
        adata_sum.write(output_path)
        print("Pipeline finished successfully!")

    finally:
        # Clean up temp files
        shutil.rmtree(temp_dir, ignore_errors=True)


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
