import os
import sys
import tempfile
import shutil
import gc
import uuid
import argparse
import numpy as np
import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import h5py
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set thread limits to 1 BEFORE importing numpy/scipy to prevent thread thrashing
# since we will be using process-level parallelism (multiprocessing)
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"


def log_mem(step):
    """Log current peak memory usage using resource module."""
    try:
        import resource

        # ru_maxrss is in KB on Linux
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
        print(f"[MEM] {step} - Peak RSS: {peak:.2f} GB")
    except Exception:
        pass


def get_sparse_x_slice(input_path, start, end):
    """Read a slice of the sparse matrix X directly from HDF5 to avoid scanpy overhead in workers."""
    with h5py.File(input_path, "r", libver="latest", swmr=True) as f:
        if "X" not in f:
            raise KeyError(f"Could not find 'X' in {input_path}")

        X = f["X"]
        if not isinstance(X, h5py.Group):
            # Fallback for dense X dataset
            total_obs = X.shape[0]
            actual_end = min(end, total_obs)
            if start >= actual_end:
                return None, None, None
            return X[start:actual_end, :], None, None

        # Sparse format (assuming CSR which is standard for anndata)
        indptr = X["indptr"]
        indices = X["indices"]
        data = X["data"]

        total_obs = len(indptr) - 1
        actual_end = min(end, total_obs)
        if start >= actual_end:
            return None, None, None

        # Read indptr for the slice
        s_indptr = indptr[start : actual_end + 1]
        nnz_start = s_indptr[0]
        nnz_end = s_indptr[-1]

        # Read data and indices for the specific range of non-zeros
        s_indices = indices[nnz_start:nnz_end]
        s_data = data[nnz_start:nnz_end]

        # Re-base indptr for the chunk (offsets within the chunk)
        s_indptr = (s_indptr - nnz_start).astype(np.int32)

        return s_data, s_indices, s_indptr


def process_chunk_batch(args):
    """Worker function for Phase 1: Reads chunks and partitions them into bin files."""
    (
        input_path,
        start_indices,
        chunk_size,
        group_indices_path,
        barcodes_per_bin,
        n_bins,
        temp_dir,
        n_vars,
    ) = args
    batch_id = uuid.uuid4().hex

    # Load shared group_indices array using memory-mapping to save memory
    group_indices = np.load(group_indices_path, mmap_mode="r")

    # Open all bin files for this batch
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

    try:
        for start in start_indices:
            end = start + chunk_size
            s_data, s_indices, s_indptr = get_sparse_x_slice(input_path, start, end)

            if s_data is None:
                continue

            if s_indptr is not None:
                # Sparse case: convert to COO-like arrays for binning
                n_rows = len(s_indptr) - 1
                r = np.repeat(np.arange(n_rows, dtype=np.int32), np.diff(s_indptr))
                c = s_indices.astype(np.int32)
                v = s_data.astype(np.float32)
            else:
                # Dense case
                X_chunk = s_data.astype(np.float32)
                r, c = X_chunk.nonzero()
                r = r.astype(np.int32)
                c = c.astype(np.int32)
                v = X_chunk[r, c]

            if len(r) == 0:
                continue

            # Map row to global unique barcode index
            t = group_indices[start + r]

            # Determine bin for each element
            bin_idx = t // barcodes_per_bin

            # Sort by bin_idx to group writes for efficiency
            sort_idx = np.argsort(bin_idx)
            bin_idx_sorted = bin_idx[sort_idx]
            t_sorted = t[sort_idx]
            c_sorted = c[sort_idx]
            v_sorted = v[sort_idx]

            # Find boundaries of each bin in the sorted arrays
            unique_bins, bin_starts = np.unique(bin_idx_sorted, return_index=True)
            bin_ends = np.append(bin_starts[1:], len(bin_idx_sorted))

            for ub, b_start, b_end in zip(unique_bins, bin_starts, bin_ends):
                bin_files_t[ub].write(t_sorted[b_start:b_end].tobytes())
                bin_files_c[ub].write(c_sorted[b_start:b_end].tobytes())
                bin_files_v[ub].write(v_sorted[b_start:b_end].tobytes())

    finally:
        for f in bin_files_t + bin_files_c + bin_files_v:
            f.close()

    return batch_id


def aggregate_bin(args):
    """Worker function for Phase 2: Aggregates binned data into CSR matrices."""
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
        csr = sp.csr_matrix((num_barcodes, n_vars), dtype=np.float32)
    else:
        t_bin = np.concatenate(t_list)
        c_bin = np.concatenate(c_list)
        v_bin = np.concatenate(v_list)

        # Sum duplicates (aggregation) within this bin
        coo = sp.coo_matrix((v_bin, (t_bin, c_bin)), shape=(num_barcodes, n_vars))
        coo.sum_duplicates()
        csr = coo.tocsr()

    out_path = os.path.join(temp_dir, f"bin_{bin_idx}_csr.npz")
    sp.save_npz(out_path, csr)

    # Clean up intermediate files
    for batch_id in batch_ids:
        for suffix in ["t", "c", "v"]:
            p = os.path.join(temp_dir, f"bin_{bin_idx}_{batch_id}_{suffix}.dat")
            if os.path.exists(p):
                os.remove(p)

    return bin_idx


def aggregate_reprocessed(input_path, output_path, n_cpus=None, row_chunk_size=20000):
    if n_cpus is None:
        try:
            n_cpus = len(os.sched_getaffinity(0))
        except AttributeError:
            n_cpus = os.cpu_count() or 1

    print(f"Starting aggregation with {n_cpus} CPUs...")
    log_mem("Start")

    # 1. Prepare aggregation map
    print("Reading metadata and preparing aggregation map...")
    with sc.read_h5ad(input_path, backed="r") as adata:
        n_obs, n_vars = adata.shape
        print(f"Dataset dimensions: {n_obs} cells x {n_vars} genes")

        # Read index in a memory-efficient way
        barcodes_all = np.array(adata.obs_names, dtype=str)
        log_mem("Loaded barcodes")

        # Split barcodes (e.g., 'ATGC-1' -> 'ATGC')
        barcodes_all = np.array(
            [b.split("-", 1)[0] for b in tqdm(barcodes_all, desc="Splitting barcodes")],
            dtype=str,
        )

        unique_barcodes, first_indices, group_indices = np.unique(
            barcodes_all, return_index=True, return_inverse=True
        )
        n_unique = len(unique_barcodes)
        print(f"Unique barcodes: {n_unique} (Reduction factor: {n_obs/n_unique:.2f}x)")

        # Selective metadata extraction using backed iloc
        obs_summed = adata.obs.iloc[first_indices].copy()
        var = adata.var.copy()

        del barcodes_all
        gc.collect()
        log_mem("Metadata mapping complete")

    # Setup temporary storage
    n_bins = max(1, n_unique // 50000)
    barcodes_per_bin = int(np.ceil(n_unique / n_bins))
    temp_dir = tempfile.mkdtemp(prefix="agg_temp_dir_", dir=os.getcwd())
    print(f"Using {n_bins} temporary bins in {temp_dir}")

    group_indices_path = os.path.join(temp_dir, "group_indices.npy")
    np.save(group_indices_path, group_indices.astype(np.int32))
    del group_indices
    gc.collect()

    try:
        # 2. Phase 1: Partitioning
        chunks_per_batch = 50
        chunk_starts = list(range(0, n_obs, row_chunk_size))
        batches = [
            chunk_starts[i : i + chunks_per_batch]
            for i in range(0, len(chunk_starts), chunks_per_batch)
        ]

        args_list = [
            (
                input_path,
                b,
                row_chunk_size,
                group_indices_path,
                barcodes_per_bin,
                n_bins,
                temp_dir,
                n_vars,
            )
            for b in batches
        ]

        print(f"Phase 1: Partitioning {n_obs} cells across {len(batches)} batches...")
        batch_ids = []
        with ProcessPoolExecutor(max_workers=n_cpus) as executor:
            futures = {
                executor.submit(process_chunk_batch, arg): arg for arg in args_list
            }
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Partitioning"
            ):
                batch_ids.append(future.result())

        log_mem("Phase 1 complete")

        # 3. Phase 2: Bin Aggregation
        print("Phase 2: Aggregating bins...")
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
                future.result()

        log_mem("Phase 2 complete")

        # 4. Phase 4: Final Matrix Assembly
        print("Phase 4: Merging final matrix...")
        total_nnz = 0
        for i in range(n_bins):
            with np.load(os.path.join(temp_dir, f"bin_{i}_csr.npz")) as data:
                total_nnz += len(data["data"])

        print(f"Final aggregated NNZ: {total_nnz}")

        final_data = np.empty(total_nnz, dtype=np.float32)
        final_indices = np.empty(total_nnz, dtype=np.int32)
        final_indptr = np.empty(n_unique + 1, dtype=np.int64)
        final_indptr[0] = 0

        curr_nnz = 0
        curr_row = 0
        for i in tqdm(range(n_bins), desc="Merging into final arrays"):
            out_path = os.path.join(temp_dir, f"bin_{i}_csr.npz")
            csr = sp.load_npz(out_path)

            n_rows_bin = csr.shape[0]
            nnz_bin = csr.nnz

            final_data[curr_nnz : curr_nnz + nnz_bin] = csr.data
            final_indices[curr_nnz : curr_nnz + nnz_bin] = csr.indices
            final_indptr[curr_row + 1 : curr_row + 1 + n_rows_bin] = (
                csr.indptr[1:] + curr_nnz
            )

            curr_nnz += nnz_bin
            curr_row += n_rows_bin

            del csr
            os.remove(out_path)

        summed_X_sparse = sp.csr_matrix(
            (final_data, final_indices, final_indptr), shape=(n_unique, n_vars)
        )
        del final_data, final_indices, final_indptr
        gc.collect()

        print("Constructing final AnnData object...")
        adata_sum = sc.AnnData(X=summed_X_sparse, obs=obs_summed, var=var)
        adata_sum.obs_names = unique_barcodes

        print(f"Saving aggregated data to {output_path}...")
        adata_sum.write(output_path)
        print("Pipeline finished successfully!")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Highly scalable barcode aggregation for Perturb-seq."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to input h5ad file")
    parser.add_argument(
        "-o", "--output", required=True, help="Path to output h5ad file"
    )
    parser.add_argument(
        "-c", "--cpus", type=int, default=None, help="Number of CPUs to use"
    )
    parser.add_argument(
        "--chunk-size", type=int, default=20000, help="Rows per read chunk"
    )
    args = parser.parse_args()

    # Respect Slurm/cgroups via os.sched_getaffinity if available
    cpus = args.cpus
    if cpus is None:
        try:
            cpus = len(os.sched_getaffinity(0))
        except AttributeError:
            cpus = os.cpu_count() or 1

    aggregate_reprocessed(
        args.input, args.output, n_cpus=cpus, row_chunk_size=args.chunk_size
    )
