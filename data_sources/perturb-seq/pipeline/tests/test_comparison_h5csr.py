#!/usr/bin/env python3
"""Focused parity and read-count check for comparison.H5CSRMatrix."""

import importlib.util
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

PIPELINE = Path(__file__).resolve().parents[1]
COMPARISON = PIPELINE / "comparison" / "comparison.py"


class ReadCountingArray:
    def __init__(self, dataset):
        self.dataset = dataset
        self.reads = 0
        self.dtype = dataset.dtype

    def __getitem__(self, key):
        self.reads += 1
        return self.dataset[key]


class Group:
    def __init__(self, group):
        self.attrs = group.attrs
        self.values = {
            name: ReadCountingArray(group[name])
            for name in ("data", "indices", "indptr")
        }

    def __getitem__(self, key):
        return self.values[key]


def write_group(parent, name, matrix):
    group = parent.create_group(name)
    csr = matrix.tocsr()
    group.attrs["shape"] = np.asarray(csr.shape, dtype=np.int64)
    group.create_dataset("data", data=csr.data)
    group.create_dataset("indices", data=csr.indices.astype(np.int64))
    group.create_dataset("indptr", data=csr.indptr.astype(np.int64))
    return Group(group)


def write_wide_group(parent, name, offset):
    group = parent.create_group(name)
    group.attrs["shape"] = np.asarray([3, 4], dtype=np.int64)
    data = group.create_dataset(
        "data", shape=(offset + 2,), dtype=np.float32, chunks=(2,), fillvalue=0
    )
    indices = group.create_dataset(
        "indices", shape=(offset + 2,), dtype=np.int64, chunks=(2,), fillvalue=0
    )
    data[offset : offset + 2] = [7, 9]
    indices[offset : offset + 2] = [1, 3]
    group.create_dataset(
        "indptr", data=np.asarray([0, offset, offset, offset + 2], dtype=np.int64)
    )
    return Group(group)


def load_comparison():
    old_argv = sys.argv
    sys.argv = [
        str(COMPARISON),
        "--dataset-id",
        "h5csr_test",
        "--curated-h5ad",
        "unused.h5ad",
        "--reprocessed-h5ad",
        "unused.h5ad",
        "--gtf",
        "unused.gtf",
    ]
    try:
        spec = importlib.util.spec_from_file_location(
            "comparison_h5csr_test", COMPARISON
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.argv = old_argv


def check(matrix, rows, cols, root, name, h5csr_matrix):
    requested_cols = slice(None) if cols is None else cols
    with h5py.File(Path(root) / f"{name}.h5", "w") as handle:
        group = write_group(handle, "matrix", matrix)
        result = h5csr_matrix(group)[rows, requested_cols]
        expected = matrix[np.asarray(rows, dtype=np.int64), :]
        if cols is not None:
            expected = expected[:, cols]
        np.testing.assert_array_equal(result.toarray(), expected.toarray())
        assert result.shape == expected.shape
        assert result.indptr.dtype == np.int64
        assert result.indices.dtype == np.int64
        reads = {key: dataset.reads for key, dataset in group.values.items()}
    return reads


def main():
    h5csr_matrix = load_comparison().H5CSRMatrix
    dense = np.zeros((12, 8), dtype=np.float32)
    dense[0, [1, 3]] = [7, 8]
    dense[4, [1, 6]] = [2, 5]
    dense[8, [0, 7]] = [3, 4]
    matrix = sparse.csr_matrix(dense)

    with tempfile.TemporaryDirectory(prefix="h5csr-parity-") as root:
        fast_reads = check(
            matrix,
            [4, 3, 4, 5],
            np.asarray([6, 1, 4, 0]),
            root,
            "fast",
            h5csr_matrix,
        )
        fallback_reads = check(matrix, [0, 4, 8], None, root, "scattered", h5csr_matrix)
        assert sum(fast_reads.values()) == 3, fast_reads
        assert sum(fallback_reads.values()) == 12, fallback_reads
        assert sum(fast_reads.values()) < sum(fallback_reads.values())

        dense_ratio = np.zeros((8, 8), dtype=np.float32)
        dense_ratio[2, 1] = 1
        dense_ratio[3, :] = np.arange(1, 9)
        dense_ratio[4, 2] = 1
        ratio_reads = check(
            sparse.csr_matrix(dense_ratio), [2, 4], None, root, "dense", h5csr_matrix
        )
        assert ratio_reads == {"data": 2, "indices": 2, "indptr": 5}, ratio_reads

        with h5py.File(Path(root) / "duplicate-cols.h5", "w") as handle:
            duplicate_cols = write_group(handle, "matrix", matrix)
            legacy = h5csr_matrix(duplicate_cols)[[4], [6, 6]]
            np.testing.assert_array_equal(
                legacy.toarray(), np.asarray([[0, 5]], dtype=np.float32)
            )
            assert duplicate_cols.values["data"].reads == 1
            assert duplicate_cols.values["indices"].reads == 1
            assert duplicate_cols.values["indptr"].reads == 2

        wide_path = Path(root) / "wide-offset.h5"
        offset = (1 << 31) + 17
        with h5py.File(wide_path, "w") as handle:
            wide = write_wide_group(handle, "matrix", offset)
            result = h5csr_matrix(wide)[[2], :]
            expected = sparse.csr_matrix(
                (
                    np.asarray([7, 9], dtype=np.float32),
                    np.asarray([1, 3], dtype=np.int64),
                    np.asarray([0, 2], dtype=np.int64),
                ),
                shape=(1, 4),
            )
            np.testing.assert_array_equal(result.toarray(), expected.toarray())
            assert result.indptr.tolist() == [0, 2]
            assert result.indptr.dtype == np.int64
            assert result.indices.dtype == np.int64
            reads_before_empty = {
                key: dataset.reads for key, dataset in wide.values.items()
            }
            assert sum(reads_before_empty.values()) == 3, reads_before_empty
            assert wide.values["data"].dataset.shape == (offset + 2,)
            empty = h5csr_matrix(wide)[np.asarray([], dtype=np.int64), :]
            assert empty.shape == (0, 4)
            assert empty.indptr.dtype == np.int64
            assert empty.indices.dtype == np.int64
            assert reads_before_empty == {
                key: dataset.reads for key, dataset in wide.values.items()
            }
        assert wide_path.stat().st_size < 1_000_000

    print(
        "H5CSR parity passed; fast reads="
        f"{sum(fast_reads.values())}, scattered fallback reads="
        f"{sum(fallback_reads.values())}, dense-span fallback reads="
        f"{sum(ratio_reads.values())}, logical offset={offset} verified"
    )


if __name__ == "__main__":
    main()
