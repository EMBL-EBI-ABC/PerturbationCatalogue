#!/usr/bin/env python3
"""
Convert a .h5ad to a .zarr store for out-of-core processing.

- Works with dense or sparse X.
- For sparse X, avoids 2-D chunks (uses auto or 1-D chunks).
- For huge files, prefer --counts-layer X to avoid materialising data.

Usage:
  scp_h5ad_to_zarr.py --input data.h5ad --counts-layer X --out data.zarr
"""

import argparse, sys
from pathlib import Path
import anndata as ad

p = argparse.ArgumentParser()
p.add_argument("--input", required=True, help="Path to input .h5ad")
p.add_argument(
    "--counts-layer",
    default="counts",
    help="Layer to use as .X (use 'X' to keep existing X)",
)
p.add_argument("--out", required=True, help="Output path (.zarr dir)")
# Optional: if you *know* X is dense and want to suggest chunks; ignored for sparse
p.add_argument(
    "--genes-chunk", type=int, default=2000, help="Hint for dense X; ignored for sparse"
)
p.add_argument(
    "--cells-chunk",
    type=int,
    default=20000,
    help="Hint for dense X; ignored for sparse",
)
p.add_argument(
    "--allow-full-load",
    action="store_true",
    help="If counts-layer != X, allow a full in-memory load to swap layers (NOT for 500GB files).",
)
args = p.parse_args()

in_path = Path(args.input)
if not in_path.exists():
    sys.exit(f"[H5AD->ZARR] Input not found: {in_path}")

out_dir = Path(args.out)
if out_dir.suffix != ".zarr":
    out_dir = out_dir.with_suffix(".zarr")
out_dir.parent.mkdir(parents=True, exist_ok=True)

print(f"[H5AD->ZARR] Reading {in_path} (backed='r') …")
try:
    A = ad.read_h5ad(str(in_path), backed="r")
except Exception as e:
    sys.exit(f"[H5AD->ZARR] Failed to read {in_path}: {e}")

# Layer -> X swap
if args.counts_layer and args.counts_layer != "X":
    # Swapping layers into X requires a non-backed AnnData (may fully load!)
    if not args.allow_full_load:
        sys.exit(
            "[H5AD->ZARR] counts-layer != 'X' but --allow-full-load not set. "
            "For huge files, ensure raw counts are already in .X, or rerun with --counts-layer X."
        )
    print(
        f"[H5AD->ZARR] Loading in-memory to swap layer '{args.counts_layer}' into X (may be heavy)…"
    )
    try:
        A_mem = ad.read_h5ad(str(in_path))  # full load
    except Exception as e:
        sys.exit(f"[H5AD->ZARR] Could not load in-memory: {e}")
    if A_mem.layers is None or args.counts_layer not in A_mem.layers:
        print(
            f"[H5AD->ZARR] Warning: layer '{args.counts_layer}' not found; leaving X unchanged.",
            file=sys.stderr,
        )
    else:
        A_mem.X = A_mem.layers[args.counts_layer]
        A_mem.layers = None
    # From here, write Zarr from in-memory object
    print(f"[H5AD->ZARR] Writing Zarr -> {out_dir} (auto chunks)…")
    try:
        A_mem.write_zarr(str(out_dir))
    except Exception as e:
        sys.exit(f"[H5AD->ZARR] Failed to write Zarr: {e}")
    print(f"[H5AD->ZARR] Done.")
    sys.exit(0)

# For backed AnnData: try writing with auto-chunking first
print(f"[H5AD->ZARR] Writing Zarr -> {out_dir} (auto chunks)…")
try:
    A.write_zarr(
        str(out_dir)
    )  # Let anndata pick chunks appropriate to dense/sparse encoders
    print("[H5AD->ZARR] Done.")
    sys.exit(0)
except Exception as e_auto:
    print(f"[H5AD->ZARR] Auto chunking failed: {e_auto}", file=sys.stderr)

# Retry with explicit chunks:
# - For dense X: use (genes_chunk, cells_chunk)
# - For sparse X: use 1-D chunk (e.g., genes_chunk) because sparse stores 1-D component arrays
# We can heuristically detect sparse encoding by trying dense chunks first, then fallback to 1-D.

print(
    f"[H5AD->ZARR] Retrying with dense-style chunks=({args.genes_chunk},{args.cells_chunk}) …"
)
try:
    A.write_zarr(str(out_dir), chunks=(args.genes_chunk, args.cells_chunk))
    print("[H5AD->ZARR] Done (dense-style chunks).")
    sys.exit(0)
except Exception as e_dense:
    print(f"[H5AD->ZARR] Dense-style chunks failed: {e_dense}", file=sys.stderr)

print("[H5AD->ZARR] Retrying with 1-D chunks for sparse encodings …")
try:
    A.write_zarr(str(out_dir), chunks=args.genes_chunk)  # 1-D chunk length
    print("[H5AD->ZARR] Done (sparse-style chunks).")
    sys.exit(0)
except Exception as e_sparse:
    sys.exit(f"[H5AD->ZARR] Failed to write Zarr even with 1-D chunks: {e_sparse}")
