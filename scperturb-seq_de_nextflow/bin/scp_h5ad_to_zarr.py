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

# Ensure required obs fields
obs = A.obs.copy()
if "perturbation" not in obs.columns:
    for alt in ("perturbed_target_symbol",):
        if alt in obs.columns:
            obs["perturbation"] = obs[alt].astype(str)
            break
if "perturbation" not in obs.columns:
    raise SystemExit(
        "Missing obs['perturbation']; tried perturbed_target_symbol as fallbacks."
    )

obs["perturbation"] = obs["perturbation"].str.replace(
    r"^control.*", "control", regex=True
)

if "sample_id" not in obs.columns:
    obs["sample_id"] = "S1"

A.obs = obs

# FIX: adata.var index structure.
A.var.drop(columns=["gene_symbol"], inplace=True)
A.var.index.name = "gene_symbol"
A.var.reset_index(inplace=True)
A.var.set_index("gene_symbol", inplace=True)
A.var["gene_symbol"] = A.var.index

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

# In scp_preprocess.py (bigmode branch), before write_zarr:
try:
    A.raw = None
except Exception:
    pass
try:
    if hasattr(A, "layers") and A.layers is not None:
        # prevent writer from iterating layers
        A.layers = None
except Exception:
    pass
for attr in ("obsm", "varm", "obsp", "varp", "uns"):
    try:
        setattr(A, attr, {})
    except Exception:
        pass

print(f"[H5AD->ZARR] Writing Zarr -> {out_dir} (auto chunks)…")
# decide chunking
n_genes = int(A.n_vars)
try:
    probe = A.X[:1, :1]  # tiny read; ok even in backed mode
    is_sparse = hasattr(probe, "toarray")
except Exception:
    is_sparse = False

if is_sparse:
    # Sparse is stored as 1-D component arrays -> pass a single positive int
    A.write_zarr(str(out_dir), chunks=100_000)  # rows per chunk
else:
    # Dense 2-D: pass (rows_chunk, cols_chunk), both > 0
    A.write_zarr(str(out_dir), chunks=(100_000, min(n_genes, 50_000)))
