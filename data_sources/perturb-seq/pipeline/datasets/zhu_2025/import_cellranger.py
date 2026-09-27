#!/usr/bin/env python3
"""Convert one Zhu donor/condition from the GEO Cell Ranger archive to H5AD."""

import argparse
import csv
import gc
import re
import shutil
import tarfile
from pathlib import Path

import anndata as ad
import h5py
import scanpy as sc


MATRIX_SUFFIX = "_sample_filtered_feature_bc_matrix.h5"
GUIDE_ALIASES = {"1-Jun": "JUN-1", "2-Jun": "JUN-2"}


def lane_id(name):
    match = re.search(r"CD4i_(R[12]L\d{2})_", name)
    if not match:
        raise ValueError(f"Cannot parse lane from {name}")
    return match.group(1)


def expected_lanes(sample):
    return (
        23
        if sample.startswith(("D1_", "D2_")) and not sample.endswith("Stim48hr")
        else 24
    )


def selected_members(archive, sample, limit=0):
    suffix = f"_{sample}{MATRIX_SUFFIX}"
    members = sorted(
        (
            m
            for m in archive.getmembers()
            if m.isfile() and Path(m.name).name.endswith(suffix)
        ),
        key=lambda m: lane_id(Path(m.name).name),
    )
    if limit:
        members = members[:limit]
    if not members:
        raise ValueError(f"No Cell Ranger matrices found for {sample}")
    lanes = [lane_id(Path(m.name).name) for m in members]
    if len(lanes) != len(set(lanes)):
        raise ValueError(f"Duplicate lanes for {sample}: {lanes}")
    if not limit and len(members) != expected_lanes(sample):
        raise ValueError(
            f"Expected {expected_lanes(sample)} lanes for {sample}, found {len(members)}"
        )
    return members


def selected_paths(directory, sample, limit=0):
    paths = sorted(
        Path(directory).glob(f"*_{sample}{MATRIX_SUFFIX}"),
        key=lambda p: lane_id(p.name),
    )
    if limit:
        paths = paths[:limit]
    lanes = [lane_id(path.name) for path in paths]
    if len(lanes) != len(set(lanes)):
        raise ValueError(f"Duplicate lanes for {sample}: {lanes}")
    if not limit and len(paths) != expected_lanes(sample):
        raise ValueError(
            f"Expected {expected_lanes(sample)} lanes for {sample}, found {len(paths)}"
        )
    return paths


def guide_labels(path):
    labels = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            guide = row["guide_id"]
            labels[guide] = "_".join(
                x for x in (row["target_symbol"], guide, row["target_ensg"]) if x
            )
    return labels


def convert_lane(source, output, sample, lane, targets):
    combined = sc.read_10x_h5(source, gex_only=False)
    feature_types = combined.var["feature_types"].astype(str)
    gene_mask = feature_types == "Gene Expression"
    guide_mask = (
        feature_types == "CRISPR Guide Capture"
    ) & ~combined.var_names.str.startswith("ProbeNTC-")
    if not gene_mask.any() or not guide_mask.any():
        raise ValueError(f"Missing gene or guide features in {source}")

    raw_guides = combined.var_names[guide_mask].astype(str).tolist()
    guide_ids = [GUIDE_ALIASES.get(name, name) for name in raw_guides]
    missing = sorted(set(guide_ids) - targets.keys())
    if missing:
        raise ValueError(
            f"Guide metadata missing {len(missing)} features; first: {missing[:5]}"
        )

    genes = combined[:, gene_mask].copy()
    genes.obsm["guides"] = combined[:, guide_mask].X.copy().tocsr()
    genes.uns["guide_names"] = [targets[name] for name in guide_ids]
    genes.var["gene_name"] = genes.var_names.astype(str)
    genes.var_names = genes.var["gene_ids"].astype(str)
    genes.var_names_make_unique()
    genes.obs["sample_id"] = sample
    genes.obs["lane_id"] = lane
    genes.obs_names = [f"{barcode}_{lane}_{sample}" for barcode in genes.obs_names]
    genes.write_h5ad(output, compression="gzip")


def run(args):
    targets = guide_labels(args.guide_targets)
    lane_dir = Path("lane_h5ads")
    lane_dir.mkdir()
    inputs = {}
    n_cells = 0

    def process_source(source, lane):
        nonlocal n_cells
        output = lane_dir / f"{lane}.h5ad"
        convert_lane(source, output, args.sample, lane, targets)
        lane_data = ad.read_h5ad(output, backed="r")
        n_cells += lane_data.n_obs
        lane_data.file.close()
        inputs[lane] = str(output)
        gc.collect()

    if args.tar:
        with tarfile.open(args.tar) as archive:
            for member in selected_members(archive, args.sample, args.limit_lanes):
                lane = lane_id(Path(member.name).name)
                source = Path(f"{lane}.h5")
                with archive.extractfile(member) as incoming, source.open(
                    "wb"
                ) as outgoing:
                    shutil.copyfileobj(incoming, outgoing, length=16 * 1024 * 1024)
                process_source(source, lane)
                source.unlink()
    else:
        for source in selected_paths(args.matrix_dir, args.sample, args.limit_lanes):
            process_source(source, lane_id(source.name))

    ad.experimental.concat_on_disk(
        inputs,
        args.output,
        axis=0,
        join="outer",
        merge="same",
        uns_merge="same",
        index_unique="-",
    )
    with h5py.File(next(iter(inputs.values())), "r") as source, h5py.File(
        args.output, "a"
    ) as target:
        if "uns" in target:
            del target["uns"]
        source.copy("uns", target)

    result = ad.read_h5ad(args.output, backed="r")
    try:
        if result.n_obs != n_cells or result.obsm["guides"].shape[0] != n_cells:
            raise ValueError(
                "Concatenated cell or guide count does not match lane inputs"
            )
        if result.obsm["guides"].shape[1] != len(result.uns["guide_names"]):
            raise ValueError("Guide matrix and guide-name counts differ")
        print(
            f"sample={args.sample} lanes={len(inputs)} cells={result.n_obs} "
            f"genes={result.n_vars} guides={result.obsm['guides'].shape[1]}"
        )
    finally:
        result.file.close()


def self_test():
    class Member:
        def __init__(self, name):
            self.name = name

        def isfile(self):
            return True

    class Archive:
        def getmembers(self):
            return [
                Member("GSM2_CD4i_R1L02_D1_Rest" + MATRIX_SUFFIX),
                Member("GSM1_CD4i_R1L01_D1_Rest" + MATRIX_SUFFIX),
                Member("GSM3_CD4i_R1L01_D2_Rest" + MATRIX_SUFFIX),
            ]

    assert [
        lane_id(Path(m.name).name)
        for m in selected_members(Archive(), "D1_Rest", limit=2)
    ] == [
        "R1L01",
        "R1L02",
    ]
    assert expected_lanes("D1_Rest") == 23
    assert expected_lanes("D4_Stim48hr") == 24
    assert GUIDE_ALIASES["1-Jun"] == "JUN-1"


def main():
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--tar")
    source.add_argument("--matrix-dir")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--guide-targets", required=True)
    parser.add_argument("--output", default="experiment_final_uncompressed.h5ad")
    parser.add_argument("--limit-lanes", type=int, default=0)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if not args.tar and not args.matrix_dir:
            parser.error("one of --tar or --matrix-dir is required")
        run(args)


if __name__ == "__main__":
    main()
