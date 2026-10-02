#!/usr/bin/env python3
"""Select one donor/state from reusable lane-level Flex count H5ADs."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import anndata as ad
import h5py


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def aggregate(args):
    samples = [
        row for row in read_tsv(args.samples) if row["sample_id"] == args.sample_id
    ]
    if len(samples) != 1:
        raise ValueError(
            f"Expected one sample row for {args.sample_id}, found {len(samples)}"
        )
    sample = samples[0]
    if sample["pool_id"] != args.pool_id:
        raise ValueError(
            "Selected sample does not belong to the requested physical pool"
        )
    safe_id = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
    if not safe_id.fullmatch(args.sample_id) or not safe_id.fullmatch(args.pool_id):
        raise ValueError("Pool and sample IDs must be safe non-empty identifiers")
    pairs = [tuple(pair.split(":")) for pair in sample["bc_cr_pairs"].split(";")]
    safe_alias = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}")
    if (
        not pairs
        or any(
            len(pair) != 2 or not safe_alias.fullmatch(alias)
            for pair in pairs
            for alias in pair
        )
        or len({bc for bc, _ in pairs}) != len(pairs)
        or len({cr for _, cr in pairs}) != len(pairs)
    ):
        raise ValueError("Sample metadata must select unique safe BC/CR aliases")

    expected_lanes = sample["pool_lanes"].split(";")
    if (
        not expected_lanes
        or len(expected_lanes) != len(set(expected_lanes))
        or any(not safe_id.fullmatch(lane) for lane in expected_lanes)
    ):
        raise ValueError("Sample metadata must list unique safe lane IDs")
    coverage_rows = [
        row for row in read_tsv(args.guide_coverage) if row["pool_id"] == args.pool_id
    ]
    guide_coverage = {row["lane_id"]: row["guide_status"] for row in coverage_rows}
    if len(guide_coverage) != len(coverage_rows):
        raise ValueError("Guide coverage manifest has duplicate lane rows")
    if set(guide_coverage) != set(expected_lanes):
        raise ValueError("Guide coverage manifest and sample lane list differ")
    if any(
        status not in {"available", "partial", "no_archived_guide_sra"}
        for status in guide_coverage.values()
    ):
        raise ValueError("Guide coverage manifest has an unsupported status")

    barcode_rows = read_tsv(args.probe_barcodes)
    barcode_pairs = {row["bc_alias"]: row["cr_alias"] for row in barcode_rows}
    bc_sequences = {row["bc_alias"]: row["bc_sequence"] for row in barcode_rows}
    if (
        not barcode_rows
        or len(barcode_pairs) != len(barcode_rows)
        or len({row["cr_alias"] for row in barcode_rows}) != len(barcode_rows)
        or len(set(bc_sequences.values())) != len(barcode_rows)
        or any(
            not re.fullmatch(r"[ACGT]{8}", sequence)
            for sequence in bc_sequences.values()
        )
        or any(barcode_pairs.get(bc) != cr for bc, cr in pairs)
    ):
        raise ValueError(
            "Sample BC/CR assignment disagrees with the pinned probe-barcode table"
        )

    available = {}
    for path in map(Path, args.lane_counts):
        match = re.fullmatch(r"(.+)__(.+)\.h5ad", path.name)
        if not match:
            raise ValueError(f"Unexpected lane-count filename: {path.name}")
        key = (match.group(1), match.group(2))
        if key in available:
            raise ValueError(f"Duplicate lane-count product: {key}")
        available[key] = path.resolve(strict=True)

    inputs, cell_count = {}, 0
    first_product = None
    validated_products = 0
    all_aliases = set(barcode_pairs)
    all_expected = {(lane, bc) for lane in expected_lanes for bc in all_aliases}
    if available.keys() != all_expected:
        missing = all_expected - available.keys()
        unexpected = available.keys() - all_expected
        raise ValueError(
            f"Lane products differ from manifest: missing={len(missing)}, unexpected={len(unexpected)}"
        )
    expected = {(lane, bc) for lane in expected_lanes for bc, _ in pairs}
    if expected - available.keys():
        raise ValueError(
            f"Missing {len(expected - available.keys())} selected lane/BC count products"
        )
    expected_by_bc = dict(pairs)
    var_ids, var_frame, guide_names = None, None, None
    counter_version_label = None
    counter_version = None
    umi_policy = None
    mapping_policy = None
    reference_hashes_json = None
    lane_cells = {}
    lane_metrics = {}
    for lane, bc in sorted(expected):
        path = available[(lane, bc)]
        if first_product is None:
            first_product = path
        validated_products += 1
        data = ad.read_h5ad(path, backed="r")
        product_n_obs = None
        try:
            product_n_obs = data.n_obs
            if data.uns.get("guide_names") is None:
                raise ValueError(f"Missing guide_names in {path}")
            current_guides = [str(value) for value in data.uns["guide_names"]]
            current_var = [str(value) for value in data.var_names]
            current_var_frame = data.var
            current_counter_label = str(data.uns.get("flex_counter_version_label", ""))
            current_counter = dict(data.uns.get("flex_counter_version", {}))
            current_umi_policy = str(data.uns.get("flex_umi_policy", ""))
            current_mapping_policy = str(data.uns.get("flex_mapping_policy", ""))
            current_reference_hashes = str(
                data.uns.get("flex_reference_hashes_json", "")
            )
            lane_metrics_json = str(data.uns.get("flex_lane_metrics_json", ""))
            if not all(
                (
                    current_counter_label,
                    current_counter,
                    current_umi_policy,
                    current_mapping_policy,
                    current_reference_hashes,
                    lane_metrics_json,
                )
            ):
                raise ValueError(f"Missing native Flex provenance in {path}")
            try:
                alias_metrics = json.loads(lane_metrics_json)
                json.loads(current_reference_hashes)
            except json.JSONDecodeError as error:
                raise ValueError(f"Malformed Flex provenance JSON in {path}") from error
            if (
                alias_metrics.get("lane_id") != path.name.split("__", 1)[0]
                or alias_metrics.get("bc_alias") != path.name.rsplit("__", 1)[1][:-5]
                or alias_metrics.get("guide_coverage_status")
                != guide_coverage.get(alias_metrics.get("lane_id"))
            ):
                raise ValueError(f"Lane metrics identity is inconsistent in {path}")
            if data.obsm["guides"].shape != (data.n_obs, len(current_guides)):
                raise ValueError(f"Guide matrix shape is inconsistent in {path}")
            if guide_names is not None and current_guides != guide_names:
                raise ValueError(f"Guide feature order differs in {path}")
            if var_ids is not None and current_var != var_ids:
                raise ValueError(f"Gene feature order differs in {path}")
            if var_frame is not None and not current_var_frame.equals(var_frame):
                raise ValueError(f"Gene feature annotations differ in {path}")
            if counter_version_label is not None and (
                current_counter_label != counter_version_label
                or current_counter != counter_version
                or current_umi_policy != umi_policy
                or current_mapping_policy != mapping_policy
                or current_reference_hashes != reference_hashes_json
            ):
                raise ValueError(f"Counter/reference provenance differs in {path}")
            guide_names, var_ids = current_guides, current_var
            var_frame = current_var_frame.copy()
            counter_version_label = current_counter_label
            counter_version = current_counter
            umi_policy = current_umi_policy
            mapping_policy = current_mapping_policy
            reference_hashes_json = current_reference_hashes
            expected_uns = {
                "flex_sample_id": args.sample_id,
                "flex_lane_id": lane,
                "flex_bc_alias": bc,
                "flex_cr_alias": expected_by_bc.get(bc, barcode_pairs.get(bc, "")),
                "flex_guide_coverage_status": guide_coverage[lane],
            }
            if any(
                str(data.uns.get(key, "")) != value
                for key, value in expected_uns.items()
            ):
                raise ValueError(f"Lane product provenance is inconsistent in {path}")
            if data.n_obs:
                obs = data.obs
                if set(obs["sample_id"].astype(str)) != {args.sample_id}:
                    raise ValueError(f"Wrong sample IDs in {path}")
                if set(obs["lane_id"].astype(str)) != {lane}:
                    raise ValueError(f"Wrong lane IDs in {path}")
                if set(obs["bc_alias"].astype(str)) != {bc}:
                    raise ValueError(f"Wrong BC alias in {path}")
                if set(obs["cr_alias"].astype(str)) != {expected_by_bc[bc]}:
                    raise ValueError(f"Wrong CR alias in {path}")
                expected_status = guide_coverage[lane]
                if set(obs["guide_coverage_status"].astype(str)) != {expected_status}:
                    raise ValueError(f"Guide coverage status is wrong in {path}")
                if not data.obs_names.is_unique:
                    raise ValueError(f"Duplicate cell IDs in {path}")
                suffix = f"{bc_sequences[bc]}-1_{lane}"
                for cell_id in data.obs_names.astype(str):
                    if not cell_id.endswith(suffix) or not re.fullmatch(
                        r"[ACGT]{16}", cell_id[: -len(suffix)]
                    ):
                        raise ValueError(
                            f"Cell identity does not match BC sequence/lane in {path}"
                        )
            cell_count += product_n_obs
            lane_cells[(lane, bc)] = product_n_obs
            lane_entry = lane_metrics.setdefault(
                lane,
                {
                    "gex_mapping": alias_metrics["gex_mapping"],
                    "guide_mapping": alias_metrics["guide_mapping"],
                    "equivalence_class_metrics_scope": alias_metrics[
                        "equivalence_class_metrics_scope"
                    ],
                    "barcode_metrics": {},
                },
            )
            for metric_name in ("gex_mapping", "guide_mapping"):
                if lane_entry[metric_name] != alias_metrics[metric_name]:
                    raise ValueError(
                        f"Lane-level {metric_name} differs across aliases in {path}"
                    )
            if (
                lane_entry["equivalence_class_metrics_scope"]
                != alias_metrics["equivalence_class_metrics_scope"]
            ):
                raise ValueError(f"Equivalence-class metric scope differs in {path}")
            lane_entry["barcode_metrics"][bc] = alias_metrics["barcode_metrics"]
        finally:
            data.file.close()
        key = f"{lane}_{bc}"
        if key in inputs:
            raise ValueError(f"Duplicate aggregation key: {key}")
        if product_n_obs:
            inputs[key] = str(path)

    if var_ids is None or guide_names is None:
        raise ValueError("No lane count products were selected")
    if first_product is None:
        raise ValueError("No lane count products were selected")
    if inputs:
        ad.experimental.concat_on_disk(
            inputs,
            args.output,
            max_loaded_elems=args.max_loaded_elems,
            axis=0,
            join="inner",
            # The installed anndata experimental writer cannot serialize the
            # Series-valued result of merge="same" into /var. Feature order
            # was validated above, so copy the pinned var group afterward.
            merge=None,
            uns_merge="same",
            index_unique=None,
        )
        with h5py.File(first_product, "r") as source, h5py.File(
            args.output, "a"
        ) as target:
            if "var" in target:
                del target["var"]
            source.copy("var", target)
    else:
        empty = ad.read_h5ad(first_product, backed="r")
        try:
            obs = pd.DataFrame(
                {
                    name: pd.Series(dtype=object)
                    for name in (
                        "sample_id",
                        "lane_id",
                        "bc_alias",
                        "cr_alias",
                        "guide_coverage_status",
                    )
                },
                index=pd.Index([], dtype=object),
            )
            empty_data = ad.AnnData(
                X=sparse.csr_matrix((0, len(var_ids)), dtype="int32"),
                obs=obs,
                var=var_frame.copy(),
            )
            empty_data.obsm["guides"] = sparse.csr_matrix(
                (0, len(guide_names)), dtype="int32"
            )
            empty_data.uns["guide_names"] = guide_names
            empty_data.write_h5ad(args.output, compression="gzip")
        finally:
            empty.file.close()
    coverage_records = [
        {
            "lane_id": lane,
            "guide_status": guide_coverage[lane],
            "cells_by_bc": {bc: lane_cells[(lane, bc)] for bc, _ in pairs},
        }
        for lane in expected_lanes
    ]
    provenance = {
        "assay": "10x Flex v1",
        "counter_version": counter_version,
        "counter_version_label": counter_version_label,
        "umi_policy": umi_policy,
        "mapping_policy": mapping_policy,
        "filtering_policy": (
            "Bustools allowlist (default threshold), correct, sort, then count; "
            "applied per physical lane x canonical BC alias to the GEX BUS."
        ),
        "equivalence_class_metrics_scope": (
            "Read-weighted EC metrics are measured after component barcode correction "
            "and suffix replacement, before allowlist filtering."
        ),
        "sample_id": args.sample_id,
        "pool_id": args.pool_id,
        "bc_cr_pairs": pairs,
        "lane_count_products": validated_products,
        "nonempty_lane_count_products": len(inputs),
        "cell_count": cell_count,
        "guide_coverage": coverage_records,
        "lane_metrics": [
            {"lane_id": lane, **lane_metrics[lane]} for lane in expected_lanes
        ],
        "reference_hashes": json.loads(reference_hashes_json),
        "gex_manifest_sha256": file_sha256(args.gex_sources),
        "guide_manifest_sha256": file_sha256(args.guide_sources),
        "guide_coverage_sha256": file_sha256(args.guide_coverage),
        "samples_sha256": file_sha256(args.samples),
        "gene_probes_sha256": file_sha256(args.gene_probes),
        "probe_barcodes_sha256": file_sha256(args.probe_barcodes),
        "guide_features_sha256": file_sha256(args.guide_features),
        "gene_count_features_sha256": file_sha256(args.gene_count_features),
        "guide_count_features_sha256": file_sha256(args.guide_count_features),
        "guide_targets_sha256": file_sha256(args.guide_targets),
    }
    with h5py.File(first_product, "r") as source, h5py.File(args.output, "a") as target:
        if "uns" in target:
            del target["uns"]
        if "uns" in source:
            source.copy("uns", target)
        for key in (
            "flex_lane_id",
            "flex_bc_alias",
            "flex_cr_alias",
            "flex_guide_coverage_status",
            "flex_lane_metrics_json",
        ):
            if key in target["uns"]:
                del target["uns"][key]
        if "guide_names" not in target["uns"]:
            raise ValueError("Concatenated H5AD lost the guide name vector")
        provenance_dataset = target["uns"].create_dataset(
            "flex_provenance_json",
            data=json.dumps(provenance, sort_keys=True),
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
        provenance_dataset.attrs["encoding-type"] = "string"
        provenance_dataset.attrs["encoding-version"] = "0.2.0"
    Path(args.metrics).write_text(
        json.dumps(provenance, sort_keys=True, indent=2) + "\n"
    )
    result = ad.read_h5ad(args.output, backed="r")
    try:
        if result.n_obs != cell_count or result.n_vars != len(var_ids):
            raise ValueError("On-disk Flex aggregation dimensions do not match inputs")
        if result.obsm["guides"].shape != (cell_count, len(guide_names)):
            raise ValueError("On-disk guide count matrix dimensions do not match")
    finally:
        result.file.close()
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--guide-coverage", required=True)
    parser.add_argument("--gex-sources", required=True)
    parser.add_argument("--guide-sources", required=True)
    parser.add_argument("--gene-probes", required=True)
    parser.add_argument("--probe-barcodes", required=True)
    parser.add_argument("--gene-count-features", required=True)
    parser.add_argument("--guide-features", required=True)
    parser.add_argument("--guide-count-features", required=True)
    parser.add_argument("--guide-targets", required=True)
    parser.add_argument("--lane-counts", nargs="+", required=True)
    parser.add_argument("--output", default="flex_sample_uncompressed.h5ad")
    parser.add_argument("--metrics", default="flex_sample_metrics.json")
    parser.add_argument("--max-loaded-elems", type=int, default=100_000_000)
    args = parser.parse_args()
    if args.max_loaded_elems < 1:
        parser.error("--max-loaded-elems must be positive")
    aggregate(args)


if __name__ == "__main__":
    main()
