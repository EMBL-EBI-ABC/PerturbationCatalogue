#!/usr/bin/env python3
"""Small synthetic checks for native Flex counting and H5AD assembly."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


PIPELINE = Path(__file__).resolve().parents[1]
DATA = PIPELINE / "datasets" / "zhu_2025"
sys.path.insert(0, str(PIPELINE / "bin"))
import cyto_flex_lane  # noqa: E402
import flex_aggregate  # noqa: E402


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_rows(path):
    return [
        line.rstrip("\r\n").split("\t")
        for line in Path(path).read_text().splitlines()
        if line
    ]


def write_fastq_pair(directory, name, r1, r2):
    r1_path, r2_path = (
        directory / f"{name}_R1.fastq.gz",
        directory / f"{name}_R2.fastq.gz",
    )
    for path, sequence in ((r1_path, r1), (r2_path, r2)):
        with gzip.open(path, "wt") as handle:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    return r1_path, r2_path


def run(command):
    subprocess.run(
        [str(value) for value in command], check=True, stdout=subprocess.DEVNULL
    )


def mtx_sum_for_feature(directory, feature_id):
    features_path = next(
        path
        for path in (directory / "features.txt.gz", directory / "features.tsv.gz")
        if path.is_file()
    )
    matrix_path = directory / "matrix.mtx.gz"
    with gzip.open(features_path, "rt") as handle:
        features = [line.rstrip("\r\n").split("\t", 1)[0] for line in handle]
    matrix = mmread(matrix_path).tocsr()
    return int(matrix[features.index(feature_id), :].sum())


def native_count_two_chunks(
    cyto, resources, source_features, count_features, alias, reads, expected_feature
):
    with tempfile.TemporaryDirectory(prefix="cyto-flex-test-") as temp:
        root = Path(temp)
        whitelist = Path(resources) / "737K-fixed-rna-profiling.txt.gz"
        barcode = gzip.open(whitelist, "rt").readline().strip()
        r1 = barcode + "ACGTACGTACGT"
        mapped = []
        for index, feature_read in enumerate(reads, 1):
            r1_path, r2_path = write_fastq_pair(root, f"chunk{index}", r1, feature_read)
            out = root / f"map{index}"
            run(
                [
                    cyto,
                    "map",
                    "gex" if alias.startswith("BC") else "crispr",
                    "--preset",
                    "gex-v1" if alias.startswith("BC") else "crispr-v1",
                    "--whitelist",
                    whitelist,
                    "--gex" if alias.startswith("BC") else "--guides",
                    source_features,
                    "--probes",
                    Path(resources)
                    / (
                        "probe-barcodes-fixed-rna-profiling-rna.txt"
                        if alias.startswith("BC")
                        else "probe-barcodes-fixed-rna-profiling-crispr.txt"
                    ),
                    "--probe-regex",
                    f"^(?:{alias})$",
                    "--num-threads",
                    1,
                    "--min-ibu-records",
                    0,
                    "--outdir",
                    out,
                    r1_path,
                    r2_path,
                ]
            )
            ibu = out / "ibu" / f"{alias}.ibu"
            if not ibu.is_file():
                raise AssertionError(f"Cyto did not emit expected alias IBU: {alias}")
            mapped.append(ibu)
        combined = root / "combined.ibu"
        run([cyto, "ibu", "cat", "-t", 1, "-o", combined, *mapped])
        decoded = subprocess.run(
            [cyto, "ibu", "view", "-i", str(combined)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        ).stdout
        records = [
            line for line in decoded.splitlines() if line and not line.startswith("#")
        ]
        if len(records) != 2:
            raise AssertionError(
                f"Expected two mapped source records, found {len(records)}"
            )
        feature_indices = {line.split()[-1] for line in records}
        if (alias.startswith("BC") and len(feature_indices) != 2) or (
            alias.startswith("CR") and len(feature_indices) != 1
        ):
            raise AssertionError(
                "Synthetic reads did not map to the expected one/two feature identities"
            )
        sort = root / "sorted.ibu"
        corrected = root / "corrected.ibu"
        run([cyto, "ibu", "sort", "-i", combined, "-o", sort, "-T", 1])
        run([cyto, "ibu", "umi", "-i", sort, "-o", corrected, "-T", 1])
        count = root / "counts"
        run(
            [
                cyto,
                "ibu",
                "count",
                "-i",
                corrected,
                "--mtx",
                "-f",
                count_features,
                "-C",
                1,
                "-s",
                f"{alias}_lane",
                "-o",
                count,
                "-t",
                1,
            ]
        )
        return mtx_sum_for_feature(count, expected_feature)


def test_native_dedup(cyto, resources):
    gex_rows = [
        {"gene_id": fields[0], "probe_seq": fields[2]}
        for fields in read_rows(DATA / "gene_probes.tsv")
    ]
    by_gene = {}
    for row in gex_rows:
        by_gene.setdefault(row["gene_id"], []).append(row)
    gene_id, two_probes = next(
        (key, rows[:2]) for key, rows in by_gene.items() if len(rows) >= 2
    )
    bc_sequence = next(
        row["bc_sequence"]
        for row in read_tsv(DATA / "probe_barcodes.tsv")
        if row["bc_alias"] == "BC001"
    )
    gex_reads = [
        row["probe_seq"] + "ACGTACGTACGTACGTAC" + bc_sequence for row in two_probes
    ]
    gene_count = native_count_two_chunks(
        cyto,
        resources,
        DATA / "gene_probes.tsv",
        DATA / "gene_count_features.tsv",
        "BC001",
        gex_reads,
        gene_id,
    )
    if gene_count != 1:
        raise AssertionError(
            f"Cyto 0.4.5 first-group two-probe smoke result changed: {gene_count}"
        )

    guide_id, anchor, spacer = read_rows(DATA / "guide_features.tsv")[0]
    barcode = next(
        row
        for row in read_tsv(DATA / "probe_barcodes.tsv")
        if row["cr_alias"] == "CR001"
    )
    guide_read = barcode["cr_sequence"] + anchor + spacer
    guide_count = native_count_two_chunks(
        cyto,
        resources,
        DATA / "guide_features.tsv",
        DATA / "guide_count_features.tsv",
        "CR001",
        [guide_read, guide_read],
        guide_id,
    )
    if guide_count != 1:
        raise AssertionError(
            f"Repeated guide UMI across two chunks should count once, found {guide_count}"
        )
    print(
        "Native Cyto smoke passed: first-group same-UMI GEX probe fixture="
        f"{gene_count}; duplicate guide across chunks={guide_count}"
    )


def make_lane_h5ad(path, sample_id, lane, bc, cr, cell, count):
    barcode_rows = read_tsv(DATA / "probe_barcodes.tsv")
    bc_sequence = next(
        row["bc_sequence"] for row in barcode_rows if row["bc_alias"] == bc
    )
    obs_name = f"{cell}{bc_sequence}-1_{lane}"
    obs = pd.DataFrame(
        {
            "sample_id": [sample_id],
            "lane_id": [lane],
            "bc_alias": [bc],
            "cr_alias": [cr],
            "guide_coverage_status": ["available"],
        },
        index=pd.Index([obs_name]),
    )
    data = ad.AnnData(
        X=sparse.csr_matrix([[count]], dtype=np.int32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(["ENSG_TEST"], name="gene_id")),
    )
    data.obsm["guides"] = sparse.csr_matrix([[count]], dtype=np.int32)
    data.uns["guide_names"] = ["TARGET_guide__ENSG_TEST"]
    data.uns["flex_sample_id"] = sample_id
    data.uns["flex_lane_id"] = lane
    data.uns["flex_bc_alias"] = bc
    data.uns["flex_cr_alias"] = cr
    data.uns["flex_guide_coverage_status"] = "available"
    data.uns["flex_cyto_version"] = "cyto 0.4.5"
    data.write_h5ad(path, compression="gzip")


def test_aggregation_and_cell_keys():
    with tempfile.TemporaryDirectory(prefix="flex-aggregate-test-") as temp:
        root = Path(temp)
        sample_id, pool_id = "sample_test", "pool_test"
        lanes = ["laneA", "laneB"]
        pairs = [("BC001", "CR001"), ("BC002", "CR002")]
        sample_path = root / "samples.tsv"
        sample_path.write_text(
            "sample_id\tpool_id\tbc_cr_pairs\tpool_lanes\n"
            f"{sample_id}\t{pool_id}\tBC001:CR001;BC002:CR002\t{';'.join(lanes)}\n"
        )
        coverage = root / "coverage.tsv"
        coverage.write_text(
            "pool_id\tlane_id\tguide_status\n"
            + "".join(f"{pool_id}\t{lane}\tavailable\n" for lane in lanes)
        )
        barcode_path = root / "barcodes.tsv"
        aliases = {
            row["bc_alias"]: row for row in read_tsv(DATA / "probe_barcodes.tsv")
        }
        barcode_path.write_text(
            "bc_alias\tbc_sequence\tcr_alias\tcr_sequence\n"
            f"BC001\t{aliases['BC001']['bc_sequence']}\tCR001\t{aliases['BC001']['cr_sequence']}\n"
            f"BC002\t{aliases['BC002']['bc_sequence']}\tCR002\t{aliases['BC002']['cr_sequence']}\n"
        )
        (
            gex_sources,
            guide_sources,
            probes,
            gene_counts,
            guides,
            guide_counts,
            targets,
        ) = ([], [], [], [], [], [], [])
        for name, content in (
            ("gex.tsv", "gex\n"),
            ("guide.tsv", "guide\n"),
            ("probes.tsv", "probe\n"),
            ("gene_counts.tsv", "gene_count\n"),
            ("guides.tsv", "guide\n"),
            ("guide_counts.tsv", "guide_count\n"),
            ("targets.tsv", "target\n"),
        ):
            path = root / name
            path.write_text(content)
            if name == "gex.tsv":
                gex_sources = path
            elif name == "guide.tsv":
                guide_sources = path
            elif name == "probes.tsv":
                probes = path
            elif name == "gene_counts.tsv":
                gene_counts = path
            elif name == "guides.tsv":
                guides = path
            elif name == "guide_counts.tsv":
                guide_counts = path
            else:
                targets = path
        products = []
        cell = "AAACAAGCAAACAAGA"  # same 16-base CBC for each barcode/lane combination
        cr_by_bc = dict(pairs)
        for lane in lanes:
            for bc, _ in pairs:
                path = root / f"{lane}__{bc}.h5ad"
                make_lane_h5ad(path, sample_id, lane, bc, cr_by_bc[bc], cell, 1)
                products.append(path)
        output, metrics = root / "sample.h5ad", root / "metrics.json"
        flex_aggregate.aggregate(
            argparse.Namespace(
                samples=sample_path,
                sample_id=sample_id,
                pool_id=pool_id,
                probe_barcodes=barcode_path,
                guide_coverage=coverage,
                lane_counts=products,
                gex_sources=gex_sources,
                guide_sources=guide_sources,
                gene_probes=probes,
                gene_count_features=gene_counts,
                guide_features=guides,
                guide_count_features=guide_counts,
                guide_targets=targets,
                output=output,
                metrics=metrics,
                max_loaded_elems=10_000,
            )
        )
        result = ad.read_h5ad(output)
        try:
            if result.n_obs != 4 or not result.obs_names.is_unique:
                raise AssertionError(
                    "Lane and BC8 namespaces did not preserve all four cells"
                )
            if set(result.obs["lane_id"]) != set(lanes):
                raise AssertionError("Lane IDs were lost during on-disk aggregation")
            if any(
                key in result.uns
                for key in (
                    "flex_lane_id",
                    "flex_bc_alias",
                    "flex_cr_alias",
                    "flex_guide_coverage_status",
                )
            ):
                raise AssertionError(
                    "Sample-level H5AD retained misleading single-lane metadata"
                )
            if result.uns["flex_cyto_version"] != "cyto 0.4.5":
                raise AssertionError(
                    "Sample-level H5AD lost the native counter version"
                )
            if int(result.obsm["guides"].sum()) != 4:
                raise AssertionError(
                    "Guide modality did not survive on-disk aggregation"
                )
            expected_barcodes = {
                f"{cell}{aliases['BC001']['bc_sequence']}-1_laneA",
                f"{cell}{aliases['BC002']['bc_sequence']}-1_laneA",
                f"{cell}{aliases['BC001']['bc_sequence']}-1_laneB",
                f"{cell}{aliases['BC002']['bc_sequence']}-1_laneB",
            }
            if set(result.obs_names) != expected_barcodes:
                raise AssertionError(
                    f"Full CBC16+BC8+lane cell identity was not preserved: "
                    f"expected={expected_barcodes}, found={set(result.obs_names)}"
                )
        finally:
            if result.isbacked:
                result.file.close()
        provenance = json.loads(metrics.read_text())
        if provenance["cell_count"] != 4:
            raise AssertionError("Aggregation provenance cell count differs")
        if (
            provenance["cyto_version"] != "cyto 0.4.5"
            or not provenance["native_umi_semantics"]
        ):
            raise AssertionError("Aggregation provenance lacks counter semantics")


def test_lane_h5ad_guide_reindex():
    with tempfile.TemporaryDirectory(prefix="flex-lane-h5ad-test-") as temp:
        root = Path(temp)
        gex_cells = ["A" * 16, "C" * 16]
        guide_cells = ["C" * 16]
        guide_matrix = sparse.csr_matrix([[0, 1]], dtype=np.int32)
        sample = {"sample_id": "sample_test"}
        output = root / "lane.h5ad"
        cyto_flex_lane.write_alias_h5ad(
            output,
            "BC001",
            "CR001",
            "laneA",
            sample,
            {"BC001": "AACCGGTT"},
            gex_cells,
            sparse.csr_matrix([[1], [2]], dtype=np.int32),
            (guide_cells, guide_matrix),
            ["ENSG_TEST"],
            ["TEST"],
            ["guideA", "guideB"],
            {"guideA": "A_guideA", "guideB": "B_guideB"},
            "available",
            "cyto 0.4.5",
        )
        data = ad.read_h5ad(output)
        try:
            if data.obsm["guides"].shape != (2, 2):
                raise AssertionError(
                    "Guide rows were not expanded to the full expression-cell shape"
                )
            if (
                int(data.obsm["guides"].sum()) != 1
                or int(data.obsm["guides"][1, 1]) != 1
            ):
                raise AssertionError(
                    "Guide counts were not assigned to the matching expression cell"
                )
        finally:
            if data.isbacked:
                data.file.close()

        empty = root / "empty.h5ad"
        cyto_flex_lane.create_empty_h5ad(
            empty,
            sample,
            "laneA",
            "BC001",
            "CR001",
            "no_archived_guide_sra",
            ["ENSG_TEST"],
            ["TEST"],
            ["A_guideA"],
            "cyto 0.4.5",
        )
        empty_data = ad.read_h5ad(empty)
        try:
            if empty_data.n_obs != 0 or empty_data.obsm["guides"].shape != (0, 1):
                raise AssertionError("Empty alias product has an invalid H5AD schema")
        finally:
            if empty_data.isbacked:
                empty_data.file.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cyto", default="/usr/local/bin/cyto")
    parser.add_argument("--resources", default="/opt/cyto/resources")
    args = parser.parse_args()
    test_native_dedup(args.cyto, args.resources)
    test_lane_h5ad_guide_reindex()
    test_aggregation_and_cell_keys()
    print(
        "Flex H5AD shape, BC8 identity, lane namespace and on-disk aggregation passed"
    )


if __name__ == "__main__":
    main()
