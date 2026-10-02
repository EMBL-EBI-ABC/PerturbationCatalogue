#!/usr/bin/env python3
"""Exercise lane-wide EC remapping, alias partitioning, and guide barcode joins."""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "bin"))
import kb_flex_lane as lane
import flex_aggregate


def run(args):
    result = subprocess.run(
        [str(value) for value in args], capture_output=True, text=True, check=False
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {args}\n{result.stderr}"
        )
    return result


def write_bus(path, rows):
    records = np.zeros(len(rows), dtype=lane.BUS_RECORD)
    for index, (barcode, umi, ec_id) in enumerate(rows):
        records[index] = (
            lane.encode_dna(barcode),
            lane.encode_dna(umi),
            ec_id,
            1,
            0,
            0,
        )
    path = Path(path)
    path.write_bytes(lane.BUS_HEADER.pack(b"BUS\0", 1, 24, 12, 0) + records.tobytes())
    return path


def native_count(bustools, bus_path, ec_path, tx_path, t2g_path, targets, root, name):
    sorted_bus = root / f"{name}.sorted.bus"
    prefix = root / name
    run([bustools, "sort", "-t", 1, "-o", sorted_bus, bus_path])
    run(
        [
            bustools,
            "count",
            "-o",
            prefix,
            "-g",
            t2g_path,
            "-e",
            ec_path,
            "-t",
            tx_path,
            "--genecounts",
            "--umi-gene",
            sorted_bus,
        ]
    )
    return lane.read_count_matrix(prefix, targets)


def test_lane_native():
    from kb_python.utils import get_bustools_binary_path

    bustools = str(get_bustools_binary_path())
    if run([bustools, "version"]).stdout.strip() != "bustools, version 0.45.1":
        raise AssertionError("Run this regression inside the pinned Flex SIF")

    cbc = "ACGTACGTACGTACGT"
    bc1, bc2 = "AACCGGTT", "CCGGAATT"
    cr1 = "TTCCAAGG"
    barcode1, barcode2 = cbc + bc1, cbc + bc2
    alias_sequences = {"BC001": bc1, "BC002": bc2}
    # Deliberately nonalphabetical: native feature order must follow the
    # reference/T2G files, not a lexical sort.
    targets = ["probe_Z", "probe_A"]

    with tempfile.TemporaryDirectory(prefix="kb-flex-lane-native-") as temp:
        root = Path(temp)
        gex_writer = lane.AliasBusWriter(
            root / "gex_aliases", sorted(alias_sequences), alias_sequences
        )
        gex_counter = lane.CanonicalEC(
            targets,
            "gex",
            root / "gex_reference",
            gex_writer,
            {"probe_Z": "ENSG1", "probe_A": "ENSG1"},
        )
        tx = root / "transcripts.txt"
        tx.write_text("probe_Z\nprobe_A\n")

        first_ec = root / "source1.ec"
        first_ec.write_text("0\t0\n1\t1\n")
        first_bus = write_bus(
            root / "source1.bus",
            [
                (barcode1, "AAAAAAAAAAAA", 0),
                (barcode1, "CCCCCCCCCCCC", 1),
            ],
        )
        gex_counter.add_source(first_bus, first_ec, tx)

        # Source two uses local EC 0 for probe_B and EC 1 for probe_A.
        # Its probe_A/UMI record duplicates source one's molecule.
        second_ec = root / "source2.ec"
        second_ec.write_text("0\t1\n1\t0\n")
        second_bus = write_bus(
            root / "source2.bus",
            [
                (barcode1, "AAAAAAAAAAAA", 1),
                (barcode2, "GGGGGGGGGGGG", 0),
            ],
        )
        gex_counter.add_source(second_bus, second_ec, tx)
        gex_ec, gex_tx = gex_counter.write_reference()

        # Equivalent target sets must converge to lane-wide IDs despite source
        # order; different BC suffixes from the same CBC remain separate cells.
        if gex_ec.read_text() != "0\t0\n1\t1\n":
            raise AssertionError(f"Unexpected canonical EC table: {gex_ec.read_text()}")
        gex_t2g = root / "gex.t2g"
        gex_t2g.write_text("probe_Z\tprobe_Z\nprobe_A\tprobe_A\n")
        bc001_barcodes, bc001_matrix = native_count(
            bustools,
            gex_writer.paths["BC001"],
            gex_ec,
            gex_tx,
            gex_t2g,
            targets,
            root,
            "gex_BC001",
        )
        bc002_barcodes, bc002_matrix = native_count(
            bustools,
            gex_writer.paths["BC002"],
            gex_ec,
            gex_tx,
            gex_t2g,
            targets,
            root,
            "gex_BC002",
        )
        if bc001_barcodes != [barcode1] or bc002_barcodes != [barcode2]:
            raise AssertionError(
                "Alias partition did not preserve full 24-base barcodes"
            )
        if bc001_matrix.toarray().tolist() != [[1, 1]]:
            raise AssertionError(
                "Source-local EC reordering or cross-source UMI dedup failed"
            )
        if bc002_matrix.toarray().tolist() != [[0, 1]]:
            raise AssertionError("The second BC alias was mixed with the first alias")

        invalid_ec = root / "negative_target.ec"
        invalid_ec.write_text("0\t-1\n")
        try:
            lane.parse_ec_map(invalid_ec, len(targets))
        except ValueError as error:
            if "outside the target list" not in str(error):
                raise
        else:
            raise AssertionError("A negative target index was accepted")

        # A corrected BUS barcode outside the pinned aliases must fail closed.
        invalid_writer = lane.AliasBusWriter(
            root / "invalid_alias", sorted(alias_sequences), alias_sequences
        )
        invalid_record = np.zeros(1, dtype=lane.BUS_RECORD)
        invalid_record[0]["barcode"] = lane.encode_dna(cbc + "AAAAAAAA")
        try:
            invalid_writer.append(invalid_record, lane.bus_header(first_bus))
        except ValueError as error:
            if "unrecognized alias suffix" not in str(error):
                raise
        else:
            raise AssertionError("An unknown alias suffix was silently accepted")
        finally:
            invalid_writer.close()

        # Guide CR is replaced with the paired canonical GEX BC suffix; the
        # unchanged CBC plus translated suffix must join only the matching cell.
        guide_target = "guide_1"
        guide_writer = lane.AliasBusWriter(
            root / "guide_aliases", sorted(alias_sequences), alias_sequences
        )
        guide_counter = lane.CanonicalEC(
            [guide_target], "guide", root / "guide_reference", guide_writer
        )
        guide_tx = root / "guide_transcripts.txt"
        guide_tx.write_text(guide_target + "\n")
        guide_ec_local = root / "guide_local.ec"
        guide_ec_local.write_text("0\t0\n")
        guide_raw = write_bus(root / "guide_raw.bus", [(cbc + cr1, "TTTTTTTTTTTT", 0)])
        replacement = root / "cr_to_bc.tsv"
        replacement.write_text(f"{cr1}\t*{bc1}\n")
        guide_replaced = root / "guide_replaced.bus"
        run(
            [
                bustools,
                "correct",
                "-r",
                "-w",
                replacement,
                "-o",
                guide_replaced,
                guide_raw,
            ]
        )
        guide_counter.add_source(guide_replaced, guide_ec_local, guide_tx)
        guide_ec, guide_tx = guide_counter.write_reference()
        guide_t2g = root / "guide.t2g"
        guide_t2g.write_text(f"{guide_target}\t{guide_target}\n")
        guide_barcodes, guide_matrix = native_count(
            bustools,
            guide_writer.paths["BC001"],
            guide_ec,
            guide_tx,
            guide_t2g,
            [guide_target],
            root,
            "guide_BC001",
        )
        if guide_barcodes != [barcode1] or guide_matrix.toarray().tolist() != [[1]]:
            raise AssertionError(
                "Guide CR replacement did not yield the paired GEX barcode"
            )

        joined, matched = lane.join_guide_counts(
            [barcode1, barcode2], (guide_barcodes, guide_matrix), 1
        )
        if matched != 1 or joined.toarray().tolist() != [[1], [0]]:
            raise AssertionError(
                "Guide counts did not join by the full CBC+GEX BC identity"
            )

    print(
        "PASS: source EC reorder, cross-source UMI dedup, BC alias partition, "
        "CR-to-paired-BC replacement, and 24-base guide join"
    )


def test_aggregate_missing_guides():
    with tempfile.TemporaryDirectory(prefix="kb-flex-aggregate-") as temp:
        root = Path(temp)
        samples = root / "samples.tsv"
        samples.write_text(
            "sample_id\tpool_id\tbc_cr_pairs\tpool_lanes\n"
            "S1\tP1\tBC001:CR001;BC002:CR002\tlaneA;laneB\n"
        )
        coverage = root / "coverage.tsv"
        coverage.write_text(
            "pool_id\tlane_id\tguide_status\n"
            "P1\tlaneA\tavailable\n"
            "P1\tlaneB\tno_archived_guide_sra\n"
        )
        barcode_file = root / "probe_barcodes.tsv"
        barcode_file.write_text(
            "bc_alias\tbc_sequence\tcr_alias\tcr_sequence\n"
            "BC001\tAACCGGTT\tCR001\tTTCCAAGG\n"
            "BC002\tCCGGAATT\tCR002\tGGTTAACC\n"
            "BC003\tACGTACGT\tCR003\tTGCATGCA\n"
        )
        metadata_files = {}
        for name in (
            "gex_sources",
            "guide_sources",
            "gene_probes",
            "gene_count_features",
            "guide_features",
            "guide_count_features",
            "guide_targets",
        ):
            path = root / f"{name}.tsv"
            path.write_text(name + "\n")
            metadata_files[name] = path

        counter_version = {
            "kb_python": "0.30.2",
            "kallisto": "kallisto, version 0.52.0",
            "bustools": "bustools, version 0.45.1",
        }
        reference_hashes = json.dumps(
            {"gex_index": "a" * 64, "guide_index": "b" * 64}, sort_keys=True
        )
        sequences = {
            "BC001": "AACCGGTT",
            "BC002": "CCGGAATT",
            "BC003": "ACGTACGT",
        }
        cr_aliases = {"BC001": "CR001", "BC002": "CR002", "BC003": "CR003"}
        cbc = "ACGTACGTACGTACGT"
        products = []
        for lane_id, guide_status in (
            ("laneA", "available"),
            ("laneB", "no_archived_guide_sra"),
        ):
            for alias, suffix in sequences.items():
                # Same CBC across both BC8 aliases and both lanes. BC8 and lane
                # identity must keep all nonempty cells distinct.
                empty_alias = lane_id == "laneB" and alias == "BC002"
                cell_ids = [] if empty_alias else [f"{cbc}{suffix}-1_{lane_id}"]
                obs = pd.DataFrame(
                    {
                        column: pd.Series(values, index=cell_ids, dtype=object)
                        for column, values in {
                            "sample_id": ["S1"] * len(cell_ids),
                            "lane_id": [lane_id] * len(cell_ids),
                            "bc_alias": [alias] * len(cell_ids),
                            "cr_alias": [cr_aliases[alias]] * len(cell_ids),
                            "guide_coverage_status": [guide_status] * len(cell_ids),
                        }.items()
                    },
                    index=pd.Index(cell_ids, dtype=object),
                )
                data = ad.AnnData(
                    X=sparse.csr_matrix(np.ones((len(cell_ids), 2), dtype=np.int32)),
                    obs=obs,
                    var=pd.DataFrame(
                        {"gene_name": ["G1", "G2"]},
                        index=pd.Index(["ENSG1", "ENSG2"], name="gene_id"),
                    ),
                )
                guide_matrix = sparse.lil_matrix((len(cell_ids), 1), dtype=np.int32)
                if cell_ids and lane_id == "laneA":
                    guide_matrix[0, 0] = 1
                data.obsm["guides"] = guide_matrix.tocsr()
                data.uns.update(
                    {
                        "guide_names": ["NTC_guide1"],
                        "flex_counter_version": counter_version,
                        "flex_counter_version_label": (
                            "kb-python 0.30.2; kallisto, version 0.52.0; "
                            "bustools, version 0.45.1"
                        ),
                        "flex_umi_policy": "probe IDs counted before stable ENSG sum",
                        "flex_mapping_policy": "ambiguous ECs are not allocated",
                        "flex_reference_hashes_json": reference_hashes,
                        "flex_sample_id": "S1",
                        "flex_lane_id": lane_id,
                        "flex_bc_alias": alias,
                        "flex_cr_alias": cr_aliases[alias],
                        "flex_guide_coverage_status": guide_status,
                        "flex_lane_metrics_json": json.dumps(
                            {
                                "pool_id": "P1",
                                "lane_id": lane_id,
                                "bc_alias": alias,
                                "cr_alias": cr_aliases[alias],
                                "guide_coverage_status": guide_status,
                                "barcode_metrics": {
                                    "cells": len(cell_ids),
                                    "observed_composite_barcodes": len(cell_ids),
                                    "retained_cells": len(cell_ids),
                                    "expression_umis": len(cell_ids),
                                    "genes_detected": 2 if cell_ids else 0,
                                    "guide_positive_cells": int(guide_matrix.sum()),
                                },
                                "gex_mapping": {"bus_records": 10},
                                "guide_mapping": {"bus_records": 2},
                                "equivalence_class_metrics_scope": (
                                    "post-correction, before allowlist"
                                ),
                            },
                            sort_keys=True,
                        ),
                    }
                )
                product = root / f"{lane_id}__{alias}.h5ad"
                data.write_h5ad(product, compression="gzip")
                products.append(product)

        output, metrics = root / "sample.h5ad", root / "metrics.json"
        flex_aggregate.aggregate(
            argparse.Namespace(
                samples=samples,
                sample_id="S1",
                pool_id="P1",
                probe_barcodes=barcode_file,
                guide_coverage=coverage,
                lane_counts=products,
                gex_sources=metadata_files["gex_sources"],
                guide_sources=metadata_files["guide_sources"],
                gene_probes=metadata_files["gene_probes"],
                gene_count_features=metadata_files["gene_count_features"],
                guide_features=metadata_files["guide_features"],
                guide_count_features=metadata_files["guide_count_features"],
                guide_targets=metadata_files["guide_targets"],
                output=output,
                metrics=metrics,
                max_loaded_elems=1000,
            )
        )
        result = ad.read_h5ad(output)
        try:
            expected_ids = {
                f"{cbc}{sequences[alias]}-1_{lane_id}"
                for lane_id in ("laneA", "laneB")
                for alias in ("BC001", "BC002")
                if not (lane_id == "laneB" and alias == "BC002")
            }
            if set(result.obs_names.astype(str)) != expected_ids:
                raise AssertionError(
                    "CBC16+BC8+lane cell identities were not preserved"
                )
            if result.n_obs != 3 or not result.obs_names.is_unique:
                raise AssertionError(
                    "Duplicate or missing sample cells after aggregation"
                )
            if result.obsm["guides"].shape != (3, 1):
                raise AssertionError("The sparse guide matrix has the wrong dimensions")
            lane_b = result.obs["lane_id"].astype(str) == "laneB"
            lane_a = result.obs["lane_id"].astype(str) == "laneA"
            if (
                not result.obs.loc[lane_b, "guide_coverage_status"]
                .astype(str)
                .eq("no_archived_guide_sra")
                .all()
            ):
                raise AssertionError("Missing-guide archive status was lost")
            if result.obsm["guides"][lane_b.to_numpy()].nnz != 0:
                raise AssertionError("Missing-guide lane has nonzero guide counts")
            if result.obsm["guides"][lane_a.to_numpy()].toarray().tolist() != [
                [1],
                [1],
            ]:
                raise AssertionError(
                    "Available guide counts did not align to GEX cells"
                )
            if result.var["gene_name"].tolist() != ["G1", "G2"]:
                raise AssertionError("Pinned gene annotations were lost in disk concat")
            if any(
                key in result.uns
                for key in (
                    "flex_lane_id",
                    "flex_bc_alias",
                    "flex_cr_alias",
                    "flex_guide_coverage_status",
                    "flex_lane_metrics_json",
                )
            ):
                raise AssertionError("Sample H5AD retained lane-specific provenance")
            provenance = json.loads(result.uns["flex_provenance_json"])
            if len(provenance["lane_metrics"]) != 2:
                raise AssertionError("Lane metrics were not retained in provenance")
            bc2_lane_b = provenance["lane_metrics"][1]["barcode_metrics"]["BC002"]
            if bc2_lane_b["retained_cells"] != 0:
                raise AssertionError("Empty lane/alias products were not preserved")
            if provenance["counter_version"]["bustools"] != "bustools, version 0.45.1":
                raise AssertionError("Native counter version was not retained")
        finally:
            if result.isbacked:
                result.file.close()
    print(
        "PASS: sparse on-disk aggregation, exact composite cell IDs, empty alias, "
        "missing-guide coverage, feature annotations, and native provenance"
    )


if __name__ == "__main__":
    test_lane_native()
    test_aggregate_missing_guides()
