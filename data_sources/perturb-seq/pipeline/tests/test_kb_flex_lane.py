#!/usr/bin/env python3
"""Exercise lane-wide EC remapping, alias partitioning, and guide barcode joins."""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE / "bin"))
import kb_flex_lane as lane


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


if __name__ == "__main__":
    test_lane_native()
