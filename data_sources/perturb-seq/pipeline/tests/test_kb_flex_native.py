#!/usr/bin/env python3
"""Exercise Flex read geometry, whitelist correction, and probe-level UMI counting."""

import csv
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


PIPELINE = Path(__file__).resolve().parents[1]
DATA = PIPELINE / "datasets" / "zhu_2025"
TECHNOLOGY = "0,0,16,1,68,76:0,16,28:1,0,50"


def run(args, *, stdout=subprocess.PIPE):
    result = subprocess.run(
        [str(value) for value in args],
        check=False,
        stdout=stdout,
        stderr=subprocess.PIPE,
        text=stdout is not subprocess.DEVNULL,
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {args}\n{result.stderr}"
        )
    return result


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta(path):
    name = None
    sequence = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name, sequence = line[1:], []
            elif line:
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence)


def hamming_one_neighbors(sequence):
    for index, base in enumerate(sequence):
        for replacement in "ACGT":
            if replacement != base:
                yield sequence[:index] + replacement + sequence[index + 1 :]


def unique_error(sequence, whitelist):
    for observed in hamming_one_neighbors(sequence):
        candidates = [
            candidate
            for candidate in hamming_one_neighbors(observed)
            if candidate in whitelist
        ]
        if observed not in whitelist and candidates == [sequence]:
            return observed
    raise AssertionError(f"No uniquely correctable single-base error for {sequence}")


def _has_unique_error(sequence, whitelist):
    try:
        unique_error(sequence, whitelist)
    except AssertionError:
        return False
    return True


def write_fastq_pair(root, name, records):
    r1_path, r2_path = root / f"{name}_R1.fastq", root / f"{name}_R2.fastq"
    with r1_path.open("w") as r1, r2_path.open("w") as r2:
        for index, (cbc, umi, probe, raw_bc) in enumerate(records):
            read_name = f"{name}_{index}"
            seq1 = cbc + umi
            seq2 = probe + "ACGTACGTACGTACGTAC" + raw_bc
            r1.write(f"@{read_name}\n{seq1}\n+\n{'I' * len(seq1)}\n")
            r2.write(f"@{read_name}\n{seq2}\n+\n{'I' * len(seq2)}\n")
    return r1_path, r2_path


def bus_text(bustools, bus_path, output_path):
    run([bustools, "text", "-o", output_path, bus_path], stdout=subprocess.DEVNULL)
    return [line.split("\t") for line in Path(output_path).read_text().splitlines()]


def matrix_counts(directory):
    genes = (directory.parent / f"{directory.name}.genes.txt").read_text().splitlines()
    barcodes = (
        (directory.parent / f"{directory.name}.barcodes.txt").read_text().splitlines()
    )
    with (directory.parent / f"{directory.name}.mtx").open() as handle:
        dimensions = next(line for line in handle if not line.startswith("%"))
        n_barcodes, n_genes, _ = map(int, dimensions.split())
        if (n_barcodes, n_genes) != (len(barcodes), len(genes)):
            raise AssertionError(
                "Bustools MatrixMarket dimensions do not match its row/column names"
            )
        counts = {}
        for line in handle:
            row, column, count = map(int, line.split())
            key = (genes[column - 1], barcodes[row - 1])
            counts[key] = counts.get(key, 0) + count
    return counts


def make_onlist(path, cbcs, raw_bcs):
    with Path(path).open("w") as handle:
        handle.write(f"{cbcs[0]}\t{raw_bcs[0]}\n")
        for cbc in cbcs[1:]:
            handle.write(f"{cbc}\t-\n")
        for raw_bc in raw_bcs[1:]:
            handle.write(f"-\t{raw_bc}\n")


def test_native():
    import kb_python
    from kb_python.utils import get_bustools_binary_path, get_kallisto_binary_path

    kb = shutil.which("kb") or "/usr/local/bin/kb"
    kallisto = get_kallisto_binary_path()
    bustools = get_bustools_binary_path()
    if not (kb and kallisto and bustools):
        raise SystemExit(
            "Run inside the pinned Flex SIF with kb, kallisto, bustools on PATH"
        )

    versions = (
        kb_python.__version__,
        run([kallisto, "version"]).stdout.strip(),
        run([bustools, "version"]).stdout.strip(),
    )
    if versions != (
        "0.30.2",
        "kallisto, version 0.52.0",
        "bustools, version 0.45.1",
    ):
        raise AssertionError(f"Unexpected native tool versions: {versions}")

    targets = dict(fasta(DATA / "gex_probe_targets.fa"))
    probe_to_gene = {
        row[0]: row[1]
        for row in (
            line.rstrip("\r\n").split("\t")
            for line in (DATA / "gex_probe_to_gene.tsv").open()
        )
    }
    probes_by_gene = {}
    for probe_id, gene_id in probe_to_gene.items():
        probes_by_gene.setdefault(gene_id, []).append(probe_id)
    pair = None
    for gene_id, probe_ids in probes_by_gene.items():
        for first, second in itertools.combinations(probe_ids, 2):
            kmers_a = {targets[first][i : i + 31] for i in range(20)}
            kmers_b = {targets[second][i : i + 31] for i in range(20)}
            if kmers_a.isdisjoint(kmers_b):
                pair = gene_id, first, second
                break
        if pair:
            break
    if pair is None:
        raise AssertionError(
            "Could not find two same-gene probes with distinct 31-mers"
        )
    gene_id, probe_a_id, probe_b_id = pair
    probe_a, probe_b = targets[probe_a_id], targets[probe_b_id]
    if len(probe_a) != 50 or len(probe_b) != 50:
        raise AssertionError(
            "Flex GEX targets must be the observed 50-base probe sequences"
        )

    cbc_path = DATA / "cbc_whitelist.txt"
    cbcs = [line.strip() for line in cbc_path.open() if line.strip()]
    raw_to_canonical = {
        row["raw_sequence"]: row["canonical_sequence"]
        for row in read_tsv(DATA / "bc_barcode_variants.tsv")
    }
    raw_to_alias = {
        row["raw_sequence"]: row["alias"]
        for row in read_tsv(DATA / "bc_barcode_variants.tsv")
    }
    bc001 = next(raw for raw, alias in raw_to_alias.items() if alias == "BC001")
    bc002 = next(raw for raw, alias in raw_to_alias.items() if alias == "BC002")
    if raw_to_canonical[bc001] == raw_to_canonical[bc002]:
        raise AssertionError("BC001 and BC002 must have different canonical suffixes")
    cbc_set = set(cbcs)
    cbc = next(value for value in cbcs if _has_unique_error(value, cbc_set))
    cbc_error = unique_error(cbc, cbc_set)
    raw_bc_error = unique_error(bc001, set(raw_to_canonical))
    raw_bc_error_target = next(
        candidate
        for candidate in hamming_one_neighbors(raw_bc_error)
        if candidate == bc001
    )
    if raw_bc_error_target != bc001:
        raise AssertionError("Raw BC single-base error does not resolve to BC001")

    with tempfile.TemporaryDirectory(prefix="kb-flex-native-") as temp:
        root = Path(temp)
        reference = root / "probes.fa"
        reference.write_text(f">{probe_a_id}\n{probe_a}\n>{probe_b_id}\n{probe_b}\n")
        index = root / "probes.idx"
        run(
            [kb, "ref", "--workflow", "custom", "-i", index, "-k", 31, reference],
            stdout=subprocess.DEVNULL,
        )
        self_t2g = root / "self.t2g"
        self_t2g.write_text(f"{probe_a_id}\t{probe_a_id}\n{probe_b_id}\t{probe_b_id}\n")

        onlist = root / "cbc_bc.onlist"
        make_onlist(onlist, cbcs, list(raw_to_canonical))
        replacement = root / "bc_replace.tsv"
        replacement.write_text(
            "".join(
                f"{raw}\t*{canonical}\n" for raw, canonical in raw_to_canonical.items()
            )
        )

        # The same physical CBC is deliberately assigned both BC aliases.
        # The repeated probe-A/UMI record crosses the independent BUS chunks.
        chunk_records = (
            [
                (cbc_error, "AAAAAAAAAAAA", probe_a, raw_bc_error),
                (cbc_error, "GGGGGGGGGGGG", probe_a, bc001),
                (cbc, "AAAAAAAAAAAA", probe_b, bc001),
            ],
            [
                (cbc, "AAAAAAAAAAAA", probe_a, bc001),
                (cbc, "CCCCCCCCCCCC", probe_a, raw_bc_error),
                (cbc, "AAAAAAAAAAAA", probe_a, bc002),
            ],
        )
        buses = []
        ecmap = txnames = None
        for index_number, records in enumerate(chunk_records, 1):
            r1, r2 = write_fastq_pair(root, f"chunk{index_number}", records)
            outdir = root / f"map{index_number}"
            run(
                [
                    kallisto,
                    "bus",
                    "-i",
                    index,
                    "-o",
                    outdir,
                    "-x",
                    TECHNOLOGY,
                    "-t",
                    1,
                    r1,
                    r2,
                ],
                stdout=subprocess.DEVNULL,
            )
            bus = outdir / "output.bus"
            buses.append(bus)
            if ecmap is None:
                ecmap, txnames = outdir / "matrix.ec", outdir / "transcripts.txt"

        corrected = root / "corrected.bus"
        run(
            [bustools, "correct", "-w", onlist, "-o", corrected, *buses],
            stdout=subprocess.DEVNULL,
        )
        corrected_rows = bus_text(bustools, corrected, root / "corrected.tsv")
        if len(corrected_rows) != sum(map(len, chunk_records)):
            raise AssertionError("CBC/BC component correction lost a valid record")
        if not any(row[0] == cbc + bc001 for row in corrected_rows):
            raise AssertionError("CBC and raw BC single-base errors were not corrected")

        replaced = root / "replaced.bus"
        run(
            [bustools, "correct", "-r", "-w", replacement, "-o", replaced, corrected],
            stdout=subprocess.DEVNULL,
        )
        replaced_rows = bus_text(bustools, replaced, root / "replaced.tsv")
        expected_barcodes = {
            cbc + raw_to_canonical[bc001],
            cbc + raw_to_canonical[bc002],
        }
        observed_barcodes = {row[0] for row in replaced_rows}
        if observed_barcodes != expected_barcodes:
            raise AssertionError(
                f"Suffix replacement did not preserve CBC and distinct BC aliases: {observed_barcodes}"
            )

        sorted_bus = root / "sorted.bus"
        run(
            [bustools, "sort", "-t", 1, "-o", sorted_bus, replaced],
            stdout=subprocess.DEVNULL,
        )
        counts = root / "counts"
        run(
            [
                bustools,
                "count",
                "-o",
                counts,
                "-g",
                self_t2g,
                "-e",
                ecmap,
                "-t",
                txnames,
                "--umi-gene",
                "--genecounts",
                sorted_bus,
            ],
            stdout=subprocess.DEVNULL,
        )
        matrix = matrix_counts(counts)
        bc001_cell = cbc + raw_to_canonical[bc001]
        bc002_cell = cbc + raw_to_canonical[bc002]
        observed_counts = {
            (probe_a_id, bc001_cell): matrix.get((probe_a_id, bc001_cell), 0),
            (probe_b_id, bc001_cell): matrix.get((probe_b_id, bc001_cell), 0),
            (probe_a_id, bc002_cell): matrix.get((probe_a_id, bc002_cell), 0),
        }
        expected_counts = {
            (probe_a_id, bc001_cell): 3,
            (probe_b_id, bc001_cell): 1,
            (probe_a_id, bc002_cell): 1,
        }
        if observed_counts != expected_counts:
            raise AssertionError(
                f"Probe-local UMI counts differ, expected {expected_counts}, got {observed_counts}"
            )
        gene_sum = sum(
            count
            for (probe_id, barcode), count in matrix.items()
            if barcode == bc001_cell and probe_to_gene[probe_id] == gene_id
        )
        if gene_sum != 4:
            raise AssertionError(
                f"Post-count probe-to-gene sum should be four, found {gene_sum}"
            )

        # A barcode one base from two CBC whitelist entries is rejected.
        ambiguous = root / "ambiguous.onlist"
        ambiguous.write_text(f"AAAAAAAAAAAAAAAA\t{bc001}\nCAAAAAAAAAAAAAAA\t-\n")
        query_cbc = "GAAAAAAAAAAAAAAA"
        r1, r2 = write_fastq_pair(
            root,
            "ambiguous",
            [(query_cbc, "ACGTACGTACGT", probe_a, bc001)],
        )
        ambiguous_map = root / "ambiguous-map"
        run(
            [
                kallisto,
                "bus",
                "-i",
                index,
                "-o",
                ambiguous_map,
                "-x",
                TECHNOLOGY,
                "-t",
                1,
                r1,
                r2,
            ],
            stdout=subprocess.DEVNULL,
        )
        ambiguous_bus = root / "ambiguous-corrected.bus"
        run(
            [
                bustools,
                "correct",
                "-w",
                ambiguous,
                "-o",
                ambiguous_bus,
                ambiguous_map / "output.bus",
            ],
            stdout=subprocess.DEVNULL,
        )
        ambiguous_rows = bus_text(bustools, ambiguous_bus, root / "ambiguous.tsv")
        if ambiguous_rows:
            raise AssertionError(
                "Bustools must reject a CBC with two equally close whitelist entries"
            )

    print(
        json.dumps(
            {
                "versions": "kb_python 0.30.2 / kallisto 0.52.0 / bustools 0.45.1",
                "technology": TECHNOLOGY,
                "probe_pair_gene": gene_id,
                "probe_a_umi_count_bc001": observed_counts[(probe_a_id, bc001_cell)],
                "probe_b_umi_count_bc001": observed_counts[(probe_b_id, bc001_cell)],
                "probe_gene_sum_bc001": gene_sum,
                "bc001_bc002_same_cbc_distinct": len(expected_barcodes) == 2,
                "cross_chunk_duplicate_collapsed": True,
                "single_base_cbc_and_bc_corrected": True,
                "ambiguous_cbc_rejected": True,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    test_native()
