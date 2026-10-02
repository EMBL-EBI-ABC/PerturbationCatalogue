#!/usr/bin/env python3
"""Count one physical Flex/Ultima lane with pinned Kallisto and Bustools inputs."""

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as scipy_io
from scipy import sparse

from bam_to_fastq import download_ranges
from stream_count import download_sra


BUS_HEADER = struct.Struct("<4sIIII")
BUS_RECORD = np.dtype(
    [
        ("barcode", "<u8"),
        ("umi", "<u8"),
        ("ec", "<u4"),
        ("count", "<u4"),
        ("flags", "<u4"),
        ("pad", "<u4"),
    ],
    align=False,
)
BUS_CHUNK_RECORDS = 1_000_000
BC_COMPONENTS = "0,0,16,1,68,76:0,16,28:1,0,50"
GUIDE_COMPONENTS = "0,0,16,1,0,8:0,16,28:1,8,0"
ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
ALIAS_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}")
DNA = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def file_hash(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_rows(path):
    return [
        line.rstrip("\r\n")
        for line in Path(path).open(encoding="utf-8")
        if line.strip()
    ]


def read_two_columns(path):
    rows = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise ValueError(
                    f"{path}:{line_number}: expected two non-empty columns"
                )
            rows.append((fields[0], fields[1]))
    return rows


def run_logged(command, log_path, timeout=24 * 3600):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as log:
        result = subprocess.run(
            [str(part) for part in command],
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(
            f"Command exited {result.returncode}: {command[0]} (see {log_path})"
        )
    return log_path.read_text(errors="replace")


def encode_dna(sequence):
    value = 0
    for base in sequence:
        try:
            value = (value << 2) | DNA[base]
        except KeyError as error:
            raise ValueError(f"Invalid DNA barcode: {sequence!r}") from error
    return value


def make_component_onlist(first, second, output):
    """Write Bustools' independent component lists without pairing unrelated rows."""
    if (
        not first
        or not second
        or len(set(first)) != len(first)
        or len(set(second)) != len(second)
    ):
        raise ValueError("Barcode component lists must be non-empty and unique")
    if any(not re.fullmatch(r"[ACGT]+", value) for value in (*first, *second)):
        raise ValueError("Barcode component lists must contain only A/C/G/T")
    with Path(output).open("w", encoding="ascii", newline="\n") as handle:
        handle.write(first[0] + "\t" + second[0] + "\n")
        for value in first[1:]:
            handle.write(value + "\t-\n")
        for value in second[1:]:
            handle.write("-\t" + value + "\n")


def parse_ec_map(path, target_count):
    rows = {}
    with Path(path).open(encoding="ascii") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed equivalence class at {path}:{line_number}")
            try:
                ec_id = int(fields[0])
                targets = tuple(sorted({int(value) for value in fields[1].split(",")}))
            except ValueError as error:
                raise ValueError(
                    f"Invalid equivalence class at {path}:{line_number}"
                ) from error
            if (
                ec_id < 0
                or not targets
                or targets[0] < 0
                or targets[-1] >= target_count
            ):
                raise ValueError(
                    f"Equivalence class is outside the target list in {path}"
                )
            if ec_id in rows:
                raise ValueError(f"Duplicate equivalence class ID {ec_id} in {path}")
            rows[ec_id] = targets
    if rows and set(rows) != set(range(len(rows))):
        raise ValueError(f"Equivalence class IDs are not contiguous in {path}")
    return rows


def bus_header(path):
    with Path(path).open("rb") as handle:
        prefix = handle.read(BUS_HEADER.size)
        if len(prefix) != BUS_HEADER.size:
            raise ValueError(f"Truncated BUS header: {path}")
        magic, version, barcode_length, umi_length, text_length = BUS_HEADER.unpack(
            prefix
        )
        if magic != b"BUS\0":
            raise ValueError(f"Invalid BUS magic in {path}")
        text = handle.read(text_length)
        if len(text) != text_length:
            raise ValueError(f"Truncated BUS header text in {path}")
    offset = BUS_HEADER.size + text_length
    size = Path(path).stat().st_size
    if size < offset or (size - offset) % BUS_RECORD.itemsize:
        raise ValueError(f"BUS records are truncated or malformed: {path}")
    return {
        "bytes": prefix + text,
        "version": version,
        "barcode_length": barcode_length,
        "umi_length": umi_length,
        "text_length": text_length,
        "record_offset": offset,
        "records": (size - offset) // BUS_RECORD.itemsize,
    }


def count_sorted_barcodes(path):
    """Count composite barcodes in a Bustools-sorted BUS without retaining them all."""
    info = bus_header(path)
    groups = 0
    previous = None
    with Path(path).open("rb") as handle:
        handle.seek(info["record_offset"])
        while True:
            block = handle.read(BUS_RECORD.itemsize * BUS_CHUNK_RECORDS)
            if not block:
                break
            if len(block) % BUS_RECORD.itemsize:
                raise ValueError(f"Partial BUS record in {path}")
            barcodes = np.frombuffer(block, dtype=BUS_RECORD)["barcode"]
            if not len(barcodes):
                continue
            if np.any(barcodes[1:] < barcodes[:-1]) or (
                previous is not None and int(barcodes[0]) < previous
            ):
                raise ValueError(f"BUS barcode order is not sorted in {path}")
            groups += int(np.count_nonzero(barcodes[1:] != barcodes[:-1]))
            if previous is None or int(barcodes[0]) != previous:
                groups += 1
            previous = int(barcodes[-1])
    return groups


class AliasBusWriter:
    """Append canonicalized BUS records to one file per barcode suffix."""

    def __init__(self, directory, aliases, sequence_by_alias):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.aliases = list(aliases)
        if not self.aliases:
            raise ValueError("At least one alias is required")
        codes = [encode_dna(sequence_by_alias[alias]) for alias in self.aliases]
        if len(set(codes)) != len(codes):
            raise ValueError("Barcode aliases do not have unique sequences")
        self.paths = {alias: self.directory / f"{alias}.bus" for alias in self.aliases}
        self.handles = {}
        self.codes = np.full(1 << 16, -1, dtype=np.int16)
        for index, code in enumerate(codes):
            if code >= len(self.codes):
                raise ValueError("Alias barcode segment must be eight bases")
            self.codes[code] = index
        self.records_by_alias = {alias: 0 for alias in self.aliases}
        self.header = None

    def append(self, records, header):
        if self.header is None:
            self.header = header
        elif (
            self.header["bytes"] != header["bytes"]
            or self.header["version"] != header["version"]
            or self.header["barcode_length"] != header["barcode_length"]
            or self.header["umi_length"] != header["umi_length"]
            or self.header["text_length"] != header["text_length"]
        ):
            raise ValueError("BUS header changed between sources in one lane")
        suffixes = (records["barcode"] & 0xFFFF).astype(np.uint16, copy=False)
        alias_ids = self.codes[suffixes]
        if np.any(alias_ids < 0):
            bad = np.unique(suffixes[alias_ids < 0])[:8].tolist()
            raise ValueError(
                f"Corrected BUS has unrecognized alias suffix codes: {bad}"
            )
        if not len(records):
            return
        order = np.argsort(alias_ids, kind="stable")
        ordered = records[order]
        ordered_ids = alias_ids[order]
        boundaries = np.searchsorted(ordered_ids, np.arange(len(self.aliases) + 1))
        for index, alias in enumerate(self.aliases):
            start, end = int(boundaries[index]), int(boundaries[index + 1])
            if start == end:
                continue
            handle = self.handles.get(alias)
            if handle is None:
                path = self.paths[alias]
                handle = path.open("wb")
                handle.write(header["bytes"])
                self.handles[alias] = handle
            ordered[start:end].tofile(handle)
            self.records_by_alias[alias] += end - start

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


class CanonicalEC:
    """Remap per-run equivalence classes to one lane-wide target-set table."""

    def __init__(
        self, target_ids, modality, directory, alias_writer, gene_by_target=None
    ):
        self.target_ids = list(target_ids)
        self.target_count = len(self.target_ids)
        self.modality = modality
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.alias_writer = alias_writer
        self.gene_by_target = gene_by_target or {}
        self.ec_id_by_targets = {}
        self.targets_by_ec_id = []
        self.metrics = {
            "bus_records": 0,
            "unique_target_reads": 0,
            "same_gene_multi_probe_reads": 0,
            "cross_gene_ambiguous_reads": 0,
            "multi_feature_ambiguous_reads": 0,
            "bus_read_count": 0,
        }
        self.transcript_sha256 = None

    def add_source(self, bus_path, ec_path, transcript_path):
        transcript_path = Path(transcript_path)
        transcript_sha = file_hash(transcript_path)
        transcripts = file_rows(transcript_path)
        if transcripts != self.target_ids:
            raise ValueError(
                f"{self.modality} index transcript order differs from its pinned feature map"
            )
        if self.transcript_sha256 is None:
            self.transcript_sha256 = transcript_sha
        elif transcript_sha != self.transcript_sha256:
            raise ValueError(
                f"{self.modality} index transcript list changed between sources"
            )

        local_ecs = parse_ec_map(ec_path, self.target_count)
        local_to_global = np.empty(len(local_ecs), dtype=np.uint32)
        for local_id, targets in local_ecs.items():
            canonical_id = self.ec_id_by_targets.get(targets)
            if canonical_id is None:
                canonical_id = len(self.targets_by_ec_id)
                if canonical_id >= np.iinfo(np.uint32).max:
                    raise ValueError("Too many equivalence classes for BUS")
                self.ec_id_by_targets[targets] = canonical_id
                self.targets_by_ec_id.append(targets)
            local_to_global[local_id] = canonical_id

        info = bus_header(bus_path)
        if (
            info["barcode_length"] != 24
            or info["umi_length"] != 12
            or info["version"] != 1
        ):
            raise ValueError(
                f"Unexpected {self.modality} BUS geometry: "
                f"version={info['version']} barcode={info['barcode_length']} UMI={info['umi_length']}"
            )
        if self.alias_writer.header is not None and (
            self.alias_writer.header["bytes"] != info["bytes"]
            or self.alias_writer.header["version"] != info["version"]
            or self.alias_writer.header["barcode_length"] != info["barcode_length"]
            or self.alias_writer.header["umi_length"] != info["umi_length"]
            or self.alias_writer.header["text_length"] != info["text_length"]
        ):
            raise ValueError("BUS header changed between sources in one lane")

        read_counts = {"bus_records": 0, "bus_read_count": 0}
        with Path(bus_path).open("rb") as handle:
            handle.seek(info["record_offset"])
            while True:
                block = handle.read(BUS_RECORD.itemsize * BUS_CHUNK_RECORDS)
                if not block:
                    break
                if len(block) % BUS_RECORD.itemsize:
                    raise ValueError(f"Partial BUS record in {bus_path}")
                records = np.frombuffer(block, dtype=BUS_RECORD).copy()
                old_ids = records["ec"].astype(np.uint64, copy=False)
                if old_ids.size and int(old_ids.max()) >= len(local_to_global):
                    raise ValueError(
                        f"BUS equivalence class ID is absent from {ec_path}"
                    )
                global_ids = local_to_global[old_ids]
                ec_weights = np.bincount(
                    old_ids.astype(np.int64, copy=False),
                    weights=records["count"].astype(np.int64, copy=False),
                    minlength=len(local_to_global),
                )
                for old_id, weight in enumerate(ec_weights):
                    if not weight:
                        continue
                    targets = local_ecs[old_id]
                    if len(targets) == 1:
                        self.metrics["unique_target_reads"] += int(weight)
                    elif self.modality == "gex":
                        genes = {
                            self.gene_by_target[self.target_ids[index]]
                            for index in targets
                        }
                        if len(genes) == 1:
                            self.metrics["same_gene_multi_probe_reads"] += int(weight)
                        else:
                            self.metrics["cross_gene_ambiguous_reads"] += int(weight)
                    else:
                        self.metrics["multi_feature_ambiguous_reads"] += int(weight)
                records["ec"] = global_ids
                self.alias_writer.append(records, info)
                read_counts["bus_records"] += len(records)
                read_counts["bus_read_count"] += int(records["count"].sum())
        if read_counts["bus_records"] != info["records"]:
            raise ValueError(
                f"Read an unexpected number of BUS records from {bus_path}"
            )
        self.metrics["bus_records"] += read_counts["bus_records"]
        self.metrics["bus_read_count"] += read_counts["bus_read_count"]

    def write_reference(self):
        if not self.transcript_sha256:
            raise ValueError(f"No {self.modality} Kallisto runs were completed")
        ec_path = self.directory / "matrix.ec"
        with ec_path.open("w", encoding="ascii", newline="\n") as handle:
            for ec_id, targets in enumerate(self.targets_by_ec_id):
                handle.write(str(ec_id) + "\t" + ",".join(map(str, targets)) + "\n")
        tx_path = self.directory / "transcripts.txt"
        tx_path.write_text("\n".join(self.target_ids) + "\n", encoding="utf-8")
        self.alias_writer.close()
        return ec_path, tx_path


def correction_summary(log_path):
    text = Path(log_path).read_text(errors="replace")
    result = {}
    for label, key in (
        ("Processed", "processed"),
        ("In on-list", "in_onlist"),
        ("Corrected", "corrected"),
        ("Uncorrected", "uncorrected"),
        ("Replaced", "replaced"),
        ("Not replaced", "not_replaced"),
    ):
        if key == "processed":
            pattern = r"^\s*Processed\s*(?:=\s*)?([\d,]+)(?:\s+BUS records)?\s*$"
        else:
            pattern = rf"^\s*{re.escape(label)}\s*=\s*([\d,]+)\s*$"
        match = re.search(pattern, text, re.M)
        if match:
            result[key] = int(match.group(1).replace(",", ""))
    return result


def parse_allowlist_metrics(log_path):
    text = Path(log_path).read_text(errors="replace")
    match = re.search(
        r"Read in ([\d,]+) BUS records, wrote ([\d,]+) barcodes to on-list with threshold ([\d,]+)",
        text,
    )
    if not match:
        raise ValueError(f"Could not read Bustools allowlist summary from {log_path}")
    return {
        "bus_records": int(match.group(1).replace(",", "")),
        "whitelist_barcodes": int(match.group(2).replace(",", "")),
        "threshold": int(match.group(3).replace(",", "")),
        "policy": "bustools allowlist default threshold per lane x BC alias",
    }


def load_manifest_features(args):
    probes = []
    gene_name_by_id = {}
    with Path(args.gene_probes).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"Malformed gene probe row {line_number}")
            gene_id, gene_name, sequence = fields
            if not re.fullmatch(r"ENSG\d+", gene_id) or not gene_name:
                raise ValueError(f"Invalid gene feature row {line_number}")
            if not re.fullmatch(r"[ACGT]+", sequence):
                raise ValueError(f"Invalid probe sequence on row {line_number}")
            if gene_id in gene_name_by_id and gene_name_by_id[gene_id] != gene_name:
                raise ValueError(f"Conflicting gene names for {gene_id}")
            gene_name_by_id[gene_id] = gene_name
            probes.append((gene_id, gene_name, sequence))
    probe_t2g = read_two_columns(args.gex_t2g)
    probe_to_gene_rows = read_two_columns(args.gex_probe_to_gene)
    if [source for source, _ in probe_t2g] != [
        source for source, _ in probe_to_gene_rows
    ]:
        raise ValueError("GEX self-T2G and probe-to-gene target orders differ")
    probe_ids = [source for source, _ in probe_t2g]
    if len(probe_ids) != len(set(probe_ids)):
        raise ValueError("GEX probe IDs are not unique")
    if any(source != target for source, target in probe_t2g):
        raise ValueError("GEX BUS T2G must map each probe target to itself")
    probe_to_gene = dict(probe_to_gene_rows)
    if set(probe_to_gene.values()) - set(gene_name_by_id):
        raise ValueError("Probe-to-gene table contains an unknown stable gene ID")
    if len(probes) != len(probe_ids):
        raise ValueError("Probe feature table and Kallisto reference sizes differ")
    gene_count_rows = [
        line.split("\t", 1)[0] for line in file_rows(args.gene_count_features)
    ]
    if gene_count_rows != [probe_to_gene[probe] for probe in probe_ids]:
        raise ValueError(
            "Legacy gene feature table disagrees with probe-to-gene reference"
        )
    gene_ids = list(dict.fromkeys(probe_to_gene[probe] for probe in probe_ids))
    gene_names = [gene_name_by_id[gene_id] for gene_id in gene_ids]

    guide_t2g = read_two_columns(args.guide_t2g)
    if any(source != target for source, target in guide_t2g):
        raise ValueError("Guide BUS T2G must map each guide target to itself")
    guide_ids = [source for source, _ in guide_t2g]
    if len(guide_ids) != len(set(guide_ids)):
        raise ValueError("Guide target IDs are not unique")
    guide_feature_ids = [
        line.split("\t", 1)[0] for line in file_rows(args.guide_features)
    ]
    guide_count_ids = [
        line.split("\t", 1)[0] for line in file_rows(args.guide_count_features)
    ]
    guide_targets = read_tsv(args.guide_targets)
    guide_by_id = {row["guide_id"]: row for row in guide_targets}
    if (
        guide_feature_ids != guide_ids
        or guide_count_ids != guide_ids
        or set(guide_by_id) != set(guide_ids)
    ):
        raise ValueError("Guide references, count features and target labels differ")
    label_by_guide = {}
    for guide_id in guide_ids:
        row = guide_by_id[guide_id]
        symbol = row["target_symbol"].strip()
        ensg = row["target_ensg"].strip()
        label_by_guide[guide_id] = (
            f"{symbol}_{guide_id}__{ensg}" if ensg else f"{symbol}_{guide_id}"
        )
    return {
        "probe_ids": probe_ids,
        "probe_to_gene": probe_to_gene,
        "gene_ids": gene_ids,
        "gene_names": gene_names,
        "gene_name_by_id": gene_name_by_id,
        "guide_ids": guide_ids,
        "guide_labels": [label_by_guide[guide_id] for guide_id in guide_ids],
    }


def load_pool_mapping(args):
    sample_rows = [
        row for row in read_tsv(args.samples) if row["pool_id"] == args.pool_id
    ]
    samples_by_bc = {}
    sample_pairs = {}
    seen_cr = set()
    expected_lanes = set()
    for sample in sample_rows:
        if not ID_RE.fullmatch(sample["sample_id"]):
            raise ValueError(f"Unsafe sample ID: {sample['sample_id']!r}")
        lanes = sample["pool_lanes"].split(";")
        if any(not ID_RE.fullmatch(lane) for lane in lanes) or len(lanes) != len(
            set(lanes)
        ):
            raise ValueError(f"Invalid lane list in {sample['sample_id']}")
        expected_lanes.update(lanes)
        pairs = [pair.split(":") for pair in sample["bc_cr_pairs"].split(";") if pair]
        if not pairs:
            raise ValueError(f"Sample {sample['sample_id']} has no BC/CR aliases")
        for pair in pairs:
            if len(pair) != 2 or not all(ALIAS_RE.fullmatch(value) for value in pair):
                raise ValueError(f"Invalid BC/CR alias pair in {sample['sample_id']}")
            bc_alias, cr_alias = pair
            if bc_alias in samples_by_bc or cr_alias in seen_cr:
                raise ValueError("Pool sample aliases are not one-to-one")
            seen_cr.add(cr_alias)
            samples_by_bc[bc_alias] = {
                "sample_id": sample["sample_id"],
                "cr_alias": cr_alias,
            }
            sample_pairs[bc_alias] = cr_alias
    if not sample_rows or not expected_lanes:
        raise ValueError(f"No sample/lane metadata for pool {args.pool_id}")
    if args.lane_id not in expected_lanes:
        raise ValueError(f"Lane {args.lane_id} is not declared by pool {args.pool_id}")
    aliases = read_tsv(args.barcode_aliases)
    bc_sequences = {row["bc_alias"]: row["bc_sequence"] for row in aliases}
    cr_sequences = {row["cr_alias"]: row["cr_sequence"] for row in aliases}
    if (
        len(bc_sequences) != len(aliases)
        or len(cr_sequences) != len(aliases)
        or len(set(bc_sequences.values())) != len(aliases)
        or len(set(cr_sequences.values())) != len(aliases)
        or set(samples_by_bc) != set(bc_sequences)
        or set(sample_pairs.values()) != set(cr_sequences)
    ):
        raise ValueError("Sample BC/CR mapping and probe barcode reference disagree")
    for alias, sequence in (*bc_sequences.items(), *cr_sequences.items()):
        if not ALIAS_RE.fullmatch(alias) or not re.fullmatch(r"[ACGT]{8}", sequence):
            raise ValueError(f"Invalid alias or sequence: {alias}={sequence}")
    return samples_by_bc, bc_sequences, cr_sequences, sample_pairs


def load_bc_variants(path, bc_sequences):
    variants = read_tsv(path)
    required = {"raw_sequence", "canonical_sequence", "alias"}
    if not variants or not required.issubset(variants[0]):
        raise ValueError(
            "BC variant table must contain raw_sequence, canonical_sequence, alias"
        )
    by_raw = {}
    for row in variants:
        raw, canonical, alias = (
            row["raw_sequence"],
            row["canonical_sequence"],
            row["alias"],
        )
        if not re.fullmatch(r"[ACGT]{8}", raw) or bc_sequences.get(alias) != canonical:
            raise ValueError(f"Invalid BC variant mapping: {row}")
        if raw in by_raw:
            raise ValueError(f"Ambiguous raw BC sequence: {raw}")
        by_raw[raw] = (canonical, alias)
    expected_aliases = set(bc_sequences)
    if {alias for _, alias in by_raw.values()} != expected_aliases:
        raise ValueError("BC variants do not cover every canonical alias")
    return variants, by_raw


def verify_cached_file(path, expected_bytes, expected_md5):
    path = Path(path)
    if (
        path.stat().st_size != int(expected_bytes)
        or file_hash(path, "md5") != expected_md5.lower()
    ):
        raise RuntimeError(f"Cached SRA differs from its pinned size/MD5: {path}")


def read_count_matrix(prefix, expected_features):
    prefix = Path(prefix)
    barcodes = file_rows(str(prefix) + ".barcodes.txt")
    features = file_rows(str(prefix) + ".genes.txt")
    if features != expected_features:
        raise ValueError(f"Bustools feature order differs from reference: {prefix}")
    matrix = scipy_io.mmread(str(prefix) + ".mtx").tocsr()
    if matrix.shape != (len(barcodes), len(features)):
        raise ValueError(
            f"Bustools matrix shape does not match barcodes/features: {prefix}"
        )
    if matrix.nnz and (
        not np.all(np.isfinite(matrix.data))
        or np.any(matrix.data < 0)
        or not np.all(matrix.data == np.floor(matrix.data))
        or np.any(matrix.data > np.iinfo(np.int32).max)
    ):
        raise ValueError(f"Bustools matrix contains invalid counts: {prefix}")
    matrix.data = matrix.data.astype(np.int32, copy=False)
    return barcodes, matrix


def aggregate_probe_matrix(matrix, probe_ids, probe_to_gene, gene_ids):
    if matrix.shape[1] != len(probe_ids):
        raise ValueError("Probe matrix columns do not match the probe reference")
    gene_index = {gene_id: index for index, gene_id in enumerate(gene_ids)}
    target_gene_rows = np.fromiter(
        (gene_index[probe_to_gene[probe]] for probe in probe_ids),
        dtype=np.int32,
        count=len(probe_ids),
    )
    coo = matrix.tocoo()
    result = sparse.coo_matrix(
        (coo.data, (coo.row, target_gene_rows[coo.col])),
        shape=(matrix.shape[0], len(gene_ids)),
        dtype=np.int32,
    ).tocsr()
    result.sum_duplicates()
    result.eliminate_zeros()
    return result


def write_alias_h5ad(
    output,
    lane_id,
    bc_alias,
    cr_alias,
    sample_id,
    guide_status,
    gene_ids,
    gene_names,
    guide_labels,
    gex_barcodes,
    gex_matrix,
    guide_matrix,
    lane_metrics,
):
    cell_ids = [f"{barcode}-1_{lane_id}" for barcode in gex_barcodes]
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError(f"Duplicate composite cells in {lane_id}/{bc_alias}")
    obs = pd.DataFrame(
        {
            "sample_id": pd.Series(
                [sample_id] * len(cell_ids), index=cell_ids, dtype=object
            ),
            "lane_id": pd.Series(
                [lane_id] * len(cell_ids), index=cell_ids, dtype=object
            ),
            "bc_alias": pd.Series(
                [bc_alias] * len(cell_ids), index=cell_ids, dtype=object
            ),
            "cr_alias": pd.Series(
                [cr_alias] * len(cell_ids), index=cell_ids, dtype=object
            ),
            "guide_coverage_status": pd.Series(
                [guide_status] * len(cell_ids), index=cell_ids, dtype=object
            ),
        },
        index=pd.Index(cell_ids, dtype=object),
    )
    var = pd.DataFrame(
        {"gene_name": gene_names},
        index=pd.Index(gene_ids, name="gene_id", dtype=object),
    )
    if guide_matrix.shape != (len(cell_ids), len(guide_labels)):
        raise ValueError("Guide matrix shape does not match the alias cells/features")
    adata = ad.AnnData(X=gex_matrix, obs=obs, var=var)
    adata.obsm["guides"] = guide_matrix
    adata.uns["guide_names"] = list(guide_labels)
    adata.uns["flex_counter_version"] = lane_metrics["counter_version"]
    adata.uns["flex_counter_version_label"] = lane_metrics["counter_version_label"]
    adata.uns["flex_umi_policy"] = (
        "Bustools --umi-gene deduplication on unique probe/guide target IDs; "
        "GEX probe counts are summed to stable ENSG after probe-level counting."
    )
    adata.uns["flex_mapping_policy"] = (
        "Kallisto BUS pseudoalignment; multi-target equivalence classes are not "
        "allocated with an EM or multimapping option."
    )
    adata.uns["flex_sample_id"] = sample_id
    adata.uns["flex_lane_id"] = lane_id
    adata.uns["flex_bc_alias"] = bc_alias
    adata.uns["flex_cr_alias"] = cr_alias
    adata.uns["flex_guide_coverage_status"] = guide_status
    adata.uns["flex_reference_hashes_json"] = json.dumps(
        lane_metrics["reference_sha256"], sort_keys=True
    )
    adata.uns["flex_lane_metrics_json"] = json.dumps(
        {
            "pool_id": lane_metrics["pool_id"],
            "lane_id": lane_id,
            "bc_alias": bc_alias,
            "cr_alias": cr_alias,
            "guide_coverage_status": guide_status,
            "barcode_metrics": lane_metrics["barcode_counts"].get(bc_alias, {}),
            "gex_mapping": lane_metrics["gex_mapping"],
            "guide_mapping": lane_metrics["guide_mapping"],
            "equivalence_class_metrics_scope": lane_metrics[
                "equivalence_class_metrics_scope"
            ],
        },
        sort_keys=True,
    )
    output = Path(output)
    temporary = output.with_suffix(output.suffix + ".partial")
    adata.write_h5ad(temporary, compression="gzip")
    temporary.replace(output)


def join_guide_counts(gex_barcodes, guide_data, guide_count):
    joined = sparse.csr_matrix((len(gex_barcodes), guide_count), dtype=np.int32)
    if guide_data is None or not gex_barcodes:
        return joined, 0
    guide_barcodes, guide_counts = guide_data
    target_row = {barcode: row for row, barcode in enumerate(gex_barcodes)}
    pairs = [
        (source_row, target_row[barcode])
        for source_row, barcode in enumerate(guide_barcodes)
        if barcode in target_row
    ]
    if not pairs:
        return joined, 0
    source_rows = np.fromiter((pair[0] for pair in pairs), dtype=np.int32)
    target_rows = np.fromiter((pair[1] for pair in pairs), dtype=np.int32)
    subset = guide_counts[source_rows].tocoo()
    joined = sparse.coo_matrix(
        (subset.data, (target_rows[subset.row], subset.col)),
        shape=(len(gex_barcodes), guide_count),
        dtype=np.int32,
    ).tocsr()
    joined.sum_duplicates()
    return joined, len(pairs)


def run_bus_source(
    source_label,
    r1,
    r2,
    modality,
    index,
    technology,
    ec_counter,
    onlist,
    replacement,
    bustools,
    work,
    logs,
    threads,
    kallisto,
):
    chunk_dir = Path(work) / f"{modality}_{source_label}"
    corrected = Path(work) / f"{modality}_{source_label}.corrected.bus"
    replaced = Path(work) / f"{modality}_{source_label}.replaced.bus"
    chunk_dir.mkdir()
    log = logs / f"{modality}_{source_label}.kallisto.log"
    started = time.monotonic()
    run_logged(
        [
            kallisto,
            "bus",
            "-i",
            index,
            "-x",
            technology,
            "-o",
            chunk_dir,
            "-t",
            threads,
            "--verbose",
            r1,
            r2,
        ],
        log,
    )
    run_info_path = chunk_dir / "run_info.json"
    if not run_info_path.is_file():
        raise RuntimeError(f"Kallisto did not create run_info.json for {source_label}")
    run_info = json.loads(run_info_path.read_text())
    if run_info.get("kallisto_version") != "0.52.0":
        raise RuntimeError(f"Unexpected Kallisto version in {source_label}")
    local_bus = chunk_dir / "output.bus"
    if not local_bus.is_file():
        raise RuntimeError(f"Kallisto did not create output.bus for {source_label}")
    raw_bus_records = bus_header(local_bus)["records"]
    if raw_bus_records != run_info.get("n_pseudoaligned"):
        raise ValueError(
            f"Kallisto BUS records differ from its mapping summary for {source_label}"
        )

    correction_log = logs / f"{modality}_{source_label}.correct.log"
    run_logged(
        [bustools, "correct", "-w", onlist, "-o", corrected, local_bus],
        correction_log,
    )
    final_bus = corrected
    replacement_metrics = {}
    if replacement:
        replace_log = logs / f"{modality}_{source_label}.replace.log"
        run_logged(
            [bustools, "correct", "-r", "-w", replacement, "-o", replaced, corrected],
            replace_log,
        )
        replacement_metrics = correction_summary(replace_log)
        final_bus = replaced
    ec_counter.add_source(
        final_bus, chunk_dir / "matrix.ec", chunk_dir / "transcripts.txt"
    )
    final_bus.unlink(missing_ok=True)
    corrected.unlink(missing_ok=True)
    shutil.rmtree(chunk_dir)
    return {
        "source": source_label,
        "raw_bus_records": raw_bus_records,
        "kallisto": {
            key: run_info.get(key)
            for key in (
                "n_processed",
                "n_pseudoaligned",
                "n_unique",
                "kallisto_version",
                "index_version",
                "k-mer length",
            )
        },
        "barcode_correction": correction_summary(correction_log),
        "alias_replacement": replacement_metrics,
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }


def filter_and_count_alias(
    bus_path,
    ec_path,
    transcript_path,
    t2g_path,
    modality,
    alias,
    work,
    logs,
    bustools,
    threads,
):
    bus_path = Path(bus_path)
    if not bus_path.is_file() or bus_header(bus_path)["records"] == 0:
        return None, {"input_bus_records": 0, "retained_cells": 0}
    prefix = Path(work) / f"count_{modality}_{alias}"
    temp = Path(work) / f"sort_{modality}_{alias}"
    temp.mkdir()
    first_sorted = Path(work) / f"{modality}_{alias}.sorted.bus"
    run_logged(
        [
            bustools,
            "sort",
            "-t",
            threads,
            "-m",
            "2G",
            "-T",
            temp,
            "-o",
            first_sorted,
            bus_path,
        ],
        logs / f"{modality}_{alias}.sort.log",
    )
    observed_barcodes = count_sorted_barcodes(first_sorted)
    whitelist_path = Path(work) / f"{modality}_{alias}.allowlist.txt"
    allow_log = logs / f"{modality}_{alias}.allowlist.log"
    if modality == "gex":
        run_logged(
            [bustools, "allowlist", "-o", whitelist_path, first_sorted],
            allow_log,
        )
        allow_metrics = parse_allowlist_metrics(allow_log)
        if allow_metrics["whitelist_barcodes"] == 0:
            first_sorted.unlink()
            shutil.rmtree(temp)
            return None, {
                "input_bus_records": bus_header(bus_path)["records"],
                "observed_composite_barcodes": observed_barcodes,
                **allow_metrics,
                "retained_cells": 0,
            }
        filtered = Path(work) / f"{modality}_{alias}.filtered.bus"
        run_logged(
            [bustools, "correct", "-w", whitelist_path, "-o", filtered, first_sorted],
            logs / f"{modality}_{alias}.filter_correct.log",
        )
    else:
        allow_metrics = {
            "policy": "guide counts are joined to GEX-filtered cells; no independent guide barcode filter"
        }
        filtered = first_sorted
    final_sorted = Path(work) / f"{modality}_{alias}.count.sorted.bus"
    run_logged(
        [
            bustools,
            "sort",
            "-t",
            threads,
            "-m",
            "2G",
            "-T",
            temp,
            "-o",
            final_sorted,
            filtered,
        ],
        logs / f"{modality}_{alias}.count_sort.log",
    )
    if filtered != first_sorted:
        filtered.unlink(missing_ok=True)
    first_sorted.unlink(missing_ok=True)
    shutil.rmtree(temp)
    run_logged(
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
            transcript_path,
            "--genecounts",
            "--umi-gene",
            final_sorted,
        ],
        logs / f"{modality}_{alias}.count.log",
    )
    final_sorted.unlink(missing_ok=True)
    if not Path(str(prefix) + ".mtx").is_file():
        raise RuntimeError(f"Bustools count did not create a matrix for {alias}")
    return prefix, {
        "input_bus_records": bus_header(bus_path)["records"],
        "observed_composite_barcodes": observed_barcodes,
        **allow_metrics,
    }


def execute(args):
    if not ID_RE.fullmatch(args.pool_id) or not ID_RE.fullmatch(args.lane_id):
        raise ValueError("Pool and lane IDs must be safe non-empty identifiers")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_counts = outdir / "lane_counts"
    if output_counts.exists() and any(output_counts.iterdir()):
        raise FileExistsError(
            f"Lane count output directory is not empty: {output_counts}"
        )
    output_counts.mkdir(exist_ok=True)
    logs = outdir / "logs"
    work = outdir / "work"
    downloads = work / "downloads"
    fastqs = work / "fastqs"
    bus_chunks = work / "bus_chunks"
    for path in (logs, downloads, fastqs, bus_chunks):
        path.mkdir(parents=True, exist_ok=True)

    samples_by_bc, bc_sequences, cr_sequences, sample_pairs = load_pool_mapping(args)
    bc_variants, variants_by_raw = load_bc_variants(
        args.bc_barcode_variants, bc_sequences
    )
    reference = load_manifest_features(args)
    pool_gex = [
        row
        for row in read_tsv(args.gex_sources)
        if row["pool_id"] == args.pool_id and row["lane_id"] == args.lane_id
    ]
    pool_guides = [
        row
        for row in read_tsv(args.guide_sources)
        if row["pool_id"] == args.pool_id and row["lane_id"] == args.lane_id
    ]
    coverage_rows = [
        row
        for row in read_tsv(args.guide_coverage)
        if row["pool_id"] == args.pool_id and row["lane_id"] == args.lane_id
    ]
    if not pool_gex:
        raise ValueError(f"No pinned GEX sources for {args.pool_id}/{args.lane_id}")
    if len(coverage_rows) != 1:
        raise ValueError(f"Expected one guide coverage row for {args.lane_id}")
    coverage = coverage_rows[0]
    guide_status = coverage["guide_status"]
    if guide_status not in {"available", "partial", "no_archived_guide_sra"}:
        raise ValueError(f"Unsupported guide coverage status: {guide_status}")
    if (
        int(coverage["gex_source_pairs"]) != len(pool_gex)
        or int(coverage["guide_sra_runs"]) != len(pool_guides)
        or ((guide_status == "no_archived_guide_sra") != (not pool_guides))
    ):
        raise ValueError("Source manifests disagree with guide coverage metadata")
    if any(row["lane_id"] != args.lane_id for row in (*pool_gex, *pool_guides)):
        raise ValueError("A source row crossed the selected physical lane boundary")

    cbc_whitelist = file_rows(args.cbc_whitelist)
    if (
        not cbc_whitelist
        or len(set(cbc_whitelist)) != len(cbc_whitelist)
        or any(not re.fullmatch(r"[ACGT]{16}", value) for value in cbc_whitelist)
    ):
        raise ValueError("CBC whitelist must contain unique 16-base DNA barcodes")
    variants_dir = work / "barcodes"
    variants_dir.mkdir()
    bc_onlist = variants_dir / "cbc_bc_onlist.tsv"
    cr_onlist = variants_dir / "cbc_cr_onlist.tsv"
    bc_replace = variants_dir / "bc_to_canonical.tsv"
    cr_to_bc_replace = variants_dir / "cr_to_bc.tsv"
    make_component_onlist(cbc_whitelist, sorted(variants_by_raw), bc_onlist)
    make_component_onlist(cbc_whitelist, list(cr_sequences.values()), cr_onlist)
    with bc_replace.open("w", encoding="ascii", newline="\n") as handle:
        for row in bc_variants:
            handle.write(row["raw_sequence"] + "\t*" + row["canonical_sequence"] + "\n")
    with cr_to_bc_replace.open("w", encoding="ascii", newline="\n") as handle:
        for bc_alias, cr_alias in sample_pairs.items():
            handle.write(cr_sequences[cr_alias] + "\t*" + bc_sequences[bc_alias] + "\n")

    import kb_python.utils

    kallisto = str(kb_python.utils.get_kallisto_binary_path())
    bustools = str(kb_python.utils.get_bustools_binary_path())
    kb_python_version = importlib.metadata.version("kb-python")
    kallisto_version = subprocess.check_output([kallisto, "version"], text=True).strip()
    bustools_version = subprocess.check_output([bustools, "version"], text=True).strip()
    if (
        kb_python_version != "0.30.2"
        or kallisto_version != "kallisto, version 0.52.0"
        or bustools_version != "bustools, version 0.45.1"
    ):
        raise RuntimeError(
            "Unexpected native counter versions: "
            f"kb-python={kb_python_version}; {kallisto_version}; {bustools_version}"
        )
    bc_aliases = sorted(samples_by_bc)
    bc_writer = AliasBusWriter(bus_chunks / "gex_aliases", bc_aliases, bc_sequences)
    cr_writer = AliasBusWriter(bus_chunks / "guide_aliases", bc_aliases, bc_sequences)
    gex_counter = CanonicalEC(
        reference["probe_ids"],
        "gex",
        bus_chunks / "gex_reference",
        bc_writer,
        reference["probe_to_gene"],
    )
    guide_counter = CanonicalEC(
        reference["guide_ids"],
        "guide",
        bus_chunks / "guide_reference",
        cr_writer,
    )
    gex_metrics, guide_metrics = [], []

    for index, source in enumerate(pool_gex):
        if not re.fullmatch(r"[0-9a-fA-F]{32}", source["r1_md5"]) or not re.fullmatch(
            r"[0-9a-fA-F]{32}", source["r2_md5"]
        ):
            raise ValueError(f"Invalid GEX checksum pin for {source['source_stem']}")
        if not source["r1_url"].startswith("https://") or not source[
            "r2_url"
        ].startswith("https://"):
            raise ValueError(
                f"GEX pair must use pinned HTTPS URLs: {source['source_stem']}"
            )
        r1 = fastqs / f"gex_{index:03d}_R1.fastq.gz"
        r2 = fastqs / f"gex_{index:03d}_R2.fastq.gz"
        try:
            source_metrics = {
                "source_stem": source["source_stem"],
                "run_accession": source["run_accession"],
                "r1_download": download_ranges(
                    source["r1_url"], int(source["r1_bytes"]), source["r1_md5"], r1
                ),
                "r2_download": download_ranges(
                    source["r2_url"], int(source["r2_bytes"]), source["r2_md5"], r2
                ),
            }
            mapped = run_bus_source(
                f"{index:03d}",
                r1,
                r2,
                "gex",
                args.gex_index,
                BC_COMPONENTS,
                gex_counter,
                bc_onlist,
                bc_replace,
                bustools,
                work,
                logs,
                args.threads,
                kallisto,
            )
            source_metrics.update(mapped)
            gex_metrics.append(source_metrics)
        finally:
            r1.unlink(missing_ok=True)
            r2.unlink(missing_ok=True)
    ec_path_gex, tx_path_gex = gex_counter.write_reference()

    local_archives = downloads / "sra"
    local_archives.mkdir(exist_ok=True)
    source_cache = (
        Path(args.source_cache_dir).resolve() if args.source_cache_dir else None
    )
    sra_bin = Path(args.sra_bin).resolve(strict=True)
    if not (sra_bin / "vdb-dump").is_file() or not (sra_bin / "fastq-dump").is_file():
        raise FileNotFoundError(
            "SRA toolkit directory must contain vdb-dump and fastq-dump"
        )
    for index, source in enumerate(pool_guides):
        accession = source["run_accession"]
        expected_bytes, expected_md5 = (
            int(source["sra_bytes"]),
            source["sra_md5"].lower(),
        )
        if not re.fullmatch(r"[0-9a-f]{32}", expected_md5):
            raise ValueError(f"Invalid SRA checksum pin for {accession}")
        candidates = (
            [
                source_cache / accession / f"{accession}.sra",
                source_cache / f"{accession}.sra",
            ]
            if source_cache
            else []
        )
        cached = next((path for path in candidates if path.is_file()), None)
        if cached:
            verify_cached_file(cached, expected_bytes, expected_md5)
            archive = cached
            download_metrics = {
                "source": "verified_source_cache",
                "archive_bytes": expected_bytes,
                "archive_md5": expected_md5,
            }
        else:
            download_metrics = download_sra(
                accession,
                local_archives,
                logs,
                run_logged,
                expected_bytes=expected_bytes,
                expected_md5=expected_md5,
            )
            archive = local_archives / accession / f"{accession}.sra"
            download_metrics["source"] = "ncbi_sdl"
        r1 = fastqs / f"guide_{index:03d}_R1.fastq.gz"
        r2 = fastqs / f"guide_{index:03d}_R2.fastq.gz"
        pairing_path = fastqs / f"guide_{index:03d}_pairs.json"
        try:
            pair_log = logs / f"guide_{index:03d}_pair.log"
            run_logged(
                [
                    sys.executable,
                    Path(__file__).with_name("stream_sra_pairs.py"),
                    "--sra",
                    archive,
                    "--accession",
                    accession,
                    "--sra-bin",
                    sra_bin,
                    "--r1",
                    r1,
                    "--r2",
                    r2,
                    "--metrics",
                    pairing_path,
                ],
                pair_log,
            )
            pairing = json.loads(pairing_path.read_text())
            if (
                pairing["is_partial"]
                or pairing["paired_records"] * 2 != pairing["sra_sequence_rows"]
            ):
                raise RuntimeError(
                    f"Guide SRA reconstruction was incomplete: {accession}"
                )
            mapped = run_bus_source(
                f"{index:03d}",
                r1,
                r2,
                "guide",
                args.guide_index,
                GUIDE_COMPONENTS,
                guide_counter,
                cr_onlist,
                cr_to_bc_replace,
                bustools,
                work,
                logs,
                args.threads,
                kallisto,
            )
            if mapped["kallisto"]["n_processed"] != pairing["paired_records"]:
                raise ValueError(
                    f"Kallisto did not process every reconstructed pair: {accession}"
                )
            guide_metrics.append(
                {
                    "run_accession": accession,
                    "download": download_metrics,
                    "pairing": pairing,
                    **mapped,
                }
            )
        finally:
            r1.unlink(missing_ok=True)
            r2.unlink(missing_ok=True)
            pairing_path.unlink(missing_ok=True)
            if cached is None:
                shutil.rmtree(local_archives / accession, ignore_errors=True)
    if pool_guides:
        ec_path_guide, tx_path_guide = guide_counter.write_reference()
    else:
        cr_writer.close()
        ec_path_guide = tx_path_guide = None

    reference_sha256 = {
        key: file_hash(path)
        for key, path in {
            "gex_index": args.gex_index,
            "gex_t2g": args.gex_t2g,
            "gex_probe_to_gene": args.gex_probe_to_gene,
            "guide_index": args.guide_index,
            "guide_t2g": args.guide_t2g,
            "cbc_whitelist": args.cbc_whitelist,
            "bc_barcode_variants": args.bc_barcode_variants,
            "gene_probes": args.gene_probes,
            "gene_count_features": args.gene_count_features,
            "guide_features": args.guide_features,
            "guide_count_features": args.guide_count_features,
            "guide_targets": args.guide_targets,
            "probe_barcodes": args.barcode_aliases,
        }.items()
    }
    metrics = {
        "pool_id": args.pool_id,
        "lane_id": args.lane_id,
        "counter_version": {
            "kb_python": kb_python_version,
            "kallisto": kallisto_version,
            "bustools": bustools_version,
        },
        "counter_version_label": (
            f"kb-python {kb_python_version}; {kallisto_version}; {bustools_version}"
        ),
        "geometry": {
            "gex_technology": BC_COMPONENTS,
            "guide_technology": GUIDE_COMPONENTS,
            "gex_cell_barcode": "CBC16+raw GEX BC8; BC suffix corrected to canonical GEX BC8",
            "guide_cell_barcode": "CBC16+CR8; CR suffix corrected then replaced by paired canonical GEX BC8",
            "umi": "R1 bases 16:28 (12bp)",
            "guide_variable_r2": "R2 bases 8:end; no fixed-length truncation",
        },
        "umi_policy": "Bustools --umi-gene on unique probe IDs, then sum probe counts by stable ENSG",
        "mapping_policy": "Kallisto BUS pseudoalignment; multi-target ECs are not allocated with EM or --multimapping",
        "equivalence_class_metrics_scope": (
            "Read-weighted equivalence-class metrics count BUS records after component "
            "barcode correction and BC/CR suffix replacement, before allowlist filtering."
        ),
        "guide_coverage_status": guide_status,
        "gex_sources": gex_metrics,
        "guide_sources": guide_metrics,
        "gex_mapping": gex_counter.metrics,
        "guide_mapping": guide_counter.metrics,
        "barcode_counts": {},
        "reference_sha256": reference_sha256,
        "source_manifest_sha256": {
            "gex": file_hash(args.gex_sources),
            "guide": file_hash(args.guide_sources),
            "guide_coverage": file_hash(args.guide_coverage),
            "samples": file_hash(args.samples),
        },
        "helper_sha256": file_hash(__file__),
    }
    for mapping, runs in (
        (metrics["gex_mapping"], gex_metrics),
        (metrics["guide_mapping"], guide_metrics),
    ):
        processed = sum(int(row["kallisto"]["n_processed"] or 0) for row in runs)
        pseudoaligned = sum(
            int(row["kallisto"]["n_pseudoaligned"] or 0) for row in runs
        )
        mapping.update(
            pseudoaligned_reads=pseudoaligned,
            unmapped_reads=processed - pseudoaligned,
            processed_reads=processed,
            pseudoalignment_rate=(pseudoaligned / processed if processed else 0.0),
        )
    metrics["gex_mapping"]["canonical_equivalence_classes"] = len(
        gex_counter.targets_by_ec_id
    )
    metrics["guide_mapping"]["canonical_equivalence_classes"] = len(
        guide_counter.targets_by_ec_id
    )

    for alias in bc_aliases:
        sample = samples_by_bc[alias]
        cr_alias = sample["cr_alias"]
        bc_bus = bc_writer.paths[alias]
        guide_bus = cr_writer.paths[alias]
        gex_result = None
        filter_metrics = {"input_bus_records": 0, "retained_cells": 0}
        if bc_bus.is_file():
            gex_prefix, filter_metrics = filter_and_count_alias(
                bc_bus,
                ec_path_gex,
                tx_path_gex,
                args.gex_t2g,
                "gex",
                alias,
                work,
                logs,
                bustools,
                args.threads,
            )
            if gex_prefix:
                gex_barcodes, probe_matrix = read_count_matrix(
                    gex_prefix, reference["probe_ids"]
                )
                if any(
                    not re.fullmatch(r"[ACGT]{24}", barcode)
                    or barcode[16:] != bc_sequences[alias]
                    for barcode in gex_barcodes
                ):
                    raise ValueError(
                        f"Filtered GEX barcodes do not match alias {alias}"
                    )
                if len(gex_barcodes) != len(set(gex_barcodes)):
                    raise ValueError(f"Duplicate GEX composite barcode in {alias}")
                gex_matrix = aggregate_probe_matrix(
                    probe_matrix,
                    reference["probe_ids"],
                    reference["probe_to_gene"],
                    reference["gene_ids"],
                )
                filter_metrics["retained_cells"] = len(gex_barcodes)
                gex_result = (gex_barcodes, gex_matrix)
        guide_result = None
        guide_filter_metrics = {"input_bus_records": 0}
        if guide_bus.is_file() and ec_path_guide is not None:
            guide_prefix, guide_filter_metrics = filter_and_count_alias(
                guide_bus,
                ec_path_guide,
                tx_path_guide,
                args.guide_t2g,
                "guide",
                alias,
                work,
                logs,
                bustools,
                args.threads,
            )
            if guide_prefix:
                guide_barcodes, guide_matrix = read_count_matrix(
                    guide_prefix, reference["guide_ids"]
                )
                if any(
                    not re.fullmatch(r"[ACGT]{24}", barcode)
                    or barcode[16:] != bc_sequences[alias]
                    for barcode in guide_barcodes
                ):
                    raise ValueError(
                        f"Guide barcode replacement did not map CR to GEX alias {alias}"
                    )
                if len(guide_barcodes) != len(set(guide_barcodes)):
                    raise ValueError(f"Duplicate guide composite barcode in {alias}")
                guide_result = (guide_barcodes, guide_matrix)
        if gex_result is None:
            gex_barcodes = []
            gex_matrix = sparse.csr_matrix(
                (0, len(reference["gene_ids"])), dtype=np.int32
            )
        else:
            gex_barcodes, gex_matrix = gex_result
        output = output_counts / f"{args.lane_id}__{alias}.h5ad"
        guide_matrix_joined, guide_joined = join_guide_counts(
            gex_barcodes,
            guide_result,
            len(reference["guide_ids"]),
        )
        guide_matrix_before_join = (
            guide_result[1]
            if guide_result is not None
            else sparse.csr_matrix((0, len(reference["guide_ids"])), dtype=np.int32)
        )
        metrics["barcode_counts"][alias] = {
            "sample_id": sample["sample_id"],
            "cr_alias": cr_alias,
            "gex_bus_records": gex_counter.alias_writer.records_by_alias[alias],
            "guide_bus_records": cr_writer.records_by_alias[alias],
            "filter": filter_metrics,
            "guide_filter": guide_filter_metrics,
            "cells": len(gex_barcodes),
            "observed_composite_barcodes": filter_metrics.get(
                "observed_composite_barcodes", 0
            ),
            "retained_cells": len(gex_barcodes),
            "expression_umis": int(gex_matrix.sum()),
            "expression_nnz": int(gex_matrix.nnz),
            "genes_detected": int(gex_matrix.getnnz(axis=0).astype(bool).sum()),
            "guide_umis_before_gex_join": int(guide_matrix_before_join.sum()),
            "guide_nnz_before_gex_join": int(guide_matrix_before_join.nnz),
            "guide_barcodes_joined": guide_joined,
            "guide_umis_joined": int(guide_matrix_joined.sum()),
            "guide_nnz_joined": int(guide_matrix_joined.nnz),
            "guide_positive_cells": int(
                np.count_nonzero(guide_matrix_joined.getnnz(axis=1))
            ),
        }
        write_alias_h5ad(
            output,
            args.lane_id,
            alias,
            cr_alias,
            sample["sample_id"],
            guide_status,
            reference["gene_ids"],
            reference["gene_names"],
            reference["guide_labels"],
            gex_barcodes,
            gex_matrix,
            guide_matrix_joined,
            metrics,
        )
        bc_bus.unlink(missing_ok=True)
        guide_bus.unlink(missing_ok=True)

    metrics_path = outdir / "lane_metrics.json"
    metrics_path.write_text(json.dumps(metrics, sort_keys=True, indent=2) + "\n")
    shutil.rmtree(work)
    return metrics_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-id", required=True)
    parser.add_argument("--lane-id", required=True)
    parser.add_argument("--gex-sources", required=True)
    parser.add_argument("--guide-sources", required=True)
    parser.add_argument("--guide-coverage", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--barcode-aliases", required=True)
    parser.add_argument("--gene-probes", required=True)
    parser.add_argument("--gene-count-features", required=True)
    parser.add_argument("--guide-features", required=True)
    parser.add_argument("--guide-count-features", required=True)
    parser.add_argument("--guide-targets", required=True)
    parser.add_argument("--gex-index", required=True)
    parser.add_argument("--gex-t2g", required=True)
    parser.add_argument("--gex-probe-to-gene", required=True)
    parser.add_argument("--guide-index", required=True)
    parser.add_argument("--guide-t2g", required=True)
    parser.add_argument("--cbc-whitelist", required=True)
    parser.add_argument("--bc-barcode-variants", required=True)
    parser.add_argument("--sra-bin", required=True)
    parser.add_argument("--source-cache-dir", default="")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    execute(args)


if __name__ == "__main__":
    main()
