#!/usr/bin/env python3
"""Count one physical 10x Flex lane from verified GEX FASTQ and guide SRA."""

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as scipy_io
from scipy import sparse

from bam_to_fastq import download_ranges
from stream_count import download_sra


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def md5_file(path):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_logged(command, log_path, timeout=12 * 3600):
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


def download_verified(url, expected_bytes, expected_md5, output):
    if not url.startswith("https://") or not re.fullmatch(
        r"[0-9a-f]{32}", expected_md5
    ):
        raise ValueError("GEX source needs an HTTPS URL and a pinned MD5")
    expected_bytes = int(expected_bytes)
    output = Path(output)
    return download_ranges(url, expected_bytes, expected_md5, output)


def verify_file(path, expected_bytes, expected_md5):
    path = Path(path)
    if (
        path.stat().st_size != int(expected_bytes)
        or md5_file(path) != expected_md5.lower()
    ):
        raise RuntimeError(f"Cached SRA does not match its pinned size/MD5: {path}")


def accumulate_ibu(cyto, source, target, logdir, threads):
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        shutil.move(str(source), target)
        return
    merged = target.with_name(target.name + ".merge")
    run_logged(
        [cyto, "ibu", "cat", "-t", threads, "-o", merged, target, source],
        logdir / (target.stem + ".cat.log"),
    )
    merged.replace(target)
    Path(source).unlink()


def load_mtx(directory, expected_suffix, expected_features):
    directory = Path(directory)
    matrix_path = directory / "matrix.mtx.gz"

    def output_file(stem):
        candidates = [directory / f"{stem}.txt.gz", directory / f"{stem}.tsv.gz"]
        found = [path for path in candidates if path.is_file()]
        if len(found) != 1:
            raise ValueError(f"Expected exactly one Cyto {stem} output in {directory}")
        return found[0]

    barcode_path = output_file("barcodes")
    feature_path = output_file("features")
    with gzip.open(feature_path, "rt", encoding="ascii") as handle:
        features = [line.rstrip("\r\n").split("\t", 1)[0] for line in handle]
    if features != expected_features:
        raise ValueError(
            "Cyto count feature order differs from its pinned feature table"
        )
    with gzip.open(barcode_path, "rt", encoding="ascii") as handle:
        barcodes = [line.rstrip("\r\n") for line in handle]
    suffix = "-" + expected_suffix
    if any(not barcode.endswith(suffix) for barcode in barcodes):
        raise ValueError("Cyto output barcode suffix differs from the lane/alias")
    cells = [barcode[: -len(suffix)] for barcode in barcodes]
    if len(cells) != len(set(cells)) or any(len(cell) != 16 for cell in cells):
        raise ValueError("Cyto output has duplicate or malformed 16-base cell barcodes")
    with gzip.open(matrix_path, "rb") as handle:
        feature_by_cell = scipy_io.mmread(handle).tocsr()
    if feature_by_cell.shape != (len(features), len(cells)):
        raise ValueError(
            "Matrix dimensions do not match the Cyto barcode/feature files"
        )
    matrix = feature_by_cell.T.tocsr()
    if matrix.nnz and (
        not np.all(np.isfinite(matrix.data))
        or np.any(matrix.data < 0)
        or not np.all(matrix.data == np.floor(matrix.data))
        or np.any(matrix.data > np.iinfo(np.int32).max)
    ):
        raise ValueError("Cyto count matrix contains invalid or non-integer UMI counts")
    matrix.data = matrix.data.astype(np.int32, copy=False)
    return cells, matrix


def guide_labels(target_rows):
    labels = {}
    for row in target_rows:
        guide_id = row["guide_id"]
        symbol = row["target_symbol"].strip()
        ensg = row["target_ensg"].strip()
        labels[guide_id] = (
            f"{symbol}_{guide_id}__{ensg}" if ensg else f"{symbol}_{guide_id}"
        )
    return labels


def matrix_features(path):
    return [
        line.split("\t", 1)[0] for line in Path(path).read_text().splitlines() if line
    ]


def create_empty_h5ad(
    output,
    sample,
    lane_id,
    bc_alias,
    cr_alias,
    guide_status,
    gene_ids,
    gene_names,
    guide_names,
    cyto_version,
):
    obs = pd.DataFrame(index=pd.Index([], dtype="object"))
    for column in (
        "sample_id",
        "lane_id",
        "bc_alias",
        "cr_alias",
        "guide_coverage_status",
    ):
        obs[column] = pd.Series(index=obs.index, dtype=object)
    var = pd.DataFrame(
        {"gene_name": gene_names}, index=pd.Index(gene_ids, name="gene_id")
    )
    adata = ad.AnnData(
        X=sparse.csr_matrix((0, len(gene_ids)), dtype=np.int32), obs=obs, var=var
    )
    adata.obsm["guides"] = sparse.csr_matrix((0, len(guide_names)), dtype=np.int32)
    adata.uns["guide_names"] = guide_names
    adata.uns["flex_cyto_version"] = cyto_version
    adata.uns["flex_sample_id"] = sample["sample_id"]
    adata.uns["flex_lane_id"] = lane_id
    adata.uns["flex_bc_alias"] = bc_alias
    adata.uns["flex_cr_alias"] = cr_alias
    adata.uns["flex_guide_coverage_status"] = guide_status
    adata.write_h5ad(output, compression="gzip")


def write_alias_h5ad(
    output,
    bc_alias,
    cr_alias,
    lane_id,
    sample,
    alias_sequences,
    gex_cells,
    gex_matrix,
    guide_data,
    gene_ids,
    gene_names,
    guide_ids,
    guide_labels_by_id,
    guide_status,
    cyto_version,
):
    obs_names = [f"{cell}{alias_sequences[bc_alias]}-1_{lane_id}" for cell in gex_cells]
    if len(obs_names) != len(set(obs_names)):
        raise ValueError(f"Duplicate author-style cell IDs for {lane_id} {bc_alias}")
    obs = pd.DataFrame(
        {
            "sample_id": sample["sample_id"],
            "lane_id": lane_id,
            "bc_alias": bc_alias,
            "cr_alias": cr_alias,
            "guide_coverage_status": guide_status,
        },
        index=pd.Index(obs_names),
    )
    var = pd.DataFrame(
        {"gene_name": gene_names}, index=pd.Index(gene_ids, name="gene_id")
    )
    guide_matrix = sparse.csr_matrix((len(gex_cells), len(guide_ids)), dtype=np.int32)
    matched_guide_barcodes = 0
    if guide_data is not None and gex_cells:
        guide_cells, guide_matrix_by_cell = guide_data
        row_by_cell = {cell: row for row, cell in enumerate(gex_cells)}
        source_rows, target_rows = [], []
        for guide_row, cell in enumerate(guide_cells):
            target_row = row_by_cell.get(cell)
            if target_row is not None:
                source_rows.append(guide_row)
                target_rows.append(target_row)
        if source_rows:
            subset = guide_matrix_by_cell[source_rows].tocoo()
            mapped_rows = np.asarray(target_rows, dtype=np.int32)[subset.row]
            guide_matrix = sparse.coo_matrix(
                (subset.data, (mapped_rows, subset.col)),
                shape=(len(gex_cells), len(guide_ids)),
                dtype=np.int32,
            ).tocsr()
        matched_guide_barcodes = len(source_rows)
    adata = ad.AnnData(X=gex_matrix, obs=obs, var=var)
    adata.obsm["guides"] = guide_matrix
    adata.uns["guide_names"] = [guide_labels_by_id[guide_id] for guide_id in guide_ids]
    adata.uns["flex_cyto_version"] = cyto_version
    adata.uns["flex_sample_id"] = sample["sample_id"]
    adata.uns["flex_lane_id"] = lane_id
    adata.uns["flex_bc_alias"] = bc_alias
    adata.uns["flex_cr_alias"] = cr_alias
    adata.uns["flex_guide_coverage_status"] = guide_status
    adata.write_h5ad(output, compression="gzip")
    return {"cells": len(gex_cells), "guide_barcodes_joined": matched_guide_barcodes}


def execute(args):
    identifier = r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}"
    if not re.fullmatch(identifier, args.pool_id) or not re.fullmatch(
        identifier, args.lane_id
    ):
        raise ValueError("Pool and lane IDs must be safe non-empty identifiers")
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=False)
    logs = outdir / "logs"
    work = outdir / "work"
    downloads = work / "downloads"
    fastqs = work / "fastqs"
    gex_maps = work / "gex_maps"
    guide_maps = work / "guide_maps"
    accum_gex = work / "ibu_gex"
    accum_guide = work / "ibu_guide"
    for path in (logs, downloads, fastqs, gex_maps, guide_maps, accum_gex, accum_guide):
        path.mkdir(parents=True)

    gex_rows = [
        row
        for row in read_tsv(args.gex_sources)
        if row["pool_id"] == args.pool_id and row["lane_id"] == args.lane_id
    ]
    guide_rows = [
        row
        for row in read_tsv(args.guide_sources)
        if row["pool_id"] == args.pool_id and row["lane_id"] == args.lane_id
    ]
    coverage_rows = [
        row
        for row in read_tsv(args.guide_coverage)
        if row["pool_id"] == args.pool_id and row["lane_id"] == args.lane_id
    ]
    if not gex_rows:
        raise ValueError(f"No GEX sources pinned for {args.pool_id} {args.lane_id}")
    if len(coverage_rows) != 1:
        raise ValueError(f"Expected one guide coverage row for {args.lane_id}")
    coverage = coverage_rows[0]
    if int(coverage["gex_source_pairs"]) != len(gex_rows) or int(
        coverage["guide_sra_runs"]
    ) != len(guide_rows):
        raise ValueError("Source manifest and guide coverage table disagree")
    guide_status = coverage["guide_status"]
    if guide_status not in {"available", "no_archived_guide_sra"} or (
        (guide_status == "available") != bool(guide_rows)
    ):
        raise ValueError("Guide-coverage status does not match pinned SRA rows")

    samples_by_bc, samples_by_cr = {}, {}
    for sample in read_tsv(args.samples):
        if sample["pool_id"] != args.pool_id:
            continue
        for field in ("sample_id", "pool_id"):
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", sample[field]):
                raise ValueError(
                    f"Unsafe {field} in sample metadata: {sample[field]!r}"
                )
        for pair in sample["bc_cr_pairs"].split(";"):
            bc_alias, cr_alias = pair.split(":")
            if not re.fullmatch(
                r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", bc_alias
            ) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", cr_alias):
                raise ValueError(f"Unsafe BC/CR alias in sample metadata: {pair!r}")
            if args.lane_id not in sample["pool_lanes"].split(";"):
                raise ValueError(
                    f"Sample {sample['sample_id']} does not include lane {args.lane_id}"
                )
            if bc_alias in samples_by_bc or cr_alias in samples_by_cr:
                raise ValueError("Duplicate BC/CR alias in sample metadata")
            samples_by_bc[bc_alias] = (cr_alias, sample)
            samples_by_cr[cr_alias] = (bc_alias, sample)
    aliases = read_tsv(args.barcode_aliases)
    bc_sequences = {row["bc_alias"]: row["bc_sequence"] for row in aliases}
    cr_sequences = {row["cr_alias"]: row["cr_sequence"] for row in aliases}
    cr_aliases = {row["cr_alias"] for row in aliases}
    if (
        len(bc_sequences) != len(aliases)
        or len(cr_sequences) != len(aliases)
        or len(cr_aliases) != len(aliases)
        or set(samples_by_bc) != set(bc_sequences)
        or set(samples_by_cr) != cr_aliases
        or set(samples_by_cr) != set(cr_sequences)
    ):
        raise ValueError(
            "Sample BC/CR mapping does not cover every barcode-reference alias"
        )
    for alias, sequence in (*bc_sequences.items(), *cr_sequences.items()):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", alias):
            raise ValueError(f"Unsafe barcode-reference alias: {alias!r}")
        if not re.fullmatch(r"[ACGT]{8}", sequence):
            raise ValueError(f"Invalid 8-base probe barcode sequence for {alias}")

    probes = []
    with Path(args.gene_probes).open() as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise ValueError(
                    "Gene probe reference must have three tab-separated fields"
                )
            probes.append({"gene_id": fields[0], "gene_name": fields[1]})
    gene_names_by_id = {}
    for row in probes:
        gene_id, gene_name = row["gene_id"], row["gene_name"]
        if gene_id in gene_names_by_id and gene_names_by_id[gene_id] != gene_name:
            raise ValueError(f"Conflicting gene names for {gene_id}")
        gene_names_by_id[gene_id] = gene_name
    gene_ids = matrix_features(args.gene_count_features)
    if len(gene_ids) != len(probes) or set(gene_ids) != set(gene_names_by_id):
        raise ValueError("Gene count feature table does not match the probe reference")
    unique_gene_ids = list(dict.fromkeys(gene_ids))
    gene_names = [gene_names_by_id[gene_id] for gene_id in unique_gene_ids]

    guide_targets = read_tsv(args.guide_targets)
    label_by_guide = guide_labels(guide_targets)
    guide_ids = matrix_features(args.guide_count_features)
    guide_feature_ids = matrix_features(args.guide_features)
    if len(guide_ids) != len(guide_feature_ids) or set(guide_ids) != set(
        label_by_guide
    ):
        raise ValueError(
            "Guide count feature table does not match guide target metadata"
        )
    if len(set(guide_ids)) != len(guide_ids):
        raise ValueError("Guide feature IDs are not unique")

    cyto = str(Path(args.cyto).resolve(strict=True))
    resources = Path(args.cyto_resources).resolve(strict=True)
    whitelist = resources / "737K-fixed-rna-profiling.txt.gz"
    rna_probes = resources / "probe-barcodes-fixed-rna-profiling-rna.txt"
    crispr_probes = resources / "probe-barcodes-fixed-rna-profiling-crispr.txt"
    if not all(path.is_file() for path in (whitelist, rna_probes, crispr_probes)):
        raise FileNotFoundError("Cyto reference resources are incomplete")
    regex = "^(?:" + "|".join(re.escape(alias) for alias in sorted(bc_sequences)) + ")$"
    crispr_regex = (
        "^(?:" + "|".join(re.escape(alias) for alias in sorted(cr_aliases)) + ")$"
    )
    nthreads = max(1, args.threads)
    metrics = {
        "pool_id": args.pool_id,
        "lane_id": args.lane_id,
        "cyto_version": subprocess.check_output([cyto, "--version"], text=True).strip(),
        "geometry": {"gex": "gex-v1", "guide": "crispr-v1"},
        "guide_coverage_status": guide_status,
        "gex_sources": [],
        "guide_sources": [],
        "barcode_counts": {},
    }

    # Each original pair is verified, mapped, then removed before its next pair starts.
    for index, source in enumerate(gex_rows):
        stem = source["source_stem"]
        r1 = fastqs / f"gex_{index:03d}_R1.fastq.gz"
        r2 = fastqs / f"gex_{index:03d}_R2.fastq.gz"
        pair = {"source_stem": stem, "run_accession": source["run_accession"]}
        pair["r1"] = download_verified(
            source["r1_url"], source["r1_bytes"], source["r1_md5"], r1
        )
        pair["r2"] = download_verified(
            source["r2_url"], source["r2_bytes"], source["r2_md5"], r2
        )
        map_dir = gex_maps / f"map_{index:03d}"
        run_logged(
            [
                cyto,
                "map",
                "gex",
                "--preset",
                "gex-v1",
                "--whitelist",
                whitelist,
                "--gex",
                args.gene_probes,
                "--probes",
                rna_probes,
                "--probe-regex",
                regex,
                "--num-threads",
                nthreads,
                "--min-ibu-records",
                0,
                "--outdir",
                map_dir,
                r1,
                r2,
            ],
            logs / f"gex_map_{index:03d}.log",
            timeout=24 * 3600,
        )
        map_metrics = map_dir / "stats" / "mapping_run.json"
        if map_metrics.is_file():
            pair["mapping"] = json.loads(map_metrics.read_text())
        for bc_alias in sorted(samples_by_bc):
            ibu = map_dir / "ibu" / f"{bc_alias}.ibu"
            if ibu.is_file():
                accumulate_ibu(cyto, ibu, accum_gex / f"{bc_alias}.ibu", logs, nthreads)
        shutil.rmtree(map_dir)
        r1.unlink()
        r2.unlink()
        metrics["gex_sources"].append(pair)

    # Guide SRA files are downloaded one at a time; the helper validates ordered read-block mates.
    sra_bin = Path(args.sra_bin).resolve(strict=True)
    local_archives = downloads / "sra"
    local_archives.mkdir()
    source_cache = (
        Path(args.source_cache_dir).resolve() if args.source_cache_dir else None
    )
    for index, source in enumerate(guide_rows):
        accession = source["run_accession"]
        expected_bytes, expected_md5 = (
            int(source["sra_bytes"]),
            source["sra_md5"].lower(),
        )
        candidates = []
        if source_cache:
            candidates = [
                source_cache / accession / f"{accession}.sra",
                source_cache / f"{accession}.sra",
            ]
        cached = next((path for path in candidates if path.is_file()), None)
        if cached:
            verify_file(cached, expected_bytes, expected_md5)
            archive = cached
            download_metrics = {
                "archive_bytes": expected_bytes,
                "archive_md5": expected_md5,
                "source": "verified_source_cache",
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
        pair_metrics = fastqs / f"guide_{index:03d}_pairs.json"
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
                pair_metrics,
            ],
            logs / f"guide_pair_{index:03d}.log",
            timeout=24 * 3600,
        )
        pairing = json.loads(pair_metrics.read_text())
        if (
            pairing["is_partial"]
            or pairing["paired_records"] * 2 != pairing["sra_sequence_rows"]
        ):
            raise RuntimeError(f"Guide SRA extraction was incomplete: {accession}")
        map_dir = guide_maps / f"map_{index:03d}"
        run_logged(
            [
                cyto,
                "map",
                "crispr",
                "--preset",
                "crispr-v1",
                "--whitelist",
                whitelist,
                "--guides",
                args.guide_features,
                "--probes",
                crispr_probes,
                "--probe-regex",
                crispr_regex,
                "--num-threads",
                nthreads,
                "--min-ibu-records",
                0,
                "--outdir",
                map_dir,
                r1,
                r2,
            ],
            logs / f"guide_map_{index:03d}.log",
            timeout=24 * 3600,
        )
        map_metrics = map_dir / "stats" / "mapping_run.json"
        run = {
            "run_accession": accession,
            "download": download_metrics,
            "pairing": pairing,
        }
        if map_metrics.is_file():
            run["mapping"] = json.loads(map_metrics.read_text())
        for cr_alias in sorted(samples_by_cr):
            ibu = map_dir / "ibu" / f"{cr_alias}.ibu"
            if ibu.is_file():
                accumulate_ibu(
                    cyto, ibu, accum_guide / f"{cr_alias}.ibu", logs, nthreads
                )
        shutil.rmtree(map_dir)
        r1.unlink()
        r2.unlink()
        pair_metrics.unlink()
        if cached is None:
            shutil.rmtree(local_archives / accession)
        metrics["guide_sources"].append(run)

    output_counts = outdir / "lane_counts"
    output_counts.mkdir()
    final_gene_ids = unique_gene_ids
    gex_id_order = final_gene_ids
    # Count all BC/CR aliases once per lane. Missing guide lanes get an explicit status in every H5AD.
    for bc_alias in sorted(samples_by_bc):
        cr_alias, sample = samples_by_bc[bc_alias]
        gex_ibu = accum_gex / f"{bc_alias}.ibu"
        guide_ibu = accum_guide / f"{cr_alias}.ibu"
        gex_data = None
        guide_data = None
        if gex_ibu.exists():
            sorted_ibu = gex_ibu.with_suffix(".sort.ibu")
            corrected_ibu = gex_ibu.with_suffix(".umi.ibu")
            run_logged(
                [cyto, "ibu", "sort", "-i", gex_ibu, "-o", sorted_ibu, "-T", nthreads],
                logs / f"{bc_alias}.sort.log",
                timeout=24 * 3600,
            )
            run_logged(
                [
                    cyto,
                    "ibu",
                    "umi",
                    "-i",
                    sorted_ibu,
                    "-o",
                    corrected_ibu,
                    "-T",
                    nthreads,
                    "-l",
                    logs / f"{bc_alias}.umi.json",
                ],
                logs / f"{bc_alias}.umi.log",
                timeout=24 * 3600,
            )
            sorted_ibu.unlink()
            gex_ibu.unlink()
            count_dir = work / f"count_gex_{bc_alias}"
            run_logged(
                [
                    cyto,
                    "ibu",
                    "count",
                    "-i",
                    corrected_ibu,
                    "--mtx",
                    "-f",
                    args.gene_count_features,
                    "-C",
                    1,
                    "-s",
                    f"{bc_alias}_{args.lane_id}",
                    "-o",
                    count_dir,
                    "-t",
                    nthreads,
                ],
                logs / f"{bc_alias}.count.log",
                timeout=24 * 3600,
            )
            gex_data = load_mtx(count_dir, f"{bc_alias}_{args.lane_id}", gex_id_order)
            shutil.rmtree(count_dir)
            corrected_ibu.unlink()
        if guide_ibu.exists():
            sorted_ibu = guide_ibu.with_suffix(".sort.ibu")
            corrected_ibu = guide_ibu.with_suffix(".umi.ibu")
            run_logged(
                [
                    cyto,
                    "ibu",
                    "sort",
                    "-i",
                    guide_ibu,
                    "-o",
                    sorted_ibu,
                    "-T",
                    nthreads,
                ],
                logs / f"{cr_alias}.sort.log",
                timeout=24 * 3600,
            )
            run_logged(
                [
                    cyto,
                    "ibu",
                    "umi",
                    "-i",
                    sorted_ibu,
                    "-o",
                    corrected_ibu,
                    "-T",
                    nthreads,
                    "-l",
                    logs / f"{cr_alias}.umi.json",
                ],
                logs / f"{cr_alias}.umi.log",
                timeout=24 * 3600,
            )
            sorted_ibu.unlink()
            guide_ibu.unlink()
            count_dir = work / f"count_guide_{cr_alias}"
            run_logged(
                [
                    cyto,
                    "ibu",
                    "count",
                    "-i",
                    corrected_ibu,
                    "--mtx",
                    "-f",
                    args.guide_count_features,
                    "-C",
                    1,
                    "-s",
                    f"{bc_alias}_{args.lane_id}",
                    "-o",
                    count_dir,
                    "-t",
                    nthreads,
                ],
                logs / f"{cr_alias}.count.log",
                timeout=24 * 3600,
            )
            guide_data = load_mtx(count_dir, f"{bc_alias}_{args.lane_id}", guide_ids)
            shutil.rmtree(count_dir)
            corrected_ibu.unlink()

        if gex_data is None:
            gex_cells = []
            expression_matrix = sparse.csr_matrix(
                (0, len(gex_id_order)), dtype=np.int32
            )
        else:
            gex_cells, expression_matrix = gex_data
        path = output_counts / f"{args.lane_id}__{bc_alias}.h5ad"
        if not gex_cells:
            create_empty_h5ad(
                path,
                sample,
                args.lane_id,
                bc_alias,
                cr_alias,
                guide_status,
                gex_id_order,
                gene_names,
                [label_by_guide[guide_id] for guide_id in guide_ids],
                metrics["cyto_version"],
            )
            metrics["barcode_counts"][bc_alias] = {
                "cr_alias": cr_alias,
                "cells": 0,
                "expression_umis": 0,
                "expression_nnz": 0,
                "guide_umis": 0,
                "guide_nnz": 0,
                "guide_barcodes_joined": 0,
            }
        else:
            summary = write_alias_h5ad(
                path,
                bc_alias,
                cr_alias,
                args.lane_id,
                sample,
                bc_sequences,
                gex_cells,
                expression_matrix,
                guide_data,
                gex_id_order,
                gene_names,
                guide_ids,
                label_by_guide,
                guide_status,
                metrics["cyto_version"],
            )
            metrics["barcode_counts"][bc_alias] = {
                "cr_alias": cr_alias,
                "cells": summary["cells"],
                "expression_umis": int(expression_matrix.sum()),
                "expression_nnz": int(expression_matrix.nnz),
                "guide_umis": int(guide_data[1].sum()) if guide_data else 0,
                "guide_nnz": int(guide_data[1].nnz) if guide_data else 0,
                "guide_barcodes_joined": summary["guide_barcodes_joined"],
            }

    metrics_path = outdir / "lane_metrics.json"
    metrics_path.write_text(json.dumps(metrics, sort_keys=True, indent=2) + "\n")
    # Keep only reusable lane products and provenance; raw sources, FASTQs and IBU scratch are transient.
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
    parser.add_argument("--cyto", default="/usr/local/bin/cyto")
    parser.add_argument("--cyto-resources", default="/opt/cyto/resources")
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
