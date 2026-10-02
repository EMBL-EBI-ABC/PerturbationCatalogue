#!/usr/bin/env python3
"""Validate pinned Zhu sources and regenerate deterministic Flex references."""

import argparse
import csv
import hashlib
import io
import re
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parent
POOL_RANGES = {
    "CD4i_R1_L01-23": ("R1", 1, 23),
    "CD4i_R2_L01-24": ("R2", 1, 24),
    "CD4i_R2_L25-48": ("R2", 25, 48),
}
AUTHOR_REPOSITORY_COMMIT = "aa5c84a973c0e1a090b0072dc5b080bf7fbbed38"
AUTHOR_SAMPLE_SHA256 = (
    "766134d11dabb5d63388e00b4d687a809c0254eb3c19d77f439a6bf9dad7abd4"
)
AUTHOR_GUIDE_LIBRARY_SHA256 = (
    "00a1bec2afc2082fc79765531696d7e22672a8ba904ea54c035858f425a657a8"
)
GUIDE_TARGETS_SHA256 = (
    "fa5fd9c8c7aae7ff2860c88ce961f0b0a36a5401db14a8dc289fb3c2f7c1243d"
)
OUTPUTS = {
    "gene_probes.tsv",
    "gene_count_features.tsv",
    "guide_features.tsv",
    "guide_count_features.tsv",
    "samples.tsv",
    "guide_coverage.tsv",
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def lane_id(pool_id, lane):
    pool, _, _ = POOL_RANGES[pool_id]
    return f"CD4i_{pool}L{lane:02d}"


def author_samples():
    sample_path = ROOT / "author_sample_metadata.suppl_table.csv"
    if hashlib.sha256(sample_path.read_bytes()).hexdigest() != AUTHOR_SAMPLE_SHA256:
        raise ValueError(
            "Author sample metadata does not match the pinned repository revision"
        )
    rows = list(csv.DictReader(sample_path.open(newline="")))
    if len(rows) != 12:
        raise ValueError(f"Expected 12 author sample rows, found {len(rows)}")
    result = []
    condition_map = {"Rest": "rest", "Stim8hr": "stim8hr", "Stim48hr": "stim48hr"}
    pattern = re.compile(r"^CD4i_(R1|R2)_(D[1-4])_(Rest|Stim8hr|Stim48hr)$")
    for row in rows:
        match = pattern.fullmatch(row["cell_sample_id"])
        if not match or row["10xrun_id"] != f"CD4i_{match[1]}":
            raise ValueError(f"Unexpected author sample row: {row['cell_sample_id']}")
        instrument, donor, author_condition = match.groups()
        condition = condition_map[author_condition]
        if (
            row["library_prep_kit"] != "GEMX_flex_v1"
            or row["sequencing_platform"] != "Ultima"
        ):
            raise ValueError(
                f"Unexpected assay or sequencing platform: {row['cell_sample_id']}"
            )
        pool = (
            "CD4i_R1_L01-23"
            if instrument == "R1"
            else "CD4i_R2_L25-48" if condition == "stim48hr" else "CD4i_R2_L01-24"
        )
        barcode_text = row["probe_hyb_loading"]
        bc_match = re.search(r"\bBC(\d{3})-(\d{3})\b", barcode_text)
        cr_match = re.search(r"\bCR(\d{3})-(\d{3})\b", barcode_text)
        if not bc_match or not cr_match:
            raise ValueError(
                f"Missing BC/CR ranges in author row: {row['cell_sample_id']}"
            )
        bc_start, bc_end = map(int, bc_match.groups())
        cr_start, cr_end = map(int, cr_match.groups())
        if bc_end - bc_start != cr_end - cr_start or (bc_start, bc_end) != (
            cr_start,
            cr_end,
        ):
            raise ValueError(
                f"Author BC and CR ranges do not pair: {row['cell_sample_id']}"
            )
        pairs = [f"BC{i:03d}:CR{i:03d}" for i in range(bc_start, bc_end + 1)]
        result.append(
            {
                "sample_id": f"zhu_2025_{donor}_{condition}_cl",
                "pool_id": pool,
                "donor": donor,
                "condition": condition,
                "donor_id": row["donor_id"],
                "author_sample": row["cell_sample_id"],
                "author_library": row["library_id"],
                "library_prep_kit": row["library_prep_kit"],
                "sequencing_platform": row["sequencing_platform"],
                "probe_hyb_loading": barcode_text,
                "bc_cr_pairs": ";".join(pairs),
            }
        )
    if len({row["sample_id"] for row in result}) != 12:
        raise ValueError(
            "Author sample metadata does not produce twelve unique donor/condition outputs"
        )
    return result


def check_lane(pool_id, lane, name):
    if pool_id not in POOL_RANGES:
        raise ValueError(f"Unknown physical pool {pool_id}")
    instrument_pool, first, last = POOL_RANGES[pool_id]
    if not first <= lane <= last:
        raise ValueError(f"{name} lane {lane} is outside {pool_id}")
    if f"CD4i_{instrument_pool}L{lane:02d}" not in name:
        raise ValueError(f"{name} does not identify {pool_id} lane {lane}")


def validate_sources():
    gex = read_tsv(ROOT / "gex_sources.tsv")
    guide = read_tsv(ROOT / "guide_sources.tsv")
    if len(gex) != 463 or len(guide) != 253:
        raise ValueError(
            f"Pinned source count changed: GEX={len(gex)}, guide={len(guide)}"
        )

    seen_stems, seen_urls = set(), set()
    gex_lanes = Counter()
    for row in gex:
        lane = int(row["lane"])
        check_lane(row["pool_id"], lane, row["source_stem"])
        if row["lane_id"] != lane_id(row["pool_id"], lane):
            raise ValueError(f"Incorrect lane_id for {row['source_stem']}")
        if row["source_stem"] in seen_stems:
            raise ValueError(f"Duplicate GEX source stem: {row['source_stem']}")
        seen_stems.add(row["source_stem"])
        for mate in ("r1", "r2"):
            url, size, digest = (
                row[mate + "_url"],
                row[mate + "_bytes"],
                row[mate + "_md5"],
            )
            if urlsplit(url).scheme != "https" or not urlsplit(url).hostname:
                raise ValueError(f"Invalid {mate.upper()} URL for {row['source_stem']}")
            if (
                not size.isdigit()
                or int(size) <= 0
                or not re.fullmatch(r"[0-9a-f]{32}", digest)
            ):
                raise ValueError(
                    f"Invalid {mate.upper()} size/checksum for {row['source_stem']}"
                )
            if url in seen_urls:
                raise ValueError(f"Duplicate pinned source URL: {url}")
            seen_urls.add(url)
        if row["r1_url"].replace("_R1_001", "_R2_001") != row["r2_url"]:
            raise ValueError(
                f"R1/R2 filenames do not form a pair: {row['source_stem']}"
            )
        gex_lanes[(row["pool_id"], lane)] += 1
    expected_gex_lanes = {
        (pool, lane)
        for pool, (_, first, last) in POOL_RANGES.items()
        for lane in range(first, last + 1)
    }
    if set(gex_lanes) != expected_gex_lanes:
        raise ValueError("GEX manifest does not cover all 71 expected lanes")

    seen_runs = set()
    guide_lanes = Counter()
    for row in guide:
        lane = int(row["lane"])
        check_lane(row["pool_id"], lane, row["library_name"])
        if row["lane_id"] != lane_id(row["pool_id"], lane):
            raise ValueError(f"Incorrect guide lane_id for {row['run_accession']}")
        if row["run_accession"] in seen_runs:
            raise ValueError(f"Duplicate guide SRA accession: {row['run_accession']}")
        seen_runs.add(row["run_accession"])
        if not re.fullmatch(r"(?:SRR|ERR|DRR)\d+", row["run_accession"]):
            raise ValueError(f"Invalid SRA accession: {row['run_accession']}")
        if not row["sra_bytes"].isdigit() or int(row["sra_bytes"]) <= 0:
            raise ValueError(f"Invalid SRA size for {row['run_accession']}")
        if not re.fullmatch(r"[0-9a-f]{32}", row["sra_md5"]):
            raise ValueError(f"Invalid SRA MD5 for {row['run_accession']}")
        guide_lanes[(row["pool_id"], lane)] += 1
    return gex, guide, gex_lanes, guide_lanes


def feature_tables():
    probe_path = (
        ROOT / "Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"
    )
    if (
        hashlib.md5(probe_path.read_bytes()).hexdigest()
        != "8d071b87b07a98cc7aabd6dcad526fef"
    ):
        raise ValueError(
            "10x 2024-A probe reference MD5 does not match the pinned source"
        )
    with probe_path.open(newline="") as handle:
        lines = iter(handle)
        metadata = [next(lines).rstrip("\r\n") for _ in range(5)]
        if metadata != [
            "#probe_set_file_format=3.0",
            "#panel_name=Chromium Human Transcriptome Probe Set v1.1.0",
            "#panel_type=predesigned",
            "#reference_genome=GRCh38",
            "#reference_version=2024-A",
        ]:
            raise ValueError(
                "Probe-set metadata/version is not the pinned 2024-A reference"
            )
        rows = csv.DictReader(lines)
        if rows.fieldnames != [
            "gene_id",
            "probe_seq",
            "probe_id",
            "included",
            "region",
            "gene_name",
        ]:
            raise ValueError("Unexpected 10x probe-set columns")
        all_rows = list(rows)
    statuses = {}
    for row in all_rows:
        statuses.setdefault(row["gene_id"], set()).add(row["included"])
    if any(len(value) != 1 for value in statuses.values()):
        raise ValueError("Included status differs among probes for one gene")
    gex_lines, seen_probes = [], set()
    for row in all_rows:
        if statuses[row["gene_id"]] != {"TRUE"}:
            continue
        sequence = row["probe_seq"]
        if len(sequence) != 50 or set(sequence) - set("ACGT"):
            raise ValueError(f"Invalid 10x probe sequence: {row['probe_id']}")
        if sequence in seen_probes:
            raise ValueError(f"Duplicate included probe sequence: {row['probe_id']}")
        seen_probes.add(sequence)
        gex_lines.append("\t".join((row["gene_id"], row["gene_name"], sequence)))

    author_guide_path = ROOT / "author_sgrna_library_metadata.suppl_table.csv"
    if (
        hashlib.sha256(author_guide_path.read_bytes()).hexdigest()
        != AUTHOR_GUIDE_LIBRARY_SHA256
    ):
        raise ValueError(
            "Author sgRNA library does not match the pinned repository revision"
        )
    author_guides = {
        row["sgRNA"]: row["seq"].upper()
        for row in csv.DictReader(author_guide_path.open(newline=""))
    }
    if len(author_guides) != 26_504:
        raise ValueError(
            "Pinned author guide library must contain 26,504 unique guide IDs"
        )
    targets_path = ROOT / "guide_targets.tsv"
    if hashlib.sha256(targets_path.read_bytes()).hexdigest() != GUIDE_TARGETS_SHA256:
        raise ValueError("Guide target metadata does not match its pinned source")
    targets = read_tsv(targets_path)
    sequences = read_tsv(ROOT / "guide_sequences.tsv")
    target_map = {row["guide_id"]: row for row in targets}
    if len(target_map) != len(targets) or len(sequences) != 26_504:
        raise ValueError("Guide target/source rows must contain 26,504 unique IDs")
    complement = str.maketrans("ACGT", "TGCA")
    guide_lines, seen_ids, seen_sequences = [], set(), set()
    for row in sequences:
        name, seq = row["guide_id"], row["author_sequence"].upper()
        if author_guides.get(name) != seq:
            raise ValueError(
                f"Guide sequence differs from the pinned author table: {name}"
            )
        if name in seen_ids or name not in target_map:
            raise ValueError(f"Duplicate or unmatched author guide ID: {name}")
        if len(seq) != 20 or set(seq) - set("ACGT"):
            raise ValueError(f"Invalid author guide sequence: {name}")
        if seq in seen_sequences:
            raise ValueError(f"Duplicate author guide sequence: {name}")
        seen_ids.add(name)
        seen_sequences.add(seq)
        guide_lines.append(
            "\t".join(
                (
                    name,
                    "GCTATGCTGTTTCCAGCTTAGCTCTTAAAC",
                    seq.translate(complement)[::-1],
                )
            )
        )
    if seen_ids != target_map.keys() or seen_ids != author_guides.keys():
        raise ValueError("Author guide sequence and target metadata IDs do not match")
    return gex_lines, guide_lines


def barcode_aliases():
    rows = read_tsv(ROOT / "probe_barcodes.tsv")
    expected = [
        f"{prefix}{number:03d}" for prefix in ("BC", "CR") for number in range(1, 17)
    ]
    aliases = [row["bc_alias"] for row in rows] + [row["cr_alias"] for row in rows]
    if aliases != expected:
        raise ValueError(
            "Probe barcode aliases must be ordered BC001-BC016, CR001-CR016"
        )
    for row in rows:
        if (
            len(row["bc_sequence"]) != 8
            or len(row["cr_sequence"]) != 8
            or set(row["bc_sequence"] + row["cr_sequence"]) - set("ACGT")
        ):
            raise ValueError(f"Invalid probe barcode sequence: {row}")
    if len({row["bc_sequence"] for row in rows}) != 16:
        raise ValueError("RNA probe barcode sequences must be unique")
    if len({row["cr_sequence"] for row in rows}) != 16:
        raise ValueError("CRISPR probe barcode sequences must be unique")
    return rows


def generated_files(gex, guide, gex_lanes, guide_lanes, gex_features, guide_features):
    barcode_aliases()
    author = author_samples()
    coverage = []
    for pool, (_, first, last) in POOL_RANGES.items():
        for lane in range(first, last + 1):
            key = (pool, lane)
            gex_n, guide_n = gex_lanes[key], guide_lanes[key]
            coverage.append(
                {
                    "pool_id": pool,
                    "lane_id": lane_id(pool, lane),
                    "gex_source_pairs": gex_n,
                    "guide_sra_runs": guide_n,
                    "guide_status": "available" if guide_n else "no_archived_guide_sra",
                }
            )
    coverage_by_pool = {}
    for pool in {row["pool_id"] for row in author}:
        rows = [row for row in coverage if row["pool_id"] == pool]
        coverage_by_pool[pool] = (
            ";".join(row["lane_id"] for row in rows if row["guide_sra_runs"]),
            ";".join(row["lane_id"] for row in rows if not row["guide_sra_runs"]),
        )
    samples = []
    for sample in author:
        pool = sample["pool_id"]
        available, missing = coverage_by_pool[pool]
        samples.append(
            {
                **sample,
                "pool_lanes": ";".join(
                    row["lane_id"] for row in coverage if row["pool_id"] == pool
                ),
                "guide_available_lanes": available,
                "guide_missing_lanes": missing,
            }
        )
    result = {
        "gene_probes.tsv": ("\n".join(gex_features) + "\n").encode(),
        "gene_count_features.tsv": (
            "\n".join(
                "\t".join((line.split("\t", 1)[0], line.split("\t", 1)[0]))
                for line in gex_features
            )
            + "\n"
        ).encode(),
        "guide_features.tsv": ("\n".join(guide_features) + "\n").encode(),
        "guide_count_features.tsv": (
            "\n".join(
                "\t".join((line.split("\t", 1)[0], line.split("\t", 1)[0]))
                for line in guide_features
            )
            + "\n"
        ).encode(),
        "samples.tsv": render_tsv(samples, list(samples[0])),
        "guide_coverage.tsv": render_tsv(coverage, list(coverage[0])),
    }
    return result


def render_tsv(rows, fields):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="verify generated files without writing"
    )
    args = parser.parse_args()
    gex, guide, gex_lanes, guide_lanes = validate_sources()
    gex_features, guide_features = feature_tables()
    outputs = generated_files(
        gex, guide, gex_lanes, guide_lanes, gex_features, guide_features
    )
    mismatches = []
    for name, data in outputs.items():
        path = ROOT / name
        if args.check:
            if not path.is_file() or path.read_bytes() != data:
                mismatches.append(name)
        else:
            path.write_bytes(data)
            print(
                f"{name}: {len(data):,} bytes, sha256={hashlib.sha256(data).hexdigest()}"
            )
    if mismatches:
        raise SystemExit("Generated files differ: " + ", ".join(mismatches))
    if args.check:
        print(
            f"Validated {len(gex)} pinned GEX pairs, {len(guide)} guide SRA sources, and {len(outputs)} generated files"
        )


if __name__ == "__main__":
    main()
