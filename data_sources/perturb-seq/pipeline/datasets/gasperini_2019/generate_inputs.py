#!/usr/bin/env python3
"""Generate Gasperini 2019 BAM sample and gene-level guide inputs."""

import csv
import gzip
import io
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


DATASET_DIR = Path(__file__).resolve().parent
FEATURES_TSV = DATASET_DIR / "gasperini_2019_atscale_features.tsv"
SAMPLES_TSV = DATASET_DIR / "gasperini_2019_atscale_samples.tsv"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
GUIDE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/"
    "GSE120861_grna_groups.at_scale.txt.gz"
)
RUN_KEY = re.compile(r"at_scale_screen\.([12][AB]_[1-8])_")


def read_ena_runs():
    query = urlencode(
        {
            "accession": "PRJNA494734",
            "result": "read_run",
            "fields": "run_accession,submitted_ftp",
            "format": "tsv",
        }
    )
    with urlopen(f"{ENA_URL}?{query}", timeout=60) as response:
        rows = list(csv.DictReader(io.TextIOWrapper(response), delimiter="\t"))

    pairs = {}
    for row in rows:
        match = RUN_KEY.search(row["submitted_ftp"])
        if not match:
            continue
        key = match.group(1)
        kind = (
            "guide"
            if "_gRNA_" in row["submitted_ftp"]
            else "mrna" if "_SI_" in row["submitted_ftp"] else ""
        )
        if not kind:
            continue
        pairs.setdefault(key, {}).setdefault(kind, []).append(row["run_accession"])

    expected = {
        f"{plate}_{well}" for plate in ("1A", "1B", "2A", "2B") for well in range(1, 9)
    }
    if set(pairs) != expected or any(
        set(value) != {"mrna", "guide"} for value in pairs.values()
    ):
        raise ValueError(f"Expected 32 complete GEX/guide pairs, found {pairs}")
    return pairs


def write_samples():
    pairs = read_ena_runs()
    with SAMPLES_TSV.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "mRNA_srrs",
                "sgRNA_srrs",
                "guide_feature_offset",
                "guide_feature_length",
            ]
        )
        for sample_id in sorted(pairs):
            pair = pairs[sample_id]
            writer.writerow(
                [
                    sample_id,
                    ";".join(sorted(f"BAM:{run}" for run in pair["mrna"])),
                    ";".join(sorted(f"BAM:{run}" for run in pair["guide"])),
                    23,
                    20,
                ]
            )


def write_features():
    raw = gzip.decompress(urlopen(GUIDE_URL, timeout=60).read()).decode()
    rows = []
    seen = set()
    for line in raw.splitlines():
        if not line.strip():
            continue
        target, sequence = line.split("\t")[:2]
        sequence = sequence.strip().upper()
        if not re.fullmatch(r"[ACGT]{20}", sequence):
            raise ValueError(f"Invalid guide sequence: {sequence!r}")
        if sequence in seen:
            raise ValueError(f"Duplicate guide sequence: {sequence}")
        seen.add(sequence)

        if target.endswith("_TSS"):
            label = f"{target[:-4]}__{sequence}"
        elif target == "bassik_mch" or target.startswith(
            ("random_", "scrambled_", "pos_control_")
        ):
            label = f"non-targeting_{target}"
        else:
            continue
        rows.append((sequence, label))

    with FEATURES_TSV.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
    print(f"Wrote {len(rows)} features to {FEATURES_TSV}")


if __name__ == "__main__":
    write_features()
    write_samples()
