#!/usr/bin/env python3

import csv
import requests
import sys
from pathlib import Path
from collections import defaultdict

ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def bytes_to_tb(num_bytes):
    """Convert bytes to terabytes (TB)."""
    return num_bytes / (1024**4)


def get_study_summary_with_samples(study_accession):
    """
    Query ENA API and return:
    - total_bytes
    - total_files
    - sample_stats: dict {sample_accession: (bytes, file_count)}
    """
    params = {
        "accession": study_accession,
        "result": "read_run",
        "fields": "fastq_bytes,sample_accession",
        "format": "tsv",
    }

    response = requests.get(ENA_API_URL, params=params)
    response.raise_for_status()

    lines = response.text.strip().split("\n")
    if len(lines) <= 1:
        return 0, 0, {}

    reader = csv.DictReader(lines, delimiter="\t")

    total_bytes = 0
    total_files = 0
    sample_stats = defaultdict(lambda: [0, 0])  # {sample: [bytes, file_count]}

    for row in reader:
        sample = row.get("sample_accession") or "UNKNOWN"

        if row.get("fastq_bytes"):
            sizes = row["fastq_bytes"].split(";")
            for size in sizes:
                size = size.strip()
                if size:
                    size_int = int(size)
                    total_bytes += size_int
                    total_files += 1
                    sample_stats[sample][0] += size_int
                    sample_stats[sample][1] += 1

    return total_bytes, total_files, sample_stats


def main(datasets_tsv):
    with open(datasets_tsv) as f:
        reader = csv.reader(f, delimiter="\t")

        print("study_accession\ttotal_size_TB\ttotal_fastq_files")

        for row in reader:
            if not row:
                continue

            study_accession = row[0].strip()
            if not study_accession:
                continue

            try:
                total_bytes, total_files, sample_stats = get_study_summary_with_samples(
                    study_accession
                )

                total_tb = bytes_to_tb(total_bytes)

                print(f"{study_accession}\t{total_tb:.3f}\t{total_files}")

                # Per-sample breakdown
                for idx, (sample, stats) in enumerate(sorted(sample_stats.items()), 1):
                    sample_bytes, sample_files = stats
                    sample_tb = bytes_to_tb(sample_bytes)
                    print(f"    #{idx} {sample} {sample_tb:.3f} {sample_files}")

            except Exception as e:
                print(f"{study_accession}\tERROR\tERROR ({e})", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} datasets.tsv")
        sys.exit(1)

    datasets_file = Path(sys.argv[1])
    if not datasets_file.exists():
        print(f"File not found: {datasets_file}")
        sys.exit(1)

    main(datasets_file)
