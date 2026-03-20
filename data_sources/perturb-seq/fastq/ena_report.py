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
    - sample_stats: dict {sample_accession: [bytes, file_count, variable_metadata_dict]}
    """
    # Fetch all available read_run fields
    url_rr = (
        "https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run&format=tsv"
    )
    res = requests.get(url_rr)
    res.raise_for_status()
    all_fields = [
        line.split("\t")[0] for line in res.text.strip().split("\n")[1:] if line
    ]

    # Exclude technical file-level prefixes and run/experiment-level redundancies
    exclude_prefixes = [
        "fastq_",
        "bam_",
        "cram_",
        "run_",
        "read_",
        "sra_",
        "submitted_",
    ]
    exclude_fields = [
        "base_count",
        "experiment_accession",
        "experiment_alias",
        "experiment_title",
        "description",
        "library_name",
        "secondary_sample_accession",
        "sample_alias",
        "study_accession",
        "secondary_study_accession",
        "submission_accession",
        "study_alias",
        "study_title",
        "first_created",
        "first_public",
        "last_updated",
        "status",
    ]

    sample_fields = [
        f
        for f in all_fields
        if not any(f.startswith(p) for p in exclude_prefixes)
        and f not in exclude_fields
    ]

    req_fields = (
        ["fastq_bytes"] + sample_fields
        if "fastq_bytes" not in sample_fields
        else sample_fields
    )

    params = {
        "accession": study_accession,
        "result": "read_run",
        "fields": ",".join(req_fields),
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
    sample_stats = defaultdict(lambda: [0, 0])
    sample_metadata = defaultdict(dict)

    for row in reader:
        sample = row.get("sample_accession") or "UNKNOWN"

        for f in sample_fields:
            if f != "sample_accession" and row.get(f):
                sample_metadata[sample][f] = row[f]

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

    # Find common keys across all samples
    all_keys = set()
    for meta in sample_metadata.values():
        all_keys.update(meta.keys())

    # Determine which fields have variable values
    variable_fields = []
    for f in all_keys:
        values = set(sample_metadata[s].get(f, "") for s in sample_metadata)
        if len(values) > 1:
            variable_fields.append(f)

    # Add sorted variable metadata to final stats
    variable_fields.sort()
    final_stats = {}
    for sample, stats in sample_stats.items():
        var_meta = {
            f: sample_metadata[sample].get(f)
            for f in variable_fields
            if sample_metadata[sample].get(f)
        }
        final_stats[sample] = [stats[0], stats[1], var_meta]

    return total_bytes, total_files, final_stats


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
                    sample_bytes, sample_files, var_meta = stats
                    sample_tb = bytes_to_tb(sample_bytes)

                    meta_str = ", ".join(f"{k}: {v}" for k, v in var_meta.items())
                    if meta_str:
                        meta_str = f"[{meta_str}]"

                    print(
                        f"    #{idx} {sample} {sample_tb:.3f} {sample_files} {meta_str}"
                    )

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
