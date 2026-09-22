#!/usr/bin/env python3
"""Generate Norman 2019 KITE features and paired GEM-group inputs."""

from io import StringIO
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd


DATASET_DIR = Path(__file__).resolve().parent
GUIDE_XLSX = DATASET_DIR / "norman_2019_guide_info.xlsx"
FEATURES_TSV = DATASET_DIR / "features.tsv"
SAMPLES_TSV = DATASET_DIR / "norman_2019_raw_samples.tsv"
ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
PROJECT = "PRJNA551220"
GUIDE_SOURCE_URL = (
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/6746554/bin/"
    "NIHMS1045467-supplement-Table_S2.xlsx"
)


def target_label(row):
    genes = [
        str(row[column]).strip()
        for column in ("gene_A", "gene_B")
        if not str(row[column]).strip().startswith("NegCtrl")
    ]
    return ";".join(genes) if genes else "non-targeting"


def generate_features():
    if not GUIDE_XLSX.is_file():
        raise FileNotFoundError(GUIDE_XLSX)

    table = pd.read_excel(GUIDE_XLSX, sheet_name="Perturbseq_sgRNA_info")
    required = {
        "gene_A",
        "gene_B",
        "GBC",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing guide-table columns: {sorted(missing)}")

    features = []
    for _, row in table.iterrows():
        label = target_label(row)
        for sequence in str(row["GBC"]).split(";"):
            sequence = sequence.strip().upper()
            if not re.fullmatch(r"[ACGT]{18}", sequence):
                raise ValueError(f"Invalid Norman guide barcode: {sequence!r}")
            features.append((sequence, label))

    features = pd.DataFrame(features, columns=["sequence", "label"])
    if features["sequence"].duplicated().any():
        raise ValueError("Norman guide barcode sequences are not unique")
    features.to_csv(FEATURES_TSV, sep="\t", index=False, header=False)
    print(f"Wrote {len(features)} features to {FEATURES_TSV}")


def fetch_metadata():
    query = urlencode(
        {
            "accession": PROJECT,
            "result": "read_run",
            "fields": "run_accession,sample_title",
            "format": "tsv",
        }
    )
    with urlopen(f"{ENA_API_URL}?{query}", timeout=60) as response:
        return pd.read_csv(StringIO(response.read().decode()), sep="\t")


def parse_library(title):
    title = str(title).strip().lower()
    match = re.search(r"gemgroup\s+(\d+)", title)
    if not match:
        return None
    if title.startswith("sgrna perturb-seq experiment"):
        modality = "mRNA"
    elif title.startswith("barcodes identifying perturbations"):
        modality = "sgRNA"
    else:
        return None
    return f"gemgroup_{match.group(1)}", modality


def generate_samples():
    metadata = fetch_metadata()
    samples = {}
    for _, row in metadata.iterrows():
        parsed = parse_library(row["sample_title"])
        if parsed is None:
            continue
        sample_id, modality = parsed
        samples.setdefault(sample_id, {"mRNA": [], "sgRNA": []})[modality].append(
            str(row["run_accession"])
        )

    expected = {f"gemgroup_{number}" for number in range(1, 9)}
    if set(samples) != expected or any(
        not values["mRNA"] or not values["sgRNA"] for values in samples.values()
    ):
        raise ValueError(f"Expected paired gemgroups 1-8, found {samples}")

    rows = [
        {
            "sample_id": sample_id,
            "mRNA_srrs": ";".join(sorted(values["mRNA"])),
            "sgRNA_srrs": ";".join(sorted(values["sgRNA"])),
        }
        for sample_id, values in sorted(
            samples.items(), key=lambda item: int(item[0].split("_")[-1])
        )
    ]
    pd.DataFrame(rows).to_csv(SAMPLES_TSV, sep="\t", index=False)
    print(f"Wrote {len(rows)} samples to {SAMPLES_TSV}")


if __name__ == "__main__":
    generate_features()
    generate_samples()
