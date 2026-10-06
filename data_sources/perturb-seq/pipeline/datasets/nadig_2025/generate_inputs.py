#!/usr/bin/env python3

from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import re
import sys

import pandas as pd


DATASET_DIR = Path(__file__).resolve().parent
GUIDE_XLSX = (
    DATASET_DIR
    / "../../../../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx"
).resolve()
FEATURES_TSV = DATASET_DIR / "features.tsv"
ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DATASETS = [
    {
        "name": "jurkat",
        "accession": "SAMN40972597",
        "library_prefix": "jurkat",
        "samples_tsv": DATASET_DIR / "jurkat_samples.tsv",
    },
    {
        "name": "hepg2",
        "accession": "SAMN40972598",
        "library_prefix": "hepg2",
        "samples_tsv": DATASET_DIR / "hepg2_samples.tsv",
    },
]


def guide_id_with_ensg(guide_id, ensg):
    guide_id = str(guide_id).strip().replace(",", "-")
    ensg = str(ensg).strip().split(".", 1)[0]
    return f"{guide_id}__{ensg}" if ensg.startswith("ENSG") else guide_id


def generate_features_tsv(xlsx_path, output_tsv):
    if not xlsx_path.exists():
        print(f"Error: could not find guide XLSX at {xlsx_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading sheet 'ST20' from {xlsx_path}")
    df = pd.read_excel(xlsx_path, sheet_name="ST20")

    df_a = df[["targeting sequence A", "sgID_A", "ensembl gene id"]].rename(
        columns={
            "targeting sequence A": "seq",
            "sgID_A": "id",
            "ensembl gene id": "ensg",
        }
    )
    df_b = df[["targeting sequence B", "sgID_B", "ensembl gene id"]].rename(
        columns={
            "targeting sequence B": "seq",
            "sgID_B": "id",
            "ensembl gene id": "ensg",
        }
    )

    features_df = pd.concat([df_a, df_b]).dropna(subset=["seq", "id"])
    features_df["seq"] = features_df["seq"].astype(str).str.strip().str.upper()
    features_df["id"] = (
        features_df["id"].astype(str).str.strip().str.replace(",", "-", regex=False)
    )
    features_df["ensg"] = (
        features_df["ensg"].astype(str).str.strip().str.split(".", n=1).str[0]
    )
    features_df = features_df[
        features_df["id"].str.startswith("non-targeting")
        | features_df["ensg"].str.startswith("ENSG")
    ]
    features_df["id"] = [
        guide_id_with_ensg(guide_id, ensg)
        for guide_id, ensg in zip(features_df["id"], features_df["ensg"])
    ]

    invalid = features_df[~features_df["seq"].str.fullmatch(r"[ACGT]{20}")]
    if not invalid.empty:
        examples = ", ".join(invalid["seq"].head(5).tolist())
        print(
            f"Error: expected 20 bp A/C/G/T guide sequences; examples: {examples}",
            file=sys.stderr,
        )
        sys.exit(1)

    features_df = features_df.groupby("seq", as_index=False).agg(
        {"id": lambda values: ";".join(sorted(set(values)))}
    )
    features_df = features_df[["seq", "id"]]
    features_df.to_csv(output_tsv, sep="\t", index=False, header=False)
    print(f"Wrote {len(features_df)} unique guides to {output_tsv}")


def ena_url(accession):
    query = urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,library_name,fastq_ftp",
            "format": "tsv",
        }
    )
    return f"{ENA_API_URL}?{query}"


def fetch_ena_metadata(accession):
    print(f"Fetching ENA metadata for {accession}")
    url = ena_url(accession)
    with urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")
    return pd.read_csv(StringIO(text), sep="\t")


def parse_library(library_name, library_prefix):
    pattern = rf"{re.escape(library_prefix)}_(mRNA|sgRNA)_(\d+)(?:_|$)"
    match = re.search(pattern, str(library_name), flags=re.IGNORECASE)
    if not match:
        return None, None
    modality = {"mrna": "mRNA", "sgrna": "sgRNA"}[match.group(1).lower()]
    return modality, match.group(2)


def generate_samples_tsv(dataset):
    metadata = fetch_ena_metadata(dataset["accession"])
    metadata[["modality", "sample_id"]] = metadata["library_name"].apply(
        lambda value: pd.Series(parse_library(value, dataset["library_prefix"]))
    )
    metadata = metadata.dropna(subset=["sample_id"])

    samples = {}
    for _, row in metadata.iterrows():
        sample_id = str(row["sample_id"])
        modality = row["modality"]
        run_accession = str(row["run_accession"])

        samples.setdefault(sample_id, {"mRNA": set(), "sgRNA": set()})
        samples[sample_id][modality].add(run_accession)

    rows = []
    for sample_id in sorted(samples, key=lambda value: int(value)):
        rows.append(
            {
                "sample_id": sample_id,
                "mRNA_srrs": ";".join(sorted(samples[sample_id]["mRNA"])),
                "sgRNA_srrs": ";".join(sorted(samples[sample_id]["sgRNA"])),
            }
        )

    output_tsv = dataset["samples_tsv"]
    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)
    print(f"Wrote {len(rows)} {dataset['name']} samples to {output_tsv}")


def main():
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    generate_features_tsv(GUIDE_XLSX, FEATURES_TSV)
    for dataset in DATASETS:
        generate_samples_tsv(dataset)


if __name__ == "__main__":
    main()
