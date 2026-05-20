#!/usr/bin/env python3

from io import StringIO
from pathlib import Path
from urllib.request import urlopen
import re
import sys

import pandas as pd


DATASET_DIR = Path(__file__).resolve().parent
GUIDE_XLSX = (
    DATASET_DIR
    / "../../../../../data_exploration/Perturbseq/supplementary/nadig_2025_guide_info.xlsx"
).resolve()
ENA_ACCESSION = "SAMN40972597"
ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    f"?accession={ENA_ACCESSION}"
    "&result=read_run"
    "&fields=run_accession,library_name,fastq_ftp"
    "&format=tsv"
)
FEATURES_TSV = DATASET_DIR / "jurkat_features.tsv"
SAMPLES_TSV = DATASET_DIR / "jurkat_samples.tsv"


def generate_features_tsv(xlsx_path, output_tsv):
    if not xlsx_path.exists():
        print(f"Error: could not find guide XLSX at {xlsx_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading sheet 'ST20' from {xlsx_path}")
    df = pd.read_excel(xlsx_path, sheet_name="ST20")

    df_a = df[["targeting sequence A", "sgID_A"]].rename(
        columns={"targeting sequence A": "seq", "sgID_A": "id"}
    )
    df_b = df[["targeting sequence B", "sgID_B"]].rename(
        columns={"targeting sequence B": "seq", "sgID_B": "id"}
    )

    features_df = pd.concat([df_a, df_b]).dropna()
    features_df["seq"] = features_df["seq"].astype(str).str.strip().str.upper()
    features_df["id"] = (
        features_df["id"].astype(str).str.strip().str.replace(",", "-", regex=False)
    )

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


def fetch_ena_metadata(url):
    print(f"Fetching ENA metadata for {ENA_ACCESSION}")
    with urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")
    return pd.read_csv(StringIO(text), sep="\t")


def parse_jurkat_library(library_name):
    match = re.search(r"jurkat_(mRNA|sgRNA)_(\d+)_", str(library_name))
    if not match:
        return None, None
    return match.group(1), match.group(2)


def generate_samples_tsv(output_tsv):
    metadata = fetch_ena_metadata(ENA_URL)
    metadata[["modality", "sample_id"]] = metadata["library_name"].apply(
        lambda value: pd.Series(parse_jurkat_library(value))
    )
    metadata = metadata.dropna(subset=["sample_id"])

    samples = {}
    for _, row in metadata.iterrows():
        sample_id = str(row["sample_id"])
        modality = row["modality"]
        run_accession = row["run_accession"]

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

    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)
    print(f"Wrote {len(rows)} samples to {output_tsv}")


def main():
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    generate_features_tsv(GUIDE_XLSX, FEATURES_TSV)
    generate_samples_tsv(SAMPLES_TSV)


if __name__ == "__main__":
    main()
