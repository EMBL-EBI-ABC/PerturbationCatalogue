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
    / "../../../../../data_exploration/Perturbseq/supplementary/replogle_2022_guide_info.xlsx"
).resolve()
ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DATASETS = [
    {
        "name": "replogle_2022_k562_essential_normalized",
        "accession": "SAMN28561243",
        "guide_sheet": "TabB_K562_day6_library",
    },
    {
        "name": "replogle_2022_rpe1_essential_normalized",
        "accession": "SAMN28561244",
        "guide_sheet": "TabC_RPE1_day7_library",
    },
    {
        "name": "replogle_2022_k562_gw_normalized",
        "accession": "SAMN28561242",
        "guide_sheet": "TabA_K562_day8_library",
    },
]


def output_path(dataset_name, suffix):
    return DATASET_DIR / f"{dataset_name}_{suffix}.tsv"


def natural_key(value):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]


def read_feature_table(xlsx_path, sheet_name):
    if not xlsx_path.exists():
        print(f"Error: could not find guide XLSX at {xlsx_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading sheet '{sheet_name}' from {xlsx_path}")
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

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
    return features_df[["seq", "id"]]


def write_feature_tables():
    feature_tables = {
        dataset["name"]: read_feature_table(GUIDE_XLSX, dataset["guide_sheet"])
        for dataset in DATASETS
    }

    signatures = {
        name: set(table.itertuples(index=False, name=None))
        for name, table in feature_tables.items()
    }
    unique_signatures = {frozenset(signature) for signature in signatures.values()}

    if len(unique_signatures) == 1:
        output_tsv = DATASET_DIR / "features.tsv"
        next(iter(feature_tables.values())).to_csv(
            output_tsv, sep="\t", index=False, header=False
        )
        print(
            f"Wrote {len(next(iter(feature_tables.values())))} guides to {output_tsv}"
        )
        return

    print(
        "Guide whitelists differ across Replogle 2022 datasets; writing one per dataset"
    )
    for dataset_name, features_df in feature_tables.items():
        output_tsv = output_path(dataset_name, "features")
        features_df.to_csv(output_tsv, sep="\t", index=False, header=False)
        print(f"Wrote {len(features_df)} guides to {output_tsv}")


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
    with urlopen(ena_url(accession), timeout=60) as response:
        text = response.read().decode("utf-8")
    return pd.read_csv(StringIO(text), sep="\t")


def parse_replogle_library(library_name):
    name = str(library_name).strip()
    match = re.match(
        r"^(?P<prefix>.+?)_(?P<modality>mRNA|sgRNA)_(?P<sample>.+?)"
        r"(?:_S\d+)?(?:_L\d+)?$",
        name,
        flags=re.IGNORECASE,
    )
    if not match:
        return None, None

    modality = {"mrna": "mRNA", "sgrna": "sgRNA"}[match.group("modality").lower()]
    sample_id = f"{match.group('prefix')}_{match.group('sample')}"
    return modality, sample_id


def generate_samples_tsv(dataset):
    metadata = fetch_ena_metadata(dataset["accession"])
    metadata[["modality", "sample_id"]] = metadata["library_name"].apply(
        lambda value: pd.Series(parse_replogle_library(value))
    )

    unmatched_with_modality = metadata[
        metadata["sample_id"].isna()
        & metadata["library_name"]
        .astype(str)
        .str.contains(r"mRNA|sgRNA", case=False, regex=True)
    ]
    if not unmatched_with_modality.empty:
        examples = ", ".join(
            unmatched_with_modality["library_name"].astype(str).head(5).tolist()
        )
        print(
            f"Error: found library names with mRNA/sgRNA that did not parse: {examples}",
            file=sys.stderr,
        )
        sys.exit(1)

    ignored_count = int(metadata["sample_id"].isna().sum())
    metadata = metadata.dropna(subset=["sample_id"])

    samples = {}
    for _, row in metadata.iterrows():
        sample_id = str(row["sample_id"])
        modality = row["modality"]
        run_accession = str(row["run_accession"])

        samples.setdefault(sample_id, {"mRNA": set(), "sgRNA": set()})
        samples[sample_id][modality].add(run_accession)

    rows = []
    skipped_incomplete = []
    for sample_id in sorted(samples, key=natural_key):
        modalities = samples[sample_id]
        if not modalities["mRNA"] or not modalities["sgRNA"]:
            skipped_incomplete.append(sample_id)
            continue

        rows.append(
            {
                "sample_id": sample_id,
                "mRNA_srrs": ";".join(sorted(modalities["mRNA"])),
                "sgRNA_srrs": ";".join(sorted(modalities["sgRNA"])),
            }
        )

    output_tsv = output_path(dataset["name"], "samples")
    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)
    print(
        f"Wrote {len(rows)} {dataset['name']} samples to {output_tsv}; "
        f"ignored {ignored_count} non-mRNA/sgRNA runs"
    )
    if skipped_incomplete:
        print(
            f"Skipped {len(skipped_incomplete)} incomplete samples without both modalities",
            file=sys.stderr,
        )


def main():
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    write_feature_tables()
    for dataset in DATASETS:
        generate_samples_tsv(dataset)


if __name__ == "__main__":
    main()
