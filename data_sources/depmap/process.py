import argparse
import pandas as pd
import json


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process gene dependency data.")
    parser.add_argument(
        "--gene-dependency", required=True, help="Path to the gene dependency CSV file."
    )
    parser.add_argument(
        "--common-essentials",
        required=True,
        help="Path to the common essentials text file.",
    )
    parser.add_argument(
        "--model-metadata", required=True, help="Path to the model metadata CSV file."
    )
    parser.add_argument(
        "--mavedb-metadata",
        required=True,
        help="Path to the MaveDB metadata JSONL file.",
    )
    parser.add_argument(
        "--output-filename", required=True, help="Path to the output JSONL file."
    )
    args = parser.parse_args()

    # Load MaveDB metadata and extract the set of normalisedGeneName
    mavedb_genes = set()
    with open(args.mavedb_metadata, "r") as f:
        for line in f:
            record = json.loads(line)
            if "normalisedGeneName" in record:
                mavedb_genes.add(record["normalisedGeneName"])

    # Load gene dependency CSV
    gene_df = pd.read_csv(args.gene_dependency)

    # Set index to the first column (cancer screen ID)
    gene_df.set_index(gene_df.columns[0], inplace=True)

    # Remove rows and columns with missing values
    # Order is important, because we need to drop the small subset of genes which have a lot of missing values.
    # (Otherwise, we will lose a lot of good screens.)
    gene_df.dropna(
        axis=1, how="any", inplace=True
    )  # Drop columns with any missing values
    gene_df.dropna(axis=0, how="any", inplace=True)  # Drop rows with any missing values

    # Remove columns present in common essentials
    with open(args.common_essentials, "r") as f:
        common_essentials = set(line.strip() for line in f)
    gene_df = gene_df.drop(
        columns=[col for col in gene_df.columns if col in common_essentials]
    )

    # Rename columns from "GENE (123)" to "GENE"
    gene_df.columns = [col.split(" ")[0] for col in gene_df.columns]

    # Create a new column with genes >= 0.95, sorted alphabetically, and add xrefs to MaveDB.
    gene_df["high_dependency_genes"] = gene_df.apply(
        lambda row: sorted(
            [
                {"name": gene, "xref": "MaveDB" if gene in mavedb_genes else ""}
                for gene in row.index[row >= 0.95].tolist()
            ],
            key=lambda x: x["name"],
        ),
        axis=1,
    )
    gene_df = gene_df[["high_dependency_genes"]]  # Drop all original columns

    # Load model metadata and join with gene_df
    model_df = pd.read_csv(args.model_metadata)
    model_df.set_index("ModelID", inplace=True)
    df = gene_df.join(model_df, how="inner")

    # Reset index to make ModelID a column
    df = df.reset_index().rename(columns={"index": "ModelID"})

    # Keep only specified columns
    columns_to_keep = [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "Age",
        "AgeCategory",
        "Sex",
        "PrimaryOrMetastasis",
        "SampleCollectionSite",
        "CatalogNumber",
        "ModelType",
        "high_dependency_genes",
    ]
    df = df[columns_to_keep]

    # Add the "xref" column
    df["xref"] = df["high_dependency_genes"].apply(
        lambda genes: (
            "In MaveDB"
            if any(gene["xref"] == "MaveDB" for gene in genes)
            else "Not in MaveDB"
        )
    )

    # Replace null values with "Unknown"
    # In source data, missing values are sometimes null and sometimes "Unknown". This
    # step fixes the inconsistency.
    na_values = {
        # Categorical fields; in alignment with other values in the field.
        "SampleCollectionSite": "Unavailable",
        "PrimaryOrMetastasis": "Unavailable",
        # This is not categorical, so better to simply keep it empty when null.
        "CatalogNumber": "",
    }
    df.fillna(value=na_values, inplace=True)

    # Save to JSONL format
    df.to_json(args.output_filename, orient="records", lines=True)


if __name__ == "__main__":
    main()
