import pandas as pd
import sys
import os


def generate_features_tsv(xlsx_path, output_tsv):
    if not os.path.exists(xlsx_path):
        print(f"Error: Could not find file {xlsx_path}")
        sys.exit(1)

    print(f"Reading sheet 'ST20' from {xlsx_path}...")
    df = pd.read_excel(xlsx_path, sheet_name="ST20")

    # Extract Guide A
    df_A = df[["sgID_A", "targeting sequence A"]].rename(
        columns={"sgID_A": "id", "targeting sequence A": "seq"}
    )

    # Extract Guide B
    df_B = df[["sgID_B", "targeting sequence B"]].rename(
        columns={"sgID_B": "id", "targeting sequence B": "seq"}
    )

    # Combine, drop missing and duplicates
    features_df = pd.concat([df_A, df_B]).dropna().drop_duplicates()

    # Save as headerless TSV
    features_df.to_csv(output_tsv, sep="\t", index=False, header=False)
    print(f"Successfully wrote {len(features_df)} unique guides to {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_features_nadig.py <path_to_xlsx> <output_tsv>")
        sys.exit(1)
    generate_features_tsv(sys.argv[1], sys.argv[2])
