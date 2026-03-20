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
    df_A = df[["targeting sequence A", "sgID_A"]].rename(
        columns={"targeting sequence A": "seq", "sgID_A": "id"}
    )

    # Extract Guide B
    df_B = df[["targeting sequence B", "sgID_B"]].rename(
        columns={"targeting sequence B": "seq", "sgID_B": "id"}
    )

    # Combine all guides
    features_df = pd.concat([df_A, df_B]).dropna()

    # Resolve duplicate guide sequences by grouping by 'seq' and joining the 'id's
    features_df = features_df.groupby("seq", as_index=False).agg(
        {"id": lambda x: ";".join(sorted(set(x)))}
    )

    # Ensure sequence is first column, id is second
    features_df = features_df[["seq", "id"]]

    # Save as headerless TSV
    features_df.to_csv(output_tsv, sep="\t", index=False, header=False)
    print(f"Successfully wrote {len(features_df)} unique guides to {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_features_nadig.py <path_to_xlsx> <output_tsv>")
        sys.exit(1)
    generate_features_tsv(sys.argv[1], sys.argv[2])
