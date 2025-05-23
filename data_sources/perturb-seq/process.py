#!/usr/bin/env python3

import argparse
import csv
import json
import sys


def main():
    parser = argparse.ArgumentParser(description="Process Perturb-Seq CSV data")
    parser.add_argument("--input-filename", required=True, help="Input CSV file")
    parser.add_argument("--output-filename", required=True, help="Output JSONL file")
    args = parser.parse_args()

    try:
        with open(args.input_filename, "r", newline="", encoding="utf-8") as infile:
            with open(args.output_filename, "w", encoding="utf-8") as outfile:
                reader = csv.DictReader(infile)

                # Validate required columns exist
                required_columns = [
                    "study_id",
                    "perturbation",
                    "gene",
                    "log2fc",
                    "pvalue",
                    "padj",
                    "mean_control",
                    "mean_perturbed",
                ]
                missing_columns = [
                    col for col in required_columns if col not in reader.fieldnames
                ]
                if missing_columns:
                    print(
                        f"Error: Missing required columns: {missing_columns}",
                        file=sys.stderr,
                    )
                    sys.exit(1)

                records_written = 0
                records_skipped = 0

                for row_num, row in enumerate(
                    reader, start=2
                ):  # Start at 2 because of header
                    try:
                        # Filter: only keep rows where padj <= 0.05
                        padj = float(row["padj"])
                        if padj > 0.05:
                            records_skipped += 1
                            continue

                        # Create record with type conversion
                        record = {
                            "record_id": f"{row['study_id']}_{row['perturbation']}_{row['gene']}",
                            "study_id": row["study_id"],
                            "perturbation": row["perturbation"],
                            "gene": row["gene"],
                            "log2fc": float(row["log2fc"]),
                            "pvalue": float(row["pvalue"]),
                            "padj": padj,
                            "mean_control": float(row["mean_control"]),
                            "mean_perturbed": float(row["mean_perturbed"]),
                        }

                        # Write as JSON line
                        outfile.write(json.dumps(record) + "\n")
                        records_written += 1

                    except ValueError as e:
                        print(f"Error processing row {row_num}: {e}", file=sys.stderr)
                        print(f"Row data: {row}", file=sys.stderr)
                        continue
                    except Exception as e:
                        print(
                            f"Unexpected error processing row {row_num}: {e}",
                            file=sys.stderr,
                        )
                        continue

                print(f"Processing complete:")
                print(f"  Records written: {records_written}")
                print(f"  Records skipped (padj > 0.05): {records_skipped}")

    except FileNotFoundError:
        print(f"Error: Input file '{args.input_filename}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
