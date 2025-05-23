#!/usr/bin/env python3

import argparse
import csv
import json
import sys


def main():
    parser = argparse.ArgumentParser(description="Process Perturb-Seq CSV data")
    parser.add_argument("--input-filename", required=True, help="Input CSV file")
    parser.add_argument("--output-filename", required=True, help="Output JSONL file")
    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="Maximum padj value to include (default: 0.05)",
    )
    parser.add_argument(
        "--log2fc-lower",
        type=float,
        default=-1.0,
        help="Lower log2fc threshold (default: -1.0)",
    )
    parser.add_argument(
        "--log2fc-higher",
        type=float,
        default=1.0,
        help="Higher log2fc threshold (default: 1.0)",
    )
    args = parser.parse_args()

    # Validate thresholds
    if args.padj_threshold <= 0 or args.padj_threshold > 1:
        print("Error: padj_threshold must be between 0 and 1", file=sys.stderr)
        sys.exit(1)

    if args.log2fc_lower >= args.log2fc_higher:
        print("Error: log2fc_lower must be less than log2fc_higher", file=sys.stderr)
        sys.exit(1)

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
                        # Filter: only keep rows where padj <= threshold
                        padj = float(row["padj"])
                        if padj > args.padj_threshold:
                            records_skipped += 1
                            continue

                        # Check log2fc thresholds
                        log2fc = float(row["log2fc"])
                        if not (
                            log2fc < args.log2fc_lower or log2fc > args.log2fc_higher
                        ):
                            records_skipped += 1
                            continue

                        # Create record with type conversion
                        record = {
                            "record_id": f"{row['study_id']}_{row['perturbation']}_{row['gene']}",
                            "study_id": row["study_id"],
                            "perturbation": row["perturbation"],
                            "gene": row["gene"],
                            "log2fc": log2fc,
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
                print(f"  Records skipped: {records_skipped}")
                print(f"  Filter criteria:")
                print(f"    padj <= {args.padj_threshold}")
                print(f"    log2fc < {args.log2fc_lower} or > {args.log2fc_higher}")

    except FileNotFoundError:
        print(f"Error: Input file '{args.input_filename}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
