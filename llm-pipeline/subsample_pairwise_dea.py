"""
Subsample pairwise DEA records proportionally across all datasets.

Takes the full pairwise DEA corpus (1.1M+ records across 7 datasets) and
samples down to a target total size, preserving each dataset's relative
share of the total rather than favoring whichever dataset happens to be
largest or listed first.

Usage:

    python3 subsample_pairwise_dea.py \
        --input_dir data/dea_pairwise_q4 \
        --output data/dea_pairwise_q4_subsampled.jsonl \
        --target_total 25000
"""
import json
import argparse
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target_total", type=int, default=25000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    input_files = sorted(Path(args.input_dir).glob("*.jsonl"))
    print(f"Found {len(input_files)} dataset files")

    dataset_records = {}
    total_available = 0
    for f in input_files:
        records = [json.loads(l) for l in open(f) if l.strip()]
        dataset_records[f.stem] = records
        total_available += len(records)
        print(f"  {f.stem}: {len(records)} records")

    print(f"\nTotal available: {total_available}")
    print(f"Target total: {args.target_total}")

    sampled = []
    for dataset_id, records in dataset_records.items():
        share = len(records) / total_available
        n_sample = round(args.target_total * share)
        n_sample = min(n_sample, len(records))
        chosen = rng.sample(records, n_sample)
        sampled.extend(chosen)
        print(f"  {dataset_id}: sampled {n_sample} of {len(records)} "
              f"({100*share:.1f}% of total)")

    rng.shuffle(sampled)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for r in sampled:
            f.write(json.dumps(r) + "\n")

    print(f"\nSaved {len(sampled)} subsampled records to {args.output}")


if __name__ == "__main__":
    main()