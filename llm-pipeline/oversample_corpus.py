"""
Oversample essential CRISPR records before gene-level split.
Increases essential:neutral ratio to improve essential recall in fine-tuning.
"""
import json
import random
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input corpus JSONL")
    parser.add_argument("--output", required=True, help="Output oversampled corpus JSONL")
    parser.add_argument("--oversample_factor", type=int, default=5, help="How many times to repeat essential records (default: 5)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)

    records = [json.loads(l) for l in open(args.input) if l.strip()]

    essential = [r for r in records if r.get("metadata", {}).get("fitness_class") == "essential"]
    non_essential = [r for r in records if r.get("metadata", {}).get("fitness_class") != "essential"]

    print(f"Original corpus: {len(records)} records")
    print(f"  Essential: {len(essential)}")
    print(f"  Non-essential: {len(non_essential)}")

    oversampled = records + essential * (args.oversample_factor - 1)
    random.shuffle(oversampled)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for r in oversampled:
            f.write(json.dumps(r) + "\n")

    essential_after = sum(1 for r in oversampled if r.get("metadata", {}).get("fitness_class") == "essential")
    print(f"\nOversampled corpus: {len(oversampled)} records")
    print(f"  Essential: {essential_after} ({100*essential_after/len(oversampled):.1f}%)")
    print(f"  Non-essential: {len(oversampled) - essential_after} ({100*(len(oversampled)-essential_after)/len(oversampled):.1f}%)")
    print(f"\nSaved to {args.output}")

if __name__ == "__main__":
    main()
