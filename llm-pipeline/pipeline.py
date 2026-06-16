import json
import logging
import argparse
from pathlib import Path

from prepare_training_data import fetch_and_process_crispr

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


def load_corpus(jsonl_path):
    """
    Load training records from a JSONL file.

    Parameters
    ----------
    jsonl_path : str or Path

    Returns
    -------
    list of training record dicts
    """
    records = []
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    log.info(f"Loaded {len(records)} records from {jsonl_path}")
    return records


def run_demo(dataset_id, output_path):
    """
    Fetch real data and show example training records.

    Parameters
    ----------
    dataset_id : str
    output_path : str
    """
    log.info(f"Running demo on {dataset_id}...")

    records, df = fetch_and_process_crispr(
        dataset_id=dataset_id,
        output_path=output_path,
        max_records=200,
    )

    if not records:
        log.error("No records returned — check API connectivity")
        return

    print("\n" + "=" * 60)
    print("PERTURBATION CATALOGUE — PIPELINE DEMO")
    print(f"Dataset: {dataset_id}")
    print("=" * 60)

    print(f"\nDataset shape: {df.shape}")
    print(f"\nFitness classification:")
    print(df["fitness_class"].value_counts().to_string())

    print(f"\nTop 5 essential genes (strongest depletion):")
    essential = df[df["fitness_class"] == "essential"].nsmallest(5, "effect_score")[
        ["gene", "effect_score", "effect_score_zscore", "cell_line"]
    ]
    print(essential.to_string(index=False))

    print(f"\nExample training record:")
    r = records[0]
    print(f"\nINSTRUCTION:\n{r['instruction']}")
    print(f"\nINPUT:\n{r['input']}")
    print(f"\nOUTPUT:\n{r['output']}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Load training data and run analysis pipeline"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run demo fetching real data from Catalogue API"
    )
    parser.add_argument(
        "--dataset_id",
        type=str,
        default="biogrid_5",
        help="Dataset ID to use for demo"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/demo_output.jsonl",
        help="Output path for demo records"
    )
    parser.add_argument(
        "--load",
        type=str,
        default=None,
        help="Path to existing JSONL corpus to load"
    )
    args = parser.parse_args()

    if args.demo:
        run_demo(args.dataset_id, args.output)
    elif args.load:
        records = load_corpus(args.load)
        print(f"Loaded {len(records)} training records")
        print(f"Example record keys: {list(records[0].keys())}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
