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


def gene_level_split(records, train_ratio=0.8, val_ratio=0.1, seed=42):
    """
    Split records by gene to prevent evaluation leakage.

    Records for the same gene are kept together in the same split.
    This forces the model to generalise to unseen genes rather than
    memorising gene-specific patterns.

    Parameters
    ----------
    records : list of dict
        Training records with metadata containing gene field.
    train_ratio : float
        Fraction of genes for training. Default 0.8.
    val_ratio : float
        Fraction of genes for validation. Default 0.1.
        Test gets the remainder (1 - train_ratio - val_ratio).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    tuple of (train_records, val_records, test_records, split_manifest)
    split_manifest is a dict mapping gene -> split label
    """
    import numpy as np

    np.random.seed(seed)

    # Get all unique genes
    all_genes = list({r["metadata"]["gene"] for r in records})
    np.random.shuffle(all_genes)

    n_total = len(all_genes)
    n_train = int(n_total * train_ratio)
    n_val = int(n_total * val_ratio)

    train_genes = set(all_genes[:n_train])
    val_genes = set(all_genes[n_train:n_train + n_val])
    test_genes = set(all_genes[n_train + n_val:])

    train_records = [r for r in records if r["metadata"]["gene"] in train_genes]
    val_records = [r for r in records if r["metadata"]["gene"] in val_genes]
    test_records = [r for r in records if r["metadata"]["gene"] in test_genes]

    split_manifest = {}
    for gene in train_genes:
        split_manifest[gene] = "train"
    for gene in val_genes:
        split_manifest[gene] = "val"
    for gene in test_genes:
        split_manifest[gene] = "test"

    log.info(
        f"Gene-level split: {len(train_genes)} train genes "
        f"({len(train_records)} records), "
        f"{len(val_genes)} val genes ({len(val_records)} records), "
        f"{len(test_genes)} test genes ({len(test_records)} records)"
    )

    return train_records, val_records, test_records, split_manifest


def save_splits(train_records, val_records, test_records, split_manifest, output_dir):
    """
    Save gene-level splits to disk.

    Parameters
    ----------
    train_records, val_records, test_records : list of dict
    split_manifest : dict — gene -> split label
    output_dir : str or Path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for split_name, split_records in [
        ("train", train_records),
        ("val", val_records),
        ("test", test_records),
    ]:
        path = output_dir / f"{split_name}.jsonl"
        with open(path, "w") as f:
            for record in split_records:
                f.write(json.dumps(record) + "\n")
        log.info(f"Saved {len(split_records)} records to {path}")

    manifest_path = output_dir / "split_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(split_manifest, f, indent=2)
    log.info(f"Saved split manifest to {manifest_path}")


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

    subparsers = parser.add_subparsers(dest="command")

    demo_parser = subparsers.add_parser(
        "demo",
        help="Fetch real data from Catalogue API and show example records"
    )
    demo_parser.add_argument(
        "--dataset_id",
        type=str,
        default="biogrid_5",
        help="Dataset ID to use for demo"
    )
    demo_parser.add_argument(
        "--output",
        type=str,
        default="data/demo_output.jsonl",
        help="Output path for demo records"
    )

    split_parser = subparsers.add_parser(
        "split",
        help="Split a training corpus by gene into train/val/test sets"
    )
    split_parser.add_argument(
        "--corpus",
        type=str,
        required=True,
        help="Path to full corpus JSONL file"
    )
    split_parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save train.jsonl, val.jsonl, test.jsonl, split_manifest.json"
    )
    split_parser.add_argument(
        "--train_ratio",
        type=float,
        default=0.8,
        help="Fraction of genes for training (default: 0.8)"
    )
    split_parser.add_argument(
        "--val_ratio",
        type=float,
        default=0.1,
        help="Fraction of genes for validation (default: 0.1)"
    )
    split_parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )

    load_parser = subparsers.add_parser(
        "load",
        help="Load and inspect an existing corpus"
    )
    load_parser.add_argument(
        "--corpus",
        type=str,
        required=True,
        help="Path to JSONL corpus to load"
    )

    args = parser.parse_args()

    if args.command == "demo":
        run_demo(args.dataset_id, args.output)

    elif args.command == "split":
        records = load_corpus(args.corpus)
        train_records, val_records, test_records, split_manifest = gene_level_split(
            records,
            train_ratio=args.train_ratio,
            val_ratio=args.val_ratio,
            seed=args.seed,
        )
        save_splits(train_records, val_records, test_records, split_manifest, args.output_dir)
        print(f"Splits saved to {args.output_dir}")
        print(f"Train: {len(train_records)} records")
        print(f"Val:   {len(val_records)} records")
        print(f"Test:  {len(test_records)} records")

    elif args.command == "load":
        records = load_corpus(args.corpus)
        print(f"Loaded {len(records)} training records")
        print(f"Example record keys: {list(records[0].keys())}")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()

