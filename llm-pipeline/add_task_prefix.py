"""
Add explicit modality prefix tokens to instruction text.

Single-variable experiment (Q2): tests whether giving the model an explicit
task signal ([CRISPR] / [DEA] / [GSEA]) reduces multi-task interference,
compared to exp_002 where the model infers modality only from instruction
phrasing.

This is a pure transform on an already-generated corpus -- no API calls, no
GPU needed. Input is exp_002's own corpus (full_corpus_oversampled.jsonl) so
this branches directly from the same baseline exp_003 did, not from exp_003
itself. Only the instruction field changes; input, output, and metadata are
untouched.

--dropout_rate controls what fraction of records get the tag stripped back
out after being loaded (default 0.0 = every record tagged, for the initial
test of whether tagging helps at all). If results are good and this gets
adopted going forward, rerun with e.g. --dropout_rate 0.3 so the resulting
model isn't dependent on the tag always being present at inference time.

Usage:

    python3 add_task_prefix.py \
        --input data/full_corpus_oversampled.jsonl \
        --output data/full_corpus_q2_prefixed.jsonl \
        --dropout_rate 0.0
"""
import json
import argparse
import random
from pathlib import Path

MODALITY_PREFIX = {
    "CRISPR_screen": "[CRISPR] ",
    "scPerturb-seq": "[DEA] ",
    "scPerturb-seq_GSEA": "[GSEA] ",
}


def add_prefix(records, dropout_rate=0.0, seed=42):
    """
    Prepend a modality tag to each record's instruction.

    dropout_rate : float in [0, 1]
        Fraction of records that keep NO tag (randomly chosen), so the
        model doesn't become fully dependent on the tag being present.
        0.0 = every record tagged. 0.3 = 70% tagged, 30% untagged.
    """
    rng = random.Random(seed)
    tagged = 0
    untagged_by_dropout = 0
    skipped = 0

    for r in records:
        modality = r["metadata"].get("modality", "")
        prefix = MODALITY_PREFIX.get(modality)
        if prefix is None:
            skipped += 1
            continue

        if dropout_rate > 0 and rng.random() < dropout_rate:
            untagged_by_dropout += 1
            continue

        r["instruction"] = prefix + r["instruction"]
        tagged += 1

    return records, tagged, untagged_by_dropout, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--dropout_rate", type=float, default=0.0,
        help="Fraction of records left untagged (default 0.0 = tag everything)"
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    records = [json.loads(l) for l in open(args.input) if l.strip()]
    print(f"Loaded {len(records)} records from {args.input}")

    records, tagged, untagged, skipped = add_prefix(
        records, dropout_rate=args.dropout_rate, seed=args.seed
    )
    print(f"Tagged: {tagged}, left untagged by dropout: {untagged}, "
          f"skipped (unrecognized modality): {skipped}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")
    print(f"Saved {len(records)} records to {args.output}")

    for modality, prefix in MODALITY_PREFIX.items():
        sample = next((r for r in records if r["metadata"].get("modality") == modality), None)
        if sample:
            print(f"\n{modality}: {sample['instruction'][:90]}")


if __name__ == "__main__":
    main()