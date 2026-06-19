import json
import logging
import argparse
import torch
from pathlib import Path
from transformers import AutoModelForCausalLM, AutoTokenizer
from peft import PeftModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)

PROMPT_TEMPLATE = (
    "### Instruction:\n{instruction}\n\n" "### Input:\n{input}\n\n" "### Response:\n"
)


def load_model(adapter_dir, model_name):
    """
    Load base model and apply saved LoRA adapter.

    Parameters
    ----------
    adapter_dir : str
        Path to saved adapter directory.
    model_name : str
        Base model HuggingFace ID.

    Returns
    -------
    tuple of (model, tokenizer)
    """
    log.info(f"Loading tokenizer from {adapter_dir}...")
    tokenizer = AutoTokenizer.from_pretrained(adapter_dir)

    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    log.info(f"Loading base model {model_name}...")
    base_model = AutoModelForCausalLM.from_pretrained(
        model_name,
        torch_dtype=torch.float32,
        trust_remote_code=True,
    )

    log.info(f"Applying LoRA adapter from {adapter_dir}...")
    model = PeftModel.from_pretrained(base_model, adapter_dir)
    model.eval()

    return model, tokenizer


def generate_response(model, tokenizer, instruction, input_text, max_new_tokens=200):
    """
    Generate model response for a single record.

    Parameters
    ----------
    model : PeftModel
    tokenizer : PreTrainedTokenizer
    instruction : str
    input_text : str
    max_new_tokens : int

    Returns
    -------
    str — generated response text
    """
    prompt = PROMPT_TEMPLATE.format(
        instruction=instruction,
        input=input_text,
    )

    inputs = tokenizer(
        prompt,
        return_tensors="pt",
        truncation=True,
        max_length=512,
    )

    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=False,
            temperature=1.0,
            pad_token_id=tokenizer.eos_token_id,
        )

    generated = tokenizer.decode(
        outputs[0][inputs["input_ids"].shape[1] :], skip_special_tokens=True
    )

    return generated.strip()


def extract_fitness_class(text):
    """
    Extract predicted fitness class from generated text.

    Parameters
    ----------
    text : str

    Returns
    -------
    str — one of: essential, anti_essential, neutral, unknown
    """
    text_lower = text.lower()

    if any(
        phrase in text_lower
        for phrase in [
            "essential for survival",
            "required for cell fitness",
            "causes significant cell depletion",
        ]
    ):
        return "essential"
    elif any(
        phrase in text_lower
        for phrase in [
            "fitness suppressor",
            "cells without this gene grow faster",
            "causes cell enrichment",
        ]
    ):
        return "anti_essential"
    elif any(
        phrase in text_lower
        for phrase in [
            "no significant fitness effect",
            "does not significantly alter",
            "does not substantially alter",
            "does not alter cell survival",
        ]
    ):
        return "neutral"
    else:
        return "unknown"


def extract_genes_from_output(text, direction="up"):
    """
    Extract gene names from DEA model output text.

    Parameters
    ----------
    text : str
    direction : str — up or down

    Returns
    -------
    list of gene name strings
    """
    import re

    if direction == "up":
        pattern = r"upregulation of[:\s]+([^;]+?)(?:\.|$)"
    else:
        pattern = r"downregulation of[:\s]+([^;]+?)(?:\.|$)"

    match = re.search(pattern, text, re.IGNORECASE)
    if not match:
        return []

    section = match.group(1)
    genes = re.findall(r"\b([A-Z][A-Z0-9\-]{1,10})\b", section)
    return genes


def evaluate_crispr(predictions, ground_truth_records):
    """
    Evaluate CRISPR fitness class predictions.

    Parameters
    ----------
    predictions : list of dict
    ground_truth_records : list of dict

    Returns
    -------
    tuple of (metrics dict, results list)
    """
    gt = {
        r["metadata"]["gene"]: r["metadata"]["fitness_class"]
        for r in ground_truth_records
    }

    correct = 0
    total = 0
    results = []

    for pred in predictions:
        gene = pred["gene"]
        if gene not in gt:
            continue

        true_class = gt[gene]
        pred_class = pred["predicted_class"]
        is_correct = pred_class == true_class

        correct += int(is_correct)
        total += 1

        results.append(
            {
                "gene": gene,
                "true_class": true_class,
                "predicted_class": pred_class,
                "correct": is_correct,
            }
        )

    accuracy = correct / total if total > 0 else 0.0
    log.info(f"Evaluated {total} genes, accuracy: {accuracy:.4f}")

    return {
        "n_evaluated": total,
        "accuracy": round(accuracy, 4),
        "correct": correct,
    }, results


def evaluate_dea(predictions, ground_truth_records, k=10):
    """
    Evaluate scPerturb-seq DEA gene predictions.

    Parameters
    ----------
    predictions : list of dict
    ground_truth_records : list of dict
    k : int

    Returns
    -------
    tuple of (metrics dict, results list)
    """
    gt = {r["metadata"]["gene"]: r["metadata"] for r in ground_truth_records}

    results = []

    for pred in predictions:
        gene = pred["gene"]
        if gene not in gt:
            continue

        true_up = gt[gene].get("top_up_genes", [])
        true_down = gt[gene].get("top_down_genes", [])

        pred_up = extract_genes_from_output(pred["predicted_text"], "up")
        pred_down = extract_genes_from_output(pred["predicted_text"], "down")

        pred_up_set = set(pred_up[:k])
        true_up_set = set(true_up[:k])
        pred_down_set = set(pred_down[:k])
        true_down_set = set(true_down[:k])

        overlap_up = len(pred_up_set & true_up_set) / max(len(true_up_set), 1)
        overlap_down = len(pred_down_set & true_down_set) / max(len(true_down_set), 1)

        results.append(
            {
                "gene": gene,
                "overlap_up": round(overlap_up, 4),
                "overlap_down": round(overlap_down, 4),
                "mean_overlap": round((overlap_up + overlap_down) / 2, 4),
            }
        )

    if not results:
        return {"n_evaluated": 0, "mean_overlap": 0.0}, results

    mean_overlap = sum(r["mean_overlap"] for r in results) / len(results)

    return {
        "n_evaluated": len(results),
        "mean_overlap_at_k": round(mean_overlap, 4),
        "k": k,
    }, results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate fine-tuned LoRA adapter on test split"
    )
    parser.add_argument(
        "--adapter_dir",
        type=str,
        required=True,
        help="Path to saved LoRA adapter directory",
    )
    parser.add_argument(
        "--splits_dir",
        type=str,
        required=True,
        help="Path to splits directory containing test.jsonl",
    )
    parser.add_argument(
        "--model_name",
        type=str,
        default="stanford-crfm/BioMedLM",
        help="Base model HuggingFace ID (default: stanford-crfm/BioMedLM)",
    )
    parser.add_argument(
        "--max_new_tokens",
        type=int,
        default=200,
        help="Maximum tokens to generate per record (default: 200)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional path to save per-gene results as JSONL",
    )
    args = parser.parse_args()

    model, tokenizer = load_model(args.adapter_dir, args.model_name)

    test_path = Path(args.splits_dir) / "test.jsonl"
    test_records = []
    with open(test_path) as f:
        for line in f:
            line = line.strip()
            if line:
                test_records.append(json.loads(line))

    log.info(f"Loaded {len(test_records)} test records")

    first_modality = test_records[0]["metadata"].get("modality", "unknown")
    log.info(f"Detected modality: {first_modality}")

    predictions = []
    for i, record in enumerate(test_records):
        gene = record["metadata"]["gene"]
        modality = record["metadata"].get("modality", "unknown")

        log.info(f"Evaluating {i+1}/{len(test_records)}: {gene}")

        generated = generate_response(
            model,
            tokenizer,
            record["instruction"],
            record["input"],
            max_new_tokens=args.max_new_tokens,
        )

        pred_class = extract_fitness_class(generated)

        predictions.append(
            {
                "gene": gene,
                "modality": modality,
                "predicted_text": generated,
                "predicted_class": pred_class,
                "true_output": record["output"],
            }
        )

        log.info(f"  Predicted class: {pred_class}")
        log.info(f"  Generated: {generated[:100]}...")

    if "CRISPR" in first_modality:
        metrics, results = evaluate_crispr(predictions, test_records)
        metric_display = f"Accuracy: {metrics['accuracy']:.4f} ({metrics['correct']}/{metrics['n_evaluated']})"
    else:
        metrics, results = evaluate_dea(predictions, test_records)
        metric_display = f"Mean overlap@k: {metrics.get('mean_overlap_at_k', 0):.4f}"

    print("\n" + "=" * 60)
    print("EVALUATION RESULTS")
    print("=" * 60)
    print(f"Model:        {args.model_name}")
    print(f"Adapter:      {args.adapter_dir}")
    print(f"Modality:     {first_modality}")
    print(f"Test records: {metrics['n_evaluated']}")
    print(f"{metric_display}")
    print("=" * 60)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            for r in results:
                f.write(json.dumps(r) + "\n")
        log.info(f"Saved per-gene results to {args.output}")


if __name__ == "__main__":
    main()
