import json
import logging
import argparse
import torch
from pathlib import Path
from transformers import AutoModelForCausalLM, AutoTokenizer
from peft import PeftModel
from benchmark import parse_genes_from_output, evaluate as benchmark_evaluate, parse_pathways_from_output, pathway_name_overlap
from rouge_score import rouge_scorer

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
    tokenizer_source = adapter_dir if (adapter_dir and adapter_dir.lower() != "none") else model_name
    log.info(f"Loading tokenizer from {tokenizer_source}...")
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_source)

    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    log.info(f"Loading base model {model_name}...")
    base_model = AutoModelForCausalLM.from_pretrained(
        model_name,
        torch_dtype=torch.bfloat16 if torch.cuda.is_available() else torch.float32,
        device_map="auto" if torch.cuda.is_available() else None,
        trust_remote_code=True,
    )

    if adapter_dir and adapter_dir.lower() != "none":
        log.info(f"Applying LoRA adapter from {adapter_dir}...")
        model = PeftModel.from_pretrained(base_model, adapter_dir)
    else:
        log.info("No adapter — running zero-shot baseline")
        model = base_model
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
        return "unparseable"


def extract_genes_from_output(text, direction="up"):
    """
    Extract gene names from DEA model output text.
    Delegates to benchmark.parse_genes_from_output for consistency.
    """
    return parse_genes_from_output(text, direction)


def compute_rouge_l(predictions, ground_truth_records):
    """
    Compute ROUGE-L between predicted and ground truth outputs.
    Format-agnostic — works for all modalities.

    Parameters
    ----------
    predictions : list of dict — gene, predicted_text
    ground_truth_records : list of dict — training records with output field

    Returns
    -------
    dict of metrics
    """
    scorer = rouge_scorer.RougeScorer(["rougeL"], use_stemmer=True)
    gt = {r["metadata"]["gene"]: r["output"] for r in ground_truth_records}

    scores = []
    for pred in predictions:
        gene = pred["gene"]
        if gene not in gt:
            continue
        score = scorer.score(gt[gene], pred["predicted_text"])
        scores.append(score["rougeL"].fmeasure)

    if not scores:
        return {"rouge_l_mean": 0.0, "rouge_l_median": 0.0, "n_evaluated": 0}

    import statistics
    return {
        "rouge_l_mean": round(sum(scores) / len(scores), 4),
        "rouge_l_median": round(statistics.median(scores), 4),
        "n_evaluated": len(scores),
    }


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
        if "fitness_class" in r["metadata"]
    }

    correct = 0
    unparseable = 0
    total = 0
    results = []

    for pred in predictions:
        gene = pred["gene"]
        if gene not in gt:
            continue

        true_class = gt[gene]
        pred_class = pred["predicted_class"]
        if pred_class == "unparseable":
            unparseable += 1
        is_correct = pred_class == true_class

        correct += int(is_correct)
        total += 1

        results.append(
            {
                "gene": gene,
                "true_class": true_class,
                "predicted_class": pred_class,
                "correct": is_correct,
                "predicted_text": pred["predicted_text"],
                "true_output": pred["true_output"],
            }
        )

    accuracy = correct / total if total > 0 else 0.0
    log.info(f"Evaluated {total} genes, accuracy: {accuracy:.4f}")

    return {
        "n_evaluated": total,
        "accuracy": round(accuracy, 4),
        "correct": correct,
        "unparseable": unparseable,
        "unparseable_pct": round(100 * unparseable / total, 1) if total > 0 else 0,
    }, results


def evaluate_dea(predictions, ground_truth_records, k=10):
    """
    Evaluate scPerturb-seq DEA gene predictions.
    Delegates to benchmark.evaluate for full metrics including
    gene set overlap, direction accuracy, and pathway enrichment.
    """
    return benchmark_evaluate(predictions, ground_truth_records, k=k)



def evaluate_gsea(predictions, ground_truth_records):
    """
    Evaluate GSEA pathway predictions using pathway name overlap.
    Uses parse_pathways_from_output() to extract pathway names from
    generated text, then compares to ground truth from metadata.
    """
    gt = {
        r["metadata"]["gene"]: {
            "activated": r["metadata"].get("activated_pathways", []),
            "suppressed": r["metadata"].get("suppressed_pathways", []),
        }
        for r in ground_truth_records
        if r["metadata"].get("modality") == "scPerturb-seq_GSEA"
    }

    results = []
    total = 0
    unparseable = 0
    f1_scores = []

    for pred in predictions:
        gene = pred["gene"]
        if gene not in gt:
            continue
        predicted_text = pred.get("predicted_text", "")
        true_activated = gt[gene]["activated"]
        true_suppressed = gt[gene]["suppressed"]
        parsed = parse_pathways_from_output(predicted_text)
        if parsed["unparseable"]:
            unparseable += 1
        act_score = pathway_name_overlap(parsed["activated"], true_activated)
        sup_score = pathway_name_overlap(parsed["suppressed"], true_suppressed)
        mean_f1 = (act_score["f1"] + sup_score["f1"]) / 2
        f1_scores.append(mean_f1)
        total += 1
        results.append({
            "gene": gene,
            "modality": "scPerturb-seq_GSEA",
            "predicted_text": predicted_text,
            "predicted_activated": parsed["activated"],
            "predicted_suppressed": parsed["suppressed"],
            "true_activated": true_activated,
            "true_suppressed": true_suppressed,
            "f1_activated": act_score["f1"],
            "f1_suppressed": sup_score["f1"],
            "f1_mean": mean_f1,
            "recall_activated": act_score["recall"],
            "recall_suppressed": sup_score["recall"],
            "unparseable": parsed["unparseable"],
        })

    mean_f1_all = sum(f1_scores) / len(f1_scores) if f1_scores else 0.0
    log.info(f"Evaluated {total} GSEA genes, mean F1: {mean_f1_all:.4f}")

    return {
        "n_evaluated": total,
        "mean_f1": round(mean_f1_all, 4),
        "unparseable": unparseable,
        "unparseable_pct": round(100 * unparseable / total, 1) if total > 0 else 0,
    }, results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate fine-tuned LoRA adapter on test split"
    )
    parser.add_argument(
        "--adapter_dir",
        type=str,
        required=False,
        default=None,
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

    # Split predictions by modality
    crispr_preds = [p for p in predictions if "CRISPR" in p["modality"]]
    dea_preds    = [p for p in predictions if p["modality"] == "scPerturb-seq"]
    gsea_preds   = [p for p in predictions if p["modality"] == "scPerturb-seq_GSEA"]

    crispr_records = [r for r in test_records if "CRISPR" in r["metadata"].get("modality","")]
    dea_records    = [r for r in test_records if r["metadata"].get("modality","") == "scPerturb-seq"]
    gsea_records   = [r for r in test_records if r["metadata"].get("modality","") == "scPerturb-seq_GSEA"]

    rouge_metrics = compute_rouge_l(predictions, test_records)

    print("\n" + "=" * 60)
    print("EVALUATION RESULTS")
    print("=" * 60)
    print(f"Model:        {args.model_name}")
    print(f"Adapter:      {args.adapter_dir}")
    print(f"Total test records: {len(predictions)}")
    print(f"ROUGE-L mean: {rouge_metrics['rouge_l_mean']:.4f}")
    print(f"ROUGE-L median: {rouge_metrics['rouge_l_median']:.4f}")

    results = []
    metrics = {"n_evaluated": len(predictions)}

    if crispr_preds:
        crispr_metrics, crispr_results = evaluate_crispr(crispr_preds, crispr_records)
        print(f"\nCRISPR ({len(crispr_preds)} records):")
        print(f"  Accuracy: {crispr_metrics['accuracy']:.4f} ({crispr_metrics['correct']}/{crispr_metrics['n_evaluated']})")
        print(f"  Unparseable: {crispr_metrics.get('unparseable', 0)} ({crispr_metrics.get('unparseable_pct', 0):.1f}%)")
        results.extend(crispr_results)

    if dea_preds:
        dea_metrics, dea_results = evaluate_dea(dea_preds, dea_records)
        print(f"\nDEA ({len(dea_preds)} records):")
        print(f"  Mean overlap@k: {dea_metrics.get('mean_overlap_at_k', dea_metrics.get('mean_overlap_both', 0)):.4f}")
        dea_text = {p["gene"]: p["predicted_text"] for p in dea_preds}
        dea_list = dea_results.to_dict("records") if hasattr(dea_results, "to_dict") else dea_results
        for r in dea_list:
            r["predicted_text"] = dea_text.get(r.get("gene"), "")
            r["modality"] = "scPerturb-seq"
        results.extend(dea_list)

    if gsea_preds:
        gsea_metrics, gsea_results = evaluate_gsea(gsea_preds, gsea_records)
        print(f"\nGSEA ({len(gsea_preds)} records):")
        print(f"  Mean F1: {gsea_metrics.get('mean_f1', 0):.4f}")
        print(f"  Unparseable: {gsea_metrics.get('unparseable', 0)} ({gsea_metrics.get('unparseable_pct', 0):.1f}%)")
        results.extend(gsea_results)

    print("=" * 60)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            for r in results:
                f.write(json.dumps(r) + "\n")
        log.info(f"Saved per-gene results to {args.output}")
        # Save summary metrics
        summary_path = Path(args.output).parent / "eval_summary.json"
        import json as json_module
        summary = {
            "model": args.model_name,
            "adapter": args.adapter_dir,
            "modality": first_modality,
            "n_evaluated": metrics["n_evaluated"],
            "task_metrics": metrics,
            "rouge_l": rouge_metrics,
        }
        with open(summary_path, "w") as f:
            f.write(json_module.dumps(summary, indent=2))
        log.info(f"Saved eval summary to {summary_path}")



if __name__ == "__main__":
    main()