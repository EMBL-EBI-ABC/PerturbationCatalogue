"""
Evaluate pairwise DEA predictions: one record per (perturbed_gene,
affected_gene) pair, output format "gene X is {up/down}regulated
(log2fc = {value})".

Two things measured per record:
  1. Direction accuracy -- did the model correctly predict up vs down?
  2. log2fc closeness -- how far off was the model's predicted magnitude
     from the true value, for records where direction was predicted
     correctly (a wrong-direction prediction makes the magnitude
     comparison meaningless, so it's excluded from that part).

Usage (as a library, called from evaluate.py):

    from evaluate_pairwise_dea import extract_pairwise_prediction, evaluate_pairwise_dea
"""
import re
import statistics


def extract_pairwise_prediction(text):
    """
    Parse a generated pairwise DEA output string into (direction, log2fc).
    Returns (None, None) if the text doesn't match the expected pattern.
    """
    match = re.search(
        r"is\s+(upregulated|downregulated)\s*\(log2fc\s*=\s*([+-]?\d+\.?\d*)\)",
        text,
        re.IGNORECASE,
    )
    if not match:
        return None, None
    direction = match.group(1).lower()
    try:
        log2fc = float(match.group(2))
    except ValueError:
        return direction, None
    return direction, log2fc


def evaluate_pairwise_dea(predictions, ground_truth_records):
    """
    predictions : list of dict -- must have 'gene' (perturbed_gene value used
        as the record key upstream) is not sufficient alone since pairwise
        records are keyed by (perturbed_gene, affected_gene); predictions
        must instead carry 'perturbed_gene', 'affected_gene', and
        'predicted_text'.
    ground_truth_records : list of dict -- the pairwise test records, with
        metadata containing perturbed_gene, affected_gene, direction, log2fc.

    Returns (metrics dict, results list).
    """
    gt = {
        (r["metadata"]["perturbed_gene"], r["metadata"]["affected_gene"]): {
            "direction": r["metadata"]["direction"],
            "log2fc": r["metadata"]["log2fc"],
        }
        for r in ground_truth_records
        if r["metadata"].get("modality") == "scPerturb-seq_pairwise"
    }

    results = []
    n_total = 0
    n_unparseable = 0
    n_direction_correct = 0
    log2fc_abs_errors = []

    for pred in predictions:
        key = (pred.get("perturbed_gene"), pred.get("affected_gene"))
        if key not in gt:
            continue

        true_direction = gt[key]["direction"]
        true_log2fc = gt[key]["log2fc"]

        pred_direction, pred_log2fc = extract_pairwise_prediction(
            pred.get("predicted_text", "")
        )

        n_total += 1
        if pred_direction is None:
            n_unparseable += 1
            results.append({
                "perturbed_gene": key[0],
                "affected_gene": key[1],
                "true_direction": true_direction,
                "predicted_direction": None,
                "true_log2fc": true_log2fc,
                "predicted_log2fc": None,
                "direction_correct": False,
                "log2fc_abs_error": None,
                "unparseable": True,
            })
            continue

        direction_correct = pred_direction == true_direction
        if direction_correct:
            n_direction_correct += 1

        log2fc_error = None
        if pred_log2fc is not None and direction_correct:
            log2fc_error = abs(pred_log2fc - true_log2fc)
            log2fc_abs_errors.append(log2fc_error)

        results.append({
            "perturbed_gene": key[0],
            "affected_gene": key[1],
            "true_direction": true_direction,
            "predicted_direction": pred_direction,
            "true_log2fc": true_log2fc,
            "predicted_log2fc": pred_log2fc,
            "direction_correct": direction_correct,
            "log2fc_abs_error": log2fc_error,
            "unparseable": False,
        })

    direction_accuracy = n_direction_correct / n_total if n_total > 0 else 0.0
    mean_log2fc_error = (
        statistics.mean(log2fc_abs_errors) if log2fc_abs_errors else None
    )
    median_log2fc_error = (
        statistics.median(log2fc_abs_errors) if log2fc_abs_errors else None
    )

    metrics = {
        "n_evaluated": n_total,
        "direction_accuracy": round(direction_accuracy, 4),
        "n_unparseable": n_unparseable,
        "unparseable_pct": round(100 * n_unparseable / n_total, 1) if n_total > 0 else 0,
        "mean_log2fc_abs_error": round(mean_log2fc_error, 4) if mean_log2fc_error is not None else None,
        "median_log2fc_abs_error": round(median_log2fc_error, 4) if median_log2fc_error is not None else None,
        "n_log2fc_comparisons": len(log2fc_abs_errors),
    }
    return metrics, results