import json
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

from catalogue_api_new import (
    query_crispr_screen,
    query_perturb_seq,
    get_dataset_metadata,
    identify_primary_score,
    pivot_gene_records,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


def fitness_class_to_text(gene, fitness_class, lfc, cell_line, condition, perturbation_type="knockout"):
    """
    Convert CRISPR fitness classification to natural language.
    """
    descriptions = {
        "essential": (
            f"{gene} is essential for survival of {cell_line} "
            f"under {condition} conditions (LFC: {lfc:.2f}). "
            f"CRISPR-mediated {perturbation_type} causes significant cell depletion, indicating "
            f"this gene is required for cell fitness under these "
            f"experimental conditions."
        ),
        "anti_essential": (
            f"{gene} acts as a fitness suppressor in {cell_line} "
            f"under {condition} conditions (LFC: {lfc:.2f}). "
            f"CRISPR-mediated {perturbation_type} causes cell enrichment, meaning cells without "
            f"this gene grow faster under these experimental conditions."
        ),
        "neutral": (
            f"{gene} shows no significant fitness effect in {cell_line} "
            f"under {condition} conditions (LFC: {lfc:.2f}). "
            f"CRISPR-mediated {perturbation_type} does not substantially alter cell survival or "
            f"proliferation under these specific conditions."
        ),
    }
    return descriptions.get(fitness_class, "Fitness effect unknown.")


def normalise_within_dataset(df, score_col="effect_score"):
    """Z-score normalise effect scores within each dataset."""
    df = df.copy()
    df[f"{score_col}_zscore"] = df.groupby("dataset_id")[score_col].transform(
        lambda x: stats.zscore(x, nan_policy="omit")
    )
    log.info(
        f"Z-score normalised {score_col} within {df['dataset_id'].nunique()} datasets"
    )
    return df


def classify_from_catalogue(df, zscore_col="effect_score_zscore", zscore_threshold=1.5, use_significant_only=False):
    """Classify genes into essential / anti_essential / neutral."""
    if use_significant_only:
        conditions = [
            df["significant"],
            pd.Series([False] * len(df), index=df.index),
        ]
    else:
        conditions = [
            df["significant"] & (df[zscore_col] < -zscore_threshold),
            df["significant"] & (df[zscore_col] > zscore_threshold),
        ]
    choices = ["essential", "anti_essential"]
    df["fitness_class"] = np.select(conditions, choices, default="neutral")
    counts = df["fitness_class"].value_counts()
    log.info(f"Classification: {counts.to_dict()}")
    return df


def catalogue_records_to_training(df, dataset_id, modality="CRISPR_screen", include_condition=True, perturbation_type="knockout"):
    """Convert harmonised CRISPR Catalogue records into training record format."""
    records = []

    for _, row in df.iterrows():
        gene = row["gene"]
        fitness_class = row.get("fitness_class", "neutral")
        effect_score = row.get("effect_score", 0.0)
        zscore = row.get("effect_score_zscore", np.nan)
        cell_line = row.get("cell_line", "unknown")
        disease = row.get("disease", "unknown")

        display_score = zscore if not np.isnan(zscore) else effect_score
        condition = disease if disease != "unknown" else "standard growth"
        output_text = fitness_class_to_text(
            gene, fitness_class, display_score, cell_line, condition, perturbation_type
        )

        record = {
            "instruction": (
                f"What is the fitness effect of CRISPR-mediated {perturbation_type} of gene {gene} "
                f"in {cell_line} under {condition}? "
                f"Describe the phenotype and its biological interpretation."
            ),
            "input": (
                f"Gene: {gene}. "
                f"Cell line: {cell_line}."
                + (f" Condition: {condition}." if include_condition and condition != "standard growth" else "")
            ),
            "output": output_text,
            "metadata": {
                "gene": gene,
                "dataset_id": dataset_id,
                "cell_line": cell_line,
                "disease": disease,
                "effect_score": float(effect_score),
                "effect_score_zscore": float(zscore) if not np.isnan(zscore) else None,
                "fitness_class": fitness_class,
                "modality": modality,
                "source": "perturbation_catalogue_api",
            },
        }
        records.append(record)

    log.info(f"Built {len(records)} training records from {dataset_id}")
    return records


def fetch_and_process_crispr(dataset_id, output_path=None, max_records=5000, include_condition=True):
    """Fetch CRISPR screen data from Catalogue API and save as training records."""
    log.info(f"Fetching CRISPR data for {dataset_id}...")

    raw = query_crispr_screen(dataset_id=dataset_id, max_records=max_records)
    if not raw:
        log.warning(f"No data returned for {dataset_id}")
        return [], pd.DataFrame()

    metadata = get_dataset_metadata(dataset_id)
    log.info(f"Dataset metadata: {metadata}")

    df = pivot_gene_records(raw)
    df["dataset_id"] = dataset_id
    for key, value in metadata.items():
        df[key] = value

    df = normalise_within_dataset(df)
    use_sig_only = (dataset_id == "biogrid_2373")
    df = classify_from_catalogue(df, use_significant_only=use_sig_only)
    records = catalogue_records_to_training(df, dataset_id, include_condition=include_condition, perturbation_type=metadata.get("library_perturbation_type", "knockout") or "knockout")

    if output_path and records:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            for record in records:
                f.write(json.dumps(record) + "\n")
        log.info(f"Saved {len(records)} records to {output_path}")

    return records, df


def fetch_and_process_perturb_seq(dataset_id, output_path=None, max_records=5000, include_condition=True):
    """Fetch scPerturb-seq DEA data from Catalogue API and save as training records."""
    log.info(f"Fetching DEA data for {dataset_id}...")

    raw = query_perturb_seq(dataset_id=dataset_id, max_records=max_records)
    if not raw:
        log.warning(f"No data returned for {dataset_id}")
        return [], pd.DataFrame()

    metadata = get_dataset_metadata(dataset_id)
    log.info(f"Dataset metadata: {metadata}")

    perturbation_effects = {}
    for r in raw:
        # FIX (Aug 18): API renamed "gene_name" to "gene_symbol".
        raw_gene_name = r.get("perturbation", {}).get("gene_symbol", "unknown")
        if "|" in str(raw_gene_name):
            parts = raw_gene_name.split("|")
            if "control_nontargeting" in parts[1].lower():
                perturbed_gene = parts[0]
            else:
                continue
        else:
            perturbed_gene = raw_gene_name
        effect = r.get("effect", {})

        record_cell_type = effect.get("cell_type") or "unknown"
        group_key = (perturbed_gene, record_cell_type)
        if group_key not in perturbation_effects:
            perturbation_effects[group_key] = {
                "up_genes": [],
                "down_genes": [],
                "cell_type": record_cell_type,
                "n_total": r.get("perturbation", {}).get("n_total", 0),
                "n_up": r.get("perturbation", {}).get("n_up", 0),
                "n_down": r.get("perturbation", {}).get("n_down", 0),
            }

        direction = effect.get("direction", "")
        log2fc = effect.get("log2fc", 0.0)
        # FIX (Aug 18): same rename applies to the affected-gene field.
        affected_gene = effect.get("gene_symbol", "unknown")
        padj = effect.get("padj", 1.0)

        if padj < 0.05:
            if not affected_gene.startswith("ENSG"):
                if direction == "increased":
                    perturbation_effects[group_key]["up_genes"].append(
                        (affected_gene, round(log2fc, 3))
                    )
                elif direction == "decreased":
                    perturbation_effects[group_key]["down_genes"].append(
                        (affected_gene, round(log2fc, 3))
                    )

    log.info(f"Grouped into {len(perturbation_effects)} unique perturbations")

    records = []
    rows = []
    cell_line = metadata.get("cell_line", "unknown")
    disease = metadata.get("disease", "unknown")
    condition = disease if disease != "unknown" else "standard growth"
    perturbation_type = metadata.get("library_perturbation_type", "knockout") or "knockout"

    for (perturbed_gene, record_cell_type), effects in perturbation_effects.items():
        cell_line_display = cell_line if cell_line != "unknown" else record_cell_type
        up_genes = sorted(effects["up_genes"], key=lambda x: abs(x[1]), reverse=True)[:10]
        down_genes = sorted(effects["down_genes"], key=lambda x: abs(x[1]), reverse=True)[:10]

        up_str = (
            ", ".join([f"{g} ({fc:+.2f})" for g, fc in up_genes])
            if up_genes else "none detected"
        )
        down_str = (
            ", ".join([f"{g} ({fc:+.2f})" for g, fc in down_genes])
            if down_genes else "none detected"
        )

        n_shown_up = len(up_genes)
        n_shown_down = len(down_genes)

        output_text = (
            f"CRISPR-mediated {perturbation_type} of {perturbed_gene} in {cell_line_display} causes "
            f"upregulation of: {up_str}; "
            f"and downregulation of: {down_str}. "
            f"Top differentially expressed genes shown: "
            f"{n_shown_up} upregulated, {n_shown_down} downregulated "
            f"(filtered by padj < 0.05, ranked by absolute log2fc)."
        )

        record = {
            "instruction": (
                f"Predict the transcriptional response to CRISPR-mediated {perturbation_type} "
                f"of gene {perturbed_gene} in {cell_line_display} under {condition}. "
                f"Describe the key upregulated and downregulated genes."
            ),
            "input": (
                f"Gene: {perturbed_gene}. "
                f"Cell line: {cell_line_display}."
                + (f" Condition: {condition}." if include_condition and condition != "standard growth" else "")
            ),
            "output": output_text,
            "metadata": {
                "gene": perturbed_gene,
                "dataset_id": dataset_id,
                "cell_line": cell_line_display,
                "disease": disease,
                "top_up_genes": [g for g, _ in up_genes],
                "top_down_genes": [g for g, _ in down_genes],
                "n_shown_up": n_shown_up,
                "n_shown_down": n_shown_down,
                "n_total_api": effects["n_total"],
                "n_up_api": effects["n_up"],
                "n_down_api": effects["n_down"],
                "modality": "scPerturb-seq",
                "source": "perturbation_catalogue_api",
            },
        }
        records.append(record)
        rows.append({
            "gene": perturbed_gene,
            "cell_line": cell_line_display,
            "n_up": effects["n_up"],
            "n_down": effects["n_down"],
            "n_total": effects["n_total"],
        })

    if output_path and records:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            for record in records:
                f.write(json.dumps(record) + "\n")
        log.info(f"Saved {len(records)} DEA records to {output_path}")

    df = pd.DataFrame(rows)
    if output_path and rows:
        genes_path = Path(output_path).parent / "genes.txt"
        with open(genes_path, "w") as f_genes:
            for row in rows:
                f_genes.write(row["gene"] + "\n")
        log.info(f"Saved gene list to {genes_path}")

    log.info(f"Built {len(records)} training records from {dataset_id}")
    return records, df


def fetch_and_process_perturb_seq_gsea(
    dataset_id, gene_names, output_path=None, fdr_threshold=0.05, top_n_pathways=5, include_condition=True
):
    """Fetch scPerturb-seq GSEA data from Catalogue API and save as training records."""
    import requests
    import time

    BASE_URL = "https://perturbation-catalogue-be-328296435987.europe-west2.run.app"

    log.info(f"Fetching GSEA data for {len(gene_names)} genes in {dataset_id}...")

    metadata = get_dataset_metadata(dataset_id)
    cell_line = metadata.get("cell_line", "unknown")
    disease = metadata.get("disease", "unknown")
    condition = disease if disease != "unknown" else "standard growth"
    perturbation_type = metadata.get("library_perturbation_type", "knockout") or "knockout"

    records = []
    rows = []
    failed = []

    for i, gene in enumerate(gene_names):
        endpoint = f"{BASE_URL}/v1/perturb-seq-gsea"
        # FIX (Aug 18): API now requires "perturbation_gene_name" — the old
        # param name "perturbed_gene_name" returns a 422 error.
        params = {"dataset_id": dataset_id, "perturbation_gene_name": gene}

        try:
            response = requests.get(endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
        except Exception as e:
            log.warning(f"GSEA query failed for {gene}: {e}")
            failed.append(gene)
            continue

        if not data:
            continue

        item = data[0] if isinstance(data, list) else data
        effects = item.get("effects", [])
        if not effects:
            continue

        activated = []
        suppressed = []

        for pathway in effects:
            fdr = pathway.get("fdr", 1.0)
            nes = pathway.get("nes", 0.0)
            term = pathway.get("term", "unknown")

            if fdr < fdr_threshold:
                clean_term = term.replace("HALLMARK_", "").replace("_", " ").title()
                if nes > 0:
                    activated.append((clean_term, round(nes, 3), fdr))
                else:
                    suppressed.append((clean_term, round(nes, 3), fdr))

        activated = sorted(activated, key=lambda x: abs(x[1]), reverse=True)[:top_n_pathways]
        suppressed = sorted(suppressed, key=lambda x: abs(x[1]), reverse=True)[:top_n_pathways]

        if not activated and not suppressed:
            continue

        act_str = (
            ", ".join([f"{term} (NES: {nes:+.2f})" for term, nes, _ in activated])
            if activated else "none detected"
        )
        sup_str = (
            ", ".join([f"{term} (NES: {nes:+.2f})" for term, nes, _ in suppressed])
            if suppressed else "none detected"
        )

        output_text = (
            f"CRISPR-mediated {perturbation_type} of {gene} in {cell_line} activates pathways: {act_str}. "
            f"Suppressed pathways: {sup_str}."
        )

        record = {
            "instruction": (
                f"What biological pathways are affected by CRISPR-mediated {perturbation_type} "
                f"of gene {gene} in {cell_line} under {condition}? "
                f"Describe the activated and suppressed pathways."
            ),
            "input": (
                f"Gene: {gene}. "
                f"Cell line: {cell_line}."
                + (f" Condition: {condition}." if include_condition and condition != "standard growth" else "")
            ),
            "output": output_text,
            "metadata": {
                "gene": gene,
                "dataset_id": dataset_id,
                "cell_line": cell_line,
                "disease": disease,
                "activated_pathways": [t for t, _, _ in activated],
                "suppressed_pathways": [t for t, _, _ in suppressed],
                "n_activated": len(activated),
                "n_suppressed": len(suppressed),
                "modality": "scPerturb-seq_GSEA",
                "source": "perturbation_catalogue_api",
            },
        }
        records.append(record)
        rows.append({
            "gene": gene,
            "cell_line": cell_line,
            "n_activated": len(activated),
            "n_suppressed": len(suppressed),
        })

        if (i + 1) % 10 == 0:
            log.info(f"Processed {i + 1}/{len(gene_names)} genes")

        time.sleep(0.2)

    log.info(f"Built {len(records)} GSEA training records from {dataset_id}")
    if failed:
        log.warning(f"Failed for {len(failed)} genes: {failed[:5]}...")

    if output_path and records:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            for record in records:
                f.write(json.dumps(record) + "\n")
        log.info(f"Saved {len(records)} GSEA records to {output_path}")

    df = pd.DataFrame(rows) if rows else pd.DataFrame()
    return records, df


def main():
    parser = argparse.ArgumentParser(
        description="Fetch perturbation data from Catalogue API and save as training records"
    )
    parser.add_argument("--modality", choices=["crispr", "perturb_seq", "gsea"], required=True, help="Data modality to fetch")
    parser.add_argument("--dataset_id", type=str, required=True, help="Catalogue dataset ID")
    parser.add_argument("--output", type=str, required=True, help="Output JSONL file path")
    parser.add_argument("--max_records", type=int, default=5000, help="Maximum records to fetch (CRISPR and DEA only)")
    parser.add_argument("--no_condition", action="store_true", help="Exclude disease/condition from input fields.")
    parser.add_argument("--genes_file", type=str, default=None, help="Text file with gene names one per line (GSEA only)")
    args = parser.parse_args()

    if args.modality == "crispr":
        records, _ = fetch_and_process_crispr(
            dataset_id=args.dataset_id,
            output_path=args.output,
            max_records=args.max_records,
            include_condition=not args.no_condition,
        )
        print(f"Saved {len(records)} CRISPR training records to {args.output}")

    elif args.modality == "perturb_seq":
        records, df = fetch_and_process_perturb_seq(
            dataset_id=args.dataset_id,
            output_path=args.output,
            max_records=args.max_records,
            include_condition=not args.no_condition,
        )
        print(f"Saved {len(records)} DEA training records to {args.output}")

    elif args.modality == "gsea":
        if not args.genes_file:
            print("Error: --genes_file required for GSEA modality")
            return
        gene_names = Path(args.genes_file).read_text().strip().splitlines()
        records, _ = fetch_and_process_perturb_seq_gsea(
            dataset_id=args.dataset_id,
            gene_names=gene_names,
            output_path=args.output,
            include_condition=not args.no_condition,
        )
        print(f"Saved {len(records)} GSEA training records to {args.output}")


if __name__ == "__main__":
    main()