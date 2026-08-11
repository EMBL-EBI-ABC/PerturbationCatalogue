"""Shared release-artifact definitions and BigQuery query builders."""

DATASETS = {
    "crispr": {
        "data": "crispr.data",
        "fields": [
            ("perturbed_target_ensg", "Perturbed Target ENSG"),
            ("perturbed_target_name", "Perturbed Target Name"),
            ("score_name", "Score Name"),
            ("score_value", "Score Value"),
            ("significant", "Significant"),
            ("significance_criteria", "Significance Criteria"),
        ],
    },
    "perturb-seq": {
        "data": "perturb_seq.pertpy_dea",
        "fields": [
            ("perturbed_target_ensg", "Perturbed Target ENSG"),
            ("perturbed_target_name", "Perturbed Target Name"),
            ("effect_gene_ensg", "Effect Gene ENSG"),
            ("effect_gene_name", "Effect Gene Name"),
            ("log2foldchange", "Log2FC"),
            ("padj", "Padj"),
            ("score_name", "Score Name"),
            ("score_value", "Score Value"),
            ("cell_type", "Cell Type"),
        ],
    },
    "mave": {
        "data": "mavedb.data",
        "fields": [
            ("perturbed_target_ensg", "Perturbed Target ENSG"),
            ("perturbed_target_name", "Perturbed Target Name"),
            ("perturbation_name", "Perturbation Name"),
            ("perturbation_position", "Position"),
            ("perturbation_aa_wt", "AA WT"),
            ("perturbation_aa_change", "AA Change"),
            ("score_name", "Score Name"),
            ("score_value", "Score Value"),
        ],
    },
}

DATASET_METADATA = "dataset_summary"


def table(project, dataset, name):
    return f"`{project}.{name}`" if "." in name else f"`{project}.{dataset}.{name}`"


def data_query(project, dataset, modality):
    source = table(project, dataset, DATASETS[modality]["data"])
    fields = []
    for field, _ in DATASETS[modality]["fields"]:
        if field == "perturbed_target_name":
            fields.append("d.perturbed_target_symbol AS perturbed_target_name")
        elif field == "effect_gene_name":
            fields.append("d.effect_gene_symbol AS effect_gene_name")
        elif field == "perturbation_position":
            fields.append(
                "SAFE_CAST(REGEXP_EXTRACT(d.perturbation_name, "
                "r'p\\.[A-Za-z]+(\\d+)') AS INT64) AS perturbation_position"
            )
        elif field == "perturbation_aa_wt":
            fields.append(
                "REGEXP_EXTRACT(d.perturbation_name, "
                "r'p\\.([A-Za-z]+)\\d+') AS perturbation_aa_wt"
            )
        elif field == "perturbation_aa_change":
            fields.append(
                "REGEXP_EXTRACT(d.perturbation_name, "
                "r'p\\.[A-Za-z]+\\d+([A-Za-z=]+)') AS perturbation_aa_change"
            )
        else:
            fields.append(f"d.{field}")
    return f"SELECT d.dataset_id, {', '.join(fields)} FROM {source} AS d"
