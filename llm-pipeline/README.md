# LLM Pipeline — Perturbation-Aware Language Model

GSoC 2026 project: Building a Perturbation-Aware LLM for Multimodal In Silico Perturbation Modelling.

## Overview

This pipeline connects to the EMBL-EBI Perturbation Catalogue REST API and converts perturbation experiment data into instruction-tuning training records for fine-tuning a biomedical language model (BioMedLM).

Three data modalities are supported:
- CRISPR screens — gene fitness effects
- scPerturb-seq DEA — differential expression responses
- scPerturb-seq GSEA — pathway-level enrichment responses

## Setup

```bash
pip install requests pandas numpy scipy
```

## How to get training data

```python
from catalogue_api import fetch_and_process_crispr
from catalogue_api import fetch_and_process_perturb_seq
from catalogue_api import fetch_and_process_perturb_seq_gsea

# CRISPR screen
records, df = fetch_and_process_crispr(
    dataset_id="biogrid_5",
    output_path="output/crispr_biogrid5.jsonl"
)

# scPerturb-seq DEA
records, df = fetch_and_process_perturb_seq(
    dataset_id="orion_2025_hct116",
    output_path="output/perturb_seq_dea.jsonl"
)

# scPerturb-seq GSEA
gsea_records, gsea_df = fetch_and_process_perturb_seq_gsea(
    dataset_id="orion_2025_hct116",
    gene_names=df["gene"].tolist(),
    output_path="output/perturb_seq_gsea.jsonl"
)
```

## How to run scripts

```bash
# Run CRISPR demo
python catalogue_api.py --demo

# Run on specific dataset
python catalogue_api.py --dataset_id biogrid_5 --output output/records.jsonl
```

## Training record format

Each record follows instruction-tuning format:

```json
{
  "instruction": "What is the fitness effect of knocking out gene X?",
  "input": "Gene: X. Cell line: Y. Condition: Z.",
  "output": "Gene X is essential for survival...",
  "metadata": {}
}
```

## Verified datasets

| Dataset | Modality | Records | Cell line |
|---|---|---|---|
| biogrid_5 | CRISPR | 4,256 | K562 (CML) |
| biogrid_2373 | CRISPR | 1,649 | SK-N-DZ (neuroblastoma) |
| orion_2025_hct116 | scPerturb-seq DEA | 427 | HCT116 (colon carcinoma) |

## Project structure

- `catalogue_api.py` — Perturbation Catalogue API client and full pipelines
- `preprocess_crispr.py` — CRISPR screen preprocessing from local MAGeCK files
- `preprocess_scrna.py` — scPerturb-seq preprocessing from local h5ad files
- `benchmark.py` — Evaluation framework with gene-level splits
- `finetune.py` — LoRA fine-tuning scaffold (requires Codon cluster)
