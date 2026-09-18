#!/bin/bash
#SBATCH --job-name=biomedlm_finetune
#SBATCH --time=08:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --gres=gpu:a100:1
#SBATCH --output=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/logs/finetune_%j.out
#SBATCH --error=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/logs/finetune_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=esranurersan@gmail.com

source /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/venv/bin/activate
cd /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/PerturbationCatalogue/llm-pipeline

pip install torch transformers peft trl accelerate --quiet

python3 finetune.py \
    --splits_dir /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/splits/ \
    --output_dir /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/runs/exp_001/ \
    --model_name stanford-crfm/BioMedLM \
    --epochs 3 \
    --batch_size 4 \
    --grad_accum 4 \
    --lr 2e-4

echo "Fine-tuning complete at $(date)"
