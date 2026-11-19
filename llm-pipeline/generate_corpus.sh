#!/bin/bash
#SBATCH --job-name=perturbation_corpus
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --output=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/logs/corpus_%j.out
#SBATCH --error=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/logs/corpus_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=esranurersan@gmail.com

source /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/venv/bin/activate
cd /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/PerturbationCatalogue/llm-pipeline

DATA=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/data

echo "Starting corpus generation at $(date)"

# CRISPR datasets
echo "=== CRISPR: biogrid_5 ==="
python3 prepare_training_data.py --modality crispr --dataset_id biogrid_5 --output $DATA/crispr_biogrid5.jsonl

echo "=== CRISPR: biogrid_2373 ==="
python3 prepare_training_data.py --modality crispr --dataset_id biogrid_2373 --output $DATA/crispr_biogrid2373.jsonl

# DEA + GSEA per dataset (run together so genes.txt is fresh for each)
echo "=== DEA + GSEA: nadig_2025_hepg2 ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id nadig_2025_hepg2 --output $DATA/dea_nadig_hepg2.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id nadig_2025_hepg2 --genes_file $DATA/genes.txt --output $DATA/gsea_nadig_hepg2.jsonl

echo "=== DEA + GSEA: nadig_2025_jurkat ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id nadig_2025_jurkat --output $DATA/dea_nadig_jurkat.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id nadig_2025_jurkat --genes_file $DATA/genes.txt --output $DATA/gsea_nadig_jurkat.jsonl

echo "=== DEA + GSEA: replogle_2022_k562_gw_normalized ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id replogle_2022_k562_gw_normalized --output $DATA/dea_replogle_k562_gw.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id replogle_2022_k562_gw_normalized --genes_file $DATA/genes.txt --output $DATA/gsea_replogle_k562_gw.jsonl

echo "=== DEA + GSEA: replogle_2022_k562_essential_normalized ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id replogle_2022_k562_essential_normalized --output $DATA/dea_replogle_k562_ess.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id replogle_2022_k562_essential_normalized --genes_file $DATA/genes.txt --output $DATA/gsea_replogle_k562_ess.jsonl

echo "=== DEA + GSEA: replogle_2022_rpe1_essential_normalized ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id replogle_2022_rpe1_essential_normalized --output $DATA/dea_replogle_rpe1.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id replogle_2022_rpe1_essential_normalized --genes_file $DATA/genes.txt --output $DATA/gsea_replogle_rpe1.jsonl

# DEA only (no GSEA) datasets
echo "=== DEA only: arce_2025 ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id arce_2025 --output $DATA/dea_arce.jsonl

echo "=== DEA only: norman_2019_raw ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id norman_2019_raw --output $DATA/dea_norman.jsonl

# Merge all
echo "=== Merging corpus ==="
cat $DATA/crispr_biogrid5.jsonl $DATA/crispr_biogrid2373.jsonl \
    $DATA/dea_nadig_hepg2.jsonl $DATA/dea_nadig_jurkat.jsonl \
    $DATA/dea_replogle_k562_gw.jsonl $DATA/dea_replogle_k562_ess.jsonl \
    $DATA/dea_replogle_rpe1.jsonl $DATA/dea_arce.jsonl $DATA/dea_norman.jsonl \
    $DATA/gsea_nadig_hepg2.jsonl $DATA/gsea_nadig_jurkat.jsonl \
    $DATA/gsea_replogle_k562_gw.jsonl $DATA/gsea_replogle_k562_ess.jsonl \
    $DATA/gsea_replogle_rpe1.jsonl > $DATA/full_corpus.jsonl

echo "=== Corpus merged ==="
wc -l $DATA/full_corpus.jsonl

# Split
echo "=== Creating gene-level splits ==="
python3 pipeline.py split --corpus $DATA/full_corpus.jsonl --output_dir /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/splits/

echo "Corpus generation complete at $(date)"
ENDOFFILEcat > /workspaces/PerturbationCatalogue/llm-pipeline/generate_corpus.sh << 'ENDOFFILE'
#!/bin/bash
#SBATCH --job-name=perturbation_corpus
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --output=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/logs/corpus_%j.out
#SBATCH --error=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/logs/corpus_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=esranurersan@gmail.com

source /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/venv/bin/activate
cd /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/PerturbationCatalogue/llm-pipeline

DATA=/hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/data

echo "Starting corpus generation at $(date)"

# CRISPR datasets
echo "=== CRISPR: biogrid_5 ==="
python3 prepare_training_data.py --modality crispr --dataset_id biogrid_5 --output $DATA/crispr_biogrid5.jsonl

echo "=== CRISPR: biogrid_2373 ==="
python3 prepare_training_data.py --modality crispr --dataset_id biogrid_2373 --output $DATA/crispr_biogrid2373.jsonl

# DEA + GSEA per dataset (run together so genes.txt is fresh for each)
echo "=== DEA + GSEA: nadig_2025_hepg2 ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id nadig_2025_hepg2 --output $DATA/dea_nadig_hepg2.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id nadig_2025_hepg2 --genes_file $DATA/genes.txt --output $DATA/gsea_nadig_hepg2.jsonl

echo "=== DEA + GSEA: nadig_2025_jurkat ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id nadig_2025_jurkat --output $DATA/dea_nadig_jurkat.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id nadig_2025_jurkat --genes_file $DATA/genes.txt --output $DATA/gsea_nadig_jurkat.jsonl

echo "=== DEA + GSEA: replogle_2022_k562_gw_normalized ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id replogle_2022_k562_gw_normalized --output $DATA/dea_replogle_k562_gw.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id replogle_2022_k562_gw_normalized --genes_file $DATA/genes.txt --output $DATA/gsea_replogle_k562_gw.jsonl

echo "=== DEA + GSEA: replogle_2022_k562_essential_normalized ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id replogle_2022_k562_essential_normalized --output $DATA/dea_replogle_k562_ess.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id replogle_2022_k562_essential_normalized --genes_file $DATA/genes.txt --output $DATA/gsea_replogle_k562_ess.jsonl

echo "=== DEA + GSEA: replogle_2022_rpe1_essential_normalized ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id replogle_2022_rpe1_essential_normalized --output $DATA/dea_replogle_rpe1.jsonl
python3 prepare_training_data.py --modality gsea --dataset_id replogle_2022_rpe1_essential_normalized --genes_file $DATA/genes.txt --output $DATA/gsea_replogle_rpe1.jsonl

# DEA only (no GSEA) datasets
echo "=== DEA only: arce_2025 ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id arce_2025 --output $DATA/dea_arce.jsonl

echo "=== DEA only: norman_2019_raw ==="
python3 prepare_training_data.py --modality perturb_seq --dataset_id norman_2019_raw --output $DATA/dea_norman.jsonl

# Merge all
echo "=== Merging corpus ==="
cat $DATA/crispr_biogrid5.jsonl $DATA/crispr_biogrid2373.jsonl \
    $DATA/dea_nadig_hepg2.jsonl $DATA/dea_nadig_jurkat.jsonl \
    $DATA/dea_replogle_k562_gw.jsonl $DATA/dea_replogle_k562_ess.jsonl \
    $DATA/dea_replogle_rpe1.jsonl $DATA/dea_arce.jsonl $DATA/dea_norman.jsonl \
    $DATA/gsea_nadig_hepg2.jsonl $DATA/gsea_nadig_jurkat.jsonl \
    $DATA/gsea_replogle_k562_gw.jsonl $DATA/gsea_replogle_k562_ess.jsonl \
    $DATA/gsea_replogle_rpe1.jsonl > $DATA/full_corpus.jsonl

echo "=== Corpus merged ==="
wc -l $DATA/full_corpus.jsonl

# Split
echo "=== Creating gene-level splits ==="
python3 pipeline.py split --corpus $DATA/full_corpus.jsonl --output_dir /hps/nobackup/mfreeberg/perturb_seq_fastq/gsoc-2026/splits/

echo "Corpus generation complete at $(date)"
