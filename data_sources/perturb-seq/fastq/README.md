# Perturb-Seq data reprocessing

## SRA toolkit configuration on the cluster (one time only)
Start an interactive job using `sinteractive`. Then:

```bash
cd $HPS_PATH
mkdir -p cache/sra
mkdir -p software
cd software
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
rm sratoolkit.current-ubuntu64.tar.gz
```

In your `~/.bashrc` / `~/.bash_profile`, add `$HPS_PATH/software/sratoolkit.3.4.1-ubuntu64/bin` to your PATH.

Configure cache settings:
1. Run vdb-config -i
2. Navigate to the CACHE tab.
3. Ensure "enable local file-caching" is checked.
4. Change the Location of user-repository to `$HPS_PATH/cache/sra`
5. Press S to save and X to exit.

## Dataset information

List of datasets raw data projects is stored in `datasets.tsv`.

To view dataset information, run `python3 ena_report.py datasets.tsv`

## Downloading example dataset

Information from the script above:

```
PRJNA1100571	2.618	2688
    #1 SAMN40972597 0.829 1792 [cell_line: Jurkat, center_name: SUB14380315, isolate: Jurkat, tissue_type: acute T cell leukemia]
    #2 SAMN40972598 1.789 896 [cell_line: HepG2, center_name: SUB14382917, isolate: HepG2, tissue_type: Hepatocellular carcinoma]
```

Here, SAMN40972597 corresponds to Perturbation Catalogue dataset accession nadig_2025_jurkat, and SAMN40972598 corresponds to nadig_2025_hepg2.

To download data for a dataset:

1. Log in to cluster
2. Set up cluster secrets from `cluster.sh`
3. `cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb-seq/fastq`
4. Start data download, example:
```bash
while read -r SAMPLE_ID; do
  echo "Starting ${SAMPLE_ID}"

  time srun --cpus-per-task 64 --mem-per-cpu 2G --time 7-00:00:00 --unbuffered \
    python3 ena_download.py \
    --sample-id "${SAMPLE_ID}" \
    --out-dir "$HPS_PATH/perturb_seq_fastq" \
    --jobs 64

done < samples.txt
```
