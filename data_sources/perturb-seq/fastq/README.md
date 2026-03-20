# Perturb-Seq data reprocessing

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

Here, SAMN40972597 corresponds to Pertur`bation Catalogue dataset accession nadig_2025_jurkat, and SAMN40972598 corresponds to nadig_2025_hepg2.

## Downloading data for the dataset

1. Log in to cluster
2. Set up cluster secrets from `cluster.sh`
3. `cd $HPS_PATH/PerturbationCatalogue/data_sources/perturb_seq/fastq`
4. Start data download, example:
```bash
srun --mem 2G --time 7-00:00:00 \
  python3 ena_download.py \
  --sample-id SAMN40972597 \
  --out-dir $HPS_PATH/perturb_seq_fastq
```
