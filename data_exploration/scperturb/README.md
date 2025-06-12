# scPerturb curated data differential expression

## adamson_2016_pilot_curated.h5ad
* Size: 117M
* Total perturbed genes: 7
* Number of cells per perturbed gene: 500ish
* Number of cell types: 1, lymphoblast

## adamson_2016_upr_epistasis_curated.h5ad
* Size: 479M
* Total perturbed genes: 15
* Number of cells per perturbed gene: 8 very very low (multiple perturbations?), most 1500ish
* Number of cell types: 1, lymphoblast

## adamson_2016_upr_perturb_seq_curated.h5ad
* Size: 1.8 GB
* Total perturbed genes: 90
* Number of cells per perturbed gene: 250-750ish, some very high, also a cluster of very low
* Number of cell types: 1, lymphoblast

## datlinger_2017_curated.h5ad
* Size: 132M
* Total perturbed genes: 32
* Number of cells per perturbed gene: 50-250ish
* Number of cell types: 1, T cell

# Data generation

Once the individual *.csv files have been produced by the Jupyter notebook, unite them into one:

```bash
# 1. Create the header with 'study_id' added
(head -n 1 $(ls *.csv | head -n 1) | sed 's/^/study_id,/' > /tmp/perturb-seq.csv)

# 2. Process all files, skipping their headers and adding basename
for f in *.csv; do bn=$(basename "$f" .csv); tail -n +2 "$f" | sed "s/^/$bn,/"; done >> /tmp/perturb-seq.csv
```