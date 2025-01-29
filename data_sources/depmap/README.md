# DepMap

## Background

[DepMap](https://depmap.org/portal/data_page/?tab=overview) gene dependency data is hosted on FigShare and is licensed under CC Attribution 4.0.

Current release being used is [24Q4](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248).

We currently use three files from them:
* `CRISPRGeneDependency.csv`: matrix where rows are cancer screens and columns are genes. Values contain probabilities which indicate how likely is a given cancer line depend on a given gene.
* `CRISPRInferredCommonEssentials.csv`: list of genes which DepMap analysis has deeemed to be essential to all cancer cell lines. These are filtered out during the processing.
* `Model.csv`: cancer cell line metadata.

Before processing, set the following environment secrets:
* LAKE_BUCKET
* WAREHOUSE_BUCKET
* ELASTIC_ENDPOINT

## Ingest gene dependency data

```bash
# Mirror.
mirror_depmap () { wget -q -O - "$1" | gsutil cp - "gs://${LAKE_BUCKET}/depmap/$2"; }
mirror_depmap https://ndownloader.figshare.com/files/51064631 CRISPRGeneDependency.csv
mirror_depmap https://ndownloader.figshare.com/files/51064916 CRISPRInferredCommonEssentials.csv
mirror_depmap https://ndownloader.figshare.com/files/51065828 Model.csv
```
