# SRA dataset discovery

`datasets.tsv` lists raw-data projects. Run `python3 ena_report.py datasets.tsv`
to inspect ENA project/sample metadata. A sample sheet for counting is generated
by the dataset-specific scripts under `../pipeline/datasets/`.

SRA retrieval and read extraction are part of the [main pipeline](../pipeline/README.md).
It streams technical and biological reads directly into kb-python without
storing FASTQ files.
