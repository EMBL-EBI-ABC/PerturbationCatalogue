# Development data inventory deck

The eight-slide PDF summarizes the development inventory across CRISPR screen, MAVE and Perturb-seq, including the 20 reprocessed single-condition Perturb-seq datasets and the exported result schemas. The Beamer source is handoff.tex; it has no external image assets.

The modality dataset counts are the current development Elasticsearch inventory (1,197 CRISPR screen, 795 MAVE, 30 Perturb-seq) and are rounded in the PDF. Perturb-seq DEA/GSEA row totals and example rows are from the development BigQuery snapshot queried on 2026-10-13; evidence.sql contains the corresponding read-only SELECTs. The exported fields are defined in dwh/release/release.py and dwh/release/release_worker.py; API mappings and download naming are in be/data_query.py.

To reproduce the PDF from this directory, run:

    pdflatex -interaction=nonstopmode -halt-on-error handoff.tex
