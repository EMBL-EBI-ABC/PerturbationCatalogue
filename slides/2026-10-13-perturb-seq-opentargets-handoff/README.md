# Perturb-Seq handoff deck

The PDF is a six-slide summary of the handoff interface for the 20 single-condition Perturb-Seq datasets. The Beamer source is handoff.tex; it has no external image assets.

The coverage counts and example rows are a read-only snapshot from dev BigQuery queried on 2026-10-13. evidence.sql contains the corresponding SELECT queries. The file and column contract is defined in dwh/release/release.py, dwh/release/release_worker.py, and the download naming in be/data_query.py.

To reproduce the PDF from this directory, run:

    pdflatex -interaction=nonstopmode -halt-on-error handoff.tex
