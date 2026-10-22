# 07: Import reviewed curation

**What to build:** A curator can import a final human-reviewed CSV or Parquet export into authoritative BigQuery metadata and correct or withdraw a score set without rerunning scientific processing.

**Blocked by:** 02: Ingest one eligible score set.

**Status:** ready-for-agent

- [ ] `IMPORT_CURATION` validates required fields, full score/experiment/experiment-set identifiers, ontology label-ID pairs and duplicate `score_urn` values before changing BigQuery.
- [ ] Curators can use parsed source metadata from GCS with the full score, experiment and experiment-set URNs intact when preparing the reviewed export.
- [ ] Only a declared final human-reviewed export is accepted; `Yes` becomes boolean `is_curated=true`, while blank or false remains ineligible and validation never grants eligibility automatically.
- [ ] Present rows replace curated values, including blanks that clear old values; absent rows remain unchanged, and the report distinguishes inserted, replaced and withdrawn datasets.
- [ ] The import identifies changes to analysis-relevant fields that require reruns, while description and ontology corrections can proceed through import and republication alone.
