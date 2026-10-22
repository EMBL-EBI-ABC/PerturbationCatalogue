# 06: Load internal BigQuery products

**What to build:** An analyst can load selected parsed and mapped products into development BigQuery tables for querying without publishing the score sets.

**Blocked by:** 03: Report exclusions and ingestion limitations; 05: Map variants and report mapping QC.

**Status:** ready-for-agent

- [ ] `LOAD_DATA` accepts selected existing local or GCS inputs and loads observations, primary scores, target mappings, variant mappings, variant QC and dataset QC at their declared grains.
- [ ] Complete keys preserve full `score_urn` and original observation identity; the primary score has no measurement-name dimension.
- [ ] Repeating a load replaces or merges the selected development scope without duplicate rows, and the load report reconciles source and BigQuery row counts.
- [ ] Internal data remains unavailable through public queries and downloads solely because it was loaded.
