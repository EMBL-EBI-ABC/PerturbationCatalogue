# 09: Serve curated MAVE data

**What to build:** Portal users can query curated MAVE score sets and their observations through Postgres-backed API searches, with source identity and QC intact.

**Blocked by:** 08: Prepare curated public projections.

**Status:** ready-for-agent

- [ ] The MAVE publication entrypoint loads the curated projection into the existing Postgres serving model using the full score-set URN as `dataset_id`.
- [ ] Direct dataset searches and filtered data queries return prepared variant fields, original HGVS and relevant QC without multiplying a source score across mapping candidates.
- [ ] An uncurated score set present in internal BigQuery tables cannot be queried through the MAVE API.
- [ ] A focused integration check exercises one curated and one uncurated score set through the serving query path.
