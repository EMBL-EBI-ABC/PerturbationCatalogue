# 08: Prepare curated public projections

**What to build:** One queryable curated-only MAVE projection supplies the portal and release with consistent metadata, scores and prepared variant fields.

**Blocked by:** 06: Load internal BigQuery products; 07: Import reviewed curation.

**Status:** ready-for-agent

- [ ] Public projections join internal products to authoritative metadata by full `score_urn` and include only rows with boolean `is_curated=true`.
- [ ] Projection keys and joins preserve one catalogue dataset per score set, source observation identity and primary-score cardinality despite multiple mapping candidates.
- [ ] Prepared position, reference and alternate fields are derived upstream once; raw HGVS and applicable QC remain available when a scalar mapping is unsupported.
- [ ] A fixture with curated and uncurated score sets proves both public data and public metadata projections contain exactly the eligible population.
