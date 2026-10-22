# 13: Compose preparation and full workflows

**What to build:** Operators can run a composed preparation workflow or a complete first-release workflow from explicit inputs, while each component remains independently runnable.

**Blocked by:** 06: Load internal BigQuery products; 07: Import reviewed curation; 12: Withdraw a public score set.

**Status:** ready-for-agent

- [ ] Preparation composes ingestion, mapping and internal loading; the full workflow additionally accepts an already available final human-reviewed curation export and publishes the curated population.
- [ ] Each standalone workflow accepts existing local or GCS prerequisites and clear output locations; missing prerequisites fail with a message naming the required input.
- [ ] A curation-only correction can run import and publication without repeating mapping or analysis, and `MAP` can run without downloading or parsing again.
- [ ] A small end-to-end run checks the same dataset identities and counts through preparation, curation and publication.
