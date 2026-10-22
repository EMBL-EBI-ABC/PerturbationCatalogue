# 14: Produce basic assay QC

**What to build:** An analyst can assess whether a parsed score set has enough usable source data for later analyses, without requiring canonical protein mapping.

**Blocked by:** 03: Report exclusions and ingestion limitations.

**Status:** ready-for-agent

- [ ] Basic assay QC runs from existing parsed metadata, observations and primary scores and writes a reusable result with reasons for unusable or insufficient data.
- [ ] The report distinguishes retained source observations from observations eligible for each analysis and gives explicit exclusion counts and denominators.
- [ ] Raw scores remain separate from any derived effect or normalization result, and unsupported mapping does not erase a source score.
- [ ] The repository-owned command-line implementation and its `ANALYSE` entrypoint agree on a small fixture.
