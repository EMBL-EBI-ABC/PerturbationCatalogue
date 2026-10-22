# 10: Serve curated search and summaries

**What to build:** Portal search, target summaries and landing counts discover only curated MAVE score sets from the shared public projection.

**Blocked by:** 08: Prepare curated public projections.

**Status:** ready-for-agent

- [ ] MAVE search documents are built only from curated metadata, keyed by full `score_urn`, and a direct search cannot reveal an uncurated dataset.
- [ ] Dataset, target and landing summaries count the same curated population, with no counts inferred from all internal MAVE rows.
- [ ] Repeated publication updates changed search documents and summaries without duplicate documents or inflated counts.
- [ ] A mixed curated/uncurated fixture demonstrates parity between search results and summary counts.
