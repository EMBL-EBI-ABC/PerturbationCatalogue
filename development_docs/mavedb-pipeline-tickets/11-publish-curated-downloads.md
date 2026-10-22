# 11: Publish curated downloads

**What to build:** Portal users can download complete or filtered MAVE data only for score sets that are currently eligible in the curated projection.

**Blocked by:** 08: Prepare curated public projections.

**Status:** ready-for-agent

- [ ] Complete dataset exports and filtered downloads draw from the same curated population as portal queries and preserve full dataset identity.
- [ ] Download columns expose prepared variant fields, source identifiers, original HGVS and applicable QC consistently with the serving model.
- [ ] Requests for an uncurated score set cannot receive a current export or a signed link to a previously generated export.
- [ ] A fixture compares downloadable rows with the curated projection, including a score set whose source observations have unsupported mappings.
