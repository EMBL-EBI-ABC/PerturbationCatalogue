# 12: Withdraw a public score set

**What to build:** When a formerly public score set loses `is_curated=true`, republication removes its public presence everywhere while retaining internal scientific products.

**Blocked by:** 09: Serve curated MAVE data; 10: Serve curated search and summaries; 11: Publish curated downloads.

**Status:** ready-for-agent

- [ ] A reviewed curation update that clears or falsifies `is_curated` causes publication to reconcile removals from Postgres rows, search documents and summary counts.
- [ ] Existing complete exports are withdrawn or blocked, and direct dataset searches, filtered downloads and signed-link requests no longer expose the score set.
- [ ] Republishing the same eligible population is idempotent; reapproving the score set restores only its current curated projection.
- [ ] An end-to-end fixture verifies withdrawal and restoration across all public paths while internal observations, mappings and QC remain intact.
