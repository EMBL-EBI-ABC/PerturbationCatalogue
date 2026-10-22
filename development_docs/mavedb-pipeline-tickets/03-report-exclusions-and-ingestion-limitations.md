# 03: Report exclusions and ingestion limitations

**What to build:** Ingestion covers the planned score-set population while explaining exclusions and retaining observations whose score or HGVS representation is limited.

**Blocked by:** 02: Ingest one eligible score set.

**Status:** ready-for-agent

- [ ] Score sets with zero targets, multiple targets or a nonhuman target are absent from eligible outputs and appear in a separate report with full `score_urn` and reason.
- [ ] A missing required `score` column or duplicate observation identity stops the affected ingestion stage; an invalid or missing score in one row retains that row with an explicit QC outcome.
- [ ] Complex, malformed and unsupported HGVS remain in source observations with their original text, affected coordinate type and reason; scalar fields remain null when they cannot represent the variant.
- [ ] Eligible dataset totals reconcile with retained observations, while excluded score sets are reported outside record-level denominators.
