# 04: Map target candidates to UniProt

**What to build:** An analyst can run `MAP` from existing parsed files and immutable UniProt references to see candidate target mappings and source-sequence validation outcomes.

**Blocked by:** 02: Ingest one eligible score set.

**Status:** ready-for-agent

- [ ] `MAP` accepts explicit local or GCS parsed inputs, UniProt references and output locations without rerunning ingestion.
- [ ] Each target candidate records its source identity, UniProt accession or isoform, alignment or offset evidence, and a mapping status; ambiguous candidates remain visible.
- [ ] Source reference symbols are checked against the source target sequence before any mapped UniProt check, with a separate result for a check not performed.
- [ ] Required mapping logic runs from repository-owned command-line code and under Nextflow, with required attribution and licences preserved.
