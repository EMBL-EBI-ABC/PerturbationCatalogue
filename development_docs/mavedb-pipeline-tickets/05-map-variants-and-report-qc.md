# 05: Map variants and report mapping QC

**What to build:** A mapped score set retains each original observation and score while reporting usable, ambiguous, mismatched and unmapped variant mappings with meaningful dataset-level rates.

**Blocked by:** 03: Report exclusions and ingestion limitations; 04: Map target candidates to UniProt.

**Status:** ready-for-agent

- [ ] Candidate variant mappings are separate from source scores and retain original accession, HGVS, source position, mapped position, reference identity and validation evidence.
- [ ] Offsets are applied exactly once; source and UniProt reference mismatches, ambiguous mappings, unmapped records and checks not performed have explicit, independent outcomes and reasons.
- [ ] Complex or unsupported variants remain as observations with null unsupported mapped fields; mapping issues do not silently exclude their curated score set.
- [ ] Dataset QC reports counts and denominators for parsing, mapping and reference checks, and a fixture proves ambiguous candidates do not inflate score counts.
