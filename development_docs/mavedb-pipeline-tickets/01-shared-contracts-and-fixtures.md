# 01: Shared contracts and fixtures

**What to build:** A checkable shared contract for MAVE source identities, observations, scores, mappings, coordinates and QC, backed by representative fixtures that future pipeline slices can use independently.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] The contracts define complete row keys and grains, including full `score_urn` as `dataset_id`, original variant accessions, and one primary `score` per observation.
- [ ] Coordinate conventions distinguish source target positions, UniProt positions and original nucleotide reference types; parsing, mapping and sequence validation outcomes have separate states and reasons.
- [ ] Fixtures cover substitutions, nonzero offsets, complex HGVS, ambiguity, reference mismatches, unmapped targets, extra score columns, and zero-target, multi-target and nonhuman score sets.
- [ ] A focused check validates fixture identities, score cardinality and the expected QC outcomes before later stages use them.
