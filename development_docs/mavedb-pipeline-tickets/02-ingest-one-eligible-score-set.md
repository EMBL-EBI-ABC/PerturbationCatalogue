# 02: Ingest one eligible score set

**What to build:** A curator or analyst can run `INGEST` on a representative single-human-target score set and obtain reusable source metadata, observations and primary scores in GCS.

**Blocked by:** 01: Shared contracts and fixtures.

**Status:** ready-for-agent

- [ ] The workflow accepts either a Zenodo dump URL or an existing dump, plus explicit output locations, and preserves the original dump and extracted source files while outputs depend on them.
- [ ] An eligible score set produces GCS parsed metadata, observations and scores with full score, experiment and experiment-set URNs, original accessions and original HGVS strings.
- [ ] Only the required `score` column enters parsed measurements; other source columns remain available in the immutable originals.
- [ ] The repository-owned parser is callable from a normal command line, and the same implementation runs under `INGEST` on a small fixture.
