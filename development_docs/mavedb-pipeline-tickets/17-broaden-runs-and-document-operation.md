# 17: Broaden runs and document operation

**What to build:** Operators can repeat the first-release pipeline across a broader set of eligible MAVE score sets and understand how to run and troubleshoot its components.

**Blocked by:** 13: Compose preparation and full workflows.

**Status:** ready-for-agent

- [ ] A broader run stays within the single-human-target scope and reconciles eligible, excluded, parsed, loaded, curated and published score-set counts.
- [ ] Resource settings account for substantial variation in observation counts, and failed prerequisites or stages identify the affected score set and next action.
- [ ] Documentation shows independent and composed commands, input ownership, immutable GCS sources, development BigQuery use and the human-reviewed curation boundary.
- [ ] The architecture description reflects parsed metadata in GCS and authoritative curated metadata in BigQuery; deferred versioning, checksums and run manifests are clearly identified as later work.
