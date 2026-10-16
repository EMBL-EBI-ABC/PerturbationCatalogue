# Production release checklist

This checklist covers the 2026.10 release from the verified development
deployment to production. It is a review document, not authorization to run a
production command. Stop at any failed gate or unresolved decision below.

The current release candidate uses the shared bucket
`gs://perturbation-catalogue-release` under `2026.10/`. The prefix contains
complete artifact sets for 2,022 datasets (1,197 CRISPR, 795 MAVE and 30
Perturb-seq), 6,126 objects in total. Development BE/FE checks passed against
that prefix. The unprefixed objects were verified unchanged and remain the
current production source until the cutover is approved. Recheck this inventory
before release; do not rely on these counts if the candidate changes.

## 1. Approve and freeze the release candidate

- [ ] Agree on the production release identifier (`2026.10` for this
  candidate), source commit and release date.
- [ ] Merge the approved stacked code changes through `dev` to `main` using the
  normal PR process. Confirm the merged BE supports `RELEASE_VERSION_PREFIX`
  and dataset-level Perturb-seq GSEA; confirm the FE includes the matching
  dataset GSEA display and download controls.
- [ ] Confirm all main-branch Cloud Build checks and the BE/FE production
  deployment builds succeeded. Record the commit, build IDs, service revisions
  and rollback revisions.
- [ ] Freeze data-changing work for the snapshot window. Compare the current
  development and production dataset inventories and agree that replacing the
  production PostgreSQL contents with the development snapshot will not discard
  production-only or newer data.
- [ ] Save the exact dataset/artifact inventory and object generations,
  sizes and checksums for the candidate. Confirm every dataset has its
  modality's expected artifact set:
  - CRISPR and MAVE: metadata JSON, CSV.GZ and Parquet.
  - Perturb-seq: metadata JSON, DEA CSV.GZ and Parquet, and GSEA CSV.GZ and
    Parquet.
- [ ] Treat known schema-valid, zero-row artifacts as present when their
  dataset result is genuinely empty. For the current candidate these include
  MAVE `urn:mavedb:00000105-a-1` and GSEA for `frangieh_2021_raw`,
  `gasperini_2019_atscale`, `jiang_2025`, `song_2025_jurkat_hiv` and
  `zhu_2025_pseudobulk`.

## 2. Prepare the production PostgreSQL database

- [ ] Source the intended environment and verify the project, instance,
  database version, region, tier, storage size/use, private network, current
  Cloud Run connection settings, and backup/PITR policy. Record the current
  production instance and its Cloud Run revision before changing anything.
- [ ] Choose the production tier and storage from measured workload and the
  current snapshot size. The source database is about 1.48 TB and its disk is
  about 2,002 GB; the old 500 GB example is obsolete. The latest operator
  report says production uses 4 vCPU/16 GiB; recheck it before release. On
  2026-10-16, development was resized from 16/64 to 4/16; after restart, SQL
  reconnected and local BE health/search/GSEA and FE checks all returned HTTP
  200, with the 2,002 GB disk unchanged. This smoke test does not establish
  production concurrency or performance. Confirm storage headroom, auto-resize,
  HA/SLA, backups, PITR, maintenance window and deletion protection with the
  owner.
- [ ] Create and wait for a manual backup of the development instance. Record
  its ID, source project/instance, database version, disk size and completion
  time. Keep it through the rollback window.
- [ ] Create a *new* production Cloud SQL instance with the approved version,
  tier, storage, network and data-protection settings. Do not restore over the
  current production instance. Confirm cross-project backup/restore access
  before the change window. The restore target must be `RUNNABLE` and use the
  same PostgreSQL version as the backup.
- [ ] Restore the development backup to the new production instance. Cloud SQL
  restore can cross projects; its resulting disk is at least the larger of the
  configured target disk and the backup disk. Instance flags are not restored,
  so compare and explicitly set the required production flags after restore.
  [Cloud SQL restore behavior](https://docs.cloud.google.com/sql/docs/postgres/backup-recovery/restore).
- [ ] Verify the restored databases, schemas, indexes, extensions, roles,
  BE service-user credentials, expected row counts and 20 reprocessed dataset
  IDs against the release receipt. Check SQL connectivity from the intended
  production VPC path and confirm disk/CPU/connection headroom. Do not point
  the live BE at this instance until these checks pass.

## 3. Refresh and verify production Elasticsearch

- [ ] Resolve the BigQuery source for the production metadata refresh before
  running the projector. The existing docs do not say whether production
  Elasticsearch should read production BigQuery or the verified dev snapshot.
  Whichever source is approved must contain the release's full metadata and
  `perturb_seq_reprocessed=true` for the 20 reprocessed datasets. If it does
  not, stop and agree on a dev-to-prod metadata handoff; do not run the full
  DWH pipeline in production.
- [ ] Record the current `dataset-summary`, `target-summary` and
  `landing-page-summary` alias targets and confirm the previous indexes will
  remain available for rollback.
- [ ] With the approved BigQuery source and production Elasticsearch
  credentials, run only the documented BigQuery-to-Elasticsearch projector
  (`dwh/bq_to_elastic/bq_to_es_projector.py`). Confirm all three fresh indexes
  completed without item errors before accepting the alias switch. The
  projector moves aliases only after all three indexes load, then prunes
  indexes beyond the live plus two previous versions; verify rollback targets
  before running it.
- [ ] Verify the aliases resolve to the new indexes; dataset, target and
  landing-page counts match the approved inventory; all 20 selected datasets
  have the reprocessed marker; and representative searches work. Record the
  new and prior index names.

## 4. Prepare the versioned artifact serving path

- [ ] Recheck that the shared bucket's `2026.10/` inventory is complete and
  matches the frozen manifest. No new artifact-generation build is needed if
  the inventory still matches. Do not modify the unprefixed objects.
- [ ] Confirm the production BE revision includes the versioned-path support
  and its runtime identity can read the bucket and sign URLs. Set
  `RELEASE_BUCKET=perturbation-catalogue-release` and
  `RELEASE_VERSION_PREFIX=2026.10` on the new production BE revision. An empty
  prefix means the legacy unprefixed path; do not leave it empty for this
  release.
- [ ] Keep the versioned prefix immutable during cutover. Record the exact
  object-generation snapshot so any later correction can be reviewed and
  scoped.

## 5. Deploy, cut over and smoke-test

- [ ] Verify the production BE and FE Cloud Run service names, region, URLs,
  service identities, VPC egress and current traffic/revisions. The documented
  main services are `perturbation-catalogue-be` and
  `perturbation-catalogue-dash`; confirm these against Cloud Run before use.
- [ ] Deploy the approved main commit to a tagged or otherwise non-serving
  revision where available. Configure the BE to use the restored database,
  approved production Elasticsearch endpoint, shared bucket and `2026.10`
  prefix. Keep the prior production revision and database intact. If the
  deployment path immediately shifts all traffic, use the agreed maintenance
  window instead of assuming a canary.
- [ ] Test the candidate revision before serving normal traffic:
  - API health and dataset discovery for all three modalities.
  - Dataset metadata and representative DEA/CRISPR/MAVE result queries.
  - Perturb-seq DEA and GSEA result queries, including the reprocessed badge
    for representative Adamson, Norman and Zhu datasets.
  - Full-download redirects for metadata, CSV.GZ and Parquet, plus GSEA CSV.GZ
    and Parquet. Confirm destinations are under `/2026.10/`; open representative
    downloads and verify non-empty files. Do not paste signed URL query strings
    into logs or chat.
  - FE landing/search pages and representative dataset pages, including the
    GSEA table and both GSEA download controls.
- [ ] Compare production API totals and result availability with the approved
  release baseline using targeted per-dataset checks; avoid an unbounded full
  result-table count on the 1.48 TB database. Check application errors, latency,
  database CPU, disk, connections and Elasticsearch health.
- [ ] After owner approval, move production traffic to the tested BE revision;
  verify the production FE still targets that BE and repeat the smoke checks
  against the public URLs. Record the exact serving revisions and traffic
  allocation.

## 6. Rollback and closeout

- [ ] If any gate fails, route Cloud Run traffic back to the recorded previous
  revision. Its previous database and empty `RELEASE_VERSION_PREFIX` continue
  to serve the existing unprefixed artifacts. Cloud Run supports routing all
  traffic back to an existing revision:
  [Cloud Run rollback](https://docs.cloud.google.com/run/docs/rollouts-rollbacks-traffic-migration).
- [ ] If Elasticsearch is the failure, restore the recorded aliases to the
  prior indexes. Do not delete either database, the backup, prior indexes, or
  release objects during rollback.
- [ ] Monitor errors, latency, database saturation and download failures for
  the agreed observation window. Record the release commit, builds, database
  instance, aliases, artifact prefix, checks and any incidents.
- [ ] Publish the GitHub release only after production checks pass. Retain the
  old production SQL instance and migration backup until the owner approves
  cleanup separately.

## Decisions required before execution

- [ ] Confirm which BigQuery project is the approved source for the production
  Elasticsearch refresh and how its 20 reprocessed markers are supplied.
- [ ] Confirm that the development snapshot may replace all current production
  PostgreSQL data, or specify how production-only data must be preserved.
- [ ] Approve production database tier, disk headroom, HA/SLA, backup/PITR and
  rollback observation window. The existing 4/16 setting and a dev smoke test
  alone do not settle this.
- [ ] Confirm the production traffic cutover window and the release/tag naming
  convention.
