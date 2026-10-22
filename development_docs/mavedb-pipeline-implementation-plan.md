# MAVE pipeline implementation plan

Implement a modular Nextflow pipeline under `data_sources/mavedb/` that brings MaveDB parsing, UniProt mapping, quality assessment, meta-analysis, metadata curation import and portal publication into the Perturbation Catalogue (PC) repository.

Build a working pipeline on representative eligible score sets first, then expand coverage and polish the interfaces. The first public release contains curated source scores, mappings and QC; analysis results follow later. Defer the versioning system, checksums and run manifests until the processing steps and outputs have settled. Initial development uses explicit input/output paths, ordinary configuration files and development cloud tables. Keep source dumps and UniProt reference objects immutable while pipeline or public outputs depend on them.

**Agreed requirements**

- Parsed source metadata lives in GCS and supplies the metadata curation workflow.
- Curated metadata lives in BigQuery and determines which entries are publicly visible.
- A catalogue dataset is identified by one full MaveDB `score_urn`, with its `experiment_urn` and `experiment_set_urn`. Include only score sets with exactly one human target (`taxId=9606`). Record the score URNs and reasons for excluding zero-target, multi-target and nonhuman score sets.
- Only datasets with `is_curated=true` in authoritative BigQuery metadata appear in portal results, search summaries, counts or public downloads. Import this state only from a final human-reviewed curation export, normalizing `Yes` to true; no separate approval field is required.
- Use the required `score` column as the only source measurement for now. Other score-file columns need not enter parsed outputs; retain the original source files in GCS.
- Parsing, mapping and meta-analysis implementations become part of the PC repository. The running pipeline must use code maintained here.
- Preserve original variant accessions and HGVS strings, including complex or ambiguous perturbations.
- Unsupported parsing, ambiguous mapping, reference mismatch and unmapped records remain explicit QC outcomes. Retain variant details and aggregate relevant dataset metrics.
- Individual modules and useful combinations of modules can run independently from existing files in GCS or local storage.

**Storage and ownership**

| Product | Primary location | Purpose |
| --- | --- | --- |
| Original dump and extracted source files | GCS | Preserve immutable source inputs, including columns omitted from parsed outputs |
| Parsed source metadata | GCS, preferably Parquet with raw source JSON retained | Source material for curation and analysis context |
| Parsed observations and primary scores | GCS Parquet; internal BigQuery tables | Scientific input files and queryable observations using the `score` column |
| UniProt references and mapping outputs | GCS; internal BigQuery mapping tables | Immutable reference inputs and explicit mapping results |
| Variant and dataset QC | GCS; internal BigQuery QC tables | Diagnostics, coverage reports and analysis eligibility |
| Meta-analysis outputs | GCS; internal BigQuery result tables | Derived effects, uncertainty and model diagnostics |
| Curated metadata | BigQuery | Authoritative curated context and public eligibility |
| Portal data and search documents | Postgres and Elasticsearch, populated from BigQuery | Serving projections for curated entries |
| Public dataset downloads | Existing release storage, populated from curated BigQuery projections | The same curated population as the portal |

The initial implementation does not need a separate parsed-metadata table in BigQuery. Its authoritative parsed metadata files remain in GCS. Internal observation tables may carry identifiers needed for joins without duplicating the full source metadata. No input checksums or run manifest are required in this iteration.

**Proposed data flow**

```mermaid
flowchart TD
    Z[Zenodo dump] --> I[Select eligible score sets and parse]
    I --> PD[Parsed observations and scores in GCS]
    I --> PM[Parsed metadata in GCS]
    I --> EX[Excluded score sets and reasons in GCS]
    PD --> M[UniProt mapping and sequence validation]
    M --> MG[Mappings and QC in GCS]
    PD --> L[Load internal BigQuery tables]
    MG --> L
    PD --> A[Meta-analysis on Slurm]
    MG --> A
    PM --> A
    A --> AG[Analysis results and diagnostics in GCS]
    AG --> L
    PM --> C[Metadata curation]
    C --> CF[Final human-reviewed curation CSV or Parquet]
    CF --> CI[Validate and import curation]
    CI --> CM[Curated metadata in BigQuery]
    L --> P[Prepare public projections using curated eligibility]
    CM --> P
    P --> PG[Postgres]
    P --> ES[Elasticsearch]
    P --> DL[Public downloads]
```

The first public release uses the ingestion, mapping, internal loading, curation import and publication branches; the `ANALYSE` branch is added later. Analysis modules should declare the inputs they actually use. Basic assay QC can operate on parsed data; analyses requiring canonical protein positions also consume mapping outputs. If analysis needs curated fields such as score interpretation, supply an explicit extract of those fields from BigQuery as another input.

**Nextflow modules and independent execution**

| Proposed workflow | Main inputs | Outputs |
| --- | --- | --- |
| `INGEST` | Zenodo file URL or existing dump, output directory | Eligible parsed metadata, observations, primary scores, ingestion QC and score-set exclusion report |
| `MAP` | Parsed files, UniProt reference files, mapping configuration | Target mappings, variant mappings and mapping QC |
| `LOAD_DATA` | Selected parsed, mapped, QC or analysis files | Internal BigQuery tables and load report |
| `ANALYSE` | Declared parsed/mapped inputs and analysis configuration | QC, normalized effects, optional combined effects and fit diagnostics |
| `IMPORT_CURATION` | Final human-reviewed curated CSV/Parquet and source metadata identifiers | Validated curated metadata in BigQuery and import report |
| `PUBLISH` | Internal BigQuery products and curated metadata | Curated public projections, Postgres, Elasticsearch and downloads |

Provide a composed preparation workflow for ingestion, mapping and internal loading. Provide a full automated workflow that can consume an already available final human-reviewed curation file. Human curation remains an activity between runs; the automated workflow imports its completed output, never an intermediate LLM-generated file.

Each standalone workflow accepts explicit input paths. Running `MAP` against existing parsed files must not require downloading or parsing again. Importing a metadata correction must not require repeating mapping or fitting. Missing prerequisites should produce a clear error naming the required input.

Separate Nextflow orchestration from the implementation scripts. Keep transformations callable from a normal command line so they can be developed and checked on small files before they are wired into Nextflow.

**Repository layout and migration**

Suggested layout; exact filenames can evolve during implementation:

```text
data_sources/mavedb/
  README.md
  pipeline/
    main.nf
    nextflow.config
    conf/
    modules/
    subworkflows/
  parsing/
  mapping/
  analysis/
    R/
    stan/
  schemas/
  containers/
  tests/
    fixtures/
```

Bring the relevant extraction, identifier handling and residue validation logic from `protvar-mave-import` into `mapping/` and `parsing/`. Bring the required R/Stan methods, analysis entrypoints and supporting resources from `ot_perturbation_catalog_mave_experiments` into `analysis/`. Preserve required attribution and licences when transferring code.

Create a single shared parsing interface. Adapt mapping and analysis to consume its outputs so each workflow uses the same accession handling, HGVS representation and position conventions. Remove dependencies on sibling repository paths and scheduler-specific batch wrappers. Retain only the analysis methods and resources required by the selected scientific products; migrate additional experiments when they become relevant.

**Data contracts to establish first**

| Product | One row represents | Essential fields |
| --- | --- | --- |
| Source score-set metadata | One eligible score set and its single human target | Full `score_urn`, `experiment_urn`, `experiment_set_urn`; study references; target and assay context |
| Source observations | One original MaveDB variant record from an eligible score set | Original accession, `score_urn`, variant number, original nucleotide/protein/splice HGVS and parsing outcomes |
| Source scores | One primary `score` value for an observation | Original accession, `score_urn` and raw `score` value |
| Target mappings | One candidate target-to-reference mapping | Score-set URN, source target identity, UniProt accession/isoform, offset or alignment details and mapping status |
| Variant mappings | One candidate mapping for an observation and target | Accession, source target, reference identity, source position, mapped position, reference/alternate residue and validation outcomes |
| Derived effects | One declared analysis result at its stated grain | Source or derived identity, method, value, uncertainty where available and contribution links |
| Curated metadata | One catalogue dataset, keyed by full `score_urn` | `score_urn`, `experiment_urn`, `experiment_set_urn`, curated experimental context and explicit public eligibility |

Keep full score-set URNs and original accessions throughout processing. Use the full `score_urn` as the catalogue `dataset_id`; retain source identity when creating compatibility fields such as `sample_id`. Do not split a score set into multiple catalogue datasets. Reject it from this pipeline's eligible set if it does not have exactly one human target, and report its `score_urn` and exclusion reason. This general rule covers the two multi-target score sets `urn:mavedb:00001237-a-1` and `urn:mavedb:00001237-a-2` previously excluded by the analysis repository.

Require the source `score` column and parse it as the sole measurement. Do not infer measurements from column position, name patterns or numeric type. Ignore other score-file columns in parsed outputs for now; they remain available in the original dump. A missing `score` column is an ingestion contract failure; missing or invalid values in individual rows produce explicit QC outcomes while retaining those rows. If a later analysis needs another source column, extend the contract and reparse the immutable source files.

Use explicit names for target-sequence positions and UniProt positions. Preserve the original nucleotide coordinate type and its reference sequence. Canonical protein mapping must account for offsets exactly once. Genomic variant mapping is a separate extension requiring an assembly, transcript and allele contract.

**Parsing, mapping and QC behavior**

Keep parsing, mapping and sequence validation outcomes separate so a record can carry several issues. For example, an observation can have valid nucleotide parsing while its protein HGVS is unsupported, or have a candidate UniProt mapping whose reference residue does not match.

| Required outcome | Variant-level record | Dataset-level summary |
| --- | --- | --- |
| Unsupported parsing | Original HGVS, affected coordinate type and unsupported construct/reason | Counts and rates by coordinate type and variant class |
| Ambiguous mapping | Candidate references/positions and reason no unique mapping was chosen | Ambiguous observations and affected score sets/targets |
| Reference mismatch | Expected and observed reference symbol, reference sequence and checked position | Mismatch counts and rates among observations actually checked |
| Unmapped | Mapping attempted and reason no usable mapping was found | Counts and rates among eligible observations |

Also distinguish malformed input, missing required information, successful checks, checks that were not performed and checks that do not apply. Missing protein HGVS does not automatically make a nucleotide observation invalid.

Validate reference symbols against the source target sequence before checking mapped symbols against UniProt. Define coordinate indexing and permitted residue/nucleotide symbols in the schema and parser behavior. Never replace a source position with a mapped position in place.

Retain complex observations even when scalar position/ref/alt fields cannot represent them. Initially those fields can remain null with an explicit unsupported status. Candidate mappings belong in a separate relation so ambiguous mappings cannot silently multiply source score rows.

Variant QC should identify the observation, check, outcome and reason, with relevant reference details. Dataset QC should include total observations, parseable observations, mapping outcomes, reference checks and metric denominators. Counts for different checks may overlap and must not be presented as mutually exclusive categories unless defined that way. Report score-set eligibility exclusions separately from record-level QC so source selection does not distort dataset denominators.

Mapping or parsing issues do not automatically remove a curated dataset or its observations. Preserve raw HGVS, expose the applicable QC status, and leave unsupported mapped fields null. Analyses that require a particular representation should select eligible observations explicitly and report excluded counts and reasons.

Structural contract failures such as missing required columns or duplicate observation identities should stop the affected stage. Expected record-level limitations should produce QC outcomes and retained source records.

**Meta-analysis execution**

Implement basic score-set QC first, then per-assay normalization and model diagnostics. Define the intended cross-study combined-effect product separately: its grouping rules, comparability assumptions and contributing observations must be explicit before public use.

Run analysis scripts inside Nextflow tasks submitted through the Slurm executor. Use containers or an explicitly defined reproducible environment for R, Stan and supporting tools. Start with score-set-level tasks; group-level aggregation follows the applicable assay tasks. Resource settings should account for large differences in observation counts.

Stage selected GCS inputs to cluster storage for reuse across tasks. Write analysis products and diagnostics to GCS, then let `LOAD_DATA` load the selected products into BigQuery. Analysis tasks should not write directly to Postgres or Elasticsearch.

Keep raw scores, normalized effects, aggregated effects and fit diagnostics distinct. Preserve links between derived substitutions and the original observations that contributed to them. Record insufficient-data, convergence or method-eligibility failures as analysis QC with reasons. Training-data filters from the existing analysis repository require review before becoming rules for catalogue analyses; its explicit two-score-set exclusion is superseded by the general single-human-target eligibility rule.

**Curation import and public eligibility**

Supply curators with parsed metadata from GCS, retaining the full `score_urn`, `experiment_urn` and `experiment_set_urn`. Accept only a final human-reviewed curation CSV or Parquet file. Validate required fields, identifiers, ontology field pairs and duplicate `score_urn` values before importing it into BigQuery. The current curation workflow sets `is_curated=Yes` on LLM rows before review; that intermediate output must not feed `IMPORT_CURATION`.

Use one canonical boolean `is_curated` field in BigQuery. Normalize `Yes` in the reviewed export to true. A row in the reviewed export with `is_curated=Yes` is approved for publication, whether its curation began with a human or an LLM; no second approval field is required. A present row with blank or false `is_curated` is ineligible for publication. Validation must not mark rows curated automatically.

Keep curated values separate from parsed source metadata. Import by full `score_urn`: insert new rows and replace the curated fields of rows present in the reviewed export. A blank field in a present row clears its previous value, including `is_curated`; clearing or setting `is_curated` to false removes public eligibility. Rows absent from an import remain unchanged. Report inserted, replaced and withdrawn datasets. This permits a curation session containing only newly curated datasets without deleting earlier ones.

Description or ontology corrections should permit curation import and republication alone. Changes to analysis-relevant fields, such as score direction or target identity, should identify the analyses requiring reruns.

**BigQuery loading and publication**

Use internal BigQuery tables for observations, primary `score` values, mappings, variant QC, dataset QC and derived results. Curated metadata has its own authoritative table keyed by full `score_urn`. Declare each table's complete key: the primary score is keyed by observation identity and has no measurement-name dimension in this iteration. Make repeat loads replace or merge the intended development scope without appending duplicates.

Prepare public models by joining internal data to eligible curated datasets keyed by `score_urn`. These models are the common input for Postgres, Elasticsearch and downloadable files. Existing `mavedb.data` and `mavedb.metadata` contracts can become compatible curated projections if that is the least disruptive integration with the warehouse.

Update the warehouse projection and release code to consume the prepared position/ref/alt fields. Parse HGVS once upstream. Review the corresponding Postgres schema, API fields and downloadable columns together.

Reuse the existing warehouse loading implementations through a MAVE-specific publication entrypoint. All public paths must enforce the same curated population: dataset search, direct dataset searches, target summaries, counts, filtered downloads and complete dataset exports.

When a previously public score set loses curated eligibility, remove its serving rows and search documents and withdraw or block its old downloads. Publication must reconcile removals as well as additions and updates. Parsed or analyzed data being present internally never grants public eligibility.

**Implementation sequence and completion criteria**

| Phase | Work | Completion criterion |
| --- | --- | --- |
| 1. Contracts and fixtures | Define initial tables, single-human-target eligibility, coordinate conventions, QC outcomes and representative inputs | Fixtures cover simple substitutions, nonzero offsets, complex HGVS, ambiguity, mismatches, unmapped targets, extra score-file columns and excluded zero-target, multi-target and nonhuman score sets |
| 2. Ingestion | Implement dump download/reuse, eligibility selection, metadata extraction, observation parsing and the primary `score` | Original identifiers/HGVS survive; eligible parsed metadata and data are written to GCS, and excluded score URNs and reasons appear in the ingestion report |
| 3. Mapping and QC | Migrate mapping logic and implement source/UniProt sequence checks | Candidate mappings are explicit; required QC outcomes exist at variant and dataset levels |
| 4. Internal loading | Load observations, primary scores, mappings and QC into development BigQuery tables | Repeated loads are safe; eligible source and loaded row counts reconcile; no public visibility is introduced |
| 5. Curation and publication | Import a small final human-reviewed curation file and connect the curated projections to serving stores and downloads | `is_curated=Yes` entries are discoverable; uncurated entries are inaccessible through every public path; present rows replace values including blanks, absent rows stay unchanged and explicit withdrawal works |
| 6. Analysis | Migrate required R/Stan code, adapt inputs and execute through Slurm | Analysis consumes shared cloud inputs and produces loadable outputs with diagnostics and lineage |
| 7. Broader runs and polish | Expand coverage within the single-human-target scope, check resource use and document independent commands | Individual modules and composed workflows run against existing inputs with clear prerequisite checks |

Use focused checks for biological parsing/mapping, score-set eligibility, reviewed-file curation import and visibility enforcement, plus a small end-to-end fixture. Verify that multiple nucleotide observations contributing to one protein substitution remain traceable, and that ambiguous mappings do not inflate score counts. Check that dataset QC totals reconcile with retained observations, while excluded score sets are reported separately.

After the pipeline is useful and its contracts are stable, revisit versioning, incremental processing and more elaborate release coordination using the behavior observed during development.

**Existing integration points**

- [MAVE source directory](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/data_sources/mavedb/README.md).
- [Current MAVE curation helpers](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/data_exploration/curation_tools/mavedb_curation_tools.py).
- [Warehouse orchestration](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/dwh/cloudbuild.yaml).
- [Postgres projections](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/dwh/bq_to_postgres/bq_to_postgres.py).
- [Release projections](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/dwh/release/release.py).
- [Portal query implementation](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/be/data_query.py).
- [Existing MAVE architecture diagram](/Users/zakirov/Documents/GitHub/PerturbationCatalogue/development_docs/mavedb-pipeline-architecture.mmd). Align it with this plan when updating the architecture documentation, including keeping parsed metadata in GCS.
