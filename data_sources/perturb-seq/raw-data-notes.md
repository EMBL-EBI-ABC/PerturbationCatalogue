# Raw sequencing data audit — 2026-09-10

Original raw-availability audit scope: the 23 rows then labelled `Conditions = Single condition` in [datasets.tsv](datasets.tsv). No sequencing data were downloaded or reanalysed in that initial audit. The four raw-data preparation columns were updated, then `Raw data access` and `Data health check` were added. The seven other rows retain their previous raw-data fields and are marked `Not checked` in the new columns. Later BAM sampling and the schema-compatibility review are documented below; Gasperini's raw-data audit remains valid despite its subsequent complex-design classification.

## Reading the preparation table

- `Raw data access`: `FASTQ` means direct FASTQ files; `BAM only` means submitted BAMs requiring reconstruction; `SRA (FASTQ extraction)` means use SRA archives to retain technical reads omitted from the ENA FASTQs; `FASTQ + SRA` means source GEX FASTQs plus guide SRA archives. `Unavailable` means no public raw reads were found in this audit, not proof of permanent absence. `Not checked` is outside the audit scope.
- `Data health check`: `OK` means no unresolved data/metadata anomaly was found in this audit, not proof of integrity or current pipeline compatibility. `Issues` flags missing cell type (Norman and Replogle K562 essential), the inconsistent Adamson UPR catalogue accession, missing Zhu R2 guide libraries / conflicting stimulation labels, unavailable Orion reads, or Gasperini's confirmed shared-control workflow incompatibility. Details are below; missing cell types remain a minor deferred curation issue. `Not checked` must not be interpreted as `OK`.
- BAM/SRA extraction and Flex/Ultima support are processing requirements, not by themselves data-health anomalies. Thus the four Zhu R1 output rows have `OK` health but still require demultiplexing and platform support before running; the other eight Zhu rows have `Issues`.
- `ENA project` is the BioProject accession shared by ENA/SRA. `ENA sample` contains the exact, semicolon-separated BioSample accessions for both gene-expression and guide libraries, including shared pooled samples where necessary.
- `Data size TB` is the size of the selected compressed archive files in decimal TB (bytes / 10^12), rounded to three decimals. The format being counted is stated in each row. These are download sizes, not extracted FASTQ sizes or processing-space estimates. BAM indexes, processed count products and alternative representations of the same reads are excluded.
- Blank accession/size means not found or unknown, not zero. An archive record alone is not evidence that the current pipeline can process it. BAM conversion, technical barcode reads, Flex probe demultiplexing and Ultima read layouts need attention where noted.
- Zhu sizes repeat the whole shared input pool on each of its four output rows. Count each pool once, not once per donor/condition. Missing guide libraries are not included in the size.
- The size table uses decimal TB (bytes / 10^12) and states which source representation is counted; it does not estimate extracted FASTQ or processing scratch space.

## Sources and accounting

Project/sample mappings were checked against GEO, ENA run metadata and, for Zhu, the authors' sample metadata. ENA reports were fetched from `https://www.ebi.ac.uk/ena/portal/api/filereport` with `result=read_run`, the BioProject accession, and these fields:

```text
study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,sample_title,experiment_title,library_name,library_strategy,library_layout,instrument_model,fastq_ftp,fastq_bytes,submitted_ftp,submitted_bytes,sra_bytes
```

Where ENA lacks a usable representation, NCBI's locator `https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=RUN_ACCESSION` supplies file sizes and public download locations. Sum one representation per run: `type=sra` for Norman, Nadig Jurkat and Zhu guides; original source FASTQs for Zhu GEX; submitted BAMs for Adamson/Gasperini. Duplicate source checksums were checked within the Zhu pools. Representative HTTP HEAD requests returned 200 with matching content lengths for SRA/cloud and ENA files; this is an availability audit, not a complete read-integrity or barcode-recovery validation.

| Dataset / shared pool | Project | Selected runs | Exact bytes counted | Format |
| --- | --- | ---: | ---: | --- |
| Adamson pilot | PRJNA354963 | 2 | 15,134,506,196 | Submitted BAM |
| Adamson UPR | PRJNA354963 | 11 | 445,802,420,411 | Submitted BAM |
| Norman | PRJNA551220 | 32 | 252,666,450,019 | SRA archive |
| Gasperini at scale | PRJNA494734 | 64 | 705,745,934,759 | Submitted BAM |
| Nadig Jurkat | PRJNA1100571 | 1,792 | 1,073,535,459,912 | SRA archive |
| Nadig HepG2 | PRJNA1100571 | 448 | 1,967,111,796,029 | ENA FASTQ.gz |
| Replogle K562 essential | PRJNA831566 | 576 | 1,216,128,877,499 | ENA FASTQ.gz |
| Replogle RPE1 essential | PRJNA831566 | 1,792 | 1,346,749,573,628 | ENA FASTQ.gz |
| Replogle K562 genome-wide | PRJNA831566 | 4,592 | 7,742,050,722,043 | ENA FASTQ.gz |
| Zhu R1 L01–L23 | PRJNA1359008 | 246 | 23,228,534,757,708 | Source GEX FASTQ.gz + guide SRA |
| Zhu R2 L01–L24 | PRJNA1359008 | 40 | 18,794,069,701,051 | Source GEX FASTQ.gz + available guide SRA |
| Zhu R2 L25–L48 | PRJNA1359008 | 38 | 18,033,011,184,252 | Source GEX FASTQ.gz + available guide SRA |

### Adamson

[GSE90546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546) links to [PRJNA354963](https://www.ebi.ac.uk/ena/browser/view/PRJNA354963). Pilot uses GSM2406675/76 (GEX and guide barcodes). UPR uses GSM2406681 (GEX, 10X010) and GSM2406682–91 (guide gemgroups 1–10). Exclude the epistasis samples GSM2406677–80.

The catalogue's UPR associated-data entry has an inconsistent `GSM2406677` accession but links to GSM2406681 and the 10X010 H5AD. The archive's UPR title and linked GEO page identify GSM2406681 as the correct raw-data match. This audit does not modify catalogue metadata.

UPR GEX run [SRR5082094](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5082094) has no ENA file links and a zero SRA RunInfo size, but the NCBI locator exposes `10X010.bam` (442,357,276,591 bytes); its public download responded to HEAD. Do not interpret missing ENA FASTQs or a zero RunInfo size as no raw data. All selected Adamson files are BAMs; reconstruction of barcode/UMI reads for the v1 chemistry remains to be validated.

### Norman and Gasperini

[Norman GSE133344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344) / [PRJNA551220](https://www.ebi.ac.uk/ena/browser/view/PRJNA551220) contains 16 samples, two runs each: GEX gemgroups 1–8 and matching guide-barcode libraries. Match by the complete sample title, not accession order. The [SRR9602535 read specification](https://www.ebi.ac.uk/ena/browser/api/xml/SRR9602535) classifies the cell-barcode/UMI read as technical; its ENA FASTQ contains the application read only. SRA extraction must retain technical reads, so the table counts SRA archive sizes.

[Gasperini GSE120861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861) / [PRJNA494734](https://www.ebi.ac.uk/ena/browser/view/PRJNA494734) has two at-scale samples: SAMN10179954 (GEX) and SAMN10179953 (gRNA enrichment), 32 BAMs each. Pair using `at_scale_screen.<1A|1B|2A|2B>_<1-8>` in the submitted filenames. Pilot and bulk-validation samples are excluded. BAM barcode/UMI recovery needs validation before use with the FASTQ pipeline.

### Previously reprocessed datasets

Confirmed the existing sample mappings in [Nadig PRJNA1100571](https://www.ebi.ac.uk/ena/browser/view/PRJNA1100571) and [Replogle PRJNA831566](https://www.ebi.ac.uk/ena/browser/view/PRJNA831566). Sizes were recomputed from the selected raw libraries, not copied from the previous estimates.

For Replogle, select `mRNA`/`sgRNA` library names; exclude 48, 56 and 273 other runs for K562 essential, RPE1 and K562 genome-wide respectively. SAMN28561245 is the separate Ultima submission and is excluded. Existing pipeline sample TSVs contain the actual run groupings; run counts, sequencing lanes and GEM groups are distinct. Nadig Jurkat also has single-file ENA FASTQs with technical reads excluded, so its size uses SRA archives instead.

### Zhu: three shared pools, twelve output datasets

The [authors' repository](https://github.com/emdann/GWT_perturbseq_analysis_2025) links raw reads to SRP643211 / [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342), corresponding to [PRJNA1359008](https://www.ebi.ac.uk/ena/browser/view/PRJNA1359008). Pool and barcode assignments are documented in GEO's `sample probe_barcode` metadata and the [authors' sample table](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/suppl_tables/sample_metadata.suppl_table.csv).

| Library-name pool | Output rows and paired GEX/CRISPR probe barcodes | Coverage found |
| --- | --- | --- |
| CD4i_R1L01–L23 | D1 Rest: 001–004; D2 Rest: 005–008; D1 Stim8hr: 009–012; D2 Stim8hr: 013–016 | 23 GEX + 23 guide BioSamples; 246 runs |
| CD4i_R2L01–L24 | D3 Rest: 001–004; D4 Rest: 005–008; D3 Stim8hr: 009–012; D4 Stim8hr: 013–016 | 24 GEX + 16 guide BioSamples; 40 runs |
| CD4i_R2L25–L48 | D1/D2/D3/D4 Stim48hr: 001–004 / 005–008 / 009–012 / 013–016 | 24 GEX + 14 guide BioSamples; 38 runs |

Use `BC` for gene-expression probe barcodes and `CR` for the corresponding CRISPR barcodes. Retain full lane identifiers when pairing libraries; R1 guide sublibraries include suffixes such as `.1`–`.10`, including `R1L22.1_CRI_lib`. Aggregate those within their GEM lane, then demultiplex donor/condition.

All 71 GEX runs are source-only in the checked ENA/SRA run reports. The NCBI locator exposes 926 original compressed FASTQs (54,488,666,614,427 bytes); the 253 available guide runs expose SRA archives. The three pools together occupy 60,055,615,643,011 bytes, about **60.056 TB**, counted once each. Do not discard source-only GEX runs based on zero normalized SRA spots/bytes.

Known limitations before reanalysis:

- No matching R2 guide runs were found in either ENA or SRA RunInfo for L07–L14 or L25, L26, L29, L30, L32, L36, L39, L41, L47, L48. All eight output rows drawing from R2 therefore have incomplete archived guide-library coverage. Their listed sizes cover the available files, not a hypothetical complete deposition.
- GEO titles for R2 L25–L48 say **24 hr**, whereas their probe mapping names and the authors' sample table say **Stim48hr**. The table follows the explicit donor/probe mapping to the existing Stim48hr datasets and flags the conflict for resolution before running.
- These are Flex libraries sequenced on Ultima UG 100. Raw inputs pool conditions even though each intended output is a single condition. The current 10x 3-prime v3 workflow must not be assumed to support them; demultiplexing and platform-specific processing remain prerequisite work, outside this discovery task.
- Other experiments in this BioProject (Tact, Th1Th2, IL10IL21 and arrayed validation) are excluded.

### Orion: no public raw-read accession found

The [author Figshare release](https://doi.org/10.25452/figshare.plus.29190726) and [author Hugging Face release](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion) provide processed H5AD/Parquet and guide metadata. Neither provides a sequencing-read accession. NCBI SRA searches for `Xaira`, `X-Atlas/Orion` and the manuscript DOI, and ENA searches for `X-Atlas`, `Fix-Cryopreserve` and the Xaira centre, returned no matching raw-read records on the audit date. Both Orion rows retain blank project, sample and size fields with an explicit not-found note. This is a dated search result, not proof that no raw reads exist anywhere.

## Follow-up: BAM recoverability and Orion publication status, 2026-09-10

Unlike the metadata-only audit above, this follow-up inspected small HTTP byte ranges of the original BAMs: the first 256 KiB of all 77 selected files (13 Adamson; 64 Gasperini), plus 2 MiB prefixes of seven representative files. No full BAM downloads or reanalysis were performed. `samtools view -h` decoded 298,815 alignment records from the 77 smaller prefixes. All decoded records had sequence/quality and the required raw barcode/UMI tags below. Truncation errors at the end of each deliberately partial download are expected, not evidence of corrupt source files. These coordinate-sorted prefixes are not random samples and do not establish whole-file integrity, complete guide recovery, or final matrix agreement.

| Dataset | Direct observations | Remaining qualification |
| --- | --- | --- |
| Adamson pilot | Both RNA SRR5082088 and guide-barcode SRR5082089 preserve `CR/CQ` (14-base raw cell barcode/quality), `UR/UQ` (10-base UMI/quality), `BC/QT` and sequence/quality. Corrected `CB/UB` also occur. Headers identify Cell Ranger 1.0 processing. | Molecular information supports both count tracks, but this assay reads synthetic guide barcodes (GBCs), not the sgRNA protospacers directly. The GBC-to-vector/target reference still needs verification. The local supplementary CSV contains protospacers/vector IDs, not GBC sequences. |
| Adamson UPR | RNA SRR5082094 and all ten guide BAMs SRR5082095–SRR5082104 preserve the same raw tags. RNA corrected `CB` suffixes cover groups 1–10. Every separate guide BAM uses local `CB` suffix `-1`. No Adamson BAM has `@RG` headers. | Preserve/reconcile GEM identities before merging: use the guide filename's group and a validated RNA grouping rule. RNA `CB` suffixes retain useful identity, but records without `CB` need a separate rule if reconstructing all raw reads. GBC reference verification is also required. Do not run a blind pooled conversion. |
| Gasperini at-scale | All 32 RNA and 32 guide BAMs preserve `CR/CY`, `UR/UY`, sequence/quality and `RG`. Every header explicitly provides `I1(BC:QT)`, `R1(CR:CY,UR:UY)`, `R2(SEQ:QUAL)` reconstruction recipes. | Strong evidence for reconstructing both tracks with the original 10x BAMs, keeping each RNA/guide pair's `1A/1B/2A/2B_1–8` identity. This is not an end-to-end conversion test. |

The [10x converter documentation](https://github.com/10XGenomics/bamtofastq) describes `--cr11` for Cell Ranger 1.0–1.1 and warns about missing read groups in older multi-GEM BAMs. Adamson GEO sample metadata also explicitly describe raw reads, cell identities, UMIs and GEM groups in the deposited BAMs. The [Adamson paper](https://doi.org/10.1016/j.cell.2016.11.048) explains the separate GBC assay. GEO cell-identity CSVs provide assigned guide names and dominant read/UMI counts, but are not substitutes for an independently rebuilt full cell-by-probe matrix.

Scientific prioritisation is separate from recoverability. The Adamson pilot is a small historical/validation dataset (5,768 profiled cells, eight GBCs); UPR offers broader focused pathway biology. Gasperini is not a small pilot: the [paper](https://doi.org/10.1016/j.cell.2018.11.029) reports 207,324 at-scale cells, 5,779 candidate enhancers and a median of 28 guides per cell. Multi-gene perturbations are allowed by the catalogue schema. The in-depth decision below classifies Gasperini as complex because its supported comparisons use target-dependent reference populations, not because it contains multi-guide cells. BAM reconstruction remains feasible at the format level irrespective of this downstream schema limitation.

Orion: the [bioRxiv API record](https://api.biorxiv.org/details/biorxiv/10.1101/2025.06.11.659105) returns version 1, dated 2025-06-16, and `published: NA`. A title/author search found no journal version of the study. Treat both Orion datasets as from the same preprint as of this check; news coverage in GEN Edge is not a journal publication of the underlying study.

## Experimental-design label review — 2026-09-10

Reviewed all 30 rows at the level of the named dataset, not every experiment in its parent publication. This is a study/sample-metadata review, not a complete inspection of every H5AD's cell assignments. The interpretation below incorporates the user's clarification and supersedes the earlier single-target restriction. TSV annotations now reflect this interpretation; no reanalysis was started.

`Conditions` answers only: **Can this dataset be meaningfully analysed and represented in the current Perturbation Catalogue data model, APIs and UI as one shared control population and one experiment containing a flat list of perturbations, each compared against that control by DEA/GSEA?**

`Single condition` means this conceptual organisation is supported. Perturbations may be single-gene, multi-gene or another type; deliberate combinations are valid flat-list entries. Several non-targeting guides may define one shared control population. Technical lanes/GEMs and donor replication do not inherently require separate catalogue datasets. However, distinct contexts requiring different matched control populations must be separated or need schema support. Multiple targets alone are never a reason for exclusion.

This field is not a raw-data availability, BAM/FASTQ compatibility, assignment-quality or statistical-power flag. For uncertain cases, identify the missing control definition explicitly. A target-specific complement (all cells lacking target X) changes with X and is not automatically one global control. A cell carrying a non-targeting guide alongside active perturbations is not automatically unperturbed. No additional multiplicity column or single-target-only eligibility gate is needed.

| Dataset IDs (grouped only where the verdict is identical) | Review of existing annotation | Consequence for the current simple-comparison scope |
| --- | --- | --- |
| `adamson_2016_pilot` | One K562 context; individually perturbed cells pooled, eight guide-barcode identities. Label supported. | Single-target comparisons in principle; retain BAM/GBC preparation caveats above. |
| `adamson_2016_upr_perturb_seq` | One K562 context in 10X010; focused UPR-regulator screen. Do not import the drug conditions of the separate 10X005 epistasis dataset. Label supported. | Single-target comparisons in principle; resolve GBC mapping, GEM identity and the catalogue accession discrepancy. |
| `norman_2019_raw` | `Single condition`: single- and double-gene CRISPRa perturbations can be flat entries against the shared non-targeting controls. | Include combinations conceptually; no single-target-only subset is required. Correct combination labels/assignments remain a processing requirement. |
| `gasperini_2019_atscale` | Complex: supported comparisons use target-dependent reference populations, not one shared control. | In-depth assessment below resolves the earlier uncertainty. Only eight deposited cells have exclusively negative-control calls; the 50k reference is a mixed perturbed cohort. Supporting the study faithfully requires target-specific comparators or a suitable model, not just context splitting. |
| `nadig_2025_jurkat`, `nadig_2025_hepg2` | Separate low-infection-rate screens, each harvested on day 7. Dual-guide constructs do not by themselves imply different target genes. Both labels supported. | Keep the two cell lines separate and retain validated assignment/control rules. |
| `replogle_2022_k562_essential_normalized`, `replogle_2022_rpe1_essential_normalized`, `replogle_2022_k562_gw_normalized` | Separate screens at K562 day 6, RPE1 day 7 and K562 day 8 respectively, not a pooled time-course row. All three labels supported. | Same-gene dual guides are compatible with single-target comparisons; retain library and batch identity. |
| `orion_2025_hct116`, `orion_2025_hek293t` | Each row is a separate cell-line screen. Intended paired guides target the same gene; the released filtered cells have valid guide pairs. Both labels supported. | Raw-data absence remains the blocker. Validate assignments during any future processing; multiple targets alone do not make an assignment invalid or schema-incompatible. |
| All twelve `zhu_2025_D[1-4]_{rest,stim8hr,stim48hr}_cl` rows | Each named output is one donor and one stimulation state/timepoint. Low-MOI library delivery and the author sample table support these labels at output level. | Raw pools contain multiple outputs and must be demultiplexed. Compare perturbations with controls from the same donor/state, not stimulated perturbations with resting controls. The four 48-hour outputs retain the unresolved GEO 24-hour versus author 48-hour label conflict. |
| `adamson_2016_upr_epistasis` | Requires splitting DMSO/tunicamycin/thapsigargin contexts, each with its own matched genetic controls. | Combinatorial UPR perturbations are allowed. Once split by treatment, each flat perturbation list can conceptually fit the schema; combinations are not an additional exclusion. |
| `song_2025_jurkat_hiv` | Complex label supported: unstimulated DMSO, PMA/ionomycin followed by GFP-positive sorting, and PMA/ionomycin followed by GFP-negative sorting. | GFP positivity is an outcome-based selection, not a separately administered treatment. Need context-matched genetic controls and care interpreting sorted-population effects. |
| `frangieh_2021_raw` | Complex label supported: untreated, IFN-gamma-treated, and IFN-gamma-pretreated/autologous-TIL co-culture survivor populations. | Stratify contexts before genetic comparisons. Co-culture selection and different collection timing are additional interpretation concerns. The model is patient-derived cultured melanoma, not freshly isolated cells; broader model-system curation is outside this review. |
| `arce_2025` | Complex label supported: two donors, Treg/Teff, resting/48-hour restimulation with CD3/CD28/CD2 in the Perturb-CITE-seq experiment. | Eight donor-by-cell-type-by-state strata in principle. Use experimental/sample metadata for identity; surface proteins alone need not define the split. |
| `zhu_2025_pseudobulk` | Complex label supported: four donors and three states, aggregated from the same underlying screen represented by the twelve single-cell outputs. | Not an independent additional raw screen; preserve donor/state aggregation and avoid double-counting. |
| `datlinger_2017` | Complex label supported for the named Jurkat dataset: untreated and anti-CD3/anti-CD28-stimulated cells. | Two potential single-context outputs, each requiring matched genetic controls. Exclude unrelated HEK293T/3T3 assay-validation samples in the parent GEO series. |
| `jiang_2025` | Complex label supported: six cell lines crossed with five pathway-specific cytokine/stimulus experiments, 24-hour stimulation. | Thirty line-by-stimulus contexts in principle, with pathway-specific guide libraries and non-targeting controls. These are not five cytokines administered together; do not assume an unstimulated arm is required or deposited for every context. |

Summary: **22 rows retain `Single condition`, including Norman's combinations; eight are complex: seven require context-specific controls/splitting and Gasperini requires target-dependent comparators.** Gasperini's classification is not a ban on high MOI or multiple targets. Multi-context datasets may become compatible after splitting; splitting Gasperini by culture context does not resolve its comparator structure. Shared-control eligibility does not remove the raw-data and processing caveats elsewhere in these notes.

### Shared-control follow-up

The follow-up review found no additional Gasperini-like comparator exception among the 22 `Single condition` rows. Adamson, Norman, Nadig, Replogle, Orion and the donor/state-specific Zhu outputs support deliberate negative-control-only populations at the study-design level. This is not an independent audit of final control purity, batch representation or post-filtering counts in every output; those remain processing checks. Control-only means assigned exclusively negative-control guides, not untouched by the experimental machinery.

For Zhu, summing `NTC single sgRNA` across lanes in the [author QC table](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/suppl_tables/QC_summaries_per_sample_lane.csv) gives 68,441–88,354 cells per donor/state output across all twelve outputs. The table distinguishes this category from multi-guide and unassigned cells. These are author-reported assignment counts, not counts after our future pipeline filtering.

## Gasperini shared-control decision

**Classification: complex under the current catalogue schema.** `gasperini_2019_atscale` is not suitable for presentation as a flat perturbation list against one shared control population using straightforward DEA/GSEA. This is a resolved eligibility decision for the deposited at-scale dataset, not a claim that multi-gene perturbations are forbidden, that the study is unusable, or that no conceivable redesigned analysis could use a fixed reference.

### Comparison actually supported by the study

The published differential-expression method labels cells according to presence or absence of each tested gRNA group. Thus the reference population for perturbation A differs from that for B. The at-scale regression also accounts for guide count, mitochondrial fraction and preparation batch. Negative-control guide groups are themselves tested for effects; their existence does not establish a cohort carrying only control guides.[1]

The structural difference is:

| Current catalogue organisation | Gasperini study comparison |
| --- | --- |
| Perturbation A versus shared control C | A-bearing cells versus cells lacking A |
| Perturbation B versus the same C | B-bearing cells versus cells lacking B |

This remains a distinction even when A and B denote multi-gene perturbations. A comparison against a mixed, fixed background would be a new estimand, not simply another encoding of the published target-specific comparison.

### Full deposited cell-assignment audit

The entire deposited `GSE120861_at_scale_screen.phenoData.txt.gz` was inspected, not a cell subsample. It has 207,324 unique cell IDs and 18 whitespace-delimited fields, with no header. Guide sequences in field 7 were matched to `GSE120861_grna_groups.at_scale.txt.gz`; their number was checked against field 11 wherever that field was not `NA`. All assigned sequences matched the 13,189-entry guide dictionary.[2][3]

Negative-control membership used the dictionary's `random_*`, `scrambled_*` and `bassik_mch` groups: 101 guides across 51 groups. This naming interpretation is independently documented in the Katsevich laboratory's reanalysis import code/report, which describes the paired negative controls and the single-guide `bassik_mch` group.[5] All remaining guides were treated as targeting; “targeting” here means an assigned targeting construct, not proven biological knockdown.

| Measured category | Cells |
| --- | ---: |
| All deposited cells | 207,324 |
| At least one assigned targeting guide | 205,789 |
| At least one negative-control guide, with or without targeting guides | 43,302 |
| Negative-control guide(s) **and** targeting guide(s) | 43,294 |
| Exclusively negative-control guide calls | **8** |
| Missing guide assignment/count (`NA`) | 1,527 |

The disjoint categories reconcile exactly: 205,789 + 8 + 1,527 = 207,324. There were no explicit zero-guide rows among non-missing assignments. Each of the eight negative-control-only cells had just one called guide; they span six of the 32 GEM libraries (`1A_3`, `1B_3`, `1A_4`, `1B_6`, `2B_5`, `2B_8`). These are candidate control-only *calls*, not independently validated unperturbed cells. Missing assignments are not evidence of unperturbed status, and undetected co-perturbations cannot be ruled out by these metadata.

**Assessment:** eight such candidates do not establish a defensible shared reference for broad, batch-aware DEA/GSEA of this entire screen. This is a scientific adequacy judgment, not a universal numerical rule that eight controls can never support any experiment. Almost all negative-control-bearing cells also have targeting guides, so treating all 43,302 as unperturbed controls would be incorrect.

As a separate check on representing complete multi-target combinations as flat entries, negative-control guides were removed and `_top_two`/`_second_two` guide groups were collapsed to their common target label. The remaining 205,789 cells had 203,468 distinct called target sets; 201,398 sets occurred once and the largest occurred six times. This does not ban combinations, but demonstrates that relabelling the random combinations does not produce a well-replicated flat experiment with a validated common control. These are deposited calls, not a fresh raw-read assignment.

### The 50k reference file is not an unperturbed control cohort

GEO describes `GSE120861_50k_reference_cells.rds.gz` as the IDs used for at-scale differential-expression testing to reduce computational cost.[4] The doubly gzip-compressed RDS contains exactly 50,000 distinct cell IDs, all present in the at-scale phenotype file. Joining them to the deposited assignments gives:

| Reference-file membership | Cells |
| --- | ---: |
| At least one targeting guide | **49,643** |
| Exclusively negative-control guide calls | 2 |
| Missing guide assignment | 355 |

For example, the file includes 538 cells assigned `FKBP2_TSS` and 533 assigned `TOP1_TSS`. Consequently it overlaps tested perturbation populations. Its filename must not be interpreted as an independently sampled, unperturbed control arm. The GEO description establishes its computational purpose but does not fully specify the selection algorithm; no claim is made here that its precise selection or exclusion logic was recovered from original author code. The paper's presence/absence comparison and the measured composition are sufficient to reject the file as evidence for an ordinary shared control cohort.

### Alternatives considered and decision boundary

- **Use every negative-control-bearing cell:** this includes 43,294 co-perturbed cells and does not define an unperturbed population. A mixed-background interpretation would require an explicitly justified new analysis, including treatment/control overlap handling.
- **Use only the eight negative-control-only calls:** too little validated, library-representative evidence to establish an adequate common reference for the complete screen. This is not a supported drop-in reanalysis route.
- **Use the 1,527 unassigned cells:** confuses missing perturbation measurement with absence of perturbation.
- **Use the 50k file:** supplies a mixed cohort overlapping target-positive groups, not an unperturbed arm. Removing the current target's carriers makes the reference target-dependent again.
- **Borrow the low-MOI pilot or bulk validation controls:** these belong to different experiments/library designs and are not an established matched single-cell control arm for the at-scale screen.
- **Support target-specific reference sets or fitted model contrasts:** scientifically appropriate directions for later investigation, but require representing the actual comparator/estimand rather than a single global control in the current catalogue.

The dataset therefore moves from “shared control unresolved” to **“Complex: target-dependent reference populations; no adequate shared control established for catalogue DEA/GSEA.”** `BAM only`, accessions and download sizes are unchanged. `Data health check = Issues` continues to flag the shared-control workflow incompatibility, not damaged BAMs; raw reconstruction and source-file integrity are separate questions. No raw-read conversion, DEA/GSEA, new control cohort or author contact was undertaken.

### Sources and reproducibility

Assessment performed 2026-09-10. Counts above are direct calculations from the deposited files; methodological conclusions are distinguished from those measurements. The scope is the named at-scale screen, not the separate pilot or bulk validation experiments. The cell metadata and reference list were downloaded; complete sequencing reads and expression matrices were not needed.

1. Gasperini et al., *A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens*, Cell (2019), [author-hosted manuscript](https://andrewjohnhill.com/documents/publications/enhancer_screen_cell.pdf), STAR Methods “Differential expression tests”, PDF page 21. [Published DOI](https://doi.org/10.1016/j.cell.2018.11.029).
2. Gasperini et al./NCBI GEO, [at-scale cell phenotype metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.phenoData.txt.gz), GSE120861. SHA-256: `16a9cd67b440808f127ae15223378b1594a6d84c6cedefacc6a288ba29c0600d`.
3. Gasperini et al./NCBI GEO, [at-scale guide dictionary](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_grna_groups.at_scale.txt.gz). SHA-256: `40e996f072c0ca3664d6aea92f2b85d996f6ce0c17f0808ef09ab4ff8953da82`.
4. Gasperini et al./NCBI GEO, [GSM3417255 metadata and file descriptions](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3417255), and [50k reference IDs](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_50k_reference_cells.rds.gz). SHA-256 of the downloaded reference file: `cec29ac8260d92e4b9cbdad90596f2763c1c308912a352030eeffc75708dcacc`.
5. Katsevich Lab, [Gasperini 2019 v2 import report](https://github.com/Katsevich-Lab/import-gasperini-2019-v2) and [processing code](https://github.com/Katsevich-Lab/import-gasperini-2019-v2/blob/main/at-scale/process_data_2.R). This is reanalysis code, not original Gasperini author code; used to cross-check guide-group interpretation, not to infer the original reference-selection algorithm.

### Sources for the broader design review

- Adamson: [study](https://doi.org/10.1016/j.cell.2016.11.048), [GSE90546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546), with pilot GSM2406675, UPR GSM2406681 and epistasis GSM2406677 explicitly separated.
- Norman: [GSE133344 overall design](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344) explicitly describes overexpressing genes alone or in combination; [study](https://doi.org/10.1126/science.aax4438).
- Gasperini: [study](https://doi.org/10.1016/j.cell.2018.11.029), including its at-scale high-MOI design.
- Nadig: [GSE264667 overall design](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667) identifies parallel low-infection-rate Jurkat/HepG2 screens and day-7 collection.
- Replogle: [author release](https://doi.org/10.25452/figshare.plus.20029387) identifies all three line/timepoint combinations; [study](https://doi.org/10.1016/j.cell.2022.05.013).
- Orion: [preprint](https://doi.org/10.1101/2025.06.11.659105) and [author dataset card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion) describe guide pairs and filtered outputs.
- Zhu: [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342), [representative sample protocol](https://www.ncbi.nlm.nih.gov/sra/SRX31479981), and [author sample metadata](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/suppl_tables/sample_metadata.suppl_table.csv). Existing barcode-pool evidence is above.
- Song: [study](https://doi.org/10.1038/s41556-025-01626-9), Figure 5 and [GSE247599](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247599).
- Frangieh: [author manuscript](https://ccsp.hms.harvard.edu/wp-content/uploads/2021/10/Frangieh-2021.pdf), especially Methods, “Perturb-CITE-Seq co-culture experiment”.
- Arce: [study](https://doi.org/10.1038/s41586-024-08314-y), Perturb-seq methods, and [GSE278572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278572). Do not substitute donor counts from other assays in this paper.
- Datlinger: [GSE92872 overall design](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92872).
- Jiang: [GSE281048 overall design](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281048) and [study](https://doi.org/10.1038/s41556-025-01622-z).

Orion release expectation: publication may strengthen incentives or requirements for raw-data deposition, but it does not establish that public FASTQs will follow. The destination journal and any author commitment are unknown. [Nature Portfolio policy](https://www.nature.com/ng/editorial-policies/reporting-standards) illustrates sequencing-data deposition requirements; it is not evidence that Orion will publish there or that a specific format will be released. Continue treating raw-read availability as unconfirmed until a release/accession is found.
