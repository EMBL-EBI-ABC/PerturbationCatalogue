# Zhu 2025 Flex/Ultima read-level inputs

The twelve `zhu_2025_D[1-4]_{rest,stim8hr,stim48hr}_cl` outputs are donor/state
subsets of three multiplexed physical pools. The authors' sample table identifies
the kit as `GEMX_flex_v1` and sequencing platform as `Ultima`; the pipeline counts
the original reads and does not use author gene-count or guide-count matrices.
Author H5ADs remain independent comparison references only.

## Read-level route

Each physical lane is counted once for all 16 paired GEX/CRISPR probe-barcode
aliases with kb-python 0.30.2, Kallisto 0.52.0 and Bustools 0.45.1. Each original
GEX FASTQ pair is size/MD5-verified, pseudoaligned and removed before the next
pair. Guide archives are downloaded one at a time; `stream_sra_pairs.py`
reconstructs their two equal sequence blocks, checking every pair's spot name
and the complete record counts. R1 contains CBC16 followed by UMI12. Guide R2
contains CR8 followed by the capture anchor and reverse-complement guide;
its observed 47–91-base lengths are preserved. GEX R2 contains the 50-base
probe target at bases 0:50 and BC8 at bases 68:76.

The included 10x probe pairs are indexed as unique probe IDs. Guide targets
contain the 30-base capture anchor plus the reverse-complement 20-base guide.
Native component-wise barcode correction uses the pinned CBC whitelist and
all 128 raw GEX BC variants; corrected BCs become one of 16 canonical aliases.
Guide CR aliases are translated to the paired canonical GEX BC sequence.
This preserves `CBC16+BC8+lane_id` as the cell identity for both modalities.

Source-specific Kallisto equivalence classes are reconciled by their target
sets before BUS records are combined. Bustools deduplicates each probe's UMIs
across all sources in a lane, then sparse probe counts are summed to stable
ENSG genes. The same UMI on two distinct probe targets can contribute one
count per probe. Multi-target classes receive no EM or multimapping allocation;
ambiguity is recorded separately from pseudoalignment rates.

Before common QC, the native GEX allowlist/correct/sort/count policy runs per
lane and canonical BC alias. Guide counts have no independent abundance
filter and are joined to RNA cells by their full composite barcode. The raw
H5AD therefore follows the standard branch's native barcode prefilter policy.
`../../tests/test_kb_flex_native.py` validates barcode correction, alias
separation, ambiguous-barcode rejection, cross-source UMI deduplication and
probe-to-gene summation inside the pinned image. Full production-pool
validation remains a separate check.

## Coverage and donor assignment

`samples.tsv` is generated from the author's pinned sample metadata and maps each
donor/state to its exact physical pool, lane set and BC/CR pairs. R1 L01–L23 has
one currently unavailable guide source for lane 13 (see
[unavailable_guide_sources.tsv](unavailable_guide_sources.tsv)); other R1 guide
source rows remain available. R2 L01–L24 and L25–L48 have missing guide lanes;
`guide_coverage.tsv` records every lane, and outputs preserve the status in
`obs` and H5AD provenance. For a lane with no archived guide source, its guide
matrix is zero-filled and explicitly marked `no_archived_guide_sra`; that matrix
must not be read as observed absence of guide counts. For lanes marked `partial`,
the available guide sources do not cover the full expected lane, so missing
guide counts are unknown rather than measured zero. Missing reads are never
replaced with author assignments. For the R2 L25–L48 conflict, inputs
follow the explicit author donor/probe mapping while the GEO title discrepancy
remains recorded for later resolution.

The raw source manifests contain 463 original paired GEX FASTQ sources and 253
guide SRA archive pins (accession, source bytes and MD5). Every GEX file is
validated against its pinned byte length and MD5. The NCBI SRA locator's size
and MD5 are checked before guide archive transfer, and the completed archive is
verified before read extraction. The first full invocation processes the smallest pool, `CD4i_R2_L25-48`,
producing `zhu_2025_D1_stim48hr_cl`, `zhu_2025_D2_stim48hr_cl`,
`zhu_2025_D3_stim48hr_cl` and `zhu_2025_D4_stim48hr_cl` together. The other
two pools follow after its four outputs validate and intermediates are removed.

## Pinned references and regeneration

`generate_inputs.py --check` validates the committed manifests and deterministic
feature/sample tables and native references. Run it without `--check` to
regenerate the derived files.
The inputs pin:

- [Author repository revision `aa5c84a973c0e1a090b0072dc5b080bf7fbbed38`](https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/aa5c84a973c0e1a090b0072dc5b080bf7fbbed38), including `sample_metadata.suppl_table.csv` (SHA-256 `766134d11dabb5d63388e00b4d687a809c0254eb3c19d77f439a6bf9dad7abd4`) and `sgRNA_library_metadata.suppl_table.csv` (SHA-256 `00a1bec2afc2082fc79765531696d7e22672a8ba904ea54c035858f425a657a8`).
- The [10x Flex v1 probe-set reference](https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-probe-set-files), distributed here as `Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv` (MD5 `8d071b87b07a98cc7aabd6dcad526fef`).
- `guide_sequences.tsv`, checked row-by-row against the pinned author guide table; `guide_targets.tsv` SHA-256 `fa5fd9c8c7aae7ff2860c88ce961f0b0a36a5401db14a8dc289fb3c2f7c1243d`.
- `gex_sources.tsv`, `guide_sources.tsv`, and `guide_coverage.tsv`, which pin the actual raw sources and lane mapping; no expiring signed download URLs are stored.
- `probe_barcodes.tsv`, pairing the 16 canonical BC/CR aliases; `bc_barcode_variants.tsv` retains all 128 GEX barcode variants.
- The independent CBC/BC/CR resource snapshots from the [Cyto 0.4.5 resource archive](https://github.com/ArcInstitute/cyto/releases/download/cyto-0.4.5/cyto-resources.tar.gz), SHA-256 `f004e5eb5e2020c3d5c894f3c19838bbfbdad548be722f875c6a3afceace13ca`. These are reference data and require no Cyto executable. `cbc_whitelist.txt.gz` has SHA-256 `6b18ff5b43f51665496a09388cca4a8cae58b22f05d8aca5fc264fa8e1600b06`; it contains 737,280 CBC16 sequences and differs from the generic 737K August 2016 list.
- Generated `gex_probe_targets.fa`, self-mapping `gex_probe_t2g.tsv`, `gex_probe_to_gene.tsv`, anchored `guide_feature_targets.fa` and self-mapping `guide_features_t2g.tsv`. There are 53,459 included GEX probes for 18,129 ENSG genes and 26,504 guides.

The libraries are identified as [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342) / [PRJNA1359008](https://www.ebi.ac.uk/ena/browser/view/PRJNA1359008). The GEO titles for R2 L25–L48 say 24 hr while the author table and donor/probe mapping say `Stim48hr`; this pipeline follows the explicit author mapping and retains that unresolved label conflict in [the raw-data notes](../../../raw-data-notes.md).
