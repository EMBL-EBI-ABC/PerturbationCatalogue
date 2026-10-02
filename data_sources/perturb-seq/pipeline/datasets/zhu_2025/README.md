# Zhu 2025 Flex/Ultima read-level inputs

The twelve `zhu_2025_D[1-4]_{rest,stim8hr,stim48hr}_cl` outputs are donor/state
subsets of three multiplexed physical pools. The authors' sample table identifies
the kit as `GEMX_flex_v1` and sequencing platform as `Ultima`; the pipeline counts
the original reads and does not use author gene-count or guide-count matrices.
Author H5ADs remain independent comparison references only.

## Read-level route

Each physical lane is counted once for all 16 paired GEX/CRISPR probe-barcode
aliases. The lane task streams each verified original GEX FASTQ pair through
the native Cyto 0.4.5 Flex probe counter, then downloads each pinned guide SRA
archive one at a time. The guide archive's `SEQUENCE` records are stored in two
equal blocks: the first block contains 28-base cell-barcode/UMI reads and the
second contains variable-length guide reads. `stream_sra_pairs.py` reconstructs
mates in lockstep and fails on unequal row counts or mismatched spot names.
R1 carries the 16-base cell barcode and 12-base UMI. In a representative
998-pair check, guide R2 lengths varied from 47 to 91 bases (mostly 58); 893
reads contained the expected 30-base anchor and 926 contained an exact
reverse-complement match to an author guide. The nominal R2 layout is the 8-base
CR barcode, capture anchor and 20-base guide sequence. Variable lengths are
retained for Cyto's anchor-aware mapping rather than forced to 58 bases.

GEX counts use the pinned 10x Human Transcriptome Probe Set v1.1.0 for GRCh38
2024-A, restricted to probes marked `included=TRUE`. Separate barcode maps
demultiplex the 8-base `BC` GEX barcode and 8-base `CR` guide barcode. All GEX
and guide IBU records are merged within a physical lane before UMI correction
and counting. Lane products preserve the full 16-base cell barcode plus GEX
barcode and `lane_id`; this prevents two probe-barcode aliquots or different
lanes from collapsing during author comparison.

Cyto 0.4.5 is the pinned read-level counter. Its probe/UMI collision behavior
is not bit-for-bit Cell Ranger Flex behavior: tied competing probe identities
for the same CBC and UMI can be discarded, a unique dominant identity may be
retained, and different UMI sequences remain distinct. The compact native
fixture in `../../tests/test_flex_native.py` is the first group in the Cyto
UMI stream and returns one for its same-UMI/two-probe case; a broader native
characterization found later tied groups are dropped. This fixture-specific
first-group result is not a rule for every cross-probe collision. Cell Ranger
documents summing counts across probe pairs targeting the same gene; see
[Cell Ranger's Flex algorithm](https://www.10xgenomics.com/support/cn/software/cell-ranger/latest/algorithms-overview/cr-flex-frp-algorithm)
for the 10x counting semantics.

## Coverage and donor assignment

`samples.tsv` is generated from the author's pinned sample metadata and maps each
donor/state to its exact physical pool, lane set and BC/CR pairs. R1 L01–L23 has
complete guide archive coverage. R2 L01–L24 and L25–L48 have missing guide
lanes; `guide_coverage.tsv` records every lane, and outputs preserve the status
in `obs` and H5AD provenance. For a lane with no archived guide source, its
guide matrix is zero-filled and explicitly marked `no_archived_guide_sra`; that
matrix must not be read as observed absence of guide counts. Missing reads are
never replaced with author assignments. For the R2 L25–L48 conflict, inputs
follow the explicit author donor/probe mapping while the GEO title discrepancy
remains recorded for later resolution.

The raw source manifests contain 463 original paired GEX FASTQ sources and 253
guide SRA archive pins (accession, source bytes and MD5). Every GEX file is
validated against its pinned byte length and MD5. The NCBI SRA locator's size
and MD5 are checked before guide archive transfer, and the completed archive is
verified before read extraction. The first complete output to validate is
`zhu_2025_D1_rest_cl`; sibling donor/state outputs reuse the same per-lane
Nextflow cache.

## Pinned references and regeneration

`generate_inputs.py --check` validates the committed manifests and deterministic
feature/sample tables. Run it without `--check` to regenerate the derived TSVs.
The inputs pin:

- [Author repository revision `aa5c84a973c0e1a090b0072dc5b080bf7fbbed38`](https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/aa5c84a973c0e1a090b0072dc5b080bf7fbbed38), including `sample_metadata.suppl_table.csv` (SHA-256 `766134d11dabb5d63388e00b4d687a809c0254eb3c19d77f439a6bf9dad7abd4`) and `sgRNA_library_metadata.suppl_table.csv` (SHA-256 `00a1bec2afc2082fc79765531696d7e22672a8ba904ea54c035858f425a657a8`).
- The [10x Flex v1 probe-set reference](https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-probe-set-files), distributed here as `Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv` (MD5 `8d071b87b07a98cc7aabd6dcad526fef`).
- `guide_sequences.tsv`, checked row-by-row against the pinned author guide table; `guide_targets.tsv` SHA-256 `fa5fd9c8c7aae7ff2860c88ce961f0b0a36a5401db14a8dc289fb3c2f7c1243d`.
- `gex_sources.tsv`, `guide_sources.tsv`, and `guide_coverage.tsv`, which pin the actual raw sources and lane mapping; no expiring signed download URLs are stored.
- `probe_barcodes.tsv`, which pairs each 8-base BC and CR sequence with its named alias and matches the barcode sequences in the pinned Cyto resources.

The libraries are identified as [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342) / [PRJNA1359008](https://www.ebi.ac.uk/ena/browser/view/PRJNA1359008). The GEO titles for R2 L25–L48 say 24 hr while the author table and donor/probe mapping say `Stim48hr`; this pipeline follows the explicit author mapping and retains that unresolved label conflict in [the raw-data notes](../../../raw-data-notes.md).
