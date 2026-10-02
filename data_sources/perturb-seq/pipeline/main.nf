#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// =============================================================================
// PIPELINE PARAMETERS
// =============================================================================
params.dataset_id = null
params.assay_mode = "standard"
params.curated_h5ad = null
params.gmt = null
params.comparison_outdir = null
params.batch_size = 50
params.min_cells_per_perturbation = 10
params.limit_perturbations = 0
params.target_sum = 10000
params.matrix_key = "X"
params.gene_map = ""
params.gsea_permutations = 1000
params.gsea_min_size = 15
params.gsea_max_size = 500
params.gsea_seed = 1
params.tie_correct = false
params.sra_bin = ""
params.sample_sheet = null
params.outdir = "results"
params.chemistry = "10xv3"
params.limit = 0
params.concat_max_loaded_elems = 100000000
params.h5repack_filter = "GZIP=4"
params.split_bam_source = ""
params.split_bam_groups = 0
params.cell_id_columns = ""
params.flex_pool_id = ""
params.flex_curated_h5ads = ""
params.flex_samples = ""
params.flex_gex_sources = ""
params.flex_guide_sources = ""
params.flex_guide_coverage = ""
params.flex_gene_probes = ""
params.flex_gene_count_features = ""
params.flex_guide_features = ""
params.flex_guide_count_features = ""
params.flex_guide_targets = ""
params.flex_probe_barcodes = ""
params.flex_source_cache = ""
params.flex_max_forks = 2
params.flex_gex_probe_targets = ""
params.flex_gex_t2g = ""
params.flex_gex_probe_to_gene = ""
params.flex_guide_feature_targets = ""
params.flex_guide_t2g = ""
params.flex_cbc_whitelist = ""
params.flex_bc_barcode_variants = ""

// cDNA Reference parameters
params.transcriptome_fa = null
params.gtf = null

// KITE Guide Reference parameters
params.features_tsv = null

// =============================================================================
// PROCESSES
// =============================================================================

/**
 * Builds the Kallisto index for the standard cDNA/mRNA workflow.
 */
process BUILD_INDEX_STANDARD {
    tag "cDNA_index"
    publishDir "${params.outdir}/reference/standard", mode: 'copy'

    input:
    path transcriptome_fa
    path gtf

    output:
    path "index.idx", emit: index
    path "t2g.txt", emit: t2g

    script:
    """
    kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa ${transcriptome_fa} ${gtf}
    """
}

/**
 * Builds the Kallisto index for the KITE (Guide RNA) workflow.
 */
process BUILD_INDEX_KITE {
    tag "KITE_index"
    publishDir "${params.outdir}/reference/kite", mode: 'copy'

    input:
    path features_tsv

    output:
    path "mismatch.idx", emit: index
    path "mismatch_t2g.txt", emit: t2g
    path "mismatch.fa", emit: fa

    script:
    """
    kb ref -i mismatch.idx -g mismatch_t2g.txt -f1 mismatch.fa --workflow kite ${features_tsv}
    """
}

/** Downloads and counts all runs for one sample's standard library. */
process KB_COUNT_STANDARD {
    tag "standard_${sample_id}"
    publishDir "${params.outdir}/counts_standard/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(accessions)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_filtered/adata.h5ad"), emit: h5ad
    path "stream_metrics.json", emit: metrics

    script:
    def sraArg = params.sra_bin ? "--sra-bin '${params.sra_bin}'" : ""
    def sourceArgs = accessions.collect { "'${it}'" }.join(' ')
    """
    python ${projectDir}/bin/stream_count.py \
      --accessions ${sourceArgs} \
      --index ${index} --t2g ${t2g} --chemistry ${chemistry} \
      --workflow standard --cpus ${task.cpus} ${sraArg}
    """
}

/** Count one staged BAM subset emitted by the generic GEM-group splitter. */
process KB_COUNT_STANDARD_BAM {
    tag "standard_${sample_id}"
    publishDir "${params.outdir}/counts_standard/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_filtered/adata.h5ad"), emit: h5ad
    path "stream_metrics.json", emit: metrics

    script:
    def sraArg = params.sra_bin ? "--sra-bin '${params.sra_bin}'" : ""
    """
    python ${projectDir}/bin/stream_count.py \
      --accessions 'BAMFILE:${bam}' \
      --index ${index} --t2g ${t2g} --chemistry ${chemistry} \
      --workflow standard --cpus ${task.cpus} ${sraArg}
    """
}

/** Split one shared BAM in one pass into independently countable GEM groups. */
process SPLIT_BAM_GEM_GROUPS {
    tag "split_${source}"
    cpus { Math.max(16, groups + 6) }

    input:
    tuple val(source), val(groups)

    output:
    path "groups/group_*.bam", emit: bams
    path "groups/split_metrics.json", emit: metrics

    script:
    if (task.cpus < groups + 3)
        error "BAM splitting needs at least ${groups + 3} CPUs for ${groups} groups"
    def samtoolsThreads = task.cpus - groups - 2
    """
    python ${projectDir}/bin/bam_to_fastq.py \
      --source '${source}' --split-groups ${groups} \
      --download-source \
      --samtools-threads ${samtoolsThreads} --output-dir groups
    """
}

/** Downloads and counts all runs for one sample's kite library. */
process KB_COUNT_KITE {
    tag "kite_${sample_id}"
    publishDir "${params.outdir}/counts_kite/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(accessions), val(feature_offset), val(feature_length)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_unfiltered/adata.h5ad"), emit: h5ad
    path "stream_metrics.json", emit: metrics

    script:
    def sraArg = params.sra_bin ? "--sra-bin '${params.sra_bin}'" : ""
    def sourceArgs = accessions.collect { "'${it}'" }.join(' ')
    """
    python ${projectDir}/bin/stream_count.py \
      --accessions ${sourceArgs} \
      --index ${index} --t2g ${t2g} --chemistry ${chemistry} \
      --workflow kite --feature-offset ${feature_offset} \
      --feature-length ${feature_length} --cpus ${task.cpus} ${sraArg}
    """
}

/** Index the pinned probe-pair and anchored guide targets once for all pool lanes. */
process BUILD_INDEX_FLEX {
    tag "Flex_probe_and_guide_indexes"
    cpus 8
    memory '32 GB'
    time '4h'

    input:
    path gex_probe_targets
    path guide_feature_targets

    output:
    path "gex_probe.idx", emit: gex_index
    path "guide_features.idx", emit: guide_index

    script:
    """
    kb ref --workflow custom -i gex_probe.idx -k 31 ${gex_probe_targets}
    kb ref --workflow custom -i guide_features.idx -k 31 ${guide_feature_targets}
    """
}

/** Count every source run for one physical Flex lane, then deduplicate within its barcode aliases. */
process KB_FLEX_LANE {
    cache 'deep'
    stageInMode 'copy'
    tag "${lane_id}"
    cpus 16
    memory { 128.GB * task.attempt }
    time '120h'
    maxForks params.flex_max_forks.toInteger()

    input:
    tuple val(pool_id), val(lane_id)
    path gex_sources
    path guide_sources
    path guide_coverage
    path samples
    path barcode_aliases
    path gene_probes
    path gene_count_features
    path guide_features
    path guide_count_features
    path guide_targets
    val source_cache_dir
    path gex_index
    path gex_t2g
    path gex_probe_to_gene
    path guide_index
    path guide_t2g
    path cbc_whitelist
    path bc_barcode_variants
    val sra_bin
    path helper_code, stageAs: "scripts/*"

    output:
    tuple val(pool_id), val(lane_id), path("flex_lane/lane_counts/*.h5ad"), path("flex_lane/lane_metrics.json"), emit: counts

    script:
    def sraRoot = sra_bin ? sra_bin : "/opt/sratoolkit.3.4.1-ubuntu64/bin"
    def cacheArg = source_cache_dir ? "--source-cache-dir '${source_cache_dir}'" : ""
    """
    rm -rf -- flex_lane
    python scripts/kb_flex_lane.py \
      --pool-id '${pool_id}' --lane-id '${lane_id}' \
      --gex-sources ${gex_sources} --guide-sources ${guide_sources} \
      --guide-coverage ${guide_coverage} --samples ${samples} \
      --barcode-aliases ${barcode_aliases} \
      --gene-probes ${gene_probes} --gene-count-features ${gene_count_features} \
      --guide-features ${guide_features} --guide-count-features ${guide_count_features} \
      --guide-targets ${guide_targets} \
      --gex-index ${gex_index} --gex-t2g ${gex_t2g} \
      --gex-probe-to-gene ${gex_probe_to_gene} \
      --guide-index ${guide_index} --guide-t2g ${guide_t2g} \
      --cbc-whitelist ${cbc_whitelist} --bc-barcode-variants ${bc_barcode_variants} --sra-bin '${sraRoot}' \
      ${cacheArg} --output-dir flex_lane --threads ${task.cpus}
    """
}

/** Reassemble one donor/state from reusable lane H5ADs using bounded on-disk concatenation. */
process AGGREGATE_FLEX_SAMPLE {
    cache 'deep'
    tag "${sample_id}"
    cpus 8
    memory { 128.GB * task.attempt }
    time '48h'

    input:
    tuple val(pool_id), val(sample_id), path(lane_counts)
    path samples
    path guide_coverage
    path gex_sources
    path guide_sources
    path gene_probes
    path probe_barcodes
    path gene_count_features
    path guide_features
    path guide_count_features
    path guide_targets
    val max_loaded_elems
    path helper_code, stageAs: "scripts/*"

    output:
    tuple val(sample_id), path("flex_sample_uncompressed.h5ad"), emit: h5ad
    path "flex_sample_metrics.json", emit: metrics

    script:
    """
    mkdir task_code
    cp -L scripts/* task_code/
    python task_code/flex_aggregate.py \
      --pool-id '${pool_id}' --sample-id '${sample_id}' \
      --samples ${samples} --guide-coverage ${guide_coverage} \
      --gex-sources ${gex_sources} --guide-sources ${guide_sources} \
      --gene-probes ${gene_probes} --probe-barcodes ${probe_barcodes} \
      --gene-count-features ${gene_count_features} \
      --guide-features ${guide_features} --guide-count-features ${guide_count_features} \
      --guide-targets ${guide_targets} \
      --lane-counts ${lane_counts.join(' ')} \
      --max-loaded-elems ${max_loaded_elems} \
      --output flex_sample_uncompressed.h5ad --metrics flex_sample_metrics.json
    """
}

/**
 * Merges mRNA and Guide RNA matrices for a single sample.
 * Aligns guides to the detected mRNA cell barcodes.
 */
process MERGE_MODALITIES {
    tag "${sample_id}"
    publishDir "${params.outdir}/merged_samples", mode: 'copy'
    
    input:
    tuple val(sample_id), path("std_adata.h5ad"), path("kite_adata.h5ad")
    
    output:
    path "${sample_id}_merged.h5ad", emit: h5ad
    path "${sample_id}_guide_diagnostics.json", emit: diagnostics
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scipy.sparse as sp
    import json

    adata_std = ad.read_h5ad("std_adata.h5ad")
    adata_kite = ad.read_h5ad("kite_adata.h5ad")

    def complement_kite_barcode_positions_8_9(barcode):
        barcode = str(barcode)
        if len(barcode) < 9:
            return barcode
        comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return barcode[:7] + barcode[7:9].translate(comp) + barcode[9:]

    def row_sums(matrix):
        if sp.issparse(matrix):
            return np.asarray(matrix.sum(axis=1)).ravel()
        return np.asarray(matrix.sum(axis=1)).ravel()

    def row_nnz(matrix):
        if sp.issparse(matrix):
            return np.diff(matrix.tocsr().indptr)
        return np.count_nonzero(matrix, axis=1)

    raw_kite_obs_names = pd.Index(adata_kite.obs_names.astype(str))
    corrected_kite_obs_names = raw_kite_obs_names.map(complement_kite_barcode_positions_8_9)
    std_obs_names = pd.Index(adata_std.obs_names.astype(str))
    raw_overlap = int(raw_kite_obs_names.isin(std_obs_names).sum())
    corrected_overlap = int(pd.Index(corrected_kite_obs_names).isin(std_obs_names).sum())
    candidates = [
        ("raw", raw_kite_obs_names),
        ("complement_positions_8_9", pd.Index(corrected_kite_obs_names)),
    ]
    candidates = [names for names in candidates if names[1].is_unique]
    if not candidates:
        raise ValueError("No unique KITE barcode representation is available")
    correction, selected_kite_obs_names = max(
        candidates,
        key=lambda item: int(item[1].isin(std_obs_names).sum()),
    )
    selected_overlap = int(selected_kite_obs_names.isin(std_obs_names).sum())
    adata_kite.obs_names = selected_kite_obs_names

    # Align kite to std barcodes
    kite_obs_map = pd.Series(np.arange(adata_kite.n_obs), index=adata_kite.obs_names)
    target_indices = kite_obs_map.reindex(adata_std.obs_names).values
    mask = ~pd.isna(target_indices)
    valid_indices = target_indices[mask].astype(int)
    
    X_found = adata_kite.X[valid_indices, :]
    row_indices = np.where(mask)[0]
    
    if sp.issparse(X_found):
        coo = X_found.tocoo()
        new_row_indices = row_indices[coo.row]
        guides_sparse = sp.csr_matrix(
            (coo.data, (new_row_indices, coo.col)),
            shape=(adata_std.n_obs, adata_kite.n_vars),
            dtype=X_found.dtype
        )
    else:
        guides_sparse = np.zeros((adata_std.n_obs, adata_kite.n_vars), dtype=X_found.dtype)
        guides_sparse[mask, :] = X_found

    adata_std.obsm['guides'] = guides_sparse
    adata_std.uns['guide_names'] = adata_kite.var_names.tolist()

    guide_umis_all_kite = row_sums(adata_kite.X)
    guide_umis_std = row_sums(guides_sparse)
    guide_features_std = row_nnz(guides_sparse)
    guide_positive_umis = guide_umis_std[guide_umis_std > 0]

    diagnostics = {
        "sample_id": "${sample_id}",
        "kite_barcode_correction": correction,
        "n_mrna_cells": int(adata_std.n_obs),
        "n_kite_barcodes": int(adata_kite.n_obs),
        "n_raw_kite_barcodes_overlapping_mrna_barcodes": raw_overlap,
        "n_corrected_kite_barcodes_overlapping_mrna_barcodes": corrected_overlap,
        "n_selected_kite_barcodes_overlapping_mrna_barcodes": selected_overlap,
        "n_barcode_overlap": int(mask.sum()),
        "pct_mrna_barcodes_with_kite_barcode": float(mask.mean()) if adata_std.n_obs else 0.0,
        "total_kite_umis_all_barcodes": float(guide_umis_all_kite.sum()),
        "total_kite_umis_overlapping_mrna_barcodes": float(guide_umis_std.sum()),
        "mrna_cells_with_any_guide_umi": int((guide_umis_std > 0).sum()),
        "mrna_cells_with_at_least_two_nonzero_guides": int((guide_features_std >= 2).sum()),
        "median_guide_umis_in_guide_positive_mrna_cells": float(np.median(guide_positive_umis)) if guide_positive_umis.size else 0.0,
        "max_guide_umis_in_mrna_cells": float(guide_umis_std.max()) if guide_umis_std.size else 0.0,
    }
    print("KITE merge diagnostics:")
    print(json.dumps(diagnostics, indent=2))
    adata_std.uns['kite_merge_diagnostics'] = diagnostics
    
    # Attach sample metadata
    adata_std.obs['sample_id'] = "${sample_id}"

    with open("${sample_id}_guide_diagnostics.json", "w") as handle:
        json.dump(diagnostics, handle, indent=2)

    adata_std.write_h5ad("${sample_id}_merged.h5ad")
    """
}

/**
 * Concatenates all samples into a final unified matrix.
 * Appends sample-specific suffixes to barcodes to prevent collisions.
 * Uses on-disk concatenation so large experiments do not require all samples in RAM.
 */
process CONCATENATE_SAMPLES {
    input:
    path "h5ads/*"
    
    output:
    path "experiment_final_uncompressed.h5ad", emit: h5ad
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    from pathlib import Path

    files = sorted(Path("h5ads").glob("*.h5ad"))
    if not files:
        raise RuntimeError("No sample H5AD files found in h5ads/")

    inputs = {}
    for path in files:
        adata = ad.read_h5ad(path, backed="r")
        try:
            sample_ids = adata.obs["sample_id"].astype(str).unique()
        finally:
            adata.file.close()

        if len(sample_ids) != 1:
            raise RuntimeError(
                f"Expected exactly one sample_id in {path}, found {sample_ids.tolist()}"
            )

        sample_id = sample_ids[0]
        if sample_id in inputs:
            raise RuntimeError(f"Duplicate sample_id across merged H5AD files: {sample_id}")

        inputs[sample_id] = str(path)

    print(f"Concatenating {len(inputs)} samples on disk...")
    ad.experimental.concat_on_disk(
        inputs,
        "experiment_final_uncompressed.h5ad",
        max_loaded_elems=${params.concat_max_loaded_elems},
        axis=0,
        join="outer",
        merge="same",
        uns_merge="same",
        index_unique="-",
    )

    # Workaround: concat_on_disk with uns_merge="same" often results in an empty uns
    # even when all inputs have the same keys. We manually copy uns from the first sample.
    import h5py
    first_sample_path = list(inputs.values())[0]
    with h5py.File(first_sample_path, "r") as f_src:
        if "uns" in f_src:
            print(f"Copying 'uns' from {first_sample_path} to final H5AD...")
            with h5py.File("experiment_final_uncompressed.h5ad", "a") as f_dst:
                if "uns" in f_dst:
                    del f_dst["uns"]
                f_src.copy("uns", f_dst)
    """
}

/**
 * Re-packs the final H5AD with HDF5 gzip compression.
 */
process COMPRESS_FINAL_H5AD {
    publishDir { params.assay_mode == "flex" ? "${params.outdir}/${dataset_id}" : params.outdir }, mode: 'copy'

    input:
    tuple val(dataset_id), path(uncompressed_h5ad)

    output:
    tuple val(dataset_id), path("experiment_final.h5ad"), emit: h5ad

    script:
    """
    h5repack -f ${params.h5repack_filter} ${uncompressed_h5ad} experiment_final.h5ad
    """
}

process QC_COMPARISON {
    tag "${dataset_id}"
    publishDir { params.assay_mode == "flex" ? "${params.outdir}/${dataset_id}" : params.outdir }, mode: 'copy', pattern: '*.filtered.h5ad'
    publishDir { params.comparison_outdir ?: "${params.outdir}/comparison_results" },
        mode: 'copy', pattern: "comparison_results/*",
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(dataset_id), path(reprocessed_h5ad, stageAs: 'reprocessed.h5ad'), path(curated_h5ad, stageAs: 'curated.h5ad')
    path reference_gtf
    path comparison_script

    output:
    tuple val(dataset_id), path('experiment_final.filtered.h5ad'), emit: h5ad
    tuple val(dataset_id), path("comparison_results/${dataset_id}"), emit: reports

    script:
    def cellIdArgs = params.cell_id_columns ? "--cell-id-columns '${params.cell_id_columns}'" : ""
    """
    python ${comparison_script} \
      --dataset-id ${dataset_id} \
      --curated-h5ad ${curated_h5ad} \
      --reprocessed-h5ad ${reprocessed_h5ad} \
      --gtf ${reference_gtf} \
      ${cellIdArgs}
    """
}

process PREPARE_INPUTS {
    tag "${dataset_id}"
    publishDir { (params.assay_mode == "flex" ? "${params.outdir}/${dataset_id}" : params.outdir) + "/dea_gsea/prep" }, mode: "copy"

    input:
    tuple val(dataset_id), path(h5ad)
    val gene_map_path
    path gtf_path

    output:
    tuple val(dataset_id), path("analysis_inputs"), path("analysis_inputs/batches/*.json"), emit: batches
    tuple val(dataset_id), path("analysis_inputs/manifest.json"), emit: manifest

    script:
    def geneMapArg = gene_map_path ? "--gene-map ${gene_map_path}" : ""
    def gtfArg = gtf_path ? "--gtf ${gtf_path}" : ""
    """
    python ${projectDir}/dea-gsea/prepare_inputs.py \
      --h5ad ${h5ad} \
      --outdir analysis_inputs \
      --dataset-id ${dataset_id} \
      --batch-size ${params.batch_size} \
      --min-cells-per-perturbation ${params.min_cells_per_perturbation} \
      --limit-perturbations ${params.limit_perturbations} \
      ${geneMapArg} \
      ${gtfArg}
    """
}


process ANALYZE_BATCH {
    tag "${batch_json.baseName}"
    publishDir { (params.assay_mode == "flex" ? "${params.outdir}/${dataset_id}" : params.outdir) + "/dea_gsea/batch_results" }, mode: "copy"

    input:
    tuple val(dataset_id), path(batch_json), path(analysis_dir), path(h5ad)
    path gmt_path

    output:
    tuple val(dataset_id), path("*.dea.parquet"), path("*.gsea.parquet"), path("*.metrics.json"), emit: results

    script:
    def tieCorrectArg = params.tie_correct ? "--tie-correct" : ""
    def gmtArg = gmt_path ? "--gmt ${gmt_path}" : ""
    """
    python ${projectDir}/dea-gsea/analyze_batch.py \
      --h5ad ${h5ad} \
      --batch-json ${batch_json} \
      --control-indices ${analysis_dir}/control_indices.npy \
      --gene-metadata ${analysis_dir}/gene_metadata.parquet \
      --outdir . \
      --dataset-id ${dataset_id} \
      --matrix-key ${params.matrix_key} \
      --target-sum ${params.target_sum} \
      --threads ${task.cpus} \
      --gsea-permutations ${params.gsea_permutations} \
      --gsea-min-size ${params.gsea_min_size} \
      --gsea-max-size ${params.gsea_max_size} \
      --gsea-seed ${params.gsea_seed} \
      ${gmtArg} \
      ${tieCorrectArg}
    """
}


process MERGE_RESULTS {
    tag "${dataset_id}"
    publishDir { (params.assay_mode == "flex" ? "${params.outdir}/${dataset_id}" : params.outdir) + "/dea_gsea" }, mode: "copy"

    input:
    tuple val(dataset_id), path(dea_files), path(gsea_files), path(metrics_files), path(manifest)

    output:
    tuple val(dataset_id), path("${dataset_id}.dea.parquet"), emit: dea
    tuple val(dataset_id), path("${dataset_id}.gsea.parquet"), emit: gsea
    tuple val(dataset_id), path("${dataset_id}.summary.json"), emit: summary

    script:
    """
    python ${projectDir}/dea-gsea/merge_results.py \
      --dataset-id ${dataset_id} \
      --outdir . \
      --manifest ${manifest} \
      --dea-files ${dea_files.join(" ")} \
      --gsea-files ${gsea_files.join(" ")} \
      --metrics-files ${metrics_files.join(" ")}
    """
}


// =============================================================================
// WORKFLOW
// =============================================================================

workflow {
    if (params.assay_mode == "standard" && (!params.dataset_id || !(params.dataset_id ==~ /[A-Za-z0-9_-]+/)))
        error "Please provide a valid --dataset_id"
    if (!params.gmt || !params.gtf)
        error "Please provide --gmt and --gtf"
    if (!(params.assay_mode in ["standard", "flex"]))
        error "--assay_mode must be standard or flex"

    gtf = file(params.gtf, checkIfExists: true)
    gene_sets = file(params.gmt, checkIfExists: true)

    if (params.assay_mode == "flex") {
        def requiredFlexInputs = [
            "flex_pool_id": params.flex_pool_id,
            "flex_curated_h5ads": params.flex_curated_h5ads,
            "flex_samples": params.flex_samples,
            "flex_gex_sources": params.flex_gex_sources,
            "flex_guide_sources": params.flex_guide_sources,
            "flex_guide_coverage": params.flex_guide_coverage,
            "flex_gene_probes": params.flex_gene_probes,
            "flex_gene_count_features": params.flex_gene_count_features,
            "flex_guide_features": params.flex_guide_features,
            "flex_guide_count_features": params.flex_guide_count_features,
            "flex_guide_targets": params.flex_guide_targets,
            "flex_probe_barcodes": params.flex_probe_barcodes,
            "flex_gex_probe_targets": params.flex_gex_probe_targets,
            "flex_gex_t2g": params.flex_gex_t2g,
            "flex_gex_probe_to_gene": params.flex_gex_probe_to_gene,
            "flex_guide_feature_targets": params.flex_guide_feature_targets,
            "flex_guide_t2g": params.flex_guide_t2g,
            "flex_cbc_whitelist": params.flex_cbc_whitelist,
            "flex_bc_barcode_variants": params.flex_bc_barcode_variants
        ]
        def missingFlexInputs = requiredFlexInputs.findAll { key, value -> !value }
        if (missingFlexInputs)
            error "Missing Flex inputs: ${missingFlexInputs.keySet().join(', ')}"
        if (!(params.flex_pool_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,127}/))
            error "Invalid --flex_pool_id"
        if (!(params.flex_max_forks.toString() ==~ /[1-9][0-9]*/) ||
            params.flex_max_forks.toInteger() > 24)
            error "--flex_max_forks must be between 1 and 24"
        if (params.cell_id_columns != "lane_id")
            error "Flex comparison requires --cell_id_columns lane_id"

        def samples = file(params.flex_samples, checkIfExists: true)
        def gexManifest = file(params.flex_gex_sources, checkIfExists: true)
        def guideManifest = file(params.flex_guide_sources, checkIfExists: true)
        def guideCoverage = file(params.flex_guide_coverage, checkIfExists: true)
        def geneProbes = file(params.flex_gene_probes, checkIfExists: true)
        def geneCountFeatures = file(params.flex_gene_count_features, checkIfExists: true)
        def guideFeatures = file(params.flex_guide_features, checkIfExists: true)
        def guideCountFeatures = file(params.flex_guide_count_features, checkIfExists: true)
        def guideTargets = file(params.flex_guide_targets, checkIfExists: true)
        def probeBarcodes = file(params.flex_probe_barcodes, checkIfExists: true)
        def gexTargets = file(params.flex_gex_probe_targets, checkIfExists: true)
        def gexT2g = file(params.flex_gex_t2g, checkIfExists: true)
        def probeToGene = file(params.flex_gex_probe_to_gene, checkIfExists: true)
        def guideTargetsFasta = file(params.flex_guide_feature_targets, checkIfExists: true)
        def guideT2g = file(params.flex_guide_t2g, checkIfExists: true)
        def cbcWhitelist = file(params.flex_cbc_whitelist, checkIfExists: true)
        def bcVariants = file(params.flex_bc_barcode_variants, checkIfExists: true)
        def flexIndexes = BUILD_INDEX_FLEX(gexTargets, guideTargetsFasta)
        def flexSraBin = params.sra_bin ?: ""
        def flexCode = ["kb_flex_lane.py", "flex_aggregate.py", "bam_to_fastq.py", "stream_count.py", "stream_sra_pairs.py", "read_router.cpp"]
            .collect { name -> file("${projectDir}/bin/${name}", checkIfExists: true) }

        def flexLanes = Channel
            .fromPath(guideCoverage)
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.pool_id == params.flex_pool_id }
            .map { row ->
                if (!(row.lane_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,127}/) || row.gex_source_pairs.toInteger() < 1)
                    error "Invalid Flex lane coverage row: ${row}"
                return [row.pool_id, row.lane_id]
            }

        def flexLaneCounts = KB_FLEX_LANE(
            flexLanes,
            Channel.value(gexManifest), Channel.value(guideManifest),
            Channel.value(guideCoverage), Channel.value(samples),
            Channel.value(probeBarcodes), Channel.value(geneProbes),
            Channel.value(geneCountFeatures), Channel.value(guideFeatures),
            Channel.value(guideCountFeatures), Channel.value(guideTargets),
            params.flex_source_cache ?: "", flexIndexes.gex_index.collect(),
            Channel.value(gexT2g), Channel.value(probeToGene),
            flexIndexes.guide_index.collect(), Channel.value(guideT2g),
            Channel.value(cbcWhitelist), Channel.value(bcVariants), flexSraBin, Channel.value(flexCode)
        )
        def laneFiles = flexLaneCounts.counts
            .map { pool, lane, h5ads, metrics -> h5ads }
            .flatten()
            .collect()
        def poolSamples = Channel.fromPath(samples)
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.pool_id == params.flex_pool_id }
            .map { row ->
                if (!(row.sample_id ==~ /[A-Za-z0-9_-]+/))
                    error "Invalid Flex sample ID: ${row.sample_id}"
                return row.sample_id
            }
            .collect()
            .map { ids ->
                if (!ids || ids.toSet().size() != ids.size())
                    error "Flex pool must declare unique samples"
                return ids
            }
        def aggregationInput = poolSamples.map { ids -> [params.flex_pool_id, ids] }
            .join(laneFiles.map { h5ads -> [params.flex_pool_id, h5ads] })
            .flatMap { pool, ids, h5ads ->
            ids.collect { sid -> [params.flex_pool_id, sid, h5ads] }
        }
        def curatedRows = Channel.fromPath(params.flex_curated_h5ads)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                if (!(row.sample_id ==~ /[A-Za-z0-9_-]+/) || !row.curated_h5ad)
                    error "Invalid Flex comparison manifest row: ${row}"
                return [row.sample_id, file(row.curated_h5ad, checkIfExists: true)]
            }
            .collect(flat: false)
        curated_counts = poolSamples.map { ids -> [params.flex_pool_id, ids] }
            .join(curatedRows.map { rows -> [params.flex_pool_id, rows] })
            .flatMap { pool, ids, rows ->
            if (rows.collect { it[0] }.toSet() != ids.toSet() || rows.size() != ids.size())
                error "Flex comparison manifest must contain exactly one author H5AD per pool sample"
            return rows
        }
        def flexRaw = AGGREGATE_FLEX_SAMPLE(
            aggregationInput,
            Channel.value(samples), Channel.value(guideCoverage),
            Channel.value(gexManifest), Channel.value(guideManifest),
            Channel.value(geneProbes), Channel.value(probeBarcodes),
            Channel.value(geneCountFeatures), Channel.value(guideFeatures),
            Channel.value(guideCountFeatures),
            Channel.value(guideTargets), params.concat_max_loaded_elems, Channel.value(flexCode)
        )
        uncompressed_final = flexRaw.h5ad
    } else {
        if (!params.sample_sheet || !params.transcriptome_fa || !params.features_tsv || !params.curated_h5ad)
            error "Please provide --sample_sheet, --transcriptome_fa, --features_tsv, and --curated_h5ad"
        if (params.split_bam_source &&
            !(params.split_bam_source ==~ /BAM:(SRR|ERR|DRR)[0-9]+/))
            error "--split_bam_source must be BAM:<run accession>"
        if (params.split_bam_source) {
            if (!(params.split_bam_groups.toString() ==~ /[0-9]+/) ||
                params.split_bam_groups.toInteger() < 2 ||
                params.split_bam_groups.toInteger() > 64)
                error "--split_bam_source must be BAM:<accession> and --split_bam_groups must be 2-64"
        } else if (params.split_bam_groups) {
            error "--split_bam_groups requires --split_bam_source"
        }
        def fa = file(params.transcriptome_fa, checkIfExists: true)
        def features = file(params.features_tsv, checkIfExists: true)
        def samplesCh = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header:true, sep:'\t')
            .map { row ->
                def sid = row.sample_id
                def mrna_srrs = row.mRNA_srrs.tokenize(';')
                def sgrna_srrs = row.sgRNA_srrs.tokenize(';')
                def valid_source = { source ->
                    source ==~ /(SRR|ERR|DRR)[0-9]+/ ||
                        (source.startsWith('BAM:') && source.size() > 4) ||
                        (source.startsWith('BAMFILE:') && source.size() > 8)
                }

                if (!(sid ==~ /[A-Za-z0-9_-]+/)) error "Invalid sample_id: ${sid}"
                def feature_offset = (row.guide_feature_offset ?: "0") as Integer
                def feature_length = (row.guide_feature_length ?: "0") as Integer
                if (feature_offset < 0 || feature_length < 0 || feature_length == 1)
                    error "Invalid guide feature trim for sample ${sid}"
                for (runs in [mrna_srrs, sgrna_srrs]) {
                    if (!runs || runs.toSet().size() != runs.size() || runs.any { !valid_source(it) })
                        error "Invalid or duplicate sequencing sources for sample ${sid}"
                }
                return [sid, mrna_srrs, sgrna_srrs, feature_offset, feature_length]
            }
        if (params.limit > 0)
            samplesCh = samplesCh.take(params.limit)

        def stdIdx = BUILD_INDEX_STANDARD(fa, gtf)
        def kiteIdx = BUILD_INDEX_KITE(features)
        if (params.split_bam_source) {
            // Each selected sample maps to exactly one numbered subset of this source.
            def groupSelections = samplesCh.map { sid, mrna, sgrna, offset, length ->
                if (mrna.size() != 1)
                    error "Split BAM samples must have one mRNA source: ${sid}"
                def prefix = params.split_bam_source + "#"
                if (!mrna[0].startsWith(prefix))
                    error "mRNA source for ${sid} must select a group from ${params.split_bam_source}"
                def groupText = mrna[0].substring(prefix.length())
                if (!(groupText ==~ /[1-9][0-9]*/))
                    error "Invalid GEM group selector for ${sid}: ${mrna[0]}"
                def group = groupText as Integer
                if (group > params.split_bam_groups.toInteger())
                    error "GEM group ${group} exceeds --split_bam_groups for ${sid}"
                return [group, sid]
            }.groupTuple().map { group, sample_ids ->
                if (sample_ids.size() != 1)
                    error "More than one sample maps to split BAM group ${group}: ${sample_ids}"
                return [group, sample_ids[0]]
            }
            def splitBams = SPLIT_BAM_GEM_GROUPS(
                Channel.of([params.split_bam_source, params.split_bam_groups.toInteger()])
            ).bams.flatten().map { bam ->
                def match = (bam.name =~ /^group_([1-9][0-9]*)\.bam$/)
                if (!match.matches()) error "Unexpected split BAM name: ${bam.name}"
                return [match[0][1] as Integer, bam]
            }
            def stdCounts = KB_COUNT_STANDARD_BAM(
                groupSelections.join(splitBams).map { group, sid, bam -> [sid, bam] },
                stdIdx.index.collect(), stdIdx.t2g.collect(), params.chemistry
            )
            def kiteCounts = KB_COUNT_KITE(
                samplesCh.map { sid, mrna, sgrna, offset, length -> [sid, sgrna, offset, length] },
                kiteIdx.index.collect(), kiteIdx.t2g.collect(), params.chemistry
            )
            def mergeChannel = stdCounts.h5ad.join(kiteCounts.h5ad)
            def merged = MERGE_MODALITIES(mergeChannel)
            def uncompressed = CONCATENATE_SAMPLES(merged.h5ad.collect())
            uncompressed_final = uncompressed.h5ad.map { h5ad -> [params.dataset_id, h5ad] }
        } else {
            def stdCounts = KB_COUNT_STANDARD(
                samplesCh.map { sid, mrna, sgrna, offset, length -> [sid, mrna] },
                stdIdx.index.collect(), stdIdx.t2g.collect(), params.chemistry
            )
            def kiteCounts = KB_COUNT_KITE(
                samplesCh.map { sid, mrna, sgrna, offset, length -> [sid, sgrna, offset, length] },
                kiteIdx.index.collect(), kiteIdx.t2g.collect(), params.chemistry
            )
            def mergeChannel = stdCounts.h5ad.join(kiteCounts.h5ad)
            def merged = MERGE_MODALITIES(mergeChannel)
            def uncompressed = CONCATENATE_SAMPLES(merged.h5ad.collect())
            uncompressed_final = uncompressed.h5ad.map { h5ad -> [params.dataset_id, h5ad] }
        }
        curated_counts = Channel.value([params.dataset_id, file(params.curated_h5ad, checkIfExists: true)])
    }

    // Step 5: Final HDF5 Compression
    raw_counts = COMPRESS_FINAL_H5AD(uncompressed_final)
    comparison_script = file("${projectDir}/comparison/comparison.py")
    filtered = QC_COMPARISON(raw_counts.h5ad.join(curated_counts), gtf, comparison_script)
    gene_map_path = params.gene_map ? file(params.gene_map, checkIfExists: true).toString() : ""
    prep = PREPARE_INPUTS(filtered.h5ad, gene_map_path, gtf)
    analysis_inputs = prep.batches.join(filtered.h5ad).flatMap { id, directory, batches, h5ad ->
        (batches instanceof List ? batches : [batches]).collect { batch -> [id, batch, directory, h5ad] }
    }
    analyzed = ANALYZE_BATCH(analysis_inputs, gene_sets)
    merged_inputs = analyzed.results.groupTuple()
        .map { id, dea, gsea, metrics -> [id, dea.flatten(), gsea.flatten(), metrics.flatten()] }
        .join(prep.manifest)
    MERGE_RESULTS(merged_inputs)
}
