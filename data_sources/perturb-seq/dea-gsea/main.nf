#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.dataset_id = "nadig_2025_jurkat"
params.h5ad = null
params.gmt = null
params.outdir = "results"
params.batch_size = 50
params.min_cells_per_perturbation = 10
params.limit_perturbations = 0
params.target_sum = 10000
params.matrix_key = "X"
params.gene_map = ""
params.gtf = ""
params.gsea_permutations = 1000
params.gsea_min_size = 15
params.gsea_max_size = 500
params.gsea_seed = 1
params.tie_correct = false
params.target_col = "perturbed_target_symbol"
params.gene_count_col = "called_knockout_gene_count"
params.control_probe_count_col = "called_control_probe_count"
params.call_type_col = "perturbation_call_type"


process PREPARE_INPUTS {
    tag "${params.dataset_id}"
    publishDir "${params.outdir}/prep", mode: "copy"

    input:
    path h5ad
    val gene_map_path
    val gtf_path

    output:
    path "analysis_inputs", emit: analysis_dir
    path "analysis_inputs/batches/*.json", emit: batches
    path "analysis_inputs/manifest.json", emit: manifest

    script:
    def geneMapArg = gene_map_path ? "--gene-map ${gene_map_path}" : ""
    def gtfArg = gtf_path ? "--gtf ${gtf_path}" : ""
    """
    python ${baseDir}/bin/prepare_inputs.py \
      --h5ad ${h5ad} \
      --outdir analysis_inputs \
      --dataset-id ${params.dataset_id} \
      --target-col ${params.target_col} \
      --gene-count-col ${params.gene_count_col} \
      --control-probe-count-col ${params.control_probe_count_col} \
      --call-type-col ${params.call_type_col} \
      --batch-size ${params.batch_size} \
      --min-cells-per-perturbation ${params.min_cells_per_perturbation} \
      --limit-perturbations ${params.limit_perturbations} \
      ${geneMapArg} \
      ${gtfArg}
    """
}


process ANALYZE_BATCH {
    tag "${batch_json.baseName}"
    publishDir "${params.outdir}/batch_results", mode: "copy"

    input:
    tuple path(batch_json), path(analysis_dir)
    path h5ad
    val gmt_path

    output:
    path "*.dea.parquet", emit: dea
    path "*.gsea.parquet", emit: gsea
    path "*.metrics.json", emit: metrics

    script:
    def tieCorrectArg = params.tie_correct ? "--tie-correct" : ""
    def gmtArg = gmt_path ? "--gmt ${gmt_path}" : ""
    """
    python ${baseDir}/bin/analyze_batch.py \
      --h5ad ${h5ad} \
      --batch-json ${batch_json} \
      --control-indices ${analysis_dir}/control_indices.npy \
      --gene-metadata ${analysis_dir}/gene_metadata.parquet \
      --outdir . \
      --dataset-id ${params.dataset_id} \
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
    publishDir "${params.outdir}", mode: "copy"

    input:
    path dea_files
    path gsea_files
    path metrics_files
    path manifest
    val dataset_id

    output:
    path "${dataset_id}.dea.parquet", emit: dea
    path "${dataset_id}.gsea.parquet", emit: gsea
    path "${dataset_id}.summary.json", emit: summary

    script:
    """
    python ${baseDir}/bin/merge_results.py \
      --dataset-id ${dataset_id} \
      --outdir . \
      --manifest ${manifest} \
      --dea-files ${dea_files.join(" ")} \
      --gsea-files ${gsea_files.join(" ")} \
      --metrics-files ${metrics_files.join(" ")}
    """
}


workflow {
    if (!params.h5ad) {
        error "Please provide --h5ad"
    }
    if (!params.gmt) {
        error "Please provide --gmt"
    }

    gene_map_path = params.gene_map ? file(params.gene_map).toString() : ""
    gtf_path = params.gtf ? file(params.gtf).toString() : ""
    gmt_path = params.gmt ? file(params.gmt).toString() : ""

    prep = PREPARE_INPUTS(file(params.h5ad), gene_map_path, gtf_path)

    analysis_inputs = prep.batches.flatten().combine(prep.analysis_dir)
    analyzed = ANALYZE_BATCH(analysis_inputs, file(params.h5ad), gmt_path)

    MERGE_RESULTS(
        analyzed.dea.collect(),
        analyzed.gsea.collect(),
        analyzed.metrics.collect(),
        prep.manifest,
        params.dataset_id
    )
}
