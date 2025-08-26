nextflow.enable.dsl=2

params.input           = params.input ?: "data/*.h5ad"
params.outdir          = params.outdir ?: "results"
params.counts_layer    = params.counts_layer ?: "counts"
params.min_reps        = params.min_reps ?: 3
params.covariates      = params.covariates ?: "cell_type_label"
params.sample_id      = params.sample_id ?: "undefined"

// auto‑bagging knobs
params.replicate_mode          = params.replicate_mode ?: "auto"
params.auto_min_cells_per_bag  = params.auto_min_cells_per_bag ?: 20
params.auto_min_bags           = params.auto_min_bags ?: 5
params.auto_max_bags           = params.auto_max_bags ?: 10
params.random_seed             = params.random_seed ?: 1

workflow {

  h5ads = Channel.fromPath(params.input)
                 .ifEmpty { error "No .h5ad matched: ${params.input}" }
  prepped = PREPROCESS(params.sample_id, h5ads)
  pb = PSEUDOBULK(params.sample_id, prepped)
  de = DE_PYDESEQ2(params.sample_id, pb.pb_counts, pb.pb_meta)
  emit:
    de.parquet_files
    de.summary
}

process PREPROCESS {
  label 'bigmem'
  tag "preprocess"
  publishDir "${params.outdir}/preprocessed", mode: 'copy'
  input:
    val sample_id
    path h5ad
  output:
    path "${sample_id}.preprocessed.h5ad"
  script:
  """
  scp_preprocess.py \
    --input ${h5ad} \
    --counts-layer ${params.counts_layer} \
    --out ${params.sample_id}.preprocessed.h5ad
  """
}

process PSEUDOBULK {
  label 'bigmem'
  tag "pseudobulk"
  publishDir "${params.outdir}/pseudobulk", mode: 'copy'
  input:
    val sample_id
    path ad
  output:
    path "${sample_id}.pb_counts.parquet", emit: pb_counts
    path "${sample_id}.pb_meta.parquet", emit: pb_meta
    path "${sample_id}.pb_summary.tsv", emit: pb_summary
  script:
  """
  scp_pseudobulk.py \
    --input ${ad} \
    --min-reps ${params.min_reps} \
    --replicate-mode ${params.replicate_mode} \
    --auto_min_cells_per_bag ${params.auto_min_cells_per_bag} \
    --auto_min_bags ${params.auto_min_bags} \
    --auto_max_bags ${params.auto_max_bags} \
    --random-seed ${params.random_seed} \
    --out-counts ${params.sample_id}.pb_counts.parquet \
    --out-meta ${params.sample_id}.pb_meta.parquet \
    --out-summary ${params.sample_id}.pb_summary.tsv
  """
}

process DE_PYDESEQ2 {
  label 'bigmemcpu'
  tag "DESeq2like"
  publishDir "${params.outdir}/de_pseudobulk", mode: 'copy'
  input:
    val sample_id
    path counts
    path meta
  output:
    path "${sample_id}.de_*.parquet", emit: parquet_files
    path "${sample_id}.summary.tsv", emit: summary
  script:
  """
  scp_de_pydeseq2.py \
    --counts ${counts} \
    --meta ${meta} \
    --covariates ${params.covariates} \
    --out-prefix ${params.sample_id}.de_ \
    --summary ${params.sample_id}.summary.tsv
  """
}