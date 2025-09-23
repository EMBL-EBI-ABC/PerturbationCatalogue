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

  // Get batch metadata from preprocessing
  batch_metadata = PREPROCESS(h5ads, params.sample_id)

  // Create a separate channel for h5ad files to use in PSEUDOBULK
  h5ads_for_pseudobulk = Channel.fromPath(params.input)
                                .ifEmpty { error "No .h5ad matched: ${params.input}" }

  // Combine h5ad files with their corresponding batch metadata
  combined_data = h5ads_for_pseudobulk
    .map { h5ad -> tuple(file(h5ad).name, h5ad) }
    .cross(
      batch_metadata.map { h5ad_staged, tsv ->
        tuple(file(h5ad_staged).name, tsv)
      }
    ) { it[0] }
    .map { h5ad_tuple, metadata_tuple ->
      def h5ad_file = h5ad_tuple[1]
      def tsv_file = metadata_tuple[1]
      tuple(h5ad_file, tsv_file)
    }

  // Expand into batch tuples
  batch_tuples = combined_data
    .flatMap { h5ad, tsv ->
        def lines = file(tsv).text.readLines()
        def header = lines[0].tokenize('\t')
        def idx_id = header.indexOf('batch_id')
        def idx_ctl = header.indexOf('control')
        def idx_per = header.indexOf('perts')

        lines.drop(1).collect { ln ->
            def f = ln.split(/\t/)
            tuple(h5ad, f[idx_id], f[idx_ctl], f[idx_per])
        }
    }

  pb = PSEUDOBULK(batch_tuples, params.sample_id)
  de = DE_PYDESEQ2(params.sample_id, pb.pb_counts, pb.pb_meta, pb.batch_id)

  emit:
    de.csv_files
    de.summary
}

process PREPROCESS {
  tag { file(h5ad).baseName }
  label 'normal'
  publishDir "${params.outdir}/meta", mode: 'copy'

  input:
    path h5ad
    val sample_id

  output:
    tuple path(h5ad), path("${sample_id}.pert_batches.tsv")

  script:
  """
  # Make a single TSV describing the batches (control + up to 999 perts)
  scp_preprocess.py \
    --input ${h5ad} \
    --others-per-file 99 \
    --out ${params.sample_id}.pert_batches.tsv
  """
}

process PSEUDOBULK {
  tag "pseudobulk_${batch_id}"
  label 'normal'
  publishDir "${params.outdir}/pseudobulk", mode: 'copy'

  input:
    tuple path(h5ad), val(batch_id), val(control), val(perts_csv)
    val sample_id

  output:
    path "${sample_id}_${batch_id}.pb_counts.parquet", emit: pb_counts
    path "${sample_id}_${batch_id}.pb_meta.parquet", emit: pb_meta
    path "${sample_id}_${batch_id}.pb_summary.tsv", emit: pb_summary
    val(batch_id), emit: batch_id

  script:
  """
  scp_pseudobulk.py \
    --input ${h5ad} \
    --pert-list '${perts_csv}' \
    --min-reps ${params.min_reps} \
    --replicate-mode ${params.replicate_mode} \
    --auto_min_cells_per_bag ${params.auto_min_cells_per_bag} \
    --auto_min_bags ${params.auto_min_bags} \
    --auto_max_bags ${params.auto_max_bags} \
    --random-seed ${params.random_seed} \
    --out-counts ${params.sample_id}_${batch_id}.pb_counts.parquet \
    --out-meta ${params.sample_id}_${batch_id}.pb_meta.parquet \
    --out-summary ${params.sample_id}_${batch_id}.pb_summary.tsv
  """
}

process DE_PYDESEQ2 {
  tag "DESeq2like"
  label 'big'
  publishDir "${params.outdir}/de_pseudobulk", mode: 'copy'

  input:
    val sample_id
    path counts
    path meta
    val(batch_id)

  output:
    path "${sample_id}.de_*.csv", emit: csv_files
    path "${sample_id}_${batch_id}.summary.tsv", emit: summary

  script:
  """
  scp_de_pydeseq2.py \
    --counts ${counts} \
    --meta ${meta} \
    --covariates ${params.covariates} \
    --out-prefix ${params.sample_id}.de_ \
    --summary ${params.sample_id}_${batch_id}.summary.tsv
  """
}