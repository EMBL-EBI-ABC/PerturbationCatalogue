#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline Parameters
params.fastq_dir = null
params.outdir = "results"
params.workflow = "standard" // "standard" for Gene Expression, "kite" for Feature Barcodes/CRISPR
params.chemistry = "10x_v3"  // e.g. 10x_v2, 10x_v3

// Reference parameters for standard workflow
params.transcriptome_fa = null
params.gtf = null

// Reference parameters for KITE workflow
params.features_tsv = null

// FASTQ read selection. Often SRA yields 3 files (e.g., _1=I1, _2=R1, _3=R2). 
// You can customize the read pattern here to select the actual R1 and R2.
params.reads_pattern = "*_{1,2}.fastq.gz" // Default assumes paired-end _1 (R1) and _2 (R2)

process BUILD_INDEX_STANDARD {
    tag "cDNA_index"
    publishDir "${params.outdir}/reference", mode: 'copy'

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

process BUILD_INDEX_KITE {
    tag "KITE_index"
    publishDir "${params.outdir}/reference", mode: 'copy'

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

process KB_COUNT {
    tag "${sample_id}"
    publishDir "${params.outdir}/counts/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path t2g
    val workflow
    val chemistry

    output:
    path "out/counts_unfiltered/*.h5ad", emit: h5ad, optional: true
    path "out/**", emit: all_outputs

    script:
    def wf_flag = workflow == 'kite' ? '--workflow kite' : ''
    """
    mkdir -p out
    kb count -i ${index} \\
             -g ${t2g} \\
             -x ${chemistry} \\
             -o out \\
             ${wf_flag} \\
             --h5ad \\
             -t ${task.cpus} \\
             ${reads.join(' ')}
    """
}

workflow {
    if (!params.fastq_dir) {
        error "Please provide --fastq_dir"
    }

    read_pairs_ch = Channel.fromFilePairs("${params.fastq_dir}/${params.reads_pattern}", size: -1)

    if (params.workflow == "standard") {
        if (!params.transcriptome_fa || !params.gtf) {
            error "Standard workflow requires --transcriptome_fa and --gtf"
        }
        fa = file(params.transcriptome_fa, checkIfExists: true)
        gtf = file(params.gtf, checkIfExists: true)
        
        index_out = BUILD_INDEX_STANDARD(fa, gtf)
        KB_COUNT(read_pairs_ch, index_out.index, index_out.t2g, params.workflow, params.chemistry)
        
    } else if (params.workflow == "kite") {
        if (!params.features_tsv) {
            error "KITE workflow requires --features_tsv (a tab-separated file mapping feature IDs to sequences)"
        }
        features = file(params.features_tsv, checkIfExists: true)
        
        index_out = BUILD_INDEX_KITE(features)
        KB_COUNT(read_pairs_ch, index_out.index, index_out.t2g, params.workflow, params.chemistry)
        
    } else {
        error "Unknown workflow: ${params.workflow}. Use 'standard' or 'kite'."
    }
}
