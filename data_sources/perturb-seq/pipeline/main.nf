#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline Parameters
params.fastq_dir = null
params.outdir = "results"
params.chemistry = "10x_v3"  // e.g. 10x_v2, 10x_v3

// Reference parameters for standard workflow (Gene Expression)
params.transcriptome_fa = null
params.gtf = null

// Reference parameters for KITE workflow (Guides/CRISPR)
params.features_tsv = null

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

process KB_COUNT_STANDARD {
    tag "${sample_id}"
    publishDir "${params.outdir}/counts_standard/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_unfiltered/adata.h5ad"), emit: h5ad, optional: true
    path "out/**", emit: all_outputs

    script:
    """
    mkdir -p out

    # Auto-detect R1 (barcode+UMI) and R2 (biological read) based on sequence length
    R1=""
    R2=""
    for fq in ${reads.join(' ')}; do
        # Extract the length of the first sequence read
        seq_len=\$(zcat \$fq | head -n 2 | tail -n 1 | tr -d '\\n' | wc -c)
        if [ "\$seq_len" -ge 20 ] && [ "\$seq_len" -le 40 ]; then
            R1=\$fq
        elif [ "\$seq_len" -gt 40 ]; then
            R2=\$fq
        fi
    done

    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
        echo "Error: Could not auto-detect R1 and R2 from fastq lengths for ${sample_id}."
        exit 1
    fi

    echo "Auto-detected R1: \$R1"
    echo "Auto-detected R2: \$R2"

    kb count -i ${index} \\
             -g ${t2g} \\
             -x ${chemistry} \\
             -o out \\
             --h5ad \\
             -t ${task.cpus} \\
             \$R1 \$R2
    """
}

process KB_COUNT_KITE {
    tag "${sample_id}"
    publishDir "${params.outdir}/counts_kite/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_unfiltered/adata.h5ad"), emit: h5ad, optional: true
    path "out/**", emit: all_outputs

    script:
    """
    mkdir -p out

    # Auto-detect R1 (barcode+UMI) and R2 (biological read) based on sequence length
    R1=""
    R2=""
    for fq in ${reads.join(' ')}; do
        # Extract the length of the first sequence read
        seq_len=\$(zcat \$fq | head -n 2 | tail -n 1 | tr -d '\\n' | wc -c)
        if [ "\$seq_len" -ge 20 ] && [ "\$seq_len" -le 40 ]; then
            R1=\$fq
        elif [ "\$seq_len" -gt 40 ]; then
            R2=\$fq
        fi
    done

    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
        echo "Error: Could not auto-detect R1 and R2 from fastq lengths for ${sample_id}."
        exit 1
    fi

    echo "Auto-detected R1: \$R1"
    echo "Auto-detected R2: \$R2"

    kb count -i ${index} \\
             -g ${t2g} \\
             -x ${chemistry} \\
             -o out \\
             --workflow kite \\
             --h5ad \\
             -t ${task.cpus} \\
             \$R1 \$R2
    """
}

process MERGE_SRR {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(std_h5ad), path(kite_h5ad)
    
    output:
    path "${sample_id}_merged.h5ad"
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    import sys
    
    srr_id = "${sample_id}"
    adata_std = ad.read_h5ad("${std_h5ad}")
    adata_kite = ad.read_h5ad("${kite_h5ad}")
    
    kite_df = adata_kite.to_df()
    # Align the kite (guide) matrix rows to the standard (cDNA) cell barcodes
    kite_df_aligned = kite_df.reindex(adata_std.obs_names, fill_value=0)
    
    # Store guide counts as a multi-dimensional array in .obsm
    adata_std.obsm['guides'] = kite_df_aligned.values
    # Save the guide names corresponding to the columns of the array
    adata_std.uns['guide_names'] = kite_df_aligned.columns.tolist()
    
    adata_std.write_h5ad(f"{srr_id}_merged.h5ad")
    """
}

process CONCAT_ALL {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path h5ad_files
    
    output:
    path "experiment_final.h5ad"
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    import glob
    import sys

    files = glob.glob("*_merged.h5ad")
    if not files:
        sys.exit("No merged h5ad files found.")
        
    keys = [f.replace("_merged.h5ad", "") for f in files]

    # Use ad.concat to combine all SRRs into one massive dataset
    # index_unique="-" will append the SRR id to the cell barcodes (e.g. AAACCC...-SRR123)
    adata = ad.concat(
        {k: ad.read_h5ad(f) for k, f in zip(keys, files)},
        label="SRR_run",
        index_unique="-"
    )

    adata.write_h5ad("experiment_final.h5ad")
    """
}

workflow {
    if (!params.fastq_dir || !params.transcriptome_fa || !params.gtf || !params.features_tsv) {
        error "Please provide --fastq_dir, --transcriptome_fa, --gtf, and --features_tsv"
    }
    
    fa = file(params.transcriptome_fa, checkIfExists: true)
    gtf = file(params.gtf, checkIfExists: true)
    features = file(params.features_tsv, checkIfExists: true)

    // Match all relevant fastq files for the given SRR prefixes (_1, _2, _3, etc.)
    read_pairs_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2,3}.fastq.gz", size: -1)

    std_idx = BUILD_INDEX_STANDARD(fa, gtf)
    kite_idx = BUILD_INDEX_KITE(features)

    std_counts = KB_COUNT_STANDARD(read_pairs_ch, std_idx.index, std_idx.t2g, params.chemistry)
    kite_counts = KB_COUNT_KITE(read_pairs_ch, kite_idx.index, kite_idx.t2g, params.chemistry)

    // Join the resulting h5ad files exactly by SRR ID
    joined_ch = std_counts.h5ad.join(kite_counts.h5ad)
    
    // Merge standard and kite matrices per-SRR
    merged_ch = MERGE_SRR(joined_ch)
    
    // Concatenate all merged SRRs into the final H5AD
    CONCAT_ALL(merged_ch.collect())
}
