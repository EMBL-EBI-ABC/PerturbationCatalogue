#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline Parameters
params.fastq_dir = null
params.outdir = "results"
params.chemistry = "10xv3"  // e.g. 10xv2, 10xv3
params.limit = 0 // Limit number of FASTQs processed (for debugging). 0 = no limit.

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
    publishDir "${params.outdir}/counts_standard", mode: 'copy'

    input:
    path reads
    path index
    path t2g
    val chemistry

    output:
    path "out/counts_filtered/adata.h5ad", emit: h5ad, optional: true
    path "out/**", emit: all_outputs

    script:
    """
    mkdir -p out

    # Create batch file for kb count (Format: sample_id \\t R1 \\t R2)
    > batch.txt
    
    # Group FASTQ files by their run prefix by stripping the _[0-9].fastq.gz suffix
    for fq in ${reads.join(' ')}; do
        base=\$(echo \$fq | sed -E 's/_[0-9]+\\.fastq\\.gz\$//')
        echo "\$base \$fq" >> file_map.txt
    done
    
    # For each run, auto-detect R1 (barcode+UMI) and R2 (biological read) based on sequence length
    awk '{print \$1}' file_map.txt | sort | uniq | while read base; do
        r1=""
        r2=""
        for fq in \$(grep "^\$base " file_map.txt | awk '{print \$2}'); do
            seq_len=\$(zcat \$fq | head -n 2 | tail -n 1 | tr -d '\\n' | wc -c)
            if [ "\$seq_len" -ge 20 ] && [ "\$seq_len" -le 40 ]; then
                r1=\$fq
            elif [ "\$seq_len" -gt 40 ]; then
                r2=\$fq
            fi
        done
        
        if [ -n "\$r1" ] && [ -n "\$r2" ]; then
            echo -e "ALL_CELLS\t\$r1\t\$r2" >> batch.txt
        else
            echo "Warning: Could not pair R1 and R2 for run \$base" >&2
        fi
    done

    echo "Batch file generated:"
    cat batch.txt

    kb count -i ${index} \\
             -g ${t2g} \\
             -x ${chemistry} \\
             -o out \\
             --h5ad \\
             --filter bustools \\
             -t ${task.cpus} \\
             batch.txt
    """
}

process KB_COUNT_KITE {
    publishDir "${params.outdir}/counts_kite", mode: 'copy'

    input:
    path reads
    path index
    path t2g
    val chemistry

    output:
    path "out/counts_filtered/adata.h5ad", emit: h5ad, optional: true
    path "out/**", emit: all_outputs

    script:
    """
    mkdir -p out

    # Create batch file for kb count (Format: sample_id \\t R1 \\t R2)
    > batch.txt
    
    # Group FASTQ files by their run prefix by stripping the _[0-9].fastq.gz suffix
    for fq in ${reads.join(' ')}; do
        base=\$(echo \$fq | sed -E 's/_[0-9]+\\.fastq\\.gz\$//')
        echo "\$base \$fq" >> file_map.txt
    done
    
    # For each run, auto-detect R1 (barcode+UMI) and R2 (biological read) based on sequence length
    awk '{print \$1}' file_map.txt | sort | uniq | while read base; do
        r1=""
        r2=""
        for fq in \$(grep "^\$base " file_map.txt | awk '{print \$2}'); do
            seq_len=\$(zcat \$fq | head -n 2 | tail -n 1 | tr -d '\\n' | wc -c)
            if [ "\$seq_len" -ge 20 ] && [ "\$seq_len" -le 40 ]; then
                r1=\$fq
            elif [ "\$seq_len" -gt 40 ]; then
                r2=\$fq
            fi
        done
        
        if [ -n "\$r1" ] && [ -n "\$r2" ]; then
            echo -e "ALL_CELLS\t\$r1\t\$r2" >> batch.txt
        else
            echo "Warning: Could not pair R1 and R2 for run \$base" >&2
        fi
    done

    echo "Batch file generated:"
    cat batch.txt

    kb count -i ${index} \\
             -g ${t2g} \\
             -x ${chemistry} \\
             -o out \\
             --workflow kite \\
             --h5ad \\
             --filter bustools \\
             -t ${task.cpus} \\
             batch.txt
    """
}

process MERGE_MODALITIES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path std_h5ad
    path kite_h5ad
    
    output:
    path "experiment_final.h5ad"
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scipy.sparse as sp

    print("Loading standard expression matrix...")
    adata_std = ad.read_h5ad("${std_h5ad}")
    
    print("Loading KITE guides matrix...")
    adata_kite = ad.read_h5ad("${kite_h5ad}")

    # Align the kite (guide) matrix rows to the standard (cDNA) cell barcodes
    kite_obs_map = pd.Series(np.arange(adata_kite.n_obs), index=adata_kite.obs_names)
    target_indices = kite_obs_map.reindex(adata_std.obs_names).values
    
    mask = ~pd.isna(target_indices)
    valid_indices = target_indices[mask].astype(int)
    
    print(f"Matched {mask.sum()} out of {adata_std.n_obs} cells with guide counts.")
    
    # Construct the aligned sparse matrix for guides
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

    print("Saving final combined matrix...")
    adata_std.write_h5ad("experiment_final.h5ad")
    print("Done!")
    """
}

workflow {
    if (!params.fastq_dir || !params.transcriptome_fa || !params.gtf || !params.features_tsv) {
        error "Please provide --fastq_dir, --transcriptome_fa, --gtf, and --features_tsv"
    }
    
    fa = file(params.transcriptome_fa, checkIfExists: true)
    gtf = file(params.gtf, checkIfExists: true)
    features = file(params.features_tsv, checkIfExists: true)

    // Collect all FASTQ files
    fastq_files = Channel.fromPath("${params.fastq_dir}/*_{1,2,3}.fastq.gz")
    if (params.limit > 0) {
        fastq_files = fastq_files.take(params.limit)
    }
    reads_ch = fastq_files.collect()

    std_idx = BUILD_INDEX_STANDARD(fa, gtf)
    kite_idx = BUILD_INDEX_KITE(features)

    std_counts = KB_COUNT_STANDARD(reads_ch, std_idx.index, std_idx.t2g, params.chemistry)
    kite_counts = KB_COUNT_KITE(reads_ch, kite_idx.index, kite_idx.t2g, params.chemistry)

    MERGE_MODALITIES(std_counts.h5ad, kite_counts.h5ad)
}
