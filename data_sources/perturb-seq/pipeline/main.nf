#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline Parameters
params.fastq_dir = null
params.metadata_tsv = null
params.sample_sheet = "samples.csv"
params.outdir = "results"
params.chemistry = "10xv3"
params.limit = 0

// Reference parameters
params.transcriptome_fa = null
params.gtf = null
params.features_tsv = null

process PREPARE_SAMPLES {
    executor 'local'
    
    input:
    path metadata_tsv
    path fastq_dir

    output:
    path "samples.csv"

    script:
    """
    python3 ${baseDir}/prepare_samples.py ${metadata_tsv} ${fastq_dir} samples.csv
    """
}

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
    tag "std_${sample_id}"
    publishDir "${params.outdir}/counts_standard/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_filtered/adata.h5ad"), emit: h5ad

    script:
    """
    mkdir -p out
    > batch.txt
    
    # Group FASTQ files by their run prefix
    for fq in ${reads.join(' ')}; do
        base=\$(echo \$fq | sed -E 's/_[0-9]+\\.fastq\\.gz\$//')
        echo "\$base \$fq" >> file_map.txt
    done
    
    # Auto-detect R1/R2 per run
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
            echo -e "${sample_id}\t\$r1\t\$r2" >> batch.txt
        fi
    done

    kb count -i ${index} -g ${t2g} -x ${chemistry} -o out --h5ad --filter bustools -t ${task.cpus} batch.txt
    """
}

process KB_COUNT_KITE {
    tag "kite_${sample_id}"
    publishDir "${params.outdir}/counts_kite/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_filtered/adata.h5ad"), emit: h5ad

    script:
    """
    mkdir -p out
    > batch.txt
    
    for fq in ${reads.join(' ')}; do
        base=\$(echo \$fq | sed -E 's/_[0-9]+\\.fastq\\.gz\$//')
        echo "\$base \$fq" >> file_map.txt
    done
    
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
            echo -e "${sample_id}\t\$r1\t\$r2" >> batch.txt
        fi
    done

    kb count -i ${index} -g ${t2g} -x ${chemistry} -o out --workflow kite --h5ad --filter bustools -t ${task.cpus} batch.txt
    """
}

process MERGE_MODALITIES {
    tag "${sample_id}"
    publishDir "${params.outdir}/merged_samples", mode: 'copy'
    
    input:
    tuple val(sample_id), path("std_adata.h5ad"), path("kite_adata.h5ad")
    
    output:
    path "${sample_id}_merged.h5ad", emit: h5ad
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scipy.sparse as sp

    adata_std = ad.read_h5ad("std_adata.h5ad")
    adata_kite = ad.read_h5ad("kite_adata.h5ad")

    # Align kite to std
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
    
    # Store sample_id in obs
    adata_std.obs['sample_id'] = "${sample_id}"

    adata_std.write_h5ad("${sample_id}_merged.h5ad")
    """
}

process CONCATENATE_SAMPLES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path "h5ads/*"
    
    output:
    path "experiment_final.h5ad"
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    import os
    import glob

    files = sorted(glob.glob("h5ads/*.h5ad"))
    adatas = []
    for f in files:
        a = ad.read_h5ad(f)
        # Suffix barcodes with sample_id to prevent collisions
        sample_id = a.obs['sample_id'].iloc[0]
        a.obs_names = a.obs_names + "-" + str(sample_id)
        adatas.append(a)

    print(f"Concatenating {len(adatas)} samples...")
    # join='outer' to ensure we keep all genes, but they should be aligned already
    merged = ad.concat(adatas, join='outer', index_unique=None, merge='same')
    
    # Re-verify guides consistency
    merged.uns['guide_names'] = adatas[0].uns['guide_names']
    
    merged.write_h5ad("experiment_final.h5ad", compression="gzip")
    """
}

workflow {
    if (!params.fastq_dir || !params.transcriptome_fa || !params.gtf || !params.features_tsv || !params.metadata_tsv) {
        error "Please provide --fastq_dir, --metadata_tsv, --transcriptome_fa, --gtf, and --features_tsv"
    }
    
    fa = file(params.transcriptome_fa)
    gtf = file(params.gtf)
    features = file(params.features_tsv)
    metadata = file(params.metadata_tsv)
    fastq_dir = file(params.fastq_dir)

    sample_sheet = PREPARE_SAMPLES(metadata, fastq_dir)
    
    samples_ch = sample_sheet
        .splitCsv(header:true)
        .map { row -> 
            def sid = row.sample_id
            def mrna_srrs = row.mRNA_srrs.split(';')
            def sgrna_srrs = row.sgRNA_srrs.split(';')
            
            def mrna_files = mrna_srrs.collect { srr -> file("${params.fastq_dir}/${srr}_{1,2,3}.fastq.gz") }.flatten()
            def sgrna_files = sgrna_srrs.collect { srr -> file("${params.fastq_dir}/${srr}_{1,2,3}.fastq.gz") }.flatten()
            
            return [sid, mrna_files, sgrna_files]
        }

    if (params.limit > 0) {
        samples_ch = samples_ch.take(params.limit)
    }

    std_idx = BUILD_INDEX_STANDARD(fa, gtf)
    kite_idx = BUILD_INDEX_KITE(features)

    // Parallel processing per sample
    std_counts = KB_COUNT_STANDARD(samples_ch.map { it[0], it[1] }, std_idx.index, std_idx.t2g, params.chemistry)
    kite_counts = KB_COUNT_KITE(samples_ch.map { it[0], it[2] }, kite_idx.index, kite_idx.t2g, params.chemistry)

    // Join cDNA and Guide results by sample_id
    merge_ch = std_counts.h5ad.join(kite_counts.h5ad)
    
    merged_samples = MERGE_MODALITIES(merge_ch)
    
    CONCATENATE_SAMPLES(merged_samples.h5ad.collect())
}
