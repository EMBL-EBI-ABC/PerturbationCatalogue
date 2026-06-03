#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// =============================================================================
// PIPELINE PARAMETERS
// =============================================================================
params.fastq_dir = null
params.sample_sheet = null
params.outdir = "results"
params.chemistry = "10xv3"
params.limit = 0
params.concat_max_loaded_elems = 100000000
params.h5repack_filter = "GZIP=4"

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

/**
 * Standard cDNA quantification per sample using Kallisto-Bustools.
 * Groups all lanes (L001-L004) for a single physical well.
 */
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
    > file_map.txt
    
    # Map SRRs to their constituent FASTQ files
    for fq in ${reads.join(' ')}; do
        base=\$(echo \$fq | sed -E 's/_[0-9]+\\.fastq\\.gz\$//')
        echo "\$base \$fq" >> file_map.txt
    done
    
    # Pair R1 (Barcode/UMI) and R2 (Transcript) based on read length
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

/**
 * KITE Guide RNA quantification per sample.
 */
process KB_COUNT_KITE {
    tag "kite_${sample_id}"
    publishDir "${params.outdir}/counts_kite/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path t2g
    val chemistry

    output:
    tuple val(sample_id), path("out/counts_unfiltered/adata.h5ad"), emit: h5ad

    script:
    """
    mkdir -p out
    > file_map.txt
    > paired_fastqs.tsv
    > fastq_args.txt
    
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
            echo -e "\$r1\t\$r2" >> paired_fastqs.tsv
            echo "\$r1" >> fastq_args.txt
            echo "\$r2" >> fastq_args.txt
        fi
    done

    if [ ! -s paired_fastqs.tsv ]; then
        echo "No valid R1/R2 sgRNA FASTQ pairs found for sample ${sample_id}" >&2
        exit 1
    fi

    echo "KITE sample ${sample_id}: using \$(wc -l < paired_fastqs.tsv) R1/R2 FASTQ pairs"
    fastq_args=\$(tr '\\n' ' ' < fastq_args.txt)

    kb count -i ${index} -g ${t2g} -x ${chemistry} -o out --workflow kite --h5ad -t ${task.cpus} \$fastq_args
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

    # The deposited sgRNA FASTQs encode barcode bases 8-9 as complements relative to mRNA.
    raw_kite_obs_names = pd.Index(adata_kite.obs_names.astype(str))
    corrected_kite_obs_names = raw_kite_obs_names.map(complement_kite_barcode_positions_8_9)
    std_obs_names = pd.Index(adata_std.obs_names.astype(str))
    raw_overlap = int(raw_kite_obs_names.isin(std_obs_names).sum())
    corrected_overlap = int(pd.Index(corrected_kite_obs_names).isin(std_obs_names).sum())
    adata_kite.obs_names = corrected_kite_obs_names

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
        "kite_barcode_correction": "complement_positions_8_9",
        "n_mrna_cells": int(adata_std.n_obs),
        "n_kite_barcodes": int(adata_kite.n_obs),
        "n_raw_kite_barcodes_overlapping_mrna_barcodes": raw_overlap,
        "n_corrected_kite_barcodes_overlapping_mrna_barcodes": corrected_overlap,
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
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path uncompressed_h5ad

    output:
    path "experiment_final.h5ad"

    script:
    """
    h5repack -f ${params.h5repack_filter} ${uncompressed_h5ad} experiment_final.h5ad
    """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {
    if (!params.fastq_dir || !params.sample_sheet || !params.transcriptome_fa || !params.gtf || !params.features_tsv) {
        error "Please provide --fastq_dir, --sample_sheet, --transcriptome_fa, --gtf, and --features_tsv"
    }
    
    fa = file(params.transcriptome_fa)
    gtf = file(params.gtf)
    features = file(params.features_tsv)
    
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header:true, sep:'\t')
        .map { row -> 
            def sid = row.sample_id
            def mrna_srrs = row.mRNA_srrs.tokenize(';')
            def sgrna_srrs = row.sgRNA_srrs.tokenize(';')
            
            def mrna_files = mrna_srrs.collect { srr -> file("${params.fastq_dir}/${srr}_{1,2,3}.fastq.gz") }.flatten()
            def sgrna_files = sgrna_srrs.collect { srr -> file("${params.fastq_dir}/${srr}_{1,2,3}.fastq.gz") }.flatten()
            
            return [sid, mrna_files, sgrna_files]
        }

    if (params.limit > 0) {
        samples_ch = samples_ch.take(params.limit)
    }

    // Step 1: Build Indices
    std_idx = BUILD_INDEX_STANDARD(fa, gtf)
    kite_idx = BUILD_INDEX_KITE(features)

    // Step 2: Quantify cDNA and Guides in parallel per sample
    std_counts = KB_COUNT_STANDARD(samples_ch.map { sid, mrna, sgrna -> [sid, mrna] }, std_idx.index.collect(), std_idx.t2g.collect(), params.chemistry)
    kite_counts = KB_COUNT_KITE(samples_ch.map { sid, mrna, sgrna -> [sid, sgrna] }, kite_idx.index.collect(), kite_idx.t2g.collect(), params.chemistry)

    // Step 3: Merge modalities per sample
    merge_ch = std_counts.h5ad.join(kite_counts.h5ad)
    merged_samples = MERGE_MODALITIES(merge_ch)
    
    // Step 4: Final Global Concatenation
    uncompressed_final = CONCATENATE_SAMPLES(merged_samples.h5ad.collect())

    // Step 5: Final HDF5 Compression
    COMPRESS_FINAL_H5AD(uncompressed_final.h5ad)
}
