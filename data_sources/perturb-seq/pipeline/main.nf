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
    tuple val(sample_id), path(std_h5ad, stageAs: 'std_adata.h5ad'), path(kite_h5ad, stageAs: 'kite_adata.h5ad')
    
    output:
    path "${sample_id}_merged.h5ad"
    
    script:
    """
    #!/usr/bin/env python3
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scipy.sparse as sp
    import sys
    
    srr_id = "${sample_id}"
    
    # Read in backed mode to avoid loading the full gene matrix into memory
    adata_std = ad.read_h5ad("${std_h5ad}", backed='r')
    adata_kite = ad.read_h5ad("${kite_h5ad}", backed='r')
    
    # Align the kite (guide) matrix rows to the standard (cDNA) cell barcodes
    # We use sparse indexing to avoid densifying the guide matrix
    kite_obs_map = pd.Series(np.arange(adata_kite.n_obs), index=adata_kite.obs_names)
    target_indices = kite_obs_map.reindex(adata_std.obs_names).values
    
    mask = ~pd.isna(target_indices)
    valid_indices = target_indices[mask].astype(int)
    
    # Construct the aligned sparse matrix for guides
    # We use the same dtype as the source matrix
    X_found = adata_kite.X[valid_indices, :]
    row_indices = np.where(mask)[0]
    
    # Build a new CSR matrix of the correct shape (n_cells_std x n_guides)
    # If the kite matrix is sparse, we build it from COO components for efficiency
    if sp.issparse(X_found):
        coo = X_found.tocoo()
        new_row_indices = row_indices[coo.row]
        guides_sparse = sp.csr_matrix(
            (coo.data, (new_row_indices, coo.col)),
            shape=(adata_std.n_obs, adata_kite.n_vars),
            dtype=X_found.dtype
        )
    else:
        # Fallback if kite.X is dense (unlikely but possible)
        guides_sparse = np.zeros((adata_std.n_obs, adata_kite.n_vars), dtype=X_found.dtype)
        guides_sparse[mask, :] = X_found
    
    # Create a new AnnData object for the merged SRR
    # We copy the obs/var/X from the standard matrix
    # Note: adata_std is backed, so we read it into memory for this step
    # or better, we create a new one and let it stream during write.
    
    # For a single SRR, we can afford to load .obs and .var
    adata_merged = ad.AnnData(
        X=adata_std.X,
        obs=adata_std.obs.copy(),
        var=adata_std.var.copy(),
        uns=adata_std.uns.copy()
    )
    
    # Store guide counts and metadata
    adata_merged.obsm['guides'] = guides_sparse
    adata_merged.uns['guide_names'] = adata_kite.var_names.tolist()
    
    # Add SRR label and ensure barcodes are globally unique across all SRRs
    adata_merged.obs['SRR_run'] = srr_id
    adata_merged.obs_names = [f"{barcode}-{srr_id}" for barcode in adata_merged.obs_names]
    
    adata_merged.write_h5ad(f"{srr_id}_merged.h5ad")
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
    import h5py
    import scipy.sparse as sp
    
    # Try to import experimental features
    try:
        from anndata.experimental import concat_on_disk, write_elem
    except ImportError:
        sys.exit("Error: Your anndata version is too old. Please use anndata >= 0.10.0 for concat_on_disk support.")

    files = sorted(glob.glob("*_merged.h5ad"))
    if not files:
        sys.exit("No merged h5ad files found.")

    # Use concat_on_disk to combine all SRRs into one massive dataset without loading matrices into memory.
    # This handles .X, .obs, and .var. 
    # join="inner" ensures we only keep genes present in all files, which is memory-efficient.
    print(f"Concatenating {len(files)} files on disk...")
    concat_on_disk(files, "experiment_final.h5ad", join="inner")

    # Manually append .obsm['guides'] and .uns['guide_names'] to the final h5ad file.
    # We load guide counts as sparse matrices. Since the guide matrix (n_cells x n_guides) 
    # is significantly smaller than the gene matrix, vstacking it in memory is generally safe.
    print("Concatenating guide counts (obsm['guides'])...")
    guide_mats = []
    for f in files:
        # Load in backed mode to only read obsm
        a = ad.read_h5ad(f, backed='r')
        # Ensure we have a sparse matrix
        if 'guides' in a.obsm:
            guide_mats.append(sp.csr_matrix(a.obsm['guides']))
        else:
            sys.exit(f"Error: .obsm['guides'] missing in {f}")

    if guide_mats:
        full_guides = sp.vstack(guide_mats)
        
        # Take guide names from the first file (assuming they are identical across all SRRs)
        a0 = ad.read_h5ad(files[0], backed='r')
        guide_names = a0.uns['guide_names']
        
        # Write to the existing HDF5 file
        with h5py.File("experiment_final.h5ad", "a") as f_out:
            write_elem(f_out, "obsm/guides", full_guides)
            write_elem(f_out, "uns/guide_names", guide_names)
            
    print("Final experiment_final.h5ad created successfully.")
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
    if (params.limit > 0) {
        read_pairs_ch = read_pairs_ch.take(params.limit)
    }

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
