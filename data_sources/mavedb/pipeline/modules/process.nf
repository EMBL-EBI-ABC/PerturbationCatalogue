process PROCESS_MAVEDB {
    tag 'Canonical metadata and scores'
    label 'mavedb_processing'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path raw_bundle, stageAs: 'raw_bundle'
    path taxonomy, stageAs: 'taxcat.dmp'
    env 'taxonomy_source'
    path processor, stageAs: 'process_mavedb.py'
    env 'python_executable'

    output:
    path 'processed', emit: bundle

    script:
    '''
    "$python_executable" process_mavedb.py \
        --input raw_bundle/mavedb_dump \
        --ingest-manifest raw_bundle/manifest.json \
        --taxonomy taxcat.dmp \
        --taxonomy-source "$taxonomy_source" \
        --output processed
    '''
}
