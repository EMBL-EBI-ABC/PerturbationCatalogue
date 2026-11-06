process FETCH_MAVEDB_TAXONOMY {
    tag 'NCBI taxonomy categories'
    label 'reference_download'
    publishDir "${params.outdir}/references", mode: 'copy', overwrite: true

    input:
    path downloader, stageAs: 'download_taxonomy.py'
    env 'taxonomy_url'
    env 'python_executable'

    output:
    path 'taxonomy', emit: bundle

    script:
    '''
    "$python_executable" download_taxonomy.py \
        --url "$taxonomy_url" \
        --output taxonomy
    '''
}

workflow MAVEDB_TAXONOMY {
    main:
    taxonomy_ch = null
    taxonomy_source = null
    if (params.taxonomy) {
        taxonomy_ch = Channel.value(file(params.taxonomy, checkIfExists: true))
        taxonomy_source = file(params.taxonomy).toUri().toString()
    } else {
        FETCH_MAVEDB_TAXONOMY(
            file("${projectDir}/../scripts/download_taxonomy.py"),
            params.taxonomy_url,
            params.python
        )
        taxonomy_ch = FETCH_MAVEDB_TAXONOMY.out.bundle.map { it.resolve('taxcat.dmp') }
        taxonomy_source = params.taxonomy_url
    }

    emit:
    taxonomy = taxonomy_ch
    source = taxonomy_source
}
