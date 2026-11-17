process FETCH_MAVEDB_TAXONOMY {
    tag 'NCBI taxonomy categories'
    label 'reference_download'
    cache false

    input:
    val job

    output:
    path 'taxonomy', emit: bundle
    path 'result.json', emit: result

    script:
    def p = MavedbCache.request(job)
    """
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/download_taxonomy.py')} \
        --url ${MavedbCache.quote(p.taxonomy_url)} --output taxonomy
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output taxonomy --require manifest.json --require taxcat.dmp > result.json
    """
}

process IMPORT_MAVEDB_TAXONOMY {
    tag 'Supplied taxonomy snapshot'
    label 'reference_download'
    cache false

    input:
    tuple val(job), path(taxonomy, stageAs: 'taxcat.dmp')

    output:
    path 'taxonomy', emit: bundle
    path 'result.json', emit: result

    script:
    def p = MavedbCache.request(job)
    def manifest = groovy.json.JsonOutput.toJson(
        [schema_version: 1, source: p.taxonomy, path: 'taxcat.dmp'] + MavedbCache.specification(job).identities.source
    )
    """
    mkdir taxonomy
    cp taxcat.dmp taxonomy/taxcat.dmp
    printf '%s\\n' ${MavedbCache.quote(manifest)} > taxonomy/manifest.json
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output taxonomy --require manifest.json --require taxcat.dmp > result.json
    """
}

workflow MAVEDB_TAXONOMY {
    take:
    engine

    main:
    decisions = Channel.of([engine: engine, stage: 'taxonomy', inputs: [:]]).map { MavedbCache.resolve(it) }
    missing = decisions.filter { it.spec != null }
    FETCH_MAVEDB_TAXONOMY(missing.filter { MavedbCache.request(it).taxonomy == null })
    IMPORT_MAVEDB_TAXONOMY(missing.filter { MavedbCache.request(it).taxonomy != null }.map { job ->
        [job, file(MavedbCache.request(job).taxonomy)]
    })
    ready = decisions.filter { it.path != null }.mix(
        FETCH_MAVEDB_TAXONOMY.out.result.map { MavedbCache.result(it) },
        IMPORT_MAVEDB_TAXONOMY.out.result.map { MavedbCache.result(it) }
    )

    emit:
    result = ready
    taxonomy = ready.map { file(it.path).resolve('taxcat.dmp') }
    source = ready.map { MavedbCache.request([engine: engine]).with { it.taxonomy ?: it.taxonomy_url } }
}
