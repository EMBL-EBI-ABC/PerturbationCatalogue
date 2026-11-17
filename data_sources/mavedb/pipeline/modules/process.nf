process PROCESS_MAVEDB {
    tag 'Canonical metadata and scores'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(raw_bundle, stageAs: 'raw_bundle'), path(taxonomy, stageAs: 'taxcat.dmp')

    output:
    path 'processed', emit: bundle
    path 'result.json', emit: result

    script:
    def p = MavedbCache.request(job)
    """
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/process_mavedb.py')} \
        --input raw_bundle/mavedb_dump --ingest-manifest raw_bundle/manifest.json \
        --taxonomy taxcat.dmp --taxonomy-source ${MavedbCache.quote(p.taxonomy ?: p.taxonomy_url)} \
        --no-tables --output processed
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output processed --require processing_manifest.json --require scores.tsv --require metadata.tsv --require dataset_processing.tsv > result.json
    """
}

workflow MAVEDB_PROCESS {
    take:
    engine
    raw
    taxonomy

    main:
    decisions = raw.combine(taxonomy).map { source, reference ->
        MavedbCache.resolve([engine: engine, stage: 'process', inputs: [ingest: source.path, taxonomy: reference.path]])
    }
    PROCESS_MAVEDB(decisions.filter { it.spec != null }.map { job ->
        [job, file(job.inputs.ingest), file(job.inputs.taxonomy).resolve('taxcat.dmp')]
    })
    ready = decisions.filter { it.path != null }.mix(PROCESS_MAVEDB.out.result.map { MavedbCache.result(it) })

    emit:
    result = ready
    bundle = ready.map { file(it.path) }
    scores = ready.map { file(it.path).resolve('scores.tsv') }
    metadata = ready.map { file(it.path).resolve('metadata.tsv') }
    audit = ready.map { file(it.path).resolve('dataset_processing.tsv') }
    manifest = ready.map { file(it.path).resolve('processing_manifest.json') }
}
