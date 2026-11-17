process GENERATE_COMPATIBILITY_TABLES {
    tag 'Raw-score compatibility tables'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(raw_bundle, stageAs: 'raw_bundle'),
          path(taxonomy, stageAs: 'taxcat.dmp'), path(processed, stageAs: 'processed')

    output:
    path 'tables', emit: bundle
    path 'result.json', emit: result

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote("${job.engine.scripts}/generate_compatibility_tables.py")} \
        --input raw_bundle/mavedb_dump --processed processed \
        --taxonomy taxcat.dmp --output tables
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output tables --require completion.json --require maves/singles.tsv --require maves/doubles.tsv --require maves/wide.tsv > result.json
    """
}

workflow MAVEDB_TABLES {
    take:
    engine
    raw
    taxonomy
    processed

    main:
    decisions = raw.combine(taxonomy).combine(processed).map { source, reference, canonical ->
        MavedbCache.resolve([
            engine: engine, stage: 'tables',
            inputs: [ingest: source.path, taxonomy: reference.path, process: canonical.path]
        ])
    }
    reused = decisions.filter { it.path != null }
    jobs = decisions.filter { it.spec != null }.map { job ->
        [job, file(job.inputs.ingest, checkIfExists: true),
         file("${job.inputs.taxonomy}/taxcat.dmp", checkIfExists: true),
         file(job.inputs.process, checkIfExists: true)]
    }
    GENERATE_COMPATIBILITY_TABLES(jobs)
    computed = GENERATE_COMPATIBILITY_TABLES.out.result.map { MavedbCache.result(it) }
    ready = reused.mix(computed)

    emit:
    result = ready
    bundle = ready.map { file(it.path, checkIfExists: true) }
}
