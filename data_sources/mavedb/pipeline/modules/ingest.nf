process DOWNLOAD_MAVEDB {
    tag 'Zenodo dump'
    label 'ingest_download'
    cache false

    input:
    val job

    output:
    tuple val(job), path('download/*', arity: '1'), path('zenodo_file_metadata.json'), emit: archive

    script:
    def p = MavedbCache.request(job)
    """
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/download_zenodo.py')} \
        --url ${MavedbCache.quote(p.zenodo_url)} --output download \
        --cache-dir ${MavedbCache.quote(p.download_dir ?: p.outdir + '/downloads')} \
        --metadata-output zenodo_file_metadata.json --no-unpack
    """
}

process INGEST_MAVEDB_ARCHIVE {
    tag 'Validate and unpack dump'
    label 'ingest_unpack'
    cache false

    input:
    tuple val(job), path(archive, stageAs: 'source.zip'), path(zenodo_metadata, stageAs: 'zenodo_file_metadata.json')

    output:
    path 'ingest/sha256-*', emit: bundle
    path 'result.json', emit: result

    script:
    def p = MavedbCache.request(job)
    """
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/ingest_mavedb_archive.py')} \
        --archive source.zip --zenodo-metadata zenodo_file_metadata.json \
        --source ${MavedbCache.quote(p.archive ?: p.zenodo_url)} \
        --expected-sha256 ${MavedbCache.quote(p.archive_sha256 ?: '')} --output ingest
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output ingest/sha256-* --require manifest.json --require mavedb_dump/main.json > result.json
    """
}

workflow MAVEDB_INGEST {
    take:
    engine

    main:
    decisions = Channel.of([engine: engine, stage: 'ingest', inputs: [:]]).map { MavedbCache.resolve(it) }
    missing = decisions.filter { it.spec != null }
    local = missing.filter { MavedbCache.request(it).archive != null }.map { job ->
        [job, file(MavedbCache.request(job).archive), file("${projectDir}/assets/no_zenodo_metadata.json")]
    }
    DOWNLOAD_MAVEDB(missing.filter { MavedbCache.request(it).archive == null })
    INGEST_MAVEDB_ARCHIVE(local.mix(DOWNLOAD_MAVEDB.out.archive))
    ready = decisions.filter { it.path != null }.mix(INGEST_MAVEDB_ARCHIVE.out.result.map { MavedbCache.result(it) })

    emit:
    result = ready
    bundle = ready.map { file(it.path) }
    dump = ready.map { file(it.path).resolve('mavedb_dump') }
    manifest = ready.map { file(it.path).resolve('manifest.json') }
}
