#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DOWNLOAD_MAVEDB; INGEST_MAVEDB_ARCHIVE } from './modules/ingest'

workflow MAVEDB_INGEST {
    main:
        if (params.archive_sha256 && !(params.archive_sha256 ==~ /[a-fA-F0-9]{64}/)) {
            error '--archive_sha256 must contain exactly 64 hexadecimal characters'
        }
        archive_ch = null
        source = null
        if (params.archive) {
            archive_ch = Channel.value([
                file(params.archive, checkIfExists: true),
                file("${projectDir}/assets/no_zenodo_metadata.json", checkIfExists: true)
            ])
            source = file(params.archive).toUri().toString()
        } else {
            if (!params.zenodo_url || !params.zenodo_url.startsWith('https://')) {
                error 'Provide an HTTPS --zenodo_url or a local ZIP with --archive'
            }
            download_dir = file(
                params.download_dir ?: "${params.outdir}/downloads"
            ).toAbsolutePath().toString()
            DOWNLOAD_MAVEDB(
                params.zenodo_url,
                file("${projectDir}/../scripts/download_zenodo.py"),
                params.python,
                download_dir
            )
            archive_ch = DOWNLOAD_MAVEDB.out.archive.combine(
                DOWNLOAD_MAVEDB.out.zenodo_metadata
            )
            source = params.zenodo_url
        }

        INGEST_MAVEDB_ARCHIVE(
            archive_ch,
            source,
            params.archive_sha256 ?: '',
            file("${projectDir}/../scripts/ingest_mavedb_archive.py"),
            params.python
        )
        bundle_ch = INGEST_MAVEDB_ARCHIVE.out.bundle
        dump_ch = bundle_ch.map { it.resolve('mavedb_dump') }
        manifest_ch = bundle_ch.map { it.resolve('manifest.json') }

        emit:
        bundle = bundle_ch
        dump = dump_ch
        manifest = manifest_ch
}

workflow {
    if (params.help) {
        log.info '''
MaveDB raw-data ingest

nextflow run data_sources/mavedb/pipeline -profile local [options]

  --archive PATH          Use an existing ZIP instead of downloading
  --zenodo_url URL        Override the pinned HTTPS Zenodo archive URL
  --archive_sha256 HEX    Require an expected SHA-256 before extraction
  --outdir PATH           Published results (default: output/nextflow)
  --download_dir PATH     Persistent ZIP cache (default: <outdir>/downloads)
  --python PATH           Python executable (default: repository .venv/bin/python)
  -resume                 Reuse successful tasks and their work directories

Only raw-data ingest is implemented. See pipeline/README.md for setup and outputs.
'''
    } else {
        MAVEDB_INGEST()
    }
}
