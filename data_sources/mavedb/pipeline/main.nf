#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DOWNLOAD_MAVEDB; INGEST_MAVEDB_ARCHIVE } from './modules/ingest'
include { MAVEDB_TAXONOMY } from './modules/taxonomy'
include { PROCESS_MAVEDB } from './modules/process'
include { MAVEDB_GMM } from './modules/gmm'

workflow MAVEDB_PROCESS {
    take:
    raw_bundle
    taxonomy
    taxonomy_source

    main:
    processor = file("${projectDir}/../scripts/process_mavedb.py")
    PROCESS_MAVEDB(
        raw_bundle, taxonomy, taxonomy_source, processor, params.python
    )
    processed_ch = PROCESS_MAVEDB.out.bundle
    scores_ch = processed_ch.map { it.resolve('scores.tsv') }
    metadata_ch = processed_ch.map { it.resolve('metadata.tsv') }
    audit_ch = processed_ch.map { it.resolve('dataset_processing.tsv') }
    manifest_ch = processed_ch.map { it.resolve('processing_manifest.json') }

    emit:
    bundle = processed_ch
    scores = scores_ch
    metadata = metadata_ch
    audit = audit_ch
    manifest = manifest_ch
}

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
MaveDB ingest and canonical processing

nextflow run data_sources/mavedb/pipeline -profile local [options]

  --archive PATH          Use an existing ZIP instead of downloading
  --zenodo_url URL        Override the pinned HTTPS Zenodo archive URL
  --archive_sha256 HEX    Require an expected SHA-256 before extraction
  --outdir PATH           Published results (default: output/nextflow)
  --download_dir PATH     Persistent ZIP cache (default: <outdir>/downloads)
  --python PATH           Python executable (default: repository .venv/bin/python)
  --taxonomy PATH         Use an existing NCBI categories.dmp (otherwise fetched)
  --taxonomy_url URL      Override the NCBI taxonomy archive source
  --ingest_only           Stop after raw-data ingest; no taxonomy acquisition
  --run_gmm               Fit GMMs after processing (opt-in)
  --processed_dir PATH    Start GMM fitting from existing canonical tables
  --gmm_profile NAME      production (default) or bounded smoke settings
  --gmm_model NAME        modes (default) or flexible
  --gmm_study_ids IDS     Comma-separated score-set IDs, or a params-file list
  -resume                 Reuse successful tasks and their work directories

See pipeline/README.md for setup, offline smoke testing, and outputs.
'''
    } else {
        if (params.ingest_only && params.run_gmm) {
            error '--ingest_only cannot be combined with --run_gmm'
        }
        if (params.processed_dir) {
            if (!params.run_gmm) error '--processed_dir requires --run_gmm'
            MAVEDB_GMM(Channel.value(file(params.processed_dir, checkIfExists: true)))
        } else {
            MAVEDB_INGEST()
            if (!params.ingest_only) {
                MAVEDB_TAXONOMY()
                MAVEDB_PROCESS(
                    MAVEDB_INGEST.out.bundle,
                    MAVEDB_TAXONOMY.out.taxonomy,
                    MAVEDB_TAXONOMY.out.source
                )
                if (params.run_gmm) MAVEDB_GMM(MAVEDB_PROCESS.out.bundle)
            }
        }
    }
}
