process DOWNLOAD_MAVEDB {
    tag 'Zenodo dump'
    label 'ingest_download'

    input:
    env 'zenodo_url'
    path download_script, stageAs: 'download_zenodo.py'
    env 'python_executable'
    env 'cache_dir'

    output:
    path 'download/*', arity: '1', emit: archive
    path 'zenodo_file_metadata.json', emit: zenodo_metadata

    script:
    '''
    "$python_executable" download_zenodo.py \
        --url "$zenodo_url" \
        --output download \
        --cache-dir "$cache_dir" \
        --metadata-output zenodo_file_metadata.json \
        --no-unpack
    '''
}

process INGEST_MAVEDB_ARCHIVE {
    tag 'Validate and unpack dump'
    label 'ingest_unpack'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple path(archive, stageAs: 'source.zip'),
          path(zenodo_metadata, stageAs: 'zenodo_file_metadata.json')
    env 'source_uri'
    env 'expected_sha256'
    path ingest_script, stageAs: 'ingest_mavedb_archive.py'
    env 'python_executable'

    output:
    path 'ingest/sha256-*', arity: '1', emit: bundle

    script:
    '''
    "$python_executable" ingest_mavedb_archive.py \
        --archive source.zip \
        --zenodo-metadata zenodo_file_metadata.json \
        --source "$source_uri" \
        --expected-sha256 "$expected_sha256" \
        --output ingest
    '''
}
