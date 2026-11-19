process FETCH_MAVEDB_TAXONOMY {
    tag 'NCBI taxonomy categories'
    label 'reference_download'
    cache false

    input:
    val job

    output:
    tuple val(job), path('taxonomy'), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind taxonomy --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory taxonomy; then
        :
    else
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/download_taxonomy.py')} \
            --url ${MavedbCache.quote(p.taxonomy_url)} --output taxonomy
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind taxonomy --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory taxonomy
    fi
    """
}

process IMPORT_MAVEDB_TAXONOMY {
    tag 'Supplied taxonomy snapshot'
    label 'reference_preparation'
    cache false

    input:
    tuple val(job), path(taxonomy_input, stageAs: 'taxcat.dmp')

    output:
    tuple val(job), path('taxonomy'), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind taxonomy --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory taxonomy; then
        :
    else
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} import-taxonomy \
            --source taxcat.dmp --source-label ${MavedbCache.quote(p.taxonomy)} --out taxonomy
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind taxonomy --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory taxonomy
    fi
    """
}

process REUSE_MAVEDB_TAXONOMY {
    tag 'Saved taxonomy snapshot'
    label 'reference_preparation'
    cache false

    input:
    tuple val(job), path(saved_taxonomy, stageAs: 'saved_taxonomy')

    output:
    tuple val(job), path('taxonomy'), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind taxonomy --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory taxonomy; then
        :
    else
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} reuse-taxonomy \
            --source saved_taxonomy --out taxonomy
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind taxonomy --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory taxonomy
    fi
    """
}
