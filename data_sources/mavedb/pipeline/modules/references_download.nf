process DOWNLOAD_MAVEDB_VARIANT_SOURCES {
    tag 'ClinVar and GWAS source bundle'
    label 'reference_download'
    cache false

    input:
    val job

    output:
    tuple val(job), path('sources'), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def downloader = job.engine.scripts + '/download_reference_data.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    def sourceUrlOptions = [
        [key: 'coding_variants_url', flag: 'coding-variants-url'],
        [key: 'evidence_url', flag: 'evidence-url'],
        [key: 'credible_sets_url', flag: 'credible-sets-url'],
    ].findAll { p[it.key] }.collect {
        "--${it.flag} ${MavedbCache.quote(p[it.key])}"
    }.join(' ')
    def refresh = p.refresh_references ? '--refresh' : ''
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind sources --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory sources; then
        :
    else
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(downloader)} \
            --only clinvar,gwas --out sources \
            --cache-dir ${MavedbCache.quote(cacheRoot + '/download_cache')} \
            ${sourceUrlOptions} ${refresh}
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind sources --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory sources
    fi
    """
}
