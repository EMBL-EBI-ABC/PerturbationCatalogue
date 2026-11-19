process PREPARE_MAVEDB_VARIANTS {
    tag 'Prepare ClinVar and GWAS variants'
    label 'reference_preparation'
    cache false

    input:
    tuple val(job), path(source_root, stageAs: 'source_root'),
          val(manifest_name), val(source_mode), val(source_locator)

    output:
    tuple val(job), path('variants'), val(manifest_name), val(source_mode),
          val(source_locator), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def processor = job.engine.scripts + '/process_variant_data.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind variants --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory variants; then
        :
    else
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(processor)} \
            --sources ${MavedbCache.quote('source_root/' + manifest_name)} --only clinvar,gwas \
            --out variants --manifest manifest.json
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind variants --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory variants
    fi
    """
}

process IMPORT_MAVEDB_VARIANTS {
    tag 'Paired ClinVar and GWAS imports'
    label 'reference_preparation'
    cache false

    input:
    tuple val(job), val(manifest_name), val(source_mode), val(source_locator),
          path(clinvar, stageAs: 'clinvar.tsv'), path(gwas, stageAs: 'gwas.tsv')

    output:
    tuple val(job), path('variants'), path('source_stub'),
          val(manifest_name), val(source_mode), val(source_locator), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind variants --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory variants; then
        :
    else
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} import-variants \
            --clinvar clinvar.tsv --gwas gwas.tsv --out variants
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind variants --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory variants
    fi
    mkdir -p source_stub
    """
}

process REUSE_MAVEDB_VARIANTS {
    tag 'Saved imported validation variants'
    label 'reference_preparation'
    cache false

    input:
    tuple val(job), val(manifest_name), val(source_mode), val(source_locator),
          path(saved_variants, stageAs: 'saved_variants')

    output:
    tuple val(job), path('variants'), path('source_stub'),
          val(manifest_name), val(source_mode), val(source_locator), emit: bundle

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def cacheRoot = "${p.references_cache_dir}/${p.references_attempt}"
    """
    if ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
        --action restore --kind variants --root ${MavedbCache.quote(cacheRoot)} \
        --spec ${MavedbCache.quote(job.spec)} --directory variants; then
        :
    else
        cp -a saved_variants variants
        ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} checkpoint \
            --action save --kind variants --root ${MavedbCache.quote(cacheRoot)} \
            --spec ${MavedbCache.quote(job.spec)} --directory variants
    fi
    mkdir -p source_stub
    """
}

process PACKAGE_MAVEDB_REFERENCES {
    tag 'Verify and publish complete reference bundle'
    label 'reference_package'
    cache false

    input:
    tuple val(job), path(taxonomy, stageAs: 'taxonomy'), path(variants, stageAs: 'variants'),
          path(source_root, stageAs: 'source_root'), val(manifest_name),
          val(source_mode), val(source_locator)

    output:
    path 'references', emit: bundle
    path 'result.json', emit: result

    script:
    def p = MavedbCache.request(job)
    def helper = job.engine.scripts + '/package_reference_outputs.py'
    def sourceManifestOption = manifest_name
        ? "--source-manifest ${MavedbCache.quote('source_root/' + manifest_name)}"
        : ''
    def sourceLocatorOption = source_mode == 'external' && source_locator
        ? "--source-manifest-locator ${MavedbCache.quote(source_locator)}"
        : ''
    def replaceableRoots = [
        p.outdir + '/references',
        p.references_cache_dir,
        job.engine.workDir,
    ].unique().collect { "--replaceable-root ${MavedbCache.quote(it)}" }.join(' ')
    """
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(helper)} assemble \
        --taxonomy taxonomy --variants variants --source-mode ${MavedbCache.quote(source_mode)} \
        ${sourceManifestOption} ${sourceLocatorOption} ${replaceableRoots} --out references
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output references \
        --require manifest.json --require taxonomy/manifest.json --require taxonomy/taxcat.dmp \
        --require variants/manifest.json --require variants/mave_clinvar_variants.tsv \
        --require variants/mave_gwas_variants.tsv > result.json
    """
}
