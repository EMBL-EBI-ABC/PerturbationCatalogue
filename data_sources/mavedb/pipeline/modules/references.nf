include { FETCH_MAVEDB_TAXONOMY; IMPORT_MAVEDB_TAXONOMY; REUSE_MAVEDB_TAXONOMY } from './taxonomy'
include { DOWNLOAD_MAVEDB_VARIANT_SOURCES } from './references_download'
include { PREPARE_MAVEDB_VARIANTS; IMPORT_MAVEDB_VARIANTS; REUSE_MAVEDB_VARIANTS; PACKAGE_MAVEDB_REFERENCES } from './references_prepare'

workflow MAVEDB_REFERENCES {
    take:
    engine

    main:
    decisions = Channel.of([engine: engine, stage: 'references', inputs: [:]])
        .map { MavedbCache.resolve(it) }
    reused = decisions.filter { it.path != null }
    missing = decisions.filter { it.spec != null }

    taxonomyImports = missing.filter { job ->
        MavedbCache.request(job).taxonomy != null
    }.map { job ->
        tuple(job, file(MavedbCache.request(job).taxonomy, checkIfExists: true))
    }
    taxonomyReuse = missing.filter { job ->
        def p = MavedbCache.request(job)
        def spec = MavedbCache.specification(job)
        p.taxonomy == null && spec.taxonomy_reuse != null
    }.map { job ->
        def reuse = MavedbCache.specification(job).taxonomy_reuse
        def source = reuse instanceof Map ? (reuse.path ?: reuse.directory) : reuse
        tuple(job, file(source.toString(), checkIfExists: true))
    }
    taxonomyFetch = missing.filter { job ->
        def p = MavedbCache.request(job)
        def spec = MavedbCache.specification(job)
        p.taxonomy == null && spec.taxonomy_reuse == null
    }
    IMPORT_MAVEDB_TAXONOMY(taxonomyImports)
    REUSE_MAVEDB_TAXONOMY(taxonomyReuse)
    FETCH_MAVEDB_TAXONOMY(taxonomyFetch)
    taxonomyComponents = IMPORT_MAVEDB_TAXONOMY.out.bundle
        .mix(REUSE_MAVEDB_TAXONOMY.out.bundle, FETCH_MAVEDB_TAXONOMY.out.bundle)

    variantImports = missing.filter { job ->
        def p = MavedbCache.request(job)
        p.validation_clinvar != null && p.validation_gwas != null
    }.map { job ->
        def p = MavedbCache.request(job)
        tuple(job, '', 'imported', '',
              file(p.validation_clinvar, checkIfExists: true),
              file(p.validation_gwas, checkIfExists: true))
    }
    variantReuse = missing.filter { job ->
        def p = MavedbCache.request(job)
        def spec = MavedbCache.specification(job)
        !(p.validation_clinvar && p.validation_gwas) && spec.variants_reuse != null
    }.map { job ->
        def spec = MavedbCache.specification(job)
        def reuse = spec.variants_reuse
        def root = reuse instanceof Map ? (reuse.path ?: reuse.directory) : reuse
        def source = root.toString().endsWith('/variants') ? root : "${root}/variants"
        tuple(job, '', 'imported', '', file(source.toString(), checkIfExists: true))
    }
    externalSourceManifests = missing.filter { job ->
        def p = MavedbCache.request(job)
        def spec = MavedbCache.specification(job)
        !(p.validation_clinvar && p.validation_gwas) && spec.variants_reuse == null &&
            p.reference_sources_manifest != null
    }.map { job ->
        def p = MavedbCache.request(job)
        def manifest = file(p.reference_sources_manifest.toString(), checkIfExists: true)
        tuple(job, file(manifest.parent.toString(), checkIfExists: true),
              manifest.fileName.toString(), 'external', manifest.toString())
    }
    retainedSourceManifests = missing.filter { job ->
        def p = MavedbCache.request(job)
        def spec = MavedbCache.specification(job)
        !(p.validation_clinvar && p.validation_gwas) && spec.variants_reuse == null &&
            p.reference_sources_manifest == null && !p.refresh_references &&
            spec.retained_source_manifest != null
    }.map { job ->
        def source = file(MavedbCache.specification(job).retained_source_manifest.toString(),
                          checkIfExists: true)
        tuple(job, file(source.parent.toString(), checkIfExists: true),
              source.fileName.toString(), 'retained', '')
    }
    freshSourceJobs = missing.filter { job ->
        def p = MavedbCache.request(job)
        def spec = MavedbCache.specification(job)
        !(p.validation_clinvar && p.validation_gwas) && spec.variants_reuse == null &&
            p.reference_sources_manifest == null &&
            (p.refresh_references || spec.retained_source_manifest == null)
    }
    DOWNLOAD_MAVEDB_VARIANT_SOURCES(freshSourceJobs)
    generatedSourceManifests = DOWNLOAD_MAVEDB_VARIANT_SOURCES.out.bundle.map { job, sources ->
        tuple(job, sources, 'manifest.json', 'generated', '')
    }
    sourceSelections = externalSourceManifests
        .mix(retainedSourceManifests, generatedSourceManifests)
    PREPARE_MAVEDB_VARIANTS(sourceSelections)

    IMPORT_MAVEDB_VARIANTS(variantImports)
    REUSE_MAVEDB_VARIANTS(variantReuse)
    importedVariants = IMPORT_MAVEDB_VARIANTS.out.bundle.map {
        job, variants, sourceRoot, manifestName, sourceMode, sourceLocator ->
        [
            components: tuple(job, variants, manifestName, sourceMode, sourceLocator),
            sources: tuple(job, sourceRoot, manifestName, sourceMode, sourceLocator),
        ]
    }
    reusedVariants = REUSE_MAVEDB_VARIANTS.out.bundle.map {
        job, variants, sourceRoot, manifestName, sourceMode, sourceLocator ->
        [
            components: tuple(job, variants, manifestName, sourceMode, sourceLocator),
            sources: tuple(job, sourceRoot, manifestName, sourceMode, sourceLocator),
        ]
    }
    preparedVariants = PREPARE_MAVEDB_VARIANTS.out.bundle.map {
        job, variants, manifestName, sourceMode, sourceLocator ->
        [
            components: tuple(job, variants, manifestName, sourceMode, sourceLocator),
            sources: null,
        ]
    }
    variantBundles = importedVariants.mix(reusedVariants, preparedVariants)
    variantComponents = variantBundles.map { it.components }
    importedSourceRoots = importedVariants.mix(reusedVariants).map { it.sources }
    sourceComponents = importedSourceRoots.mix(sourceSelections)

    packages = taxonomyComponents.map { [taxonomy: it] }
        .combine(variantComponents.map { [variant: it] })
        .combine(sourceComponents.map { [source: it] })
        .map { taxonomyRecord, variantRecord, sourceRecord ->
            def taxonomy = taxonomyRecord.taxonomy
            def variant = variantRecord.variant
            def source = sourceRecord.source
            def job = taxonomy[0]
            if (job.spec != variant[0].spec || job.spec != source[0].spec) {
                error 'Reference component tasks received different stage specifications'
            }
            if (variant[3] != source[3] || variant[2] != source[2]) {
                error 'Reference source and variant component modes differ'
            }
            tuple(job, taxonomy[1], variant[1], source[1], source[2], source[3], source[4])
        }
    PACKAGE_MAVEDB_REFERENCES(packages)
    computed = PACKAGE_MAVEDB_REFERENCES.out.result.map { MavedbCache.result(it) }
    ready = reused.mix(computed).map { reference ->
        def bundle = file(reference.path, checkIfExists: true)
        def taxonomyManifest = new groovy.json.JsonSlurper()
            .parseText(bundle.resolve('taxonomy/manifest.json').text)
        reference + [
            stage: 'references',
            taxonomy_source: taxonomyManifest.source,
            clinvar: bundle.resolve('variants/mave_clinvar_variants.tsv'),
            gwas: bundle.resolve('variants/mave_gwas_variants.tsv'),
        ]
    }

    emit:
    result = ready
    bundle = ready.map { file(it.path, checkIfExists: true) }
    taxonomy = ready.map { file(it.path, checkIfExists: true).resolve('taxonomy/taxcat.dmp') }
    clinvar = ready.map { file(it.path, checkIfExists: true).resolve('variants/mave_clinvar_variants.tsv') }
    gwas = ready.map { file(it.path, checkIfExists: true).resolve('variants/mave_gwas_variants.tsv') }
}
