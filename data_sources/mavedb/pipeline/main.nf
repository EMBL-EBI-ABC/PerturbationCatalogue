#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAVEDB_INGEST } from './modules/ingest'
include { MAVEDB_REFERENCES } from './modules/references'
include { MAVEDB_PROCESS } from './modules/process'
include { MAVEDB_TABLES } from './modules/tables'
include { MAVEDB_GMM } from './modules/gmm'
include { MAVEDB_MERGE } from './modules/merge'
include { MAVEDB_QC } from './modules/qc'
include { MAVEDB_REPORTS } from './modules/reports'
include { MAVEDB_VALIDATION } from './modules/validation'
include { MAVEDB_RESIDUE_CHECK } from './modules/residue_check'

workflow {
    if (params.help) {
        log.info '''
MaveDB pipeline: generate normal output tables; reuse completed results by default.

nextflow run data_sources/mavedb/pipeline -profile local [options]

  --outdir PATH            Saved analysis and durable result cache
  --archive PATH           Local ZIP instead of the pinned Zenodo URL
  --taxonomy PATH          Import a local categories.dmp into references
  --reference_sources_manifest PATH
                           Verified raw Open Targets source bundle manifest
  --references_cache_dir PATH
                           Resumable reference task components (default under tmp)
  --refresh_references     Refresh sources; requires --rerun references
  --processed_dir PATH     Existing canonical scores.tsv/metadata.tsv
  --ingest_only            Stop after ingest
  --run_gmm                Enable GMM fitting and merging
  --run_qc                 Enable basic QC
  --run_reports            Enable per-study reports and GMM dependencies
  --run_validation         Validate GMM fits against prepared ClinVar/GWAS tables
  --uniprot_residues PATH   Check reference residues against a local ProtVar TSV
  --validation_clinvar PATH Imported ClinVar TSV (requires paired GWAS path)
  --validation_gwas PATH   Imported GWAS TSV (requires paired ClinVar path)
  --rerun STAGES           Comma-separated stage names, or all enabled stages
                          ingest,references,process,tables,gmm,merge,qc,reports,validation,residue_check
  --gmm_profile NAME       production (default) or bounded smoke
  --gmm_model NAME         modes (default) or flexible
  --validation_formats FORMATS
                           Comma-separated plot formats: png,pdf (default)
  --validation_dpi N       Raster plot resolution (default 150)
  --validation_no_plots    Write validation tables without figures
  --gmm_study_ids IDS      Comma-separated study IDs or params-file list
  -resume                  Optional native task resume; --rerun still forces work

See pipeline/README.md for all fitting/resource settings and migration.
'''
    } else {
        def root = file(params.outdir).toAbsolutePath()
        java.nio.file.Files.createDirectories(root)
        def engine = [python: file(params.python).toAbsolutePath().toString(),
                      runtime: "${projectDir}/../scripts/pipeline_runtime.py",
                      scripts: "${projectDir}/../scripts",
                      workDir: workflow.workDir.toString()]
        def writerLock = MavedbCache.acquire(engine, root)
        workflow.onComplete { writerLock.outputStream.close(); writerLock.waitFor() }
        def prepared = MavedbCache.call(engine, ['prepare'], params as Map)
        engine.request = prepared.request
        def request = new groovy.json.JsonSlurper().parseText(file(prepared.request).text)
        def enabled = request.enabled as Set
        def results = Channel.empty()
        def processed
        def raw
        def references
        if (params.processed_dir) {
            def resolved = MavedbCache.resolve([engine: engine, stage: 'process', inputs: [:]])
            processed = Channel.value(resolved)
        } else if (request.saved_processing) {
            processed = Channel.value([path: request.saved_processing, stage: 'process'])
        } else {
            MAVEDB_INGEST(engine)
            raw = MAVEDB_INGEST.out.result
            results = results.mix(raw.map { [ingest: it.path] })
        }
        if (!params.ingest_only && enabled.contains('references')) {
            MAVEDB_REFERENCES(engine)
            references = MAVEDB_REFERENCES.out.result
            results = results.mix(references.map { [references: it.path] })
        }
        if (!params.ingest_only && !params.processed_dir && !request.saved_processing) {
            MAVEDB_PROCESS(engine, raw, references)
            processed = MAVEDB_PROCESS.out.result
        }
        if (!params.ingest_only) {
            results = results.mix(processed.map { [process: it.path] })
            if (enabled.contains('residue_check')) {
                MAVEDB_RESIDUE_CHECK(engine, processed, request.uniprot_residues)
                results = results.mix(MAVEDB_RESIDUE_CHECK.out.result.map { [residue_check: it.path] })
            }
            def tables
            if (enabled.contains('tables')) {
                MAVEDB_TABLES(engine, raw, references, processed)
                tables = MAVEDB_TABLES.out.result
                results = results.mix(tables.map { [tables: it.path] })
            }
            if (enabled.contains('gmm')) {
                def gmmInput = enabled.contains('validation')
                    ? processed.combine(references).map { canonical, reference -> canonical }
                    : processed
                MAVEDB_GMM(engine, gmmInput)
                def fitted = MAVEDB_GMM.out.result
                results = results.mix(fitted.map { [gmm: it.path] })
                MAVEDB_MERGE(engine, fitted)
                def merged = MAVEDB_MERGE.out.result
                results = results.mix(merged.map { [merge: it.path] })
                if (enabled.contains('reports')) {
                    MAVEDB_REPORTS(engine, processed, merged)
                    results = results.mix(MAVEDB_REPORTS.out.result.map { [reports: it.path] })
                }
                if (enabled.contains('validation')) {
                    MAVEDB_VALIDATION(engine, processed, merged, references)
                    results = results.mix(MAVEDB_VALIDATION.out.result.map { [validation: it.path] })
                }
            }
            if (enabled.contains('qc')) {
                MAVEDB_QC(engine, processed, tables)
                results = results.mix(MAVEDB_QC.out.result.map { [qc: it.path] })
            }
        }
        results.collect().map { values ->
            def paths = [:]
            values.each { paths.putAll(it) }
            MavedbCache.call(engine, ['finish', '--request', engine.request], paths)
        }.view { 'MaveDB results verified and published.' }
    }
}
