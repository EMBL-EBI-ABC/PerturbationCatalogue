process MAVEDB_GMM_VALIDATION {
    tag 'GMM validation against controls, ClinVar and GWAS'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(merged, stageAs: 'gmm'),
          path(scores, stageAs: 'scores.tsv'),
          path(metadata, stageAs: 'metadata.tsv'),
          path(clinvar, stageAs: 'clinvar.tsv'),
          path(gwas, stageAs: 'gwas.tsv')

    output:
    path 'validation', emit: bundle
    path 'result.json', emit: result

    script:
    def request = MavedbCache.request(job)
    def modelRoot = (request.gmm_profile == 'smoke' ? 'smoke/' : '') + "maves/${request.gmm_model}"
    def formats = request.validation_formats instanceof List
        ? request.validation_formats
        : request.validation_formats.toString().tokenize(',')
    def formatArguments = MavedbCache.arguments(formats.unique())
    def dpi = request.validation_dpi ?: 150
    def noPlots = request.validation_no_plots ? '--no-plots' : ''
    def requiredOutputs = [
        'controls.tsv',
        'clinvar_matches.tsv',
        'gwas_matches.tsv',
        'ctrl_roc.tsv',
        'clinvar_roc.tsv',
        'gwas_roc.tsv',
        'auc_summary.tsv',
        'counts.tsv',
        'manifest.json',
    ].collect { "--require ${MavedbCache.quote(it)}" }.join(' ')
    """
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export MPLBACKEND=Agg
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    mkdir -p .mplconfig
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote("${job.engine.scripts}/../analysis/mave_gmm_validation.py")} \
        --gmm-dir ${MavedbCache.quote("gmm/${modelRoot}")} \
        --metadata metadata.tsv --source-scores scores.tsv \
        --clinvar clinvar.tsv --gwas gwas.tsv --out validation \
        --formats ${formatArguments} --dpi ${MavedbCache.quote(dpi)} ${noPlots}
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output validation ${requiredOutputs} > result.json
    """
}

workflow MAVEDB_VALIDATION {
    take:
    engine
    processed
    merged
    references

    main:
    decisions = processed.combine(merged).combine(references).map { canonical, model, reference ->
        MavedbCache.resolve([
            engine: engine, stage: 'validation',
            inputs: [
                process: canonical.path,
                merge: model.path,
                references: reference.path,
            ]
        ])
    }
    reused = decisions.filter { it.path != null }
    jobs = decisions.filter { it.spec != null }.map { job ->
        tuple(job,
              file(job.inputs.merge, checkIfExists: true),
              file("${job.inputs.process}/scores.tsv", checkIfExists: true),
              file("${job.inputs.process}/metadata.tsv", checkIfExists: true),
              file("${job.inputs.references}/variants/mave_clinvar_variants.tsv", checkIfExists: true),
              file("${job.inputs.references}/variants/mave_gwas_variants.tsv", checkIfExists: true))
    }
    MAVEDB_GMM_VALIDATION(jobs)
    computed = MAVEDB_GMM_VALIDATION.out.result.map { MavedbCache.result(it) }
    ready = reused.mix(computed)

    emit:
    result = ready
    bundle = ready.map { file(it.path, checkIfExists: true) }
}
