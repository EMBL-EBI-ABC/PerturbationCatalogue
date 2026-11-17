process MAVEDB_SUMMARY_REPORTS {
    tag 'Per-study GMM summary reports'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(merged, stageAs: 'gmm'), path(metadata, stageAs: 'metadata.tsv')

    output:
    path 'reports', emit: bundle
    path 'result.json', emit: result

    script:
    def p = MavedbCache.request(job)
    def modelRoot = (p.gmm_profile == 'smoke' ? 'smoke/' : '') + "maves/${p.gmm_model}"
    """
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote("${job.engine.scripts}/summary_plot_batch.py")} \
        --gmm-dir ${MavedbCache.quote("gmm/${modelRoot}")} --meta metadata.tsv \
        --model ${MavedbCache.quote(p.gmm_model)} --out reports
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output reports > result.json
    """
}

workflow MAVEDB_REPORTS {
    take:
    engine
    processed
    merged

    main:
    decisions = processed.combine(merged).map { canonical, model ->
        MavedbCache.resolve([
            engine: engine, stage: 'reports',
            inputs: [process: canonical.path, merge: model.path]
        ])
    }
    reused = decisions.filter { it.path != null }
    jobs = decisions.filter { it.spec != null }.map { job ->
        [job, file(job.inputs.merge, checkIfExists: true),
         file("${job.inputs.process}/metadata.tsv", checkIfExists: true)]
    }
    MAVEDB_SUMMARY_REPORTS(jobs)
    computed = MAVEDB_SUMMARY_REPORTS.out.result.map { MavedbCache.result(it) }
    ready = reused.mix(computed)

    emit:
    result = ready
    bundle = ready.map { file(it.path, checkIfExists: true) }
}
