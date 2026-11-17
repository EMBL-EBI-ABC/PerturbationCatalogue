process PACKAGE_GMM_OUTPUTS {
    tag 'Package and merge GMM outputs'
    cache false

    input:
    tuple val(job), path(collected, stageAs: 'collected')

    output:
    path 'gmm', emit: bundle
    path 'result.json', emit: result

    script:
    """
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/package_gmm_outputs.py')} \
        --bundle collected --out gmm
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output gmm --require completion.json --require batch_plan.json > result.json
    """
}

workflow MAVEDB_MERGE {
    take:
    engine
    fitted

    main:
    decisions = fitted.map {
        MavedbCache.resolve([engine: engine, stage: 'merge', inputs: [gmm: it.path]])
    }
    PACKAGE_GMM_OUTPUTS(decisions.filter { it.spec != null }.map { job -> [job, file(job.inputs.gmm)] })
    ready = decisions.filter { it.path != null }.mix(PACKAGE_GMM_OUTPUTS.out.result.map { MavedbCache.result(it) })

    emit:
    result = ready
    bundle = ready.map { file(it.path) }
}
