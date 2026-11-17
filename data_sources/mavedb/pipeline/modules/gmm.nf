/** CLI options for the frozen GMM attempt; shared by planning and fitting. */
def gmmArguments(Map request) {
    def options = [
        profile: 'profile', model: 'model', variants_per_batch: 'variants_per_batch',
        seed: 'seed', chains: 'chains', parallel_chains: 'parallel-chains',
        iter_warmup: 'iter-warmup', iter_sampling: 'iter-sampling',
        max_studies: 'max-studies', max_study_rows: 'max-study-rows',
        max_fit_rows: 'max-fit-rows', memory_budget_mb: 'memory-budget-mb',
    ]
    def argv = []
    options.each { key, option ->
        def value = request["gmm_${key}"]
        if (value != null) argv.addAll(["--${option}", value])
    }
    def studies = request.gmm_study_ids instanceof List
        ? request.gmm_study_ids
        : (request.gmm_study_ids ?: '').toString().tokenize(',').collect { it.trim() }.findAll { it }
    studies.each { study -> argv.addAll(['--study-id', study]) }
    return MavedbCache.arguments(argv)
}

process PLAN_GMM_BATCHES {
    tag 'Freeze GMM batch plan'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(scores, stageAs: 'scores.tsv'),
          path(metadata, stageAs: 'metadata.tsv'),
          path(fitter, stageAs: 'fit_gmm_batch.py'),
          path(modules, stageAs: 'modules/*'),
          path(stan_sources, stageAs: 'stan_sources')

    output:
    path 'plan_bundle', emit: bundle
    path 'result.json', emit: result

    script:
    def arguments = gmmArguments(MavedbCache.request(job))
    """
    export PYTHONDONTWRITEBYTECODE=1
    mkdir -p plan_bundle/code/data_sources/mavedb/scripts plan_bundle/code/data_sources/mavedb/analysis/stan
    cp fit_gmm_batch.py plan_bundle/code/data_sources/mavedb/scripts/
    cp modules/*.py plan_bundle/code/data_sources/mavedb/analysis/
    cp stan_sources/*.stan plan_bundle/code/data_sources/mavedb/analysis/stan/
    ${MavedbCache.quote(job.engine.python)} plan_bundle/code/data_sources/mavedb/scripts/fit_gmm_batch.py \
        --scores scores.tsv --meta metadata.tsv --list-batches \
        --summary-json plan_bundle/batch_plan.json ${arguments}
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output plan_bundle --require batch_plan.json --require code/data_sources/mavedb/scripts/fit_gmm_batch.py > result.json
    """
}

process FIT_GMM_BATCH {
    tag { "GMM batch ${job.index}" }
    label 'gmm_fit'
    cache false

    input:
    tuple val(job), path(frozen, stageAs: 'frozen'),
          path(scores, stageAs: 'scores.tsv'), path(metadata, stageAs: 'metadata.tsv')

    output:
    path 'batch', emit: bundle
    path 'result.json', emit: result

    script:
    def request = MavedbCache.request(job)
    def arguments = gmmArguments(request)
    def failFlag = request.gmm_fail_on_study_error ? '--fail-on-study-error' : ''
    """
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export MAVEDB_GMM_CACHE_DIR=${MavedbCache.quote(request.gmm_cache_dir)}
    mkdir -p tmp
    export TMPDIR="\$PWD/tmp"
    ${MavedbCache.quote(job.engine.python)} frozen/code/data_sources/mavedb/scripts/fit_gmm_batch.py \
        --scores scores.tsv --meta metadata.tsv --index ${job.index} \
        --out batch --expected-plan frozen/batch_plan.json \
        --summary-json batch/summary.json ${arguments} ${failFlag}
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output batch --require summary.json > result.json
    """
}

process COLLECT_GMM_BATCHES {
    tag 'Verify all planned GMM batches'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(frozen, stageAs: 'frozen'),
          path(batches, stageAs: 'batch-inputs/*')

    output:
    path 'gmm', emit: bundle
    path 'result.json', emit: result

    script:
    def request = MavedbCache.request(job)
    def failFlag = request.gmm_fail_on_study_error ? '--fail-on-study-error' : ''
    def batchArguments = MavedbCache.arguments(batches instanceof List ? batches : [batches])
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${MavedbCache.quote(job.engine.python)} frozen/code/data_sources/mavedb/scripts/fit_gmm_batch.py \
        --collect-plan frozen/batch_plan.json --batch-results ${batchArguments} \
        --out gmm ${failFlag}
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output gmm --require completion.json --require batch_plan.json > result.json
    """
}

workflow MAVEDB_GMM {
    take:
    engine
    processed

    main:
    planning = processed.map { canonical ->
        MavedbCache.resolve([engine: engine, stage: 'gmm-plan', inputs: [process: canonical.path]])
    }
    reusedPlans = planning.filter { it.path != null }
    PLAN_GMM_BATCHES(planning.filter { it.spec != null }.map { job ->
        def scripts = job.engine.scripts
        def modules = ['gmm.py', '_gmm_batch.py', '_gmm_inputs.py', '_gmm_stan.py', '_gmm_pipeline.py'].collect {
            file("${scripts}/../analysis/${it}", checkIfExists: true)
        }
        tuple(job,
              file("${job.inputs.process}/scores.tsv", checkIfExists: true),
              file("${job.inputs.process}/metadata.tsv", checkIfExists: true),
              file("${scripts}/fit_gmm_batch.py", checkIfExists: true),
              modules, file("${scripts}/../analysis/stan", checkIfExists: true))
    })
    plans = reusedPlans.mix(PLAN_GMM_BATCHES.out.result.map { MavedbCache.result(it) })

    fitting = plans.combine(processed).flatMap { plan, canonical ->
        if (plan.collected) return []
        def frozen = new groovy.json.JsonSlurper().parseText(file("${plan.path}/batch_plan.json").text)
        if (!frozen.batch_count) error 'GMM selection contains no eligible MAVE studies'
        frozen.assignments.collect { it[0] }.unique().sort().collect { index ->
            MavedbCache.resolve([engine: engine, stage: 'gmm-batch', index: index,
                                 inputs: ['gmm-plan': plan.path, process: canonical.path]])
        }
    }
    reusedBatches = fitting.filter { it.path != null }
    FIT_GMM_BATCH(fitting.filter { it.spec != null }.map { job ->
        tuple(job, file(job.inputs['gmm-plan'], checkIfExists: true),
              file("${job.inputs.process}/scores.tsv", checkIfExists: true),
              file("${job.inputs.process}/metadata.tsv", checkIfExists: true))
    })
    batches = reusedBatches.mix(FIT_GMM_BATCH.out.result.map { MavedbCache.result(it) })

    collecting = plans.filter { !it.collected }.combine(processed)
        .combine(batches.collect().map { [batches: it] }).map { plan, canonical, completed ->
            def inputs = ['gmm-plan': plan.path, process: canonical.path]
            completed.batches.eachWithIndex { batch, index -> inputs["batch-${index}"] = batch.path }
            MavedbCache.resolve([engine: engine, stage: 'gmm', inputs: inputs])
        }
    reusedCollections = collecting.filter { it.path != null }
    COLLECT_GMM_BATCHES(collecting.filter { it.spec != null }.map { job ->
        def batchPaths = job.inputs.findAll { key, value -> key.startsWith('batch-') }
            .collect { key, value -> file(value, checkIfExists: true) }
        tuple(job, file(job.inputs['gmm-plan'], checkIfExists: true), batchPaths)
    })
    result_ch = plans.filter { it.collected }.map { [stage: 'gmm', path: it.path] }
        .mix(reusedCollections, COLLECT_GMM_BATCHES.out.result.map { MavedbCache.result(it) })

    emit:
    result = result_ch
    bundle = result_ch.map { file(it.path) }
    plan = plans.map { file(it.path).resolve('batch_plan.json') }
}
