def quoteGmmArgument(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process PLAN_GMM_BATCHES {
    tag 'Freeze GMM batch plan'
    publishDir "${params.outdir}/gmm-plan", mode: 'copy', pattern: 'batch_plan.json', overwrite: true

    input:
    path processed, stageAs: 'processed'
    path fitter, stageAs: 'fit_gmm_batch.py'
    path modules, stageAs: 'modules/*'
    path stan_sources, stageAs: 'stan_sources'
    val arguments
    env 'python_executable'

    output:
    path 'batch_plan.json', emit: plan
    path 'code', emit: code

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    mkdir -p code/data_sources/mavedb/scripts code/data_sources/mavedb/analysis/stan
    cp fit_gmm_batch.py code/data_sources/mavedb/scripts/
    cp modules/*.py code/data_sources/mavedb/analysis/
    cp stan_sources/*.stan code/data_sources/mavedb/analysis/stan/
    "\$python_executable" code/data_sources/mavedb/scripts/fit_gmm_batch.py \
        --scores processed/scores.tsv --meta processed/metadata.tsv \
        --list-batches --summary-json batch_plan.json ${arguments}
    """
}

process FIT_GMM_BATCH {
    tag "GMM batch ${index}"
    label 'gmm_fit'

    input:
    tuple val(index), path(plan, stageAs: 'batch_plan.json'),
          path(code, stageAs: 'code'), path(processed, stageAs: 'processed')
    val arguments
    val fail_flag
    env 'python_executable'
    env 'MAVEDB_GMM_CACHE_DIR'

    output:
    path "batch-${index}", emit: bundle

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    mkdir -p tmp
    export TMPDIR="\$PWD/tmp"
    "\$python_executable" code/data_sources/mavedb/scripts/fit_gmm_batch.py \
        --scores processed/scores.tsv --meta processed/metadata.tsv \
        --index ${index} --out batch-${index} --expected-plan batch_plan.json \
        --summary-json batch-${index}/summary.json ${fail_flag}${arguments}
    """
}

process COLLECT_GMM_BATCHES {
    tag 'Verify all planned GMM batches'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path plan, stageAs: 'batch_plan.json'
    path batches
    path code, stageAs: 'code'
    val fail_flag
    env 'python_executable'

    output:
    path 'gmm', emit: bundle

    script:
    def batchArguments = batches.collect { quoteGmmArgument(it) }.join(' ')
    """
    export PYTHONDONTWRITEBYTECODE=1
    "\$python_executable" code/data_sources/mavedb/scripts/fit_gmm_batch.py \
        --collect-plan batch_plan.json --batch-results ${batchArguments} ${fail_flag}--out gmm
    """
}

workflow MAVEDB_GMM {
    take:
    processed

    main:
    if (params.gmm_max_parallel.toInteger() < 1) error '--gmm_max_parallel must be positive'
    options = [
        profile: params.gmm_profile, model: params.gmm_model,
        variants_per_batch: params.gmm_variants_per_batch, seed: params.gmm_seed,
        chains: params.gmm_chains, 'parallel-chains': params.gmm_parallel_chains,
        'iter-warmup': params.gmm_iter_warmup, 'iter-sampling': params.gmm_iter_sampling,
        'max-studies': params.gmm_max_studies, 'max-study-rows': params.gmm_max_study_rows,
        'max-fit-rows': params.gmm_max_fit_rows, 'memory-budget-mb': params.gmm_memory_budget_mb,
    ]
    argv = options.findAll { key, value -> value != null }.collectMany { key, value ->
        ["--${key}", value]
    }
    studies = params.gmm_study_ids instanceof List ? params.gmm_study_ids : params.gmm_study_ids.toString().tokenize(',')
    studies.each { study -> argv.addAll(['--study-id', study]) }
    arguments = argv.collect { quoteGmmArgument(it) }.join(' ')
    failFlag = params.gmm_fail_on_study_error ? '--fail-on-study-error ' : ''
    sourceFiles = ['gmm.py', '_gmm_batch.py', '_gmm_inputs.py', '_gmm_stan.py', '_gmm_pipeline.py'].collect {
        file("${projectDir}/../analysis/${it}", checkIfExists: true)
    }
    PLAN_GMM_BATCHES(
        processed, file("${projectDir}/../scripts/fit_gmm_batch.py"), sourceFiles,
        file("${projectDir}/../analysis/stan"), arguments, params.python
    )
    jobs = PLAN_GMM_BATCHES.out.plan.combine(PLAN_GMM_BATCHES.out.code).combine(processed).flatMap {
        plan, code, inputs ->
        def parsed = new groovy.json.JsonSlurper().parseText(plan.text)
        if (!parsed.batch_count) error 'GMM selection contains no eligible MAVE studies'
        def indices = parsed.assignments.collect { it[0] }.unique().sort()
        indices.collect { index -> [index, plan, code, inputs] }
    }
    FIT_GMM_BATCH(jobs, arguments, failFlag, params.python, file(params.gmm_cache_dir).toAbsolutePath().toString())
    COLLECT_GMM_BATCHES(
        PLAN_GMM_BATCHES.out.plan, FIT_GMM_BATCH.out.bundle.collect(), PLAN_GMM_BATCHES.out.code, failFlag, params.python
    )

    emit:
    bundle = COLLECT_GMM_BATCHES.out.bundle
    plan = PLAN_GMM_BATCHES.out.plan
}
