process MAVEDB_BASIC_QC {
    tag 'Basic source-score QC'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(scores, stageAs: 'scores.tsv'),
          path(metadata, stageAs: 'metadata.tsv'), path(singles, stageAs: 'singles.tsv')

    output:
    path 'qc', emit: bundle
    path 'result.json', emit: result

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote("${job.engine.scripts}/../analysis/mave_qc.py")} \
        --scores scores.tsv --metadata metadata.tsv --singles singles.tsv --output qc
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output qc > result.json
    """
}

workflow MAVEDB_QC {
    take:
    engine
    processed
    tables

    main:
    decisions = processed.combine(tables).map { canonical, compatibility ->
        MavedbCache.resolve([
            engine: engine, stage: 'qc',
            inputs: [process: canonical.path, tables: compatibility.path]
        ])
    }
    reused = decisions.filter { it.path != null }
    jobs = decisions.filter { it.spec != null }.map { job ->
        [job, file("${job.inputs.process}/scores.tsv", checkIfExists: true),
         file("${job.inputs.process}/metadata.tsv", checkIfExists: true),
         file("${job.inputs.tables}/maves/singles.tsv", checkIfExists: true)]
    }
    MAVEDB_BASIC_QC(jobs)
    computed = MAVEDB_BASIC_QC.out.result.map { MavedbCache.result(it) }
    ready = reused.mix(computed)

    emit:
    result = ready
    bundle = ready.map { file(it.path, checkIfExists: true) }
}
