process CHECK_MAVEDB_UNIPROT_RESIDUES {
    tag 'UniProt reference residue check'
    label 'mavedb_processing'
    cache false

    input:
    tuple val(job), path(scores, stageAs: 'scores.tsv'),
          path(metadata, stageAs: 'metadata.tsv'), path(residues, stageAs: 'residues.tsv')

    output:
    path 'result.json', emit: result

    script:
    """
    export PYTHONDONTWRITEBYTECODE=1
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.scripts + '/check_uniprot_residues.py')} \
        --scores scores.tsv --metadata metadata.tsv --residues residues.tsv --output residue_check
    ${MavedbCache.quote(job.engine.python)} ${MavedbCache.quote(job.engine.runtime)} register \
        --spec ${MavedbCache.quote(job.spec)} --output residue_check \
        --require residue_checks.tsv --require summary.tsv \
        --require dataset_summary.tsv --require applied_offset_corrections.tsv > result.json
    """
}

workflow MAVEDB_RESIDUE_CHECK {
    take:
    engine
    processed
    residues

    main:
    decisions = processed.map { canonical ->
        MavedbCache.resolve([
            engine: engine, stage: 'residue_check',
            inputs: [process: canonical.path, residues: residues]
        ])
    }
    reused = decisions.filter { it.path != null }
    jobs = decisions.filter { it.spec != null }.map { job ->
        [job, file("${job.inputs.process}/scores.tsv", checkIfExists: true),
         file("${job.inputs.process}/metadata.tsv", checkIfExists: true),
         file(job.inputs.residues, checkIfExists: true)]
    }
    CHECK_MAVEDB_UNIPROT_RESIDUES(jobs)
    computed = CHECK_MAVEDB_UNIPROT_RESIDUES.out.result.map { MavedbCache.result(it) }
    ready = reused.mix(computed)

    emit:
    result = ready
}
