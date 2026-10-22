# 15: Produce normalized assay effects on Slurm

**What to build:** An analyst can run the required repository-owned R/Stan assay methods on Slurm and obtain reproducible normalized effects, diagnostics and traceable BigQuery results.

**Blocked by:** 05: Map variants and report mapping QC; 06: Load internal BigQuery products; 14: Produce basic assay QC.

**Status:** ready-for-agent

- [ ] Required analysis methods and resources are maintained in this repository, with attribution and licences preserved and no sibling-repository or scheduler-wrapper dependency.
- [ ] Nextflow stages declared GCS inputs for reuse, executes score-set-level tasks through Slurm in a defined environment, and writes effects and diagnostics to GCS before selected results enter BigQuery through `LOAD_DATA`.
- [ ] Analysis eligibility uses explicit parsing and mapping requirements; insufficient data, model convergence and method failures produce analysis QC with excluded counts and reasons.
- [ ] Normalized effects, raw scores and fit diagnostics remain distinct, and derived substitutions link back to every contributing original observation, including multiple nucleotide variants for one protein substitution.
