# Norman 2019

`generate_inputs.py` reads the committed author Table S2 copy and creates the
18 bp GBC feature reference plus the eight gemgroup sample rows. The feature
labels are already in the pipeline's standard form: target pairs use `;`, and
negative-control-only pairs use `non-targeting`.

The author H5AD has paired identities such as `AHR_NegCtrl0__AHR_NegCtrl0`
and uses `NegCtrl` names or blank guide IDs. `prepare_curated.py` derives the
standard `obs["perturbation"]` column (`AHR`, `AHR;FEV`, or `non-targeting`)
without loading the expression matrix, so the shared comparison code receives
the same label vocabulary as the pipeline output.

`gene_ensg.tsv` records the four historical Norman symbols whose current
reference annotation uses updated names; those labels are emitted with their
Ensembl IDs.
