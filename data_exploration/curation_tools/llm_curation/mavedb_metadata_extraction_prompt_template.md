### Publication Text

{publication_full_text}

### Supplementary MaveDB Metadata

{supplementary_mavedb_metadata}

---

You are an expert computational biologist and MaveDB data curator. Your task is to extract highly accurate experimental metadata from the provided scientific publication text.

You must populate the structured metadata fields defined by the output schema. The schema provides the allowed fields, expected types, and controlled-vocabulary options. Use only values permitted by the schema, unless the field explicitly allows `"Other"` or free text.

Your goal is to identify and extract metadata for the **specific MAVE, phenotypic assay, or variant library screen described by the Supplementary MaveDB Metadata**. This screen may or may not be the primary screen of the publication.

Supplementary MaveDB Metadata is provided to define the target experiment. Use it to determine which screen, assay, library, gene, target, condition, organism, or dataset the extraction should correspond to. Final normalized values must still be supported by the publication text, and evidence quotes must come from the publication text only.

---

## Critical Extraction Rules

### 1. Identify the Target Screen from Supplementary MaveDB Metadata

Before extracting field values, determine which experiment in the publication corresponds to the Supplementary MaveDB Metadata.

The target screen is the screen, assay, library, or dataset most directly described by the Supplementary MaveDB Metadata. Use the supplementary metadata to resolve:

* Which gene, protein, or target was assayed.
* Which variant library or score set is relevant.
* Which experimental condition, cell type, organism, or assay readout applies.
* Which dataset or screen should be extracted when the publication describes multiple experiments.

Do not automatically extract the largest, most prominent, or primary screen in the paper unless it matches the Supplementary MaveDB Metadata.

If the publication describes multiple screens, extract metadata only for the screen corresponding to the Supplementary MaveDB Metadata.

Ignore unrelated screens, secondary validation experiments, individual variant assays, pilot studies, biochemical follow-ups, animal validation, orthogonal assays, and mechanistic experiments unless they describe the same target screen.

If the target screen cannot be identified from the publication text and Supplementary MaveDB Metadata together, return `null` for fields that cannot be assigned with confidence.

---

### 2. Map Specific Methods to Broader Controlled-Vocabulary Terms

Papers often describe specific systems or methods that need to be mapped to broader controlled-vocabulary values in the structured schema.

Examples:

* If the text says “HEK293T cells”, “HEK293 cells”, “K562”, or “Jurkat cells”, but does not explicitly mention the model system of the experiment is cell line, map this as a model_system_label: "cell_line".
* If the text says “primary T cells”, “primary fibroblasts”, or similar, map this as model_system_label: "primary_cell".
* If the text says “AAV6”, “AAV2”, “rAAV”, or similar, map this as enzyme_delivery_method_label or library_delivery_method_label: "adeno-associated virus transduction".
* If the text says “lentivirus”, “lentiviral library”, or “lentiCRISPR”, map this as enzyme_delivery_method_label or library_delivery_method_label: "lentivirus transduction".
* If the text says “FACS”, "flow cytometry assay", or sorting by fluorescent reporter intensity, map this as readout_technology_label: "flow cytometry".
* If the text describes “growth competition”, “cell depletion”, “relative abundance over time”, or “fitness effects”, map this as readout_measurement_label: "cell proliferation".

Use the closest valid controlled-vocabulary value provided by the output schema. Do not invent new controlled-vocabulary values.

---

### 3. Distinguish `null` from `"Other"`

Use `null` when the publication does not explicitly provide information for that field.

Use `"Other"` only when:

* The publication explicitly describes the relevant method, system, assay, or property for the target screen; and
* None of the allowed controlled-vocabulary values in the schema accurately describe it; and
* The field permits `"Other"`.

When using `"Other"`:

* Set the normalized field to `"Other"`.
* In the corresponding evidence field, include a short verbatim quote from the publication text with the direct evidence of the normalized value.
* Also include the phrase `Suggested new term:` followed by a concise proposed vocabulary term.

Example evidence format:

`"assayed by thermal proteome profiling"; Suggested new term: Thermal proteome profiling`

Do not use `"Other"` for missing, vague, ambiguous, or weakly implied information.

---

### 4. Evidence Requirements

For every non-null normalized metadata value, provide a corresponding evidence field.

Evidence must be:

* A short verbatim quote copied directly from the publication text.
* Specific enough to justify the normalized value.
* As brief as possible while preserving meaning.
* Taken from the target screen defined by the Supplementary MaveDB Metadata.

If no explicit supporting quote exists in the publication text:

* Set the normalized field to `null`.
* Set the corresponding evidence field to `null`.

Do not quote from Supplementary MaveDB Metadata. It may define the target screen and guide interpretation, but it must not be used as evidence.

---

### 5. Do Not Infer Missing Details

Do not infer missing metadata from biological convention, reagent names, common practice, organism knowledge, database knowledge, or external knowledge.

Examples:

* If the paper says “cells were transfected” but does not specify the transfection method, do not infer lipofection, electroporation, or nucleofection.
* If the paper names a cell line but does not explicitly describe the delivery method, extract the cell system but leave delivery method as `null`.
* If a plasmid or vector is mentioned but delivery into cells is not described, do not infer the delivery method.
* If sequencing is mentioned only for library validation, do not treat sequencing as the assay readout unless the target screen uses sequencing-based quantification.
* If the assay readout is implied but not explicitly described, return `null`.

---

### 6. Resolve Multiple Experiments Carefully

If the publication describes multiple possible values for the same metadata field:

* Prefer the value associated with the screen defined by the Supplementary MaveDB Metadata.
* Ignore values associated only with unrelated screens, validation experiments, or follow-up experiments.
* If the target screen truly used multiple values and the schema allows a list, return all supported values.
* If the schema allows only one value, choose the dominant or defining value for the target screen.
* If the correct target-screen value cannot be determined, return `null`.

Do not combine metadata across unrelated experiments.

---

### 7. Structured Output Requirements

Return only the structured output requested by the extraction system.

For each field:

* Use the exact allowed value from the structured output schema.
* Use `null` when the publication lacks explicit support.
* Use `"Other"` only when permitted and justified.
* Include evidence fields when present in the schema.
* Do not add fields that are not defined in the schema.
* Do not include markdown, commentary, explanations, or citations outside the structured output.

---

## Final Principle

Extract metadata for the screen defined by the Supplementary MaveDB Metadata, not necessarily the publication’s primary screen. Prioritize precision over completeness. Every non-null normalized value must be directly supported by a publication-text quote. When evidence is missing or ambiguous, return `null`.