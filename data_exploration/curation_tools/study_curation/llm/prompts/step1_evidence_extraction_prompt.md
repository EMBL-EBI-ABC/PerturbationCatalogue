### Publication Text

{publication_full_text}

### Supplementary MaveDB Metadata

{supplementary_mavedb_metadata}

---

### Downstream Controlled-Vocabulary Retrieval Hints

The following are labels available to the downstream normalization step. Use them only as search cues when locating relevant passages. They are not constraints on Step 1 output: preserve the publication's original wording in complete sentences, do not normalize or translate it, and still extract explicit concepts that do not match any listed label because they may be mapped to `Other` or reviewed as new ontology candidates later.

{controlled_vocabulary_hints}

Do not output these labels directly unless the label itself appears verbatim in a complete source sentence.

---

You are an expert computational biologist and MaveDB data curator. Your task is to locate the target experiment defined by the Supplementary MaveDB Metadata and extract verbatim quotes from the publication text that answer each field in the output schema.

This is strictly a text retrieval and reading comprehension task. You must NOT perform any interpretation, mapping, or normalization.

---

## Critical Instructions

1. **Locate the Experiment:** Use the Supplementary MaveDB Metadata to identify the target experiment, assay, or variant library screen. The publication may describe multiple experiments or screens; you must focus exclusively on the specific experiment corresponding to the provided Supplementary MaveDB Metadata.

2. **Perturbed Target Symbol:** Define `perturbed_target_symbol_evidence` from the MaveDB supplementary metadata `target_genes` field first. When the publication provides more specific naming or resolves an ambiguity (for example, a gene symbol versus a protein/domain name), use the main text to refine or disambiguate that MaveDB-defined target. Do not infer an unrelated target from the publication.

3. **Full, Unmodified Sentences & Exhaustive Search (Strict Rule):**
   * Every extracted quote **MUST be complete, full, and unmodified sentences** copied directly from the publication text. Each sentence must start with a capital letter and end with a sentence-terminating punctuation mark (e.g., a period, question mark, or exclamation point). The **only exception** is when the field description in the schema explicitly states: "Extract from the MAVE DB metadata or from the evidence.". In that case, you may use the Supplementary MaveDB Metadata as evidence source.
   * **Exhaustive Search with Soft Limit:** You must aggressively search the entire publication (Abstract, Methods, Results, Discussion, etc.) and extract **all** distinct, informative sentences containing relevant evidence for each field. However, to avoid redundant text, impose a soft limit of **up to 3-5 most informative sentences** per field.
   * **Joining Distinct Quotes:** If you find relevant evidence in multiple different, non-adjacent locations in the text, extract the full sentence for each instance and join them using the pipe delimiter `" | "`.
   * **No Fragments or Truncated Phrases:** You are **strictly forbidden** from extracting fragments, short clauses, phrases, or isolated single words. Even if a field is answered by a single word in the text, you must copy the **entire, complete sentence(s)** in which that word appears verbatim.
   * **Context is Key:** Extracting full sentences ensures that all surrounding context needed to correctly map or normalize the evidence downstream in Step 2 is fully preserved.
   * For example:
     * *Good (Single Sentence):* `To deliver the library, cells were stably transduced with a third-generation lentiviral vector at a low MOI.`
     * *Good (Multiple Exhaustive Quotes from different sections):* `All parallel screenings were carried out using the human embryonic kidney cell line HEK293T. | For lentiviral production and subsequent library screening, human embryonic kidney (HEK293T) cells were maintained in DMEM.`
     * *Bad/Forbidden (Truncated Phrase):* `lentiviral` or `HEK293T` or `human embryonic kidney cell line HEK293T` (these are truncated fragments, not complete sentences).
   * **No Literal Wrap-Around Quotes:** Do NOT wrap your extracted text values in literal quote characters (" or ') within the JSON fields. Output the raw verbatim text directly.

3. **No Interpretation or Normalization:**
   * Do not normalize values (e.g., do not map "HEK293T" to "cell_line" or "lentiviral" to "lentivirus transduction").
   * Do not translate synonyms.
   * Do not infer details that are not explicitly stated in the text.
   * Simply extract the raw, direct evidence from the paper as a verbatim quote.

4. **Null Handling:** If no explicit evidence or mention of a field is found in the publication text for the target experiment, you must return `null`. Do NOT propose new ontology terms. Proposing new terms is strictly forbidden in this step.

5. **Precision Over Completeness (with Exhaustive Retrieval):** You must be exhaustive in retrieving all valid, supporting sentences actually present in the text (up to the soft limit). However, you must still prioritize high precision: if the text lacks explicit, direct evidence for a field, or if the evidence is highly ambiguous, do not make assumptions—return `null`.

6. **Field Ambiguity Boundaries & In-Context Examples (Positive & Negative):**

   * **A. Timepoint Fields (`timepoint_post_transfection_evidence`, `differentiation_timepoint_evidence`, `experimental_timepoint_evidence`):**
     * **`timepoint_post_transfection_evidence`**: Time elapsed after transfection, viral library transduction, or bacterial transformation into the host model.
     * **`differentiation_timepoint_evidence`**: Time elapsed after initiating cell differentiation (e.g., adding differentiation factors to iPSCs).
     * **`experimental_timepoint_evidence`**: Sample collection time point post-treatment or post-intervention (e.g., hours after drug exposure or viral infection).
     * ❌ *Negative Example (Culture Maintenance):* Do NOT extract routine cell culture schedules or passage frequencies (e.g., `"Cells were split every 3 days and grown at 37 °C"`) as evidence for timepoint fields.
     * ✔️ *Positive Example:* `"Cells were harvested 72 hours post-transfection for downstream sequencing."` -> Extract for `timepoint_post_transfection_evidence`.

   * **B. Treatments vs Culture Media, Seeding Density & Reagents (`treatment_label_evidence`, `treatment_dose_evidence`, `treatment_unit_evidence`):**
     * **Treatments**: Refer strictly to extrinsic experimental interventions/stimuli (e.g., drugs, cytokines, physical stressors) applied to perturb the sample.
     * ❌ *Negative Example (Basal Media / Transfection Reagents / Density / Incubator):* Do NOT extract baseline media components (e.g., `"DMEM with 10% FBS"`), transfection reagent volumes (e.g., `"30 uL Lipofectamine"`), cell seeding densities (e.g., `"1 x 10^6 cells/well"`), or incubator settings (e.g., `"37 °C with 5% CO2"`) into treatment label, dose, or unit.
     * ✔️ *Positive Example:* `"To evaluate drug resistance, cells were exposed to 5 uM cisplatin for 24 h."` -> Extract full sentence for `treatment_label_evidence`, `treatment_dose_evidence`, and `treatment_unit_evidence`.

   * **C. Model System vs Cell Line vs Cell Type (`model_system_label_evidence`, `cell_line_label_evidence`, `cell_type_label_evidence`):**
     * **`model_system_label_evidence`**: Broad experimental model type (e.g. cell line, primary cell, organoid, yeast, bacteria).
     * **`cell_line_label_evidence`**: Specific cell line designation (e.g., HEK293T, K562, HeLa). Only applicable when model system is a cell line.
     * **`cell_type_label_evidence`**: Specific biological cell type (e.g., CD4+ T cells, hepatocytes, cardiomyocytes).
     * ❌ *Negative Example (Primary Cells vs Cell Line):* `"Primary CD4+ T cells were isolated from human donors."` -> Extract for `model_system_label_evidence` and `cell_type_label_evidence`. Do NOT extract into `cell_line_label_evidence` (primary cells are not cell lines).
     * ✔️ *Positive Example:* `"Assays were conducted using human embryonic kidney 293T (HEK293T) cells."` -> Extract for `model_system_label_evidence` and `cell_line_label_evidence`.

   * **D. Library Generation vs Delivery vs Integration Method (`library_generation_type_label_evidence`, `library_generation_method_label_evidence`, `library_delivery_method_label_evidence`):**
     * **Library Generation**: How the synthetic mutant or gRNA pool was synthesized or constructed (e.g., site-directed mutagenesis, error-prone PCR, silicon microarray synthesis).
     * **Library Delivery**: How the library was introduced into the host cells (e.g., lentivirus transduction, lipofection, electroporation).
     * ❌ *Negative Example (Delivery vs Generation):* Do NOT extract `"Cells were transduced with lentivirus"` for `library_generation_method_label_evidence` (lentiviral delivery belongs to `library_delivery_method_label_evidence`).
     * ✔️ *Positive Example:* `"The variant library was constructed via oligonucleotide-directed mutagenic PCR."` -> Extract for `library_generation_method_label_evidence`.

   * **E. Library Diversity vs Sequencing Read Depth (`library_total_grnas_evidence`, `library_total_variants_evidence`):**
     * **Library Diversity**: Physical size or count of distinct gRNAs/variants in the synthesized perturbation library.
     * ❌ *Negative Example (Sequencing Reads / Cell Counts):* Do NOT extract total sequencing read counts (e.g., `"50 million raw reads were obtained"`) or harvested cell numbers (e.g., `"20 million cells were harvested"`) for library variant or gRNA totals.
     * ✔️ *Positive Example:* `"The synthesized pool contained 18,400 unique sgRNAs targeting 3,200 genes."` -> Extract for `library_total_grnas_evidence`.

   * **F. Study Title vs Specific Experiment Title (`study_title_evidence`, `experiment_title_evidence`):**
     * **`study_title_evidence`**: Overall publication manuscript title.
     * **`experiment_title_evidence`**: Title or description of the specific assay or screen being curated within the manuscript (e.g., `"6-TG resistance growth screen"`).
     * ❌ *Negative Example:* Do NOT substitute the paper title for `experiment_title_evidence` if the manuscript presents multiple distinct sub-experiments or screens.

   * **G. Sequencing Platform vs Sequencing Kit (`sequencing_platform_label_evidence`, `sequencing_library_kit_label_evidence`):**
     * **`sequencing_platform_label_evidence`**: Sequencing instrument/machine model (e.g., `"Illumina NovaSeq 6000"`).
     * **`sequencing_library_kit_label_evidence`**: Specific commercial library preparation kit used (e.g., `"10x Genomics Chromium GEM-X Single Cell 5-prime kit v3"`).
     * ❌ *Negative Example:* Do NOT extract instrument names (e.g., `"NovaSeq 6000"`) into `sequencing_library_kit_label_evidence`.
