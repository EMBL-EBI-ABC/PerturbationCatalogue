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

2. **Full, Unmodified Sentences & Exhaustive Search (Strict Rule):**
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
