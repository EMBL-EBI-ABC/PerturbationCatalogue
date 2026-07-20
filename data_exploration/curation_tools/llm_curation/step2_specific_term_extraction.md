### Step 1 Evidence (JSON)

{step1_evidence}

### Supplementary MaveDB Metadata

{supplementary_mavedb_metadata}

---

You are an expert computational biologist and MaveDB data curator. Your task is to review the verbatim publication evidence extracted in Step 1, cross-reference it with any Supplementary MaveDB Metadata, and normalize the evidence into the controlled vocabularies defined by the target output schema.

This is strictly an evidence mapping and classification task. You must normalize the verbatim quotes into standard ontology terms / labels.

---

## Critical Instructions

1. **Mapping Verbatim Evidence to Schema:**
   - Look at the `*_evidence` fields in the input Step 1 JSON.
   - Map each piece of evidence to its corresponding normalized field in the `SpecificTermExtractionSchema` (e.g. map `model_system_label_evidence` or any related context to `model_system_label`).
   - If `dataset_id` or basic metadata (like `study_title`, `study_uri`, etc.) is present in the Step 1 JSON or Supplementary MaveDB Metadata, format and populate them.

2. **Synonym Resolution & Case Insensitivity:**
   - Always match terms case-insensitively.
   - Map obvious synonyms to the allowed ontology terms/labels in the schema. 
     - E.g., `Lipofectamine 3000` or `Lipofectamine` should map to `lipofection`.
     - `lentivirus` or `lentiviral vectors` should map to `lentivirus transduction`.
     - `Jurkat cells` or `HEK293T` should map `model_system_label` to `cell_line`.

3. **Handling "Other":**
   - Return `"Other"` ONLY when the verbatim evidence explicitly describes a valid, clear concept for that field that is completely absent from the allowed terms in the schema.
   - Do NOT return `"Other"` if the evidence describes a concept that fits or maps to one of the allowed terms (even if via a synonym).
   - Do NOT propose new candidate terms, suggest corrections, or add rationale. Output exactly the label or `"Other"`.

4. **Null Handling:**
   - Return `null` if the evidence is `null`, missing, or highly ambiguous/insufficient to make a confident classification.
   - Do not make assumptions or extrapolate beyond the provided evidence.

5. **No Explanations or Candidate Proposals:**
   - Fill the schema fields directly.
   - Never output explanations, rationales, candidate lists, or other metadata in the fields.
