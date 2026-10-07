### Target Metadata Field
Name: {field_name}

### Existing Controlled Vocabulary (Allowed Terms)
{controlled_vocabulary}

### Unmapped Verbatim Evidence Strings (Classified as "Other") and Their Source Files
{evidence_list}

---

You are an Expert Ontology Curator. Your task is to analyze the provided list of verbatim evidence strings (which could not be mapped to any existing terms in the controlled vocabulary for `{field_name}` and were thus classified as `"Other"`) and identify distinct, recurring, and missing ontological concepts.

For this field, synthesize the evidence strings into concise, reusable ontology candidate terms. Do not omit an evidence source merely because it appears once: when a singleton has a defensible label, return it as a candidate too.

### Key Guidelines:
1. **Identify Distinct, Recurring Concepts:** 
   - Group synonym evidence strings together. For example, if multiple publications mention "Adeno-associated virus serotype 9 transduction" or "AAV9 gene delivery", they represent a single recurring concept: "adeno-associated virus transduction".
2. **Propose Reusable, General-Purpose Labels:**
   - Propose labels that are standard, reusable, and general-purpose ontology terms (following lowercase/naming conventions of the existing vocabulary).
   - Avoid highly specific protocol names, manufacturer details, or raw numeric parameters unless they form a distinct ontological class.
3. **Draft Quality Rationales and Evidence:**
   - For each proposed term, provide a clear rationale explaining why this concept should be added to the ontology.
   - List the direct quotes from the input evidence supporting this concept.
   - Provide the source file for each piece of supporting evidence.
4. **Cover Every Source:**
   - Include every supplied source file in the supporting evidence for at least one candidate whenever a defensible term can be proposed.
   - Do not invent a label when the evidence is insufficient; the application will explicitly route uncovered evidence to manual curator review.

Your output must strictly conform to the FieldCandidates schema.
