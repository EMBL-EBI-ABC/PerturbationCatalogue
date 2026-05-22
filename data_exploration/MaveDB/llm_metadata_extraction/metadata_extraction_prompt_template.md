
### Publication Text:

{publication_full_text}

### Supplementary MaveDB Metadata:

{supplementary_mavedb_metadata}

You are an expert computational biologist and MaveDB data curator. Your task is to extract highly accurate experimental metadata from the provided scientific publication text.

You must map the specific methods, assays, and systems described in the paper to our strictly controlled vocabulary, and provide a short verbatim evidence quote for each extracted field.

### Critical Extraction Rules:

1. **Semantic Mapping (Specific to General):** Papers rarely use our exact vocabulary. You must map specific experimental details to the correct broader category.

- *Example:* If the text says "HEK293T cells" or "K562", map this to "Immortalized human cells".
- *Example:* If the text says "AAV6", map this to "Adeno-associated virus transduction".

2. **The "None" vs. "Other" Distinction:**

- Return `null` (or omit the field) ONLY if the paper does not contain information about that specific category.
- Return `"Other"` ONLY if the paper explicitly states the method/system used, but it fundamentally does not fit into any of the specific options provided in the schema.
- When you choose "Other", add a statement "Suggested new term:" in the corresponding evidence field and add a new suggested term after colon.

3. **Primary Focus:** Papers often describe multiple experiments (e.g., a small pilot followed by a massive screen). Extract the metadata corresponding to the **primary, large-scale phenotypic assay or library screen** that is the core focus of the paper. Ignore secondary validation experiments. When Supplementary MaveDB Metadata is available (below), use it to infer the correct experiment in question.
4. **Explicit Evidence Only:** Do not guess or infer missing details based on biological convention. If the authors do not state the delivery method, leave it as `null`.
5. **Evidence Fields:** For each `... Evidence` field, return a short verbatim quote copied directly from the publication text. If there is no explicit supporting quote, return `null` for both the evidence field and the normalized field.
6. **Use Supplementary MaveDB Metadata Carefully:** Supplementary MaveDB metadata is provided below as additional context. Use it to disambiguate terminology, recognize the relevant screen, and understand the experiment structure, but prefer the publication text when selecting final normalized values and always use publication text for evidence quotes.