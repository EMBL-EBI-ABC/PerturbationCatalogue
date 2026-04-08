# Ensembl VEP Integration — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add two Ensembl VEP tools (`predict_variant_consequence`, `batch_variant_consequences`) to the AI Explorer for variant consequence prediction including SpliceAI, LOFTEE, CADD, AlphaMissense, regulatory features, and all variant types.

**Architecture:** Direct REST integration with `https://rest.ensembl.org/vep/human/` following the existing ProtVar/STRING handler pattern. Auto-detect input format (rsID, HGVS, VCF-style region) and route to the correct endpoint. Full plugin suite enabled. Results returned as dicts to Gemini, which presents them via the existing `create_visualization` table. No frontend changes needed.

**Tech Stack:** Python (async), aiohttp, Google Gemini function declarations, Ensembl VEP REST API

---

### Task 1: Add format detection helper and VEP constants

**Files:**
- Modify: `be/ai_chat.py` — insert after line 2065 (`PROTVAR_BASE = ...`)

**Step 1: Add VEP base URL constant and format detection function**

Insert after the `PROTVAR_BASE` line (line 2065):

```python
# --- Ensembl VEP ---

VEP_BASE = "https://rest.ensembl.org/vep/human"
VEP_PARAMS = {
    "CADD": "1",
    "SpliceAI": "1",
    "AlphaMissense": "1",
    "LoF": "1",
    "Conservation": "1",
    "domains": "1",
    "hgvs": "1",
    "canonical": "1",
    "uniprot": "1",
    "protein": "1",
}


def _detect_variant_format(variant: str) -> str:
    """Detect variant input format: 'rsid', 'hgvs', or 'region'.

    - rsID: starts with 'rs' followed by digits (e.g. rs56116432)
    - Region/VCF-style: chr:pos:ref:alt or chr:start-end:strand/allele
    - HGVS: everything else (transcript:c.X, genomic g.X, protein p.X)
    """
    v = variant.strip()
    if v.lower().startswith("rs") and v[2:].isdigit():
        return "rsid"
    # VCF-style region: e.g. "9:22125504:G:C" or "9:22125503-22125504:1/C"
    parts = v.split(":")
    if len(parts) >= 3 and parts[0].replace("chr", "").replace("X", "1").replace("Y", "1").replace("MT", "1").isdigit():
        return "region"
    return "hgvs"
```

**Step 2: Verify no syntax errors**

Run: `cd /Users/alexey/PerturbationCatalogue/be && python -c "import ai_chat; print('OK')"`
Expected: `OK`

**Step 3: Commit**

```
feat: add Ensembl VEP constants and format detection helper
```

---

### Task 2: Add `predict_variant_consequence` declaration

**Files:**
- Modify: `be/ai_chat.py` — insert declaration after `GET_VARIANT_STRUCTURAL_CONTEXT_DECLARATION` (after line 527), add to `INTERNAL_TOOL_DECLARATIONS` list

**Step 1: Add the FunctionDeclaration**

Insert after line 527 (end of `GET_VARIANT_STRUCTURAL_CONTEXT_DECLARATION`):

```python
PREDICT_VARIANT_CONSEQUENCE_DECLARATION = types.FunctionDeclaration(
    name="predict_variant_consequence",
    description=(
        "Predict the functional consequence of a genetic variant using Ensembl VEP (Variant Effect Predictor). "
        "Returns consequence type (missense, frameshift, splice, regulatory, etc.), impact severity, "
        "affected gene/transcript, protein change, and in-silico predictions: SIFT, PolyPhen, CADD, "
        "SpliceAI (splicing impact), AlphaMissense (missense pathogenicity), LOFTEE (loss-of-function), "
        "and conservation scores. Also reports regulatory feature consequences (enhancer/promoter disruption) "
        "and colocated known variants with population frequencies. "
        "Handles ALL variant types (SNPs, indels, frameshifts) — not limited to missense like ProtVar. "
        "Use for: 'What is the consequence of variant X?', 'Predict effect of rs123', "
        "'Is this a splice variant?', 'Does this variant affect regulatory regions?', "
        "'VEP annotation for variant X'. "
        "For protein stability (FoldX) and structural context (pockets, interfaces), use ProtVar tools instead."
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "variant": types.Schema(
                type="STRING",
                description=(
                    "Variant identifier in any format: "
                    "rsID (e.g. 'rs56116432'), "
                    "HGVS (e.g. 'ENST00000366667:c.803C>T', '9:g.22125504G>C'), "
                    "or VCF-style region (e.g. '9:22125504:G:C')"
                ),
            ),
        },
        required=["variant"],
    ),
)

BATCH_VARIANT_CONSEQUENCES_DECLARATION = types.FunctionDeclaration(
    name="batch_variant_consequences",
    description=(
        "Predict functional consequences for a batch of variants (up to 200) using Ensembl VEP. "
        "Returns a summary per variant: most severe consequence, impact, gene, protein change, "
        "and key scores (SIFT, PolyPhen, CADD, AlphaMissense). Efficient for annotating variant lists "
        "from MAVE, CRISPR, or other screening results. "
        "Use for: 'Annotate these variants', 'VEP for this list of variants', "
        "'What are the consequences of these SNPs?'"
    ),
    parameters=types.Schema(
        type="OBJECT",
        properties={
            "variants": types.Schema(
                type="ARRAY",
                items=types.Schema(
                    type="STRING",
                    description="Variant identifier (rsID, HGVS, or VCF-style region)",
                ),
                description=(
                    "List of variant identifiers (max 200). "
                    "All formats accepted: rsIDs, HGVS, VCF-style regions. Can be mixed."
                ),
            ),
        },
        required=["variants"],
    ),
)
```

**Step 2: Add both declarations to INTERNAL_TOOL_DECLARATIONS**

In the `INTERNAL_TOOL_DECLARATIONS` list (line ~653), add before the closing `]`:

```python
    PREDICT_VARIANT_CONSEQUENCE_DECLARATION,
    BATCH_VARIANT_CONSEQUENCES_DECLARATION,
```

**Step 3: Verify no syntax errors**

Run: `cd /Users/alexey/PerturbationCatalogue/be && python -c "import ai_chat; print(len(ai_chat.INTERNAL_TOOL_DECLARATIONS))"`
Expected: `21` (was 19, added 2)

**Step 4: Commit**

```
feat: add VEP tool declarations (predict_variant_consequence, batch_variant_consequences)
```

---

### Task 3: Implement `_tool_predict_variant_consequence` handler

**Files:**
- Modify: `be/ai_chat.py` — insert handler after VEP constants/helper (after Task 1 code), add to `TOOL_HANDLERS` dict

**Step 1: Add the handler function**

Insert after the `_detect_variant_format` function:

```python
def _extract_canonical_consequence(data: list) -> dict:
    """Extract the most relevant consequence from a VEP response.

    Prioritises the canonical transcript, then falls back to the transcript
    with the most severe impact.
    """
    if not data:
        return {}

    variant = data[0]
    result = {
        "input": variant.get("input", variant.get("id", "")),
        "allele_string": variant.get("allele_string", ""),
        "location": f"{variant.get('seq_region_name', '')}:{variant.get('start', '')}-{variant.get('end', '')}",
        "most_severe_consequence": variant.get("most_severe_consequence", ""),
        "strand": variant.get("strand"),
    }

    # Pick canonical transcript, or most severe
    tx_cons = variant.get("transcript_consequences", [])
    chosen = None
    for tc in tx_cons:
        if tc.get("canonical") == 1:
            chosen = tc
            break
    if not chosen and tx_cons:
        impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
        tx_cons_sorted = sorted(tx_cons, key=lambda t: impact_order.get(t.get("impact", "MODIFIER"), 3))
        chosen = tx_cons_sorted[0]

    if chosen:
        result["gene_symbol"] = chosen.get("gene_symbol", "")
        result["gene_id"] = chosen.get("gene_id", "")
        result["transcript_id"] = chosen.get("transcript_id", "")
        result["biotype"] = chosen.get("biotype", "")
        result["impact"] = chosen.get("impact", "")
        result["consequence_terms"] = chosen.get("consequence_terms", [])
        result["amino_acids"] = chosen.get("amino_acids")
        result["codons"] = chosen.get("codons")
        result["protein_start"] = chosen.get("protein_start")
        result["hgvsc"] = chosen.get("hgvsc")
        result["hgvsp"] = chosen.get("hgvsp")

        # In-silico predictions
        predictions = {}
        if chosen.get("sift_prediction"):
            predictions["sift"] = f"{chosen['sift_prediction']}({chosen.get('sift_score', '')})"
        if chosen.get("polyphen_prediction"):
            predictions["polyphen"] = f"{chosen['polyphen_prediction']}({chosen.get('polyphen_score', '')})"
        if chosen.get("cadd_phred") is not None:
            predictions["cadd_phred"] = chosen["cadd_phred"]
        if chosen.get("cadd_raw") is not None:
            predictions["cadd_raw"] = chosen["cadd_raw"]
        # AlphaMissense
        if chosen.get("am_class"):
            predictions["alphamissense"] = f"{chosen['am_class']}({chosen.get('am_pathogenicity', '')})"
        # SpliceAI
        splice_keys = [k for k in chosen if k.startswith("spliceai_pred_ds_")]
        if splice_keys:
            splice_scores = {k.replace("spliceai_pred_", ""): chosen[k] for k in sorted(chosen.keys()) if k.startswith("spliceai_pred_")}
            predictions["spliceai"] = splice_scores
        # LOFTEE
        if chosen.get("lof"):
            predictions["loftee"] = chosen["lof"]
            if chosen.get("lof_filter"):
                predictions["loftee_filter"] = chosen["lof_filter"]
        # Conservation
        if chosen.get("conservation") is not None:
            predictions["conservation"] = chosen["conservation"]

        if predictions:
            result["predictions"] = predictions

        # Protein domains
        domains = chosen.get("domains")
        if domains:
            result["protein_domains"] = [
                f"{d.get('db', '')}:{d.get('name', '')}" for d in domains[:5]
            ]

        # UniProt
        if chosen.get("swissprot"):
            result["uniprot_id"] = chosen["swissprot"][0] if isinstance(chosen["swissprot"], list) else chosen["swissprot"]

    # Regulatory consequences
    reg_cons = variant.get("regulatory_feature_consequences", [])
    if reg_cons:
        result["regulatory_consequences"] = [
            {
                "regulatory_feature_id": rc.get("regulatory_feature_id", ""),
                "biotype": rc.get("biotype", ""),
                "consequence_terms": rc.get("consequence_terms", []),
                "impact": rc.get("impact", ""),
            }
            for rc in reg_cons[:5]
        ]

    # Motif consequences
    motif_cons = variant.get("motif_feature_consequences", [])
    if motif_cons:
        result["motif_consequences"] = [
            {
                "motif_feature_id": mc.get("motif_feature_id", ""),
                "consequence_terms": mc.get("consequence_terms", []),
                "impact": mc.get("impact", ""),
            }
            for mc in motif_cons[:5]
        ]

    # Colocated variants (population frequencies)
    colocated = variant.get("colocated_variants", [])
    if colocated:
        coloc_summary = []
        for cv in colocated[:3]:
            entry = {"id": cv.get("id", "")}
            freqs = cv.get("frequencies", {})
            if freqs:
                # Get global frequency from first allele
                for allele, allele_freqs in freqs.items():
                    if "gnomade" in allele_freqs:
                        entry["gnomad_global"] = allele_freqs["gnomade"]
                    break
            coloc_summary.append(entry)
        result["colocated_variants"] = coloc_summary

    return result


async def _tool_predict_variant_consequence(args: dict) -> dict:
    """Predict variant consequence using Ensembl VEP REST API."""
    variant = args.get("variant", "").strip()
    if not variant:
        return {"error": "variant is required"}

    fmt = _detect_variant_format(variant)

    # Build URL based on format
    if fmt == "rsid":
        url = f"{VEP_BASE}/id/{variant}"
    elif fmt == "region":
        # Convert "9:22125504:G:C" to "9:22125504-22125504:1/C" for the region endpoint
        parts = variant.split(":")
        if len(parts) == 4:
            chrom, pos, ref, alt = parts
            url = f"{VEP_BASE}/region/{chrom}:{pos}-{int(pos) + len(ref) - 1}:1/{alt}"
        else:
            # Already in region format or close enough
            url = f"{VEP_BASE}/region/{variant}"
    else:
        url = f"{VEP_BASE}/hgvs/{variant}"

    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = dict(VEP_PARAMS)

    try:
        async with aiohttp.ClientSession() as http:
            async with http.get(
                url, params=params, headers=headers, timeout=aiohttp.ClientTimeout(total=30)
            ) as resp:
                if resp.status == 400:
                    body = await resp.text()
                    return {"error": f"Ensembl VEP error (invalid input): {body[:300]}"}
                if resp.status == 429:
                    return {"error": "Ensembl VEP rate limit reached. Please try again in a moment."}
                if resp.status != 200:
                    body = await resp.text()
                    return {"error": f"Ensembl VEP returned status {resp.status}: {body[:300]}"}
                data = await resp.json()
    except Exception as exc:
        return {"error": f"Ensembl VEP API error: {str(exc)}"}

    result = _extract_canonical_consequence(data)
    result["ensembl_vep_url"] = f"https://www.ensembl.org/Homo_sapiens/Tools/VEP"
    return result
```

**Step 2: Add to TOOL_HANDLERS**

In the `TOOL_HANDLERS` dict (~line 2504), add:

```python
    "predict_variant_consequence": _tool_predict_variant_consequence,
```

**Step 3: Verify no syntax errors**

Run: `cd /Users/alexey/PerturbationCatalogue/be && python -c "import ai_chat; print('predict_variant_consequence' in ai_chat.TOOL_HANDLERS)"`
Expected: `True`

**Step 4: Commit**

```
feat: implement predict_variant_consequence handler (Ensembl VEP)
```

---

### Task 4: Implement `_tool_batch_variant_consequences` handler

**Files:**
- Modify: `be/ai_chat.py` — insert handler after `_tool_predict_variant_consequence`, add to `TOOL_HANDLERS`

**Step 1: Add the batch handler function**

Insert after `_tool_predict_variant_consequence`:

```python
async def _tool_batch_variant_consequences(args: dict) -> dict:
    """Batch-predict variant consequences using Ensembl VEP POST endpoints."""
    variants = args.get("variants", [])
    if not variants:
        return {"error": "variants list is required"}
    if len(variants) > 200:
        return {"error": f"Maximum 200 variants per batch, got {len(variants)}"}

    # Group variants by format for batch POST
    groups: dict[str, list[str]] = {"rsid": [], "hgvs": [], "region": []}
    for v in variants:
        fmt = _detect_variant_format(str(v).strip())
        groups[fmt].append(str(v).strip())

    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = dict(VEP_PARAMS)

    all_results = []

    async with aiohttp.ClientSession() as http:
        # POST rsIDs
        if groups["rsid"]:
            url = f"{VEP_BASE}/id"
            body = {"ids": groups["rsid"]}
            try:
                async with http.post(
                    url, json=body, params=params, headers=headers,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        all_results.extend(data)
                    else:
                        body_text = await resp.text()
                        all_results.append({"error": f"VEP rsID batch error ({resp.status}): {body_text[:200]}"})
            except Exception as exc:
                all_results.append({"error": f"VEP rsID batch error: {str(exc)}"})

        # POST HGVS
        if groups["hgvs"]:
            url = f"{VEP_BASE}/hgvs"
            body = {"hgvs_notations": groups["hgvs"]}
            try:
                async with http.post(
                    url, json=body, params=params, headers=headers,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        all_results.extend(data)
                    else:
                        body_text = await resp.text()
                        all_results.append({"error": f"VEP HGVS batch error ({resp.status}): {body_text[:200]}"})
            except Exception as exc:
                all_results.append({"error": f"VEP HGVS batch error: {str(exc)}"})

        # POST regions — need to convert "9:22125504:G:C" to VCF-like format "9 22125504 . G C . . ."
        if groups["region"]:
            vcf_lines = []
            for v in groups["region"]:
                parts = v.split(":")
                if len(parts) == 4:
                    chrom, pos, ref, alt = parts
                    vcf_lines.append(f"{chrom} {pos} . {ref} {alt} . . .")
                else:
                    vcf_lines.append(v)
            url = f"{VEP_BASE}/region"
            body = {"variants": vcf_lines}
            try:
                async with http.post(
                    url, json=body, params=params, headers=headers,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        all_results.extend(data)
                    else:
                        body_text = await resp.text()
                        all_results.append({"error": f"VEP region batch error ({resp.status}): {body_text[:200]}"})
            except Exception as exc:
                all_results.append({"error": f"VEP region batch error: {str(exc)}"})

    # Build summary table
    summaries = []
    errors = []
    for item in all_results:
        if "error" in item:
            errors.append(item["error"])
            continue
        row = _extract_canonical_consequence([item])
        summaries.append(row)

    output = {
        "variant_count": len(summaries),
        "annotations": summaries,
    }
    if errors:
        output["errors"] = errors

    return output
```

**Step 2: Add to TOOL_HANDLERS**

```python
    "batch_variant_consequences": _tool_batch_variant_consequences,
```

**Step 3: Verify no syntax errors**

Run: `cd /Users/alexey/PerturbationCatalogue/be && python -c "import ai_chat; print('batch_variant_consequences' in ai_chat.TOOL_HANDLERS)"`
Expected: `True`

**Step 4: Commit**

```
feat: implement batch_variant_consequences handler (Ensembl VEP)
```

---

### Task 5: Update system prompt with VEP routing guidance

**Files:**
- Modify: `be/ai_chat.py` — update `SYSTEM_INSTRUCTION` string

**Step 1: Add VEP routing to CRITICAL TOOL ROUTING section**

In the CRITICAL TOOL ROUTING block (after line 2680, the ProtVar routing lines), add:

```
- "What is the consequence of variant X?" / "VEP for X" / "Predict effect of rs123" / "Is this a splice variant?" / "Regulatory impact of variant" → call predict_variant_consequence (Ensembl VEP)
- "Annotate these variants" / "VEP for this list" / "Consequences of these SNPs" → call batch_variant_consequences (Ensembl VEP)
```

**Step 2: Add ENSEMBL VEP section to system prompt**

After the existing ProtVar workflow section (after line 2791), add:

```
ENSEMBL VEP (VARIANT CONSEQUENCE PREDICTION):
You have access to Ensembl VEP for comprehensive variant consequence prediction.

- predict_variant_consequence: Predict the functional consequence of ANY variant type (SNPs, indels, frameshifts, splice variants, regulatory). Returns consequence terms, impact severity (HIGH/MODERATE/LOW/MODIFIER), affected gene/transcript, protein change, and in-silico predictions: SIFT, PolyPhen, CADD, SpliceAI (splicing), AlphaMissense (missense pathogenicity), LOFTEE (loss-of-function), conservation. Also reports regulatory feature consequences (enhancer/promoter/CTCF disruption) and colocated known variants with gnomAD frequencies. Accepts rsIDs (rs56116432), HGVS (ENST00000366667:c.803C>T), or VCF-style (9:22125504:G:C).

- batch_variant_consequences: Annotate up to 200 variants in one call. Efficient for processing variant lists from screening results. Returns a summary per variant. Present results as a table using create_visualization.

VEP vs ProtVar guidance:
- Use VEP (predict_variant_consequence) for: consequence type prediction, splicing impact (SpliceAI), loss-of-function (LOFTEE), regulatory/motif consequences, any variant type (indels, frameshifts, not just missense), novel variants
- Use ProtVar (annotate_variant) for: protein stability (FoldX ddG), structural context (binding pockets, PPI interfaces), EVE/ESM-1b scores, PTM disruption
- For missense variants, both tools are complementary — VEP gives consequence + CADD + SpliceAI + regulatory context, ProtVar gives stability + structural context
```

**Step 3: Verify no syntax errors**

Run: `cd /Users/alexey/PerturbationCatalogue/be && python -c "import ai_chat; print('predict_variant_consequence' in ai_chat.SYSTEM_INSTRUCTION)"`
Expected: `True`

**Step 4: Commit**

```
feat: add Ensembl VEP routing guidance to system prompt
```

---

### Task 6: Update plan.md status

**Files:**
- Modify: `plan.md` — update task 11 status and Phase 4 row

**Step 1: Update status**

Change task 11 status from `**Status: TODO**` to `**Status: DONE**` and update the Phase 4 row in the phased rollout table.

**Step 2: Commit**

```
docs: mark Ensembl VEP integration (task 11) as done
```
