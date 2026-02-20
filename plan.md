# AI Explorer Integration Plan

## Current State

The AI Explorer uses **Gemini 2.5 Flash** with function calling, querying internal Elasticsearch + PostgreSQL and **Open Targets** via MCP. It renders tables, pie charts, and bar charts in a Data Portal panel via SSE streaming. There are 5 internal tools + 4 Open Targets MCP tools.

---

## Tier 1 -- High Impact, Integrate First

### 1. AlphaFold + 3D Protein Structure Viewer

**Why:** When a researcher asks about a perturbed gene, showing the predicted 3D structure is the single biggest "wow moment." Overlay MAVE variant effect data onto the structure to show which residues matter.

**API:** `GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}` -- returns structure URLs, pLDDT confidence, PAE data. No auth needed. An [AlphaFold MCP server](https://github.com/Augmented-Nature/AlphaFold-MCP-Server) also exists with tools for `get_structure`, `get_confidence_scores`, batch retrieval, and export.

**Viewer:** Embed **PDBe-Molstar** as a web component -- it's what AlphaFold DB itself uses, and it's maintained by EBI (your host). Just two CDN tags + `<pdbe-molstar molecule-id="AF-P51587-F1" alphafold-view="true">`. Add a `renderProteinStructure` case in `chat.js`, new `viz_type: "protein_structure"`.

**Example flow:** User asks "Show me BRCA2 structure" -> AI resolves gene to UniProt P51587 -> fetches AlphaFold metadata -> renders interactive 3D viewer with confidence coloring in the Data Portal.

**Status: DONE** (Phase 1 implemented)

---

### 2. UniProt -- Protein Function & Context

**Why:** Every perturbed gene encodes a protein. UniProt provides the canonical protein function description, domain architecture, disease associations, subcellular location, GO terms, and known variants/mutagenesis data. This is the "what does this gene do?" answer.

**API:** `GET https://rest.uniprot.org/uniprotkb/{accession}` + search/ID mapping endpoints. No auth, free. [Augmented-Nature UniProt MCP](https://github.com/Augmented-Nature/Augmented-Nature-UniProt-MCP-Server) has 26 tools, or you can add direct REST function declarations (simpler, fits existing pattern).

**Key tools to add:**
- `lookup_protein(gene_name)` -- function, domains, disease, subcellular location, GO terms
- `map_identifiers(ids, from_db, to_db)` -- bridge gene symbols to UniProt/Ensembl IDs
- `get_protein_variants(uniprot_id)` -- known variants + mutagenesis (especially relevant for MAVE interpretation)

**Status: DONE** (All three tools implemented: `lookup_protein`, `map_identifiers`, `get_protein_variants`)

---

### 3. Europe PMC -- Literature Search

**Why:** "What papers describe perturbation of gene X?" is a natural follow-up question. Europe PMC indexes all PubMed content PLUS has **text-mined entities** (genes, diseases, chemicals extracted from full text), making it far more useful than raw PubMed for cross-linking.

**API:** `GET https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=...` -- no auth, free, well-documented. Returns structured JSON with annotations.

**MCP servers available:** [BioMCP](https://github.com/genomoncology/biomcp) wraps Europe PMC + PubMed/PubTator3 + bioRxiv + clinical trials + variant annotations (COSMIC, ClinVar, CIViC) in one package. This is the most efficient single integration -- adding BioMCP as a second MCP server alongside Open Targets would immediately give you literature, variant annotation, and clinical trial search.

**Status: DONE** (Direct REST integration via `search_literature` tool querying Europe PMC API)

---

### 4. Pharos -- "Is This Gene Druggable?"

**Why:** After identifying a gene with strong perturbation effects, the translational question is "can we drug it?" Pharos classifies all ~20,000 human protein targets by druggability level (Tclin/Tchem/Tbio/Tdark) and shows existing drugs, ligands, and bioactivity data. Complements (doesn't duplicate) Open Targets.

**API:** GraphQL at `https://pharos.nih.gov/api` -- modern, flexible, well-documented.

**MCP server:** [pharos-mcp-server](https://github.com/QuentinCody/pharos-mcp-server) exists and is ready to integrate using your existing MCP client pattern.

**Status: DONE** (Direct GraphQL integration via `get_druggability` tool querying Pharos API)

---

## Tier 2 -- Molecular Consequence & Visualization

### 5. ProtVar -- Variant Molecular Consequences

**Why:** The key bridge between perturbation results and molecular impact. When the Catalogue shows a variant hit, ProtVar answers "What does this variant do to the protein?" with stability predictions (FoldX ddG), pathogenicity scores (EVE, ESM-1b, CADD), binding pocket/interface disruption, and PTM context. This data is **not available through Open Targets** -- OT links out to ProtVar but does not ingest its annotations. ProtVar is an EBI tool, making it a natural fit.

**API:** REST at `https://www.ebi.ac.uk/ProtVar/api/` ([Swagger docs](https://www.ebi.ac.uk/ProtVar/api/swagger-ui/index.html)). No auth, free. Accepts VCF, HGVS, dbSNP IDs, or UniProt accession + position. Pre-computes all 28M+ possible human missense variants for near-instant retrieval.

**Key data returned:**
- Protein stability change (FoldX ddG, >200M pre-computed values)
- Pathogenicity predictions (EVE, ESM-1b, CADD v1.6, AlphaMissense)
- Structural context (AlphaFold mapping, binding pockets via AutoSite, PPI interfaces)
- PTM disruption and active site proximity
- Population frequencies (gnomAD, TOPMed, ExAC)
- ClinVar/COSMIC classifications

**Key tools added:**
- `annotate_variant(variant)` -- full molecular consequence annotation for a missense variant (mapping + pathogenicity scores + FoldX stability)
- `get_variant_structural_context(uniprot_id, position)` -- structural impact at a specific residue (PDB structures, pockets, interactions, co-located variants)

**Limitation:** Human missense variants only (not indels, splicing, or structural variants).

**MCP server:** None exists. Direct REST integration (fits existing pattern).

**OT overlap note:** Open Targets includes basic VEP consequence terms and AlphaMissense scores, but ProtVar provides deeper protein-centric context (stability, pockets, interfaces, PTM) that OT does not have.

**Status: DONE** (Both tools implemented: `annotate_variant`, `get_variant_structural_context`)

---

The following visualizations render in the existing Data Portal using Plotly.js (already loaded) or lightweight CDN scripts:

### 6. Volcano Plot (Perturb-seq DEA)
The gold-standard differential expression visualization. The `perturb_seq_dea` table already has `log2_fc` and `padj` -- just add a `renderVolcanoPlot` function in `chat.js` using Plotly scatter traces. **Zero new dependencies.**

**Status: TODO**

### 7. Needle/Lollipop Plot (MAVE Variants)
Show variant effect scores mapped along protein sequence with domain annotations. The `mave_data` table has `position` and `score`. Render with Plotly scatter + shapes. **Zero new dependencies.**

**Status: TODO**

### 8. Gene Interaction Network (Cytoscape.js)
Interactive network: perturbed gene in center, differentially expressed genes radiating out with edge width = |log2FC|. Load `cytoscape.js` from CDN (~250KB). Transforms Perturb-seq results from a table into a visual network. **One CDN dependency.**

**Status: TODO**

### 9. Gene Summary Cards
Structured HTML cards with gene name, function (from UniProt), druggability (from Pharos), key stats from the Catalogue, and quick links to AlphaFold/UniProt/Open Targets. **Zero dependencies, just HTML/CSS.**

**Status: DONE** (Phase 1 implemented as `gene_card` viz_type)

---

## Tier 3 -- Domain-Specific Data Enrichment

### 10. DepMap -- Cancer CRISPR Dependencies
**Why:** The world's largest CRISPR screen dataset (1,865+ cancer cell lines). When the Catalogue shows a CRISPR hit for gene X, DepMap answers "Is X essential across all cancers or just this lineage?"

**API:** [Sanger DepMap REST API](https://api.cellmodelpassports.sanger.ac.uk/swagger) (JSONAPI v1.0, free). Broad's API is more download-oriented. A partial MCP tool exists via [BioAgent](https://zitniklab.hms.harvard.edu/bioagent/tools/remote/depmap_24q2.html).

**Integration:** Add as Gemini function declarations wrapping the Sanger REST API.

**Status: TODO**

### 11. MaveDB -- External MAVE Score Sets
**Why:** The Catalogue already has MAVE data, but MaveDB is the canonical source with 7M+ variant measurements. Cross-referencing lets the AI pull the latest scores and additional datasets.

**API:** Excellent FastAPI at `https://api.mavedb.org/docs`. No auth. Add as direct REST function declarations.

**Status: TODO**

### 12. Reactome -- Pathway Enrichment
**Why:** "What pathways are affected by this perturbation?" is a core analysis question. Submit a gene list from Perturb-seq results and get enriched pathways back.

**API:** REST at `https://reactome.org/AnalysisService/identifiers/` (POST gene list, get enrichment). No MCP server yet, but the REST API is straightforward.

**Status: TODO**

### 13. STRING -- Protein Interaction Networks
**Why:** Known physical and functional protein associations provide mechanistic context for perturbation effects.

**API:** REST at `https://string-db.org/api/`. MCP server available via [MCPMed](https://mcpmed.org/).

**Status: TODO**

### 14. Ensembl VEP -- Advanced Variant Consequence Prediction

**Why:** While Open Targets already ingests basic VEP consequence terms, direct VEP access unlocks the full plugin ecosystem: **SpliceAI** (splicing impact), **LOFTEE** (loss-of-function), **regulatory feature consequences** (enhancer/promoter/CTCF disruption), **motif feature consequences** (TF binding disruption), and **real-time annotation of novel variants** not yet in OT's release cycle. Handles all variant types (SNPs, indels, CNVs, structural) unlike ProtVar's missense-only scope.

**API:** REST at `https://rest.ensembl.org/vep/:species/hgvs/:hgvs_notation` (GET) or batch POST for up to 200 variants. No auth, free. [Full docs](https://rest.ensembl.org/documentation/info/vep_hgvs_get).

**Key tools to add:**
- `predict_variant_consequence(variant)` -- consequence terms, affected transcripts, amino acid changes, regulatory impact
- `batch_variant_consequences(variants)` -- batch annotation for variant lists from MAVE/CRISPR results

**MCP server:** An [unofficial Ensembl MCP server](https://github.com/Augmented-Nature/Ensembl-MCP-Server) exists with 25 tools including `get_variant_consequences`. No official Ensembl MCP yet.

**OT overlap note:** OT already has VEP consequence terms for variants in its index. Direct VEP is most valuable for: (1) regulatory/splicing consequences OT doesn't expose, (2) novel variants not yet in OT, (3) the SpliceAI/LOFTEE plugin annotations.

**Priority:** Lower than ProtVar -- add after ProtVar is integrated, since OT covers basic VEP needs.

**Status: TODO**

---

## Architecture Recommendation

The `ai_chat.py` already has a clean MCP client pattern (`init_open_targets_mcp` -> `_call_open_targets_tool`). Generalize this into a **multi-MCP dispatcher**:

```
Startup:
  1. init_mcp("Open Targets", "https://mcp.platform.opentargets.org/mcp")
  2. init_mcp("BioMCP", "<biomcp-url>")           # literature + variants
  3. init_mcp("Pharos", "<pharos-mcp-url>")        # druggability
  4. init_mcp("AlphaFold", "<alphafold-mcp-url>")  # protein structures

Runtime:
  tool_name -> look up which MCP server owns it -> dispatch
```

For resources without MCP servers (UniProt, MaveDB, DepMap, Reactome), add direct REST function declarations alongside the existing internal tools (`search_datasets`, `query_perturbation_data`, etc.).

---

## Suggested Phased Rollout

| Phase | Integrations | New Viz Types | Effort |
|-------|-------------|---------------|--------|
| **Phase 1** | UniProt (REST), AlphaFold (REST), Mol* viewer, Europe PMC, Pharos | `protein_structure`, `gene_card` | **DONE** |
| **Phase 2** | ProtVar (REST) **DONE**, volcano plot, needle plot | `volcano_plot`, `needle_plot` | 1-2 weeks |
| **Phase 3** | DepMap (REST), MaveDB (REST), Cytoscape networks | `network` | 2-3 weeks |
| **Phase 4** | Reactome, STRING, Ensembl VEP | `clustergram`, pathway diagrams | 2-3 weeks |

---

## Key Ecosystem Projects to Watch

- **[BioContextAI](https://www.nature.com/articles/s41587-025-02900-9)** -- EMBL-EBI backed community registry of MCP servers for biomedical data (Nature Biotech paper)
- **[Holy Bio MCP](https://github.com/longevity-genie/holy-bio-mcp)** -- 51+ bioinformatics functions via `gget-mcp`, `biothings-mcp`, etc.
- **[MCPmed](https://mcpmed.org/)** -- Curated hub of validated bioinformatics MCP servers (GEO, STRING, EMBL-EBI)
- **[BioinfoMCP](https://arxiv.org/html/2510.02139v1)** -- Auto-converts CLI bioinformatics tools to MCP servers (38 tools, 94.7% success rate)

The bioinformatics MCP ecosystem is rapidly maturing. The existing Open Targets MCP integration puts the project ahead of the curve -- the architecture is already in place to plug in additional servers as they become available.
