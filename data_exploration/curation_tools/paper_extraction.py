"""
paper_extraction.py — LLM-based extraction of genetic perturbation experiment metadata
from scientific paper text.

Given a paper .txt file and the ObsSchema field definitions from
perturbseq_anndata_schema.py, this module uses the Google Gemini API with structured
JSON output to produce per-experiment metadata records compatible with the existing
curation pipeline.

Usage:
    from curation_tools.paper_extraction import extract_paper, batch_extract, ExtractionConfig
    from pathlib import Path

    # Single paper
    records = extract_paper(Path("pdf_text/27869803.txt"))

    # Batch
    batch_extract(
        txt_dir=Path("pdf_text"),
        output_dir=Path("curated_gemini_v2"),
    )
"""

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Fields extractable from paper text
# ---------------------------------------------------------------------------
# Excluded: per-cell fields (guide_sequence, perturbed_target_*, sample_id,
# perturbation_name, significant, significance_criteria, technical/biological_replicate,
# score_interpretation, perturbed_target_chromosome*), all *_id ontology fields
# (require ontology lookup), license_*, associated_datasets, dataset_id,
# perturbed_target_number (per-cell count).
EXTRACTABLE_FIELDS: list[str] = [
    # dataset
    "data_modality",
    # perturbation
    "perturbation_type_label",
    "timepoint",
    "treatment_label",
    # model system
    "model_system_label",
    "species",
    "tissue_label",
    "cell_type_label",
    "cell_line_label",
    "sex_label",
    "developmental_stage_label",
    "disease_label",
    # study
    "study_title",
    "study_uri",
    "study_year",
    "first_author",
    "last_author",
    # experiment
    "experiment_title",
    "experiment_summary",
    "number_of_perturbed_targets",
    "number_of_perturbed_samples",
    # library generation
    "library_generation_type_label",
    "library_generation_method_label",
    "enzyme_delivery_method_label",
    "library_delivery_method_label",
    "enzyme_integration_state_label",
    "library_integration_state_label",
    "enzyme_expression_control_label",
    "library_expression_control_label",
    # library specs
    "library_name",
    "library_uri",
    "library_format_label",
    "library_scope_label",
    "library_perturbation_type_label",
    "library_manufacturer",
    "library_lentiviral_generation",
    "library_grnas_per_target",
    "library_total_grnas",
    "library_total_variants",
    # assay
    "readout_dimensionality_label",
    "readout_type_label",
    "readout_technology_label",
    "method_name_label",
    "method_uri",
    "sequencing_library_kit_label",
    "sequencing_platform_label",
    "sequencing_strategy_label",
    "software_counts_label",
    "software_analysis_label",
    "reference_genome_label",
]

# Fields that require specific controlled vocabulary values (from isin checks in ObsSchema)
CONSTRAINED_FIELDS: dict[str, list[str]] = {
    "data_modality": ["Perturb-seq", "CRISPR screen", "MAVE"],
    "perturbation_type_label": ["CRISPRn", "CRISPRi", "CRISPRa", "DMS"],
    "model_system_label": ["cell_line", "primary_cell", "organoid", "yeast"],
    "sex_label": ["female", "male", "mixed", "unknown"],
    "developmental_stage_label": [
        "embryonic",
        "fetal",
        "neonatal",
        "child",
        "adolescent",
        "adult",
        "senior adult",
    ],
    "enzyme_delivery_method_label": [
        "lipofection",
        "nucleofection",
        "retrovirus transduction",
        "lentivirus transduction",
        "transformation",
        "nanoparticle-mediated transfection",
    ],
    "library_delivery_method_label": [
        "lipofection",
        "nucleofection",
        "retrovirus transduction",
        "lentivirus transduction",
        "transformation",
        "nanoparticle-mediated transfection",
    ],
    "enzyme_integration_state_label": [
        "random locus integration",
        "targeted locus integration",
        "native locus replacement",
        "non-integrative transgene expression",
    ],
    "library_integration_state_label": [
        "random locus integration",
        "targeted locus integration",
        "native locus replacement",
        "non-integrative transgene expression",
    ],
    "enzyme_expression_control_label": [
        "constitutive transgene expression",
        "inducible transgene expression",
        "native promoter-driven transgene expression",
        "degradation domain-based transgene control",
    ],
    "library_expression_control_label": [
        "constitutive transgene expression",
        "inducible transgene expression",
        "native promoter-driven transgene expression",
        "degradation domain-based transgene control",
    ],
    "library_format_label": ["pooled", "arrayed", "arrayed|pooled", "in vivo"],
    "library_scope_label": ["focused", "genome-wide"],
    "library_perturbation_type_label": [
        "knockout",
        "inhibition",
        "activation",
        "base editing",
        "prime editing",
        "mutagenesis",
    ],
    "readout_dimensionality_label": [
        "single-dimensional assay",
        "high-dimensional assay",
    ],
    "readout_type_label": ["transcriptomic", "proteomic", "phenotypic"],
    "readout_technology_label": [
        "single-cell rna-seq",
        "population growth assay",
        "flow cytometry",
    ],
    "method_name_label": [
        "Perturb-seq",
        "Perturb-CITE-seq",
        "scRNA-seq",
        "proliferation CRISPR screen",
        "DMS-TileSeq",
        "DMS-BarSeq",
        "Joined and refined DMS-BarSeq and DMS-TileSeq",
        "Combined DMS-BarSeq and DMS-TileSeq",
    ],
    "sequencing_library_kit_label": [
        "10x Genomics Chromium GEM-X Single Cell 5-prime kit v3",
        "10x Genomics Chromium Next GEM Single Cell 5-prime HT Kit v2",
        "10x Genomics Single Cell 3-prime",
        "10x Genomics Single Cell 3-prime v2",
        "10x Genomics Single Cell 3-prime v3",
        "Nextera XT DNA Library Preparation Kit",
        "GEM-X Flex Gene Expression Human n-plex kit",
    ],
    "sequencing_platform_label": [
        "Illumina NovaSeq X",
        "Illumina NovaSeq X Plus",
        "Illumina HiSeq 4000",
        "Illumina HiSeq 2500",
        "Illumina HiSeq 2000",
        "Illumina NovaSeq 6000",
        "Illumina NextSeq 500",
        "Ultima Genomics UG100",
    ],
    "sequencing_strategy_label": [
        "barcode sequencing",
        "direct sequencing",
        "barcode sequencing|direct sequencing",
    ],
    "software_counts_label": ["custom", "MaGeCK", "CellRanger", "Drop-seq Tools"],
    "software_analysis_label": [
        "custom",
        "MAGeCK",
        "Achilles",
        "TRADE",
        "Seurat",
        "MAST",
        "scanpy",
    ],
    "reference_genome_label": ["GRCh38", "GRCh37"],
}

# Fields whose values should be coerced to integers (or left None)
_INTEGER_FIELDS: set[str] = {
    "study_year",
    "library_total_grnas",
    "library_total_variants",
    "number_of_perturbed_targets",
    "number_of_perturbed_samples",
}

# Human-readable descriptions for each extractable field (sourced from ObsSchema)
_FIELD_DESCRIPTIONS: dict[str, str] = {
    "data_modality": "Data modality of the experiment. Perturb-seq = single-cell RNA-seq readout; CRISPR screen = population-level phenotypic readout; MAVE = deep mutational scanning.",
    "perturbation_type_label": "Perturbation type of the investigated sample.",
    "timepoint": "Duration of the experiment in ISO 8601 format. Example: P1DT12H30M15S (1 day, 12 hours, 30 minutes, 15 seconds). P16DT0H0M0S = 16 days.",
    "treatment_label": "Treatment/compound used to stimulate cells. Use ChEMBL compound label for chemical entities (e.g. 'doxycycline', 'ricin'). Use null if no treatment information is mentioned in the text. If there are treatment conditions, and there is an untreated condition, use 'untreated_control'.",
    "model_system_label": "Type of biological model system.",
    "species": "Species name. Currently only 'Homo sapiens' is supported.",
    "tissue_label": "Tissue of origin (UBERON ontology label, e.g. 'blood', 'lung').",
    "cell_type_label": "Cell type label (Cell Ontology, CL). E.g. 'erythroleukemia cell line cell'.",
    "cell_line_label": "Cell line name (Cell Line Ontology / Cellosaurus). E.g. 'K562', 'HEK293T'.",
    "sex_label": "Biological sex of the cell source.",
    "developmental_stage_label": "Developmental stage of the cell source.",
    "disease_label": "Disease associated with the model system (MONDO ontology). E.g. 'chronic myelogenous leukemia'.",
    "study_title": "Full title of the publication.",
    "study_uri": "DOI or URI of the publication. E.g. '10.1016/j.cell.2014.09.029'.",
    "study_year": "Year the paper was published. Integer.",
    "first_author": "Full name(s) of the first author(s).",
    "last_author": "Full name(s) of the last/corresponding author(s).",
    "experiment_title": "A concise, descriptive title for this specific experiment (not the paper title).",
    "experiment_summary": "1-3 sentence summary of what was done and why in this experiment.",
    "number_of_perturbed_targets": "Total number of distinct gene/genomic loci targeted. Integer.",
    "number_of_perturbed_samples": "Total number of perturbed cells/samples in the experiment. Integer or null.",
    "library_generation_type_label": "Library generation type (EFO). E.g. 'endogenous perturbation method'.",
    "library_generation_method_label": "Specific CRISPR system used. E.g. 'SpCas9', 'dCas9-KRAB', 'dCas9-SunTag', 'ABE8e', 'PE2'.",
    "enzyme_delivery_method_label": "How the Cas9/effector protein was delivered into cells.",
    "library_delivery_method_label": "How the guide RNA library was delivered into cells.",
    "enzyme_integration_state_label": "How the Cas9/effector integrates (or not) into the genome.",
    "library_integration_state_label": "How the guide RNA library integrates (or not) into the genome.",
    "enzyme_expression_control_label": "How Cas9/effector expression is controlled.",
    "library_expression_control_label": "How guide RNA expression is controlled.",
    "library_name": "Name of the guide RNA library used. E.g. 'Genome-scale CRISPRi sgRNA library v2'.",
    "library_uri": "Accession or URI for the guide RNA library, if publicly deposited. E.g. Addgene ID.",
    "library_format_label": "Physical format of the library screen.",
    "library_scope_label": "Coverage scope of the library.",
    "library_perturbation_type_label": "Effect type of the perturbation in the library.",
    "library_manufacturer": "Lab or vendor that produced the library. E.g. 'Addgene', 'Broad Institute', 'Weissman lab'.",
    "library_lentiviral_generation": "Lentiviral generation number if applicable (e.g. '2', '3').",
    "library_grnas_per_target": "Number of guide RNAs per targeted gene/locus. E.g. '5', '4-6'.",
    "library_total_grnas": "Total number of guide RNAs in the library. Integer.",
    "library_total_variants": "For MAVE studies: total number of sequence variants in the library. Integer.",
    "readout_dimensionality_label": "Whether the readout captures one measurement per cell ('single-dimensional assay') or many ('high-dimensional assay').",
    "readout_type_label": "Category of biological readout.",
    "readout_technology_label": "Specific technology used for the readout.",
    "method_name_label": "Established method name for the experimental approach.",
    "method_uri": "URI or DOI associated with the method protocol, if available.",
    "sequencing_library_kit_label": "Sequencing library preparation kit. E.g. '10x Genomics Single Cell 3-prime v3'.",
    "sequencing_platform_label": "Sequencing instrument platform.",
    "sequencing_strategy_label": "Sequencing strategy used to read out guide RNA identity.",
    "software_counts_label": "Software used to count reads/UMIs.",
    "software_analysis_label": "Software used for statistical analysis of screen hits.",
    "reference_genome_label": "Human reference genome assembly version.",
}


# ---------------------------------------------------------------------------
# FieldMeta dataclass
# ---------------------------------------------------------------------------


@dataclass
class FieldMeta:
    name: str
    description: str
    allowed_values: list[str] | None  # None means free text
    nullable: bool
    python_type: str  # "str", "int", "float"


# ---------------------------------------------------------------------------
# SchemaIntrospector
# ---------------------------------------------------------------------------


class SchemaIntrospector:
    """
    Builds FieldMeta objects for each extractable field and generates
    prompt fragments and a Gemini-compatible response schema.
    """

    @classmethod
    def get_field_meta(
        cls, field_names: list[str] = EXTRACTABLE_FIELDS
    ) -> dict[str, FieldMeta]:
        """Return a dict of field_name -> FieldMeta for each extractable field."""
        result: dict[str, FieldMeta] = {}
        for name in field_names:
            description = _FIELD_DESCRIPTIONS.get(name, f"See ObsSchema for '{name}'.")
            allowed = CONSTRAINED_FIELDS.get(name)
            nullable = name not in {
                "data_modality",
                "perturbation_type_label",
                "model_system_label",
                "study_title",
                "study_uri",
                "study_year",
                "experiment_title",
                "number_of_perturbed_targets",
                "species",
            }
            python_type = "int" if name in _INTEGER_FIELDS else "str"
            result[name] = FieldMeta(
                name=name,
                description=description,
                allowed_values=allowed,
                nullable=nullable,
                python_type=python_type,
            )
        return result

    @classmethod
    def build_field_reference_block(cls, field_meta: dict[str, FieldMeta]) -> str:
        """Return a formatted multi-line string describing each field for the prompt."""
        lines: list[str] = ["FIELD REFERENCE:", "=" * 60, ""]
        for meta in field_meta.values():
            req = "OPTIONAL" if meta.nullable else "REQUIRED"
            lines.append(f"{meta.name} [{req}]")
            lines.append(f"  Description: {meta.description}")
            if meta.allowed_values:
                lines.append(
                    f"  Allowed values (use exactly as written): "
                    + " | ".join(meta.allowed_values)
                )
            if meta.python_type == "int":
                lines.append("  Type: integer or null")
            lines.append("")
        return "\n".join(lines)

    @classmethod
    def build_response_schema(cls) -> dict[str, Any]:
        """
        Return a JSON Schema dict compatible with Gemini's response_schema parameter.
        Shape: array of objects, all fields nullable, constrained fields have enum.
        """
        properties: dict[str, Any] = {}
        for name in EXTRACTABLE_FIELDS:
            if name in _INTEGER_FIELDS:
                prop: dict[str, Any] = {"type": "integer", "nullable": True}
            else:
                prop = {"type": "string", "nullable": True}
            if name in CONSTRAINED_FIELDS:
                prop["enum"] = CONSTRAINED_FIELDS[name]
            properties[name] = prop

        return {
            "type": "array",
            "items": {
                "type": "object",
                "properties": properties,
                "required": [],  # all nullable in schema; prompt handles null instruction
            },
        }


# ---------------------------------------------------------------------------
# FewShotProvider
# ---------------------------------------------------------------------------


@dataclass
class _FewShotExample:
    paper_id: str
    paper_text_snippet: str
    extracted_records: list[dict[str, Any]]


class FewShotProvider:
    """
    Loads few-shot examples from curated_gemini/ paired with paper text from pdf_text/.

    Selects up to n_examples diverse examples prioritizing coverage of different
    method_name_label values (CRISPR screen, Perturb-seq, etc.).
    """

    MAX_SNIPPET_CHARS: int = 4_000

    def __init__(
        self,
        curated_gemini_dir: Path,
        pdf_text_dir: Path,
        n_examples: int = 3,
        seed: int = 42,
    ) -> None:
        self.curated_gemini_dir = Path(curated_gemini_dir)
        self.pdf_text_dir = Path(pdf_text_dir)
        self.n_examples = n_examples
        self.seed = seed

    def _find_paired_examples(self) -> list[tuple[Path, Path]]:
        """Find (json_path, txt_path) pairs where both files exist (PMID stem match)."""
        pairs: list[tuple[Path, Path]] = []
        for json_path in sorted(self.curated_gemini_dir.glob("*.json")):
            txt_path = self.pdf_text_dir / f"{json_path.stem}.txt"
            if txt_path.exists():
                pairs.append((json_path, txt_path))
        logger.debug(
            "Found %d paired few-shot examples in %s",
            len(pairs),
            self.curated_gemini_dir,
        )
        return pairs

    def _method_label(self, records: list[dict[str, Any]]) -> str:
        """Return the method_name_label of the first non-null record."""
        for rec in records:
            val = rec.get("method_name_label")
            if val:
                return val
        return "unknown"

    def get_examples(self) -> list[_FewShotExample]:
        """
        Select up to n_examples prioritizing diversity across method_name_label.
        Uses seed for reproducibility.
        """
        import random

        rng = random.Random(self.seed)
        pairs = self._find_paired_examples()
        if not pairs:
            logger.warning(
                "No paired few-shot examples found. Few-shot block will be empty."
            )
            return []

        # Load records for each pair and group by method label
        loaded: list[tuple[str, str, list[dict]]] = []  # (paper_id, method_label, records)
        for json_path, txt_path in pairs:
            try:
                records = json.loads(json_path.read_text(encoding="utf-8"))
                if not isinstance(records, list) or not records:
                    continue
                method = self._method_label(records)
                snippet = txt_path.read_text(encoding="utf-8", errors="replace")[
                    : self.MAX_SNIPPET_CHARS
                ]
                loaded.append((json_path.stem, method, records, snippet))  # type: ignore[arg-type]
            except Exception as exc:
                logger.debug("Skipping few-shot pair %s: %s", json_path.stem, exc)

        # Diversity selection: pick one per unique method_label, up to n_examples
        seen_methods: set[str] = set()
        selected: list[_FewShotExample] = []
        rng.shuffle(loaded)  # type: ignore[arg-type]
        for paper_id, method, records, snippet in loaded:  # type: ignore[misc]
            if method not in seen_methods:
                seen_methods.add(method)
                selected.append(
                    _FewShotExample(
                        paper_id=paper_id,
                        paper_text_snippet=snippet,
                        extracted_records=records,
                    )
                )
            if len(selected) >= self.n_examples:
                break

        # If still short, fill with remaining (any method)
        if len(selected) < self.n_examples:
            selected_ids = {ex.paper_id for ex in selected}
            for paper_id, method, records, snippet in loaded:  # type: ignore[misc]
                if paper_id not in selected_ids:
                    selected.append(
                        _FewShotExample(
                            paper_id=paper_id,
                            paper_text_snippet=snippet,
                            extracted_records=records,
                        )
                    )
                    if len(selected) >= self.n_examples:
                        break

        logger.debug("Selected %d few-shot examples.", len(selected))
        return selected

    def format_for_prompt(self, examples: list[_FewShotExample]) -> str:
        """Format examples as fenced blocks for inclusion in the user prompt."""
        if not examples:
            return ""
        blocks: list[str] = []
        for i, ex in enumerate(examples, start=1):
            records_json = json.dumps(ex.extracted_records, indent=2, ensure_ascii=False)
            blocks.append(
                f"--- FEW-SHOT EXAMPLE {i} (paper_id: {ex.paper_id}) ---\n"
                f"PAPER TEXT (first {self.MAX_SNIPPET_CHARS} characters):\n"
                f"{ex.paper_text_snippet}\n\n"
                f"EXPECTED OUTPUT (JSON array):\n"
                f"{records_json}\n"
                f"--- END EXAMPLE {i} ---"
            )
        return "\n\n".join(blocks)


# ---------------------------------------------------------------------------
# PromptBuilder
# ---------------------------------------------------------------------------


class PromptBuilder:
    """
    Assembles system and user prompts for Gemini extraction.
    """

    EXPERIMENT_DEFINITION: str = """
WHAT COUNTS AS ONE EXPERIMENT RECORD:
======================================
One record = one coherent perturbation experiment run under a fixed set of conditions.

CREATE A SEPARATE RECORD when ANY of the following differ:
- perturbation_type_label (e.g., CRISPRi vs. CRISPRa — always separate records)
- cell_line_label or cell_type_label (different cell type or cell line)
- treatment_label (different chemical treatment or stimulus)
- library_scope_label when the paper explicitly describes independent screens
  (e.g., a genome-wide screen AND a focused validation screen are separate)
- readout_type_label (e.g., transcriptomic vs. phenotypic readout in the same paper)

DO NOT create separate records for:
- Technical or biological replicates of the same experiment
- Multiple time points of the same screen (use the predominant timepoint)
- Small validation experiments (individual gRNA re-tests, clonogenic assays,
  Western blots, survival assays validating a few hits)
- Computational re-analyses of previously published data

EXAMPLES:
- CRISPRi screen + CRISPRa screen in K562 → 2 records
- Same screen repeated in 3 cell lines → 3 records
- Genome-wide screen + focused sublibrary screen (same cell line, same perturbation type) → 2 records
- Perturb-seq (scRNA-seq readout) + parallel growth screen (phenotypic readout) → 2 records
"""

    def build_system_prompt(self, field_reference: str) -> str:
        return (
            "You are a scientific literature curation expert specializing in "
            "genetic perturbation experiments (CRISPR screens, Perturb-seq, "
            "deep mutational scanning). Your task is to extract structured metadata "
            "from scientific paper text.\n\n"
            + self.EXPERIMENT_DEFINITION
            + "\n\n"
            + field_reference
            + "\n\nOUTPUT FORMAT:\n"
            "================\n"
            "Return ONLY a valid JSON array. Each element is an object with ALL "
            "of the following keys. Use null (not the string 'null') for fields "
            "that cannot be determined from the paper text. "
            "Do not include fields not listed in the FIELD REFERENCE above. "
            "Do not add prose or explanation outside the JSON.\n\n"
            "Keys (in order): "
            + ", ".join(EXTRACTABLE_FIELDS)
        )

    def build_user_prompt(
        self,
        paper_text: str,
        few_shot_block: str,
        max_paper_chars: int = 80_000,
        additional_instructions: str | None = None,
    ) -> str:
        truncated = paper_text[:max_paper_chars]
        if len(paper_text) > max_paper_chars:
            truncated += "\n\n[... paper truncated ...]"

        parts: list[str] = []
        if additional_instructions:
            parts.append("ADDITIONAL INSTRUCTIONS:\n" + additional_instructions.strip())
        if few_shot_block:
            parts.append("EXAMPLES OF CORRECT EXTRACTION:\n" + few_shot_block)
            parts.append(
                "--- NOW PROCESS THE FOLLOWING PAPER ---\n"
                "Extract ALL experiment records from the paper below. "
                "Return ONLY the JSON array.\n\n"
                "PAPER TEXT:\n" + truncated
            )
        else:
            parts.append(
                "Extract ALL experiment records from the paper below. "
                "Return ONLY the JSON array.\n\n"
                "PAPER TEXT:\n" + truncated
            )

        return "\n\n".join(parts)


# ---------------------------------------------------------------------------
# ExtractionConfig
# ---------------------------------------------------------------------------

_HERE = Path(__file__).resolve().parent
_TEXT_EXTRACTION_DIR = _HERE.parent / "text_extraction"


@dataclass
class ExtractionConfig:
    gemini_model: str = "gemini-2.0-flash"
    curated_gemini_dir: Path = field(
        default_factory=lambda: _TEXT_EXTRACTION_DIR / "curated_gemini"
    )
    pdf_text_dir: Path = field(
        default_factory=lambda: _TEXT_EXTRACTION_DIR / "pdf_text"
    )
    output_dir: Path = field(
        default_factory=lambda: _TEXT_EXTRACTION_DIR / "curated_gemini_v2"
    )
    n_few_shot: int = 3
    max_paper_chars: int = 200_000
    temperature: float = 0.1
    max_output_tokens: int = 65_536
    max_retries: int = 3
    retry_delay_seconds: float = 5.0
    skip_existing: bool = True


# ---------------------------------------------------------------------------
# PaperExtractor
# ---------------------------------------------------------------------------


class ExtractionError(RuntimeError):
    """Raised when the LLM extraction fails after all retries."""


class PaperExtractor:
    """
    Main extraction class. Wraps schema introspection, few-shot loading,
    prompt building, Gemini API calls, and output normalization.

    Requires the GEMINI_API_KEY environment variable to be set, or
    application default credentials if using Vertex AI.

    Example
    -------
    >>> extractor = PaperExtractor()
    >>> records = extractor.extract_from_file(Path("pdf_text/27869803.txt"))
    >>> print(len(records), "experiment records extracted")
    """

    def __init__(self, config: ExtractionConfig | None = None) -> None:
        self.config = config or ExtractionConfig()
        self._client = self._init_client()
        self._field_meta = SchemaIntrospector.get_field_meta()
        self._response_schema = SchemaIntrospector.build_response_schema()
        field_reference = SchemaIntrospector.build_field_reference_block(self._field_meta)
        self._prompt_builder = PromptBuilder()
        self._system_prompt = self._prompt_builder.build_system_prompt(field_reference)
        self._few_shot_provider = FewShotProvider(
            curated_gemini_dir=self.config.curated_gemini_dir,
            pdf_text_dir=self.config.pdf_text_dir,
            n_examples=self.config.n_few_shot,
        )
        self._few_shot_block: str | None = None  # lazy-loaded on first call

    def _init_client(self):
        """Initialize the Gemini client. Deferred import so the module loads without google-genai.

        Supports two auth methods (checked in order):
        1. GEMINI_API_KEY env var       — Google AI Studio API key
        2. GOOGLE_APPLICATION_CREDENTIALS — GCP service account JSON (uses Vertex AI backend)
           Also requires GOOGLE_CLOUD_PROJECT env var, or project_id is read from the JSON.
        """
        try:
            from google import genai
        except ImportError as exc:
            raise ImportError(
                "google-genai is required for PaperExtractor. "
                "Install it with: pip install google-genai>=1.0.0"
            ) from exc
        import os, json as _json
        if os.environ.get("GEMINI_API_KEY"):
            return genai.Client()
        # Service account path: use Vertex AI backend
        sa_path = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS")
        if not sa_path:
            raise EnvironmentError(
                "No credentials found. Set either GEMINI_API_KEY or "
                "GOOGLE_APPLICATION_CREDENTIALS."
            )
        project = os.environ.get("GOOGLE_CLOUD_PROJECT")
        if not project:
            # Read project_id directly from the service account JSON
            with open(sa_path) as f:
                project = _json.load(f).get("project_id")
        if not project:
            raise EnvironmentError(
                "Could not determine GCP project. Set GOOGLE_CLOUD_PROJECT env var."
            )
        return genai.Client(vertexai=True, project=project, location="us-central1")

    def _get_few_shot_block(self) -> str:
        if self._few_shot_block is None:
            examples = self._few_shot_provider.get_examples()
            self._few_shot_block = self._few_shot_provider.format_for_prompt(examples)
        return self._few_shot_block

    def extract(
        self,
        paper_text: str,
        paper_id: str = "unknown",
        additional_instructions: str | None = None,
    ) -> list[dict[str, Any]]:
        """
        Extract experiment records from paper text.

        Parameters
        ----------
        paper_text : str
            Full text of the scientific paper.
        paper_id : str
            Identifier for logging (e.g., PMID or filename stem).
        additional_instructions : str, optional
            Extra instructions prepended to the user prompt, e.g.
            "Only extract Perturb-seq experiments. Ignore CRISPR screens."

        Returns
        -------
        list[dict]
            List of experiment record dicts, one per experiment.
        """
        few_shot_block = self._get_few_shot_block()
        user_prompt = self._prompt_builder.build_user_prompt(
            paper_text=paper_text,
            few_shot_block=few_shot_block,
            max_paper_chars=self.config.max_paper_chars,
            additional_instructions=additional_instructions,
        )
        raw_records = self._call_gemini(self._system_prompt, user_prompt, paper_id)
        return self._validate_and_normalize(raw_records, paper_id)

    def extract_from_file(
        self,
        txt_path: Path | str,
        additional_instructions: str | None = None,
    ) -> list[dict[str, Any]]:
        """Read a paper .txt file and extract experiment records."""
        txt_path = Path(txt_path)
        paper_text = txt_path.read_text(encoding="utf-8", errors="replace")
        return self.extract(
            paper_text,
            paper_id=txt_path.stem,
            additional_instructions=additional_instructions,
        )

    def extract_batch(
        self,
        txt_dir: Path | str | None = None,
        output_dir: Path | str | None = None,
        paper_ids: list[str] | None = None,
    ) -> dict[str, list[dict[str, Any]]]:
        """
        Extract records for all .txt files in txt_dir.

        Parameters
        ----------
        txt_dir : Path, optional
            Directory of paper .txt files. Defaults to config.pdf_text_dir.
        output_dir : Path, optional
            Directory to write output JSON files. Defaults to config.output_dir.
        paper_ids : list[str], optional
            If given, only process these paper IDs (filename stems).

        Returns
        -------
        dict[str, list[dict]]
            Mapping of paper_id -> extracted records for newly processed papers.
        """
        src_dir = Path(txt_dir) if txt_dir else self.config.pdf_text_dir
        out_dir = Path(output_dir) if output_dir else self.config.output_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        txt_files = sorted(src_dir.glob("*.txt"))
        if paper_ids is not None:
            id_set = set(paper_ids)
            txt_files = [f for f in txt_files if f.stem in id_set]

        results: dict[str, list[dict[str, Any]]] = {}
        total = len(txt_files)

        for i, txt_path in enumerate(txt_files, start=1):
            out_path = out_dir / f"{txt_path.stem}.json"
            if self.config.skip_existing and out_path.exists():
                logger.info("[%d/%d] Skipping existing: %s", i, total, txt_path.stem)
                continue

            logger.info("[%d/%d] Extracting: %s", i, total, txt_path.stem)
            try:
                records = self.extract_from_file(txt_path)
                out_path.write_text(
                    json.dumps(records, indent=2, ensure_ascii=False),
                    encoding="utf-8",
                )
                results[txt_path.stem] = records
                logger.info(
                    "  -> %d records written to %s", len(records), out_path.name
                )
            except ExtractionError as exc:
                logger.error("  -> FAILED for %s: %s", txt_path.stem, exc)
            except Exception as exc:
                logger.exception(
                    "  -> Unexpected error for %s: %s", txt_path.stem, exc
                )

        return results

    def _call_gemini(
        self, system_prompt: str, user_prompt: str, paper_id: str = "unknown"
    ) -> list[dict[str, Any]]:
        """Call the Gemini API with retry logic. Returns parsed JSON list."""
        from google import genai
        from google.genai import types

        last_exc: Exception | None = None

        for attempt in range(self.config.max_retries):
            try:
                response = self._client.models.generate_content(
                    model=self.config.gemini_model,
                    contents=[
                        types.Content(
                            role="user",
                            parts=[types.Part(text=user_prompt)],
                        )
                    ],
                    config=types.GenerateContentConfig(
                        system_instruction=system_prompt,
                        response_mime_type="application/json",
                        response_schema=self._response_schema,
                        temperature=self.config.temperature,
                        max_output_tokens=self.config.max_output_tokens,
                    ),
                )
                raw = response.text
                parsed = json.loads(raw)
                if not isinstance(parsed, list):
                    raise ValueError(
                        f"Expected JSON array, got {type(parsed).__name__}: {raw[:200]}"
                    )
                return parsed

            except Exception as exc:
                exc_name = type(exc).__name__
                # Retry on quota/availability errors
                retryable_names = {
                    "ResourceExhausted",
                    "ServiceUnavailable",
                    "InternalServerError",
                    "DeadlineExceeded",
                    "TooManyRequests",
                }
                is_json_truncation = isinstance(exc, json.JSONDecodeError)
                if exc_name in retryable_names or "429" in str(exc) or "503" in str(exc) or is_json_truncation:
                    delay = self.config.retry_delay_seconds * (2**attempt)
                    logger.warning(
                        "[%s] Retryable error on attempt %d/%d: %s. "
                        "Waiting %.1fs before retry.",
                        paper_id,
                        attempt + 1,
                        self.config.max_retries,
                        exc,
                        delay,
                    )
                    time.sleep(delay)
                    last_exc = exc
                else:
                    raise ExtractionError(
                        f"Non-retryable Gemini error for {paper_id}: {exc}"
                    ) from exc

        raise ExtractionError(
            f"Gemini extraction failed for {paper_id} after "
            f"{self.config.max_retries} retries. Last error: {last_exc}"
        ) from last_exc

    def _validate_and_normalize(
        self, records: list[dict[str, Any]], paper_id: str = "unknown"
    ) -> list[dict[str, Any]]:
        """
        Post-process extracted records:
        - Ensure all EXTRACTABLE_FIELDS keys are present (fill missing with None).
        - Warn on constrained-field violations (do not raise).
        - Coerce integer fields to int or None.
        """
        normalized: list[dict[str, Any]] = []
        for i, rec in enumerate(records):
            if not isinstance(rec, dict):
                logger.warning(
                    "[%s] Record %d is not a dict (%s), skipping.",
                    paper_id,
                    i,
                    type(rec).__name__,
                )
                continue

            out: dict[str, Any] = {}
            for field_name in EXTRACTABLE_FIELDS:
                value = rec.get(field_name)

                # Coerce integer fields
                if field_name in _INTEGER_FIELDS and value is not None:
                    try:
                        value = int(str(value).replace(",", "").strip())
                    except (ValueError, TypeError):
                        logger.warning(
                            "[%s] Record %d: could not coerce %s=%r to int.",
                            paper_id,
                            i,
                            field_name,
                            value,
                        )
                        value = None

                # Warn on constraint violations (don't reject)
                if (
                    value is not None
                    and field_name in CONSTRAINED_FIELDS
                    and value not in CONSTRAINED_FIELDS[field_name]
                ):
                    logger.warning(
                        "[%s] Record %d: %s=%r not in allowed values %s.",
                        paper_id,
                        i,
                        field_name,
                        value,
                        CONSTRAINED_FIELDS[field_name],
                    )

                out[field_name] = value

            normalized.append(out)

        return normalized

    # ------------------------------------------------------------------
    # Introspection helpers (useful in notebooks)
    # ------------------------------------------------------------------

    def get_system_prompt(self) -> str:
        """Return the system prompt for inspection."""
        return self._system_prompt

    def get_user_prompt_preview(
        self, paper_text: str = "[PAPER TEXT WOULD GO HERE]"
    ) -> str:
        """Return a user prompt preview (truncated paper_text) for inspection."""
        return self._prompt_builder.build_user_prompt(
            paper_text=paper_text,
            few_shot_block=self._get_few_shot_block(),
            max_paper_chars=self.config.max_paper_chars,
        )

    def get_user_prompt(
        self,
        paper_text: str,
        additional_instructions: str | None = None,
    ) -> str:
        """Return the full user prompt for a given paper text, for inspection."""
        return self._prompt_builder.build_user_prompt(
            paper_text=paper_text,
            few_shot_block=self._get_few_shot_block(),
            max_paper_chars=self.config.max_paper_chars,
            additional_instructions=additional_instructions,
        )


# ---------------------------------------------------------------------------
# Module-level convenience functions
# ---------------------------------------------------------------------------


def extract_paper(
    txt_path: Path | str,
    output_path: Path | str | None = None,
    config: ExtractionConfig | None = None,
) -> list[dict[str, Any]]:
    """
    Extract experiment records from a single paper text file.

    Parameters
    ----------
    txt_path : Path or str
        Path to the paper .txt file.
    output_path : Path or str, optional
        If given, write the output JSON to this path.
    config : ExtractionConfig, optional
        Custom extraction config. Uses defaults if None.

    Returns
    -------
    list[dict]
        Extracted experiment records.
    """
    extractor = PaperExtractor(config)
    records = extractor.extract_from_file(Path(txt_path))
    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        logger.info("Wrote %d records to %s", len(records), out)
    return records


def batch_extract(
    txt_dir: Path | str,
    output_dir: Path | str,
    config: ExtractionConfig | None = None,
    paper_ids: list[str] | None = None,
) -> None:
    """
    Batch-extract experiment records for all .txt files in txt_dir.

    Parameters
    ----------
    txt_dir : Path or str
        Directory containing paper .txt files.
    output_dir : Path or str
        Directory to write output JSON files (one per paper).
    config : ExtractionConfig, optional
        Custom extraction config. Uses defaults if None.
    paper_ids : list[str], optional
        If given, only process papers with these IDs (filename stems).
    """
    extractor = PaperExtractor(config)
    extractor.extract_batch(
        txt_dir=Path(txt_dir),
        output_dir=Path(output_dir),
        paper_ids=paper_ids,
    )
