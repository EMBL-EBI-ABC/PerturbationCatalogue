"""Pipeline for discovering ontology candidate terms from 'Other' evidence across the corpus (Step 3a)."""

import argparse
import json
import os
import traceback
import types
import typing
from pathlib import Path
from typing import Any, Type, get_args, get_origin

import instructor
from instructor.cache import DiskCache
from pydantic import BaseModel

from curation_tools.llm_curation.llm_curation_schema import (
    FieldCandidates,
    SpecificTermExtractionSchema,
)
from curation_tools.llm_curation.logging_utils import (
    _ensure_log_file,
    append_log_line,
    print_status_block,
)

# Persistent disk-based instructor cache
cache = DiskCache(directory=".instructor_cache")

DEFAULT_LLM_MODEL_NAME = os.getenv("LLM_MODEL_NAME", "google/gemini-3.6-flash")
JSON_INDENT = 2


def get_controlled_vocabulary_and_description(field_name: str) -> tuple[list[str], str]:
    """Programmatically extract the allowed controlled vocabulary list and field description from the schema."""
    field = SpecificTermExtractionSchema.model_fields.get(field_name)
    if not field:
        return [], ""

    description = field.description or ""
    annotation = field.annotation

    def extract_literals(annot) -> list[str]:
        origin = get_origin(annot)
        if origin in (typing.Union, types.UnionType):
            literals = []
            for arg in get_args(annot):
                literals.extend(extract_literals(arg))
            return literals
        elif origin is typing.Literal:
            return list(get_args(annot))
        return []

    return extract_literals(annotation), description


def aggregate_unmapped_evidence(
    step1_dir: Path,
    step2_dir: Path,
    log_file: Path,
    selected_files: list[str | Path] | None = None,
) -> dict[str, list[dict[str, str]]]:
    """Iterate over Step 2 normalized artifacts, find fields with value "Other", and retrieve their original verbatim evidence string from corresponding Step 1 files."""
    unmapped_evidence_map: dict[str, list[dict[str, str]]] = {}

    step1_dir = Path(step1_dir).resolve()
    step2_dir = Path(step2_dir).resolve()

    if not step2_dir.is_dir():
        print_status_block(
            log_file,
            "Step 2 directory not found",
            f"Path: {step2_dir}",
        )
        return unmapped_evidence_map

    step2_files = [
        f
        for f in sorted(step2_dir.glob("*.json"))
        if not f.name.endswith("_audit.json")
    ]
    if not step2_files and (step2_dir / "step2_normalized").is_dir():
        step2_files = [
            f
            for f in sorted((step2_dir / "step2_normalized").glob("*.json"))
            if not f.name.endswith("_audit.json")
        ]
    if not step2_files:
        step2_files = [
            f
            for f in sorted(step2_dir.rglob("*.json"))
            if not f.name.endswith("_audit.json")
        ]

    if selected_files is not None:
        selected_names = {Path(file_name).name for file_name in selected_files}
        step2_files = [
            file_path for file_path in step2_files if file_path.name in selected_names
        ]

    append_log_line(
        log_file,
        f"Found {len(step2_files)} Step 2 files to analyze for 'Other' values.",
    )

    for s2_file in step2_files:
        s1_file = step1_dir / s2_file.name
        if (
            not s1_file.is_file()
            and (step1_dir / "step1_evidence" / s2_file.name).is_file()
        ):
            s1_file = step1_dir / "step1_evidence" / s2_file.name
        if not s1_file.is_file():
            matches = list(step1_dir.rglob(s2_file.name))
            if matches:
                s1_file = matches[0]

        if not s1_file.is_file():
            append_log_line(
                log_file,
                f"Corresponding Step 1 file not found for {s2_file.name} in {step1_dir}. Skipping.",
            )
            continue

        try:
            s2_data = json.loads(s2_file.read_text(encoding="utf-8"))
            s1_data = json.loads(s1_file.read_text(encoding="utf-8"))
        except Exception as e:
            append_log_line(
                log_file,
                f"Error parsing JSON for {s2_file.name}: {e}. Skipping.",
            )
            continue

        for field_name, val in s2_data.items():
            if val == "Other":
                # Construct step 1 evidence key (e.g., cell_line_label -> cell_line_label_evidence)
                s1_key = f"{field_name}_evidence"
                evidence_val = s1_data.get(s1_key) or s1_data.get(field_name)

                if isinstance(evidence_val, str):
                    cleaned_evidence = evidence_val.strip()
                    # Filter out empty or obvious null/none placeholders
                    if cleaned_evidence and cleaned_evidence.lower() not in (
                        "null",
                        "none",
                        "n/a",
                        "na",
                        "",
                    ):
                        unmapped_evidence_map.setdefault(field_name, []).append(
                            {
                                "evidence": cleaned_evidence,
                                "source_file": s2_file.name,
                            }
                        )

    return unmapped_evidence_map


def discover_candidates(
    step1_dir: Path,
    step2_dir: Path,
    output_dir: Path,
    log_file: Path,
    prompt_template_file: Path,
    model_name: str,
    verbose: bool = False,
    selected_files: list[str | Path] | None = None,
) -> Path:
    """Run corpus-level ontology candidate discovery on fields with 'Other' values."""
    output_dir = Path(output_dir).resolve()
    log_file = _ensure_log_file(log_file)
    prompt_template_file = Path(prompt_template_file).resolve()

    print_status_block(
        log_file,
        "Starting Step 3a Ontology Candidate Discovery",
        f"Step 1 Dir: {step1_dir}",
        f"Step 2 Dir: {step2_dir}",
        f"Output Dir: {output_dir}",
        f"Model Name: {model_name}",
    )

    # Aggregate unmapped evidence from Step 1 and Step 2 files
    unmapped_evidence = aggregate_unmapped_evidence(
        step1_dir, step2_dir, log_file, selected_files=selected_files
    )

    if not unmapped_evidence:
        print_status_block(
            log_file,
            "Ontology Candidate Discovery Complete (No candidates)",
            "No 'Other' values found across any fields in the analyzed Step 2 files.",
        )
        out_file = output_dir / "step3_ontology_candidates.json"
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_file.write_text(json.dumps({}, indent=JSON_INDENT), encoding="utf-8")
        return out_file

    prompt_template = prompt_template_file.read_text(encoding="utf-8")

    client = instructor.from_provider(
        model_name,
        location="global",
        vertexai=True,
        cache=cache,
    )

    results: dict[str, list[dict[str, Any]]] = {}

    for field_name, evidence_list in unmapped_evidence.items():
        append_log_line(
            log_file,
            f"Analyzing field '{field_name}' with {len(evidence_list)} unique unmapped evidence strings.",
        )

        vocab, desc = get_controlled_vocabulary_and_description(field_name)
        controlled_vocab_str = ""
        if desc:
            controlled_vocab_str += f"Description: {desc}\n"
        if vocab:
            controlled_vocab_str += f"Allowed Terms: {json.dumps(vocab)}"
        else:
            controlled_vocab_str += "Allowed Terms: [No controlled vocabulary list defined, open text mapping expected]"

        prompt = prompt_template.format(
            field_name=field_name,
            controlled_vocabulary=controlled_vocab_str,
            evidence_list=evidence_list,
        )

        if verbose:
            print_status_block(
                log_file,
                f"[VERBOSE] Step 3 Candidate Discovery Prompt for field: {field_name}",
                prompt,
            )

        try:
            response = client.create(
                response_model=FieldCandidates,
                messages=[{"role": "user", "content": prompt}],
                thinking_config={
                    "thinking_level": "high",
                },
                generation_config={
                    "temperature": 0.2,
                },
            )

            # Validate and store the results
            if isinstance(response, FieldCandidates):
                results[field_name] = [
                    {
                        "proposed_new_term": candidate.proposed_new_term,
                        "supporting_evidence": [
                            c.model_dump() for c in candidate.supporting_evidence
                        ],
                        "rationale": candidate.rationale,
                    }
                    for candidate in response.candidates
                ]

        except Exception as e:
            append_log_line(
                log_file,
                f"Error calling LLM for field '{field_name}': {e}\n{traceback.format_exc()}",
            )

    output_path = output_dir / "step3_ontology_candidates.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(results, indent=JSON_INDENT), encoding="utf-8")

    print_status_block(
        log_file,
        "Step 3 Ontology Candidate Discovery Complete",
        f"Consolidated candidates report saved: {output_path}",
    )

    return output_path


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for Step 3a ontology candidate discovery."""
    parser = argparse.ArgumentParser(
        description="Analyze unmapped 'Other' evidence across Step 2 outputs to discover missing ontology candidates (Step 3a).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--step1-dir",
        type=Path,
        required=True,
        help="Directory containing Step 1 evidence JSON files.",
    )
    parser.add_argument(
        "--step2-dir",
        type=Path,
        required=True,
        help="Directory containing Step 2 normalized JSON files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory used to write Step 3 candidates report.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        required=True,
        help="Log file used for progress and errors.",
    )
    parser.add_argument(
        "--prompt-template-file",
        type=Path,
        default=Path(__file__).parent / "step3_candidate_discovery_prompt.md",
        help="Prompt template file used for candidate discovery.",
    )
    parser.add_argument(
        "--llm-model",
        default=DEFAULT_LLM_MODEL_NAME,
        help="LLM model ID to use.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print verbose prompt logs.",
    )
    return parser


def main() -> None:
    """Main CLI entry point for Step 3 ontology candidate discovery."""
    args = build_parser().parse_args()
    discover_candidates(
        step1_dir=args.step1_dir,
        step2_dir=args.step2_dir,
        output_dir=args.output_dir,
        log_file=args.log_file,
        prompt_template_file=args.prompt_template_file,
        model_name=args.llm_model,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
