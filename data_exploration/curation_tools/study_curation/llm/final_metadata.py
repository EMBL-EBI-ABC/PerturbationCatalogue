"""Pure Step 5 projection from persisted artifacts to final metadata payloads."""

from __future__ import annotations

from typing import Any

from curation_tools.perturbseq_anndata_schema import ObsSchema
from curation_tools.study_curation.sources.mavedb import (
    extract_curated_mavedb_prompt_metadata,
    merge_prompt_metadata_value,
)

OBS_SCHEMA_FIELDS = tuple(ObsSchema.to_schema().columns.keys())
PROVENANCE_FIELDS = ("__source_urns", "__source_files")


def _is_missing(value: object) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _unique_strings(values: list[object]) -> list[str]:
    result: list[str] = []
    for value in values:
        if isinstance(value, str) and value and value not in result:
            result.append(value)
    return result


def _join_values(values: list[object]) -> str | None:
    unique_values = _unique_strings(values)
    return "|".join(unique_values) if unique_values else None


def _extract_mavedb_fields(entries: list[dict[str, Any]]) -> dict[str, object]:
    """Map stored MaveDB snapshots into fields defined by ``ObsSchema``."""
    target_records: list[dict[str, Any]] = []
    prompt_metadata: dict[str, object] = {}
    for entry in entries:
        target_records.extend(
            target
            for target in entry.get("targetGenes", [])
            if isinstance(target, dict)
        )
        for field_name, field_value in extract_curated_mavedb_prompt_metadata(
            entry
        ).items():
            prompt_metadata[field_name] = merge_prompt_metadata_value(
                prompt_metadata.get(field_name), field_value
            )
    target_names = _unique_strings([target.get("name") for target in target_records])
    taxonomies = [
        ((target.get("targetSequence") or {}).get("taxonomy") or {}).get("organismName")
        for target in target_records
    ]
    result: dict[str, object] = {}
    if target_names:
        result["perturbed_target_symbol"] = _join_values(target_names)
        result["perturbed_target_number"] = len(target_names)
        result["number_of_perturbed_targets"] = str(len(target_names))
    if prompt_metadata.get("total_variants") is not None:
        result["library_total_variants"] = prompt_metadata["total_variants"]
    for source_name, target_name in (
        ("primary_publication_titles", "study_title"),
        ("primary_publication_dois", "study_uri"),
    ):
        value = prompt_metadata.get(source_name)
        joined = _join_values(value if isinstance(value, list) else [value])
        if joined:
            result[target_name] = joined
    years = prompt_metadata.get("primary_publication_years")
    if isinstance(years, list) and years:
        result["study_year"] = years[0]
    for source_name, target_name in (
        ("experiment_title", "experiment_title"),
        ("license", "license_label"),
    ):
        if prompt_metadata.get(source_name):
            result[target_name] = prompt_metadata[source_name]
    summary = _join_values(
        [
            prompt_metadata.get("experiment_short_description"),
            prompt_metadata.get("score_set_short_description"),
        ]
    )
    if summary:
        result["experiment_summary"] = summary
    if any(taxonomy == "Homo sapiens" for taxonomy in taxonomies):
        result["species"] = "Homo sapiens"
    if prompt_metadata:
        result["data_modality"] = "MAVE"
    return result


def merge_final_metadata_from_entries(
    normalized_payload: dict[str, Any], source_entries: list[dict[str, Any]]
) -> dict[str, object]:
    """Project a stored normalised payload and stored MaveDB snapshots to ObsSchema."""
    final_payload = {
        field_name: normalized_payload.get(field_name)
        for field_name in OBS_SCHEMA_FIELDS
    }
    for field_name, value in _extract_mavedb_fields(source_entries).items():
        if field_name in final_payload and _is_missing(final_payload[field_name]):
            final_payload[field_name] = value
    for field_name in PROVENANCE_FIELDS:
        if field_name in normalized_payload:
            final_payload[field_name] = normalized_payload[field_name]
    return final_payload
