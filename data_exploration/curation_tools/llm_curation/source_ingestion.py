"""Adapters that snapshot selected external caches into a new curation run."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence

from curation_tools.llm_curation.mavedb.processing import format_identifier_for_lookup
from curation_tools.llm_curation.workflow import (
    CurationItemInput,
    CreateRunRequest,
    MaveDBSnapshotInput,
    SourceDocumentInput,
)

PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_SCHEMA_PATH = PACKAGE_DIR / "llm_curation_schema.py"


def _snapshot_dois(payload: Mapping[str, object]) -> list[str]:
    """Return DOI links recorded directly in a MaveDB snapshot.

    The generated URN-to-DOI map may intentionally contain only primary
    publications, while a selected Markdown source can be a secondary
    publication.  The raw snapshot is therefore the authoritative fallback.
    """
    sources: list[object] = [payload]
    experiment = payload.get("experiment")
    if isinstance(experiment, Mapping):
        sources.append(experiment)
    dois: list[str] = []
    for source in sources:
        for key in ("primaryPublicationIdentifiers", "secondaryPublicationIdentifiers"):
            identifiers = source.get(key) if isinstance(source, Mapping) else None
            if not isinstance(identifiers, list):
                continue
            for identifier in identifiers:
                if not isinstance(identifier, Mapping):
                    continue
                value = identifier.get("doi")
                if not value and str(identifier.get("dbName", "")).casefold() == "doi":
                    value = identifier.get("identifier")
                if isinstance(value, str) and value and value not in dois:
                    dois.append(value)
    return dois


def build_create_run_request(
    *,
    run_name: str,
    publication_dir: str | Path,
    mavedb_metadata_dir: str | Path,
    urn_to_dois_file: str | Path | None,
    step1_prompt: str | Path = PACKAGE_DIR / "step1_evidence_extraction_prompt.md",
    step2_prompt: str | Path = PACKAGE_DIR / "step2_specific_term_extraction.md",
    step3_prompt: str | Path = PACKAGE_DIR / "step3_candidate_discovery_prompt.md",
    schema_path: str | Path = DEFAULT_SCHEMA_PATH,
    model: str = "google/gemini-3.7-flash",
    max_workers: int = 8,
    selected_urns: Sequence[str] | None = None,
) -> CreateRunRequest:
    """Build an immutable-create request from selected source cache locations."""
    publications = Path(publication_dir).resolve()
    metadata_dir = Path(mavedb_metadata_dir).resolve()
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    if not publications.is_dir():
        raise FileNotFoundError(f"Publication directory not found: {publications}")
    if not metadata_dir.is_dir():
        raise FileNotFoundError(f"MaveDB metadata directory not found: {metadata_dir}")
    publication_paths = sorted(publications.glob("*.md"))
    if not publication_paths:
        raise ValueError(f"No Markdown publication sources found in: {publications}")
    snapshots: list[MaveDBSnapshotInput] = []
    for path in sorted(metadata_dir.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid MaveDB JSON input {path}: {exc}") from exc
        if not isinstance(payload, dict):
            raise ValueError(f"MaveDB snapshot must be an object: {path}")
        urn = str(payload.get("urn") or payload.get("urnMaveDb") or "")
        if not urn:
            raise ValueError(f"MaveDB snapshot has no URN: {path}")
        snapshots.append(MaveDBSnapshotInput(urn, path.name, payload, path))
    if not snapshots:
        raise ValueError(f"No MaveDB JSON snapshots found in: {metadata_dir}")
    selected = {urn.strip() for urn in selected_urns or () if urn.strip()}
    if selected:
        available = {snapshot.urn for snapshot in snapshots}
        unknown = selected - available
        if unknown:
            raise ValueError(f"Selected MaveDB URNs were not found: {sorted(unknown)}")
        snapshots = [snapshot for snapshot in snapshots if snapshot.urn in selected]

    mapping_path = Path(urn_to_dois_file).resolve() if urn_to_dois_file else None
    if mapping_path is None:
        raise ValueError(
            "A URN-to-DOI mapping is required to create linked curation items"
        )
    try:
        mapping: Mapping[str, object] = json.loads(
            mapping_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid URN-to-DOI mapping {mapping_path}: {exc}") from exc
    if not isinstance(mapping, dict):
        raise ValueError("URN-to-DOI mapping must be a JSON object")
    documents: list[SourceDocumentInput] = []
    document_labels_by_urn: dict[str, list[str]] = {}
    for path in publication_paths:
        label = path.name
        matching_sources = [
            (snapshot, doi)
            for snapshot in snapshots
            for doi in [
                *(
                    [str(value) for value in mapping.get(snapshot.urn, [])]
                    if isinstance(mapping.get(snapshot.urn), list)
                    else []
                ),
                *_snapshot_dois(snapshot.payload),
            ]
            if format_identifier_for_lookup(doi) == path.stem
        ]
        matching_urns = tuple(
            dict.fromkeys(snapshot.urn for snapshot, _ in matching_sources)
        )
        if not matching_urns:
            if selected:
                continue
            raise ValueError(
                f"No MaveDB snapshot is linked to publication source: {label}. "
                "The supplied mapping and selected MaveDB metadata contain no matching DOI."
            )
        doi = matching_sources[0][1]
        documents.append(
            SourceDocumentInput(label, path.read_text(encoding="utf-8"), path, doi)
        )
        for urn in matching_urns:
            document_labels_by_urn.setdefault(urn, []).append(label)

    items: list[CurationItemInput] = []
    unlinked_selected_urns: list[str] = []
    for snapshot in snapshots:
        document_labels = tuple(
            dict.fromkeys(document_labels_by_urn.get(snapshot.urn, ()))
        )
        if not document_labels:
            if selected:
                unlinked_selected_urns.append(snapshot.urn)
            continue
        items.append(
            CurationItemInput(
                item_label=snapshot.urn,
                document_labels=document_labels,
                source_urns=(snapshot.urn,),
                prompt_context={
                    "publication_source": document_labels[0],
                    "publication_sources": document_labels,
                },
                primary_dataset_id=snapshot.urn,
            )
        )
    if unlinked_selected_urns:
        raise ValueError(
            "Selected MaveDB URNs have no linked publication source: "
            f"{sorted(unlinked_selected_urns)}"
        )
    if not items:
        raise ValueError("No publication sources match the selected MaveDB URNs")
    prompt_paths = {
        "step1": Path(step1_prompt),
        "step2": Path(step2_prompt),
        "step3": Path(step3_prompt),
    }
    return CreateRunRequest(
        run_name=run_name,
        documents=tuple(documents),
        mavedb_snapshots=tuple(snapshots),
        items=tuple(items),
        prompt_templates={
            name: path.resolve().read_text(encoding="utf-8")
            for name, path in prompt_paths.items()
        },
        model_settings={"model": model, "max_workers": max_workers},
        schema_path=Path(schema_path).resolve(),
        configuration={
            "creation_sources": {
                "publication_dir": str(publications),
                "mavedb_metadata_dir": str(metadata_dir),
                "urn_to_dois_file": str(mapping_path) if mapping_path else None,
                "selected_urns": sorted(selected),
            }
        },
    )
