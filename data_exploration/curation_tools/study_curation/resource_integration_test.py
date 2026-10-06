"""Integration checks for packaged prompts and custom schema snapshots."""

import hashlib
import json
from pathlib import Path

import pytest

from curation_tools.study_curation import paths
from curation_tools.study_curation.llm.source_ingestion import build_create_run_request
from curation_tools.study_curation.llm.workflow import CurationWorkflow


def test_run_uses_packaged_prompts_and_preserves_schema_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    schema_source = paths.DEFAULT_SCHEMA_PATH.read_text(encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    publications = tmp_path / "publications"
    metadata = tmp_path / "metadata"
    publications.mkdir()
    metadata.mkdir()
    (publications / "10_1234_example.md").write_text("A MAVE publication.")
    (metadata / "record.json").write_text(json.dumps({"urn": "urn:mavedb:1"}))
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"urn:mavedb:1": ["10.1234/example"]}))

    schema_path = tmp_path / "custom_schema.py"
    schema_path.write_text(schema_source, encoding="utf-8")
    request = build_create_run_request(
        run_name="resources",
        publication_dir=publications,
        mavedb_metadata_dir=metadata,
        urn_to_dois_file=mapping,
        schema_path=schema_path,
    )
    assert all(request.prompt_templates.values())
    workflow = CurationWorkflow(tmp_path / "runs")
    workflow.create_run(request)
    original_configuration = workflow.run_configuration("resources")

    expected_hash = hashlib.sha256(schema_path.read_bytes()).hexdigest()
    assert workflow.preview_schema_update("resources")["schema_hash"] == expected_hash
    assert workflow.apply_schema_terms("resources", expected_hash).status == "applied"
    assert workflow.run_configuration("resources") == original_configuration
    assert schema_path.read_text(encoding="utf-8") == schema_source
