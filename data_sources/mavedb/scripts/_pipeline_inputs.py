"""Reuse canonical processing before acquiring any upstream source snapshots."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

from data_sources.mavedb.scripts._pipeline_adopt import adopt_processing
from data_sources.mavedb.scripts._pipeline_cache import (
    CacheError,
    ResultStore,
    file_identity,
)


def _same_source(requested: str, original: str) -> bool:
    """Compare recorded file URIs with their equivalent local CLI paths."""

    def normalized(source: str) -> str:
        parsed = urlparse(source)
        if parsed.scheme == "file":
            return str(Path(unquote(parsed.path)).expanduser().resolve())
        if not parsed.scheme:
            return str(Path(source).expanduser().resolve())
        return source

    return normalized(str(requested)) == normalized(str(original))


def _identity(value: Any, label: str) -> dict[str, Any]:
    if (
        not isinstance(value, dict)
        or not isinstance(value.get("bytes"), int)
        or isinstance(value.get("bytes"), bool)
        or value["bytes"] < 0
        or not isinstance(value.get("sha256"), str)
        or re.fullmatch(r"[a-fA-F0-9]{64}", value["sha256"]) is None
    ):
        raise CacheError(
            f"Saved processing has invalid {label} provenance; use --rerun process"
        )
    return {key: value[key] for key in ("bytes", "sha256")}


def validate_adopted_processing(
    p: dict[str, Any], path: Path, manifest: dict[str, Any]
) -> dict[str, Any]:
    """Check requested source identities without network access or source reads
    when the original files are no longer available.

    The caller verifies the completed canonical output bytes separately. Return
    portable provenance identities suitable for an adopted process result whose
    original raw-file inventory was never registered in this result store.
    """
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1:
        raise CacheError(
            f"Unsupported processing provenance in {path}; use --rerun process"
        )
    ingest = manifest.get("ingest")
    taxonomy = manifest.get("taxonomy")
    if not isinstance(ingest, dict) or not isinstance(taxonomy, dict):
        raise CacheError(
            f"Missing processing source provenance in {path}; use --rerun process"
        )
    archive_hash = ingest.get("archive_sha256")
    if (
        not isinstance(archive_hash, str)
        or re.fullmatch(r"[a-fA-F0-9]{64}", archive_hash) is None
    ):
        raise CacheError(
            "Saved processing has invalid archive identity; use --rerun process"
        )
    metadata = _identity(ingest.get("metadata"), "source metadata")
    tax_identity = _identity(taxonomy, "taxonomy")
    expected = p.get("archive_sha256")
    if expected and expected.lower() != archive_hash.lower():
        raise CacheError(
            "Requested archive checksum differs from saved processing; use matching inputs or --rerun ingest,process"
        )
    archive = p.get("archive")
    if archive and Path(archive).is_file():
        if file_identity(Path(archive))["sha256"] != archive_hash:
            raise CacheError(
                "Source archive differs from saved processing; use --rerun ingest,process"
            )
    else:
        requested = archive or p.get("zenodo_url")
        original = ingest.get("source")
        if (
            not requested
            or not isinstance(original, str)
            or not _same_source(requested, original)
        ):
            raise CacheError(
                "Requested archive source differs from saved processing; use original settings or --rerun ingest,process"
            )
    taxonomy_path = p.get("taxonomy")
    if taxonomy_path and Path(taxonomy_path).is_file():
        if file_identity(Path(taxonomy_path)) != tax_identity:
            raise CacheError(
                "Taxonomy differs from saved processing; use --rerun taxonomy,process"
            )
    else:
        requested = taxonomy_path or p.get("taxonomy_url")
        original = taxonomy.get("source")
        if (
            not requested
            or not isinstance(original, str)
            or not _same_source(requested, original)
        ):
            raise CacheError(
                "Requested taxonomy source differs from saved processing; use original settings or --rerun taxonomy,process"
            )
    return {"archive": archive_hash, "metadata": metadata, "taxonomy": tax_identity}


def _manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads((path / "processing_manifest.json").read_text())
    except (OSError, UnicodeError, ValueError) as error:
        raise CacheError(
            f"Cannot read saved processing provenance in {path}; use --rerun process"
        ) from error
    if not isinstance(value, dict):
        raise CacheError(
            f"Invalid processing provenance in {path}; use --rerun process"
        )
    return value


def validate_registered_sources(
    store: ResultStore, record: dict[str, Any], provenance: dict[str, Any]
) -> None:
    """Detect replacement upstream generations using their saved inventories.

    Reading a completion record is enough here: canonical output reuse does not
    require rehashing the complete raw dump or staging it into a new task.
    """
    ingest = store.record("ingest")
    if ingest is not None:
        if "raw" in record["inputs"]:
            current_raw = {
                name.removeprefix("mavedb_dump/"): identity
                for name, identity in ingest["outputs"].items()
                if name.startswith("mavedb_dump/")
            }
            if current_raw != record["inputs"]["raw"]:
                raise CacheError(
                    "Saved source data changed after processing; use --rerun process"
                )
        else:
            # Legacy processing knows the archive hash, rather than a complete
            # per-file raw inventory. Verify the small registered source receipt.
            receipt = store.root / ingest["result_path"] / "manifest.json"
            try:
                if file_identity(receipt) != ingest["outputs"].get("manifest.json"):
                    raise CacheError(
                        "Registered archive receipt differs; use --rerun ingest"
                    )
                source = json.loads(receipt.read_text())
            except (OSError, UnicodeError, ValueError) as error:
                raise CacheError(
                    "Cannot validate registered source archive; use --rerun ingest"
                ) from error
            if source.get("archive_sha256") != provenance["archive"]:
                raise CacheError(
                    "Saved source archive changed after processing; use --rerun process"
                )
    taxonomy = store.record("taxonomy")
    if (
        taxonomy is not None
        and taxonomy["outputs"].get("taxcat.dmp") != provenance["taxonomy"]
    ):
        raise CacheError("Saved taxonomy changed after processing; use --rerun process")


def reusable_processing(p: dict[str, Any], store: ResultStore) -> str | None:
    """Resolve completed canonical results before scheduling upstream stages.

    Derived tables require their matching raw snapshot.
    Explicit source regeneration also goes through the ordinary stage workflow.
    """
    if (
        p.get("ingest_only")
        or p.get("processed_dir")
        or "tables" in p.get("enabled", [])
        or set(p.get("selected", [])) & {"ingest", "taxonomy", "process"}
    ):
        return None
    record = store.record("process")
    if record is not None:
        if record["settings"].get("external"):
            return None
        path = store.get("process")
        manifest = _manifest(path)
        provenance = validate_adopted_processing(p, path, manifest)
        validate_registered_sources(store, record, provenance)
        return str(path)
    legacy = store.root / "processed"
    if not legacy.exists():
        if legacy.is_symlink():
            raise CacheError(
                "Saved processing publication is broken; use --rerun process"
            )
        return None
    try:
        manifest = adopt_processing(legacy)
        provenance = validate_adopted_processing(p, legacy, manifest)
        # Do not attach an old canonical result to a newly replaced source.
        validate_registered_sources(store, {"inputs": provenance}, provenance)
    except (OSError, ValueError, KeyError) as error:
        raise CacheError(
            f"Cannot adopt saved processing: {error}. Use --rerun process"
        ) from error
    return str(store.commit("process", legacy, {"adopted": True}, provenance))
