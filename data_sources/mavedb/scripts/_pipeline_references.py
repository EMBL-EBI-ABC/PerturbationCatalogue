"""Validation and provenance helpers for the unified references stage."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

from data_sources.mavedb.scripts._pipeline_cache import CacheError, file_identity

VARIANT_FILES = {
    "clinvar": "mave_clinvar_variants.tsv",
    "gwas": "mave_gwas_variants.tsv",
}
VARIANT_SOURCE_NAMES = ("coding_variants", "evidence", "credible_sets", "id_map")
VARIANT_REQUIRED = {
    "clinvar": {
        "uniprotAccessions",
        "clinicalSignificances",
        "aminoAcidPosition",
        "referenceAminoAcid",
        "alternateAminoAcid",
    },
    "gwas": {
        "studyId",
        "beta",
        "pValueMantissa",
        "pValueExponent",
        "uniprotAccessions",
        "aminoAcidPosition",
        "referenceAminoAcid",
        "alternateAminoAcid",
    },
}


def read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CacheError(f"Cannot read {label} at {path}: {error}") from error
    if not isinstance(value, dict):
        raise CacheError(f"Invalid {label} at {path}: expected a JSON object")
    return value


def validate_tsv(path: Path, name: str) -> tuple[Path, int]:
    """Verify a prepared validation table's header and count its data rows."""
    path = Path(path)
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise CacheError(f"--validation_{name} file does not exist: {path}") from error
    if not resolved.is_file():
        raise CacheError(f"--validation_{name} must be a regular TSV file: {resolved}")
    try:
        if resolved.stat().st_size == 0:
            raise CacheError(f"--validation_{name} is empty: {resolved}")
        with resolved.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header:
                raise CacheError(f"--validation_{name} has no TSV header: {resolved}")
            if len(header) != len(set(header)):
                raise CacheError(f"--validation_{name} has duplicate TSV column names")
            missing = VARIANT_REQUIRED[name] - set(header)
            if missing:
                raise CacheError(
                    f"--validation_{name} is missing required TSV columns: "
                    + ", ".join(sorted(missing))
                )
            rows = sum(1 for _ in reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise CacheError(
            f"Cannot read --validation_{name} TSV at {resolved}: {error}"
        ) from error
    return resolved, rows


def _identity_record(value: Any, label: str) -> dict[str, int | str]:
    if not isinstance(value, dict):
        raise CacheError(f"Invalid source identity for {label}")
    size = value.get("bytes")
    digest = value.get("sha256")
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise CacheError(f"Invalid byte count for source {label}")
    if not isinstance(digest, str) or re.fullmatch(r"[a-fA-F0-9]{64}", digest) is None:
        raise CacheError(f"Invalid SHA-256 digest for source {label}")
    return {"bytes": size, "sha256": digest.lower()}


def _source_entry_identity(value: Any, label: str) -> dict[str, Any]:
    """Extract the content identities of one v3 source record."""
    if not isinstance(value, dict):
        raise CacheError(f"Source manifest entry {label} must be an object")
    files = value.get("files")
    if files is not None:
        if not isinstance(files, list) or not files:
            raise CacheError(f"Source manifest entry {label} requires file records")
        identities = []
        for index, descriptor in enumerate(files):
            identity = _identity_record(descriptor, f"{label}[{index}]")
            source_path = (
                descriptor.get("path") if isinstance(descriptor, dict) else None
            )
            if not isinstance(source_path, str) or not source_path:
                raise CacheError(f"Source manifest entry {label}[{index}] needs a path")
            identities.append({"path": source_path, **identity})
        return {"files": identities}
    identity = _identity_record(value, label)
    source_path = value.get("path")
    if not isinstance(source_path, str) or not source_path:
        raise CacheError(f"Source manifest entry {label} needs a path")
    return {"file": {"path": source_path, **identity}}


def _source_content_identity(value: Any, label: str) -> dict[str, Any]:
    """Compare a source record's bytes while ignoring its storage location."""
    if not isinstance(value, dict):
        raise CacheError(f"Source manifest entry {label} must be an object")
    files = value.get("files")
    if files is not None:
        if not isinstance(files, list) or not files:
            raise CacheError(f"Source manifest entry {label} requires file records")
        contents = [_identity_record(item, f"{label} partition") for item in files]
        return {
            "files": sorted(contents, key=lambda item: (item["sha256"], item["bytes"]))
        }
    if isinstance(value.get("file"), dict):
        return _identity_record(value["file"], label)
    return _identity_record(value, label)


def source_manifest_identities(path: Path, verify_files: bool = True) -> dict[str, Any]:
    """Read and, by default, verify the selected v3 inputs for both variant tables."""
    from data_sources.mavedb.analysis._reference_io import verify_file
    from data_sources.mavedb.analysis._reference_sources import load_sources

    manifest_path = Path(path).expanduser().resolve(strict=True)
    if not manifest_path.is_file():
        raise CacheError(f"--reference_sources_manifest must be a regular file: {path}")
    raw_manifest = read_json(manifest_path, "raw reference source manifest")
    try:
        manifest = load_sources(manifest_path)
    except (OSError, UnicodeError, ValueError, KeyError) as error:
        raise CacheError(
            f"Cannot validate raw source manifest {manifest_path}: {error}"
        ) from error
    if manifest.get("schema_version") != 3:
        raise CacheError(f"Unsupported raw reference source manifest: {manifest_path}")
    sources = manifest.get("sources")
    raw_sources = raw_manifest.get("sources")
    if not isinstance(sources, dict) or not isinstance(raw_sources, dict):
        raise CacheError(
            f"Raw reference source manifest is missing sources: {manifest_path}"
        )
    missing = sorted(set(VARIANT_SOURCE_NAMES) - set(sources))
    if missing:
        raise CacheError(
            "Raw reference source manifest is missing required source(s): "
            + ", ".join(missing)
        )
    required_sources = {}
    for name in VARIANT_SOURCE_NAMES:
        required_sources[name] = _source_entry_identity(raw_sources[name], name)
        if verify_files:
            record = sources[name]
            descriptors = record.get("files")
            try:
                if isinstance(descriptors, list):
                    for descriptor in descriptors:
                        verify_file(descriptor)
                else:
                    verify_file(record)
            except (OSError, ValueError, KeyError) as error:
                raise CacheError(
                    f"Raw reference source {name} failed verification: {error}"
                ) from error
    return {
        "schema_version": 3,
        "sources": required_sources,
    }


def taxonomy_manifest_identity(taxonomy_dir: Path) -> dict[str, int | str]:
    """Verify a saved taxonomy directory against its own manifest."""
    taxonomy_dir = Path(taxonomy_dir)
    taxonomy = taxonomy_dir / "taxcat.dmp"
    manifest_path = taxonomy_dir / "manifest.json"
    if taxonomy.is_symlink() or manifest_path.is_symlink():
        raise CacheError(f"Saved taxonomy contains a symbolic link: {taxonomy_dir}")
    actual = file_identity(taxonomy)
    manifest = read_json(manifest_path, "taxonomy manifest")
    if manifest.get("schema_version") not in {None, 1} or manifest.get("path") not in {
        None,
        "taxcat.dmp",
    }:
        raise CacheError(f"Invalid saved taxonomy manifest: {manifest_path}")
    recorded = _identity_record(manifest, "taxonomy")
    if actual != recorded:
        raise CacheError(f"Saved taxonomy differs from its manifest: {taxonomy_dir}")
    return actual


def _validate_prepared_variants(bundle: Path) -> dict[str, Any]:
    manifest_path = bundle / "variants" / "manifest.json"
    manifest = read_json(manifest_path, "prepared variant manifest")
    if manifest.get("schema_version") != 1:
        raise CacheError(f"Unsupported prepared variant manifest: {manifest_path}")
    results = manifest.get("results")
    if not isinstance(results, dict):
        raise CacheError(
            f"Prepared variant manifest is missing results: {manifest_path}"
        )
    for name, filename in VARIANT_FILES.items():
        entry = results.get(name)
        if not isinstance(entry, dict):
            raise CacheError(f"Prepared variant manifest is missing {name} output")
        if entry.get("file") != filename:
            raise CacheError(f"Prepared {name} output must be named {filename}")
        expected = _identity_record(entry, f"prepared {name} output")
        actual = file_identity(bundle / "variants" / filename)
        if actual != expected:
            raise CacheError(f"Prepared {name} output differs from its manifest")
        _, rows = validate_tsv(bundle / "variants" / filename, name)
        if "rows" in entry and entry["rows"] != rows:
            raise CacheError(f"Prepared {name} row count differs from its manifest")
    return manifest


def validate_variant_component(bundle: Path) -> dict[str, Any]:
    """Verify a saved variants directory before reusing it as a component."""
    return _validate_prepared_variants(Path(bundle))


def _validate_raw_source_manifest(path: Path, verify_files: bool) -> dict[str, Any]:
    return source_manifest_identities(path, verify_files=verify_files)


def validate_reference_bundle(
    bundle: Path, verify_sources: bool = True
) -> dict[str, Any]:
    """Validate the complete public references bundle before it is committed."""
    bundle = Path(bundle)
    top = read_json(bundle / "manifest.json", "references bundle manifest")
    if top.get("schema_version") != 1:
        raise CacheError(
            f"Unsupported references bundle manifest: {bundle / 'manifest.json'}"
        )
    taxonomy_identity = taxonomy_manifest_identity(bundle / "taxonomy")
    taxonomy_manifest = read_json(
        bundle / "taxonomy" / "manifest.json", "taxonomy manifest"
    )
    variants = _validate_prepared_variants(bundle)
    top_taxonomy = top.get("taxonomy")
    if (
        not isinstance(top_taxonomy, dict)
        or top_taxonomy.get("path") != "taxonomy/taxcat.dmp"
        or top_taxonomy.get("identity") != taxonomy_identity
        or top_taxonomy.get("manifest_identity")
        != file_identity(bundle / "taxonomy" / "manifest.json")
        or top_taxonomy.get("source") != taxonomy_manifest.get("source")
    ):
        raise CacheError("References bundle taxonomy provenance is inconsistent")
    top_variants = top.get("variants")
    if (
        not isinstance(top_variants, dict)
        or top_variants.get("path") != "variants"
        or top_variants.get("imported") != bool(variants.get("imported"))
        or top_variants.get("manifest_identity")
        != file_identity(bundle / "variants" / "manifest.json")
        or top_variants.get("clinvar")
        != file_identity(bundle / "variants" / VARIANT_FILES["clinvar"])
        or top_variants.get("gwas")
        != file_identity(bundle / "variants" / VARIANT_FILES["gwas"])
    ):
        raise CacheError("References bundle variant provenance is inconsistent")
    generated = not bool(variants.get("imported"))
    source_manifest = bundle / "sources" / "manifest.json"
    external_manifest = bundle / "sources" / "external.json"
    if generated:
        if source_manifest.is_file():
            source_ids = _validate_raw_source_manifest(source_manifest, verify_sources)
        elif external_manifest.is_file():
            external = read_json(external_manifest, "external source descriptor")
            if external.get("schema_version") not in {1, 3} or not isinstance(
                external.get("sources"), dict
            ):
                raise CacheError(
                    f"Invalid external source descriptor: {external_manifest}"
                )
            missing = sorted(set(VARIANT_SOURCE_NAMES) - set(external["sources"]))
            if missing:
                raise CacheError(
                    "External source descriptor is missing required source(s): "
                    + ", ".join(missing)
                )
            source_ids = {
                "schema_version": 3,
                "sources": {
                    name: _source_entry_identity(external["sources"][name], name)
                    for name in VARIANT_SOURCE_NAMES
                },
            }
            if verify_sources:
                source_ref = external.get("source_manifest")
                if not isinstance(source_ref, dict) or not isinstance(
                    source_ref.get("path"), str
                ):
                    raise CacheError(
                        f"External source descriptor lacks its manifest locator: {external_manifest}"
                    )
                raw_manifest_path = (
                    Path(source_ref["path"]).expanduser().resolve(strict=True)
                )
                if file_identity(raw_manifest_path) != _identity_record(
                    source_ref, "external source manifest"
                ):
                    raise CacheError(
                        f"External raw source manifest changed: {raw_manifest_path}"
                    )
                actual_sources = source_manifest_identities(raw_manifest_path)
                expected_sources = {
                    "schema_version": 3,
                    "sources": {
                        name: _source_entry_identity(external["sources"][name], name)
                        for name in VARIANT_SOURCE_NAMES
                    },
                }
                if actual_sources != expected_sources:
                    raise CacheError(
                        "External source descriptor differs from its raw source manifest"
                    )
        else:
            raise CacheError(
                "Generated variant outputs require sources/manifest.json or sources/external.json"
            )
        variant_sources = variants.get("sources")
        if not isinstance(variant_sources, dict):
            raise CacheError("Generated variant manifest is missing source identities")
        missing_sources = sorted(set(VARIANT_SOURCE_NAMES) - set(variant_sources))
        if missing_sources:
            raise CacheError(
                "Generated variant manifest is missing source identities: "
                + ", ".join(missing_sources)
            )
        for name in VARIANT_SOURCE_NAMES:
            expected_content = _source_content_identity(
                source_ids["sources"][name], name
            )
            actual_content = _source_content_identity(variant_sources[name], name)
            if actual_content != expected_content:
                raise CacheError(
                    f"Prepared variant source identity differs from raw source {name}"
                )
    top_sources = top.get("sources")
    if top_sources is None:
        if generated:
            raise CacheError("Generated references bundle is missing source provenance")
    elif not isinstance(top_sources, dict):
        raise CacheError("References bundle source provenance is invalid")
    else:
        source_path = top_sources.get("manifest")
        if not isinstance(source_path, str):
            raise CacheError("References bundle source provenance has no manifest path")
        relative = Path(source_path)
        if relative.is_absolute() or ".." in relative.parts:
            raise CacheError("References bundle source manifest path is unsafe")
        source_file = bundle / relative
        if not source_file.is_file() or file_identity(source_file) != top_sources.get(
            "manifest_identity"
        ):
            raise CacheError(
                "References bundle source manifest differs from its record"
            )
    return top


def implementation_identities(repo: Path, paths: tuple[Path, ...]) -> dict[str, Any]:
    """Record reference-stage code identities without including run-specific values."""
    identities = {}
    missing = []
    for path in paths:
        if path.is_file():
            try:
                label = path.resolve().relative_to(repo.resolve()).as_posix()
            except ValueError:
                label = path.name
            identities[label] = file_identity(path)
        else:
            missing.append(path)
    if missing:
        raise CacheError(
            "Reference implementation is incomplete; missing file(s): "
            + ", ".join(str(path) for path in missing)
        )
    return identities
