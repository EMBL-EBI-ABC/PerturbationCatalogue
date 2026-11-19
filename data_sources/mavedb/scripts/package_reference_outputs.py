#!/usr/bin/env python3
"""Assemble and checkpoint reference components for the MaveDB pipeline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from data_sources.mavedb.analysis._reference_io import verify_file
from data_sources.mavedb.analysis._reference_sources import (
    DEFAULT_ID_MAP,
    SOURCE_SCHEMA_VERSION,
    SOURCE_URLS,
    load_sources,
)

CLINVAR_FILE = "mave_clinvar_variants.tsv"
GWAS_FILE = "mave_gwas_variants.tsv"
CLINVAR_REQUIRED = {
    "uniprotAccessions",
    "clinicalSignificances",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
}
GWAS_REQUIRED = {
    "studyId",
    "beta",
    "pValueMantissa",
    "pValueExponent",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
}
REQUIRED_SOURCES = ("coding_variants", "evidence", "credible_sets", "id_map")
PARQUET_SOURCES = {"coding_variants", "evidence", "credible_sets"}


def file_identity(path: Path) -> dict[str, int | str]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _tsv_summary(path: Path, required: set[str]) -> dict[str, int | str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Reference TSV is missing or empty: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) != len(set(header)):
            raise ValueError(f"Reference TSV has a missing or duplicate header: {path}")
        missing = sorted(required - set(header))
        if missing:
            raise ValueError(
                f"Reference TSV {path} is missing required columns: {', '.join(missing)}"
            )
        rows = sum(1 for _ in reader)
    return {"rows": rows, **file_identity(path)}


def _validate_taxonomy(directory: Path) -> dict[str, Any]:
    taxonomy = directory / "taxcat.dmp"
    manifest_path = directory / "manifest.json"
    manifest = _json(manifest_path)
    if not taxonomy.is_file() or not taxonomy.stat().st_size:
        raise ValueError(f"Taxonomy bundle is missing taxcat.dmp: {directory}")
    expected = manifest.get("sha256")
    actual = file_identity(taxonomy)
    if expected != actual["sha256"] or manifest.get("bytes") != actual["bytes"]:
        raise ValueError(f"Taxonomy checksum differs from its manifest: {taxonomy}")
    return manifest


def _validate_source_manifest(path: Path) -> dict[str, Any]:
    path = Path(path).expanduser().resolve(strict=True)
    manifest = _json(path)
    if manifest.get("schema_version") != SOURCE_SCHEMA_VERSION:
        raise ValueError(f"Unsupported source manifest: {path}")
    raw_sources = manifest.get("sources")
    if not isinstance(raw_sources, dict) or not set(REQUIRED_SOURCES).issubset(
        raw_sources
    ):
        missing = sorted(set(REQUIRED_SOURCES) - set(raw_sources or {}))
        raise ValueError(
            f"Source manifest is missing selected products: {', '.join(missing)}"
        )
    resolved_sources = load_sources(path)["sources"]
    for name in REQUIRED_SOURCES:
        record = resolved_sources[name]
        if name in PARQUET_SOURCES:
            for descriptor in record["files"]:
                verify_file(descriptor)
        else:
            verify_file(record)
    return {
        **{
            key: value
            for key, value in manifest.items()
            if key not in {"sources", "homology"}
        },
        "sources": {name: raw_sources[name] for name in REQUIRED_SOURCES},
    }


def _validate_variants(directory: Path) -> dict[str, Any]:
    manifest = _json(directory / "manifest.json")
    if manifest.get("schema_version") != 1:
        raise ValueError(
            f"Unsupported variants manifest: {directory / 'manifest.json'}"
        )
    results = manifest.get("results")
    if not isinstance(results, dict):
        raise ValueError(f"Variants manifest is missing results: {directory}")
    for name, filename, required in (
        ("clinvar", CLINVAR_FILE, CLINVAR_REQUIRED),
        ("gwas", GWAS_FILE, GWAS_REQUIRED),
    ):
        result = results.get(name)
        if not isinstance(result, dict) or result.get("file") != filename:
            raise ValueError(
                f"Variants manifest is missing the {name} table: {directory}"
            )
        identity = {
            "bytes": result.get("bytes"),
            "sha256": result.get("sha256"),
        }
        actual = file_identity(directory / filename)
        if actual != identity:
            raise ValueError(
                f"Variant {name} table differs from its manifest: {directory}"
            )
        summary = _tsv_summary(directory / filename, required)
        if result.get("rows") != summary["rows"]:
            raise ValueError(f"Variant {name} row count differs from its manifest")
    return manifest


def _source_file(root: Path, relative: str, label: str) -> Path:
    path = Path(relative)
    if not relative or path.is_absolute():
        raise ValueError(f"Source {label} must have a relative path")
    resolved = (root / path).resolve(strict=True)
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"Source {label} escapes its source bundle")
    if not resolved.is_file():
        raise ValueError(f"Source {label} is not a regular file: {resolved}")
    return resolved


def _materialize_source_bundle(source_manifest: Path, destination: Path) -> Path:
    """Copy only the selected, verified source files and normalize the manifest name."""
    source_manifest = source_manifest.expanduser().resolve(strict=True)
    manifest = _validate_source_manifest(source_manifest)
    source_root = source_manifest.parent
    destination.mkdir(parents=True, exist_ok=False)
    for name in REQUIRED_SOURCES:
        record = manifest["sources"][name]
        if name in PARQUET_SOURCES:
            descriptors = record["files"]
        else:
            descriptors = [record]
        for index, descriptor in enumerate(descriptors):
            relative = descriptor["path"]
            source = _source_file(source_root, relative, f"{name}[{index}]")
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)

    output_manifest = destination / "manifest.json"
    _write_json(output_manifest, manifest)
    _validate_source_manifest(output_manifest)
    return output_manifest


def _rewrite_variant_source_manifest(
    directory: Path, descriptor: dict[str, int | str]
) -> None:
    path = directory / "manifest.json"
    manifest = _json(path)
    if manifest.get("imported"):
        return
    manifest["source_manifest"] = descriptor
    manifest["source_manifest_schema_version"] = SOURCE_SCHEMA_VERSION
    _write_json(path, manifest)


def import_taxonomy(source: Path, source_label: str, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(source, output / "taxcat.dmp")
    _write_json(
        output / "manifest.json",
        {
            "schema_version": 1,
            "source": source_label,
            "path": "taxcat.dmp",
            **file_identity(output / "taxcat.dmp"),
        },
    )


def reuse_taxonomy(source: Path, output: Path) -> None:
    if source.is_dir():
        shutil.copytree(source, output)
        _validate_taxonomy(output)
        return
    import_taxonomy(source, str(source), output)


def import_variants(clinvar: Path, gwas: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(clinvar, output / CLINVAR_FILE)
    shutil.copyfile(gwas, output / GWAS_FILE)
    results = {
        "clinvar": {
            "file": CLINVAR_FILE,
            **_tsv_summary(output / CLINVAR_FILE, CLINVAR_REQUIRED),
        },
        "gwas": {
            "file": GWAS_FILE,
            **_tsv_summary(output / GWAS_FILE, GWAS_REQUIRED),
        },
    }
    _write_json(
        output / "manifest.json",
        {
            "schema_version": 1,
            "imported": True,
            "source": "paired validation TSV inputs",
            "inputs": {
                "clinvar": {"path": str(clinvar.resolve()), **file_identity(clinvar)},
                "gwas": {"path": str(gwas.resolve()), **file_identity(gwas)},
            },
            "stages": ["clinvar", "gwas"],
            "results": results,
        },
    )


def _copy_component(source: Path, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)
    for name in ("pipeline_manifest.json", ".pipeline-result.json"):
        (destination / name).unlink(missing_ok=True)


def assemble(
    taxonomy: Path,
    variants: Path,
    output: Path,
    source_mode: str,
    source_manifest: Path | None,
    replaceable_roots: Sequence[Path] = (),
    source_manifest_locator: Path | None = None,
) -> None:
    output = Path(output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".references-", dir=output.parent) as temp:
        staged = Path(temp) / "references"
        staged.mkdir()
        taxonomy_manifest = _validate_taxonomy(taxonomy)
        _copy_component(taxonomy, staged / "taxonomy")
        variant_manifest = _validate_variants(variants)
        _copy_component(variants, staged / "variants")

        source_record: dict[str, Any] | None = None
        local_source_manifest: Path | None = None
        external_source_manifest: Path | None = None
        if source_mode in {"generated", "retained"}:
            if source_manifest is None:
                raise ValueError(
                    f"{source_mode} source mode requires a source manifest"
                )
            source_manifest = source_manifest.expanduser().resolve(strict=True)
            local_source_manifest = _materialize_source_bundle(
                source_manifest, staged / "sources"
            )
            source_record = {
                "mode": source_mode,
                "manifest": "sources/manifest.json",
                "manifest_identity": file_identity(staged / "sources/manifest.json"),
            }
        elif source_mode == "external":
            if source_manifest is None:
                raise ValueError("External source mode requires a source manifest")
            source_manifest = source_manifest.expanduser().resolve(strict=True)
            source_data = _validate_source_manifest(source_manifest)
            external_source_manifest = (
                source_manifest_locator.expanduser().resolve(strict=True)
                if source_manifest_locator is not None
                else source_manifest
            )
            if file_identity(external_source_manifest) != file_identity(
                source_manifest
            ):
                raise ValueError(
                    "Staged external source manifest differs from its original locator"
                )
            original_data = _validate_source_manifest(external_source_manifest)
            if original_data["sources"] != source_data["sources"]:
                raise ValueError(
                    "Staged external sources differ from the original source manifest"
                )
            source_root = external_source_manifest.parent.resolve()
            candidates = [external_source_manifest]
            for name in REQUIRED_SOURCES:
                record = source_data["sources"][name]
                candidates.append(source_root / record["path"])
                candidates.extend(
                    source_root / descriptor["path"]
                    for descriptor in record.get("files", [])
                )
            replaceable = [output, *replaceable_roots]
            must_materialize = any(
                candidate.resolve().is_relative_to(root.expanduser().resolve())
                for root in replaceable
                for candidate in candidates
            )
            if must_materialize:
                local_source_manifest = _materialize_source_bundle(
                    source_manifest, staged / "sources"
                )
                source_record = {
                    "mode": "external-materialized",
                    "manifest": "sources/manifest.json",
                    "manifest_identity": file_identity(
                        staged / "sources/manifest.json"
                    ),
                    "origin_source_manifest": {
                        "path": str(external_source_manifest),
                        **file_identity(external_source_manifest),
                    },
                }
            else:
                sources_dir = staged / "sources"
                sources_dir.mkdir()
                descriptor = {
                    "schema_version": 1,
                    "source_manifest": {
                        "path": str(external_source_manifest),
                        **file_identity(external_source_manifest),
                    },
                    "source_root": str(source_root),
                    "raw_schema_version": SOURCE_SCHEMA_VERSION,
                    "sources": source_data["sources"],
                }
                _write_json(sources_dir / "external.json", descriptor)
                external_source_manifest = Path(descriptor["source_manifest"]["path"])
                source_record = {
                    "mode": source_mode,
                    "manifest": "sources/external.json",
                    "manifest_identity": file_identity(sources_dir / "external.json"),
                }
        elif source_mode != "imported":
            raise ValueError(f"Unknown source mode: {source_mode}")

        if not variant_manifest.get("imported"):
            if local_source_manifest is not None:
                source_descriptor = {
                    "path": "../sources/manifest.json",
                    **file_identity(local_source_manifest),
                }
            elif external_source_manifest is not None:
                source_descriptor = {
                    "path": str(external_source_manifest),
                    **file_identity(external_source_manifest),
                }
            else:
                raise ValueError("Generated variants require a source manifest")
            _rewrite_variant_source_manifest(staged / "variants", source_descriptor)
            _validate_variants(staged / "variants")

        manifest = {
            "schema_version": 1,
            "taxonomy": {
                "path": "taxonomy/taxcat.dmp",
                "source": taxonomy_manifest.get("source"),
                "identity": file_identity(staged / "taxonomy/taxcat.dmp"),
                "manifest_identity": file_identity(staged / "taxonomy/manifest.json"),
            },
            "sources": source_record,
            "variants": {
                "path": "variants",
                "imported": bool(variant_manifest.get("imported")),
                "manifest_identity": file_identity(staged / "variants/manifest.json"),
                "clinvar": file_identity(staged / f"variants/{CLINVAR_FILE}"),
                "gwas": file_identity(staged / f"variants/{GWAS_FILE}"),
            },
        }
        _write_json(staged / "manifest.json", manifest)
        backup = Path(temp) / "previous"
        if output.exists():
            if not output.is_dir():
                raise FileExistsError(f"Reference output is not a directory: {output}")
            output.rename(backup)
        try:
            staged.rename(output)
        except (OSError, KeyboardInterrupt, SystemExit):
            if output.exists():
                shutil.rmtree(output)
            if backup.exists():
                backup.rename(output)
            raise


def _implementation_subset(spec: dict[str, Any], kind: str) -> dict[str, Any]:
    identities = spec.get("identities", {})
    implementation = identities.get("implementation", {})
    if not isinstance(implementation, dict):
        return {}
    marker_sets = {
        "taxonomy": (
            "download_taxonomy.py",
            "package_reference_outputs.py",
            "pipeline/modules/taxonomy.nf",
        ),
        "sources": (
            "download_reference_data.py",
            "analysis/_reference_downloads.py",
            "analysis/_reference_sources.py",
            "analysis/_reference_io.py",
            "pipeline/modules/references.nf",
            "pipeline/modules/references_download.nf",
        ),
        "variant_processing": (
            "process_variant_data.py",
            "analysis/_reference_sources.py",
            "analysis/_reference_io.py",
            "analysis/_variant_data.py",
            "pipeline/modules/references.nf",
            "pipeline/modules/references_prepare.nf",
        ),
        "variant_generation": (
            "download_reference_data.py",
            "analysis/_reference_downloads.py",
            "analysis/_reference_sources.py",
            "analysis/_reference_io.py",
            "process_variant_data.py",
            "analysis/_variant_data.py",
            "pipeline/modules/references.nf",
            "pipeline/modules/references_download.nf",
            "pipeline/modules/references_prepare.nf",
        ),
        "variant_import": (
            "package_reference_outputs.py",
            "pipeline/modules/references.nf",
            "pipeline/modules/references_prepare.nf",
        ),
    }
    if kind == "variants":
        if identities.get("variants_import") is not None or spec.get("variants_reuse"):
            markers = marker_sets["variant_import"]
        elif identities.get("source_manifest") is not None:
            markers = marker_sets["variant_processing"]
        else:
            markers = marker_sets["variant_generation"]
    elif kind in {"taxonomy", "sources"}:
        markers = marker_sets[kind]
    else:
        raise ValueError(f"Unsupported checkpoint component: {kind}")
    return {
        name: value
        for name, value in implementation.items()
        if any(name.endswith(marker) for marker in markers)
    }


def _variant_reuse_identity(spec: dict[str, Any]) -> dict[str, Any] | None:
    value = spec.get("variants_reuse")
    if not value:
        return None
    directory = Path(value).expanduser().resolve(strict=True)
    if directory.name != "variants" and (directory / "variants").is_dir():
        directory = directory / "variants"
    return {
        "manifest": file_identity(directory / "manifest.json"),
        "clinvar": file_identity(directory / CLINVAR_FILE),
        "gwas": file_identity(directory / GWAS_FILE),
    }


def _taxonomy_reuse_identity(spec: dict[str, Any]) -> dict[str, Any] | None:
    value = spec.get("taxonomy_reuse")
    if not value:
        return None
    path = Path(value).expanduser().resolve(strict=True)
    if path.is_file():
        return {"taxonomy": file_identity(path)}
    _validate_taxonomy(path)
    return {
        "taxonomy": file_identity(path / "taxcat.dmp"),
        "manifest": file_identity(path / "manifest.json"),
    }


def _source_checkpoint_inputs(
    request: dict[str, Any], spec: dict[str, Any]
) -> dict[str, Any]:
    refresh = bool(request.get("refresh_references"))
    inputs: dict[str, Any] = {
        "urls": {
            name: request.get(f"{name}_url") or SOURCE_URLS[name]
            for name in ("coding_variants", "evidence", "credible_sets")
        },
        "id_map": file_identity(DEFAULT_ID_MAP),
        "refresh_references": refresh,
    }
    if refresh:
        inputs["attempt"] = (
            request.get("references_attempt")
            or Path(spec.get("references_component_cache", "")).name
        )
    return inputs


def _checkpoint_identity(spec_path: Path, kind: str) -> str:
    spec = _json(spec_path)
    identities = spec.get("identities", {})
    if not isinstance(identities, dict):
        raise ValueError(f"Invalid checkpoint identities in {spec_path}")
    request_path = spec.get("request")
    request = _json(Path(request_path)) if request_path else {}

    if kind == "taxonomy":
        inputs: dict[str, Any] = {
            "taxonomy_import": identities.get("taxonomy_import"),
            "taxonomy_url": request.get("taxonomy_url"),
            "taxonomy_reuse": _taxonomy_reuse_identity(spec),
        }
    elif kind == "sources":
        inputs = _source_checkpoint_inputs(request, spec)
    elif kind == "variants":
        inputs = {
            "source_manifest": identities.get("source_manifest"),
            "variants_import": identities.get("variants_import"),
            "variants_reuse": _variant_reuse_identity(spec),
        }
        if not any(
            inputs[name] is not None
            for name in ("source_manifest", "variants_import", "variants_reuse")
        ):
            inputs["generated_sources"] = _source_checkpoint_inputs(request, spec)
    else:
        raise ValueError(f"Unsupported checkpoint component: {kind}")

    identity = json.dumps(
        {
            "kind": kind,
            "settings": spec.get("settings"),
            "inputs": inputs,
            "implementation": _implementation_subset(spec, kind),
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(identity.encode("utf-8")).hexdigest()


def _validate_checkpoint(path: Path, kind: str) -> None:
    if not path.is_dir():
        raise ValueError(f"Missing {kind} checkpoint directory: {path}")
    if kind == "taxonomy":
        _validate_taxonomy(path)
    elif kind == "sources":
        _validate_source_manifest(path / "manifest.json")
    elif kind == "variants":
        _validate_variants(path)
    else:
        raise ValueError(f"Unsupported checkpoint component: {kind}")


def checkpoint(action: str, kind: str, root: Path, spec: Path, directory: Path) -> bool:
    """Restore/save one resumable component under the configured attempt root."""
    identity = _checkpoint_identity(spec, kind)
    checkpoint_root = root / "component_checkpoints"
    checkpoint_root.mkdir(parents=True, exist_ok=True)
    data_path = root / kind
    marker = checkpoint_root / f"{kind}.json"
    if action == "restore":
        if not data_path.is_dir() or not marker.is_file():
            return False
        record = _json(marker)
        if record.get("identity") != identity:
            return False
        _validate_checkpoint(data_path, kind)
        if directory.exists():
            shutil.rmtree(directory)
        shutil.copytree(data_path, directory)
        return True
    if action != "save":
        raise ValueError(f"Unknown checkpoint action: {action}")
    _validate_checkpoint(directory, kind)
    temporary = Path(tempfile.mkdtemp(prefix=f".{kind}-", dir=root))
    staged = temporary / kind
    backup = temporary / "previous"
    try:
        shutil.copytree(directory, staged)
        _validate_checkpoint(staged, kind)
        if data_path.exists():
            data_path.rename(backup)
        try:
            staged.rename(data_path)
            _write_json(marker, {"identity": identity, "kind": kind})
        except OSError:
            if data_path.exists():
                shutil.rmtree(data_path)
            if backup.exists():
                backup.rename(data_path)
            raise
    finally:
        shutil.rmtree(temporary, ignore_errors=True)
    return True


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    taxonomy_parser = subparsers.add_parser("import-taxonomy")
    taxonomy_parser.add_argument("--source", type=Path, required=True)
    taxonomy_parser.add_argument("--source-label", required=True)
    taxonomy_parser.add_argument("--out", type=Path, required=True)

    reuse_parser = subparsers.add_parser("reuse-taxonomy")
    reuse_parser.add_argument("--source", type=Path, required=True)
    reuse_parser.add_argument("--out", type=Path, required=True)

    variants_parser = subparsers.add_parser("import-variants")
    variants_parser.add_argument("--clinvar", type=Path, required=True)
    variants_parser.add_argument("--gwas", type=Path, required=True)
    variants_parser.add_argument("--out", type=Path, required=True)

    assemble_parser = subparsers.add_parser("assemble")
    assemble_parser.add_argument("--taxonomy", type=Path, required=True)
    assemble_parser.add_argument("--variants", type=Path, required=True)
    assemble_parser.add_argument(
        "--source-mode",
        choices=("generated", "retained", "external", "imported"),
        required=True,
    )
    assemble_parser.add_argument("--source-manifest", type=Path)
    assemble_parser.add_argument("--source-manifest-locator", type=Path)
    assemble_parser.add_argument(
        "--replaceable-root", type=Path, action="append", default=[]
    )
    assemble_parser.add_argument("--out", type=Path, required=True)

    checkpoint_parser = subparsers.add_parser("checkpoint")
    checkpoint_parser.add_argument(
        "--action", choices=("restore", "save"), required=True
    )
    checkpoint_parser.add_argument(
        "--kind", choices=("taxonomy", "sources", "variants"), required=True
    )
    checkpoint_parser.add_argument("--root", type=Path, required=True)
    checkpoint_parser.add_argument("--spec", type=Path, required=True)
    checkpoint_parser.add_argument("--directory", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "import-taxonomy":
        import_taxonomy(args.source, args.source_label, args.out)
    elif args.command == "reuse-taxonomy":
        reuse_taxonomy(args.source, args.out)
    elif args.command == "import-variants":
        import_variants(args.clinvar, args.gwas, args.out)
    elif args.command == "assemble":
        assemble(
            args.taxonomy,
            args.variants,
            args.out,
            args.source_mode,
            args.source_manifest,
            args.replaceable_root,
            args.source_manifest_locator,
        )
    else:
        restored = checkpoint(
            args.action, args.kind, args.root, args.spec, args.directory
        )
        if args.action == "restore" and not restored:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
