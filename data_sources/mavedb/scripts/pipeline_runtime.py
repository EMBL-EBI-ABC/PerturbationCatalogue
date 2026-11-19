#!/usr/bin/env python3
"""Durable stage policy used by the MaveDB Nextflow workflow.

Nextflow schedules computation. This module resolves portable completed results,
checkpoints each successful task, and publishes only a consistent result set.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import shutil
import sys
import tempfile
import uuid
from typing import Any, Sequence

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from data_sources.mavedb.scripts._pipeline_cache import (
    CacheError,
    ResultStore,
    file_identity,
)
from data_sources.mavedb.scripts._pipeline_references import (
    VARIANT_FILES,
    implementation_identities,
    read_json as read_reference_json,
    source_manifest_identities,
    taxonomy_manifest_identity,
    validate_reference_bundle,
    validate_tsv,
    validate_variant_component,
)
from data_sources.mavedb.scripts.download_taxonomy import (
    DEFAULT_URL as DEFAULT_TAXONOMY_URL,
)

LOGGER = logging.getLogger(__name__)
SCRIPTS = Path(__file__).resolve().parent
STAGES = {
    "ingest",
    "references",
    "process",
    "tables",
    "gmm",
    "merge",
    "validation",
    "qc",
    "reports",
}
PUBLIC = {"gmm-plan": "gmm", "gmm-batch": "gmm"}

# Keep these headers aligned with _validation_inputs.py. The policy checks the
# small headers with the standard library before Nextflow can schedule fitting.
VALIDATION_CLINVAR_REQUIRED = {
    "uniprotAccessions",
    "clinicalSignificances",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
}
VALIDATION_GWAS_REQUIRED = {
    "studyId",
    "beta",
    "pValueMantissa",
    "pValueExponent",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
}
VALIDATION_IMPLEMENTATION = (
    Path(__file__).resolve(),
    Path(__file__).resolve().parent / "_pipeline_cache.py",
    REPO / "data_sources/mavedb/pipeline/modules/validation.nf",
    REPO / "data_sources/mavedb/analysis/mave_gmm_validation.py",
    REPO / "data_sources/mavedb/analysis/_validation_inputs.py",
    REPO / "data_sources/mavedb/analysis/_validation_core.py",
    REPO / "data_sources/mavedb/analysis/_validation_plots.py",
    REPO / "data_sources/mavedb/analysis/_homology_inputs.py",
    REPO / "data_sources/mavedb/analysis/_homology_clinvar.py",
    REPO / "data_sources/mavedb/analysis/_reference_io.py",
    REPO / "data_sources/mavedb/analysis/mave_summary/data.py",
    REPO / "data_sources/mavedb/analysis/_mave_qc.py",
    REPO / "data_sources/mavedb/scripts/process_mavedb.py",
)
REFERENCES_IMPLEMENTATION = (
    Path(__file__).resolve(),
    Path(__file__).resolve().parent / "_pipeline_cache.py",
    Path(__file__).resolve().parent / "_pipeline_references.py",
    SCRIPTS / "download_taxonomy.py",
    SCRIPTS / "download_reference_data.py",
    SCRIPTS / "process_variant_data.py",
    REPO / "data_sources/mavedb/analysis/_reference_sources.py",
    REPO / "data_sources/mavedb/analysis/_reference_downloads.py",
    REPO / "data_sources/mavedb/analysis/_reference_io.py",
    REPO / "data_sources/mavedb/analysis/_variant_data.py",
    REPO / "data_sources/mavedb/analysis/_variant_annotations.py",
    REPO / "data_sources/mavedb/pipeline/main.nf",
    REPO / "data_sources/mavedb/pipeline/modules/references.nf",
    REPO / "data_sources/mavedb/pipeline/modules/references_download.nf",
    REPO / "data_sources/mavedb/pipeline/modules/references_prepare.nf",
    REPO / "data_sources/mavedb/pipeline/modules/taxonomy.nf",
    SCRIPTS / "package_reference_outputs.py",
)


def stage_error(
    stage: str,
    store: ResultStore,
    error: Exception,
    requested_inputs: dict[str, Any] | None = None,
) -> CacheError:
    descendants = {
        "ingest": [
            "process",
            "tables",
            "gmm",
            "merge",
            "validation",
            "qc",
            "reports",
        ],
        "taxonomy": [
            "process",
            "tables",
            "gmm",
            "merge",
            "validation",
            "qc",
            "reports",
        ],
        "references": [],
        "process": ["tables", "gmm", "merge", "validation", "qc", "reports"],
        "tables": ["qc"],
        "gmm": ["merge", "validation", "reports"],
        "merge": ["validation", "reports"],
    }
    descendants_for_error = descendants.get(stage, [])
    if stage == "references":
        record = store.record("references")
        changed = set()
        if record and requested_inputs is not None:
            previous = record.get("inputs", {})
            changed = {
                name
                for name in set(previous) | set(requested_inputs)
                if previous.get(name) != requested_inputs.get(name)
            }
        if not changed or "implementation" in changed or "taxonomy_import" in changed:
            if changed:
                descendants_for_error = [
                    "process",
                    "tables",
                    "gmm",
                    "merge",
                    "validation",
                    "qc",
                    "reports",
                ]
            elif record is None:
                descendants_for_error = []
            else:
                outputs = record.get("outputs", {})
                bundle = store.root / "references"
                try:
                    taxonomy_changed = outputs.get(
                        "taxonomy/taxcat.dmp"
                    ) != file_identity(bundle / "taxonomy/taxcat.dmp")
                except OSError:
                    taxonomy_changed = True
                variants_changed = False
                for filename in VARIANT_FILES.values():
                    name = f"variants/{filename}"
                    try:
                        if outputs.get(name) != file_identity(bundle / name):
                            variants_changed = True
                    except OSError:
                        variants_changed = True
                if taxonomy_changed:
                    descendants_for_error = [
                        "process",
                        "tables",
                        "gmm",
                        "merge",
                        "validation",
                        "qc",
                        "reports",
                    ]
                elif variants_changed:
                    descendants_for_error = ["validation"]
                else:
                    descendants_for_error = []
        elif changed & {"variants_import", "source_manifest"}:
            descendants_for_error = ["validation"]
    affected = [stage]
    for name in descendants_for_error:
        if store.record(name) is not None:
            affected.append(name)
    command = f"--rerun {','.join(affected)}"
    if stage == "references" and "source settings differ" in str(error).lower():
        command = "--rerun references --refresh_references"
    return CacheError(
        f"{error}. Affected stages: {', '.join(affected)}. Continue with {command}"
    )


def provenance_notice(stage: str, path: Path) -> None:
    if stage == "gmm":
        from data_sources.mavedb.analysis._gmm_pipeline import implementation_identity

        saved = read_json(path / "batch_plan.json").get("implementation")
        if saved is not None and saved != implementation_identity():
            LOGGER.warning(
                "REUSE gmm with its saved implementation provenance; use --rerun gmm to apply current fitting code"
            )
    elif stage == "process" and (path / "processing_manifest.json").is_file():
        saved = read_json(path / "processing_manifest.json")["implementation"][
            "processor_sha256"
        ]
        if saved != file_identity(SCRIPTS / "process_mavedb.py")["sha256"]:
            LOGGER.warning(
                "REUSE process with its saved implementation provenance; use --rerun process to apply current processing code"
            )


def read_json(path: Path | str) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def identity(path: Path) -> dict[str, Any]:
    if path.is_file():
        return file_identity(path)
    if not path.is_dir():
        raise CacheError(f"Required input does not exist: {path}")
    return {
        p.relative_to(path).as_posix(): file_identity(p)
        for p in sorted(path.rglob("*"))
        if p.is_file()
        and p.name not in {".pipeline-result.json", "pipeline_manifest.json"}
    }


def gmm_arguments(p: dict[str, Any]) -> list[str]:
    options = {
        "profile": "profile",
        "model": "model",
        "variants_per_batch": "variants_per_batch",
        "seed": "seed",
        "chains": "chains",
        "parallel_chains": "parallel-chains",
        "iter_warmup": "iter-warmup",
        "iter_sampling": "iter-sampling",
        "max_studies": "max-studies",
        "max_study_rows": "max-study-rows",
        "max_fit_rows": "max-fit-rows",
        "memory_budget_mb": "memory-budget-mb",
    }
    args = []
    for key, option in options.items():
        value = p.get(f"gmm_{key}")
        if value is not None:
            args.extend([f"--{option}", str(value)])
    studies = p.get("gmm_study_ids", [])
    if isinstance(studies, str):
        studies = [s.strip() for s in studies.split(",") if s.strip()]
    for study in studies:
        args.extend(["--study-id", study])
    return args


def gmm_configuration(p: dict[str, Any]) -> dict[str, Any]:
    # Resolve profile defaults without importing CmdStan or reading the inputs.
    from data_sources.mavedb.scripts.fit_gmm_batch import build_parser, _build_config
    from data_sources.mavedb.analysis._gmm_pipeline import configuration

    argv = gmm_arguments(p)
    parser = build_parser()
    return configuration(_build_config(parser, parser.parse_args(argv), argv))


def settings(
    stage: str, p: dict[str, Any], store: ResultStore | None = None
) -> dict[str, Any]:
    if stage == "ingest":
        return {
            "source": p.get("archive") or p["zenodo_url"],
            "expected_sha256": p.get("archive_sha256"),
        }
    if stage == "references":
        requested: dict[str, Any] = {}
        taxonomy_url = p.get("taxonomy_url")
        if (
            not p.get("taxonomy")
            and taxonomy_url
            and taxonomy_url != DEFAULT_TAXONOMY_URL
        ):
            requested["taxonomy_url"] = taxonomy_url
        for name in ("coding_variants_url", "evidence_url", "credible_sets_url"):
            if p.get(name):
                requested[name] = p[name]
        if store is not None:
            try:
                previous = store.record("references")
            except CacheError:
                if "references" in p.get("selected", []) or store.pending("references"):
                    previous = None
                else:
                    raise
            saved = previous.get("settings", {}) if previous else {}
            if (
                not p.get("refresh_references")
                and not p.get("taxonomy")
                and taxonomy_url == DEFAULT_TAXONOMY_URL
            ):
                saved_taxonomy_url = saved.get("taxonomy_url")
                if saved_taxonomy_url:
                    requested.setdefault("taxonomy_url", saved_taxonomy_url)
            if not p.get("refresh_references"):
                for name in (
                    "coding_variants_url",
                    "evidence_url",
                    "credible_sets_url",
                ):
                    if not p.get(name) and name in saved:
                        requested[name] = saved[name]
        return requested
    if stage == "taxonomy":
        return {"source": p.get("taxonomy") or p["taxonomy_url"]}
    if stage in {"gmm", "gmm-plan", "gmm-batch"}:
        return gmm_configuration(p) | {
            "fail_on_study_error": p["gmm_fail_on_study_error"]
        }
    if stage == "validation":
        no_plots = bool(p.get("validation_no_plots", False))
        formats = p.get("validation_formats", ["png", "pdf"])
        return {
            "formats": [] if no_plots else list(dict.fromkeys(formats)),
            "dpi": int(p.get("validation_dpi", 150)),
            "no_plots": no_plots,
        }
    return {}


def dependencies(stage: str, paths: dict[str, str]) -> dict[str, Any]:
    found: dict[str, Any] = {}

    def add(name: str, path: Path) -> None:
        found[name] = identity(path)

    if stage in {"process", "tables"} and "ingest" in paths:
        add("raw", Path(paths["ingest"]) / "mavedb_dump")
    if stage in {"process", "tables"}:
        if "references" in paths:
            add("taxonomy", Path(paths["references"]) / "taxonomy/taxcat.dmp")
        elif "taxonomy" in paths:
            # Compatibility for older task invocations and saved workflow requests.
            add("taxonomy", Path(paths["taxonomy"]) / "taxcat.dmp")
    if stage in {"tables", "gmm-plan", "gmm", "gmm-batch", "qc", "reports"}:
        for name in ("scores.tsv", "metadata.tsv"):
            add(name, Path(paths["process"]) / name)
    if stage == "gmm-batch" and "gmm-plan" in paths:
        add("plan", Path(paths["gmm-plan"]) / "batch_plan.json")
    if stage == "merge":
        add("gmm", Path(paths["gmm"]))
    if stage == "qc":
        add("singles", Path(paths["tables"]) / "maves/singles.tsv")
    if stage == "reports":
        # Report metadata is separate from fitting inputs; record both explicitly.
        add("merged", Path(paths["merge"]))
    if stage == "validation":
        for name in ("scores.tsv", "metadata.tsv"):
            add(name, Path(paths["process"]) / name)
        add("merge", Path(paths["merge"]))
        for name in ("clinvar", "gwas"):
            if "references" in paths:
                reference = Path(paths["references"]) / "variants" / VARIANT_FILES[name]
            else:
                reference = Path(paths[name])
            reference, _ = validate_tsv(reference, name)
            found[name] = file_identity(reference)
        found["implementation"] = {
            path.relative_to(REPO).as_posix(): file_identity(path)
            for path in VALIDATION_IMPLEMENTATION
        }
    return found


def _validate_validation_reference(path: Path, name: str, required: set[str]) -> Path:
    """Check a local reference TSV's header and return its resolved file path."""
    resolved, _ = validate_tsv(Path(path).expanduser(), name)
    return resolved


def _prepare_reference_imports(p: dict[str, Any]) -> None:
    """Validate and normalize explicitly imported taxonomy and paired TSVs."""
    clinvar = p.get("validation_clinvar")
    gwas = p.get("validation_gwas")
    has_clinvar = isinstance(clinvar, (str, os.PathLike)) and bool(str(clinvar).strip())
    has_gwas = isinstance(gwas, (str, os.PathLike)) and bool(str(gwas).strip())
    if has_clinvar != has_gwas:
        raise CacheError(
            "--validation_clinvar and --validation_gwas must be supplied together"
        )
    if not has_clinvar:
        return
    for name, required in (
        ("clinvar", VALIDATION_CLINVAR_REQUIRED),
        ("gwas", VALIDATION_GWAS_REQUIRED),
    ):
        p[f"validation_{name}"] = str(
            _validate_validation_reference(
                Path(p[f"validation_{name}"]), name, required
            )
        )


def _prepare_validation(p: dict[str, Any]) -> None:
    """Normalize validation plotting options; references come from the stage."""
    raw_formats = p.get("validation_formats", ["png", "pdf"])
    if isinstance(raw_formats, str):
        formats = [item.strip() for item in raw_formats.split(",") if item.strip()]
    elif isinstance(raw_formats, (list, tuple)):
        formats = [str(item).strip() for item in raw_formats if str(item).strip()]
    else:
        raise CacheError("--validation_formats must be a list of png and/or pdf")
    if not formats or set(formats) - {"png", "pdf"}:
        raise CacheError("--validation_formats must contain png and/or pdf")
    p["validation_formats"] = list(dict.fromkeys(formats))
    try:
        dpi = int(p.get("validation_dpi", 150))
    except (TypeError, ValueError) as error:
        raise CacheError("--validation_dpi must be a positive integer") from error
    if dpi < 1:
        raise CacheError("--validation_dpi must be a positive integer")
    p["validation_dpi"] = dpi
    p["validation_no_plots"] = _boolean_option(
        p.get("validation_no_plots", False), "--validation_no_plots"
    )


def _has_variant_imports(p: dict[str, Any]) -> bool:
    return bool(p.get("validation_clinvar") and p.get("validation_gwas"))


def _reference_attempt_options(p: dict[str, Any]) -> dict[str, Any]:
    """Options that must stay stable while a references attempt is resumed."""
    return {
        name: p.get(name)
        for name in (
            "taxonomy",
            "taxonomy_url",
            "validation_clinvar",
            "validation_gwas",
            "reference_sources_manifest",
            "coding_variants_url",
            "evidence_url",
            "credible_sets_url",
            "references_cache_dir",
            "refresh_references",
        )
    }


def _pending_reference_request(
    store: ResultStore, attempt: str
) -> dict[str, Any] | None:
    """Find the saved request that owns a pending references attempt."""
    requests = store.root / "pipeline_requests"
    if not requests.is_dir():
        return None
    candidates = sorted(
        requests.glob("*.json"), key=lambda path: path.stat().st_mtime_ns, reverse=True
    )
    for path in candidates:
        try:
            request = read_json(path)
        except (OSError, ValueError, CacheError):
            continue
        if request.get("references_attempt") == attempt:
            return request
    return None


def _reference_inputs(
    p: dict[str, Any], store: ResultStore, verify_source_files: bool = True
) -> dict[str, Any]:
    """Identify explicit references inputs and retain omitted saved imports."""
    inputs: dict[str, Any] = {
        "implementation": implementation_identities(REPO, REFERENCES_IMPLEMENTATION)
    }
    if p.get("taxonomy"):
        inputs["taxonomy_import"] = file_identity(Path(p["taxonomy"]))
    if _has_variant_imports(p):
        inputs["variants_import"] = {
            name: file_identity(Path(p[f"validation_{name}"]))
            for name in ("clinvar", "gwas")
        }
    elif p.get("reference_sources_manifest"):
        inputs["source_manifest"] = source_manifest_identities(
            Path(p["reference_sources_manifest"]), verify_files=verify_source_files
        )

    if p.get("refresh_references"):
        return inputs
    pending = store.pending("references")
    pending_request = (
        _pending_reference_request(store, pending)
        if pending is not None and not p.get("references_rerun_explicit")
        else None
    )
    try:
        previous = store.record("references")
    except CacheError:
        if "references" in p.get("selected", []) or pending is not None:
            previous = None
        else:
            raise
    if previous is None:
        return inputs
    saved = previous.get("inputs", {})
    pending_inputs = (
        pending_request.get("reference_inputs", {})
        if isinstance(pending_request, dict)
        else {}
    )
    if not isinstance(saved, dict):
        saved = {}
    if isinstance(pending_inputs, dict):
        saved = saved | pending_inputs
    if not p.get("taxonomy") and "taxonomy_import" in saved:
        inputs["taxonomy_import"] = saved["taxonomy_import"]
    if not _has_variant_imports(p) and not p.get("reference_sources_manifest"):
        if "variants_import" in saved:
            inputs["variants_import"] = saved["variants_import"]
        elif "source_manifest" in saved:
            inputs["source_manifest"] = saved["source_manifest"]
    return inputs


def _preflight_reference_cache(p: dict[str, Any], store: ResultStore) -> None:
    """Reject changed requested inputs before any consumer tasks are scheduled."""
    if "references" not in p["enabled"]:
        return
    requested = _reference_inputs(p, store)
    p["reference_inputs"] = requested
    pending = store.pending("references")
    if pending is not None and not p.get("references_rerun_explicit"):
        previous_request = _pending_reference_request(store, pending)
        current_options = _reference_attempt_options(p)
        p["references_attempt_options"] = current_options
        if previous_request is not None:
            previous_options = previous_request.get("references_attempt_options")
            previous_inputs = previous_request.get("reference_inputs")
            previous_settings = previous_request.get("references_settings")
            if previous_options is not None and previous_options != current_options:
                raise CacheError(
                    "Unfinished references attempt options differ; use --rerun references to start a new attempt"
                )
            if previous_inputs is not None and previous_inputs != requested:
                raise CacheError(
                    "Unfinished references attempt inputs differ; use --rerun references to start a new attempt"
                )
            if previous_settings is not None:
                p["references_settings"] = previous_settings
                return
        p["references_settings"] = settings("references", p, store)
        return
    try:
        record = store.record("references")
    except CacheError:
        if "references" in p["selected"]:
            return
        raise
    if record is None:
        return
    changed = {
        name
        for name in set(record.get("inputs", {})) | set(requested)
        if record.get("inputs", {}).get(name) != requested.get(name)
    }
    requested_settings = settings("references", p, store)
    p["references_settings"] = requested_settings
    p["references_attempt_options"] = _reference_attempt_options(p)
    settings_changed = record["settings"] != requested_settings
    if "references" not in p["selected"] and (changed or settings_changed):
        message = (
            "Reference source settings differ; use --rerun references --refresh_references"
            if settings_changed
            else "Reference inputs differ from the saved references bundle; use --rerun references"
        )
        error = CacheError(message)
        raise stage_error("references", store, error, requested)
    if "references" in p["selected"]:
        if settings_changed and not p.get("refresh_references"):
            raise CacheError(
                "Changing reference source URLs requires --rerun references --refresh_references"
            )
        # Variant imports and source snapshots are known validation inputs.
        # Require their saved consumer to be selected before any analysis jobs
        # can start. A code or taxonomy change is checked against the actual
        # prepared file bytes after the references component commits.
        required = (
            ["validation"]
            if changed & {"variants_import", "source_manifest"}
            and store.record("validation") is not None
            else []
        )
        missing = [name for name in required if name not in p["selected"]]
        if missing:
            all_reruns = ["references", *required]
            raise CacheError(
                "Changed reference inputs affect saved stages: "
                + ", ".join(required)
                + ". Continue with --rerun "
                + ",".join(all_reruns)
            )


def _saved_references(store: ResultStore) -> Path | None:
    if store.record("references") is None:
        return None
    bundle = store.get("references")
    if bundle is None:
        return None
    validate_reference_bundle(bundle, verify_sources=False)
    return bundle


def _taxonomy_reuse(
    p: dict[str, Any], store: ResultStore, references: Path | None
) -> Path | None:
    if p.get("taxonomy") or p.get("refresh_references"):
        return None
    if references is not None:
        return references / "taxonomy"
    legacy = store.record("taxonomy")
    if legacy is not None:
        saved = store.get("taxonomy")
        taxonomy_manifest_identity(saved)
        return saved
    candidate = store.root / "references" / "taxonomy"
    if candidate.is_dir() and not candidate.is_symlink():
        taxonomy_manifest_identity(candidate)
        return candidate
    return None


def _retained_source_manifest(
    p: dict[str, Any], references: Path | None
) -> Path | None:
    if (
        p.get("refresh_references")
        or p.get("reference_sources_manifest")
        or _has_variant_imports(p)
        or references is None
    ):
        return None
    local = references / "sources" / "manifest.json"
    if local.is_file():
        source_manifest_identities(local)
        return local
    external_descriptor = references / "sources" / "external.json"
    if external_descriptor.is_file():
        descriptor = read_reference_json(
            external_descriptor, "external source descriptor"
        )
        source = descriptor.get("source_manifest")
        if isinstance(source, dict) and isinstance(source.get("path"), str):
            manifest = Path(source["path"]).expanduser().resolve(strict=True)
            if file_identity(manifest) != {
                key: source.get(key) for key in ("bytes", "sha256")
            }:
                raise CacheError(f"External raw source manifest changed: {manifest}")
            source_manifest_identities(manifest)
            return manifest
    return None


def _references_spec_fields(
    p: dict[str, Any], store: ResultStore, references: Path | None
) -> dict[str, Any]:
    retained_source = _retained_source_manifest(p, references)
    variants_reuse = None
    if (
        references is not None
        and not p.get("refresh_references")
        and not _has_variant_imports(p)
        and not p.get("reference_sources_manifest")
    ):
        variant_manifest_path = references / "variants" / "manifest.json"
        if variant_manifest_path.is_file():
            variant_manifest = read_reference_json(
                variant_manifest_path, "prepared variant manifest"
            )
            if variant_manifest.get("imported"):
                try:
                    validate_variant_component(references)
                except CacheError as error:
                    raise CacheError(
                        "Saved imported variants are damaged; supply both validation TSVs or restore the saved bundle"
                    ) from error
                variants_reuse = references
    taxonomy = _taxonomy_reuse(p, store, references)
    return {
        "taxonomy_reuse": str(taxonomy) if taxonomy is not None else None,
        "retained_source_manifest": (
            str(retained_source) if retained_source is not None else None
        ),
        "retained_source_root": (
            str(retained_source.parent) if retained_source is not None else None
        ),
        "variants_reuse": str(variants_reuse) if variants_reuse is not None else None,
        "requested_source_manifest": p.get("reference_sources_manifest"),
    }


def _boolean_option(value: Any, option: str) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no", ""}:
            return False
    raise CacheError(f"{option} must be a boolean value")


def prepare(p: dict[str, Any]) -> Path:
    p["run_validation"] = _boolean_option(
        p.get("run_validation", False), "--run_validation"
    )
    rerun = p.get("rerun") or ""
    if not isinstance(rerun, (str, list)):
        raise CacheError("--rerun requires a stage name, comma-separated names, or all")
    selected = set(rerun.split(",")) if isinstance(rerun, str) else set(rerun)
    selected = {s.strip() for s in selected if s.strip()}
    if "taxonomy" in selected:
        raise CacheError("--rerun taxonomy was renamed to --rerun references")
    unknown = selected - STAGES - {"all"}
    if unknown or ("all" in selected and len(selected) != 1):
        raise CacheError(
            f"Invalid --rerun selection: {rerun}. Choose {', '.join(sorted(STAGES))}, or all"
        )
    p["references_rerun_explicit"] = selected == {"all"} or "references" in selected
    for name in (
        "archive",
        "taxonomy",
        "processed_dir",
        "outdir",
        "python",
        "gmm_cache_dir",
        "download_dir",
        "reference_sources_manifest",
    ):
        if p.get(name):
            # Resolving a virtualenv interpreter symlink selects the base Python
            # and loses its installed dependencies.
            path = Path(p[name]).expanduser()
            p[name] = str(path.absolute() if name == "python" else path.resolve())
    p["refresh_references"] = _boolean_option(
        p.get("refresh_references", False), "--refresh_references"
    )
    cache_dir = p.get("references_cache_dir") or (
        SCRIPTS.parent / "tmp" / "nextflow" / "references"
    )
    p["references_cache_dir"] = str(Path(cache_dir).expanduser().resolve())
    if int(p["gmm_max_parallel"]) < 1:
        raise CacheError("--gmm_max_parallel must be positive")
    if p.get("archive_sha256"):
        import re

        if not re.fullmatch(r"[a-fA-F0-9]{64}", p["archive_sha256"]):
            raise CacheError(
                "--archive_sha256 must contain exactly 64 hexadecimal characters"
            )
        p["archive_sha256"] = p["archive_sha256"].lower()
    enabled = (
        {"ingest"}
        if p["ingest_only"]
        else (
            {"process"}
            if p.get("processed_dir")
            else {"ingest", "references", "process", "tables"}
        )
    )
    for flag, names in (
        ("run_gmm", {"gmm", "merge"}),
        ("run_validation", {"references", "gmm", "merge", "validation"}),
        ("compatibility_tables", {"tables"}),
        ("run_qc", {"tables", "qc"}),
        ("run_reports", {"gmm", "merge", "reports"}),
    ):
        if p.get(flag):
            enabled |= names
    if selected != {"all"}:
        enabled |= selected
    if enabled & {"gmm", "merge", "reports", "validation"}:
        enabled |= {"gmm", "merge"}
    if "validation" in enabled:
        enabled.add("references")
    if "reports" in enabled:
        enabled.add("process")
    if "qc" in enabled:
        enabled |= {"tables", "process"}
    if p["ingest_only"] and enabled != {"ingest"}:
        raise CacheError(
            "--ingest_only cannot be combined with processing or analysis stages"
        )
    if p.get("processed_dir"):
        if p["ingest_only"]:
            raise CacheError("--processed_dir cannot be combined with --ingest_only")
        if selected & {"ingest", "tables"} or "tables" in enabled:
            raise CacheError(
                "--processed_dir cannot regenerate source stages or derived tables; supply the original archive instead"
            )
        enabled.discard("ingest")
    _prepare_reference_imports(p)
    if "references" not in enabled and (
        _has_variant_imports(p) or p.get("reference_sources_manifest")
    ):
        raise CacheError(
            "Reference imports require the references stage; use --rerun references or --run_validation"
        )
    if _has_variant_imports(p) and p.get("reference_sources_manifest"):
        raise CacheError(
            "--reference_sources_manifest cannot be combined with paired validation TSV imports"
        )
    if "validation" in enabled:
        _prepare_validation(p)
    p["enabled"] = sorted(enabled)
    p["selected"] = sorted(enabled if selected == {"all"} else selected)
    p["run_id"] = uuid.uuid4().hex
    store = ResultStore(Path(p["outdir"]))
    pending_references = store.pending("references")
    if (
        p["refresh_references"]
        and not p["references_rerun_explicit"]
        and pending_references is None
    ):
        raise CacheError("--refresh_references requires explicit --rerun references")
    from data_sources.mavedb.scripts._pipeline_inputs import reusable_processing

    try:
        p["saved_processing"] = reusable_processing(p, store)
        if p["saved_processing"]:
            provenance_notice("process", Path(p["saved_processing"]))
    except CacheError as error:
        raise stage_error("process", store, error) from error
    _preflight_reference_cache(p, store)
    if "references" in enabled:
        if p["references_rerun_explicit"]:
            p["references_attempt"] = store.begin("references", rerun=True)
        elif pending_references is not None:
            p["references_attempt"] = store.begin("references")
        elif store.record("references") is None:
            p["references_attempt"] = store.begin("references")
        else:
            p["references_attempt"] = None
        p["references_attempt_options"] = _reference_attempt_options(p)
        p.setdefault("references_settings", settings("references", p, store))
    if "gmm" in enabled:
        # Validate configuration early, but never fit in the policy process.
        settings("gmm", p)
        if (
            "gmm" in p["selected"]
            or store.pending("gmm")
            or store.record("gmm") is None
        ):
            p["gmm_attempt"] = store.begin("gmm", rerun="gmm" in p["selected"])
        else:
            p["gmm_attempt"] = store.record("gmm")["settings"].get("attempt", "adopted")
    path = Path(p["outdir"]) / "pipeline_requests" / f"{p['run_id']}.json"
    write_json(path, p)
    return path


def stage_key(stage: str, p: dict[str, Any], index: int | None = None) -> str:
    if stage in {"gmm-plan", "gmm-batch"}:
        return f"{stage}-{p['gmm_attempt']}" + (
            f"-{index}" if index is not None else ""
        )
    return stage


def adopt(
    stage: str, p: dict[str, Any], store: ResultStore, paths: dict[str, str]
) -> None:
    if stage not in {"taxonomy", "process", "gmm", "tables"}:
        return
    from data_sources.mavedb.scripts._pipeline_adopt import adopt_processing, unpack_gmm

    root = Path(p["outdir"])
    if stage in p["selected"] or store.record(stage) is not None:
        return
    if stage == "taxonomy" and (root / "references/taxonomy").exists():
        # Pre-registry runs already published the snapshot used for processing.
        # Adopt it before scheduling a fetch from NCBI's moving reference URL.
        legacy = root / "references/taxonomy"
        manifest = read_json(legacy / "manifest.json")
        actual = file_identity(legacy / "taxcat.dmp")
        if actual != {key: manifest.get(key) for key in ("bytes", "sha256")}:
            raise CacheError(
                "Saved taxonomy differs from its manifest; restore it or explicitly use --rerun taxonomy"
            )
        source_inputs = {}
        if p.get("taxonomy"):
            if file_identity(Path(p["taxonomy"])) != actual:
                raise CacheError(
                    "Requested taxonomy differs from the saved snapshot; use --rerun taxonomy"
                )
            source_inputs["source"] = actual
        elif (
            p["taxonomy_url"] != DEFAULT_TAXONOMY_URL
            and manifest.get("source") != p["taxonomy_url"]
        ):
            raise CacheError(
                "Requested taxonomy source differs from the saved snapshot; supply the original --taxonomy or use --rerun taxonomy"
            )
        store.commit(stage, legacy, settings(stage, p), source_inputs)
    if (
        stage == "process"
        and not p.get("processed_dir")
        and (root / "processed").exists()
    ):
        manifest = adopt_processing(root / "processed")
        deps = dependencies(stage, paths)
        raw_manifest = read_json(Path(paths["ingest"]) / "manifest.json")
        if manifest["ingest"]["archive_sha256"] != raw_manifest[
            "archive_sha256"
        ] or manifest["ingest"]["metadata"] != file_identity(
            Path(paths["ingest"]) / "mavedb_dump/main.json"
        ):
            raise CacheError(
                "Existing processing used a different source archive; use --rerun process"
            )
        if deps.get("taxonomy") != {
            k: manifest["taxonomy"][k] for k in ("bytes", "sha256")
        }:
            raise CacheError(
                "Existing processing used a different taxonomy; restore its matching snapshot with --taxonomy PATH --rerun taxonomy. Use --rerun process only when deliberately updating the processed data"
            )
        store.commit(stage, root / "processed", settings(stage, p), deps)
    if stage == "gmm" and (root / "gmm").exists():
        with tempfile.TemporaryDirectory(prefix="mavedb-adopt-") as temp:
            collected = Path(temp) / "collected"
            plan = unpack_gmm(root / "gmm", collected)
            expected = gmm_configuration(p)
            if plan["configuration"] != expected:
                raise CacheError(
                    "Existing GMM settings differ; use the original settings, a new --outdir, or --rerun gmm"
                )
            deps = dependencies("gmm", paths)
            if plan["inputs"] != {
                "scores": deps["scores.tsv"],
                "metadata": deps["metadata.tsv"],
            }:
                raise CacheError("Existing GMM inputs differ; use --rerun gmm,merge")
            store.commit("gmm", collected, settings("gmm", p), deps)
            if "merge" not in p["selected"] and store.record("merge") is None:
                native = store.get("gmm")
                store.commit(
                    "merge",
                    root / "gmm",
                    {},
                    dependencies("merge", {"gmm": str(native)}),
                )

    if stage == "tables" and "process" in paths:
        canonical = Path(paths["process"])
        manifest_path = canonical / "processing_manifest.json"
        if manifest_path.is_file() and read_json(manifest_path).get("settings", {}).get(
            "tables",
            read_json(manifest_path)
            .get("settings", {})
            .get("compatibility_tables", False),
        ):
            manifest = adopt_processing(canonical)
            with tempfile.TemporaryDirectory(prefix="mavedb-adopt-tables-") as temp:
                target = Path(temp)
                for name in ("singles", "doubles", "wide"):
                    destination = target / "maves" / f"{name}.tsv"
                    destination.parent.mkdir(exist_ok=True)
                    shutil.copyfile(canonical / "maves" / f"{name}.tsv", destination)
                write_json(
                    target / "completion.json",
                    {
                        "schema_version": 1,
                        "adopted_processing": file_identity(manifest_path),
                        "outputs": {
                            name: value
                            for name, value in manifest["outputs"].items()
                            if name.startswith("maves/")
                        },
                    },
                )
                store.commit(stage, target, {}, dependencies(stage, paths))


def _check_reference_upgrade(
    store: ResultStore, repairing_registered_result: bool = False
) -> None:
    """Refuse to replace unregistered siblings in a legacy references folder."""
    directory = store.root / "references"
    if repairing_registered_result:
        try:
            state = read_json(store.state_path)
        except (OSError, ValueError):
            state = {}
        results = state.get("results", {}) if isinstance(state, dict) else {}
        if isinstance(results, dict) and "references" in results:
            return
    else:
        try:
            if store.record("references") is not None:
                return
        except CacheError:
            raise
    if not directory.is_dir():
        return
    unexpected = sorted(
        path.name for path in directory.iterdir() if path.name != "taxonomy"
    )
    if unexpected:
        raise CacheError(
            "Cannot replace unregistered files in the legacy references directory: "
            + ", ".join(unexpected)
            + ". Move them to a safe location before creating the references bundle"
        )


def _resolve_references(
    request: Path, p: dict[str, Any], store: ResultStore
) -> dict[str, Any]:
    requested_inputs = p.get("reference_inputs") or _reference_inputs(
        p, store, verify_source_files=False
    )
    selected = "references" in p.get("selected", [])
    pending = store.pending("references")
    resuming = pending is not None
    existing: Path | None = None
    if not selected and not resuming:
        try:
            existing = store.get(
                "references", settings("references", p, store), requested_inputs
            )
            if existing is not None:
                validate_reference_bundle(
                    existing,
                    verify_sources=(
                        (existing / "sources" / "external.json").is_file()
                        and not p.get("reference_sources_manifest")
                    ),
                )
                LOGGER.info("REUSE references")
                return {"path": str(existing), "stage": "references"}
        except CacheError as error:
            raise stage_error("references", store, error, requested_inputs) from error
    else:
        existing = store.root / "references"
        if not existing.is_dir():
            existing = None

    _check_reference_upgrade(store, repairing_registered_result=selected or resuming)
    p["references_attempt"] = (
        p.get("references_attempt")
        or pending
        or store.begin("references", rerun=selected)
    )
    spec = {
        "request": str(request),
        "stage": "references",
        "key": "references",
        "inputs": {},
        "identities": requested_inputs,
        "settings": p.get("references_settings") or settings("references", p, store),
        "index": None,
    }
    spec.update(_references_spec_fields(p, store, existing))
    spec["references_component_cache"] = str(
        Path(p["references_cache_dir"]) / p["references_attempt"]
    )
    spec_path = (
        Path(p["outdir"]) / "pipeline_requests" / p["run_id"] / "references.json"
    )
    write_json(spec_path, spec)
    return {"spec": str(spec_path), "stage": "references"}


def resolve(
    request: Path, stage: str, paths: dict[str, str], index: int | None = None
) -> dict[str, Any]:
    p = read_json(request)
    store = ResultStore(Path(p["outdir"]))
    if stage == "references":
        return _resolve_references(request, p, store)
    public = PUBLIC.get(stage, stage)
    key = stage_key(stage, p, index)
    if (
        stage == "taxonomy"
        and not p.get("taxonomy")
        and p["taxonomy_url"] == DEFAULT_TAXONOMY_URL
        and stage not in p["selected"]
        and store.record(stage) is not None
    ):
        # An omitted override means continue the saved analysis. The default
        # acquisition URL must not replace an explicitly imported snapshot or
        # require its original task work directory to remain available.
        try:
            cached = store.get(stage)
        except CacheError as error:
            raise stage_error(stage, store, error) from error
        LOGGER.info("REUSE saved taxonomy snapshot")
        return {"path": str(cached), "stage": stage}
    if stage == "process" and p.get("processed_dir"):
        external = Path(p["processed_dir"])
        config = {"external": True}
        if (
            not external.exists()
            and store.record("process") is not None
            and stage not in p["selected"]
        ):
            return {"path": str(store.get("process", config)), "stage": stage}
        deps = {
            name: file_identity(external / name)
            for name in ("scores.tsv", "metadata.tsv")
        }
        try:
            cached = (
                None if stage in p["selected"] else store.get("process", config, deps)
            )
        except CacheError as error:
            raise stage_error("process", store, error) from error
        if cached is None:
            with tempfile.TemporaryDirectory(prefix="mavedb-canonical-") as temp:
                snapshot = Path(temp)
                for name in (
                    "scores.tsv",
                    "metadata.tsv",
                    "processing_manifest.json",
                    "dataset_processing.tsv",
                ):
                    if (external / name).is_file():
                        shutil.copyfile(external / name, snapshot / name)
                if {name: file_identity(snapshot / name) for name in deps} != deps:
                    raise CacheError(
                        "Canonical inputs changed while copying; retry with stable inputs"
                    )
                cached = store.commit("process", snapshot, config, deps)
        return {"path": str(cached), "stage": stage}
    if stage == "process" and stage not in p["selected"]:
        record = store.record(stage)
        if record and record["settings"].get("adopted"):
            from data_sources.mavedb.scripts._pipeline_inputs import (
                validate_adopted_processing,
                validate_registered_sources,
            )

            cached = store.get(stage)
            provenance = validate_adopted_processing(
                p, cached, read_json(cached / "processing_manifest.json")
            )
            validate_registered_sources(store, record, provenance)
            return {"path": str(cached), "stage": stage}
    if stage == "gmm-plan":
        try:
            adopt("gmm", p, store, paths)
        except (OSError, ValueError, KeyError) as error:
            raise stage_error("gmm", store, error) from error
        # A complete collection needs no planning, code snapshot or fit tasks.
        if "gmm" not in p["selected"] and not store.pending("gmm"):
            try:
                ready = store.get("gmm", settings("gmm", p), dependencies("gmm", paths))
            except CacheError as error:
                raise stage_error("gmm", store, error) from error
            if ready is not None:
                provenance_notice("gmm", ready)
                return {"path": str(ready), "stage": stage, "collected": True}
    adopt(stage, p, store, paths)
    deps = dependencies(stage, paths)
    if stage in {"ingest", "taxonomy"}:
        source = p.get("archive" if stage == "ingest" else "taxonomy")
        if source:
            if Path(source).is_file():
                deps["source"] = file_identity(Path(source))
            elif store.record(key) is not None:
                deps = store.record(key)["inputs"]
            else:
                raise CacheError(f"Required source file is missing: {source}")
    config = settings(stage, p)
    # Immutable attempt-specific keys already implement explicit GMM regeneration.
    forced = (
        public in p["selected"] or (stage == "gmm" and store.pending("gmm"))
    ) and stage not in {"gmm-plan", "gmm-batch"}
    try:
        cached = None if forced else store.get(key, config, deps)
    except CacheError as error:
        raise stage_error(public, store, error) from error
    if cached is not None:
        LOGGER.info("REUSE %s", key)
        return {"path": str(cached), "stage": stage}
    spec = {
        "request": str(request),
        "stage": stage,
        "key": key,
        "inputs": paths,
        "identities": deps,
        "settings": config,
        "index": index,
    }
    spec_path = Path(p["outdir"]) / "pipeline_requests" / p["run_id"] / f"{key}.json"
    write_json(spec_path, spec)
    return {"spec": str(spec_path), "stage": stage}


def register(spec_path: Path, output: Path, required: Sequence[str] = ()) -> Path:
    """Validate and checkpoint outputs produced directly by a Nextflow module."""
    spec = read_json(spec_path)
    p = read_json(spec["request"])
    stage = spec["stage"]
    if spec["key"] != stage_key(stage, p, spec.get("index")):
        raise CacheError("Stage registration key differs from its saved request")
    if stage == "references":
        store = ResultStore(Path(p["outdir"]))
        current_inputs = _reference_inputs(p, store)
        current_settings = p.get("references_settings") or settings(stage, p, store)
    else:
        current_inputs = dependencies(stage, spec["inputs"])
        current_settings = settings(stage, p)
    if stage in {"ingest", "taxonomy"}:
        source = p.get("archive" if stage == "ingest" else "taxonomy")
        if source and Path(source).is_file():
            current_inputs["source"] = file_identity(Path(source))
    if current_inputs != spec["identities"]:
        raise CacheError(
            f"Inputs changed while {stage} was running; retry with stable inputs"
        )
    if current_settings != spec["settings"]:
        raise CacheError("Stage registration settings differ from its saved request")
    output = Path(output)
    if not output.is_dir():
        raise CacheError(f"Stage {stage} produced no result directory: {output}")
    for name in required:
        relative = Path(name)
        candidate = output / relative
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or not candidate.resolve().is_relative_to(output.resolve())
            or candidate.is_symlink()
            or not candidate.is_file()
        ):
            raise CacheError(
                f"Stage {stage} is missing a required regular output file: {name}"
            )
    if stage == "references":
        validate_reference_bundle(
            output,
            verify_sources=not bool(
                p.get("reference_sources_manifest") and not _has_variant_imports(p)
            ),
        )
    LOGGER.info("REGISTER %s", spec["key"])
    return ResultStore(Path(p["outdir"])).commit(
        spec["key"], Path(output), spec["settings"], spec["identities"]
    )


def finish(request: Path, paths: dict[str, str]) -> None:
    p = read_json(request)
    store = ResultStore(Path(p["outdir"]))
    # Check existing descendants even when their execution was not requested.
    current = dict(paths)
    stale = []
    for stage in STAGES:
        try:
            exists = store.record(stage) is not None
        except CacheError:
            if stage == "references":
                stale.append("references")
                continue
            raise
        if stage not in current and exists:
            try:
                current[stage] = str(store.get(stage))
            except CacheError:
                if stage == "references":
                    stale.append("references")
                else:
                    raise
    try:
        references_record = store.record("references")
    except CacheError:
        references_record = None
        if "references" not in stale:
            stale.append("references")
    if references_record is not None:
        try:
            checked_references = store.get("references")
            if checked_references is None:
                raise CacheError("Saved references result is missing")
            validate_reference_bundle(
                checked_references,
                verify_sources=(
                    (checked_references / "sources" / "external.json").is_file()
                    and not p.get("reference_sources_manifest")
                ),
            )
            current["references"] = str(checked_references)
        except CacheError:
            stale.append("references")
            outputs = references_record.get("outputs", {})
            bundle = store.root / "references"
            taxonomy_name = "taxonomy/taxcat.dmp"
            taxonomy_expected = outputs.get(taxonomy_name)
            try:
                taxonomy_changed = (
                    taxonomy_expected is None
                    or file_identity(bundle / taxonomy_name) != taxonomy_expected
                )
            except OSError:
                taxonomy_changed = True
            variants_changed = False
            for filename in VARIANT_FILES.values():
                output_name = f"variants/{filename}"
                try:
                    if (
                        outputs.get(output_name) is None
                        or file_identity(bundle / output_name) != outputs[output_name]
                    ):
                        variants_changed = True
                except OSError:
                    variants_changed = True
            if taxonomy_changed:
                stale.extend(
                    name
                    for name in ("process", "tables")
                    if store.record(name) is not None and name not in stale
                )
            if variants_changed and store.record("validation") is not None:
                stale.append("validation")
    for stage in (
        "process",
        "tables",
        "gmm",
        "merge",
        "validation",
        "qc",
        "reports",
    ):
        record = store.record(stage)
        if record is None:
            continue
        try:
            if stage == "process" and record["settings"].get("adopted"):
                from data_sources.mavedb.scripts._pipeline_inputs import (
                    validate_adopted_processing,
                    validate_registered_sources,
                )

                actual = validate_adopted_processing(
                    p,
                    Path(current["process"]),
                    read_json(Path(current["process"]) / "processing_manifest.json"),
                )
                validate_registered_sources(store, record, actual)
            elif stage == "process" and record["settings"].get("external"):
                actual = {
                    name: file_identity(Path(current["process"]) / name)
                    for name in ("scores.tsv", "metadata.tsv")
                }
            elif stage == "validation":
                validation_paths = dict(current)
                if "references" not in validation_paths and "validation" in p.get(
                    "enabled", []
                ):
                    for name in ("clinvar", "gwas"):
                        if p.get(f"validation_{name}"):
                            validation_paths[name] = p[f"validation_{name}"]
                elif "references" not in validation_paths:
                    reference_paths = record["settings"].get("reference_paths", {})
                    for name in ("clinvar", "gwas"):
                        if name not in validation_paths and reference_paths.get(name):
                            validation_paths[name] = reference_paths[name]
                actual = dependencies(stage, validation_paths)
            else:
                actual = dependencies(stage, current)
            # The plan is an internal fitting checkpoint, not a dependency of a completed collection.
            recorded = {k: v for k, v in record["inputs"].items() if k != "plan"}
            if actual != recorded:
                stale.append(stage)
        except (KeyError, CacheError):
            stale.append(stage)
    if stale:
        # Include transitive descendants so the diagnostic gives a complete next command.
        for parent, children in (
            ("process", ["tables", "gmm", "qc", "reports", "validation"]),
            ("tables", ["qc"]),
            ("gmm", ["merge", "reports", "validation"]),
            ("merge", ["reports", "validation"]),
        ):
            if parent in stale:
                stale.extend(s for s in children if store.record(s) and s not in stale)
        if "references" in stale:
            # Consumers already present above are the only ones whose direct
            # reference files changed. A changed ClinVar/GWAS table alone does
            # not invalidate taxonomy-based processing or GMM.
            stale = list(dict.fromkeys(stale))
        raise CacheError(
            "Changed inputs affect saved results: "
            + ", ".join(stale)
            + ". Continue with --rerun "
            + ",".join(stale)
        )
    # Stages already own their ordinary output directories. Completion of the
    # workflow validates their lineage; no data-copying publication layer exists.
    root = Path(p["outdir"])
    write_json(
        root / "pipeline_publication.json",
        {
            "schema_version": 2,
            "run_id": p["run_id"],
            "status": "completed",
            "results": {
                name: os.path.relpath(path, root) for name, path in current.items()
            },
            "manifests": {
                name: file_identity(Path(path) / "pipeline_manifest.json")
                for name, path in current.items()
                if (Path(path) / "pipeline_manifest.json").is_file()
            },
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=("prepare", "resolve", "register", "finish", "lock")
    )
    parser.add_argument("--root", type=Path)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--stage")
    parser.add_argument("--index", type=int)
    parser.add_argument("--spec", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require", action="append", default=[])
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        if args.action == "lock":
            with ResultStore(args.root).lock():
                print("LOCKED", flush=True)
                sys.stdin.read()
            return
        payload = json.load(sys.stdin) if args.action != "register" else {}
        if args.action == "prepare":
            result = {"request": str(prepare(payload))}
        elif args.action == "resolve":
            result = resolve(args.request, args.stage, payload, args.index)
        elif args.action == "register":
            result = {
                "path": str(register(args.spec, args.output, args.require)),
                "stage": read_json(args.spec)["stage"],
            }
        else:
            finish(args.request, payload)
            result = {"status": "completed"}
    except (
        CacheError,
        OSError,
        ValueError,
        KeyError,
    ) as error:
        parser.exit(1, f"MaveDB pipeline: {error}\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
