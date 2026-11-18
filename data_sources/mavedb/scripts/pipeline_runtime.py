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
from data_sources.mavedb.scripts.download_taxonomy import (
    DEFAULT_URL as DEFAULT_TAXONOMY_URL,
)

LOGGER = logging.getLogger(__name__)
SCRIPTS = Path(__file__).resolve().parent
STAGES = {"ingest", "taxonomy", "process", "tables", "gmm", "merge", "qc", "reports"}
PUBLIC = {"gmm-plan": "gmm", "gmm-batch": "gmm"}


def stage_error(stage: str, store: ResultStore, error: Exception) -> CacheError:
    descendants = {
        "ingest": ["process", "tables", "gmm", "merge", "qc", "reports"],
        "taxonomy": ["process", "tables", "gmm", "merge", "qc", "reports"],
        "process": ["tables", "gmm", "merge", "qc", "reports"],
        "tables": ["qc"],
        "gmm": ["merge", "reports"],
        "merge": ["reports"],
    }
    affected = [stage]
    for name in descendants.get(stage, []):
        if store.record(name) is not None:
            affected.append(name)
    return CacheError(
        f'{error}. Affected stages: {", ".join(affected)}. Continue with --rerun {",".join(affected)}'
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


def settings(stage: str, p: dict[str, Any]) -> dict[str, Any]:
    if stage == "ingest":
        return {
            "source": p.get("archive") or p["zenodo_url"],
            "expected_sha256": p.get("archive_sha256"),
        }
    if stage == "taxonomy":
        return {"source": p.get("taxonomy") or p["taxonomy_url"]}
    if stage in {"gmm", "gmm-plan", "gmm-batch"}:
        return gmm_configuration(p) | {
            "fail_on_study_error": p["gmm_fail_on_study_error"]
        }
    return {}


def dependencies(stage: str, paths: dict[str, str]) -> dict[str, Any]:
    found: dict[str, Any] = {}

    def add(name: str, path: Path) -> None:
        found[name] = identity(path)

    if stage in {"process", "tables"} and "ingest" in paths:
        add("raw", Path(paths["ingest"]) / "mavedb_dump")
    if stage in {"process", "tables"} and "taxonomy" in paths:
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
    return found


def prepare(p: dict[str, Any]) -> Path:
    rerun = p.get("rerun") or ""
    if not isinstance(rerun, (str, list)):
        raise CacheError("--rerun requires a stage name, comma-separated names, or all")
    selected = set(rerun.split(",")) if isinstance(rerun, str) else set(rerun)
    selected = {s.strip() for s in selected if s.strip()}
    unknown = selected - STAGES - {"all"}
    if unknown or ("all" in selected and len(selected) != 1):
        raise CacheError(
            f"Invalid --rerun selection: {rerun}. Choose {', '.join(sorted(STAGES))}, or all"
        )
    for name in (
        "archive",
        "taxonomy",
        "processed_dir",
        "outdir",
        "python",
        "gmm_cache_dir",
        "download_dir",
    ):
        if p.get(name):
            # Resolving a virtualenv interpreter symlink selects the base Python
            # and loses its installed dependencies.
            path = Path(p[name]).expanduser()
            p[name] = str(path.absolute() if name == "python" else path.resolve())
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
            else {"ingest", "taxonomy", "process", "tables"}
        )
    )
    for flag, names in (
        ("run_gmm", {"gmm", "merge"}),
        ("compatibility_tables", {"tables"}),
        ("run_qc", {"tables", "qc"}),
        ("run_reports", {"gmm", "merge", "reports"}),
    ):
        if p.get(flag):
            enabled |= names
    if selected != {"all"}:
        enabled |= selected
    if enabled & {"gmm", "merge", "reports"}:
        enabled |= {"gmm", "merge"}
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
        if selected & {"ingest", "taxonomy", "tables"} or "tables" in enabled:
            raise CacheError(
                "--processed_dir cannot regenerate source stages or derived tables; supply the original archive instead"
            )
        enabled -= {"ingest", "taxonomy"}
    p["enabled"] = sorted(enabled)
    p["selected"] = sorted(enabled if selected == {"all"} else selected)
    p["run_id"] = uuid.uuid4().hex
    store = ResultStore(Path(p["outdir"]))
    from data_sources.mavedb.scripts._pipeline_inputs import reusable_processing

    try:
        p["saved_processing"] = reusable_processing(p, store)
        if p["saved_processing"]:
            provenance_notice("process", Path(p["saved_processing"]))
    except CacheError as error:
        raise stage_error("process", store, error) from error
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


def resolve(
    request: Path, stage: str, paths: dict[str, str], index: int | None = None
) -> dict[str, Any]:
    p = read_json(request)
    store = ResultStore(Path(p["outdir"]))
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
    current_inputs = dependencies(stage, spec["inputs"])
    if stage in {"ingest", "taxonomy"}:
        source = p.get("archive" if stage == "ingest" else "taxonomy")
        if source and Path(source).is_file():
            current_inputs["source"] = file_identity(Path(source))
    if current_inputs != spec["identities"]:
        raise CacheError(
            f"Inputs changed while {stage} was running; retry with stable inputs"
        )
    if settings(stage, p) != spec["settings"]:
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
    LOGGER.info("REGISTER %s", spec["key"])
    return ResultStore(Path(p["outdir"])).commit(
        spec["key"], Path(output), spec["settings"], spec["identities"]
    )


def finish(request: Path, paths: dict[str, str]) -> None:
    p = read_json(request)
    store = ResultStore(Path(p["outdir"]))
    # Check existing descendants even when their execution was not requested.
    current = dict(paths)
    for stage in STAGES:
        if stage not in current and store.record(stage) is not None:
            current[stage] = str(store.get(stage))
    stale = []
    for stage in ("process", "tables", "gmm", "merge", "qc", "reports"):
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
            ("process", ["tables", "gmm", "qc", "reports"]),
            ("tables", ["qc"]),
            ("gmm", ["merge", "reports"]),
            ("merge", ["reports"]),
        ):
            if parent in stale:
                stale.extend(s for s in children if store.record(s) and s not in stale)
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
