"""Input adapters for the MAVE homology analysis; no model fitting is required."""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from data_sources.mavedb.analysis.mave_summary.data import _load_canonical_meta
from data_sources.mavedb.scripts.process_mavedb import AA_3_TO_1, SINGLE_PATTERN

LOGGER = logging.getLogger(__name__)


@dataclass
class HomologyInputs:
    """Normalized selected-model tables and their source files."""

    meta: pd.DataFrame
    scores: pd.DataFrame
    singles: pd.DataFrame
    params: pd.DataFrame
    paths: list[Path]


def read_tsv(path: Path, required: set[str] | None = None) -> pd.DataFrame:
    """Read TSVs without interpreting amino-acid strings as missing values."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input table not found: {path}")
    frame = pd.read_csv(
        path, sep="\t", keep_default_na=False, na_values=["", "NA"], low_memory=False
    )
    if required and (missing := required - set(frame.columns)):
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    return frame


def boolean_values(values: pd.Series, name: str) -> pd.Series:
    """Read explicit R/Python boolean spellings and preserve missing values."""
    text = values.astype("string").str.strip().str.lower()
    mapped = text.map(
        {"true": True, "t": True, "1": True, "false": False, "f": False, "0": False}
    )
    if (text.notna() & mapped.isna()).any():
        raise ValueError(f"Invalid boolean value in {name}")
    return mapped.astype("boolean")


def _selected(frame: pd.DataFrame, *, required: bool = False) -> pd.DataFrame:
    if "top_model" not in frame:
        if required:
            raise ValueError("Native GMM tables require a top_model column")
        return frame
    return frame.loc[boolean_values(frame["top_model"], "top_model").fillna(False)]


def _metadata(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0)
    if "dataset_id" in header:
        meta = _load_canonical_meta(path)
    else:
        meta = read_tsv(path, {"score_urn"})
    aliases = {
        "experiment_id": "experiment_urn",
        "experiment_set_id": "experiment_set_urn",
    }
    for old, new in aliases.items():
        if new not in meta and old in meta:
            meta[new] = meta[old]
    for column in (
        "experiment_urn",
        "experiment_set_urn",
        "organism",
        "target_name",
        "identifier_uniprot",
    ):
        if column not in meta:
            meta[column] = pd.NA
    if "meta_study" not in meta:
        meta["meta_study"] = False
    meta["meta_study"] = boolean_values(meta["meta_study"], "meta_study")
    if meta["score_urn"].isna().any() or meta["score_urn"].duplicated().any():
        raise ValueError("Metadata must contain one nonmissing score_urn per score set")
    # R filter(!meta_study) drops unknown as well as true values.
    return meta.loc[meta["meta_study"].eq(False).fillna(False)].copy()


def _native_paths(gmm_dir: Path) -> tuple[list[Path], list[Path]]:
    tables: dict[str, dict[int, Path]] = {"score": {}, "fitted_params": {}}
    for prefix in tables:
        for path in gmm_dir.glob(f"{prefix}*.tsv"):
            match = re.fullmatch(rf"{prefix}(\d+)\.tsv", path.name)
            if match:
                index = int(match[1])
                if index in tables[prefix]:
                    raise ValueError(f"Duplicate {prefix} batch index {index}")
                tables[prefix][index] = path
    if not tables["score"] or tables["score"].keys() != tables["fitted_params"].keys():
        raise ValueError(
            f"{gmm_dir}: require matching score<N>.tsv and fitted_params<N>.tsv batches"
        )
    indices = sorted(tables["score"])
    return tuple([tables[prefix][i] for i in indices] for prefix in tables)


def _native_scores(paths: list[Path], urns: set[str]) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    seen: set[str] = set()
    for path in paths:
        file_urns: set[str] = set()
        for chunk in pd.read_csv(
            path,
            sep="\t",
            keep_default_na=False,
            na_values=["", "NA"],
            chunksize=100_000,
        ):
            required = {"score_urn", "hgvs_pro", "w", "top_model", "model"}
            if missing := required - set(chunk.columns):
                raise ValueError(f"{path}: missing columns: {sorted(missing)}")
            file_urns.update(chunk["score_urn"].dropna())
            chunk = _selected(chunk, required=True)
            chunk = chunk.loc[chunk["score_urn"].isin(urns)]
            parts.append(chunk)
        if overlap := seen & file_urns:
            raise ValueError(
                f"Overlapping score sets in GMM batches: {sorted(overlap)[:5]}"
            )
        seen.update(file_urns)
    if not parts:
        return pd.DataFrame(columns=["score_urn", "hgvs_pro", "w"])
    return pd.concat(parts, ignore_index=True)


def annotate_positions(scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Parse simple HGVS calls into UniProt coordinates on each original row."""
    scores = scores.copy()
    parsed = scores["hgvs_pro"].astype("string").str.extract(SINGLE_PATTERN)
    parsed.columns = ["wt", "position", "mut"]
    parsed["mut"] = parsed["mut"].mask(parsed["mut"].eq("="), parsed["wt"])
    parsed["wt"] = parsed["wt"].str.title().map(AA_3_TO_1)
    mutant = parsed["mut"].str.title()
    mutant = mutant.mask(parsed["mut"].str.lower().eq("del"), "del")
    parsed["mut"] = mutant.map(AA_3_TO_1)
    offsets = (
        pd.to_numeric(meta.set_index("score_urn")["offset_uniprot"], errors="coerce")
        if "offset_uniprot" in meta
        else pd.Series(dtype=float)
    )
    parsed["position"] = pd.to_numeric(parsed["position"], errors="coerce") + (
        scores["score_urn"].map(offsets).fillna(0)
    )
    for column in ("position", "wt", "mut"):
        scores[column] = parsed[column]
    return scores


def load_inputs(
    *,
    metadata: Path,
    gmm_dir: Path | None = None,
    data_root: Path | None = None,
    model: str = "modes",
) -> HomologyInputs:
    """Accept native batch exports or consolidated legacy R tables."""
    if (gmm_dir is None) == (data_root is None):
        raise ValueError("Choose exactly one of gmm_dir and data_root")
    meta = _metadata(metadata)
    paths = [metadata]
    legacy_singles: pd.DataFrame | None = None
    if gmm_dir is not None:
        score_paths, param_paths = _native_paths(gmm_dir)
        paths.extend([*score_paths, *param_paths])
        params = pd.concat(
            [
                _selected(read_tsv(p, {"score_urn", "sigma_neut"}), required=True)
                for p in param_paths
            ],
            ignore_index=True,
        )
        scores = _native_scores(score_paths, set(meta["score_urn"]))
    else:
        root = data_root / "maves" / f"gmm_{model}"
        score_path, singles_path, param_path = (
            root / f"{name}.tsv" for name in ("scores", "singles", "params")
        )
        paths.extend([score_path, singles_path, param_path])
        scores = _selected(read_tsv(score_path, {"score_urn", "hgvs_pro", "w"}))
        params = _selected(
            read_tsv(param_path, {"score_urn", "sigma_neut", "top_model"})
        )
        legacy_singles = read_tsv(
            singles_path, {"score_urn", "position", "wt", "mut", "w"}
        )
    params = params.loc[params["score_urn"].isin(meta["score_urn"])].copy()
    if params["score_urn"].duplicated().any():
        raise ValueError("More than one selected model per score set in parameters")
    params["sigma_neut"] = pd.to_numeric(params["sigma_neut"], errors="coerce")
    if not (np.isfinite(params["sigma_neut"]) & params["sigma_neut"].gt(0)).all():
        raise ValueError("Selected model sigma_neut must be finite and positive")
    scores = scores.loc[scores["score_urn"].isin(meta["score_urn"])].copy()
    scores["w"] = pd.to_numeric(scores["w"], errors="coerce")
    scores = scores.loc[np.isfinite(scores["w"]) & scores["hgvs_pro"].notna()].copy()
    if missing := set(scores["score_urn"]) - set(params["score_urn"]):
        raise ValueError(
            f"Finite fitted scores without selected parameters: {sorted(missing)[:5]}"
        )
    if gmm_dir is not None and not scores.empty:
        if "model" not in params:
            raise ValueError("Native parameter tables require a model column")
        expected = scores["score_urn"].map(params.set_index("score_urn")["model"])
        if (scores["model"].isna() | scores["model"].ne(expected)).any():
            raise ValueError("Selected score and parameter models do not agree")
    if gmm_dir is not None or not {"position", "wt", "mut"}.issubset(scores.columns):
        scores = annotate_positions(scores, meta)
    if legacy_singles is None:
        singles = (
            scores.dropna(subset=["position", "wt", "mut"])
            .groupby(
                ["score_urn", "position", "wt", "mut"], as_index=False, sort=False
            )["w"]
            .mean()
        )
    else:
        singles = legacy_singles.loc[
            legacy_singles["score_urn"].isin(params["score_urn"])
        ].copy()
        singles["w"] = pd.to_numeric(singles["w"], errors="coerce")
        singles = singles.loc[np.isfinite(singles["w"])].copy()
    for frame in (scores, singles):
        frame["position"] = pd.to_numeric(frame["position"], errors="coerce")
        frame.drop(columns=["identifier_uniprot"], errors="ignore", inplace=True)
        frame["identifier_uniprot"] = frame["score_urn"].map(
            meta.set_index("score_urn")["identifier_uniprot"]
        )
    LOGGER.info(
        "Loaded %d score sets, %d score rows and %d single substitutions",
        len(meta),
        len(scores),
        len(singles),
    )
    return HomologyInputs(meta, scores, singles, params, paths)


def input_fingerprint(paths: list[Path]) -> dict[str, str]:
    """Content hashes make cache reuse independent of timestamps and file sizes."""
    fingerprints = {}
    for path in sorted(set(p.resolve() for p in paths)):
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        fingerprints[str(path)] = digest.hexdigest()
    return fingerprints
