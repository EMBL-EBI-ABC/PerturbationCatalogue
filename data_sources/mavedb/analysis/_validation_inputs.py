"""Input adapters for native MAVE GMM validation tables."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from data_sources.mavedb.analysis._homology_clinvar import combine_clnsigs
from data_sources.mavedb.analysis._homology_inputs import (
    _native_paths,
    annotate_positions,
    boolean_values,
    read_metadata,
)
from data_sources.mavedb.scripts.process_mavedb import AMINO_ACIDS

LOGGER = logging.getLogger(__name__)
_CHUNK_ROWS = 100_000
_SCORE_REQUIRED = {
    "score_urn",
    "dataset_id",
    "source_row",
    "model",
    "top_model",
    "hgvs_pro",
    "s",
    "w",
    "mode_density",
}
_PARAM_REQUIRED = {
    "score_urn",
    "model",
    "top_model",
    "sigma_neut",
    "median",
    "mad",
}
_SOURCE_REQUIRED = {"dataset_id", "source_row", "hgvs_pro", "wt_mut"}
_CLINVAR_REQUIRED = {
    "uniprotAccessions",
    "clinicalSignificances",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
}
_GWAS_REQUIRED = {
    "studyId",
    "beta",
    "pValueMantissa",
    "pValueExponent",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
}
_KEY_COLUMNS = ("identifier_uniprot", "position", "wt", "mut")


@dataclass
class ValidationInputs:
    """Normalized observations, fitted parameters, reference variants and files."""

    meta: pd.DataFrame
    scores: pd.DataFrame
    params: pd.DataFrame
    clinvar: pd.DataFrame
    gwas: pd.DataFrame
    paths: list[Path]


def _header(path: Path, required: set[str]) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input table not found: {path}")
    try:
        columns = list(
            pd.read_csv(path, sep="\t", nrows=0, keep_default_na=False).columns
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise ValueError(f"{path}: could not read TSV header: {exc}") from exc
    if missing := required - set(columns):
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    return columns


def _blank_to_na(values: pd.Series) -> pd.Series:
    text = values.astype("string")
    return text.mask(text.str.strip().eq(""))


def _numeric(frame: pd.DataFrame, column: str, path: Path) -> pd.Series:
    raw = frame[column]
    result = pd.to_numeric(raw, errors="coerce")
    invalid = raw.notna() & raw.astype("string").str.strip().ne("") & result.isna()
    if invalid.any():
        examples = raw.loc[invalid].astype(str).drop_duplicates().head(5).tolist()
        raise ValueError(f"{path}: invalid numeric value in {column}: {examples}")
    return result


def _read_native_scores(score_paths: list[Path], score_urns: set[str]) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    seen_urns: set[str] = set()
    result_columns: list[str] | None = None
    dtype = {
        "score_urn": "string",
        "dataset_id": "string",
        "source_row": "string",
        "model": "string",
        "top_model": "string",
        "hgvs_pro": "string",
        "hgvs_nt": "string",
        "s": "string",
        "w": "string",
        "mode_density": "string",
    }
    na_values = {column: ["", "NA"] for column in ("top_model", "s", "w")} | {
        "hgvs_pro": [""],
        "hgvs_nt": [""],
    }

    for path in score_paths:
        columns = _header(path, _SCORE_REQUIRED)
        result_columns = columns
        file_urns: set[str] = set()
        try:
            chunks = pd.read_csv(
                path,
                sep="\t",
                dtype={key: value for key, value in dtype.items() if key in columns},
                keep_default_na=False,
                na_values=na_values,
                chunksize=_CHUNK_ROWS,
            )
            for chunk in chunks:
                chunk["score_urn"] = _blank_to_na(chunk["score_urn"])
                file_urns.update(chunk["score_urn"].dropna().astype(str))
                selected = chunk.loc[
                    boolean_values(chunk["top_model"], "top_model")
                    .fillna(False)
                    .to_numpy(dtype=bool)
                ]
                selected = selected.loc[selected["score_urn"].isin(score_urns)]
                if not selected.empty:
                    parts.append(selected)
        except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
            raise ValueError(
                f"{path}: could not read native score table: {exc}"
            ) from exc

        if overlap := seen_urns & file_urns:
            raise ValueError(
                "Overlapping score sets in GMM batches: "
                + ", ".join(sorted(overlap)[:5])
            )
        seen_urns.update(file_urns)

    if parts:
        scores = pd.concat(parts, ignore_index=True)
    else:
        scores = pd.DataFrame(columns=result_columns or sorted(_SCORE_REQUIRED))

    for column in ("score_urn", "dataset_id", "source_row", "model", "hgvs_pro"):
        if column in scores:
            scores[column] = _blank_to_na(scores[column])
    for column in ("s", "w"):
        scores[column] = _numeric(scores, column, score_paths[0])
    return scores


def _read_native_params(param_paths: list[Path], score_urns: set[str]) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    columns: list[str] | None = None
    dtype = {
        "score_urn": "string",
        "model": "string",
        "top_model": "string",
        "sigma_neut": "string",
        "median": "string",
        "mad": "string",
    }
    na_values = {
        column: ["", "NA"] for column in ("top_model", "sigma_neut", "median", "mad")
    }

    for path in param_paths:
        columns = _header(path, _PARAM_REQUIRED)
        try:
            frame = pd.read_csv(
                path,
                sep="\t",
                dtype={key: value for key, value in dtype.items() if key in columns},
                keep_default_na=False,
                na_values=na_values,
                low_memory=False,
            )
        except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
            raise ValueError(
                f"{path}: could not read native parameter table: {exc}"
            ) from exc
        frame["score_urn"] = _blank_to_na(frame["score_urn"])
        selected = frame.loc[
            boolean_values(frame["top_model"], "top_model").fillna(False)
        ].copy()
        selected = selected.loc[selected["score_urn"].isin(score_urns)]
        if not selected.empty:
            parts.append(selected)

    if parts:
        params = pd.concat(parts, ignore_index=True)
    else:
        params = pd.DataFrame(columns=columns or sorted(_PARAM_REQUIRED))

    for column in ("score_urn", "model"):
        if column in params:
            params[column] = _blank_to_na(params[column])
    for column in ("sigma_neut", "median", "mad"):
        params[column] = _numeric(params, column, param_paths[0])
    return params


def _validate_score_identities(scores: pd.DataFrame) -> None:
    if scores.empty:
        return
    identity_columns = ("score_urn", "dataset_id", "source_row", "model", "hgvs_pro")
    if scores[list(identity_columns)].isna().any().any():
        raise ValueError("Selected native scores have missing observation identities")
    if scores["dataset_id"].ne(scores["score_urn"]).any():
        raise ValueError("Native score dataset_id does not match score_urn")
    if not scores["source_row"].str.fullmatch(r"[1-9][0-9]*").all():
        raise ValueError("Native source_row values must be canonical positive integers")
    if scores.duplicated(["dataset_id", "source_row"]).any():
        raise ValueError(
            "Selected native scores contain duplicate observation identities"
        )


def _recover_wt_mut(source_scores: Path, scores: pd.DataFrame) -> pd.Series:
    """Stream source flags and verify each selected row against its source HGVS."""
    columns = _header(source_scores, _SOURCE_REQUIRED)
    optional_score_urn = "score_urn" in columns
    optional_hgvs_nt = "hgvs_nt" in columns
    usecols = ["dataset_id", "source_row", "hgvs_pro", "wt_mut"]
    if optional_score_urn:
        usecols.append("score_urn")
    optional_score = "score" in columns
    if optional_score:
        usecols.append("score")
    if optional_hgvs_nt and "hgvs_nt" in scores:
        usecols.append("hgvs_nt")

    if scores.empty:
        return pd.Series(pd.array([], dtype="boolean"), index=scores.index)

    wanted = scores[
        ["dataset_id", "source_row", "hgvs_pro", "s"]
        + (["hgvs_nt"] if optional_hgvs_nt and "hgvs_nt" in scores else [])
    ].copy()
    wanted_keys = list(zip(wanted["dataset_id"], wanted["source_row"], strict=True))
    wanted_key_set = set(wanted_keys)
    wanted_ids = set(wanted["dataset_id"].astype(str))
    if optional_hgvs_nt and "hgvs_nt" in wanted:
        wanted["hgvs_nt"] = _blank_to_na(wanted["hgvs_nt"])
    wanted = wanted.rename(columns={"hgvs_pro": "native_hgvs_pro"})
    wanted = wanted.rename(columns={"s": "native_s"})
    if "hgvs_nt" in wanted:
        wanted = wanted.rename(columns={"hgvs_nt": "native_hgvs_nt"})

    matched_flags: dict[tuple[str, str], bool] = {}
    matched_keys: set[tuple[str, str]] = set()
    dtype = {
        "dataset_id": "string",
        "source_row": "string",
        "hgvs_pro": "string",
        "wt_mut": "string",
        "score_urn": "string",
        "hgvs_nt": "string",
        "score": "string",
    }
    na_values = {column: [""] for column in ("hgvs_pro", "hgvs_nt", "wt_mut")} | {
        "score": ["", "NA"]
    }

    try:
        chunks = pd.read_csv(
            source_scores,
            sep="\t",
            usecols=usecols,
            dtype={key: value for key, value in dtype.items() if key in usecols},
            keep_default_na=False,
            na_values=na_values,
            chunksize=_CHUNK_ROWS,
        )
        for chunk in chunks:
            chunk = chunk.loc[chunk["dataset_id"].isin(wanted_ids)].copy()
            if chunk.empty:
                continue
            chunk["dataset_id"] = _blank_to_na(chunk["dataset_id"])
            chunk["source_row"] = _blank_to_na(chunk["source_row"])
            chunk = chunk.loc[chunk["source_row"].notna()]
            relevant_keys = list(
                zip(chunk["dataset_id"], chunk["source_row"], strict=True)
            )
            chunk = chunk.loc[[key in wanted_key_set for key in relevant_keys]].copy()
            if chunk.empty:
                continue
            relevant_keys = list(
                zip(chunk["dataset_id"], chunk["source_row"], strict=True)
            )
            if chunk.duplicated(["dataset_id", "source_row"]).any() or any(
                key in matched_keys for key in relevant_keys
            ):
                raise ValueError(
                    f"{source_scores}: duplicate selected (dataset_id, source_row) identity"
                )

            joined = chunk.merge(
                wanted,
                on=["dataset_id", "source_row"],
                how="inner",
                validate="one_to_one",
                sort=False,
            )
            source_hgvs = _blank_to_na(joined["hgvs_pro"])
            native_hgvs = _blank_to_na(joined["native_hgvs_pro"])
            mismatch = (
                source_hgvs.isna() | native_hgvs.isna() | source_hgvs.ne(native_hgvs)
            )
            if mismatch.any():
                examples = joined.loc[mismatch, ["dataset_id", "source_row"]].head(5)
                raise ValueError(
                    f"{source_scores}: source HGVS does not match native score identities: "
                    f"{examples.to_dict(orient='records')}"
                )
            if optional_score_urn:
                source_urn = _blank_to_na(joined["score_urn"])
                if (
                    source_urn.notna().any()
                    and source_urn.fillna(joined["dataset_id"])
                    .ne(joined["dataset_id"])
                    .any()
                ):
                    raise ValueError(
                        f"{source_scores}: score_urn does not match dataset_id"
                    )
            if "native_hgvs_nt" in joined:
                source_hgvs_nt = _blank_to_na(joined["hgvs_nt"])
                native_hgvs_nt = _blank_to_na(joined["native_hgvs_nt"])
                both_present = source_hgvs_nt.notna() & native_hgvs_nt.notna()
                mismatch_nt = both_present & source_hgvs_nt.ne(native_hgvs_nt)
                if mismatch_nt.any():
                    raise ValueError(
                        f"{source_scores}: source hgvs_nt does not match native score identities"
                    )
            if optional_score:
                source_score = _numeric(
                    joined.rename(columns={"score": "source_score"}),
                    "source_score",
                    source_scores,
                )
                native_score = pd.to_numeric(joined["native_s"], errors="coerce")
                same_score = np.isclose(
                    source_score.to_numpy(dtype=float),
                    native_score.to_numpy(dtype=float),
                    rtol=1e-12,
                    atol=1e-12,
                    equal_nan=True,
                )
                if not same_score.all():
                    raise ValueError(
                        f"{source_scores}: source score does not match native score identities"
                    )

            flags = boolean_values(joined["wt_mut"], "source wt_mut")
            if flags.isna().any():
                raise ValueError(
                    f"{source_scores}: selected source rows have missing wt_mut"
                )
            for key, flag in zip(relevant_keys, flags.astype(bool), strict=True):
                matched_flags[key] = bool(flag)
                matched_keys.add(key)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        if isinstance(exc, ValueError) and str(exc).startswith(str(source_scores)):
            raise
        raise ValueError(
            f"{source_scores}: could not validate source scores: {exc}"
        ) from exc

    unmatched = wanted_key_set - matched_keys
    if unmatched:
        raise ValueError(
            f"{source_scores}: missing selected score identities: "
            f"{sorted(unmatched)[:5]}"
        )
    return pd.Series(
        pd.array([matched_flags[key] for key in wanted_keys], dtype="boolean"),
        index=scores.index,
    )


def _read_reference(
    path: Path,
    required: set[str],
    *,
    identifiers: tuple[str, ...],
    text_columns: tuple[str, ...],
) -> pd.DataFrame:
    columns = _header(path, required)
    dtype = {column: "string" for column in (*identifiers, *text_columns)}
    dtype = {key: value for key, value in dtype.items() if key in columns}
    try:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=dtype,
            keep_default_na=False,
            na_values={column: ["", "NA"] for column in columns},
            low_memory=False,
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise ValueError(f"{path}: could not read reference TSV: {exc}") from exc


def _normalize_clinvar(path: Path) -> pd.DataFrame:
    raw = _read_reference(
        path,
        _CLINVAR_REQUIRED,
        identifiers=("uniprotAccessions",),
        text_columns=(
            "clinicalSignificances",
            "referenceAminoAcid",
            "alternateAminoAcid",
        ),
    )
    clinvar = raw.loc[:, sorted(_CLINVAR_REQUIRED)].rename(
        columns={
            "uniprotAccessions": "identifier_uniprot",
            "clinicalSignificances": "clnsig",
            "aminoAcidPosition": "position",
            "referenceAminoAcid": "wt",
            "alternateAminoAcid": "mut",
        }
    )
    clinvar["identifier_uniprot"] = _blank_to_na(clinvar["identifier_uniprot"])
    clinvar["wt"] = _blank_to_na(clinvar["wt"])
    clinvar["mut"] = _blank_to_na(clinvar["mut"])
    clinvar["position"] = _numeric(clinvar, "position", path)
    accepted = clinvar.loc[
        clinvar["wt"].isin(AMINO_ACIDS) & clinvar["mut"].isin(AMINO_ACIDS)
    ].dropna(subset=list(_KEY_COLUMNS))
    if accepted.empty:
        return pd.DataFrame(columns=[*_KEY_COLUMNS, "clnsig"])
    return (
        accepted.groupby(list(_KEY_COLUMNS), sort=False, dropna=False)["clnsig"]
        .agg(combine_clnsigs)
        .reset_index()
    )


def _normalize_gwas(path: Path) -> pd.DataFrame:
    raw = _read_reference(
        path,
        _GWAS_REQUIRED,
        identifiers=("studyId", "uniprotAccessions"),
        text_columns=("referenceAminoAcid", "alternateAminoAcid"),
    )
    gwas = raw.loc[:, sorted(_GWAS_REQUIRED)].rename(
        columns={
            "studyId": "study_id",
            "pValueMantissa": "pvalue_mantissa",
            "pValueExponent": "pvalue_exponent",
            "uniprotAccessions": "identifier_uniprot",
            "aminoAcidPosition": "position",
            "referenceAminoAcid": "wt",
            "alternateAminoAcid": "mut",
        }
    )
    gwas["study_id"] = _blank_to_na(gwas["study_id"])
    gwas["identifier_uniprot"] = _blank_to_na(gwas["identifier_uniprot"])
    gwas["wt"] = _blank_to_na(gwas["wt"])
    gwas["mut"] = _blank_to_na(gwas["mut"])
    for column in ("position", "beta", "pvalue_mantissa", "pvalue_exponent"):
        gwas[column] = _numeric(gwas, column, path)
    return gwas.dropna(subset=list(_KEY_COLUMNS)).reset_index(drop=True)


def _validate_params(params: pd.DataFrame) -> None:
    if params.empty:
        return
    if params["score_urn"].isna().any() or params["model"].isna().any():
        raise ValueError("Selected native parameters have missing score_urn or model")
    if params["score_urn"].duplicated().any():
        raise ValueError("More than one selected native model per score set")
    for column in ("sigma_neut", "median", "mad"):
        values = params[column].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Selected native parameter {column} must be finite")
    if not params["sigma_neut"].gt(0).all():
        raise ValueError("Selected native parameter sigma_neut must be positive")
    if not params["mad"].gt(0).all():
        raise ValueError("Selected native parameter mad must be positive")


def load_inputs(
    *,
    metadata: Path,
    gmm_dir: Path,
    source_scores: Path,
    clinvar: Path,
    gwas: Path,
) -> ValidationInputs:
    """Read canonical inputs and selected native-model observations for validation."""
    meta = read_metadata(metadata)
    score_paths, param_paths = _native_paths(gmm_dir)
    score_urns = set(meta["score_urn"].astype(str))
    scores = _read_native_scores(score_paths, score_urns)
    _validate_score_identities(scores)
    scores["wt_mut"] = _recover_wt_mut(source_scores, scores)
    scores = annotate_positions(scores, meta)
    scores["synonymous_mut"] = scores["wt"].eq(scores["mut"]).fillna(False).astype(bool)
    scores["nonsense_mut"] = scores["mut"].eq("stop").fillna(False).astype(bool)
    scores["identifier_uniprot"] = scores["score_urn"].map(
        meta.set_index("score_urn")["identifier_uniprot"]
    )

    params = _read_native_params(param_paths, score_urns)
    _validate_params(params)
    if not scores.empty:
        missing_params = set(scores["score_urn"]) - set(params["score_urn"])
        if missing_params:
            raise ValueError(
                "Selected native scores without selected parameters: "
                + ", ".join(sorted(missing_params)[:5])
            )
        expected_models = scores["score_urn"].map(
            params.set_index("score_urn")["model"]
        )
        if scores["model"].ne(expected_models).any():
            raise ValueError("Selected native score and parameter models do not agree")

    clinvar_frame = _normalize_clinvar(clinvar)
    gwas_frame = _normalize_gwas(gwas)
    paths = [metadata, source_scores, *score_paths, *param_paths, clinvar, gwas]
    LOGGER.info(
        "Loaded %d selected native score rows from %d score sets, %d ClinVar variants, "
        "and %d GWAS variants",
        len(scores),
        scores["score_urn"].nunique(),
        len(clinvar_frame),
        len(gwas_frame),
    )
    return ValidationInputs(meta, scores, params, clinvar_frame, gwas_frame, paths)
