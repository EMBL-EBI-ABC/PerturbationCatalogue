"""Deterministic score validation tables for fitted MAVE GMM scores."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from data_sources.mavedb.analysis._homology_clinvar import calc_roc, pr_auc

_CONTROL_SCORE_TYPES = ("s", "ns", "w", "p")
_VARIANT_SCORE_TYPES = ("s", "ns", "w", "aw", "p", "ap")
_ROC_COLUMNS = (
    "thresh",
    "tp",
    "tn",
    "fp",
    "fn",
    "tpr",
    "tnr",
    "fpr",
    "precision",
    "auc",
)
_ROC_OUTPUT_COLUMNS = (
    "benchmark",
    "score_type",
    "score_urn",
    *_ROC_COLUMNS,
    "pr_auc",
)
_SUMMARY_COLUMNS = (
    "benchmark",
    "score_type",
    "score_urn",
    "auc",
    "pr_auc",
    "n_positive",
    "n_negative",
    "n_observations",
)
_COUNTS_COLUMNS = (
    "benchmark",
    "score_urn",
    "n_observations",
    "n_positive",
    "n_negative",
)
_VARIANT_KEYS = ("identifier_uniprot", "position", "wt", "mut")
_SCORE_REQUIRED = (
    "score_urn",
    "s",
    "w",
    "wt_mut",
    "synonymous_mut",
    "nonsense_mut",
    *_VARIANT_KEYS,
)
_PARAM_REQUIRED = ("score_urn", "sigma_neut", "median", "mad")
_META_REQUIRED = ("score_urn", "taxid")
_CLINVAR_REQUIRED = (*_VARIANT_KEYS, "clnsig")
_GWAS_REQUIRED = (
    *_VARIANT_KEYS,
    "study_id",
    "beta",
    "pvalue_mantissa",
    "pvalue_exponent",
)


def _prepare_input(
    frame: pd.DataFrame, required: Sequence[str], name: str
) -> pd.DataFrame:
    """Check non-empty input schemas while allowing schema-less empty inputs."""
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame")
    missing = [column for column in required if column not in frame.columns]
    if missing and not frame.empty:
        raise ValueError(f"{name} is missing required columns: {missing}")
    if missing:
        frame = frame.copy()
        for column in missing:
            frame[column] = pd.Series(dtype="object")
    return frame.copy()


def _check_unique_keys(frame: pd.DataFrame, column: str, name: str) -> None:
    present = frame[column].notna()
    if frame.loc[present, column].duplicated().any():
        raise ValueError(f"{name} must contain at most one row per {column}")


def _safe_left_merge(
    left: pd.DataFrame, right: pd.DataFrame, keys: Sequence[str]
) -> pd.DataFrame:
    """Left-join without allowing null join keys to match each other."""
    valid_right = right.loc[right.loc[:, list(keys)].notna().all(axis=1)].copy()
    return left.merge(valid_right, how="left", on=list(keys), sort=False)


def _numeric(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").astype(float)


def _transform_scores(scores: pd.DataFrame, params: pd.DataFrame) -> pd.DataFrame:
    parameter_values = params.loc[:, list(_PARAM_REQUIRED)].copy()
    parameter_values = parameter_values.loc[parameter_values["score_urn"].notna()]
    joined = _safe_left_merge(scores, parameter_values, ("score_urn",))

    s = _numeric(joined["s"])
    w = _numeric(joined["w"])
    median = _numeric(joined["median"])
    mad = _numeric(joined["mad"])
    sigma_neut = _numeric(joined["sigma_neut"])
    with np.errstate(divide="ignore", invalid="ignore"):
        joined["ns"] = (s - median) / mad
        joined["aw"] = (w - 1.0).abs()
        joined["p"] = w / sigma_neut
        joined["ap"] = (w - 1.0).abs() / sigma_neut
    return joined


def _control_labels(frame: pd.DataFrame) -> pd.Series:
    wt_mut = frame["wt_mut"].fillna(False).astype(bool)
    synonymous_mut = frame["synonymous_mut"].fillna(False).astype(bool)
    nonsense_mut = frame["nonsense_mut"].fillna(False).astype(bool)
    is_negative = wt_mut | synonymous_mut
    labels = pd.Series(pd.NA, index=frame.index, dtype="boolean")
    labels.loc[is_negative] = False
    labels.loc[~is_negative & nonsense_mut] = True
    return labels


def _clinvar_labels(frame: pd.DataFrame) -> pd.Series:
    labels = frame["clnsig"].map(
        {
            "benign": False,
            "likely benign": False,
            "pathogenic": True,
            "likely pathogenic": True,
        }
    )
    return labels.astype("boolean")


def _matched_clinvar(scores: pd.DataFrame, clinvar: pd.DataFrame) -> pd.DataFrame:
    matchable = scores.loc[
        scores.loc[:, list(_VARIANT_KEYS)].notna().all(axis=1)
    ].copy()
    matches = matchable.merge(
        clinvar.loc[clinvar.loc[:, list(_VARIANT_KEYS)].notna().all(axis=1)].copy(),
        how="inner",
        on=list(_VARIANT_KEYS),
        sort=False,
    )
    matches["cat"] = _clinvar_labels(matches)
    return matches


def _matched_gwas(
    scores: pd.DataFrame, meta: pd.DataFrame, gwas: pd.DataFrame
) -> pd.DataFrame:
    human_urns = set(
        meta.loc[pd.to_numeric(meta["taxid"], errors="coerce").eq(9606), "score_urn"]
        .dropna()
        .tolist()
    )
    human_scores = scores.loc[scores["score_urn"].isin(human_urns)].copy()
    human_scores = human_scores.loc[
        human_scores.loc[:, list(_VARIANT_KEYS)].notna().all(axis=1)
    ].copy()
    valid_gwas = gwas.loc[gwas.loc[:, list(_VARIANT_KEYS)].notna().all(axis=1)].copy()
    matches = human_scores.merge(
        valid_gwas,
        how="left",
        on=list(_VARIANT_KEYS),
        sort=False,
    )
    matches["gwas_hit"] = matches["study_id"].notna().astype(bool)
    return matches


def _empty_roc_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "benchmark": pd.Series(dtype="object"),
            "score_type": pd.Series(dtype="object"),
            "score_urn": pd.Series(dtype="object"),
            **{column: pd.Series(dtype="float64") for column in _ROC_COLUMNS},
            "pr_auc": pd.Series(dtype="float64"),
        }
    )


def _empty_summary_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            **{column: pd.Series(dtype="object") for column in _SUMMARY_COLUMNS[:3]},
            **{column: pd.Series(dtype="float64") for column in ("auc", "pr_auc")},
            **{column: pd.Series(dtype="int64") for column in _SUMMARY_COLUMNS[-3:]},
        }
    )


def _empty_counts_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "benchmark": pd.Series(dtype="object"),
            "score_urn": pd.Series(dtype="object"),
            **{column: pd.Series(dtype="int64") for column in _COUNTS_COLUMNS[-3:]},
        }
    )


def _group_rows(
    frame: pd.DataFrame,
    score_types: Sequence[str],
    label_column: str,
    benchmark: str,
    greater_types: set[str],
    max_steps_all: int,
    max_steps_per_study: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    roc_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    counts_rows: list[dict[str, Any]] = []

    if frame.empty:
        return _empty_roc_frame(), _empty_summary_frame(), _empty_counts_frame()

    studies = [
        (score_urn, group.copy())
        for score_urn, group in frame.groupby(
            "score_urn", sort=False, dropna=False, observed=True
        )
    ]
    groups = [("all", frame.copy()), *studies]

    for score_urn, group in groups:
        labels = group[label_column]
        labels_valid = labels.notna()
        true = labels.loc[labels_valid].astype(bool)
        counts_rows.append(
            {
                "benchmark": benchmark,
                "score_urn": score_urn,
                "n_observations": int(len(true)),
                "n_positive": int(true.sum()),
                "n_negative": int((~true).sum()),
            }
        )

    for score_type in score_types:
        for score_urn, group in groups:
            labels = group[label_column]
            numeric_score = _numeric(group[score_type])
            valid = labels.notna() & numeric_score.notna()
            truth = labels.loc[valid].astype(bool)
            score = numeric_score.loc[valid]
            positive_count = int(truth.sum())
            negative_count = int((~truth).sum())
            max_steps = max_steps_all if score_urn == "all" else max_steps_per_study
            curve = calc_roc(
                truth,
                score,
                greater=score_type in greater_types,
                max_steps=max_steps,
            )
            if curve.empty:
                auc_value = float("nan")
                pr_auc_value = float("nan")
            else:
                pr_auc_value = pr_auc(curve["tpr"], curve["precision"])
                curve = curve.copy()
                curve["pr_auc"] = pr_auc_value
                curve.insert(0, "score_urn", score_urn)
                curve.insert(0, "score_type", score_type)
                curve.insert(0, "benchmark", benchmark)
                roc_frames.append(curve.loc[:, list(_ROC_OUTPUT_COLUMNS)])
                auc_value = float(curve["auc"].iloc[0])

            summary_rows.append(
                {
                    "benchmark": benchmark,
                    "score_type": score_type,
                    "score_urn": score_urn,
                    "auc": auc_value,
                    "pr_auc": pr_auc_value,
                    "n_positive": positive_count,
                    "n_negative": negative_count,
                    "n_observations": int(len(truth)),
                }
            )

    roc = pd.concat(roc_frames, ignore_index=True) if roc_frames else _empty_roc_frame()
    summary = (
        pd.DataFrame.from_records(summary_rows, columns=list(_SUMMARY_COLUMNS))
        if summary_rows
        else _empty_summary_frame()
    )
    counts = (
        pd.DataFrame.from_records(counts_rows, columns=list(_COUNTS_COLUMNS))
        if counts_rows
        else _empty_counts_frame()
    )
    return roc, summary, counts


def evaluate_validation(
    *,
    scores: pd.DataFrame,
    params: pd.DataFrame,
    meta: pd.DataFrame,
    clinvar: pd.DataFrame,
    gwas: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Build control, ClinVar, and human GWAS validation tables.

    ``scores`` contains the selected individual observations; ``params`` has
    one selected parameter row per ``score_urn``. Variant joins preserve
    duplicate matches and never match missing keys. All returned tables are
    deterministic and retain source observation columns where applicable.
    """
    score_input = _prepare_input(scores, _SCORE_REQUIRED, "scores")
    parameter_input = _prepare_input(params, _PARAM_REQUIRED, "params")
    meta_input = _prepare_input(meta, _META_REQUIRED, "meta")
    clinvar_input = _prepare_input(clinvar, _CLINVAR_REQUIRED, "clinvar")
    gwas_input = _prepare_input(gwas, _GWAS_REQUIRED, "gwas")
    _check_unique_keys(parameter_input, "score_urn", "params")
    _check_unique_keys(meta_input, "score_urn", "meta")

    score_table = _transform_scores(score_input, parameter_input)

    controls = score_table.copy()
    controls["cat"] = _control_labels(controls)
    controls = controls.loc[controls["cat"].notna()].copy()

    clinvar_matches = _matched_clinvar(score_table, clinvar_input)
    gwas_matches = _matched_gwas(score_table, meta_input, gwas_input)

    control_roc, control_summary, control_counts = _group_rows(
        controls,
        _CONTROL_SCORE_TYPES,
        "cat",
        "controls",
        set(),
        4000,
        1000,
    )
    clinvar_binary = clinvar_matches.loc[clinvar_matches["cat"].notna()].copy()
    clinvar_roc, clinvar_summary, clinvar_counts = _group_rows(
        clinvar_binary,
        _VARIANT_SCORE_TYPES,
        "cat",
        "clinvar",
        {"aw", "ap"},
        6000,
        1000,
    )

    if gwas_matches.empty:
        gwas_for_roc = gwas_matches.copy()
    else:
        hit_by_urn = gwas_matches.groupby(
            "score_urn", sort=False, dropna=False, observed=True
        )["gwas_hit"].transform("any")
        gwas_for_roc = gwas_matches.loc[hit_by_urn].copy()
    gwas_roc, gwas_summary, gwas_counts = _group_rows(
        gwas_for_roc,
        _VARIANT_SCORE_TYPES,
        "gwas_hit",
        "gwas",
        {"aw", "ap"},
        2000,
        500,
    )

    summary_frames = [
        frame
        for frame in (control_summary, clinvar_summary, gwas_summary)
        if not frame.empty
    ]
    counts_frames = [
        frame
        for frame in (control_counts, clinvar_counts, gwas_counts)
        if not frame.empty
    ]

    return {
        "scores": score_table,
        "controls": controls,
        "clinvar_matches": clinvar_matches,
        "gwas_matches": gwas_matches,
        "ctrl_roc": control_roc,
        "clinvar_roc": clinvar_roc,
        "gwas_roc": gwas_roc,
        "auc_summary": (
            pd.concat(summary_frames, ignore_index=True)
            if summary_frames
            else _empty_summary_frame()
        ),
        "counts": (
            pd.concat(counts_frames, ignore_index=True)
            if counts_frames
            else _empty_counts_frame()
        ),
    }
