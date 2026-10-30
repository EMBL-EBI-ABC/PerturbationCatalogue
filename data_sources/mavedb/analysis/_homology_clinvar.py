"""Build ClinVar comparisons for the MAVE homology analysis."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import pandas as pd

from data_sources.mavedb.analysis._mave_homology import HOMOLOGY_TYPES

_AA_SORTED = tuple("ACDEFGHIKLMNPQRSTVWY")
_RECOGNISED_CLNSIG = {
    "pathogenic",
    "likely pathogenic",
    "benign",
    "likely benign",
    "conflicting classifications of pathogenicity",
}
_SCORE_TYPES = ("w", "aw", "p", "ap")
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
_COUNT_COLUMNS = ("type", "same_wt", "n_variants", "n_studies")
_ROC_GROUP_COLUMNS = ("score_type", "type", "same_wt", "homolog_urn")
_SUMMARY_COLUMNS = ("type", "same_wt", "auc", "pr_auc", "n_variants", "n_studies")


def combine_clnsigs(x: Iterable[Any] | Any) -> str:
    """Combine ClinVar labels using the precedence in ``mave_homology.R``."""
    if isinstance(x, str):
        values = (x,)
    elif isinstance(x, Iterable):
        values = x
    else:
        values = (x,)

    labels = {
        value
        for value in values
        if isinstance(value, str) and value in _RECOGNISED_CLNSIG
    }
    if not labels:
        return "uncertain significance"
    if "conflicting classifications of pathogenicity" in labels:
        return "conflicting classifications of pathogenicity"
    if labels & {"pathogenic", "likely pathogenic"} and labels & {
        "benign",
        "likely benign",
    }:
        return "conflicting classifications of pathogenicity"
    if "pathogenic" in labels:
        return "pathogenic"
    if "likely pathogenic" in labels:
        return "likely pathogenic"
    if "benign" in labels:
        return "benign"
    return "likely benign"


def _empty_roc() -> pd.DataFrame:
    return pd.DataFrame({column: pd.Series(dtype="float64") for column in _ROC_COLUMNS})


def _step_auc(x: Sequence[float], y: Sequence[float]) -> float:
    """Return the left-step area used by the source R analysis."""
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if len(x_values) != len(y_values):
        raise ValueError("x and y must have the same length")
    if len(x_values) < 2:
        return 0.0
    return float(np.sum(np.diff(x_values) * y_values[:-1]))


def pr_auc(recall: Sequence[float], precision: Sequence[float]) -> float:
    """Compute the R-compatible left-step area under a precision-recall curve."""
    recall_values = np.asarray(recall, dtype=float)
    precision_values = np.asarray(precision, dtype=float)
    if len(recall_values) != len(precision_values):
        raise ValueError("recall and precision must have the same length")

    valid = ~np.isnan(recall_values) & ~np.isnan(precision_values)
    recall_values = recall_values[valid]
    precision_values = precision_values[valid]
    if len(recall_values) == 0:
        return float("nan")
    return _step_auc(
        np.concatenate(([0.0], recall_values)),
        np.concatenate(([precision_values[0]], precision_values)),
    )


def calc_roc(
    true: Sequence[Any] | pd.Series,
    score: Sequence[float] | pd.Series,
    greater: bool = True,
    max_steps: int = 6000,
) -> pd.DataFrame:
    """Calculate threshold counts and ROC metrics with source-R semantics.

    A threshold classifies observations as positive when ``score >= thresh``
    if ``greater`` is true, and ``score <= thresh`` otherwise. Tied values are
    evaluated together. The AUC uses the R helper's left-step rule rather than
    trapezoidal integration or average precision.
    """
    if max_steps <= 0:
        raise ValueError("max_steps must be positive")

    truth = pd.Series(true).reset_index(drop=True)
    scores = pd.to_numeric(pd.Series(score), errors="coerce").reset_index(drop=True)
    if len(truth) != len(scores):
        raise ValueError("true and score must have the same length")
    valid = truth.notna() & scores.notna()
    truth = truth.loc[valid].astype(bool).to_numpy()
    scores = scores.loc[valid].to_numpy(dtype=float)
    if len(scores) == 0:
        return _empty_roc()

    unique_scores = np.unique(scores)
    if len(unique_scores) > max_steps:
        stride = max(1, len(unique_scores) // max_steps)
        # R's seq(0, n, by=stride) is used as one-based vector indices;
        # index zero is ignored by R, leaving stride, 2*stride, ... .
        selected = np.arange(stride, len(unique_scores) + 1, stride) - 1
        thresholds = np.concatenate(([-np.inf], unique_scores[selected], [np.inf]))
    else:
        thresholds = np.concatenate(([-np.inf], unique_scores, [np.inf]))

    order = np.argsort(scores, kind="stable")
    sorted_scores = scores[order]
    sorted_truth = truth[order]
    positive_prefix = np.concatenate(([0], np.cumsum(sorted_truth, dtype=np.int64)))
    negative_prefix = np.concatenate(([0], np.cumsum(~sorted_truth, dtype=np.int64)))
    total_positive = int(positive_prefix[-1])
    total_negative = int(negative_prefix[-1])
    if greater:
        prefix_index = np.searchsorted(sorted_scores, thresholds, side="left")
        tp = total_positive - positive_prefix[prefix_index]
        fp = total_negative - negative_prefix[prefix_index]
        fn = total_positive - tp
        tn = total_negative - fp
    else:
        prefix_index = np.searchsorted(sorted_scores, thresholds, side="right")
        tp = positive_prefix[prefix_index]
        fp = negative_prefix[prefix_index]
        fn = total_positive - tp
        tn = total_negative - fp

    def divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
        result = np.full(denominator.shape, np.nan, dtype=float)
        np.divide(numerator, denominator, out=result, where=denominator != 0)
        return result

    curve = pd.DataFrame(
        {
            "thresh": thresholds,
            "tp": tp,
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "tpr": divide(tp, tp + fn),
            "tnr": divide(tn, tn + fp),
            "fpr": divide(fp, tn + fp),
            "precision": divide(tp, tp + fp),
        }
    )
    curve = curve.sort_values(
        "thresh", ascending=not greater, kind="stable"
    ).reset_index(drop=True)
    curve["auc"] = _step_auc(curve["fpr"], curve["tpr"])
    return curve.loc[:, list(_ROC_COLUMNS)]


def _empty_homologs(homology_columns: Sequence[str] = ()) -> pd.DataFrame:
    columns = [column for column in homology_columns if column != "taxonomy_level"]
    columns.extend(
        (
            "uniprot",
            "homolog_uniprot",
            "homolog_urn",
            "position",
            "homolog_position",
            "wt",
            "homolog_wt",
            "same_wt",
            "mut",
            "w",
            "clnsig",
        )
    )
    return pd.DataFrame(columns=list(dict.fromkeys(columns)))


def _empty_counts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "type": pd.Series(dtype="object"),
            "same_wt": pd.Series(dtype="object"),
            "n_variants": pd.Series(dtype="int64"),
            "n_studies": pd.Series(dtype="int64"),
        }
    )


def _empty_clinvar_roc() -> pd.DataFrame:
    return pd.DataFrame(
        {
            **{column: pd.Series(dtype="object") for column in _ROC_GROUP_COLUMNS},
            **{column: pd.Series(dtype="float64") for column in _ROC_COLUMNS},
            "pr_auc": pd.Series(dtype="float64"),
            "n_variants": pd.Series(dtype="int64"),
            "n_studies": pd.Series(dtype="int64"),
        }
    )


def _label_clinvar_variants(clinvar_variants: pd.DataFrame) -> pd.DataFrame:
    required = (
        "uniprotAccessions",
        "clinicalSignificances",
        "aminoAcidPosition",
        "referenceAminoAcid",
        "alternateAminoAcid",
    )
    missing = [column for column in required if column not in clinvar_variants]
    if missing:
        raise ValueError(f"clinvar_variants is missing required columns: {missing}")

    accepted = clinvar_variants.loc[
        clinvar_variants["referenceAminoAcid"].isin(_AA_SORTED)
        & clinvar_variants["alternateAminoAcid"].isin(_AA_SORTED),
        list(required),
    ].copy()
    accepted = accepted.rename(
        columns={
            "uniprotAccessions": "identifier_uniprot",
            "clinicalSignificances": "clnsig",
            "aminoAcidPosition": "position",
            "referenceAminoAcid": "wt",
            "alternateAminoAcid": "mut",
        }
    )
    grouped = accepted.groupby(
        ["identifier_uniprot", "position", "wt", "mut"], sort=False, dropna=False
    )["clnsig"].agg(combine_clnsigs)
    return grouped.reset_index()


def _build_clinvar_homologs(
    combined: pd.DataFrame,
    scores: pd.DataFrame,
    homology_key: pd.DataFrame,
) -> pd.DataFrame:
    required_key = (
        "target_uniprot",
        "source_uniprot",
        "target_position",
        "source_position",
        "target_wt",
        "source_wt",
        "type",
    )
    missing_key = [column for column in required_key if column not in homology_key]
    if missing_key:
        raise ValueError(f"homology_key is missing required columns: {missing_key}")
    required_scores = ("score_urn", "identifier_uniprot", "w", "position", "wt", "mut")
    missing_scores = [column for column in required_scores if column not in scores]
    if missing_scores:
        raise ValueError(f"scores is missing required columns: {missing_scores}")

    key = homology_key.loc[
        homology_key["target_uniprot"].notna()
        & homology_key["source_uniprot"].notna()
        & homology_key["target_uniprot"].ne(homology_key["source_uniprot"])
    ].copy()
    if combined.empty or key.empty:
        return _empty_homologs(homology_key.columns)

    joined = combined.merge(
        key,
        how="left",
        left_on=["identifier_uniprot", "position", "wt"],
        right_on=["target_uniprot", "target_position", "target_wt"],
        sort=False,
    )
    score_lookup = scores.loc[:, list(required_scores)].rename(
        columns={
            "identifier_uniprot": "_score_uniprot",
            "position": "_score_position",
            "wt": "_score_wt",
            "mut": "_score_mut",
            "w": "_score_w",
        }
    )
    score_join_columns = [
        "_score_uniprot",
        "_score_position",
        "_score_wt",
        "_score_mut",
    ]
    score_lookup = score_lookup.dropna(subset=score_join_columns)
    joined = joined.dropna(
        subset=["source_uniprot", "source_position", "source_wt", "mut"]
    )
    joined = joined.merge(
        score_lookup,
        how="inner",
        left_on=["source_uniprot", "source_position", "source_wt", "mut"],
        right_on=["_score_uniprot", "_score_position", "_score_wt", "_score_mut"],
        sort=False,
    )
    if joined.empty:
        return _empty_homologs(homology_key.columns)

    joined["type"] = joined["type"].map(HOMOLOGY_TYPES)
    same_wt = joined["wt"].eq(joined["source_wt"])
    joined["same_wt"] = np.where(same_wt, "Same WT", "Diff. WT")
    joined.loc[joined["wt"].isna() | joined["source_wt"].isna(), "same_wt"] = pd.NA

    metadata_columns = [
        column for column in homology_key.columns if column != "taxonomy_level"
    ]
    result = joined.loc[:, metadata_columns].copy()
    result["uniprot"] = joined["identifier_uniprot"]
    result["homolog_uniprot"] = joined["source_uniprot"]
    result["homolog_urn"] = joined["score_urn"]
    result["position"] = joined["position"]
    result["homolog_position"] = joined["source_position"]
    result["wt"] = joined["wt"]
    result["homolog_wt"] = joined["source_wt"]
    result["same_wt"] = joined["same_wt"]
    result["mut"] = joined["mut"]
    result["w"] = pd.to_numeric(joined["_score_w"], errors="coerce")
    result["clnsig"] = joined["clnsig"]
    return result.reset_index(drop=True)


def _expanded_count_rows(homologs: pd.DataFrame) -> pd.DataFrame:
    columns = ["type", "homolog_urn", "same_wt"]
    base = homologs.loc[:, columns].copy()
    expansions = (
        base,
        base.assign(same_wt="Both"),
        base.assign(type="All"),
        base.assign(type="All", same_wt="Both"),
    )
    return pd.concat(expansions, ignore_index=True)


def _clinvar_counts(homologs: pd.DataFrame) -> pd.DataFrame:
    if homologs.empty:
        return _empty_counts()
    expanded = _expanded_count_rows(homologs)
    return (
        expanded.groupby(["type", "same_wt"], sort=False, dropna=False)
        .agg(
            n_variants=("homolog_urn", "size"),
            n_studies=("homolog_urn", lambda values: values.nunique(dropna=False)),
        )
        .reset_index()
        .loc[:, _COUNT_COLUMNS]
    )


def _expanded_roc_rows(binary: pd.DataFrame) -> pd.DataFrame:
    expanded_rows: list[pd.DataFrame] = []
    for type_value, type_group in ((None, binary), ("All", binary)):
        for same_wt_value, same_wt_group in ((None, type_group), ("Both", type_group)):
            for urn_value, group in ((None, same_wt_group), ("All", same_wt_group)):
                expanded = group.copy()
                if type_value is not None:
                    expanded["type"] = type_value
                if same_wt_value is not None:
                    expanded["same_wt"] = same_wt_value
                if urn_value is not None:
                    expanded["homolog_urn"] = urn_value
                expanded_rows.append(expanded)
    return pd.concat(expanded_rows, ignore_index=True)


def _clinvar_roc(
    homologs: pd.DataFrame,
    params: pd.DataFrame,
    counts: pd.DataFrame,
) -> pd.DataFrame:
    if homologs.empty:
        return _empty_clinvar_roc()
    if "score_urn" not in params or "sigma_neut" not in params:
        raise ValueError("params must contain score_urn and sigma_neut columns")

    parameter_lookup = params.loc[:, ["score_urn", "sigma_neut"]].rename(
        columns={"score_urn": "homolog_urn"}
    )
    annotated = homologs.loc[
        :, ["type", "homolog_urn", "wt", "homolog_wt", "same_wt", "mut", "w", "clnsig"]
    ].merge(parameter_lookup, how="left", on="homolog_urn", sort=False)
    annotated["cat"] = annotated["clnsig"].map(
        {
            "benign": False,
            "likely benign": False,
            "pathogenic": True,
            "likely pathogenic": True,
        }
    )
    annotated = annotated.loc[annotated["cat"].notna()].copy()
    if annotated.empty:
        return _empty_clinvar_roc()

    w = pd.to_numeric(annotated["w"], errors="coerce")
    sigma_neut = pd.to_numeric(annotated["sigma_neut"], errors="coerce")
    annotated["aw"] = (w - 1).abs()
    with np.errstate(divide="ignore", invalid="ignore"):
        annotated["p"] = w / sigma_neut
        annotated["ap"] = (w - 1).abs() / sigma_neut

    long = annotated.melt(
        id_vars=["type", "homolog_urn", "same_wt", "cat"],
        value_vars=["w", "aw", "p", "ap"],
        var_name="score_type",
        value_name="score",
    )
    expanded = _expanded_roc_rows(long)
    result_rows: list[pd.DataFrame] = []
    for group_key, group in expanded.groupby(
        list(_ROC_GROUP_COLUMNS), sort=False, dropna=False, observed=True
    ):
        score_type, type_value, same_wt_value, homolog_urn = group_key
        curve = calc_roc(
            group["cat"],
            group["score"],
            greater=score_type in {"aw", "ap"},
            max_steps=6000,
        )
        if curve.empty:
            continue
        curve["pr_auc"] = pr_auc(curve["tpr"], curve["precision"])
        for column, value in zip(_ROC_GROUP_COLUMNS, group_key):
            curve[column] = value
        result_rows.append(curve)

    if not result_rows:
        return _empty_clinvar_roc()
    roc = pd.concat(result_rows, ignore_index=True)
    roc = roc.merge(counts, how="left", on=["type", "same_wt"], sort=False)
    roc["n_variants"] = roc["n_variants"].astype("Int64")
    roc["n_studies"] = roc["n_studies"].astype("Int64")
    return roc.loc[
        :,
        [
            *_ROC_GROUP_COLUMNS,
            *_ROC_COLUMNS[:-1],
            "auc",
            "pr_auc",
            "n_variants",
            "n_studies",
        ],
    ]


def evaluate_clinvar(
    clinvar_variants: pd.DataFrame,
    scores: pd.DataFrame,
    params: pd.DataFrame,
    homology_key: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Return the ClinVar homolog tables, counts, curves, and summary.

    ``clinvar_counts`` follows the source R scope: ``n_variants`` counts
    homolog-score rows before non-binary labels are removed, and ``n_studies``
    counts distinct homolog score URNs. Counts are expanded over each homology
    type and wild-type agreement, plus the ``All``/``Both`` rollups.
    """
    combined = _label_clinvar_variants(clinvar_variants)
    homologs = _build_clinvar_homologs(combined, scores, homology_key)
    counts = _clinvar_counts(homologs)
    roc = _clinvar_roc(homologs, params, counts)
    if roc.empty:
        summary = pd.DataFrame(columns=_SUMMARY_COLUMNS)
    else:
        summary = (
            roc.loc[
                roc["score_type"].eq("p") & roc["homolog_urn"].eq("All"),
                list(_SUMMARY_COLUMNS),
            ]
            .drop_duplicates()
            .reset_index(drop=True)
        )
    return {
        "combined_clinvar": combined,
        "clinvar_homologs": homologs,
        "clinvar_counts": counts,
        "clinvar_roc": roc,
        "clinvar_summary": summary,
    }
