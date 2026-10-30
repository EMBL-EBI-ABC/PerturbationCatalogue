"""Same-gene and homologous-gene MAVE score correlations.

This ports the correlation portion of ``analysis/mave_homology.R``. The
calculation is intentionally directed: when both datasets are available, it
retains A-to-B and B-to-A rows, matching the original per-target workflow.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from data_sources.mavedb.analysis.gmm import classify_point_sd

HOMOLOGY_TYPES = {
    "ortholog_one2one": "Ortholog (one to one)",
    "ortholog_one2many": "Ortholog (one to many)",
    "ortholog_many2many": "Ortholog (many to many)",
    "other_paralog": "Paralog (diff. species)",
    "within_species_paralog": "Paralog (same species)",
    "same_gene": "Identical",
    "alt_allele": "Alt Allele",
}

MODE_NAMES = ("ko", "del", "neut", "pos", "noise")
HOMOLOGY_KEY_COLUMNS = (
    "taxonomy_level",
    "target_species",
    "source_species",
    "type",
    "target_id",
    "target_uniprot",
    "source_id",
    "source_uniprot",
    "target_wt",
    "source_wt",
    "target_perc_pos",
    "target_perc_id",
    "target_position",
    "source_position",
)
CORRELATION_COLUMNS = (
    "relationship",
    "study_relationship",
    "gene_relationship",
    "percent_identity",
    "percent_positions",
    "experiment_set_urn1",
    "experiment_urn1",
    "score_urn1",
    "species1",
    "uniprot1",
    "name1",
    "experiment_set_urn2",
    "experiment_urn2",
    "score_urn2",
    "species2",
    "uniprot2",
    "name2",
    "shared_variants",
    "pearson",
    "pearson_low",
    "pearson_high",
    "pearson_p",
    "spearman",
    "spearman_p",
    "cohens_kappa",
    "pearson_q",
    "spearman_q",
)

_META_COLUMNS = (
    "experiment_set_urn",
    "experiment_urn",
    "organism",
    "target_name",
    "identifier_uniprot",
)
_HOMOLOGY_REQUIRED = (
    "taxonomy_level",
    "target_species",
    "source_species",
    "type",
    "target_id",
    "source_id",
    "target_align_seq",
    "source_align_seq",
    "target_perc_pos",
    "target_perc_id",
)


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def _present(value: Any) -> bool:
    """Return whether a scalar is a usable identifier/name/amino acid."""

    if pd.isna(value):
        return False
    return str(value).strip() != ""


def _versionless_gene(value: Any) -> str | None:
    if not _present(value):
        return None
    return re.sub(r"\.\d+$", "", str(value).strip())


def _empty_homology_key() -> pd.DataFrame:
    result = pd.DataFrame(columns=HOMOLOGY_KEY_COLUMNS)
    result["target_position"] = pd.Series(dtype="Int64")
    result["source_position"] = pd.Series(dtype="Int64")
    return result


def build_homology_key(homology: pd.DataFrame, mave_ids: pd.DataFrame) -> pd.DataFrame:
    """Expand aligned Ensembl homology rows into per-residue coordinates.

    ``mave_ids`` may use either the source-table names ``From``/``To`` or the
    normalized names ``uniprot_id``/``ensembl_id``. Ensembl version suffixes
    are removed before matching, as in the source R analysis.
    """

    _require_columns(homology, _HOMOLOGY_REQUIRED, "homology")
    if {"From", "To"}.issubset(mave_ids.columns):
        id_map = mave_ids.loc[:, ["From", "To"]].rename(
            columns={"From": "uniprot_id", "To": "ensembl_id"}
        )
    elif {"uniprot_id", "ensembl_id"}.issubset(mave_ids.columns):
        id_map = mave_ids.loc[:, ["uniprot_id", "ensembl_id"]].copy()
    else:
        raise ValueError(
            "mave_ids must contain From/To or uniprot_id/ensembl_id columns"
        )

    id_map["uniprot_id"] = id_map["uniprot_id"].map(
        lambda value: str(value).strip() if _present(value) else None
    )
    id_map["ensembl_gene"] = id_map["ensembl_id"].map(_versionless_gene)
    id_map = id_map.dropna(subset=["uniprot_id", "ensembl_gene"])
    if id_map.empty:
        return _empty_homology_key()

    source = homology.loc[:, _HOMOLOGY_REQUIRED].copy()
    source["target_id"] = source["target_id"].map(_versionless_gene)
    source["source_id"] = source["source_id"].map(_versionless_gene)
    source = source.dropna(subset=["target_id", "source_id"])

    target_map = id_map.rename(
        columns={"ensembl_gene": "target_id", "uniprot_id": "target_uniprot"}
    )
    source_map = id_map.rename(
        columns={"ensembl_gene": "source_id", "uniprot_id": "source_uniprot"}
    )
    pairs = source.merge(target_map, how="inner", on="target_id", sort=False)
    pairs = pairs.merge(source_map, how="inner", on="source_id", sort=False)
    pairs = pairs.loc[pairs["target_uniprot"] != pairs["source_uniprot"]]

    records: list[dict[str, Any]] = []
    for row in pairs.itertuples(index=False):
        target_sequence = row.target_align_seq
        source_sequence = row.source_align_seq
        if not isinstance(target_sequence, str) or not isinstance(source_sequence, str):
            continue
        if len(target_sequence) != len(source_sequence):
            raise ValueError(
                "homology alignment sequences must have equal lengths for "
                f"{row.target_id!r} and {row.source_id!r}"
            )
        if not target_sequence:
            continue

        target_position = 0
        source_position = 0
        for target_wt, source_wt in zip(target_sequence, source_sequence):
            if target_wt == "-":
                mapped_target_position: Any = pd.NA
            else:
                target_position += 1
                mapped_target_position = target_position
            if source_wt == "-":
                mapped_source_position: Any = pd.NA
            else:
                source_position += 1
                mapped_source_position = source_position

            records.append(
                {
                    "taxonomy_level": row.taxonomy_level,
                    "target_species": row.target_species,
                    "source_species": row.source_species,
                    "type": row.type,
                    "target_id": row.target_id,
                    "target_uniprot": row.target_uniprot,
                    "source_id": row.source_id,
                    "source_uniprot": row.source_uniprot,
                    "target_wt": target_wt,
                    "source_wt": source_wt,
                    "target_perc_pos": row.target_perc_pos,
                    "target_perc_id": row.target_perc_id,
                    "target_position": mapped_target_position,
                    "source_position": mapped_source_position,
                }
            )

    if not records:
        return _empty_homology_key()
    result = pd.DataFrame.from_records(records, columns=HOMOLOGY_KEY_COLUMNS)
    result["target_position"] = pd.array(result["target_position"], dtype="Int64")
    result["source_position"] = pd.array(result["source_position"], dtype="Int64")
    return result


def correlation_stats(
    w1: Sequence[float],
    w2: Sequence[float],
    mode1: Sequence[str],
    mode2: Sequence[str],
) -> dict[str, int | float]:
    """Return R-compatible correlation, Fisher interval, and fixed-level kappa."""

    first = np.asarray(w1, dtype=float).reshape(-1)
    second = np.asarray(w2, dtype=float).reshape(-1)
    first_modes = pd.Series(mode1, dtype="object").reset_index(drop=True)
    second_modes = pd.Series(mode2, dtype="object").reset_index(drop=True)
    if not (len(first) == len(second) == len(first_modes) == len(second_modes)):
        raise ValueError("scores and modes must have equal lengths")
    if len(first) == 0:
        return {}
    if len(first) < 3:
        raise ValueError("at least three shared variants are required")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("scores must be finite")
    if np.unique(first).size < 2 or np.unique(second).size < 2:
        raise ValueError("scores must be nonconstant")

    pearson_test = stats.pearsonr(first, second)
    pearson_interval = (
        pearson_test.confidence_interval(confidence_level=0.95)
        if len(first) > 3
        else None
    )
    spearman_test = stats.spearmanr(first, second)

    if first_modes.isna().any() or second_modes.isna().any():
        kappa = float("nan")
    else:
        first_codes = pd.Categorical(first_modes, categories=MODE_NAMES).codes
        second_codes = pd.Categorical(second_modes, categories=MODE_NAMES).codes
        if (first_codes < 0).any() or (second_codes < 0).any():
            kappa = float("nan")
        else:
            first_counts = np.bincount(first_codes, minlength=len(MODE_NAMES))
            second_counts = np.bincount(second_codes, minlength=len(MODE_NAMES))
            observed = float(np.mean(first_codes == second_codes))
            expected = float(np.dot(first_counts, second_counts) / len(first) ** 2)
            kappa = (
                (observed - expected) / (1.0 - expected)
                if expected != 1.0
                else float("nan")
            )

    return {
        "shared_variants": int(len(first)),
        "pearson": float(pearson_test.statistic),
        "pearson_low": (
            float(pearson_interval.low) if pearson_interval else float("nan")
        ),
        "pearson_high": (
            float(pearson_interval.high) if pearson_interval else float("nan")
        ),
        "pearson_p": float(pearson_test.pvalue),
        "spearman": float(spearman_test.statistic),
        "spearman_p": float(spearman_test.pvalue),
        "cohens_kappa": float(kappa),
    }


def _p_adjust_fdr(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjustment matching R's ``p.adjust(..., 'fdr')``."""

    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if finite_indices.size == 0:
        return adjusted
    order = np.argsort(values[finite_indices], kind="mergesort")
    sorted_indices = finite_indices[order]
    sorted_values = values[sorted_indices]
    count = len(sorted_values)
    scaled = sorted_values * count / np.arange(1, count + 1)
    sorted_adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    adjusted[sorted_indices] = np.minimum(sorted_adjusted, 1.0)
    return adjusted


def _is_false(value: Any) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() == "false"
    return value is False or isinstance(value, (bool, np.bool_)) and not bool(value)


def _normalize_same_gene_species(value: Any) -> Any:
    if not _present(value):
        return pd.NA
    return str(value).lower().replace(" ", "_", 1)


def _study_relationship(first: pd.Series, second: pd.Series) -> str:
    first_experiment = first["experiment_urn"]
    second_experiment = second["experiment_urn"]
    if (
        _present(first_experiment)
        and _present(second_experiment)
        and first_experiment == second_experiment
    ):
        return "Same Experiment"
    first_set = first["experiment_set_urn"]
    second_set = second["experiment_set_urn"]
    if _present(first_set) and _present(second_set) and first_set == second_set:
        return "Same Study"
    return "Different Study"


def _classified_modes(values: np.ndarray, sigma_neut: Any) -> np.ndarray:
    sigma = float(sigma_neut) if _present(sigma_neut) else float("nan")
    return classify_point_sd(values, sigma, 0.99)


def _correlation_row(
    first_scores: np.ndarray,
    second_scores: np.ndarray,
    first_sigma: Any,
    second_sigma: Any,
) -> dict[str, int | float]:
    return correlation_stats(
        first_scores,
        second_scores,
        _classified_modes(first_scores, first_sigma),
        _classified_modes(second_scores, second_sigma),
    )


def _base_output_row(
    *,
    score_urn1: Any,
    score_urn2: Any,
    first: pd.Series,
    second: pd.Series,
    gene_relationship: str,
    percent_identity: Any,
    percent_positions: Any,
    species1: Any,
    species2: Any,
) -> dict[str, Any]:
    study_relationship = _study_relationship(first, second)
    return {
        "relationship": f"{study_relationship} - {gene_relationship}",
        "study_relationship": study_relationship,
        "gene_relationship": gene_relationship,
        "percent_identity": percent_identity,
        "percent_positions": percent_positions,
        "experiment_set_urn1": first["experiment_set_urn"],
        "experiment_urn1": first["experiment_urn"],
        "score_urn1": score_urn1,
        "species1": species1,
        "uniprot1": first["identifier_uniprot"],
        "name1": first["target_name"],
        "experiment_set_urn2": second["experiment_set_urn"],
        "experiment_urn2": second["experiment_urn"],
        "score_urn2": score_urn2,
        "species2": species2,
        "uniprot2": second["identifier_uniprot"],
        "name2": second["target_name"],
    }


def _prepared_inputs(
    meta: pd.DataFrame,
    scores: pd.DataFrame,
    singles: pd.DataFrame,
    params: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[Any, pd.Series], dict[Any, float]]:
    _require_columns(
        meta,
        ["score_urn", "meta_study", *_META_COLUMNS],
        "meta",
    )
    _require_columns(scores, ["score_urn", "hgvs_pro", "w"], "scores")
    _require_columns(singles, ["score_urn", "position", "wt", "mut", "w"], "singles")
    _require_columns(params, ["score_urn", "sigma_neut"], "params")

    if meta["score_urn"].duplicated().any():
        raise ValueError("meta must contain at most one row per score_urn")
    if params["score_urn"].duplicated().any():
        raise ValueError("params must contain at most one selected row per score_urn")

    eligible_meta = meta.loc[meta["meta_study"].map(_is_false)].copy()
    eligible_meta = eligible_meta.loc[eligible_meta["score_urn"].map(_present)]
    metadata = eligible_meta.set_index("score_urn", drop=False)
    meta_lookup = {score_urn: row for score_urn, row in metadata.iterrows()}

    score_values = scores.loc[:, ["score_urn", "hgvs_pro", "w"]].copy()
    score_values["w"] = pd.to_numeric(score_values["w"], errors="coerce")
    score_w = score_values["w"].to_numpy(dtype=float, na_value=np.nan)
    score_values = score_values.loc[
        score_values["score_urn"].isin(meta_lookup) & np.isfinite(score_w)
    ]
    score_values = score_values.merge(
        eligible_meta.loc[:, ["score_urn", *_META_COLUMNS]],
        how="inner",
        on="score_urn",
        sort=False,
        validate="many_to_one",
    )

    single_values = singles.loc[:, ["score_urn", "position", "wt", "mut", "w"]].copy()
    single_values["w"] = pd.to_numeric(single_values["w"], errors="coerce")
    single_values["position"] = pd.to_numeric(
        single_values["position"], errors="coerce"
    )
    single_w = single_values["w"].to_numpy(dtype=float, na_value=np.nan)
    single_position = single_values["position"].to_numpy(dtype=float, na_value=np.nan)
    single_values = single_values.loc[
        single_values["score_urn"].isin(meta_lookup)
        & np.isfinite(single_w)
        & np.isfinite(single_position)
    ]
    single_values = single_values.merge(
        eligible_meta.loc[:, ["score_urn", *_META_COLUMNS]],
        how="inner",
        on="score_urn",
        sort=False,
        validate="many_to_one",
    )

    parameter_values = params.loc[:, ["score_urn", "sigma_neut"]].copy()
    sigma_lookup = {
        row.score_urn: row.sigma_neut
        for row in parameter_values.itertuples(index=False)
    }
    return score_values, single_values, meta_lookup, sigma_lookup


def calculate_correlations(
    meta: pd.DataFrame,
    scores: pd.DataFrame,
    singles: pd.DataFrame,
    params: pd.DataFrame,
    homology_key: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate directed same-gene and homologous-gene score correlations.

    Metadata is rejoined by ``score_urn`` inside this function so callers may
    pass score and singles tables with only their native measurement fields.
    Meta-studies, nonfinite scores, self-pairs, and pairs with fewer than three
    shared nonconstant values are omitted.
    """

    _require_columns(
        homology_key,
        HOMOLOGY_KEY_COLUMNS,
        "homology_key",
    )
    score_values, single_values, metadata, sigma_lookup = _prepared_inputs(
        meta, scores, singles, params
    )
    if score_values.empty:
        return pd.DataFrame(columns=CORRELATION_COLUMNS)

    urns = list(pd.unique(score_values["score_urn"]))
    hgvs_means = {
        urn: group.groupby("hgvs_pro", dropna=False, sort=False, as_index=False)[
            "w"
        ].mean()
        for urn, group in score_values.groupby("score_urn", sort=False)
    }
    singles_by_urn = {
        urn: group for urn, group in single_values.groupby("score_urn", sort=False)
    }
    singles_by_uniprot = {
        uniprot: group
        for uniprot, group in single_values.loc[
            single_values["identifier_uniprot"].map(_present)
        ].groupby("identifier_uniprot", sort=False)
    }
    homology_by_target = {
        uniprot: group
        for uniprot, group in homology_key.loc[
            homology_key["target_uniprot"].map(_present)
        ].groupby("target_uniprot", sort=False)
    }
    uniprot_groups: dict[str, list[Any]] = {}
    name_groups: dict[str, list[Any]] = {}
    for urn in urns:
        row = metadata[urn]
        identifier = row["identifier_uniprot"]
        if _present(identifier):
            uniprot_groups.setdefault(str(identifier).strip(), []).append(urn)
        if _present(row["target_name"]):
            name_groups.setdefault(str(row["target_name"]).strip(), []).append(urn)

    output: list[dict[str, Any]] = []
    for urn in urns:
        first = metadata[urn]
        first_scores = hgvs_means[urn].rename(columns={"w": "w1"})

        identifier = first["identifier_uniprot"]
        if _present(identifier):
            candidate_urns = uniprot_groups.get(str(identifier).strip(), [])
        elif _present(first["target_name"]):
            # R falls back by name for the current target, including comparator
            # datasets that have UniProt identifiers of their own.
            candidate_urns = name_groups.get(str(first["target_name"]).strip(), [])
        else:
            candidate_urns = []

        if candidate_urns:
            for second_urn in candidate_urns:
                if second_urn == urn:
                    continue
                shared = first_scores.merge(
                    hgvs_means[second_urn].rename(columns={"w": "w2"}),
                    how="inner",
                    on="hgvs_pro",
                    sort=False,
                )
                if len(shared) < 3:
                    continue
                w1 = shared["w1"].to_numpy(dtype=float)
                w2 = shared["w2"].to_numpy(dtype=float)
                if np.unique(w1).size < 2 or np.unique(w2).size < 2:
                    continue
                metrics = _correlation_row(
                    w1,
                    w2,
                    sigma_lookup.get(urn, np.nan),
                    sigma_lookup.get(second_urn, np.nan),
                )
                second = metadata[second_urn]
                row_out = _base_output_row(
                    score_urn1=urn,
                    score_urn2=second_urn,
                    first=first,
                    second=second,
                    gene_relationship="Identical",
                    percent_identity=100,
                    percent_positions=100,
                    species1=_normalize_same_gene_species(first["organism"]),
                    species2=_normalize_same_gene_species(second["organism"]),
                )
                row_out.update(metrics)
                output.append(row_out)

        if not _present(identifier) or single_values.empty or homology_key.empty:
            continue
        target_key = homology_by_target.get(
            str(identifier).strip(), homology_key.iloc[:0]
        )
        target_singles = singles_by_urn.get(urn, single_values.iloc[:0])
        target_key = target_key.loc[
            target_key["target_position"].notna()
            & target_key["source_position"].notna()
            & target_key["target_wt"].map(_present)
            & target_key["source_wt"].map(_present)
        ]
        target_singles = target_singles.loc[
            target_singles["position"].notna()
            & target_singles["wt"].map(_present)
            & target_singles["mut"].map(_present)
        ]
        if target_key.empty or target_singles.empty:
            continue

        target_matches = target_key.merge(
            target_singles.loc[:, ["position", "wt", "mut", "score_urn", "w"]],
            how="inner",
            left_on=["target_position", "target_wt"],
            right_on=["position", "wt"],
            sort=False,
        ).rename(columns={"score_urn": "target_urn", "w": "target_w"})
        if target_matches.empty:
            continue

        source_ids = list(pd.unique(target_key["source_uniprot"]))
        source_parts = [
            singles_by_uniprot[uniprot]
            for uniprot in source_ids
            if uniprot in singles_by_uniprot
        ]
        if not source_parts:
            continue
        source_singles = pd.concat(source_parts, ignore_index=True)
        source_singles = source_singles.loc[
            source_singles["position"].notna()
            & source_singles["wt"].map(_present)
            & source_singles["mut"].map(_present)
        ].loc[:, ["identifier_uniprot", "score_urn", "position", "wt", "mut", "w"]]
        homolog_matches = target_matches.merge(
            source_singles,
            how="inner",
            left_on=["source_uniprot", "source_position", "source_wt", "mut"],
            right_on=["identifier_uniprot", "position", "wt", "mut"],
            sort=False,
            suffixes=("_target", "_source"),
        ).rename(columns={"score_urn": "source_urn", "w": "source_w"})
        homolog_matches = homolog_matches.loc[
            homolog_matches["target_urn"] != homolog_matches["source_urn"]
        ]
        if homolog_matches.empty:
            continue

        group_columns = [
            "type",
            "target_urn",
            "source_urn",
            "target_species",
            "source_species",
            "target_uniprot",
            "source_uniprot",
            "target_perc_pos",
            "target_perc_id",
        ]
        for keys, group in homolog_matches.groupby(
            group_columns, dropna=False, sort=False
        ):
            target_w = group["target_w"].to_numpy(dtype=float)
            source_w = group["source_w"].to_numpy(dtype=float)
            if (
                len(group) < 3
                or np.unique(target_w).size < 2
                or np.unique(source_w).size < 2
            ):
                continue
            (
                raw_type,
                target_urn,
                source_urn,
                target_species,
                source_species,
                target_uniprot,
                source_uniprot,
                target_perc_pos,
                target_perc_id,
            ) = keys
            gene_relationship = HOMOLOGY_TYPES.get(raw_type)
            if gene_relationship is None:
                continue
            metrics = _correlation_row(
                source_w,
                target_w,
                sigma_lookup.get(source_urn, np.nan),
                sigma_lookup.get(target_urn, np.nan),
            )
            first_homolog = metadata[source_urn]
            second_homolog = metadata[target_urn]
            row_out = _base_output_row(
                score_urn1=source_urn,
                score_urn2=target_urn,
                first=first_homolog,
                second=second_homolog,
                gene_relationship=gene_relationship,
                percent_identity=target_perc_id,
                percent_positions=target_perc_pos,
                species1=source_species,
                species2=target_species,
            )
            # Homology metadata is authoritative for the gene endpoints and
            # species, while names come from the enriched MAVE metadata.
            row_out["uniprot1"] = source_uniprot
            row_out["uniprot2"] = target_uniprot
            row_out.update(metrics)
            output.append(row_out)

    if not output:
        return pd.DataFrame(columns=CORRELATION_COLUMNS)

    result = pd.DataFrame.from_records(output)
    result = result.loc[result["score_urn1"] != result["score_urn2"]].copy()
    result["pearson_q"] = _p_adjust_fdr(result["pearson_p"].to_numpy(dtype=float))
    result["spearman_q"] = _p_adjust_fdr(result["spearman_p"].to_numpy(dtype=float))
    return result.loc[:, CORRELATION_COLUMNS].reset_index(drop=True)
