"""Prepare the tabular statistics used by the MAVE QC figures.

The source tables are intentionally read in chunks. Eligible score rows are
staged in a temporary SQLite database so the full score and singles tables do
not need to be held in memory at once.
"""

from __future__ import annotations

import json
import logging
import math
import re
import sqlite3
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices
from scipy import stats

LOGGER = logging.getLogger(__name__)

_METADATA_REQUIRED = {
    "dataset_id",
    "experiment_id",
    "experiment_set_id",
    "primary_analysis",
    "n_analyzed_variants",
    "target_name",
    "organism",
    "taxid",
    "meta_study",
    "single_mutations",
    "multiple_mutations",
}
_METADATA_NUMERIC = (
    "n_analyzed_variants",
    "completeness",
    "completeness_missense",
    "stop_proportion",
    "del_proportion",
    "synonymous_proportion",
    "coverage_50",
    "coverage_75",
    "coverage_100",
    "coverage_nonsense",
    "coverage_synonymous",
    "mean_muts",
    "mean_backgrounds",
    "n_positions",
    "n_subs",
)
_COMPLETENESS_METRICS = (
    "completeness",
    "completeness_missense",
    "stop_proportion",
    "del_proportion",
    "synonymous_proportion",
    "wt_mutants",
)
_COVERAGE_METRICS = (
    "coverage_50",
    "coverage_75",
    "coverage_100",
    "coverage_nonsense",
    "coverage_synonymous",
)
_DATASET_COLUMNS = ["name", "mutation_type", "n"]
_ORGANISM_COLUMNS = ["taxid", "organism", "kingdom", "n"]
_SET_SUMMARY_COLUMNS = ["type", "urn", "n", "mut_types"]
_GENE_COLUMNS = [
    "target_name",
    "score_sets",
    "experiments",
    "experiment_sets",
    "meta_analyses",
    "measurements",
]
_VARIANT_COLUMNS = [
    "target_name",
    "hgvs_pro",
    "wt_mut",
    "score_sets",
    "experiments",
    "experiment_sets",
    "measurements",
]
_METRIC_COLUMNS = ["score_urn", "human", "n", "metric", "value"]
_MULTIPLE_MUTATION_COLUMNS = ["score_urn", "mean_muts", "mean_backgrounds"]
_CORRELATION_COLUMNS = [
    "target_name",
    "organism",
    "set_1",
    "set_2",
    "spearman",
    "pearson",
    "match_type",
]
_EVOLUTION_COLUMNS = [
    "score_urn",
    "organism",
    "kingdom",
    "estimate",
    "statistic",
    "pvalue",
    "n_positions",
    "n_subs",
]


def prepare_qc_data(
    scores_path: Path,
    metadata_path: Path,
    singles_path: Path,
    *,
    chunk_size: int = 250_000,
    matrix_path: Path | None = None,
) -> dict[str, pd.DataFrame]:
    """Return the tables underlying the MAVE QC plots.

    Only rows assigned to ``primary_analysis == "maves"`` are considered.
    Score-based variant summaries and score-set correlations use non-meta
    studies, finite numeric scores, and non-empty protein HGVS values. Rows
    whose target metadata is a multi-target JSON array remain in metadata-wide
    summaries but are omitted from target-attributed calculations.
    """
    scores_path = Path(scores_path)
    metadata_path = Path(metadata_path)
    singles_path = Path(singles_path)
    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    for label, path in (
        ("scores", scores_path),
        ("metadata", metadata_path),
        ("singles", singles_path),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} input file does not exist: {path}")

    metadata = _load_metadata(metadata_path)
    matrix_file = matrix_path or Path(__file__).parent / "data" / "BLOSUM100.txt"
    matrix = _load_blosum_matrix(Path(matrix_file))

    target_ids, group_ids, score_ids, score_metadata = _make_score_dimensions(metadata)
    non_meta_score_urns = set(
        metadata.loc[~metadata["meta_study"], "score_urn"].dropna().tolist()
    )
    target_score_urns = {
        score_urn
        for score_urn, target_id in target_ids.items()
        if target_id is not None
    }
    evo_score_urns = target_score_urns
    target_name_by_id = {
        target_ids[row.score_urn]: row.target_name
        for row in metadata.itertuples()
        if row.score_urn is not None and target_ids.get(row.score_urn) is not None
    }
    group_label_by_id = {
        group_ids[row.score_urn]: (row.target_name, row.organism)
        for row in metadata.itertuples()
        if row.score_urn is not None and group_ids.get(row.score_urn) is not None
    }

    outputs = _metadata_outputs(metadata)

    with tempfile.TemporaryDirectory(prefix="mave-qc-") as temp_dir:
        database_path = Path(temp_dir) / "qc.sqlite3"
        connection = sqlite3.connect(str(database_path))
        try:
            _configure_database(connection)
            _create_database_tables(connection)
            with connection:
                connection.executemany(
                    """
                    INSERT INTO score_metadata
                        (score_id, score_urn, target_id, group_id,
                         experiment_urn, experiment_set_urn)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    score_metadata,
                )

            ambiguous_urns = set(
                metadata.loc[metadata["target_is_multi"], "score_urn"].dropna().tolist()
            )
            eligible_score_rows, staged_score_rows, ambiguous_score_rows = (
                _stage_scores(
                    connection,
                    scores_path,
                    score_ids,
                    target_ids,
                    group_ids,
                    non_meta_score_urns,
                    ambiguous_urns,
                    chunk_size,
                )
            )
            _create_score_indexes(connection)
            LOGGER.info(
                "Staged %d of %d eligible non-meta score rows across target groups",
                staged_score_rows,
                eligible_score_rows,
            )
            if ambiguous_urns:
                LOGGER.warning(
                    "Omitted %d eligible score rows from %d multi-target datasets "
                    "from target-attributed QC tables; metadata summaries retain each dataset once",
                    ambiguous_score_rows,
                    len(ambiguous_urns),
                )

            outputs["data_per_variant"] = _data_per_variant(
                connection, target_name_by_id
            )
            LOGGER.info(
                "Prepared %d target-attributed variant rows from %d staged score rows",
                len(outputs["data_per_variant"]),
                staged_score_rows,
            )
            outputs["score_set_correlation"] = _score_set_correlations(
                connection,
                group_label_by_id,
                score_ids,
            )
            LOGGER.info(
                "Prepared %d score-set correlation pairs",
                len(outputs["score_set_correlation"]),
            )

            _stage_singles(
                connection,
                singles_path,
                score_ids,
                evo_score_urns,
                matrix,
                chunk_size,
            )
            _create_evolution_index(connection)
            outputs["study_evo_correlation"] = _evolution_correlations(
                connection,
                metadata,
                score_ids,
                evo_score_urns,
            )
        finally:
            connection.close()

    return outputs


def _empty_outputs() -> dict[str, pd.DataFrame]:
    """Create all result keys with stable columns for empty inputs."""
    return {
        "datasets": pd.DataFrame(columns=_DATASET_COLUMNS),
        "organisms": pd.DataFrame(columns=_ORGANISM_COLUMNS),
        "set_summary": pd.DataFrame(columns=_SET_SUMMARY_COLUMNS),
        "data_per_gene": pd.DataFrame(columns=_GENE_COLUMNS),
        "data_per_variant": pd.DataFrame(columns=_VARIANT_COLUMNS),
        "completeness": pd.DataFrame(columns=_METRIC_COLUMNS),
        "coverage": pd.DataFrame(columns=_METRIC_COLUMNS),
        "multiple_mutation": pd.DataFrame(columns=_MULTIPLE_MUTATION_COLUMNS),
        "score_set_correlation": pd.DataFrame(columns=_CORRELATION_COLUMNS),
        "study_evo_correlation": pd.DataFrame(columns=_EVOLUTION_COLUMNS),
    }


def _load_metadata(path: Path) -> pd.DataFrame:
    """Read and normalize metadata while preserving literal target ``NA``."""
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing_columns = sorted(_METADATA_REQUIRED - set(frame.columns))
    if missing_columns:
        raise ValueError(
            "metadata input is missing required columns: " + ", ".join(missing_columns)
        )

    dataset_ids = frame["dataset_id"].map(_identifier)
    duplicated = dataset_ids.notna() & dataset_ids.duplicated(keep=False)
    if duplicated.any():
        examples = sorted(dataset_ids[duplicated].dropna().unique().tolist())[:5]
        raise ValueError(
            "metadata contains duplicate dataset_id values; expected one row per "
            f"score set (examples: {', '.join(examples)})"
        )

    for column in _METADATA_NUMERIC:
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if "wt_mutants" not in frame.columns:
        frame["wt_mutants"] = ""
    frame["wt_mutants"] = frame["wt_mutants"].map(_numeric_bool)
    frame["n"] = frame["n_analyzed_variants"]

    for column in ("meta_study", "single_mutations", "multiple_mutations"):
        frame[column] = frame[column].map(_parse_bool)

    frame["score_urn"] = dataset_ids
    frame["experiment_urn"] = frame["experiment_id"].map(_identifier)
    frame["experiment_set_urn"] = frame["experiment_set_id"].map(_identifier)
    frame["primary_analysis"] = frame["primary_analysis"].map(_plain_text)
    frame["target_name_raw"] = frame["target_name"].map(
        lambda value: _global_text(value, preserve_literal_na=True)
    )
    frame["target_name"] = frame["target_name"].map(
        lambda value: _scalar_json_text(value, preserve_literal_na=True)
    )
    frame["target_is_multi"] = frame["target_name_raw"].map(_is_multi_json)
    frame["organism_raw"] = frame["organism"].map(_global_text)
    frame["organism"] = frame["organism"].map(_scalar_json_text)
    frame["taxid_raw"] = frame["taxid"].map(_global_text)
    frame["taxid"] = frame["taxid"].map(_scalar_json_text)
    if "kingdom" not in frame.columns:
        frame["kingdom"] = ""
    frame["kingdom"] = frame["kingdom"].map(_plain_text)
    frame["primary_analysis"] = frame["primary_analysis"].str.lower()
    return frame.loc[frame["primary_analysis"] == "maves"].copy()


def _make_score_dimensions(
    metadata: pd.DataFrame,
) -> tuple[
    dict[str, int | None], dict[str, int | None], dict[str, int], list[tuple[Any, ...]]
]:
    """Assign compact IDs to metadata dimensions used by SQLite staging."""
    score_urns = sorted(metadata["score_urn"].dropna().unique().tolist())
    score_ids = {score_urn: index for index, score_urn in enumerate(score_urns)}

    target_names = sorted(metadata["target_name"].dropna().unique().tolist())
    target_id_by_name = {
        target_name: index for index, target_name in enumerate(target_names)
    }

    group_keys = {
        (row.target_name, row.organism)
        for row in metadata.itertuples()
        if row.score_urn is not None and row.target_name is not None
    }
    sorted_group_keys = sorted(
        group_keys,
        key=lambda key: (key[0], key[1] is not None, key[1] or ""),
    )
    group_id_by_key = {
        group_key: index for index, group_key in enumerate(sorted_group_keys)
    }

    target_id_by_score: dict[str, int | None] = {}
    group_id_by_score: dict[str, int | None] = {}
    metadata_by_score = metadata.set_index("score_urn", drop=False)
    dimension_rows: list[tuple[Any, ...]] = []
    for score_urn in score_urns:
        row = metadata_by_score.loc[score_urn]
        target_name = row["target_name"]
        organism = row["organism"]
        target_id = target_id_by_name.get(target_name)
        group_id = (
            group_id_by_key.get((target_name, organism))
            if target_id is not None
            else None
        )
        target_id_by_score[score_urn] = target_id
        group_id_by_score[score_urn] = group_id
        dimension_rows.append(
            (
                score_ids[score_urn],
                score_urn,
                target_id,
                group_id,
                _sqlite_text(row["experiment_urn"]),
                _sqlite_text(row["experiment_set_urn"]),
            )
        )

    return target_id_by_score, group_id_by_score, score_ids, dimension_rows


def _metadata_outputs(metadata: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build metadata-only summaries and metric tables."""
    outputs = _empty_outputs()
    if metadata.empty:
        return outputs

    mutation_types = metadata.apply(_mutation_type, axis=1)
    metadata = metadata.assign(mutation_type=mutation_types)
    outputs["datasets"] = _dataset_counts(metadata)
    outputs["organisms"] = _organism_counts(metadata)

    non_meta = metadata.loc[~metadata["meta_study"]].copy()
    outputs["set_summary"] = _set_summary(non_meta)
    outputs["data_per_gene"] = _data_per_gene(metadata)
    outputs["completeness"] = _long_metrics(non_meta, _COMPLETENESS_METRICS)
    outputs["coverage"] = _long_metrics(non_meta, _COVERAGE_METRICS)
    outputs["multiple_mutation"] = _multiple_mutation(metadata)
    return outputs


def _dataset_counts(metadata: pd.DataFrame) -> pd.DataFrame:
    columns = [
        ("score_urn", "score_urn"),
        ("experiment_set_urn", "experiment_set_urn"),
        ("experiment_urn", "experiment_urn"),
        ("target_name_raw", "target_name"),
        ("organism_raw", "organism"),
    ]
    records: list[dict[str, Any]] = []
    meta_urn = re.compile(r"urn:mavedb:\d+-\d+-\d+")
    for _, row in metadata.iterrows():
        for source_column, name in columns:
            value = row[source_column]
            if isinstance(value, str) and meta_urn.search(value):
                name = "meta_score_urn"
            records.append(
                {
                    "name": name,
                    "mutation_type": row["mutation_type"],
                    "value": value,
                }
            )
    frame = pd.DataFrame.from_records(records)
    if frame.empty:
        return pd.DataFrame(columns=_DATASET_COLUMNS)
    result = (
        frame.groupby(["name", "mutation_type"], dropna=False)["value"]
        .nunique(dropna=False)
        .rename("n")
        .reset_index()
    )
    return (
        result[_DATASET_COLUMNS]
        .sort_values(["name", "mutation_type"], kind="stable")
        .reset_index(drop=True)
    )


def _organism_counts(metadata: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        metadata.groupby(
            ["taxid_raw", "organism_raw", "kingdom"], dropna=False, sort=False
        )
        .size()
        .rename("n")
        .reset_index()
        .rename(columns={"taxid_raw": "taxid", "organism_raw": "organism"})
    )
    grouped["kingdom"] = grouped["kingdom"].fillna("Other")
    return (
        grouped[_ORGANISM_COLUMNS]
        .sort_values("n", kind="stable")
        .reset_index(drop=True)
    )


def _set_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    summaries: list[pd.DataFrame] = []
    for summary_type, urn_column, count_column in (
        ("experiment_set", "experiment_set_urn", "experiment_urn"),
        ("experiment", "experiment_urn", "score_urn"),
    ):
        rows: list[dict[str, Any]] = []
        for urn, group in metadata.groupby(urn_column, dropna=False, sort=False):
            mutation_types = group["mutation_type"].drop_duplicates().tolist()
            rows.append(
                {
                    "type": summary_type,
                    "urn": urn,
                    "n": group[count_column].nunique(dropna=False),
                    "mut_types": (
                        mutation_types[0] if len(mutation_types) == 1 else "Mixed"
                    ),
                }
            )
        summaries.append(pd.DataFrame.from_records(rows, columns=_SET_SUMMARY_COLUMNS))

    score_rows = metadata[["score_urn", "n", "mutation_type"]].rename(
        columns={"score_urn": "urn", "n": "n", "mutation_type": "mut_types"}
    )
    score_rows.insert(0, "type", "score")
    summaries.append(score_rows[_SET_SUMMARY_COLUMNS])
    combined = pd.concat(summaries, ignore_index=True)
    return combined[_SET_SUMMARY_COLUMNS]


def _data_per_gene(metadata: pd.DataFrame) -> pd.DataFrame:
    target_metadata = metadata.loc[metadata["target_name"].notna()]
    records: list[dict[str, Any]] = []
    for target_name, group in target_metadata.groupby(
        "target_name", dropna=False, sort=False
    ):
        non_meta = group.loc[~group["meta_study"]]
        measurements = group["n"].sum(skipna=False)
        records.append(
            {
                "target_name": target_name,
                "score_sets": non_meta["score_urn"].nunique(dropna=False),
                "experiments": non_meta["experiment_urn"].nunique(dropna=False),
                "experiment_sets": non_meta["experiment_set_urn"].nunique(dropna=False),
                "meta_analyses": int(group["meta_study"].sum()),
                "measurements": measurements,
            }
        )
    result = pd.DataFrame.from_records(records, columns=_GENE_COLUMNS)
    return result.sort_values("target_name", kind="stable").reset_index(drop=True)


def _long_metrics(
    metadata: pd.DataFrame, metric_names: tuple[str, ...]
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in metadata.iterrows():
        human = "Human" if row["organism"] == "Homo sapiens" else "Other"
        for metric in metric_names:
            value = row[metric]
            if pd.isna(value):
                continue
            rows.append(
                {
                    "score_urn": row["score_urn"],
                    "human": human,
                    "n": row["n"],
                    "metric": metric,
                    "value": value,
                }
            )
    return pd.DataFrame.from_records(rows, columns=_METRIC_COLUMNS)


def _multiple_mutation(metadata: pd.DataFrame) -> pd.DataFrame:
    selected = metadata.loc[metadata["mean_muts"].notna()]
    result = selected[["score_urn", "mean_muts", "mean_backgrounds"]].copy()
    return result.reset_index(drop=True)[_MULTIPLE_MUTATION_COLUMNS]


def _configure_database(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = FILE")
    connection.execute("PRAGMA cache_size = -32768")


def _create_database_tables(connection: sqlite3.Connection) -> None:
    connection.executescript("""
        CREATE TABLE score_metadata (
            score_id INTEGER PRIMARY KEY,
            score_urn TEXT NOT NULL UNIQUE,
            target_id INTEGER,
            group_id INTEGER,
            experiment_urn TEXT,
            experiment_set_urn TEXT
        );
        CREATE TABLE variant_scores (
            score_id INTEGER NOT NULL,
            target_id INTEGER NOT NULL,
            group_id INTEGER NOT NULL,
            hgvs_pro TEXT NOT NULL,
            wt_mut INTEGER,
            score REAL NOT NULL
        );
        CREATE TABLE evo_values (
            score_id INTEGER NOT NULL,
            score REAL NOT NULL,
            blosum100 REAL NOT NULL
        );
        """)


def _stage_scores(
    connection: sqlite3.Connection,
    path: Path,
    score_ids: dict[str, int],
    target_ids: dict[str, int | None],
    group_ids: dict[str, int | None],
    non_meta_score_urns: set[str],
    ambiguous_score_urns: set[str],
    chunk_size: int,
) -> tuple[int, int, int]:
    required = {"dataset_id", "hgvs_pro", "score", "wt_mut"}
    _validate_table_columns(path, "scores", required)
    usecols = sorted(required)
    rows_seen = 0
    eligible_rows = 0
    inserted_rows = 0
    ambiguous_rows = 0

    LOGGER.info("Reading scores in chunks of %d rows: %s", chunk_size, path)
    for chunk_number, chunk in enumerate(
        pd.read_csv(
            path,
            sep="\t",
            usecols=usecols,
            dtype=str,
            keep_default_na=False,
            chunksize=chunk_size,
        ),
        start=1,
    ):
        chunk = chunk.reset_index(drop=True)
        rows_seen += len(chunk)
        score_values = pd.to_numeric(chunk["score"], errors="coerce").to_numpy(
            dtype=float
        )
        finite = np.isfinite(score_values)
        valid_hgvs = chunk["hgvs_pro"].map(_is_protein_hgvs).to_numpy(dtype=bool)
        is_non_meta = chunk["dataset_id"].isin(non_meta_score_urns).to_numpy(dtype=bool)
        eligible = finite & valid_hgvs & is_non_meta
        eligible_rows += int(eligible.sum())
        if not eligible.any():
            if chunk_number % 4 == 0:
                LOGGER.info("Read %d score rows", rows_seen)
            continue

        score_urn_values = chunk["dataset_id"].to_numpy(dtype=object)
        eligible_indices = np.flatnonzero(eligible)
        staged: list[tuple[Any, ...]] = []
        for row_index in eligible_indices:
            score_urn = score_urn_values[row_index]
            target_id = target_ids.get(score_urn)
            group_id = group_ids.get(score_urn)
            if target_id is None or group_id is None:
                if score_urn in ambiguous_score_urns:
                    ambiguous_rows += 1
                continue
            staged.append(
                (
                    score_ids[score_urn],
                    target_id,
                    group_id,
                    str(chunk.at[row_index, "hgvs_pro"]).strip(),
                    _bool_code(chunk.at[row_index, "wt_mut"]),
                    float(score_values[row_index]),
                )
            )

        if staged:
            with connection:
                connection.executemany(
                    """
                    INSERT INTO variant_scores
                        (score_id, target_id, group_id, hgvs_pro, wt_mut, score)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    staged,
                )
            inserted_rows += len(staged)
        if chunk_number % 4 == 0:
            LOGGER.info("Read %d score rows; staged %d", rows_seen, inserted_rows)

    LOGGER.info("Finished reading %d score rows", rows_seen)
    return eligible_rows, inserted_rows, ambiguous_rows


def _create_score_indexes(connection: sqlite3.Connection) -> None:
    with connection:
        connection.execute(
            "CREATE INDEX idx_variant_scores_target ON variant_scores(target_id, hgvs_pro, wt_mut)"
        )
        connection.execute(
            "CREATE INDEX idx_variant_scores_group ON variant_scores(group_id, score_id, hgvs_pro)"
        )


def _data_per_variant(
    connection: sqlite3.Connection,
    target_name_by_id: dict[int, str],
) -> pd.DataFrame:
    rows = connection.execute("""
        SELECT
            s.target_id,
            s.hgvs_pro,
            s.wt_mut,
            COUNT(DISTINCT s.score_id) AS score_sets,
            COUNT(DISTINCT COALESCE(m.experiment_urn, '')) AS experiments,
            COUNT(DISTINCT COALESCE(m.experiment_set_urn, '')) AS experiment_sets,
            COUNT(*) AS measurements
        FROM variant_scores AS s
        JOIN score_metadata AS m ON m.score_id = s.score_id
        GROUP BY s.target_id, s.hgvs_pro, s.wt_mut
        ORDER BY s.target_id, s.hgvs_pro, s.wt_mut
        """)
    result = pd.DataFrame.from_records(
        rows,
        columns=[
            "target_id",
            "hgvs_pro",
            "wt_mut",
            "score_sets",
            "experiments",
            "experiment_sets",
            "measurements",
        ],
    )
    if result.empty:
        return pd.DataFrame(columns=_VARIANT_COLUMNS)
    result["target_name"] = result["target_id"].map(target_name_by_id)
    result["wt_mut"] = result["wt_mut"].map(_bool_from_code)
    return result[_VARIANT_COLUMNS]


def _score_set_correlations(
    connection: sqlite3.Connection,
    group_label_by_id: dict[int, tuple[str, str | None]],
    score_ids: dict[str, int],
) -> pd.DataFrame:
    score_urn_by_id = {score_id: score_urn for score_urn, score_id in score_ids.items()}
    group_labels = connection.execute("""
        SELECT group_id, MIN(target_id) AS target_id
        FROM variant_scores
        GROUP BY group_id
        HAVING COUNT(DISTINCT score_id) > 1
        ORDER BY group_id
        """).fetchall()
    records: list[tuple[Any, ...]] = []

    for group_index, (group_id, target_id) in enumerate(group_labels, start=1):
        target_name, organism = group_label_by_id[group_id]
        frame = pd.read_sql_query(
            """
            SELECT score_id, hgvs_pro, AVG(score) AS score
            FROM variant_scores
            WHERE group_id = ?
            GROUP BY score_id, hgvs_pro
            """,
            connection,
            params=(group_id,),
        )
        if frame.empty:
            continue
        frame["score"] = frame["score"].replace([np.inf, -np.inf], np.nan)
        frame = frame.dropna(subset=["score"])
        if frame.empty:
            continue
        wide = frame.pivot(index="hgvs_pro", columns="score_id", values="score")
        ordered_ids = sorted(wide.columns.tolist())
        if len(ordered_ids) < 2:
            continue
        wide = wide[ordered_ids]
        pearson = wide.corr(method="pearson", min_periods=2).to_numpy(dtype=float)
        spearman = wide.corr(method="spearman", min_periods=2).to_numpy(dtype=float)

        for left_index in range(len(ordered_ids)):
            for right_index in range(left_index + 1, len(ordered_ids)):
                pearson_value = pearson[left_index, right_index]
                spearman_value = spearman[left_index, right_index]
                if not np.isfinite(pearson_value) or not np.isfinite(spearman_value):
                    continue
                set_1 = score_urn_by_id[ordered_ids[left_index]]
                set_2 = score_urn_by_id[ordered_ids[right_index]]
                records.append(
                    (
                        target_name,
                        organism,
                        set_1,
                        set_2,
                        float(spearman_value),
                        float(pearson_value),
                        _match_type(set_1, set_2),
                    )
                )
        if group_index % 100 == 0:
            LOGGER.info(
                "Calculated score-set correlations for %d target groups", group_index
            )

    result = pd.DataFrame.from_records(records, columns=_CORRELATION_COLUMNS)
    if result.empty:
        return result
    return result.sort_values(
        ["target_name", "organism", "set_1", "set_2"], kind="stable"
    ).reset_index(drop=True)


def _stage_singles(
    connection: sqlite3.Connection,
    path: Path,
    score_ids: dict[str, int],
    evo_score_urns: set[str],
    matrix: Any,
    chunk_size: int,
) -> None:
    required = {"score_urn", "wt", "mut", "score"}
    _validate_table_columns(path, "singles", required)
    allowed_amino_acids = set(str(matrix.alphabet))
    matrix_values = {
        (wild_type, mutant): float(matrix[wild_type, mutant])
        for wild_type in allowed_amino_acids
        for mutant in allowed_amino_acids
    }
    rows_seen = 0
    inserted_rows = 0
    LOGGER.info("Reading singles in chunks of %d rows: %s", chunk_size, path)
    for chunk_number, chunk in enumerate(
        pd.read_csv(
            path,
            sep="\t",
            usecols=sorted(required),
            dtype=str,
            keep_default_na=False,
            chunksize=chunk_size,
        ),
        start=1,
    ):
        chunk = chunk.reset_index(drop=True)
        rows_seen += len(chunk)
        score_values = pd.to_numeric(chunk["score"], errors="coerce").to_numpy(
            dtype=float
        )
        valid = np.isfinite(score_values)
        valid &= chunk["score_urn"].isin(evo_score_urns).to_numpy(dtype=bool)
        valid &= chunk["wt"].isin(allowed_amino_acids).to_numpy(dtype=bool)
        valid &= chunk["mut"].isin(allowed_amino_acids).to_numpy(dtype=bool)
        staged: list[tuple[int, float, float]] = []
        for row_index in np.flatnonzero(valid):
            score_urn = chunk.at[row_index, "score_urn"]
            wild_type = chunk.at[row_index, "wt"]
            mutant = chunk.at[row_index, "mut"]
            staged.append(
                (
                    score_ids[score_urn],
                    float(score_values[row_index]),
                    matrix_values[(wild_type, mutant)],
                )
            )
        if staged:
            with connection:
                connection.executemany(
                    "INSERT INTO evo_values (score_id, score, blosum100) VALUES (?, ?, ?)",
                    staged,
                )
            inserted_rows += len(staged)
        if chunk_number % 4 == 0:
            LOGGER.info("Read %d singles rows; staged %d", rows_seen, inserted_rows)
    LOGGER.info("Finished reading %d singles rows; staged %d", rows_seen, inserted_rows)


def _create_evolution_index(connection: sqlite3.Connection) -> None:
    with connection:
        connection.execute("CREATE INDEX idx_evo_values_score ON evo_values(score_id)")


def _evolution_correlations(
    connection: sqlite3.Connection,
    metadata: pd.DataFrame,
    score_ids: dict[str, int],
    evo_score_urns: set[str],
) -> pd.DataFrame:
    rows = connection.execute(
        "SELECT DISTINCT score_id FROM evo_values ORDER BY score_id"
    ).fetchall()
    score_urn_by_id = {score_id: score_urn for score_urn, score_id in score_ids.items()}
    metadata_by_score = metadata.set_index("score_urn", drop=False)
    records: list[dict[str, Any]] = []
    for (score_id,) in rows:
        score_urn = score_urn_by_id[score_id]
        if score_urn not in evo_score_urns:
            continue
        frame = pd.read_sql_query(
            "SELECT score, blosum100 FROM evo_values WHERE score_id = ?",
            connection,
            params=(score_id,),
        )
        if len(frame) <= 1:
            continue
        score_values = frame["score"].to_numpy(dtype=float)
        blosum_values = frame["blosum100"].to_numpy(dtype=float)
        if not np.isfinite(score_values).all() or not np.isfinite(blosum_values).all():
            continue
        if np.ptp(score_values) == 0 or np.ptp(blosum_values) == 0:
            continue
        correlation = stats.pearsonr(score_values, blosum_values)
        estimate = float(correlation.statistic)
        if len(frame) > 2:
            if abs(estimate) == 1:
                statistic = math.copysign(math.inf, estimate)
            else:
                statistic = estimate * math.sqrt((len(frame) - 2) / (1 - estimate**2))
        else:
            statistic = math.nan
        row = metadata_by_score.loc[score_urn]
        kingdom = row["kingdom"] if pd.notna(row["kingdom"]) else "Other"
        records.append(
            {
                "score_urn": score_urn,
                "organism": row["organism"],
                "kingdom": kingdom,
                "estimate": estimate,
                "statistic": statistic,
                "pvalue": float(correlation.pvalue),
                "n_positions": row["n_positions"],
                "n_subs": row["n_subs"],
            }
        )
    return pd.DataFrame.from_records(records, columns=_EVOLUTION_COLUMNS)


def _load_blosum_matrix(path: Path) -> Any:
    if not path.is_file():
        raise FileNotFoundError(f"BLOSUM100 matrix file does not exist: {path}")
    with path.open("r", encoding="utf-8") as matrix_file:
        matrix = substitution_matrices.read(matrix_file)
    if not matrix.alphabet:
        raise ValueError(f"BLOSUM100 matrix is empty: {path}")
    return matrix


def _validate_table_columns(path: Path, label: str, required: set[str]) -> None:
    columns = set(pd.read_csv(path, sep="\t", nrows=0, keep_default_na=False).columns)
    missing = sorted(required - columns)
    if missing:
        raise ValueError(
            f"{label} input is missing required columns: " + ", ".join(missing)
        )


def _mutation_type(row: pd.Series) -> str:
    has_multiple = bool(row["multiple_mutations"])
    has_single = bool(row["single_mutations"])
    if has_multiple and has_single:
        return "Both"
    if has_multiple:
        return "Multiple"
    if has_single:
        return "Single"
    return "Neither??"


def _match_type(set_1: str, set_2: str) -> str:
    experiment_1 = re.sub(r"-[a-z0-9]*$", "", set_1)
    experiment_2 = re.sub(r"-[a-z0-9]*$", "", set_2)
    experiment_set_1 = re.sub(r"-[a-z0-9]*-[a-z0-9]*$", "", set_1)
    experiment_set_2 = re.sub(r"-[a-z0-9]*-[a-z0-9]*$", "", set_2)
    if experiment_1 == experiment_2:
        return "Same Experiment"
    if experiment_set_1 == experiment_set_2:
        return "Same Study"
    return "Same Gene"


def _parse_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    normalized = str(value).strip().lower()
    return normalized in {"true", "t", "1", "yes"}


def _numeric_bool(value: Any) -> float:
    normalized = _global_text(value)
    if normalized is None:
        return math.nan
    lowered = normalized.lower()
    if lowered in {"true", "t", "yes"}:
        return 1.0
    if lowered in {"false", "f", "no"}:
        return 0.0
    try:
        return float(normalized)
    except ValueError:
        return math.nan


def _bool_code(value: Any) -> int | None:
    normalized = _global_text(value)
    if normalized is None:
        return None
    lowered = normalized.lower()
    if lowered in {"true", "t", "1", "yes"}:
        return 1
    if lowered in {"false", "f", "0", "no"}:
        return 0
    return None


def _bool_from_code(value: Any) -> bool | None:
    if value is None or pd.isna(value):
        return None
    return bool(value)


def _is_protein_hgvs(value: Any) -> bool:
    text = _global_text(value)
    return text is not None and text != "NA"


def _identifier(value: Any) -> str | None:
    return _global_text(value)


def _plain_text(value: Any) -> str | None:
    return _global_text(value)


def _global_text(value: Any, *, preserve_literal_na: bool = False) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "null":
        return None
    if text == "NA" and not preserve_literal_na:
        return None
    return text


def _scalar_json_text(value: Any, *, preserve_literal_na: bool = False) -> str | None:
    text = _global_text(value, preserve_literal_na=preserve_literal_na)
    if text is None:
        return None
    if text.startswith(("[", "{")):
        try:
            decoded = json.loads(text)
        except json.JSONDecodeError:
            return text
        if isinstance(decoded, list):
            if len(decoded) != 1:
                return None
            decoded = decoded[0]
        if isinstance(decoded, dict) or decoded is None:
            return None
        text = str(decoded).strip()
        if not text or text.lower() == "null":
            return None
        if text == "NA" and not preserve_literal_na:
            return None
    return text


def _is_multi_json(value: Any) -> bool:
    if value is None or pd.isna(value):
        return False
    text = str(value).strip()
    if not text.startswith("["):
        return False
    try:
        decoded = json.loads(text)
    except json.JSONDecodeError:
        return False
    return isinstance(decoded, list) and len(decoded) > 1


def _sqlite_text(value: Any) -> str | None:
    return None if value is None or pd.isna(value) else str(value)
