"""Chunked, identity-preserving input preparation for MaveDB GMM fitting."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import sqlite3
import tempfile
from collections import Counter
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

_CATALOGUE_SCORE_COLUMNS = (
    "dataset_id",
    "source_row",
    "hgvs_nt",
    "hgvs_pro",
    "score",
)
_CATALOGUE_META_COLUMNS = (
    "dataset_id",
    "meta_study",
    "analysis_types",
    "target_name",
    "target_category",
    "taxid",
)
_LEGACY_SCORE_COLUMNS = ("score_urn", "hgvs_nt", "hgvs_pro", "score")
_LEGACY_META_COLUMNS = ("score_urn", "meta_study")
_EXCLUSION_COLUMNS = ("study_id", "source_row", "input_row", "reason")
_ELIGIBILITY_COLUMNS = ("study_id", "reason", "details")
_MISSING_TEXT = {"", "na", "null"}


class InputDataError(ValueError):
    """Raised when an input table violates the declared fitting contract."""


@dataclass
class _MetadataStatus:
    reason: str | None
    details: str


@dataclass
class PreparedInputs:
    """Selected studies staged in scratch space and removed by ``close``."""

    study_counts: dict[str, int]
    study_qc: dict[str, dict[str, int]]
    exclusions_path: Path
    eligibility_path: Path
    _database_path: Path = field(repr=False)
    _input_format: str = field(repr=False)
    _staged_studies: tuple[str, ...] = field(repr=False)
    _temporary_directory: tempfile.TemporaryDirectory = field(repr=False)
    _closed: bool = field(default=False, init=False, repr=False)

    def __enter__(self) -> PreparedInputs:
        self._ensure_open()
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def close(self) -> None:
        """Remove staged input and QC files; safe to call more than once."""

        if not self._closed:
            self._temporary_directory.cleanup()
            self._closed = True

    def iter_studies(
        self, study_ids: Sequence[str]
    ) -> Iterator[tuple[str, pd.DataFrame]]:
        """Yield requested selected studies in caller order, one frame at a time."""

        self._ensure_open()
        requested = tuple(study_ids)
        if len(set(requested)) != len(requested):
            raise InputDataError("study_ids must not contain duplicates")
        unknown = [
            study_id for study_id in requested if study_id not in self._staged_studies
        ]
        if unknown:
            raise InputDataError(
                "requested study was not staged: "
                + ", ".join(repr(study_id) for study_id in unknown)
            )

        connection = sqlite3.connect(self._database_path)
        try:
            for study_id in requested:
                rows = connection.execute(
                    """
                    SELECT dataset_id, source_row, hgvs_nt, hgvs_pro, score_urn, score
                    FROM staged_rows
                    WHERE study_id = ?
                    ORDER BY source_order
                    """,
                    (study_id,),
                ).fetchall()
                if self._input_format == "catalogue":
                    frame = _catalogue_frame(rows)
                else:
                    frame = _legacy_frame(rows)
                yield study_id, frame
                del frame, rows
        finally:
            connection.close()

    def _ensure_open(self) -> None:
        if self._closed:
            raise RuntimeError("prepared inputs are closed")


def prepare_inputs(
    scores: Path,
    meta: Path,
    input_format: str = "catalogue",
    study_ids: tuple[str, ...] = (),
    max_studies: int | None = None,
    max_study_rows: int | None = None,
    seed: int | None = None,
    chunk_size: int = 5000,
    stage_selector: Callable[[dict[str, int]], Sequence[str]] | None = None,
    counts_only: bool = False,
) -> PreparedInputs:
    """Validate, select, and stage GMM inputs without loading full score tables.

    Catalogue rows are unique by ``(dataset_id, source_row)``. Eligibility is
    determined from metadata before row-level checks; valid rows require a
    finite numeric score and protein HGVS. Missing nucleotide HGVS is allowed.
    When a row cap and explicit seed are set, each study gets an independent
    deterministic reservoir sample. Without a seed, the random generator uses
    system entropy. Selected rows are returned in source order.
    """

    _validate_options(
        input_format,
        study_ids,
        max_studies,
        max_study_rows,
        seed,
        chunk_size,
        stage_selector,
        counts_only,
    )
    scores_path = _input_path(scores, "scores")
    meta_path = _input_path(meta, "metadata")
    score_columns = (
        _CATALOGUE_SCORE_COLUMNS
        if input_format == "catalogue"
        else _LEGACY_SCORE_COLUMNS
    )
    metadata_columns = (
        _CATALOGUE_META_COLUMNS if input_format == "catalogue" else _LEGACY_META_COLUMNS
    )
    _validate_headers(scores_path, score_columns)
    _validate_headers(meta_path, metadata_columns)

    temporary_directory = tempfile.TemporaryDirectory(prefix="mavedb-gmm-inputs-")
    succeeded = False
    connection: sqlite3.Connection | None = None
    try:
        scratch = Path(temporary_directory.name)
        database_path = scratch / "staged.sqlite3"
        exclusions_path = scratch / "exclusions.tsv"
        eligibility_path = scratch / "eligibility.tsv"
        metadata_status = _read_metadata(
            meta_path, metadata_columns, input_format, chunk_size
        )
        _write_eligibility(eligibility_path, metadata_status)
        _validate_requested_studies(study_ids, metadata_status)

        connection = _create_stage_database(database_path)
        raw_counts: Counter[str] = Counter()
        eligible_counts: Counter[str] = Counter()
        exclusion_counts: dict[str, Counter[str]] = {}
        with exclusions_path.open("w", newline="", encoding="utf-8") as handle:
            exclusion_writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            exclusion_writer.writerow(_EXCLUSION_COLUMNS)
            _count_and_validate_rows(
                scores_path,
                score_columns,
                metadata_status,
                input_format,
                chunk_size,
                connection,
                raw_counts,
                eligible_counts,
                exclusion_counts,
                exclusion_writer,
            )

            selected_studies = _select_studies(
                study_ids, max_studies, metadata_status, raw_counts, eligible_counts
            )
            study_counts = {
                study_id: eligible_counts[study_id] for study_id in selected_studies
            }
            if counts_only:
                staged_studies: list[str] = []
            elif stage_selector is None:
                staged_studies = list(selected_studies)
            else:
                staged_studies = _validate_stage_selection(
                    stage_selector(dict(study_counts)), selected_studies
                )
                if not staged_studies:
                    raise InputDataError("stage_selector selected no studies")
            if staged_studies:
                _stage_selected_rows(
                    scores_path,
                    score_columns,
                    metadata_status,
                    input_format,
                    chunk_size,
                    selected_studies,
                    staged_studies,
                    max_study_rows,
                    seed,
                    connection,
                    raw_counts,
                    eligible_counts,
                    exclusion_writer,
                )

        connection.execute(
            "CREATE INDEX staged_rows_by_study ON staged_rows(study_id, source_order)"
        )
        connection.commit()
        study_qc = _build_study_qc(
            staged_studies,
            raw_counts,
            eligible_counts,
            exclusion_counts,
            max_study_rows,
            connection,
            input_format,
        )
        for study_id, expected_rows in study_qc.items():
            actual_rows = connection.execute(
                "SELECT COUNT(*) FROM staged_rows WHERE study_id = ?", (study_id,)
            ).fetchone()[0]
            if actual_rows != expected_rows["selected_rows"]:
                raise InputDataError(
                    f"staged row count mismatch for {study_id!r}: "
                    f"expected {expected_rows['selected_rows']}, found {actual_rows}"
                )
        connection.close()
        connection = None

        prepared = PreparedInputs(
            study_counts=study_counts,
            study_qc=study_qc,
            exclusions_path=exclusions_path,
            eligibility_path=eligibility_path,
            _database_path=database_path,
            _input_format=input_format,
            _staged_studies=tuple(staged_studies),
            _temporary_directory=temporary_directory,
        )
        succeeded = True
        return prepared
    finally:
        if connection is not None:
            connection.close()
        if not succeeded:
            temporary_directory.cleanup()


def _validate_options(
    input_format: str,
    study_ids: Sequence[str],
    max_studies: int | None,
    max_study_rows: int | None,
    seed: int | None,
    chunk_size: int,
    stage_selector: Callable[[dict[str, int]], Sequence[str]] | None,
    counts_only: bool,
) -> None:
    if input_format not in {"catalogue", "legacy"}:
        raise InputDataError("input_format must be 'catalogue' or 'legacy'")
    if isinstance(study_ids, str):
        raise InputDataError("study_ids must be a sequence of study identifiers")
    if any(not isinstance(study_id, str) or study_id == "" for study_id in study_ids):
        raise InputDataError("study_ids must contain nonempty strings")
    if len(set(study_ids)) != len(study_ids):
        raise InputDataError("study_ids must not contain duplicates")
    _validate_positive_integer("max_studies", max_studies, allow_none=True)
    _validate_positive_integer("max_study_rows", max_study_rows, allow_none=True)
    _validate_positive_integer("chunk_size", chunk_size)
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise InputDataError("seed must be an integer or None")
    if stage_selector is not None and not callable(stage_selector):
        raise InputDataError("stage_selector must be callable or None")
    if not isinstance(counts_only, bool):
        raise InputDataError("counts_only must be a boolean")
    if counts_only and stage_selector is not None:
        raise InputDataError("counts_only cannot be combined with stage_selector")
    if study_ids and max_studies is not None and len(study_ids) > max_studies:
        raise InputDataError("max_studies cannot truncate explicit study_ids")


def _validate_positive_integer(
    name: str, value: int | None, allow_none: bool = False
) -> None:
    if value is None and allow_none:
        return
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise InputDataError(f"{name} must be a positive integer")


def _input_path(path: Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} input is not a file: {resolved}")
    return resolved


def _validate_headers(path: Path, required_columns: Sequence[str]) -> None:
    try:
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            header = next(csv.reader(handle, delimiter="\t"), None)
    except UnicodeDecodeError as error:
        raise InputDataError(f"{path} is not valid UTF-8 TSV") from error
    if not header:
        raise InputDataError(f"{path} is empty or has no TSV header")
    if any(column == "" for column in header):
        raise InputDataError(f"{path} has an empty column name")
    if len(set(header)) != len(header):
        raise InputDataError(f"{path} has duplicate column names")
    missing = [column for column in required_columns if column not in header]
    if missing:
        raise InputDataError(
            f"{path} is missing required columns: {', '.join(missing)}"
        )


def _iter_table_chunks(
    path: Path, columns: Sequence[str], chunk_size: int
) -> Iterator[pd.DataFrame]:
    try:
        with pd.read_csv(
            path,
            sep="\t",
            usecols=list(columns),
            dtype=str,
            keep_default_na=False,
            na_filter=False,
            chunksize=chunk_size,
            encoding="utf-8-sig",
        ) as reader:
            for chunk in reader:
                yield chunk.loc[:, list(columns)]
    except pd.errors.ParserError as error:
        raise InputDataError(f"could not parse TSV input {path}: {error}") from error
    except UnicodeDecodeError as error:
        raise InputDataError(f"{path} is not valid UTF-8 TSV") from error


def _read_metadata(
    path: Path, columns: Sequence[str], input_format: str, chunk_size: int
) -> dict[str, _MetadataStatus]:
    key_column = "dataset_id" if input_format == "catalogue" else "score_urn"
    metadata: dict[str, _MetadataStatus] = {}
    for chunk in _iter_table_chunks(path, columns, chunk_size):
        for row in chunk.to_dict(orient="records"):
            study_id = row[key_column]
            if study_id == "":
                raise InputDataError(f"empty metadata key in {path}")
            if study_id in metadata:
                raise InputDataError(f"duplicate metadata key {study_id!r} in {path}")
            status = _metadata_status(row, input_format)
            metadata[study_id] = status
    return metadata


def _metadata_status(row: dict[str, str], input_format: str) -> _MetadataStatus:
    meta_study = _parse_bool(row["meta_study"])
    if input_format == "legacy":
        if meta_study is True:
            return _MetadataStatus("meta_study", "meta_study=TRUE")
        if meta_study is None:
            return _MetadataStatus(
                "ambiguous_meta_study", "meta_study must be TRUE or FALSE"
            )
        return _MetadataStatus(None, "eligible legacy score set")

    analysis_types = _parse_analysis_types(row["analysis_types"])
    target_fields = {
        name: _parse_target_values(row[name], name)
        for name in ("target_name", "target_category", "taxid")
    }
    target_cardinalities = {name: len(values) for name, values in target_fields.items()}
    reasons: list[str] = []
    details: list[str] = []

    if meta_study is True:
        reasons.append("meta_study")
        details.append("meta_study=TRUE")
    elif meta_study is None:
        reasons.append("ambiguous_meta_study")
        details.append("meta_study must be TRUE or FALSE")

    if analysis_types is None:
        reasons.append("missing_analysis_types")
        details.append("analysis_types is missing")
    elif "maves" not in analysis_types:
        reasons.append("not_maves")
        details.append("analysis_types does not include maves")

    if any(length > 1 for length in target_cardinalities.values()):
        reasons.append("multi_target")
        details.append(f"target field cardinalities={target_cardinalities}")
    else:
        for name, values in target_fields.items():
            if not values or values[0] is None:
                reason = f"missing_{name}"
                reasons.append(reason)
                details.append(f"{name} is missing")
        categories = target_fields["target_category"]
        if categories and categories[0] is not None:
            if categories[0] != "protein_coding":
                reasons.append("unsupported_target_category")
                details.append(f"target_category={categories[0]!r}")
        taxids = target_fields["taxid"]
        if taxids and taxids[0] is not None:
            taxid = _parse_taxid(taxids[0])
            if taxid is None:
                reasons.append("ambiguous_taxid")
                details.append(f"taxid is not an unambiguous integer: {taxids[0]!r}")
            elif taxid != 9606:
                reasons.append("nonhuman_target")
                details.append(f"taxid={taxid}")

    if reasons:
        return _MetadataStatus(";".join(reasons), "; ".join(details))
    return _MetadataStatus(None, "single human MAVE target")


def _parse_bool(value: str) -> bool | None:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    return None


def _parse_analysis_types(value: str) -> list[str] | None:
    normalized = value.strip()
    if normalized.lower() in _MISSING_TEXT:
        return None
    try:
        parsed = json.loads(normalized)
    except json.JSONDecodeError as error:
        raise InputDataError(f"analysis_types is malformed JSON: {value!r}") from error
    if not isinstance(parsed, list) or any(
        not isinstance(item, str) for item in parsed
    ):
        raise InputDataError("analysis_types must be a JSON array of strings")
    return parsed


def _parse_target_values(value: str, field_name: str) -> list[str | None]:
    normalized = value.strip()
    if normalized.lower() in _MISSING_TEXT:
        return []
    if normalized.startswith("[") or normalized.lower() == "null":
        try:
            parsed = json.loads(normalized)
        except json.JSONDecodeError as error:
            raise InputDataError(
                f"{field_name} is malformed JSON: {value!r}"
            ) from error
        if parsed is None:
            return []
        if not isinstance(parsed, list):
            raise InputDataError(f"{field_name} must be a scalar or JSON array")
        values = parsed
    elif normalized.startswith("{"):
        raise InputDataError(f"{field_name} must be a scalar or JSON array")
    else:
        values = [normalized]

    converted: list[str | None] = []
    for item in values:
        if item is None:
            converted.append(None)
        elif isinstance(item, (str, int, float)) and not isinstance(item, bool):
            text = str(item)
            converted.append(None if text.strip().lower() in _MISSING_TEXT else text)
        else:
            raise InputDataError(f"{field_name} array values must be scalar")
    return converted


def _parse_taxid(value: str) -> int | None:
    try:
        parsed = int(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def _write_eligibility(path: Path, metadata_status: dict[str, _MetadataStatus]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(_ELIGIBILITY_COLUMNS)
        for study_id, status in metadata_status.items():
            if status.reason is not None:
                writer.writerow((study_id, status.reason, status.details))


def _validate_requested_studies(
    study_ids: Sequence[str], metadata_status: dict[str, _MetadataStatus]
) -> None:
    for study_id in study_ids:
        if study_id not in metadata_status:
            raise InputDataError(
                f"requested study {study_id!r} is absent from metadata"
            )
        reason = metadata_status[study_id].reason
        if reason is not None:
            raise InputDataError(
                f"requested study {study_id!r} is metadata-ineligible: {reason}"
            )


def _create_stage_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("""
        CREATE TABLE observation_keys (
            study_id TEXT NOT NULL,
            source_row TEXT NOT NULL,
            PRIMARY KEY (study_id, source_row)
        )
        """)
    connection.execute("""
        CREATE TABLE staged_rows (
            study_id TEXT NOT NULL,
            source_order INTEGER NOT NULL,
            dataset_id TEXT,
            source_row TEXT,
            score_urn TEXT,
            hgvs_nt TEXT,
            hgvs_pro TEXT,
            score REAL NOT NULL,
            PRIMARY KEY (study_id, source_order)
        )
        """)
    connection.execute("""
        CREATE TABLE reservoir_rows (
            study_id TEXT NOT NULL,
            slot INTEGER NOT NULL,
            source_order INTEGER NOT NULL,
            dataset_id TEXT,
            source_row TEXT,
            score_urn TEXT,
            hgvs_nt TEXT,
            hgvs_pro TEXT,
            score REAL NOT NULL,
            PRIMARY KEY (study_id, slot)
        )
        """)
    return connection


def _count_and_validate_rows(
    scores_path: Path,
    columns: Sequence[str],
    metadata_status: dict[str, _MetadataStatus],
    input_format: str,
    chunk_size: int,
    connection: sqlite3.Connection,
    raw_counts: Counter[str],
    eligible_counts: Counter[str],
    exclusion_counts: dict[str, Counter[str]],
    exclusion_writer: csv.writer,
) -> None:
    study_column = "dataset_id" if input_format == "catalogue" else "score_urn"
    source_order = 0
    for chunk in _iter_table_chunks(scores_path, columns, chunk_size):
        for row in chunk.to_dict(orient="records"):
            source_order += 1
            study_id = row[study_column]
            if study_id == "":
                raise InputDataError(
                    f"empty {study_column} in scores record {source_order}"
                )
            if study_id not in metadata_status:
                raise InputDataError(
                    f"score study {study_id!r} has no matching metadata record"
                )
            source_row = row.get("source_row")
            if input_format == "catalogue":
                _validate_source_row(source_row, source_order)
                try:
                    connection.execute(
                        "INSERT INTO observation_keys(study_id, source_row) VALUES (?, ?)",
                        (study_id, source_row),
                    )
                except sqlite3.IntegrityError as error:
                    raise InputDataError(
                        f"duplicate catalogue observation key ({study_id!r}, {source_row!r})"
                    ) from error
            raw_counts[study_id] += 1
            status = metadata_status[study_id]
            if status.reason is not None:
                _record_exclusion(
                    exclusion_writer,
                    input_format,
                    study_id,
                    source_row,
                    source_order,
                    f"metadata:{status.reason}",
                )
                continue
            _, reasons = _score_row_status(row)
            if reasons:
                for reason in reasons:
                    exclusion_counts.setdefault(study_id, Counter())[reason] += 1
                _record_exclusion(
                    exclusion_writer,
                    input_format,
                    study_id,
                    source_row,
                    source_order,
                    ";".join(reasons),
                )
                continue
            eligible_counts[study_id] += 1
        connection.commit()


def _validate_source_row(source_row: str | None, source_order: int) -> None:
    valid = source_row is not None and source_row.isascii() and source_row.isdigit()
    if valid:
        try:
            value = int(source_row)
        except ValueError:
            valid = False
        else:
            valid = value > 0 and str(value) == source_row
    if not valid:
        raise InputDataError(
            f"source_row must be a canonical positive integer at scores record {source_order}"
        )


def _score_row_status(row: dict[str, str]) -> tuple[float | None, list[str]]:
    score_text = row["score"].strip()
    reasons: list[str] = []
    score: float | None = None
    if score_text.lower() in _MISSING_TEXT:
        reasons.append("missing_score")
    else:
        try:
            score = float(score_text)
        except ValueError:
            reasons.append("invalid_score")
        else:
            if not math.isfinite(score):
                reasons.append("nonfinite_score")

    if row["hgvs_pro"].strip().lower() in _MISSING_TEXT:
        reasons.append("missing_hgvs_pro")
    return score, reasons


def _select_studies(
    requested: Sequence[str],
    max_studies: int | None,
    metadata_status: dict[str, _MetadataStatus],
    raw_counts: Counter[str],
    eligible_counts: Counter[str],
) -> list[str]:
    if requested:
        selected = list(requested)
        for study_id in selected:
            if raw_counts.get(study_id, 0) == 0:
                raise InputDataError(f"requested study {study_id!r} has no score rows")
        return selected

    eligible_studies = sorted(
        study_id
        for study_id, status in metadata_status.items()
        if status.reason is None and raw_counts.get(study_id, 0) > 0
    )
    eligible_studies.sort(
        key=lambda study_id: (eligible_counts.get(study_id, 0) == 0, study_id)
    )
    if max_studies is not None:
        eligible_studies = eligible_studies[:max_studies]
    if not eligible_studies:
        raise InputDataError("no eligible studies selected")
    return eligible_studies


def _validate_stage_selection(
    staged_studies: Sequence[str], planned_studies: Sequence[str]
) -> list[str]:
    if isinstance(staged_studies, str):
        raise InputDataError("stage_selector must return a sequence of study IDs")
    staged = list(staged_studies)
    if any(not isinstance(study_id, str) for study_id in staged):
        raise InputDataError("stage_selector returned a non-string study ID")
    if len(set(staged)) != len(staged):
        raise InputDataError("stage_selector returned duplicate study IDs")
    planned = set(planned_studies)
    unknown = [study_id for study_id in staged if study_id not in planned]
    if unknown:
        raise InputDataError(
            "stage_selector returned IDs outside the planned selection: "
            + ", ".join(repr(study_id) for study_id in unknown)
        )
    return staged


def _stage_selected_rows(
    scores_path: Path,
    columns: Sequence[str],
    metadata_status: dict[str, _MetadataStatus],
    input_format: str,
    chunk_size: int,
    planned_studies: Sequence[str],
    staged_studies: Sequence[str],
    max_study_rows: int | None,
    seed: int | None,
    connection: sqlite3.Connection,
    raw_counts: Counter[str],
    eligible_counts: Counter[str],
    exclusion_writer: csv.writer,
) -> None:
    study_column = "dataset_id" if input_format == "catalogue" else "score_urn"
    planned = set(planned_studies)
    staged = set(staged_studies)
    rngs = {
        study_id: random.Random(_study_seed(seed, study_id))
        for study_id in staged_studies
    }
    selected_valid_counts: Counter[str] = Counter()
    second_raw_counts: Counter[str] = Counter()
    source_order = 0
    for chunk in _iter_table_chunks(scores_path, columns, chunk_size):
        for row in chunk.to_dict(orient="records"):
            source_order += 1
            study_id = row[study_column]
            second_raw_counts[study_id] += 1
            status = metadata_status[study_id]
            if status.reason is not None:
                continue
            score, reasons = _score_row_status(row)
            if reasons:
                continue
            selected_valid_counts[study_id] += 1
            if study_id not in planned:
                _record_exclusion(
                    exclusion_writer,
                    input_format,
                    study_id,
                    row.get("source_row"),
                    source_order,
                    "study_limit",
                )
                continue
            if study_id not in staged:
                continue

            values = _stage_values(study_id, source_order, row, input_format, score)
            if max_study_rows is None:
                connection.execute(
                    """
                    INSERT INTO staged_rows(
                        study_id, source_order, dataset_id, source_row,
                        score_urn, hgvs_nt, hgvs_pro, score
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    values,
                )
                continue

            count = selected_valid_counts[study_id]
            if count <= max_study_rows:
                slot = count - 1
                connection.execute(
                    """
                    INSERT INTO reservoir_rows(
                        study_id, slot, source_order, dataset_id, source_row,
                        score_urn, hgvs_nt, hgvs_pro, score
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (study_id, slot, *values[1:]),
                )
                continue

            slot = rngs[study_id].randrange(count)
            if slot >= max_study_rows:
                _record_exclusion(
                    exclusion_writer,
                    input_format,
                    study_id,
                    row.get("source_row"),
                    source_order,
                    "study_row_cap",
                )
                continue
            replaced = connection.execute(
                """
                SELECT source_order, source_row
                FROM reservoir_rows
                WHERE study_id = ? AND slot = ?
                """,
                (study_id, slot),
            ).fetchone()
            if replaced is None:
                raise InputDataError(
                    f"missing reservoir slot {slot} for study {study_id!r}"
                )
            _record_exclusion(
                exclusion_writer,
                input_format,
                study_id,
                replaced[1],
                replaced[0],
                "study_row_cap",
            )
            connection.execute(
                """
                INSERT OR REPLACE INTO reservoir_rows(
                    study_id, slot, source_order, dataset_id, source_row,
                    score_urn, hgvs_nt, hgvs_pro, score
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (study_id, slot, *values[1:]),
            )
        connection.commit()

    if dict(second_raw_counts) != dict(raw_counts):
        raise InputDataError("scores input changed while it was being prepared")
    if dict(selected_valid_counts) != dict(eligible_counts):
        raise InputDataError(
            "eligible score rows changed while inputs were being prepared"
        )

    if max_study_rows is not None:
        for study_id in staged_studies:
            connection.execute(
                """
                INSERT INTO staged_rows(
                    study_id, source_order, dataset_id, source_row,
                    score_urn, hgvs_nt, hgvs_pro, score
                )
                SELECT study_id, source_order, dataset_id, source_row,
                       score_urn, hgvs_nt, hgvs_pro, score
                FROM reservoir_rows
                WHERE study_id = ?
                """,
                (study_id,),
            )
    connection.commit()


def _stage_values(
    study_id: str,
    source_order: int,
    row: dict[str, str],
    input_format: str,
    score: float | None,
) -> tuple[Any, ...]:
    if score is None:
        raise InputDataError(
            f"internal error: selected score is missing for {study_id!r}"
        )
    if input_format == "catalogue":
        return (
            study_id,
            source_order,
            row["dataset_id"],
            row["source_row"],
            study_id,
            row["hgvs_nt"],
            row["hgvs_pro"],
            score,
        )
    return (
        study_id,
        source_order,
        None,
        None,
        study_id,
        row["hgvs_nt"],
        row["hgvs_pro"],
        score,
    )


def _study_seed(seed: int | None, study_id: str) -> int | None:
    if seed is None:
        return None
    seed_text = str(seed)
    digest = hashlib.sha256(f"{seed_text}\0{study_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False)


def _record_exclusion(
    writer: csv.writer,
    input_format: str,
    study_id: str,
    source_row: str | None,
    source_order: int,
    reason: str,
) -> None:
    if input_format == "catalogue":
        writer.writerow((study_id, source_row, "", reason))
    else:
        writer.writerow((study_id, "", source_order, reason))


def _build_study_qc(
    selected_studies: Sequence[str],
    raw_counts: Counter[str],
    eligible_counts: Counter[str],
    reason_counts: dict[str, Counter[str]],
    max_study_rows: int | None,
    connection: sqlite3.Connection,
    input_format: str,
) -> dict[str, dict[str, int]]:
    qc: dict[str, dict[str, int]] = {}
    for study_id in selected_studies:
        raw = raw_counts[study_id]
        eligible = eligible_counts[study_id]
        selected = (
            min(eligible, max_study_rows) if max_study_rows is not None else eligible
        )
        counts = reason_counts.get(study_id, Counter())
        selected_text_bytes = _selected_text_bytes(connection, study_id, input_format)
        qc[study_id] = {
            "raw_rows": raw,
            "eligible_rows": eligible,
            "selected_rows": selected,
            "excluded_rows": raw - selected,
            "invalid_rows": raw - eligible,
            "subset_excluded_rows": eligible - selected,
            "missing_score_rows": counts["missing_score"],
            "invalid_score_rows": counts["invalid_score"],
            "nonfinite_score_rows": counts["nonfinite_score"],
            "missing_hgvs_pro_rows": counts["missing_hgvs_pro"],
            "selected_text_bytes": selected_text_bytes,
            "staged_bytes_estimate": selected_text_bytes + selected * 512,
        }
    return qc


def _selected_text_bytes(
    connection: sqlite3.Connection, study_id: str, input_format: str
) -> int:
    if input_format == "catalogue":
        query = """
            SELECT dataset_id, source_row, hgvs_nt, hgvs_pro, score_urn, score
            FROM staged_rows
            WHERE study_id = ?
        """
        columns = 6
    else:
        query = """
            SELECT score_urn, hgvs_nt, hgvs_pro, score
            FROM staged_rows
            WHERE study_id = ?
        """
        columns = 4
    byte_count = 0
    for row in connection.execute(query, (study_id,)):
        text_values = row[:-1]
        byte_count += sum(
            len(value.encode("utf-8")) for value in text_values if value is not None
        )
        byte_count += len(format(float(row[-1]), ".17g").encode("ascii"))
        byte_count += columns
    return byte_count


def _catalogue_frame(rows: list[tuple[Any, ...]]) -> pd.DataFrame:
    dataset_ids = [row[0] for row in rows]
    source_rows = [row[1] for row in rows]
    frame = pd.DataFrame(
        {
            "dataset_id": pd.array(dataset_ids, dtype="string"),
            "source_row": pd.array(source_rows, dtype="string"),
            "hgvs_nt": _nullable_hgvs([row[2] for row in rows]),
            "hgvs_pro": _nullable_hgvs([row[3] for row in rows]),
            "score": pd.Series([row[5] for row in rows], dtype="float64"),
            "score_urn": pd.array(dataset_ids, dtype="string"),
        }
    )
    return frame.loc[:, [*_CATALOGUE_SCORE_COLUMNS, "score_urn"]]


def _legacy_frame(rows: list[tuple[Any, ...]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "score_urn": pd.array([row[4] for row in rows], dtype="string"),
            "hgvs_nt": _nullable_hgvs([row[2] for row in rows]),
            "hgvs_pro": _nullable_hgvs([row[3] for row in rows]),
            "score": pd.Series([row[5] for row in rows], dtype="float64"),
        }
    )


def _nullable_hgvs(values: list[str | None]) -> pd.arrays.StringArray:
    normalized = [
        pd.NA if value is None or value.strip().lower() in _MISSING_TEXT else value
        for value in values
    ]
    return pd.array(normalized, dtype="string")
