"""Streaming helpers and Open Targets coding-variant reference exports."""

from __future__ import annotations

import csv
import json
import math
import sqlite3
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Sequence

import pandas as pd
import pyarrow.dataset as ds

CODING_COLUMNS = (
    "targetId",
    "variantId",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
)
CLINVAR_COLUMNS = (
    "targetId",
    "variantId",
    "variantHgvsId",
    "variantFunctionalConsequenceId",
    "clinicalSignificances",
    "score",
    "variantEffect",
    "directionOnTrait",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
)
CLINVAR_EVIDENCE_REQUIRED_COLUMNS = (
    "targetId",
    "variantId",
    "clinicalSignificances",
)
GWAS_INPUT_COLUMNS = (
    "studyLocusId",
    "studyId",
    "variantId",
    "chromosome",
    "position",
    "region",
    "beta",
    "zScore",
    "pValueMantissa",
    "pValueExponent",
)
GWAS_COLUMNS = (
    *GWAS_INPUT_COLUMNS,
    "targetId",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
)


def _source_paths(paths: Path | Sequence[Path]) -> str | list[str]:
    if isinstance(paths, Path):
        return str(paths)
    source = [str(path) for path in paths]
    if not source:
        raise ValueError("At least one Parquet path is required")
    return source


def validate_batch_size(batch_size: int) -> None:
    if (
        isinstance(batch_size, bool)
        or not isinstance(batch_size, int)
        or batch_size < 1
    ):
        raise ValueError("batch_size must be a positive integer")


def _dataset_batches(
    paths: Path | Sequence[Path],
    columns: Sequence[str],
    batch_size: int,
    *,
    required_columns: Sequence[str] | None = None,
    filter_expression: ds.Expression | None = None,
) -> Iterator[list[dict[str, Any]]]:
    validate_batch_size(batch_size)
    if not columns:
        raise ValueError("At least one Parquet column is required")
    if len(set(columns)) != len(columns):
        raise ValueError("Parquet columns must be unique")
    required = set(columns if required_columns is None else required_columns)
    if not required.issubset(columns):
        raise ValueError("Required Parquet columns must be selected columns")

    dataset = ds.dataset(
        _source_paths(paths), format="parquet", exclude_invalid_files=False
    )
    missing = sorted(required - set(dataset.schema.names))
    if missing:
        label = "requested" if required_columns is None else "required"
        raise ValueError(f"Parquet dataset is missing {label} columns: {missing}")
    selected = [column for column in columns if column in dataset.schema.names]
    if not selected:
        raise ValueError("Parquet dataset has none of the requested columns")

    scanner = dataset.scanner(
        columns=selected,
        batch_size=batch_size,
        filter=filter_expression,
    )
    for batch in scanner.to_batches():
        yield batch.to_pylist()


def parquet_batches(
    paths: Path | Sequence[Path],
    columns: Sequence[str],
    *,
    batch_size: int = 65_536,
    required_columns: Sequence[str] | None = None,
    filter_column: str | None = None,
    filter_values: Sequence[Any] | None = None,
) -> Iterator[list[dict[str, Any]]]:
    """Yield projected Parquet batches with optional columns and row filtering.

    All ``columns`` are required by default. Pass ``required_columns`` when some
    projected fields may be absent from a source release; those missing fields
    are omitted from each row and can be read with ``row.get``. A filter column
    and its accepted values must be supplied together.
    """
    if (filter_column is None) != (filter_values is None):
        raise ValueError("filter_column and filter_values must be supplied together")

    filter_expression = None
    if filter_column is not None:
        filter_expression = ds.field(filter_column).isin(list(filter_values))
    yield from _dataset_batches(
        paths,
        columns,
        batch_size,
        required_columns=required_columns,
        filter_expression=filter_expression,
    )


def explode(value: Any) -> list[Any]:
    """Match tidyr::unnest(): drop null/empty lists and keep scalar values."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


@contextmanager
def coding_index(
    coding: Path | Sequence[Path],
    directory: Path,
    *,
    batch_size: int = 65_536,
    genes: Sequence[str] | None = None,
) -> Iterator[sqlite3.Connection]:
    """Build a disk-backed coding-variant index and yield its open connection.

    Rows are exploded on ``uniprotAccessions`` before insertion. Duplicate
    source rows and duplicate accessions are kept, while null and empty lists
    contribute no rows. ``genes`` applies a source-side ``targetId`` filter.
    """
    validate_batch_size(batch_size)
    directory.mkdir(parents=True, exist_ok=True)
    temp_dir = tempfile.TemporaryDirectory(prefix=".coding-index-", dir=directory)
    connection: sqlite3.Connection | None = None
    try:
        db_path = Path(temp_dir.name) / "coding.sqlite"
        connection = sqlite3.connect(db_path)
        connection.execute(
            "CREATE TABLE coding (target TEXT, variant TEXT, uniprot TEXT, "
            "position INTEGER, wt TEXT, mut TEXT)"
        )

        filter_expression = None
        if genes is not None:
            filter_expression = ds.field("targetId").isin(list(genes))

        insert_sql = "INSERT INTO coding VALUES (?, ?, ?, ?, ?, ?)"
        for batch in _dataset_batches(
            coding,
            CODING_COLUMNS,
            batch_size,
            filter_expression=filter_expression,
        ):
            connection.executemany(
                insert_sql,
                (
                    (
                        row["targetId"],
                        row["variantId"],
                        accession,
                        row["aminoAcidPosition"],
                        row["referenceAminoAcid"],
                        row["alternateAminoAcid"],
                    )
                    for row in batch
                    for accession in explode(row["uniprotAccessions"])
                ),
            )

        connection.commit()
        connection.execute("CREATE INDEX coding_variant ON coding(variant)")
        connection.execute(
            "CREATE INDEX coding_target_variant ON coding(target, variant)"
        )
        connection.commit()
        yield connection
    finally:
        if connection is not None:
            connection.close()
        temp_dir.cleanup()


def tsv_value(value: Any) -> Any:
    """Serialize nulls and nested values using the MAVE reference TSV rules."""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if (
        value is None
        or value is pd.NA
        or isinstance(value, float)
        and math.isnan(value)
    ):
        return "NA"
    if isinstance(value, (dict, list)):
        return json.dumps(value, separators=(",", ":"), sort_keys=True)
    return value


def generate_clinvar(
    evidence: Path | Sequence[Path],
    connection: sqlite3.Connection,
    mapping: pd.DataFrame,
    destination: Path,
    *,
    batch_size: int = 65_536,
) -> dict[str, int]:
    """Join mapped ClinVar evidence to coding rows using an existing index.

    Evidence is left-joined on ``targetId`` and ``variantId``. Duplicate coding
    rows and accessions are retained, clinical significance lists are exploded,
    and evidence with no coding match remains in the output.
    """
    validate_batch_size(batch_size)
    if "To" not in mapping:
        raise ValueError("ID mapping requires a To column")
    genes = sorted(set(mapping["To"].dropna().astype(str)))
    destination.parent.mkdir(parents=True, exist_ok=True)
    coding_rows = evidence_rows = output_rows = 0

    if genes:
        targets = ", ".join("?" for _ in genes)
        coding_rows = connection.execute(
            f"SELECT COUNT(*) FROM coding WHERE target IN ({targets})", genes
        ).fetchone()[0]

    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(CLINVAR_COLUMNS)
        for batch in parquet_batches(
            evidence,
            CLINVAR_COLUMNS[:8],
            batch_size=batch_size,
            required_columns=CLINVAR_EVIDENCE_REQUIRED_COLUMNS,
            filter_column="targetId",
            filter_values=genes,
        ):
            for row in batch:
                evidence_rows += 1
                labels = explode(row.get("clinicalSignificances"))
                if not labels:
                    continue
                matches = connection.execute(
                    "SELECT uniprot, position, wt, mut FROM coding "
                    "WHERE target = ? AND variant = ? ORDER BY rowid",
                    (row.get("targetId"), row.get("variantId")),
                ).fetchall()
                if not matches:
                    matches = [(None, None, None, None)]
                for label in labels:
                    values = [
                        label if column == "clinicalSignificances" else row.get(column)
                        for column in CLINVAR_COLUMNS[:8]
                    ]
                    for match in matches:
                        writer.writerow(
                            [tsv_value(value) for value in (*values, *match)]
                        )
                        output_rows += 1

    return {
        "coding_rows": coding_rows,
        "evidence_rows": evidence_rows,
        "rows": output_rows,
    }


def generate_gwas(
    credible: Path | Sequence[Path],
    connection: sqlite3.Connection,
    destination: Path,
    *,
    batch_size: int = 65_536,
) -> dict[str, int]:
    """Stream credible-set rows and join coding variants by ``variantId``.

    Membership is based on all coding rows, as in the source R pipeline. The
    join preserves duplicate coding records and emits one row per UniProt
    accession. Both source data and the join index are processed in bounded
    batches without collecting either complete dataset in memory.
    """
    validate_batch_size(batch_size)
    destination.parent.mkdir(parents=True, exist_ok=True)
    credible_rows = matched_rows = output_rows = 0
    coding_rows = connection.execute("SELECT COUNT(*) FROM coding").fetchone()[0]
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(GWAS_COLUMNS)
        for batch in parquet_batches(
            credible, GWAS_INPUT_COLUMNS, batch_size=batch_size
        ):
            for row in batch:
                credible_rows += 1
                matches = connection.execute(
                    "SELECT target, uniprot, position, wt, mut "
                    "FROM coding WHERE variant IS ? ORDER BY rowid",
                    (row["variantId"],),
                ).fetchall()
                if not matches:
                    continue
                matched_rows += 1
                for match in matches:
                    writer.writerow(
                        [tsv_value(row[column]) for column in GWAS_INPUT_COLUMNS]
                        + [tsv_value(value) for value in match]
                    )
                    output_rows += 1
    return {
        "coding_rows": coding_rows,
        "credible_rows": credible_rows,
        "matched_rows": matched_rows,
        "rows": output_rows,
    }
