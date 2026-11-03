"""Streaming helpers and Open Targets GWAS variant export."""

from __future__ import annotations

import csv
import sqlite3
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Sequence

import pyarrow.dataset as ds

from data_sources.mavedb.analysis._homology_references import _tsv_value

CODING_COLUMNS = (
    "targetId",
    "variantId",
    "uniprotAccessions",
    "aminoAcidPosition",
    "referenceAminoAcid",
    "alternateAminoAcid",
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


def _validate_batch_size(batch_size: int) -> None:
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
    filter_expression: ds.Expression | None = None,
) -> Iterator[list[dict[str, Any]]]:
    _validate_batch_size(batch_size)
    if not columns:
        raise ValueError("At least one Parquet column is required")
    if len(set(columns)) != len(columns):
        raise ValueError("Parquet columns must be unique")

    dataset = ds.dataset(
        _source_paths(paths), format="parquet", exclude_invalid_files=True
    )
    missing = sorted(set(columns) - set(dataset.schema.names))
    if missing:
        raise ValueError(f"Parquet dataset is missing requested columns: {missing}")

    scanner = dataset.scanner(
        columns=list(columns),
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
) -> Iterator[list[dict[str, Any]]]:
    """Yield bounded batches of dictionaries, requiring every selected column."""
    yield from _dataset_batches(paths, columns, batch_size)


def _explode(value: Any) -> list[Any]:
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
    _validate_batch_size(batch_size)
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
            records = [
                (
                    row["targetId"],
                    row["variantId"],
                    accession,
                    row["aminoAcidPosition"],
                    row["referenceAminoAcid"],
                    row["alternateAminoAcid"],
                )
                for row in batch
                for accession in _explode(row["uniprotAccessions"])
            ]
            connection.executemany(insert_sql, records)
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
    return _tsv_value(value)


def generate_gwas(
    credible: Path | Sequence[Path],
    coding: Path | Sequence[Path],
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
    _validate_batch_size(batch_size)
    destination.parent.mkdir(parents=True, exist_ok=True)
    credible_rows = matched_rows = output_rows = 0
    with tempfile.TemporaryDirectory(
        prefix=".gwas-join-", dir=destination.parent
    ) as working:
        with coding_index(coding, Path(working), batch_size=batch_size) as connection:
            coding_rows = connection.execute("SELECT COUNT(*) FROM coding").fetchone()[
                0
            ]
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
                                [
                                    tsv_value(row[column])
                                    for column in GWAS_INPUT_COLUMNS
                                ]
                                + [tsv_value(value) for value in match]
                            )
                            output_rows += 1
    return {
        "coding_rows": coding_rows,
        "credible_rows": credible_rows,
        "matched_rows": matched_rows,
        "rows": output_rows,
    }
