"""Extract variant score rows from the CSV files in a MaveDB dump."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Iterator

SIMPLE_P_PATTERN = re.compile(
    r"^(?:[^:]+:)?p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)$"
)
HGVS_SCHEME_PATTERN = re.compile(r"^(?:[^:]+:)?([a-z]\.)")
ACCESSION_VARIANT_PATTERN = re.compile(r"^(.+?)#(\d+)$")
MISSING_SCORE_VALUES = {"", "None", "NaN", "nan", "null", "NA"}
CSV_FIELDS = (
    "accession",
    "variant_num",
    "hgvs_nt",
    "hgvs_pro",
    "is_simple_p",
    "ref_aa",
    "position",
    "alt_aa",
    "score",
)

MAVEDB_DIRECTORY = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = MAVEDB_DIRECTORY / "tmp" / "mavedb_dump" / "csv"
DEFAULT_OUTPUT = MAVEDB_DIRECTORY / "tmp" / "mave_score.csv"


@dataclass
class ScoreSummary:
    """Counts collected while score rows are streamed to the output CSV."""

    files: int = 0
    rows: int = 0
    hgvs_schemes: Counter[str] = field(default_factory=Counter)
    protein_variant_complexity: Counter[str] = field(default_factory=Counter)

    def add_row(self, scheme: str, is_simple_p: bool) -> None:
        self.rows += 1
        self.hgvs_schemes[scheme] += 1
        if scheme == "p.":
            category = "simple" if is_simple_p else "complex"
            self.protein_variant_complexity[category] += 1


def is_simple_p_scheme(hgvs_pro: str) -> bool:
    """Return whether a protein HGVS value is a single-residue substitution."""
    return SIMPLE_P_PATTERN.fullmatch(hgvs_pro) is not None


def parse_simple_p(hgvs_pro: str) -> tuple[str | None, int | None, str | None]:
    """Return reference amino acid, position, and alternate amino acid."""
    match = SIMPLE_P_PATTERN.fullmatch(hgvs_pro)
    if match is None:
        return None, None, None
    ref_aa, position, alt_aa = match.groups()
    return ref_aa, int(position), alt_aa


def detect_hgvs_scheme(value: str | None) -> str:
    """Return the HGVS scheme prefix, or ``empty`` / ``unknown``."""
    if not value:
        return "empty"
    match = HGVS_SCHEME_PATTERN.match(value)
    return match.group(1) if match else "unknown"


def _score_files(input_path: Path) -> list[Path]:
    """Find score CSVs in a dump directory or a CSV directory."""
    if not input_path.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_path}")

    csv_directory = input_path / "csv"
    if csv_directory.is_dir():
        input_path = csv_directory

    files = sorted(input_path.glob("*.scores.csv"))
    if not files:
        raise FileNotFoundError(
            f"No *.scores.csv files found in input directory: {input_path}"
        )
    return files


def _iter_score_rows(
    score_files: list[Path], summary: ScoreSummary
) -> Iterator[dict[str, object]]:
    """Read score files one at a time and yield extracted variant rows."""
    for score_file in score_files:
        try:
            with score_file.open("r", newline="", encoding="utf-8-sig") as input_file:
                reader = csv.DictReader(input_file)
                for source_row in reader:
                    yield _extract_score_row(source_row, summary)
        except (OSError, UnicodeError, csv.Error) as error:
            raise ValueError(f"Could not process {score_file.name}: {error}") from error
        summary.files += 1


def _extract_score_row(
    source_row: dict[str, str | None], summary: ScoreSummary
) -> dict[str, object]:
    """Normalize one source row and parse simple protein substitutions."""
    full_accession = (source_row.get("accession") or "").strip()
    accession_match = ACCESSION_VARIANT_PATTERN.fullmatch(full_accession)
    if accession_match:
        accession = accession_match.group(1)
        variant_num: int | str = int(accession_match.group(2))
    else:
        accession = full_accession
        variant_num = "NA"

    hgvs_pro = (source_row.get("hgvs_pro") or "").strip()
    scheme = detect_hgvs_scheme(hgvs_pro)
    is_simple_p = scheme == "p." and is_simple_p_scheme(hgvs_pro)
    if is_simple_p:
        ref_aa, position, alt_aa = parse_simple_p(hgvs_pro)
    else:
        ref_aa, position, alt_aa = None, None, None

    summary.add_row(scheme, is_simple_p)
    raw_score = (source_row.get("score") or "").strip()
    score = "NA" if raw_score in MISSING_SCORE_VALUES else raw_score

    return {
        "accession": accession,
        "variant_num": variant_num,
        "hgvs_nt": source_row.get("hgvs_nt") or "",
        "hgvs_pro": hgvs_pro,
        "is_simple_p": is_simple_p,
        "ref_aa": ref_aa if ref_aa else "NA",
        "position": position if position is not None else "NA",
        "alt_aa": alt_aa if alt_aa else "NA",
        "score": score,
    }


def extract_scores(input_path: Path, output_path: Path) -> ScoreSummary:
    """Stream extracted variant scores from a MaveDB dump into a CSV."""
    score_files = _score_files(input_path)
    output_resolved = output_path.resolve()
    if output_resolved in {score_file.resolve() for score_file in score_files}:
        raise ValueError("The output CSV path must differ from every input score file.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary = ScoreSummary()
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in _iter_score_rows(score_files, summary):
            writer.writerow(row)

    return summary


def _print_summary(summary: ScoreSummary, output_path: Path) -> None:
    print(f"Finished streaming {summary.files} files into {output_path}\n")
    print(f"Total rows processed: {summary.rows}")
    print("HGVS scheme breakdown in hgvs_pro:")
    for scheme, count in sorted(summary.hgvs_schemes.items()):
        print(f"  {scheme}: {count}")

    if "p." in summary.hgvs_schemes:
        print("\n'p.' scheme breakdown:")
        for category, count in sorted(summary.protein_variant_complexity.items()):
            print(f"  {category}: {count}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract variant scores from MaveDB score CSV files."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "Path to the dump directory or its csv/ directory "
            f"(default: {DEFAULT_INPUT})."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT}).",
    )
    args = parser.parse_args()

    try:
        summary = extract_scores(args.input, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    _print_summary(summary, args.output)


if __name__ == "__main__":
    main()
