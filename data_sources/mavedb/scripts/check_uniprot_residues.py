#!/usr/bin/env python3
"""Check MaveDB reference amino acids against an exported ProtVar residue TSV."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from collections.abc import Collection, Iterable, Iterator
import csv
from dataclasses import dataclass
import json
from operator import itemgetter
from pathlib import Path

STATUSES = (
    "match",
    "mismatch",
    "ambiguous_target",
    "missing_uniprot_accession",
    "unsupported_reference_residue",
    "invalid_uniprot_position",
    "missing_reference_accession",
    "missing_reference_position",
    "conflicting_reference_residues",
)
OFFSET_SHIFTS = (-1, 1)


@dataclass(frozen=True)
class Target:
    accession: str = ""
    offset: int | None = None
    issue: str = ""

    def uniprot_position(self, position: int) -> int:
        return position + (self.offset or 0)


@dataclass
class ShiftComparison:
    comparable_checks: int = 0
    current_matches: int = 0
    shifted_matches: int = 0
    mismatches_resolved: int = 0
    matches_lost: int = 0

    def add(self, current: bool, shifted: bool) -> None:
        self.comparable_checks += 1
        self.current_matches += current
        self.shifted_matches += shifted
        self.mismatches_resolved += shifted and not current
        self.matches_lost += current and not shifted

    def qualifies(self, total_checks: int) -> bool:
        count = self.comparable_checks
        return (
            count >= 10
            and 100 * count >= 90 * total_checks
            and 100 * self.shifted_matches >= 95 * count
            and 100 * (self.shifted_matches - self.current_matches) >= 20 * count
        )


def read_tsv(path: Path, *fields: str) -> Iterator[tuple[str, ...]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        missing = set(fields) - set(header)
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        select = itemgetter(*(header.index(field) for field in fields))
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"{path}:{reader.line_num} has the wrong column count")
            yield select(row)


def target_values(value: str) -> list[str | int | None]:
    return json.loads(value) if value.startswith("[") else [value]


def read_targets(path: Path) -> dict[str, Target | None]:
    targets: dict[str, Target | None] = {}
    for dataset, taxid, category, accession, offset in read_tsv(
        path,
        "dataset_id",
        "taxid",
        "target_category",
        "identifier_uniprot",
        "offset_uniprot",
    ):
        if not dataset or dataset in targets:
            raise ValueError(f"Missing or duplicate dataset_id in {path}: {dataset!r}")
        taxa = target_values(taxid)
        categories = target_values(category)
        eligible = any(
            str(taxon) == "9606" and kind == "protein_coding"
            for taxon, kind in zip(taxa, categories, strict=True)
        )
        if not eligible:
            targets[dataset] = None
        elif len(taxa) > 1 or accession.startswith("[") or offset.startswith("["):
            targets[dataset] = Target(issue="ambiguous_target")
        elif accession in {"", "NA"}:
            targets[dataset] = Target(issue="missing_uniprot_accession")
        else:
            targets[dataset] = Target(
                accession, None if offset in {"", "NA"} else int(offset)
            )
    return targets


def select_offset_shift(
    observations: Iterable[tuple[int, str, str]],
    target: Target,
    positions: dict[int, set[str]],
) -> tuple[int, ShiftComparison | None, int]:
    """Compare offsets on shared reference positions and accept a single winner."""
    comparisons = {shift: ShiftComparison() for shift in OFFSET_SHIFTS}
    total = 0
    for position, _, residue in observations:
        mapped = target.uniprot_position(position)
        current = positions.get(mapped, ())
        if len(residue) != 1 or len(current) != 1:
            continue
        total += 1
        for shift, comparison in comparisons.items():
            shifted = positions.get(mapped + shift, ())
            if len(shifted) == 1:
                comparison.add(residue in current, residue in shifted)
    passing = [
        (shift, comparison)
        for shift, comparison in comparisons.items()
        if comparison.qualifies(total)
    ]
    if len(passing) == 1:
        shift, comparison = passing[0]
        return shift, comparison, total
    return 0, None, total


def residue_status(
    target: Target,
    residue: str,
    position: int,
    reference: Collection[str],
    has_reference_accession: bool,
) -> str:
    if target.issue:
        return target.issue
    if len(residue) != 1:
        return "unsupported_reference_residue"
    if position < 1:
        return "invalid_uniprot_position"
    if not has_reference_accession:
        return "missing_reference_accession"
    if not reference:
        return "missing_reference_position"
    if len(reference) > 1:
        return "conflicting_reference_residues"
    return "match" if residue in reference else "mismatch"


def check_residues(scores: Path, metadata: Path, residues: Path, output: Path) -> None:
    targets = read_targets(metadata)
    observations: defaultdict[str, Counter[tuple[int, str, str]]] = defaultdict(Counter)
    requested: dict[str, dict[int, set[str]]] = {}
    for dataset, simple, position, abbreviation, residue in read_tsv(
        scores, "dataset_id", "is_simple_p", "position", "ref_aa_abbr", "ref_aa_sym"
    ):
        if dataset not in targets:
            raise ValueError(f"Score dataset is absent from metadata: {dataset}")
        target = targets[dataset]
        if target is None or simple.upper() != "TRUE":
            continue
        source_position = int(position)
        observations[dataset][source_position, abbreviation, residue] += 1

    for dataset, dataset_observations in observations.items():
        target = targets[dataset]
        assert target is not None
        if not target.accession:
            continue
        for source_position, _, residue in dataset_observations:
            mapped_position = target.uniprot_position(source_position)
            if len(residue) == 1 and mapped_position >= 0:
                positions = requested.setdefault(target.accession, {})
                for position in (
                    mapped_position - 1,
                    mapped_position,
                    mapped_position + 1,
                ):
                    if position > 0:
                        positions.setdefault(position, set())

    reference_accessions: set[str] = set()
    for accession, position, residue in read_tsv(
        residues, "accession", "protein_position", "protein_seq"
    ):
        if accession in requested:
            reference_accessions.add(accession)
            if position and residue:
                values = requested[accession].get(int(position))
                if values is not None:
                    values.add(residue)

    checks: Counter[tuple[str, bool]] = Counter()
    counts: Counter[tuple[str, bool]] = Counter()
    output.mkdir(parents=True, exist_ok=True)
    with (
        (output / "residue_checks.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as details,
        (output / "dataset_summary.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as datasets,
        (output / "applied_offset_corrections.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as corrections,
    ):
        details_writer = csv.writer(details, delimiter="\t", lineterminator="\n")
        datasets_writer = csv.writer(datasets, delimiter="\t", lineterminator="\n")
        corrections_writer = csv.writer(
            corrections, delimiter="\t", lineterminator="\n"
        )
        details_writer.writerow(
            (
                "dataset_id",
                "uniprot_accession",
                "source_position",
                "uniprot_offset",
                "offset_defaulted",
                "uniprot_position",
                "hgvs_ref_aa_abbr",
                "hgvs_ref_aa",
                "reference_residues",
                "status",
                "observations",
                "effective_uniprot_offset",
                "effective_uniprot_position",
                "effective_reference_residues",
                "effective_status",
                "offset_applied",
            )
        )
        datasets_writer.writerow(
            (
                "dataset_id",
                "uniprot_accession",
                "uniprot_offset",
                "offset_defaulted",
                "effective_uniprot_offset",
                "offset_shift",
                "offset_applied",
                "eligible_checks",
                "baseline_matches",
                "baseline_mismatches",
                "compared_checks",
                "matches",
                "mismatches",
                "uncheckable_checks",
                "agreement_fraction",
                "check_coverage_fraction",
            )
        )
        corrections_writer.writerow(
            (
                "dataset_id",
                "uniprot_accession",
                "current_offset",
                "offset_defaulted",
                "shift",
                "suggested_offset",
                "baseline_checks",
                "comparable_checks",
                "current_matches",
                "shifted_matches",
                "mismatches_resolved",
                "matches_lost",
            )
        )
        for dataset, dataset_observations in sorted(observations.items()):
            target = targets[dataset]
            assert target is not None
            positions = requested.get(target.accession, {})
            has_reference_accession = target.accession in reference_accessions
            shift, comparison, baseline_checks = select_offset_shift(
                dataset_observations, target, positions
            )
            effective_offset = (
                (target.offset or 0) + shift if target.accession else None
            )
            defaulted = str(target.offset is None).upper() if target.accession else ""
            applied = bool(shift)
            applied_text = str(applied).upper()
            baseline: Counter[str] = Counter()
            effective: Counter[str] = Counter()
            for (position, abbreviation, residue), count in sorted(
                dataset_observations.items()
            ):
                mapped = target.uniprot_position(position)
                reference = positions.get(mapped, ())
                status = residue_status(
                    target, residue, mapped, reference, has_reference_accession
                )
                effective_position = mapped + shift
                effective_reference = positions.get(effective_position, ())
                effective_status = (
                    residue_status(
                        target,
                        residue,
                        effective_position,
                        effective_reference,
                        has_reference_accession,
                    )
                    if applied
                    else status
                )
                details_writer.writerow(
                    (
                        dataset,
                        target.accession,
                        position,
                        target.offset,
                        defaulted,
                        mapped if target.accession else None,
                        abbreviation,
                        "" if residue == "NA" else residue,
                        ",".join(sorted(reference)),
                        status,
                        count,
                        effective_offset,
                        effective_position if target.accession else None,
                        ",".join(sorted(effective_reference)),
                        effective_status,
                        applied_text,
                    )
                )
                baseline[status] += 1
                effective[effective_status] += 1
                checks[effective_status, applied] += 1
                counts[effective_status, applied] += count
            eligible = len(dataset_observations)
            compared = effective["match"] + effective["mismatch"]
            datasets_writer.writerow(
                (
                    dataset,
                    target.accession,
                    target.offset,
                    defaulted,
                    effective_offset,
                    shift,
                    applied_text,
                    eligible,
                    baseline["match"],
                    baseline["mismatch"],
                    compared,
                    effective["match"],
                    effective["mismatch"],
                    eligible - compared,
                    effective["match"] / compared if compared else None,
                    compared / eligible,
                )
            )
            if comparison is not None:
                corrections_writer.writerow(
                    (
                        dataset,
                        target.accession,
                        target.offset,
                        defaulted,
                        shift,
                        effective_offset,
                        baseline_checks,
                        comparison.comparable_checks,
                        comparison.current_matches,
                        comparison.shifted_matches,
                        comparison.mismatches_resolved,
                        comparison.matches_lost,
                    )
                )

    with (output / "summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("status", "checks", "observations", "offset_applied"))
        writer.writerows(
            (
                status,
                checks[status, applied],
                counts[status, applied],
                str(applied).upper(),
            )
            for status in STATUSES
            for applied in (False, True)
            if not applied or checks[status, applied]
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("scores", "metadata", "residues", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    try:
        check_residues(args.scores, args.metadata, args.residues, args.output)
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
