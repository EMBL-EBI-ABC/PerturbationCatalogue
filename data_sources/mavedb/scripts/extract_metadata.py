"""Extract target gene identifiers from a MaveDB Zenodo dump."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

SOURCE_METADATA_FIELDS = (
    "dataset_id",
    "score_urn",
    "experiment_urn",
    "experiment_set_urn",
    "target_tax_id",
    "target_source_target_id",
    "target_sequence_id",
    "target_sequence",
    "study_references",
    "target_name",
    "target_category",
    "target_sequence_type",
    "target_sequence_label",
    "target_accession_assembly",
    "score_set_title",
    "score_set_short_description",
    "score_set_method_text",
    "score_set_num_variants",
    "experiment_title",
    "experiment_short_description",
    "experiment_method_text",
    "experiment_abstract_text",
    "experiment_creation_date",
    "experiment_doi_identifiers",
    "experiment_external_links",
    "experiment_extra_metadata",
    "experiment_modification_date",
    "experiment_primary_publication_identifiers",
    "experiment_published_date",
    "experiment_raw_read_identifiers",
    "experiment_secondary_publication_identifiers",
    "experiment_set_creation_date",
    "experiment_set_modification_date",
    "experiment_set_published_date",
    "score_set_abstract_text",
    "score_set_creation_date",
    "score_set_data_usage_policy",
    "score_set_dataset_columns",
    "score_set_doi_identifiers",
    "score_set_external_links",
    "score_set_extra_metadata",
    "score_set_license",
    "score_set_meta_analyzed_by_score_set_urns",
    "score_set_meta_analyzes_score_set_urns",
    "score_set_modification_date",
    "score_set_primary_publication_identifiers",
    "score_set_published_date",
    "score_set_secondary_publication_identifiers",
    "score_set_superseded_score_set",
    "score_set_superseding_score_set",
    "target_accession",
    "target_mapped_hgnc_name",
    "target_taxonomy_article_reference",
    "target_taxonomy_code",
    "target_taxonomy_common_name",
    "target_taxonomy_genome_id",
    "target_taxonomy_genome_identifier_id",
    "target_taxonomy_has_described_species_name",
    "target_taxonomy_id",
    "target_taxonomy_organism_name",
    "target_taxonomy_rank",
    "target_taxonomy_record_type",
    "target_taxonomy_tax_id",
    "target_taxonomy_url",
    "target_uniprot_id_from_mapped_metadata",
)
LEGACY_FIELDS = (
    "urn",
    "numVariants",
    "Gene",
    "GeneCategory",
    "GeneTaxId",
    "UniProt",
    "UniProtOffset",
    "RefSeq",
    "RefSeqOffset",
    "Ensembl",
    "EnsemblOffset",
)
CSV_FIELDS = (*SOURCE_METADATA_FIELDS, *LEGACY_FIELDS)
IDENTIFIER_DATABASES = ("UniProt", "RefSeq", "Ensembl")


@dataclass
class ExtractionSummary:
    """Counts collected while target gene rows are written."""

    rows: int = 0
    human_rows: int = 0
    nonhuman_rows: int = 0
    protein_coding_rows: int = 0
    rows_with_refseq: int = 0
    rows_with_gene_name: int = 0
    rows_with_uniprot: int = 0
    ensembl_prefix_counts: Counter[str] = field(default_factory=Counter)

    def add_row(self, row: dict[str, object]) -> None:
        """Update row-level summary counts."""
        self.rows += 1
        if row["GeneTaxId"] == "9606":
            self.human_rows += 1
        else:
            self.nonhuman_rows += 1

        if row["GeneCategory"] == "protein_coding":
            self.protein_coding_rows += 1
        if row["RefSeq"]:
            self.rows_with_refseq += 1
        if row["Gene"]:
            self.rows_with_gene_name += 1
        if row["UniProt"]:
            self.rows_with_uniprot += 1


def _as_mapping(value: object) -> dict[str, object]:
    """Return a JSON object value as a mapping, or an empty mapping."""
    return value if isinstance(value, dict) else {}


def _mappings(values: object) -> Iterator[dict[str, object]]:
    """Yield object members from a JSON array, skipping malformed members."""
    if isinstance(values, list):
        for value in values:
            if isinstance(value, dict):
                yield value


def _parse_offset(value: object) -> int | None:
    """Convert a non-empty offset to an integer, ignoring invalid values."""
    if value is None or value == "":
        return None
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return None


def _csv_value(value: object) -> object:
    """Encode structured values compactly so CSV cells retain their structure."""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return value


def _study_references(*records: dict[str, object]) -> list[str]:
    """Collect distinct DOI, identifier, or URL references from source records."""
    references: list[str] = []
    seen: set[str] = set()
    source_fields = (
        "primaryPublicationIdentifiers",
        "secondaryPublicationIdentifiers",
        "doiIdentifiers",
    )

    for record in records:
        for field_name in source_fields:
            values = record.get(field_name)
            if not isinstance(values, list):
                continue
            for value in values:
                if isinstance(value, str):
                    reference = value
                elif isinstance(value, dict):
                    reference = next(
                        (
                            value.get(key)
                            for key in ("doi", "identifier", "url")
                            if isinstance(value.get(key), str) and value.get(key)
                        ),
                        "",
                    )
                else:
                    reference = ""
                if reference and reference not in seen:
                    references.append(reference)
                    seen.add(reference)

    return references


def _experiment_sets(data: object) -> list[object]:
    """Get experiment sets from a full export object or a sliced JSON list."""
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        experiment_sets = data.get("experimentSets", [])
        if isinstance(experiment_sets, list):
            return experiment_sets
    raise ValueError(
        "Input JSON must be an object with an experimentSets list or a list."
    )


def _iter_rows(
    experiment_sets: list[object], summary: ExtractionSummary
) -> Iterator[dict[str, object]]:
    """Yield one CSV row per target gene in every score set."""
    for experiment_set in _mappings(experiment_sets):
        for experiment in _mappings(experiment_set.get("experiments")):
            for score_set in _mappings(experiment.get("scoreSets")):
                for target_gene in _mappings(score_set.get("targetGenes")):
                    target_sequence = _as_mapping(target_gene.get("targetSequence"))
                    target_accession = _as_mapping(target_gene.get("targetAccession"))
                    taxonomy = _as_mapping(target_sequence.get("taxonomy"))
                    score_urn = score_set.get("urn")
                    target_tax_id = taxonomy.get("code", taxonomy.get("taxId"))

                    source_metadata: dict[str, object] = {
                        "dataset_id": score_urn,
                        "score_urn": score_urn,
                        "experiment_urn": experiment.get("urn"),
                        "experiment_set_urn": experiment_set.get(
                            "urn", experiment.get("experimentSetUrn")
                        ),
                        "target_tax_id": target_tax_id,
                        "target_source_target_id": (
                            str(target_gene["id"])
                            if target_gene.get("id") is not None
                            else None
                        ),
                        "target_sequence_id": (
                            str(target_sequence["id"])
                            if target_sequence.get("id") is not None
                            else None
                        ),
                        "target_sequence": target_sequence.get("sequence"),
                        "study_references": _study_references(experiment, score_set),
                        "target_name": target_gene.get("name"),
                        "target_category": target_gene.get("category"),
                        "target_sequence_type": target_sequence.get("sequenceType"),
                        "target_sequence_label": target_sequence.get("label"),
                        "target_accession_assembly": target_accession.get("assembly"),
                        "score_set_title": score_set.get("title"),
                        "score_set_short_description": score_set.get(
                            "shortDescription"
                        ),
                        "score_set_method_text": score_set.get("methodText"),
                        "score_set_num_variants": score_set.get("numVariants"),
                        "experiment_title": experiment.get("title"),
                        "experiment_short_description": experiment.get(
                            "shortDescription"
                        ),
                        "experiment_method_text": experiment.get("methodText"),
                        "experiment_abstract_text": experiment.get("abstractText"),
                        "experiment_creation_date": experiment.get("creationDate"),
                        "experiment_doi_identifiers": experiment.get("doiIdentifiers"),
                        "experiment_external_links": experiment.get("externalLinks"),
                        "experiment_extra_metadata": experiment.get("extraMetadata"),
                        "experiment_modification_date": experiment.get(
                            "modificationDate"
                        ),
                        "experiment_primary_publication_identifiers": experiment.get(
                            "primaryPublicationIdentifiers"
                        ),
                        "experiment_published_date": experiment.get("publishedDate"),
                        "experiment_raw_read_identifiers": experiment.get(
                            "rawReadIdentifiers"
                        ),
                        "experiment_secondary_publication_identifiers": experiment.get(
                            "secondaryPublicationIdentifiers"
                        ),
                        "experiment_set_creation_date": experiment_set.get(
                            "creationDate"
                        ),
                        "experiment_set_modification_date": experiment_set.get(
                            "modificationDate"
                        ),
                        "experiment_set_published_date": experiment_set.get(
                            "publishedDate"
                        ),
                        "score_set_abstract_text": score_set.get("abstractText"),
                        "score_set_creation_date": score_set.get("creationDate"),
                        "score_set_data_usage_policy": score_set.get("dataUsagePolicy"),
                        "score_set_dataset_columns": score_set.get("datasetColumns"),
                        "score_set_doi_identifiers": score_set.get("doiIdentifiers"),
                        "score_set_external_links": score_set.get("externalLinks"),
                        "score_set_extra_metadata": score_set.get("extraMetadata"),
                        "score_set_license": score_set.get("license"),
                        "score_set_meta_analyzed_by_score_set_urns": score_set.get(
                            "metaAnalyzedByScoreSetUrns"
                        ),
                        "score_set_meta_analyzes_score_set_urns": score_set.get(
                            "metaAnalyzesScoreSetUrns"
                        ),
                        "score_set_modification_date": score_set.get(
                            "modificationDate"
                        ),
                        "score_set_primary_publication_identifiers": score_set.get(
                            "primaryPublicationIdentifiers"
                        ),
                        "score_set_published_date": score_set.get("publishedDate"),
                        "score_set_secondary_publication_identifiers": score_set.get(
                            "secondaryPublicationIdentifiers"
                        ),
                        "score_set_superseded_score_set": score_set.get(
                            "supersededScoreSet"
                        ),
                        "score_set_superseding_score_set": score_set.get(
                            "supersedingScoreSet"
                        ),
                        "target_accession": target_accession.get("accession"),
                        "target_mapped_hgnc_name": target_gene.get("mappedHgncName"),
                        "target_taxonomy_article_reference": taxonomy.get(
                            "articleReference"
                        ),
                        "target_taxonomy_code": taxonomy.get("code"),
                        "target_taxonomy_common_name": taxonomy.get("commonName"),
                        "target_taxonomy_genome_id": taxonomy.get("genomeId"),
                        "target_taxonomy_genome_identifier_id": taxonomy.get(
                            "genomeIdentifierId"
                        ),
                        "target_taxonomy_has_described_species_name": taxonomy.get(
                            "hasDescribedSpeciesName"
                        ),
                        "target_taxonomy_id": taxonomy.get("id"),
                        "target_taxonomy_organism_name": taxonomy.get("organismName"),
                        "target_taxonomy_rank": taxonomy.get("rank"),
                        "target_taxonomy_record_type": taxonomy.get("recordType"),
                        "target_taxonomy_tax_id": taxonomy.get("taxId", target_tax_id),
                        "target_taxonomy_url": taxonomy.get("url"),
                        "target_uniprot_id_from_mapped_metadata": target_gene.get(
                            "uniprotIdFromMappedMetadata"
                        ),
                    }

                    row = {
                        field_name: _csv_value(source_metadata[field_name])
                        for field_name in SOURCE_METADATA_FIELDS
                    }
                    row.update(
                        {
                            "urn": score_urn,
                            "numVariants": score_set.get("numVariants"),
                            "Gene": target_gene.get("name", ""),
                            "GeneCategory": target_gene.get("category", ""),
                            "GeneTaxId": (
                                str(target_tax_id) if target_tax_id is not None else ""
                            ),
                            "UniProt": "",
                            "UniProtOffset": None,
                            "RefSeq": "",
                            "RefSeqOffset": None,
                            "Ensembl": "",
                            "EnsemblOffset": None,
                        }
                    )
                    for external_id in _mappings(
                        target_gene.get("externalIdentifiers")
                    ):
                        identifier = _as_mapping(external_id.get("identifier"))
                        database_name = identifier.get("dbName")
                        identifier_value = identifier.get("identifier")

                        if database_name == "Ensembl" and identifier_value:
                            summary.ensembl_prefix_counts[
                                str(identifier_value)[:4]
                            ] += 1

                        if (
                            database_name in IDENTIFIER_DATABASES
                            and not row[str(database_name)]
                        ):
                            row[str(database_name)] = identifier_value
                            row[f"{database_name}Offset"] = _parse_offset(
                                external_id.get("offset")
                            )

                    if not row["UniProt"]:
                        row["UniProt"] = (
                            source_metadata["target_uniprot_id_from_mapped_metadata"]
                            or ""
                        )

                    summary.add_row(row)
                    yield row


def _resolve_input_file(input_path: Path) -> Path:
    """Accept either main.json itself or the directory containing it."""
    json_path = input_path / "main.json" if input_path.is_dir() else input_path
    if not json_path.is_file():
        raise FileNotFoundError(
            f"Could not find main.json or input JSON file at: {json_path}"
        )
    return json_path


def extract_score_set_data(
    input_path: Path, output_csv: Path, dry_run: bool = False
) -> ExtractionSummary:
    """Extract target gene identifiers to CSV and return summary counts."""
    json_path = _resolve_input_file(input_path)
    if json_path.resolve() == output_csv.resolve():
        raise ValueError("The output CSV path must differ from the input JSON path.")

    with json_path.open("r", encoding="utf-8") as input_file:
        data = json.load(input_file)
    experiment_sets = _experiment_sets(data)
    summary = ExtractionSummary()

    if dry_run:
        for _ in _iter_rows(experiment_sets, summary):
            pass
        print(
            f"Dry run: skipped writing {summary.rows} target gene rows to {output_csv}"
        )
        return summary

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in _iter_rows(experiment_sets, summary):
            writer.writerow(row)

    print(f"Extracted {summary.rows} target gene rows to {output_csv}")
    return summary


def _print_summary(summary: ExtractionSummary) -> None:
    print("\nSummary:")
    print(f"  Human target genes (taxId=9606): {summary.human_rows}")
    print(f"  Non-human target genes: {summary.nonhuman_rows}")
    print(f"  Protein-coding target genes: {summary.protein_coding_rows}")
    print(f"  Rows with RefSeq: {summary.rows_with_refseq}")
    print("  Ensembl identifier prefixes:")
    for prefix, count in sorted(summary.ensembl_prefix_counts.items()):
        print(f"     {prefix}: {count}")
    print(f"  Rows with gene name: {summary.rows_with_gene_name}")
    print(f"  Rows with UniProt: {summary.rows_with_uniprot}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract target gene identifiers from an unpacked MaveDB Zenodo dump."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data_sources/mavedb/tmp"),
        help=(
            "Path to main.json or the unpacked dump directory "
            "(default: data_sources/mavedb/tmp)."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data_sources/mavedb/tmp/mave_identifier.csv"),
        help=(
            "Output CSV path " "(default: data_sources/mavedb/tmp/mave_identifier.csv)."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse and print summary counts without writing the CSV.",
    )
    args = parser.parse_args()

    try:
        summary = extract_score_set_data(args.input, args.output, dry_run=args.dry_run)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    _print_summary(summary)


if __name__ == "__main__":
    main()
