"""Generate pharmacogenomics, drug, and ProtVar pocket annotations."""

from __future__ import annotations

import csv
import math
import re
import tempfile
from collections.abc import Sequence
from itertools import chain
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds

from data_sources.mavedb.analysis._homology_references import _explode
from data_sources.mavedb.analysis._variant_data import (
    coding_index,
    parquet_batches,
    tsv_value,
)

PHARMACOGENOMICS_COLUMNS = (
    "targetFromSourceId",
    "variantId",
    "drugs",
    "genotypeId",
    "genotype",
    "genotypeAnnotationText",
    "directionality",
    "evidenceLevel",
    "variantFunctionalConsequenceId",
    "studyId",
    "phenotypeFromSourceId",
    "phenotypeText",
    "pgxCategory",
    "literature",
)
PHARMACOGENOMICS_OUTPUT_COLUMNS = (
    "ensembl_id",
    "uniprot_id",
    "variant_id",
    "position",
    "ref",
    "alt",
    "genotypeId",
    "drug",
    "chembl_id",
    "direct_target",
    "genotype",
    "genotype_id",
    "genotype_desc",
    "directionality",
    "evidence_level",
    "functional_consequence",
    "study",
    "phenotype",
    "phenotype_desc",
    "pgx_category",
    "literature",
)
_DRUG_COLUMNS = (
    "linkedTargets",
    "id",
    "name",
    "canonicalSmiles",
    "inchiKey",
    "drugType",
    "blackBoxWarning",
    "yearOfFirstApproval",
    "maximumClinicalTrialPhase",
    "isApproved",
    "hasBeenWithdrawn",
    "synonyms",
    "description",
)
DRUG_OUTPUT_COLUMNS = (
    "uniprot_id",
    "ensembl_id",
    "target_count",
    "chembl_id",
    "name",
    "smiles",
    "inchi_key",
    "type",
    "black_box_warning",
    "year_approved",
    "maximum_trial_phase",
    "approved",
    "withdrawn",
    "synonyms",
    "description",
)
_POCKET_REQUIRED_COLUMNS = (
    "struct_id",
    "pocket_id",
    "pocket_resid",
    "pocket_score_combined_scaled",
    "pocket_rad_gyration",
    "pocket_plddt_mean",
)
_POCKET_PREFIX = (
    "uniprot_id",
    "pocket_id",
    "position",
    "high_confidence",
    "pocket_score",
)


def generate_pharmacogenomics(
    pharma: Path | Sequence[Path],
    coding: Path | Sequence[Path],
    mapping: pd.DataFrame,
    destination: Path,
    *,
    batch_size: int = 65_536,
) -> dict[str, int]:
    """Write the mapped Open Targets pharmacogenomics annotations as a TSV."""
    mapped_uniprots, _ = _mapping_values(mapping)
    _validate_batch_size(batch_size)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output_rows = 0
    pharma_columns = _pharmacogenomics_columns(pharma)

    with tempfile.TemporaryDirectory(
        prefix=".pharmacogenomics-join-", dir=destination.parent
    ) as directory:
        with coding_index(coding, Path(directory), batch_size=batch_size) as connection:
            batches = iter(
                parquet_batches(pharma, pharma_columns, batch_size=batch_size)
            )
            first_batch = next(batches, None)
            with destination.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(PHARMACOGENOMICS_OUTPUT_COLUMNS)
                for batch in chain(
                    () if first_batch is None else (first_batch,), batches
                ):
                    for row in batch:
                        for drug in _explode(row["drugs"]):
                            if drug is None:
                                drug = {}
                            if not isinstance(drug, dict):
                                raise ValueError(
                                    "Pharmacogenomics drugs entries must be structs"
                                )
                            matches = connection.execute(
                                "SELECT uniprot, position, wt, mut "
                                "FROM coding WHERE target IS ? AND variant IS ? "
                                "ORDER BY rowid",
                                (row["targetFromSourceId"], row["variantId"]),
                            )
                            for uniprot, position, reference, alternate in matches:
                                if uniprot not in mapped_uniprots:
                                    continue
                                output = (
                                    row["targetFromSourceId"],
                                    uniprot,
                                    row["variantId"],
                                    position,
                                    reference,
                                    alternate,
                                    row["genotypeId"],
                                    drug.get("drugFromSource"),
                                    drug.get("drugId"),
                                    _direct_target_value(row, drug),
                                    row["genotype"],
                                    row["genotypeId"],
                                    row["genotypeAnnotationText"],
                                    row["directionality"],
                                    row["evidenceLevel"],
                                    row["variantFunctionalConsequenceId"],
                                    row["studyId"],
                                    row["phenotypeFromSourceId"],
                                    row["phenotypeText"],
                                    row["pgxCategory"],
                                    _collapsed_list(row["literature"]),
                                )
                                writer.writerow([tsv_value(value) for value in output])
                                output_rows += 1

    return {"rows": output_rows}


def generate_drugs(
    drugs: Path | Sequence[Path],
    mapping: pd.DataFrame,
    destination: Path,
    *,
    batch_size: int = 65_536,
) -> dict[str, int]:
    """Write drug annotations linked to MAVE-mapped Ensembl targets."""
    _validate_batch_size(batch_size)
    pairs = _mapping_pairs(mapping)
    uniprots_by_gene: dict[str, list[str]] = {}
    for uniprot, gene in pairs:
        uniprots_by_gene.setdefault(gene, []).append(uniprot)

    destination.parent.mkdir(parents=True, exist_ok=True)
    output_rows = 0
    batches = iter(parquet_batches(drugs, _DRUG_COLUMNS, batch_size=batch_size))
    first_batch = next(batches, None)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(DRUG_OUTPUT_COLUMNS)
        for batch in chain(() if first_batch is None else (first_batch,), batches):
            for row in batch:
                linked_targets = row["linkedTargets"]
                if linked_targets is None:
                    continue
                if not isinstance(linked_targets, dict):
                    raise ValueError("Drug linkedTargets entries must be structs")
                target_count = linked_targets.get("count")
                for gene in _explode(linked_targets.get("rows")):
                    if gene is None:
                        continue
                    for uniprot in uniprots_by_gene.get(str(gene), ()):
                        output = (
                            uniprot,
                            gene,
                            target_count,
                            row["id"],
                            row["name"],
                            row["canonicalSmiles"],
                            row["inchiKey"],
                            row["drugType"],
                            row["blackBoxWarning"],
                            row["yearOfFirstApproval"],
                            row["maximumClinicalTrialPhase"],
                            row["isApproved"],
                            row["hasBeenWithdrawn"],
                            _collapsed_list(row["synonyms"]),
                            row["description"],
                        )
                        writer.writerow([tsv_value(value) for value in output])
                        output_rows += 1

    return {"rows": output_rows}


def generate_pockets(
    source: Path,
    mapping: pd.DataFrame,
    destination: Path,
    *,
    batch_size: int = 65_536,
) -> dict[str, int]:
    """Explode ProtVar residue sets and retain pockets for mapped UniProt IDs."""
    _validate_batch_size(batch_size)
    mapped_uniprots, _ = _mapping_values(mapping)
    header = list(pd.read_csv(source, sep="\t", nrows=0).columns)
    missing = set(_POCKET_REQUIRED_COLUMNS) - set(header)
    if missing:
        raise ValueError(f"Pocket TSV is missing required columns: {sorted(missing)}")
    metric_start = header.index("pocket_rad_gyration")
    metric_end = header.index("pocket_plddt_mean")
    if metric_start > metric_end:
        raise ValueError("Pocket metric columns are in an invalid order")
    metric_columns = tuple(header[metric_start : metric_end + 1])
    output_columns = (*_POCKET_PREFIX, *metric_columns)

    destination.parent.mkdir(parents=True, exist_ok=True)
    output_rows = 0
    chunks = pd.read_csv(
        source,
        sep="\t",
        chunksize=batch_size,
        dtype={"struct_id": "string", "pocket_resid": "string"},
    )
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_columns)
        for chunk in chunks:
            for row in chunk.to_dict(orient="records"):
                uniprot = _uniprot_from_struct_id(row["struct_id"])
                if uniprot not in mapped_uniprots:
                    continue
                for position in _pocket_residues(row["pocket_resid"]):
                    score = _numeric_or_missing(row["pocket_score_combined_scaled"])
                    high_confidence = None if score is None else score > 800
                    output = (
                        uniprot,
                        row["pocket_id"],
                        position,
                        high_confidence,
                        score,
                        *(row[column] for column in metric_columns),
                    )
                    writer.writerow([tsv_value(value) for value in output])
                    output_rows += 1

    return {"rows": output_rows}


def _validate_batch_size(batch_size: int) -> None:
    if (
        isinstance(batch_size, bool)
        or not isinstance(batch_size, int)
        or batch_size <= 0
    ):
        raise ValueError("Batch size must be a positive integer")


def _pharmacogenomics_columns(
    pharma: Path | Sequence[Path],
) -> tuple[str, ...]:
    if isinstance(pharma, Path):
        sources: str | list[str] = str(pharma)
    else:
        sources = [str(path) for path in pharma]
        if not sources:
            raise ValueError("At least one Parquet path is required")
    schema = ds.dataset(sources, format="parquet", exclude_invalid_files=True).schema
    missing = sorted(set(PHARMACOGENOMICS_COLUMNS) - set(schema.names))
    if missing:
        raise ValueError(f"Parquet dataset is missing requested columns: {missing}")

    drug_type = schema.field("drugs").type
    while (
        pa.types.is_list(drug_type)
        or pa.types.is_large_list(drug_type)
        or pa.types.is_fixed_size_list(drug_type)
    ):
        drug_type = drug_type.value_type
    nested_fields = (
        {field.name for field in drug_type} if pa.types.is_struct(drug_type) else set()
    )
    required_drug_fields = {"drugFromSource", "drugId"}
    if not required_drug_fields.issubset(nested_fields):
        missing_drug_fields = sorted(required_drug_fields - nested_fields)
        raise ValueError(
            "Pharmacogenomics drugs structs are missing required fields: "
            f"{missing_drug_fields}"
        )

    if "isDirectTarget" in schema.names:
        return (*PHARMACOGENOMICS_COLUMNS, "isDirectTarget")
    if "isDirectTarget" in nested_fields:
        return PHARMACOGENOMICS_COLUMNS
    raise ValueError(
        "Pharmacogenomics data requires isDirectTarget at the top level "
        "or inside drugs structs"
    )


def _mapping_values(mapping: pd.DataFrame) -> tuple[set[str], list[str]]:
    pairs = _mapping_pairs(mapping)
    return {uniprot for uniprot, _ in pairs}, list(
        dict.fromkeys(gene for _, gene in pairs)
    )


def _mapping_pairs(mapping: pd.DataFrame) -> list[tuple[str, str]]:
    if not {"From", "To"}.issubset(mapping.columns):
        raise ValueError("ID mapping requires From and To columns")
    pairs = []
    for uniprot, gene in mapping[["From", "To"]].itertuples(index=False, name=None):
        if _missing(uniprot) or _missing(gene):
            continue
        uniprot, gene = str(uniprot).strip(), str(gene).strip()
        if uniprot and gene:
            pairs.append((uniprot, re.sub(r"\.\d+$", "", gene)))
    if not pairs:
        raise ValueError("ID mapping contains no usable UniProt/Ensembl pairs")
    return pairs


def _collapsed_list(value: Any) -> Any:
    if _missing(value):
        return None
    if not isinstance(value, list):
        return value
    return ",".join("NA" if _missing(item) else str(item) for item in value)


def _direct_target_value(row: dict[str, Any], drug: dict[str, Any]) -> Any:
    value = row.get("isDirectTarget")
    return drug.get("isDirectTarget") if _missing(value) else value


def _uniprot_from_struct_id(value: Any) -> str | None:
    if _missing(value):
        return None
    match = re.fullmatch(r"([A-Z0-9]*)-.*", str(value))
    if match is None:
        raise ValueError(f"Invalid ProtVar struct_id: {value!r}")
    return match.group(1)


def _pocket_residues(value: Any) -> list[int | None]:
    if _missing(value):
        return [None]
    residues = re.sub(r"[{}]", "", str(value)).split(",")
    if len(residues) == 1 and not residues[0].strip():
        return [None]
    result = []
    for residue in residues:
        try:
            result.append(int(float(residue.strip())))
        except (OverflowError, ValueError):
            result.append(None)
    return result


def _numeric_or_missing(value: Any) -> int | float | None:
    if _missing(value):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return int(number) if number.is_integer() else number


def _missing(value: Any) -> bool:
    if value is None or value is pd.NA:
        return True
    if isinstance(value, (float, int)) and not isinstance(value, bool):
        return math.isnan(value) if isinstance(value, float) else False
    return False
