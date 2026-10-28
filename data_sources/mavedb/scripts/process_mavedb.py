#!/usr/bin/env python3
"""Export complete MaveDB scores and selected dataset metadata.

Read a Zenodo dump (main.json and csv/*.scores.csv), processing one score set at
a time. Scores are not normalized. Standard deviations use the sample formula
(n - 1), and missing output values are written as NA. If the local taxonomy
cache is absent, the NCBI category dump is downloaded to data_sources/mavedb/tmp.
Write scores.tsv and metadata.tsv, joined by dataset_id. Optional compatibility
tables retain the original R singles/doubles/wide layout. No FASTA files are made.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
from functools import lru_cache
import json
import logging
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
from typing import TextIO
from urllib.request import Request, urlopen

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
import numpy as np
import pandas as pd

LOGGER = logging.getLogger(__name__)
MAVEDB_DIRECTORY = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = MAVEDB_DIRECTORY / "tmp" / "mavedb_dump"
DEFAULT_OUTPUT = MAVEDB_DIRECTORY / "output" / "processed"
DEFAULT_TAXONOMY = MAVEDB_DIRECTORY / "tmp" / "taxcat.dmp"
TAXONOMY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz"
EXCLUDED_URNS = {"urn:mavedb:00001237-a-1", "urn:mavedb:00001237-a-2"}
SOURCE_METADATA_FIELDS = {
    "title": "title",
    "short_description": "description",
    "abstract_text": "abstract",
    "num_variants": "n_variants_reported",
}
TARGET_SOURCE_FIELDS = {
    "target_name": ("name",),
    "target_category": ("category",),
    "sequence_type": ("target_sequence", "sequence_type"),
    "sequence": ("target_sequence", "sequence"),
    "label": ("target_sequence", "label"),
    "taxid": ("target_sequence", "taxonomy", "taxid"),
    "organism": ("target_sequence", "taxonomy", "organism_name"),
    "uniprot_id_from_mapped_metadata": ("uniprot_id_from_mapped_metadata",),
}
TARGET_DERIVED_FIELDS = {
    "kingdom": ("target_sequence", "taxonomy", "kingdom"),
    "reported_taxid": ("target_sequence", "taxonomy", "reported_taxid"),
    "derived_protein_sequence": ("derived_protein_sequence",),
}
TARGET_DATABASES = ("uniprot", "ensembl", "refseq")
SOURCE_COLUMNS = (
    (
        "dataset_id experiment_id experiment_set_id title description abstract "
        "n_variants_reported"
    ).split()
    + list(TARGET_SOURCE_FIELDS)
    + [
        f"{field}_{database}"
        for database in TARGET_DATABASES
        for field in ("identifier", "offset")
    ]
)
TSV_TEXT_ESCAPES = str.maketrans({"\\": "\\\\", "\t": "\\t", "\r": "\\r", "\n": "\\n"})
# Keep the combined R token for compatibility and match its components too.
WT_HGVS_CODES = {"p.=", "p.(=)", "_wt", "_sy", "n.=", "n.(=)", "n.=, n.(=)"}
AA_3_TO_1 = dict(
    zip(
        "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Asx Glx Xaa Sec Pyl Ter del".split(),
        list("ARNDCQEGHILKMFPSTWYVBZXUO") + ["stop", "del"],
    )
)
AA_3_TO_1["*"] = "stop"
AMINO_ACIDS = sorted("ACDEFGHIKLMNPQRSTVWY")
# Match the corrected R output order within each score, n and sd block.
MUTATIONS = [*AMINO_ACIDS, "del", "stop"]
SINGLE_PATTERN = re.compile(
    r"^(?:[^:]+:)?p\.([A-Za-z]+)([0-9]+)(?:delins)?([A-Za-z]+|=|\*)$"
)
DOUBLE_PATTERN = re.compile(
    r"^p\.\[?([A-Za-z]*)([0-9]*)([A-Za-z]*);?" r"([A-Za-z]*)([0-9]*)([A-Za-z]*)"
)
SCORE_COLUMNS = (
    "dataset_id source_row variant_num source_accession hgvs_nt hgvs_pro hgvs_splice "
    "is_simple_p ref_aa_abbr ref_aa_sym position alt_aa_abbr alt_aa_sym score score_raw"
).split()
SIMPLE_P_PATTERN = re.compile(
    r"^(?:[^:]+:)?p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)$"
)
ANALYSIS_TYPES = ("maves", "mpras", "splicing")
IDENTITY_COLUMNS = (
    "score_urn experiment_set_urn experiment_urn target_name identifier_uniprot taxid organism"
).split()
SINGLE_COLUMNS = IDENTITY_COLUMNS + "position wt mut score n sd".split()
DOUBLE_COLUMNS = (
    IDENTITY_COLUMNS
    + ("hgvs_pro position1 wt1 mut1 position2 wt2 mut2 score n sd").split()
)
WIDE_COLUMNS = (
    IDENTITY_COLUMNS
    + ["position", "wt"]
    + [f"{mut}_{metric}" for metric in ("score", "n", "sd") for mut in MUTATIONS]
    + "missing synonymous_score median_score mean_score sd_score stop synonymous del".split()
)
SCORE_METRICS = "mean_score median_score sd_score min_score max_score".split()
DERIVED_COLUMNS = (
    (
        "analysis_types primary_analysis source_accession_prefix n_variants "
        "n_scored_variants n_missing_scores n_invalid_scores n_analyzed_variants"
    ).split()
    + SCORE_METRICS
    + (
        "meta_study position_reference single_mutations multiple_mutations wt_mutants "
        "n_subs n_positions n_muts mean_variant_reps max_muts min_muts mean_muts "
        "mean_backgrounds region_start region_end region_length n_outlier p_outlier "
        "mean_position_median_score completeness completeness_missense positions_50 "
        "positions_75 positions_100 coverage_50 coverage_75 coverage_100 coverage_nonsense "
        "coverage_synonymous stop_proportion del_proportion synonymous_proportion delins_mutations"
    ).split()
)
BOOLEAN_COLUMNS = {
    "is_simple_p",
    "wt_mut",
    "meta_study",
    "single_mutations",
    "multiple_mutations",
    "wt_mutants",
    "delins_mutations",
    "stop",
    "synonymous",
    "del",
}
INTEGER_COLUMNS = {
    "source_row",
    "variant_num",
    "n_variants",
    "n_variants_reported",
    "n_scored_variants",
    "n_missing_scores",
    "n_invalid_scores",
    "n_analyzed_variants",
    "position",
    "position1",
    "position2",
    "n",
    "n_subs",
    "n_positions",
    "n_muts",
    "max_muts",
    "min_muts",
    "region_start",
    "region_end",
    "region_length",
    "n_outlier",
    "positions_50",
    "positions_75",
    "positions_100",
    "missing",
    "offset_uniprot",
    "offset_ensembl",
    "offset_refseq",
} | {f"{mut}_n" for mut in MUTATIONS}


@lru_cache(maxsize=None)
def snake_case(name: str) -> str:
    """Use one spelling for source field names."""
    return re.sub(r"(?<=[a-z0-9])(?=[A-Z])", "_", name).lower()


def present(value: object) -> bool:
    """Treat zero and FALSE as information, but omit empty source fields."""
    return value is not None and value != "" and value != [] and value != {}


def value_key(value: object) -> str:
    """Compare structured values independently of key order and text whitespace."""

    def canonical(item: object) -> object:
        if isinstance(item, str):
            return " ".join(item.split())
        if isinstance(item, dict):
            return {key: canonical(child) for key, child in item.items()}
        if isinstance(item, list):
            return [canonical(child) for child in item]
        return item

    return json.dumps(canonical(value), ensure_ascii=False, sort_keys=True)


def normalized(value: object) -> object:
    """Retain nonempty source data, normalizing names and duplicate list entries."""
    if isinstance(value, dict):
        return {
            snake_case(key): normalized(item)
            for key, item in value.items()
            if present(item)
        }
    if isinstance(value, list):
        result = []
        seen = set()
        for item in value:
            item = normalized(item)
            key = value_key(item)
            if present(item) and key not in seen:
                result.append(item)
                seen.add(key)
        return result
    return value


def normalize_target(target: dict[str, object]) -> dict[str, object]:
    """Select the reference extractors' target fields and retain derived provenance."""
    result = normalized(
        {
            key: value
            for key, value in target.items()
            if snake_case(key)
            in {
                "name",
                "category",
                "target_sequence",
                "external_identifiers",
                "uniprot_id_from_mapped_metadata",
            }
        }
    )
    sequence = result.get("target_sequence", {})
    for key in list(sequence):
        if key not in {"sequence_type", "sequence", "label", "taxonomy"}:
            del sequence[key]
    taxonomy = sequence.get("taxonomy", {})
    for key in list(taxonomy):
        if key not in {"code", "tax_id", "organism_name"}:
            del taxonomy[key]
    taxid = taxonomy.pop("code", taxonomy.get("tax_id"))
    legacy = taxonomy.pop("tax_id", None)
    if taxid is not None:
        taxonomy["taxid"] = str(taxid)
    if legacy is not None and str(legacy) != str(taxid):
        taxonomy["reported_taxid"] = str(legacy)
    mapped = result.pop("uniprot_id_from_mapped_metadata", None)
    identifiers = result.get("external_identifiers", [])
    for external in identifiers:
        for key in list(external):
            if key not in {"identifier", "offset"}:
                del external[key]
        identifier = external.get("identifier", {})
        for key in list(identifier):
            if key not in {"db_name", "identifier"}:
                del identifier[key]
    matches = [
        external
        for external in identifiers
        if str(external.get("identifier", {}).get("db_name", "")).lower() == "uniprot"
    ]
    if mapped:
        same = next(
            (
                external
                for external in matches
                if external.get("identifier", {}).get("identifier") == mapped
            ),
            None,
        )
        if same is not None:
            same["sources"] = ["external_identifiers", "mapped_metadata"]
        elif not matches:
            identifiers.append(
                {
                    "identifier": {"db_name": "UniProt", "identifier": mapped},
                    "source": "mapped_metadata",
                }
            )
            result["external_identifiers"] = identifiers
        else:
            result["uniprot_id_from_mapped_metadata"] = mapped
    return result


def analysis_target(
    target: dict[str, object],
    dataset_id: str,
    experiment_id: str,
    experiment_set_id: str,
) -> dict[str, object]:
    """Build the small target view used by the R-compatible analysis functions."""
    sequence = target.get("target_sequence", {})
    taxonomy = sequence.get("taxonomy", {})
    row = {
        "score_urn": dataset_id,
        "experiment_urn": experiment_id,
        "experiment_set_urn": experiment_set_id,
        "target_name": target.get("name"),
        "target_category": target.get("category"),
        "taxid": taxonomy.get("taxid"),
        "organism": taxonomy.get("organism_name"),
        "kingdom": taxonomy.get("kingdom"),
        "sequence_type": sequence.get("sequence_type"),
        "label": sequence.get("label"),
        "dna_sequence": None,
        "protein_sequence": None,
    }
    if sequence.get("sequence_type") in {"dna", "protein"}:
        row[f"{sequence['sequence_type']}_sequence"] = sequence.get("sequence")
    for database in ("uniprot", "ensembl", "refseq"):
        row[f"identifier_{database}"] = None
        row[f"offset_{database}"] = None
    for external in target.get("external_identifiers", []):
        identifier = external.get("identifier", {})
        database = str(identifier.get("db_name", "")).lower()
        if database in {"uniprot", "ensembl", "refseq"}:
            key = f"identifier_{database}"
            if row[key] is None:
                row[key] = identifier.get("identifier")
                offset = external.get("offset")
                try:
                    row[f"offset_{database}"] = (
                        int(offset) if offset is not None else None
                    )
                except (ValueError, TypeError):
                    LOGGER.warning(
                        "%s: invalid %s offset %r", dataset_id, database, offset
                    )
    return row


def load_datasets(input_path: Path) -> dict[str, dict[str, object]]:
    """Read only source fields extracted by the MAVE experiments or ProtVar scripts."""
    with input_path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    experiment_sets = data if isinstance(data, list) else data["experimentSets"]
    result = {}
    for experiment_set in experiment_sets:
        for experiment in experiment_set.get("experiments", []):
            for score_set in experiment.get("scoreSets", []):
                dataset_id = score_set["urn"]
                if dataset_id in result:
                    raise ValueError(f"Duplicate dataset_id in main.json: {dataset_id}")
                experiment_id = experiment.get("urn") or dataset_id.rsplit("-", 1)[0]
                experiment_set_id = (
                    experiment_set.get("urn") or dataset_id.rsplit("-", 2)[0]
                )
                targets = [
                    normalize_target(target)
                    for target in score_set.get("targetGenes") or []
                ]
                metadata = {
                    "dataset_id": dataset_id,
                    "experiment_id": experiment_id,
                    "experiment_set_id": experiment_set_id,
                    **{
                        SOURCE_METADATA_FIELDS[snake_case(key)]: value
                        for key, value in score_set.items()
                        if snake_case(key) in SOURCE_METADATA_FIELDS and present(value)
                    },
                    "targets": targets,
                }
                result[dataset_id] = {
                    "metadata": metadata,
                    "targets": [
                        analysis_target(
                            target, dataset_id, experiment_id, experiment_set_id
                        )
                        for target in targets
                    ],
                }
    return result


def load_metadata(input_path: Path) -> dict[str, list[dict[str, object]]]:
    """Compatibility view of the shared dataset reader."""
    return {key: value["targets"] for key, value in load_datasets(input_path).items()}


def download_taxonomy(destination: Path) -> None:
    """Download NCBI's taxcat archive and atomically cache categories.dmp."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(
        TAXONOMY_URL,
        headers={"User-Agent": "PerturbationCatalogue MaveDB processor"},
    )
    with tempfile.TemporaryDirectory(
        prefix=".taxcat-", dir=destination.parent
    ) as temporary_directory:
        temporary_path = Path(temporary_directory)
        archive_path = temporary_path / "taxcat.tar.gz"
        taxonomy_path = temporary_path / "taxcat.dmp"
        with urlopen(request, timeout=60) as response, archive_path.open(
            "wb"
        ) as archive_file:
            shutil.copyfileobj(response, archive_file)

        try:
            with tarfile.open(archive_path, mode="r:gz") as archive:
                members = [
                    member
                    for member in archive.getmembers()
                    if Path(member.name).name == "categories.dmp" and member.isfile()
                ]
                if len(members) != 1:
                    raise ValueError(
                        "NCBI taxcat archive must contain one categories.dmp file"
                    )
                source = archive.extractfile(members[0])
                if source is None:
                    raise ValueError("Could not read categories.dmp from NCBI archive")
                with source, taxonomy_path.open("wb") as output_file:
                    shutil.copyfileobj(source, output_file)
        except tarfile.TarError as error:
            raise ValueError(f"Invalid NCBI taxcat archive: {error}") from error

        os.replace(taxonomy_path, destination)


def resolve_taxonomy(taxonomy_path: Path | None) -> Path:
    """Use an explicit taxonomy file or download it to the local tmp cache."""
    if taxonomy_path is not None:
        if not taxonomy_path.is_file():
            raise FileNotFoundError(f"Taxonomy file does not exist: {taxonomy_path}")
        return taxonomy_path

    if not DEFAULT_TAXONOMY.is_file():
        LOGGER.info("Downloading NCBI taxonomy categories to %s", DEFAULT_TAXONOMY)
        download_taxonomy(DEFAULT_TAXONOMY)
    return DEFAULT_TAXONOMY


def add_kingdoms(
    metadata: dict[str, list[dict[str, object]]], taxonomy_path: Path
) -> None:
    """Join NCBI taxcat.dmp by its third column, as in the R workflow."""
    wanted = {row["taxid"] for targets in metadata.values() for row in targets}
    kingdoms = {
        "A": "Archaea",
        "B": "Bacteria",
        "E": "Eukaryota",
        "V": "Virus",
        "U": "Unclassified",
        "O": "Other",
    }
    lookup = {}
    with taxonomy_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3 and fields[2] in wanted:
                lookup[fields[2]] = kingdoms.get(fields[0], fields[0])
    for targets in metadata.values():
        for row in targets:
            row["kingdom"] = lookup.get(row["taxid"], row["kingdom"])


def read_scores(path: Path) -> pd.DataFrame:
    """Read all variant records; retain missing/invalid scores for the main export."""
    scores = pd.read_csv(
        path,
        usecols=lambda column: column
        in {"accession", "hgvs_nt", "hgvs_pro", "hgvs_splice", "score"},
        dtype="string",
        keep_default_na=False,
        na_values=["", "NA"],
        encoding="utf-8-sig",
    )
    if "score" not in scores:
        raise ValueError(f"{path}: missing required score column")
    raw = scores["score"].str.strip()
    numeric = pd.to_numeric(raw, errors="coerce").astype(float)
    invalid = (
        numeric.isna() & raw.notna() & ~raw.isin(["None", "NaN", "nan", "null", "NULL"])
    )
    scores["score_raw"] = raw.where(invalid)
    if invalid.any():
        LOGGER.warning(
            "%s: retaining %d invalid scores outside analysis", path.name, invalid.sum()
        )
    scores["score"] = numeric
    for column in ("accession", "hgvs_nt", "hgvs_splice", "hgvs_pro"):
        if column not in scores:
            scores[column] = pd.Series(pd.NA, index=scores.index, dtype="string")
    scores["hgvs_pro"] = scores["hgvs_pro"].str.strip()
    scores["source_row"] = np.arange(1, len(scores) + 1)
    accession = scores["accession"].str.strip()
    parts = accession.str.extract(r"^(.+?)#([0-9]+)$")
    scores["variant_num"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
    scores["source_accession"] = accession.where(parts[0].isna())
    scores["source_prefix"] = parts[0]
    parsed = scores["hgvs_pro"].str.extract(SIMPLE_P_PATTERN)
    scores["is_simple_p"] = parsed[0].notna()
    scores["ref_aa_abbr"] = parsed[0]
    scores["position"] = pd.to_numeric(parsed[1], errors="coerce").astype("Int64")
    scores["alt_aa_abbr"] = parsed[2]
    scores["ref_aa_sym"] = parsed[0].map(AA_3_TO_1)
    scores["alt_aa_sym"] = parsed[2].where(parsed[2].ne("="), parsed[0]).map(AA_3_TO_1)
    return scores


def score_summary(scores: pd.Series) -> dict[str, float]:
    """Summarize unaggregated measurements using sample standard deviation."""
    return dict(
        zip(
            SCORE_METRICS,
            (scores.mean(), scores.median(), scores.std(), scores.min(), scores.max()),
        )
    )


def with_identity(
    frame: pd.DataFrame, metadata: dict[str, object], columns: list[str]
) -> pd.DataFrame:
    """Attach constant target identity columns after aggregating score rows."""
    for column in IDENTITY_COLUMNS:
        frame[column] = metadata.get(column)
    return frame.reindex(columns=columns)


def make_singles(scores: pd.DataFrame, metadata: dict[str, object]) -> pd.DataFrame:
    """Collapse replicate single substitutions, synonymous calls and deletions."""
    hgvs = scores["hgvs_pro"]
    eligible = ~hgvs.isin(WT_HGVS_CODES) & ~hgvs.str.contains(
        r"[;?_]|fs\*?[0-9]*$", na=True
    )
    parsed = hgvs[eligible].str.extract(SINGLE_PATTERN)
    parsed.columns = ["wt", "position", "mut"]
    parsed = parsed[parsed["mut"].str.len().le(3).fillna(False)].copy()
    parsed["mut"] = parsed["mut"].mask(parsed["mut"].eq("="), parsed["wt"])
    parsed["wt"] = parsed["wt"].map(AA_3_TO_1)
    parsed["mut"] = parsed["mut"].map(AA_3_TO_1)
    parsed["position"] = pd.to_numeric(parsed["position"]) + (
        metadata.get("offset_uniprot") or 0
    )
    parsed["score"] = scores.loc[parsed.index, "score"]
    singles = parsed.groupby(
        ["position", "wt", "mut"], dropna=False, as_index=False
    ).agg(score=("score", "mean"), n=("score", "size"), sd=("score", "std"))
    return with_identity(singles, metadata, SINGLE_COLUMNS)


def make_doubles(scores: pd.DataFrame, metadata: dict[str, object]) -> pd.DataFrame:
    """Keep all calls with at most two components in sets containing multi mutants.

    Like R's align_start parsing, unparseable fields remain NA. This table also
    includes single mutants and WT controls from those score sets.
    """
    eligible = scores["hgvs_pro"].str.count(";").le(1)
    doubles = (
        scores.loc[eligible]
        .groupby("hgvs_pro", as_index=False)
        .agg(score=("score", "mean"), n=("score", "size"), sd=("score", "std"))
    )
    parsed = doubles["hgvs_pro"].str.extract(DOUBLE_PATTERN)
    parsed.columns = "wt1 position1 mut1 wt2 position2 mut2".split()
    for column in ("position1", "position2"):
        parsed[column] = pd.to_numeric(parsed[column], errors="coerce") + (
            metadata.get("offset_uniprot") or 0
        )
    for column in ("wt1", "wt2", "mut1", "mut2"):
        parsed[column] = parsed[column].map(AA_3_TO_1)
    doubles = pd.concat([doubles, parsed], axis=1).sort_values(
        "position1 wt1 mut1 position2 wt2 mut2".split(), na_position="last"
    )
    return with_identity(doubles, metadata, DOUBLE_COLUMNS)


def make_wide(singles: pd.DataFrame, metadata: dict[str, object]) -> pd.DataFrame:
    """Pivot canonical single mutants and calculate per-position statistics."""
    valid = singles.dropna(subset=["position", "wt"])
    valid = valid[valid["mut"].isin(MUTATIONS)]
    if valid.empty:
        return pd.DataFrame(columns=WIDE_COLUMNS)
    wide = valid.pivot(
        index=["position", "wt"], columns="mut", values=["score", "n", "sd"]
    )
    wide = wide.reindex(
        columns=pd.MultiIndex.from_product([["score", "n", "sd"], MUTATIONS])
    )
    wide.columns = [f"{mut}_{metric}" for metric, mut in wide.columns]
    wide = wide.reset_index()
    aa_scores = wide[[f"{aa}_score" for aa in AMINO_ACIDS]]
    wide["missing"] = wide[[f"{mut}_score" for mut in MUTATIONS]].isna().sum(axis=1)
    # Look up the reference amino acid by name, independently of column order.
    wide["synonymous_score"] = [
        wide.at[index, f"{wt}_score"] if wt in AMINO_ACIDS else np.nan
        for index, wt in wide["wt"].items()
    ]
    wide["median_score"] = aa_scores.median(axis=1)
    wide["mean_score"] = aa_scores.mean(axis=1)
    wide["sd_score"] = aa_scores.std(axis=1)
    for mutation in ("stop", "synonymous", "del"):
        wide[mutation] = wide[f"{mutation}_score"].notna()
    return with_identity(wide, metadata, WIDE_COLUMNS)


def position_summary(wide: pd.DataFrame) -> dict[str, object]:
    """Compute the R completeness and coverage metrics across targeted positions."""
    if wide.empty:
        return {}
    completeness = (
        22 - wide["missing"] - wide["stop"].astype(int) - wide["del"].astype(int)
    ) / 20
    missense = (20 * completeness - wide["synonymous"].astype(int)) / 19
    start, end = wide["position"].min(), wide["position"].max()
    length = end - start + 1
    result = {
        "region_start": start,
        "region_end": end,
        "region_length": length,
        "completeness": completeness.mean(),
        "completeness_missense": missense.mean(),
        "mean_position_median_score": wide["median_score"].mean(),
        "coverage_nonsense": wide["stop"].sum() / length,
        "coverage_synonymous": wide["synonymous"].sum() / length,
    }
    for mutation in ("stop", "del", "synonymous"):
        result[f"{mutation}_proportion"] = wide[mutation].mean()
    for threshold, selected in (
        (50, missense > 0.5),
        (75, missense > 0.75),
        (100, missense == 1),
    ):
        result[f"positions_{threshold}"] = int(selected.sum())
        result[f"coverage_{threshold}"] = selected.sum() / length
    return result


def multiple_summary(hgvs: pd.Series) -> dict[str, float]:
    """Count mutations per measurement and distinct variant backgrounds."""
    parts = hgvs.str.replace(r"p.\[|\]", "", regex=True).str.split(";")
    counts = parts.str.len()
    backgrounds = pd.DataFrame({"hgvs": hgvs, "mutation": parts}).explode("mutation")
    return {
        "mean_muts": counts.mean(),
        "max_muts": counts.max(),
        "min_muts": counts.min(),
        "mean_backgrounds": backgrounds.groupby("mutation")["hgvs"].nunique().mean(),
    }


def translate_dna(sequence: object) -> str | None:
    """Translate complete in-frame DNA with the standard genetic code."""
    if not isinstance(sequence, str) or len(sequence) % 3:
        return None
    try:
        return str(Seq(sequence).translate())
    except TranslationError as error:
        LOGGER.warning("Could not translate target DNA: %s", error)
        return None


def any_with_na(values: pd.Series) -> bool | None:
    """Match R's any(): TRUE wins; otherwise retain missing logical values."""
    if values.any():
        return True
    return None if values.isna().any() else False


def nucleotide_summary(scores: pd.DataFrame, category: str) -> dict[str, object]:
    """Summarize MPRA or splicing nucleotide calls without protein parsing."""
    hgvs = scores["hgvs_nt"]
    positions = hgvs.str.findall(r"[0-9]+").explode().dropna()
    missing = hgvs.isna().any()
    start = int(positions.astype(int).min()) if len(positions) and not missing else None
    end = int(positions.astype(int).max()) if len(positions) and not missing else None
    result = {
        **score_summary(scores["score"]),
        "n_muts": hgvs.str.split(";").explode().nunique(dropna=False),
        "n_positions": positions.nunique() + int(missing),
        "region_start": start,
        "region_end": end,
        "region_length": end - start + 1 if start is not None else None,
        "wt_mutants": any_with_na(scores["wt_mut"]),
        "single_mutations": any_with_na(
            hgvs.str.contains(
                r"^n.\d+[ACTG]>[ACTG]"
                if category == "mpras"
                else r"^g.\d+[ACTG]>[ACTG]"
            )
        ),
        "multiple_mutations": any_with_na(hgvs.str.contains(";")),
    }
    if category == "mpras":
        result["delins_mutations"] = any_with_na(hgvs.str.contains("delins"))
    return result


def json_safe(value: object) -> object:
    """Make nested metadata valid JSON, including numpy scalars and missing values."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    return (
        None
        if value is pd.NA or (isinstance(value, float) and not np.isfinite(value))
        else value
    )


def write_tsv(
    frame: pd.DataFrame,
    destination: Path | TextIO,
    *,
    header: bool = True,
    quote: bool = False,
    single_line: bool = False,
    max_cell_chars: int | None = None,
) -> None:
    """Use R-compatible NA, logical values and integer output columns."""
    frame = frame.copy()
    for column in frame.columns:
        if single_line:
            # Escape scalar text before JSON serialization, so JSON cells remain
            # valid JSON while every dataset occupies exactly one physical line.
            frame[column] = frame[column].map(
                lambda value: (
                    value.translate(TSV_TEXT_ESCAPES)
                    if isinstance(value, str)
                    else value
                )
            )
        structured = (
            frame[column].map(lambda value: isinstance(value, (dict, list)))
            if frame[column].dtype == object
            else pd.Series(False, index=frame.index)
        )
        if column in BOOLEAN_COLUMNS or pd.api.types.is_bool_dtype(frame[column]):
            frame[column] = frame[column].astype(object)
            frame.loc[~structured, column] = frame.loc[~structured, column].map(
                {True: "TRUE", False: "FALSE"}
            )
        elif column in INTEGER_COLUMNS:
            scalar = frame.loc[~structured, column].astype("Int64")
            frame[column] = frame[column].astype(object)
            frame.loc[~structured, column] = scalar.astype(object)
        if structured.any():
            frame[column] = frame[column].map(
                lambda value: (
                    json.dumps(
                        json_safe(value), ensure_ascii=False, separators=(",", ":")
                    )
                    if isinstance(value, (dict, list))
                    else value
                )
            )
        if max_cell_chars is not None:
            lengths = frame[column].map(
                lambda value: len(value) if isinstance(value, str) else 0
            )
            if lengths.gt(max_cell_chars).any():
                raise ValueError(
                    f"Metadata column {column!r} exceeds {max_cell_chars:,} characters; "
                    "no values were truncated"
                )
    frame.to_csv(
        destination,
        sep="\t",
        index=False,
        header=header,
        na_rep="NA",
        quoting=csv.QUOTE_ALL if quote else csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )


def category_scores(
    scores: pd.DataFrame, targets: list[dict[str, object]]
) -> dict[str, list[tuple[int, pd.DataFrame]]]:
    """Decide dataset eligibility once using the original analysis selection rules."""
    scored = scores.dropna(subset=["score"])
    result: dict[str, list[tuple[int, pd.DataFrame]]] = {}
    for index, target in enumerate(targets):
        category = target["target_category"]
        if category == "protein_coding":
            protein = scored.dropna(subset=["hgvs_pro"]).copy()
            if not protein.empty:
                protein["wt_mut"] = protein["hgvs_pro"].isin(WT_HGVS_CODES)
                result.setdefault("maves", []).append((index, protein))
        if category == "regulatory":
            nucleotide = scored.copy()
            if not nucleotide.empty:
                nucleotide["wt_mut"] = nucleotide["hgvs_nt"].isin(WT_HGVS_CODES)
                result.setdefault("mpras", []).append((index, nucleotide))
        if category in {"protein_coding", "other_noncoding"}:
            splice = (
                scored.dropna(subset=["hgvs_splice"]).copy()
                if category == "protein_coding"
                else scored.copy()
            )
            if not splice.empty:
                splice["wt_mut"] = splice["hgvs_nt"].str.contains("A>A|C>C|G>G|T>T")
                result.setdefault("splicing", []).append((index, splice))
    return result


def mave_analysis(
    scores: pd.DataFrame, target: dict[str, object]
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    """Calculate the corrected R metrics without duplicating source variant rows."""
    hgvs = scores["hgvs_pro"]
    counts = hgvs.str.count(";") + 1
    singles = make_singles(scores, target)
    wide = make_wide(singles, target)
    summary = score_summary(scores["score"])
    sd = summary["sd_score"]
    outliers = (
        int(((scores["score"] - summary["mean_score"]).abs() > 2 * sd).sum())
        if pd.notna(sd)
        else None
    )
    summary.update(
        n_analyzed_variants=len(scores),
        n_subs=int(counts.sum()),
        n_positions=hgvs.str.findall(r"[0-9]+").explode().nunique(),
        multiple_mutations=bool((counts > 1).any()),
        single_mutations=bool((counts == 1).any()),
        wt_mutants=bool(scores["wt_mut"].any()),
        n_outlier=outliers,
        p_outlier=outliers / len(scores) if outliers is not None else None,
        mean_variant_reps=singles["n"].mean(),
    )
    summary.update(position_summary(wide))
    if summary["multiple_mutations"]:
        summary.update(multiple_summary(hgvs))
    return summary, singles, wide


def dataset_analysis(
    dataset: dict[str, object],
    scores: pd.DataFrame,
    categories: dict[str, list[tuple[int, pd.DataFrame]]],
    compatibility: dict[str, TextIO],
) -> dict[str, object]:
    """One dataset row with primary-analysis metrics and source target metadata."""
    row = dict(dataset["metadata"])
    targets = dataset["targets"]
    source_targets = row["targets"]
    primary = next(category for category in ANALYSIS_TYPES if category in categories)
    row.update(
        analysis_types=[
            category for category in ANALYSIS_TYPES if category in categories
        ],
        primary_analysis=primary,
        meta_study=bool(re.search(r"urn:mavedb:\d*-0-", row["dataset_id"])),
        n_variants=len(scores),
        n_scored_variants=int(scores["score"].notna().sum()),
        n_invalid_scores=int(scores["score_raw"].notna().sum()),
        n_missing_scores=int(
            (scores["score"].isna() & scores["score_raw"].isna()).sum()
        ),
    )
    row.update(score_summary(scores["score"].dropna()))
    row["n_analyzed_variants"] = len(
        set().union(*(set(subset.index) for _, subset in categories[primary]))
    )
    for index, target in enumerate(targets):
        source = source_targets[index]
        if target.get("kingdom"):
            source.setdefault("target_sequence", {}).setdefault("taxonomy", {})[
                "kingdom"
            ] = target["kingdom"]
        if (
            target["target_category"] == "protein_coding"
            and target["protein_sequence"] is None
            and target["dna_sequence"] is not None
        ):
            translated = translate_dna(target["dna_sequence"])
            if translated is not None:
                source["derived_protein_sequence"] = translated
    for category in ANALYSIS_TYPES:
        analyses = categories.get(category, [])
        for index, subset in analyses:
            target = targets[index]
            if category == "maves":
                metrics, singles, wide = mave_analysis(subset, target)
                if compatibility:
                    write_tsv(singles, compatibility["singles"], header=False)
                    write_tsv(wide, compatibility["wide"], header=False)
                    if metrics["multiple_mutations"]:
                        write_tsv(
                            make_doubles(subset, target),
                            compatibility["doubles"],
                            header=False,
                        )
            else:
                metrics = {
                    "n_analyzed_variants": len(subset),
                    **nucleotide_summary(subset, category),
                }
            metrics["position_reference"] = (
                "uniprot"
                if category == "maves" and target.get("offset_uniprot") is not None
                else "source"
            )
            if category == primary and len(analyses) == 1:
                row.update(metrics)
    return row


def unpack_targets(row: dict[str, object]) -> dict[str, object]:
    """Use scalar target columns, or aligned arrays when a dataset has several targets."""
    result = {key: value for key, value in row.items() if key != "targets"}
    targets = []
    for target in row.get("targets", []):
        flat = {}
        for column, path in (TARGET_SOURCE_FIELDS | TARGET_DERIVED_FIELDS).items():
            value = target
            for key in path:
                value = value.get(key) if isinstance(value, dict) else None
            if column in TARGET_SOURCE_FIELDS or value is not None:
                flat[column] = value
        for database in TARGET_DATABASES:
            identifiers = [
                external
                for external in target.get("external_identifiers", [])
                if str(external.get("identifier", {}).get("db_name", "")).lower()
                == database
            ]
            values = {
                f"identifier_{database}": [
                    external["identifier"].get("identifier") for external in identifiers
                ],
                f"offset_{database}": [
                    external.get("offset") for external in identifiers
                ],
                f"identifier_{database}_sources": [
                    external.get("sources")
                    or ([external["source"]] if external.get("source") else None)
                    for external in identifiers
                ],
            }
            for column, items in values.items():
                if not column.endswith("_sources") or any(
                    item is not None for item in items
                ):
                    flat[column] = items[0] if len(items) == 1 else items or None
        targets.append(flat)
    for column in dict.fromkeys(key for target in targets for key in target):
        values = [target.get(column) for target in targets]
        result[column] = values[0] if len(values) == 1 else values
    return result


def metadata_frame(rows: list[dict[str, object]]) -> pd.DataFrame:
    """Expand target metadata, keeping dataset-level derived fields last."""
    frame = pd.DataFrame([unpack_targets(row) for row in rows])
    target_derived = list(TARGET_DERIVED_FIELDS) + [
        f"identifier_{database}_sources" for database in TARGET_DATABASES
    ]
    columns = (
        SOURCE_COLUMNS
        + DERIVED_COLUMNS
        + [column for column in target_derived if column in frame]
    )
    return frame.reindex(columns=columns)


def resolve_metadata_path(input_path: Path) -> Path:
    """Accept main.json, the dump, its parent extraction directory or csv/."""
    if not input_path.is_dir():
        return input_path
    candidates = [input_path / "main.json", input_path / "mavedb_dump/main.json"]
    if input_path.name == "csv":
        candidates.append(input_path.parent / "main.json")
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"No main.json found under {input_path}")


def process_mavedb(
    input_path: Path,
    output_path: Path,
    *,
    scores_path: Path | None = None,
    taxonomy_path: Path | None = None,
    compatibility_tables: bool = False,
) -> dict[str, int]:
    """Publish two canonical tables, with the same dataset selection for both."""
    metadata_path = resolve_metadata_path(input_path)
    scores_path = (
        scores_path if scores_path is not None else metadata_path.parent / "csv"
    )
    score_files = sorted(scores_path.glob("*.scores.csv"))
    if not score_files:
        raise FileNotFoundError(f"No *.scores.csv files found in {scores_path}")
    datasets = load_datasets(metadata_path)
    target_metadata = {key: value["targets"] for key, value in datasets.items()}
    add_kingdoms(target_metadata, resolve_taxonomy(taxonomy_path))
    output_path.mkdir(parents=True, exist_ok=True)
    rows = []
    score_count = 0
    with tempfile.TemporaryDirectory(
        prefix=".processing-", dir=output_path
    ) as temporary:
        stage = Path(temporary)
        with ExitStack() as stack:
            score_handle = stack.enter_context(
                (stage / "scores.tsv").open("w", encoding="utf-8", newline="")
            )
            write_tsv(pd.DataFrame(columns=SCORE_COLUMNS), score_handle)
            compatibility = {}
            if compatibility_tables:
                directory = stage / "maves"
                directory.mkdir()
                for name, columns in [
                    ("singles", SINGLE_COLUMNS),
                    ("doubles", DOUBLE_COLUMNS),
                    ("wide", WIDE_COLUMNS),
                ]:
                    handle = stack.enter_context(
                        (directory / f"{name}.tsv").open(
                            "w", encoding="utf-8", newline=""
                        )
                    )
                    compatibility[name] = handle
                    write_tsv(pd.DataFrame(columns=columns), handle)
            for index, score_file in enumerate(score_files, start=1):
                dataset_id = score_file.name.removesuffix(".scores.csv").replace(
                    "urn-mavedb-", "urn:mavedb:", 1
                )
                if dataset_id in EXCLUDED_URNS:
                    continue
                dataset = datasets.get(dataset_id)
                if dataset is None or not dataset["targets"]:
                    LOGGER.warning("%s: no target metadata; skipping", dataset_id)
                    continue
                scores = read_scores(score_file)
                categories = category_scores(scores, dataset["targets"])
                if not categories:
                    continue
                scores["dataset_id"] = dataset_id
                # The legacy name is internal to the optional R compatibility functions.
                scores["score_urn"] = dataset_id
                row = dataset_analysis(dataset, scores, categories, compatibility)
                prefixes = scores["source_prefix"].dropna().unique().tolist()
                aliases = [prefix for prefix in prefixes if prefix != dataset_id]
                if len(prefixes) == 1 and aliases:
                    row["source_accession_prefix"] = aliases[0]
                elif aliases:
                    mismatched = scores["source_prefix"].notna() & scores[
                        "source_prefix"
                    ].ne(dataset_id)
                    scores.loc[mismatched, "source_accession"] = scores.loc[
                        mismatched, "accession"
                    ]
                write_tsv(
                    scores.reindex(columns=SCORE_COLUMNS), score_handle, header=False
                )
                rows.append(row)
                score_count += len(scores)
                if index % 100 == 0:
                    LOGGER.info("Processed %d/%d score files", index, len(score_files))
        write_tsv(
            metadata_frame(rows),
            stage / "metadata.tsv",
            quote=True,
            single_line=True,
            max_cell_chars=32767,
        )
        for name in ("scores.tsv", "metadata.tsv"):
            os.replace(stage / name, output_path / name)
        # Remove only files managed by the previous version of this processor.
        # Other files in the output directory are never touched.
        for category in ANALYSIS_TYPES:
            directory = output_path / category
            for name in ("full.tsv", "meta.tsv", "meta_columns.tsv"):
                (directory / name).unlink(missing_ok=True)
            if category == "maves":
                for name in ("singles.tsv", "doubles.tsv", "wide.tsv"):
                    destination = directory / name
                    if compatibility_tables:
                        directory.mkdir(exist_ok=True)
                        os.replace(stage / "maves" / name, destination)
                    else:
                        destination.unlink(missing_ok=True)
            if directory.is_dir() and not any(directory.iterdir()):
                directory.rmdir()
    return {"datasets": len(rows), "scores": score_count}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Dump directory or main.json path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Directory for scores.tsv, metadata.tsv and optional compatibility tables.",
    )
    parser.add_argument(
        "--scores", type=Path, help="Override the dump's csv/ directory."
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        help="Use this existing NCBI taxcat.dmp instead of the cached/downloaded file.",
    )
    parser.add_argument(
        "--compatibility-tables",
        action="store_true",
        help="Also generate maves/singles.tsv, doubles.tsv and wide.tsv in the R layout.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        counts = process_mavedb(
            args.input,
            args.output,
            scores_path=args.scores,
            taxonomy_path=args.taxonomy,
            compatibility_tables=args.compatibility_tables,
        )
    except (OSError, ValueError, KeyError, pd.errors.ParserError) as error:
        parser.error(str(error))
    LOGGER.info("Wrote %s to %s", counts, args.output)


if __name__ == "__main__":
    main()
