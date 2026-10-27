#!/usr/bin/env python3
"""Reproduce the MAVE, MPRA and splicing TSV analysis of process_mavedb_data.R.

Read a Zenodo dump (main.json and csv/*.scores.csv), processing one score set at
a time. Scores are not normalized. Standard deviations use the sample formula
(n - 1), and missing output values are written as NA. If the local taxonomy
cache is absent, the NCBI category dump is downloaded to data_sources/mavedb/tmp.
No FASTA files are made.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
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
FULL_COLUMNS = "score_urn hgvs_nt hgvs_splice hgvs_pro wt_mut score".split()
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
MAVE_META_COLUMNS = (
    "experiment_set_urn experiment_urn score_urn meta_study title desc abstract target_category "
    "kingdom organism taxid target_name label identifier_uniprot identifier_ensembl identifier_refseq "
    "offset_uniprot offset_ensembl offset_refseq dna_sequence protein_sequence "
    "single_mutations multiple_mutations wt_mutants stop_proportion del_proportion synonymous_proportion "
    "n n_subs n_positions mean_variant_reps max_muts min_muts mean_muts mean_backgrounds "
    "region_start region_end region_length n_outlier p_outlier mean_score median_score sd_score "
    "min_score max_score mean_position_median_score completeness completeness_missense "
    "positions_50 positions_75 positions_100 coverage_50 coverage_75 coverage_100 "
    "coverage_nonsense coverage_synonymous"
).split()
NUCLEOTIDE_META_PREFIX = (
    "experiment_set_urn experiment_urn score_urn meta_study title desc abstract target_category "
    "kingdom organism taxid target_name label identifier_ensembl identifier_refseq identifier_uniprot "
    "offset_ensembl offset_refseq offset_uniprot dna_sequence"
).split()
NUCLEOTIDE_META_SUFFIX = (
    "n n_muts n_positions region_start region_end region_length wt_mutants "
    "single_mutations multiple_mutations"
).split()
SCORE_METRICS = "mean_score median_score sd_score min_score max_score".split()
MPRA_META_COLUMNS = (
    NUCLEOTIDE_META_PREFIX
    + NUCLEOTIDE_META_SUFFIX
    + ["delins_mutations"]
    + SCORE_METRICS
)
SPLICING_META_COLUMNS = (
    NUCLEOTIDE_META_PREFIX
    + ["protein_sequence"]
    + NUCLEOTIDE_META_SUFFIX
    + SCORE_METRICS
)
META_COLUMNS = {
    "maves": MAVE_META_COLUMNS,
    "mpras": MPRA_META_COLUMNS,
    "splicing": SPLICING_META_COLUMNS,
}
BOOLEAN_COLUMNS = {
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


COMMON_COLUMN_DESCRIPTIONS = {
    "experiment_set_urn": ("str", "URN for a group of linked experiments"),
    "experiment_urn": ("str", "URN for a single assay"),
    "score_urn": ("str", "URN for a score set calculated against an experiment"),
    "meta_study": ("bool", "Whether the score set is raw data or a meta analysis set"),
    "title": ("str", "Paper title"),
    "desc": ("str", "MaveDB description of the entry"),
    "abstract": ("str", "Paper abstract"),
    "target_category": ("str", "Target category"),
    "kingdom": ("str", "Phyogenetic kingdom of target"),
    "organism": ("str", "Organism of target"),
    "taxid": ("str", "Organism NCBI TaxID"),
    "target_name": ("str", "Target gene name"),
    "label": ("str", "Target label/additional description"),
    "identifier_uniprot": ("str", "Target Uniprot ID"),
    "identifier_ensembl": ("str", "Target Ensembl ID"),
    "identifier_refseq": ("str", "Target RefSeq ID"),
    "offset_uniprot": ("int", "Offset between protein_sequence and Uniprot"),
    "offset_ensembl": ("int", "Offset between dna_sequence and Ensembl"),
    "offset_refseq": ("int", "Offset between dna_sequence and RefSeq"),
    "dna_sequence": ("str", "Target nucleotide sequence"),
    "protein_sequence": ("str", "Target AA sequence"),
    "single_mutations": ("bool", "Whether the study has single mutations"),
    "multiple_mutations": ("bool", "Whether the study has multiple mutations"),
    "wt_mutants": (
        "bool",
        "Whether the study has any variants marked as WT/fully synonymous (note this is p.= and different to synonymous proportion)",
    ),
    "stop_proportion": ("float", "Proportion of positions with a nonsense mutation"),
    "del_proportion": (
        "float",
        "Proportion of positions with a single AA deletion mutation",
    ),
    "synonymous_proportion": (
        "float",
        "Proportion of positions with an explicit p.Met1Met/p.Met1= style synonymous substitution",
    ),
    "n": ("int", "Number of variants"),
    "n_subs": ("int", "Number of distinct substitutions"),
    "n_positions": ("int", "Number of positions covered"),
    "mean_variant_reps": ("float", "Mean number of measurements of each variant"),
    "min_muts": ("int", "Smallest number of variants in one measurement"),
    "max_muts": ("int", "Largest number of variants in one measurement"),
    "mean_muts": (
        "int",
        "Mean number of substitution per variant (NA for studies without multiple mutation)",
    ),
    "mean_backgrounds": (
        "int",
        "Mean number of background each variant is found in (NA for studies without multiple mutation)",
    ),
    "region_start": ("int", "First mutated position"),
    "region_end": ("int", "Final mutated position"),
    "region_length": (
        "int",
        "Length of region covered (not all subsection will be mutated though, compare n_positions)",
    ),
    "n_outlier": (
        "int",
        "Number of measurements more than 2 s.d. from the mean score (caution - MAVEs are not normally distributed, this is a rough measurement)",
    ),
    "p_outlier": (
        "float",
        "Proportion of measurements more than 2 s.d. from the mean score",
    ),
    "mean_score": ("float", "Mean fitness score"),
    "median_score": ("float", "Median fitness score"),
    "sd_score": ("float", "SD of fitness score"),
    "min_score": ("float", "Minimum fitness score"),
    "max_score": ("float", "Maximum fitness score"),
    "mean_position_median_score": (
        "float",
        "Mean value of the median fitness score at a position",
    ),
    "completeness": (
        "float",
        "Proportion of single substitutions measured at positions with at least one variant",
    ),
    "completeness_missense": (
        "float",
        "Proportion of missense substitutions measured at positions with at least one variant",
    ),
    "positions_50": (
        "int",
        "Number of positions with at least 50% of missense variants measured",
    ),
    "positions_75": (
        "int",
        "Number of positions with at least 75% of missense variants measured",
    ),
    "positions_100": ("int", "Number of positions with all missense variants measured"),
    "coverage_50": (
        "float",
        "positions_50/region_length (i.e. coverage across the whole region spanned, not just mutated positions)",
    ),
    "coverage_75": (
        "float",
        "positions_75/region_length (i.e. coverage across the whole region spanned, not just mutated positions)",
    ),
    "coverage_100": (
        "float",
        "positions_100/region_length (i.e. coverage across the whole region spanned, not just mutated positions)",
    ),
    "coverage_nonsense": (
        "float",
        "Proportion of region spanned with nonsense variants measured",
    ),
    "coverage_synonymous": (
        "float",
        "Proportion of region spanned with synonymous variants measured",
    ),
    "n_muts": ("int", "Number of distinct mutation tested in any combination"),
    "delins_mutations": ("bool", "Study includes complex delins mutations"),
}
COLUMN_DESCRIPTIONS = {
    "maves": {
        **COMMON_COLUMN_DESCRIPTIONS,
    },
    "mpras": {
        **COMMON_COLUMN_DESCRIPTIONS,
        "n": ("int", "Total variants measured"),
        "n_positions": ("int", "Number of positions mutated"),
        "region_start": ("int", "Start of the targeted region"),
        "region_end": ("int", "End of the targeted region"),
        "region_length": ("int", "Length of the targted region"),
        "wt_mutants": ("bool", "Study includes WT alleles"),
        "single_mutations": ("bool", "Study includes single substitutions"),
        "multiple_mutations": ("bool", "Study includes multiple mutations"),
        "sd_score": ("float", "SD of fitness scores"),
    },
    "splicing": {
        **COMMON_COLUMN_DESCRIPTIONS,
        "n": ("int", "Total variants measured"),
        "n_positions": ("int", "Number of positions mutated"),
        "region_start": ("int", "Start of the targeted region"),
        "region_end": ("int", "End of the targeted region"),
        "region_length": ("int", "Length of the targted region"),
        "wt_mutants": ("bool", "Study includes WT alleles"),
        "single_mutations": ("bool", "Study includes single substitutions"),
        "multiple_mutations": ("bool", "Study includes multiple mutations"),
        "sd_score": ("float", "SD of fitness scores"),
    },
}


def load_metadata(input_path: Path) -> dict[str, list[dict[str, object]]]:
    """Flatten target metadata, retaining targets with no external identifiers."""
    with input_path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    experiment_sets = data if isinstance(data, list) else data["experimentSets"]
    result: dict[str, list[dict[str, object]]] = {}
    for experiment_set in experiment_sets:
        for experiment in experiment_set.get("experiments", []):
            for score_set in experiment.get("scoreSets", []):
                urn = score_set["urn"]
                targets = result.setdefault(urn, [])
                for target in score_set.get("targetGenes") or []:
                    sequence = target.get("targetSequence") or {}
                    taxonomy = sequence.get("taxonomy") or {}
                    taxid = taxonomy.get("taxId", taxonomy.get("code"))
                    row = {
                        "score_urn": urn,
                        "experiment_set_urn": experiment_set.get("urn")
                        or urn.rsplit("-", 2)[0],
                        "experiment_urn": experiment.get("urn")
                        or urn.rsplit("-", 1)[0],
                        "meta_study": bool(re.search(r"urn:mavedb:\d*-0-", urn)),
                        "title": score_set.get("title"),
                        "desc": score_set.get("shortDescription"),
                        "abstract": score_set.get("abstractText"),
                        "n": score_set.get("numVariants"),
                        "target_name": target.get("name"),
                        "target_category": target.get("category"),
                        "label": sequence.get("label"),
                        "sequence_type": sequence.get("sequenceType"),
                        "dna_sequence": None,
                        "protein_sequence": None,
                        "taxid": str(taxid) if taxid is not None else None,
                        "organism": taxonomy.get("organismName"),
                        "kingdom": taxonomy.get("kingdom"),
                    }
                    if sequence.get("sequenceType") in {"dna", "protein"}:
                        row[f"{sequence['sequenceType']}_sequence"] = sequence.get(
                            "sequence"
                        )
                    for database in ("uniprot", "ensembl", "refseq"):
                        row[f"identifier_{database}"] = None
                        row[f"offset_{database}"] = None
                    for external in target.get("externalIdentifiers") or []:
                        identifier = external.get("identifier") or {}
                        database = (identifier.get("dbName") or "").lower()
                        if database in {"uniprot", "ensembl", "refseq"}:
                            key = f"identifier_{database}"
                            if row[key] is not None and row[key] != identifier.get(
                                "identifier"
                            ):
                                raise ValueError(
                                    f"{urn}: multiple {database} identifiers for one target"
                                )
                            row[key] = identifier.get("identifier")
                            offset = external.get("offset")
                            row[f"offset_{database}"] = (
                                int(offset) if offset is not None else None
                            )
                    targets.append(row)
    return result


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
    """Read only relevant columns and discard missing or unparseable scores."""
    scores = pd.read_csv(
        path,
        usecols=lambda column: column in FULL_COLUMNS,
        dtype="string",
        keep_default_na=False,
        na_values=["", "NA"],
        encoding="utf-8-sig",
    )
    if "score" not in scores:
        raise ValueError(f"{path}: missing required score column")
    numeric = pd.to_numeric(scores["score"], errors="coerce")
    invalid = (
        numeric.isna()
        & scores["score"].notna()
        & ~scores["score"].isin(["None", "NaN", "nan", "null", "NULL"])
    )
    if invalid.any():
        LOGGER.warning("%s: dropping %d invalid scores", path.name, invalid.sum())
    scores["score"] = numeric.astype(float)
    for column in ("hgvs_nt", "hgvs_splice", "hgvs_pro"):
        if column not in scores:
            scores[column] = pd.Series(pd.NA, index=scores.index, dtype="string")
    return scores.dropna(subset=["score"]).reset_index(drop=True)


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


def write_tsv(
    frame: pd.DataFrame,
    destination: Path | TextIO,
    *,
    header: bool = True,
    quote: bool = False,
) -> None:
    """Use R-compatible NA, logical values and integer output columns."""
    frame = frame.copy()
    for column in frame.columns:
        if column in BOOLEAN_COLUMNS:
            frame[column] = frame[column].map({True: "TRUE", False: "FALSE"})
        elif column in INTEGER_COLUMNS:
            frame[column] = frame[column].astype("Int64")
    frame.to_csv(
        destination,
        sep="\t",
        index=False,
        header=header,
        na_rep="NA",
        quoting=csv.QUOTE_ALL if quote else csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )


def process_mavedb(
    input_path: Path,
    output_path: Path,
    *,
    scores_path: Path | None = None,
    taxonomy_path: Path | None = None,
) -> dict[str, int]:
    """Write the requested TSVs; return the number of metadata rows by category."""
    metadata_path = input_path / "main.json" if input_path.is_dir() else input_path
    scores_path = (
        scores_path if scores_path is not None else metadata_path.parent / "csv"
    )
    score_files = sorted(scores_path.glob("*.scores.csv"))
    if not score_files:
        raise FileNotFoundError(f"No *.scores.csv files found in {scores_path}")
    taxonomy_path = resolve_taxonomy(taxonomy_path)
    metadata = load_metadata(metadata_path)
    add_kingdoms(metadata, taxonomy_path)
    summaries: dict[str, list[dict[str, object]]] = {
        category: [] for category in META_COLUMNS
    }
    with ExitStack() as stack:
        outputs = {}
        for category in META_COLUMNS:
            directory = output_path / category
            directory.mkdir(parents=True, exist_ok=True)
            schemas = {"full": FULL_COLUMNS}
            if category == "maves":
                schemas.update(
                    singles=SINGLE_COLUMNS, doubles=DOUBLE_COLUMNS, wide=WIDE_COLUMNS
                )
            for name, columns in schemas.items():
                handle = stack.enter_context(
                    (directory / f"{name}.tsv").open("w", encoding="utf-8", newline="")
                )
                outputs[category, name] = handle
                write_tsv(pd.DataFrame(columns=columns), handle)
        for index, score_file in enumerate(score_files, start=1):
            urn = score_file.name.removesuffix(".scores.csv").replace(
                "urn-mavedb-", "urn:mavedb:", 1
            )
            if urn in EXCLUDED_URNS:
                continue
            targets = metadata.get(urn)
            if not targets:
                LOGGER.warning("%s: no target metadata; skipping", urn)
                continue
            scores = read_scores(score_file)
            if scores.empty:
                continue
            scores["score_urn"] = urn
            for category in META_COLUMNS:
                if category == "maves":
                    selected = [
                        target
                        for target in targets
                        if target["target_category"] == "protein_coding"
                    ]
                    if not selected:
                        continue
                    subset = scores.dropna(subset=["hgvs_pro"]).copy()
                    subset["wt_mut"] = subset["hgvs_pro"].isin(WT_HGVS_CODES)
                elif category == "mpras":
                    selected = [
                        target
                        for target in targets
                        if target["target_category"] == "regulatory"
                    ]
                    if not selected:
                        continue
                    subset = scores.copy()
                    subset["wt_mut"] = subset["hgvs_nt"].isin(WT_HGVS_CODES)
                else:
                    selected = [
                        target
                        for target in targets
                        if target["target_category"]
                        in {"protein_coding", "other_noncoding"}
                    ]
                    if not selected or (
                        all(
                            target["target_category"] == "protein_coding"
                            for target in selected
                        )
                        and scores["hgvs_splice"].isna().all()
                    ):
                        continue
                    subset = scores.copy()
                    subset["wt_mut"] = subset["hgvs_nt"].str.contains("A>A|C>C|G>G|T>T")
                category_summaries = []
                all_singles, all_wide, all_measurements = [], [], []
                for target in selected:
                    target_scores = subset
                    if (
                        category == "splicing"
                        and target["target_category"] == "protein_coding"
                    ):
                        target_scores = subset.dropna(subset=["hgvs_splice"])
                    if target_scores.empty:
                        continue
                    write_tsv(
                        target_scores[FULL_COLUMNS],
                        outputs[category, "full"],
                        header=False,
                    )
                    if category != "maves":
                        category_summaries.append(dict(target))
                        all_measurements.append(target_scores)
                        continue
                    hgvs = target_scores["hgvs_pro"]
                    counts = hgvs.str.count(";") + 1
                    row = {
                        **target,
                        "n": len(target_scores),
                        "n_subs": int(counts.sum()),
                        "n_positions": hgvs.str.findall(r"[0-9]+").explode().nunique(),
                        "multiple_mutations": bool((counts > 1).any()),
                        "single_mutations": bool((counts == 1).any()),
                        "wt_mutants": bool(target_scores["wt_mut"].any()),
                    }
                    if row["protein_sequence"] is None:
                        row["protein_sequence"] = translate_dna(row["dna_sequence"])
                    singles = make_singles(target_scores, target)
                    wide = make_wide(singles, target)
                    write_tsv(singles, outputs[category, "singles"], header=False)
                    write_tsv(wide, outputs[category, "wide"], header=False)
                    if row["multiple_mutations"]:
                        write_tsv(
                            make_doubles(target_scores, target),
                            outputs[category, "doubles"],
                            header=False,
                        )
                    category_summaries.append(row)
                    all_singles.append(singles)
                    all_wide.append(wide)
                    all_measurements.append(target_scores)
                if not category_summaries:
                    continue
                combined = pd.concat(all_measurements, ignore_index=True)
                if category == "maves":
                    summary = score_summary(combined["score"])
                    sd = summary["sd_score"]
                    outliers = (
                        int(
                            (
                                (combined["score"] - summary["mean_score"]).abs()
                                > 2 * sd
                            ).sum()
                        )
                        if pd.notna(sd)
                        else None
                    )
                    summary.update(
                        n_outlier=outliers,
                        p_outlier=(
                            outliers / len(combined) if outliers is not None else None
                        ),
                        mean_variant_reps=pd.concat(all_singles)["n"].mean(),
                    )
                    summary.update(
                        position_summary(pd.concat(all_wide, ignore_index=True))
                    )
                    if any(row["multiple_mutations"] for row in category_summaries):
                        summary.update(multiple_summary(combined["hgvs_pro"]))
                else:
                    summary = nucleotide_summary(combined, category)
                summaries[category].extend(
                    {**row, **summary} for row in category_summaries
                )
            if index % 100 == 0:
                LOGGER.info("Processed %d/%d score files", index, len(score_files))
    for category, rows in summaries.items():
        columns = META_COLUMNS[category]
        # R only emits protein_sequence for splicing when such targets exist.
        if category == "splicing" and not any(
            row["sequence_type"] == "protein" for row in rows
        ):
            columns = [column for column in columns if column != "protein_sequence"]
        frame = pd.DataFrame(rows).reindex(columns=columns)
        write_tsv(frame, output_path / category / "meta.tsv", quote=True)
        description_columns = list(META_COLUMNS[category])
        if category == "maves":
            # The R descriptor table lists min_muts before max_muts, whereas
            # meta.tsv itself lists max_muts first.
            first, second = (
                description_columns.index(name) for name in ("min_muts", "max_muts")
            )
            description_columns[first], description_columns[second] = (
                description_columns[second],
                description_columns[first],
            )
        descriptions = [
            dict(
                column=column,
                type=COLUMN_DESCRIPTIONS[category][column][0],
                description=COLUMN_DESCRIPTIONS[category][column][1],
            )
            for column in description_columns
        ]
        write_tsv(
            pd.DataFrame(descriptions), output_path / category / "meta_columns.tsv"
        )
    return {category: len(rows) for category, rows in summaries.items()}


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
        help="Parent directory for maves/, mpras/ and splicing/.",
    )
    parser.add_argument(
        "--scores", type=Path, help="Override the dump's csv/ directory."
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        help="Use this existing NCBI taxcat.dmp instead of the cached/downloaded file.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        counts = process_mavedb(
            args.input,
            args.output,
            scores_path=args.scores,
            taxonomy_path=args.taxonomy,
        )
    except (OSError, ValueError, KeyError, pd.errors.ParserError) as error:
        parser.error(str(error))
    LOGGER.info("Wrote metadata rows %s to %s", counts, args.output)


if __name__ == "__main__":
    main()
