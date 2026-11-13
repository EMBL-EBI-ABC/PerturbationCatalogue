"""Streaming readers for the MAVE summary-plot input tables."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

import pandas as pd
from pandas.api.types import is_bool_dtype, is_float_dtype, is_integer_dtype

TABLE_NAMES = (
    "scores",
    "singles",
    "density",
    "posterior",
    "fit_summary",
    "params",
    "loo",
)
_READR_NA_VALUES = ["", "NA"]
_STREAM_CHUNK_ROWS = 100_000
_INTEGER_COLUMNS = {"position", "iteration", "chain", "draw", "n", "n_positions"}
_TEXT_COLUMNS = {
    "dna_sequence",
    "hgvs_pro",
    "hgvs_nt",
    "identifier_uniprot",
    "mode",
    "mode_density",
    "mode_sd",
    "model",
    "mut",
    "name",
    "organism",
    "target_name",
    "variable",
    "wt",
}
_NATIVE_PREFIXES = {
    "scores": "score",
    "density": "density",
    "params": "fitted_params",
    "posterior": "posterior",
    "fit_summary": "fit_summary",
    "loo": "loo",
}
_NATIVE_REQUIRED_COLUMNS = {
    "scores": {"score_urn", "model", "top_model", "hgvs_pro", "s", "w"},
    "density": {
        "score_urn",
        "model",
        "top_model",
        "mode",
        "name",
        "x",
        "density",
    },
    "params": {"score_urn", "model", "top_model", "median", "mad"},
    "posterior": {
        "score_urn",
        "model",
        "top_model",
        "iteration",
        "chain",
        "draw",
    },
    "fit_summary": {"score_urn", "model", "top_model", "variable", "rhat"},
    "loo": {"score_urn", "model", "elpd_loo", "se_elpd_loo"},
}


class TableDataError(ValueError):
    """Raised when an input table is missing required structure or has bad values."""


def _table_paths(data_root: Path, model: str) -> dict[str, Path]:
    gmm_root = data_root / "maves" / f"gmm_{model}"
    return {
        "meta": data_root / "maves" / "meta.tsv",
        **{name: gmm_root / f"{name}.tsv" for name in TABLE_NAMES},
    }


def _read_header(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input table not found: {path}")
    try:
        return list(
            pd.read_csv(
                path,
                sep="\t",
                nrows=0,
                keep_default_na=False,
                na_values=_READR_NA_VALUES,
            ).columns
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise TableDataError(f"Could not read TSV header from {path}: {exc}") from exc


def _infer_dtype_map(path: Path, columns: list[str]) -> dict[str, str]:
    """Infer a stable chunk dtype map from the same first 1,000 rows as readr."""
    try:
        sample = pd.read_csv(
            path,
            sep="\t",
            nrows=1_000,
            keep_default_na=False,
            na_values=_READR_NA_VALUES,
            dtype={"score_urn": "string"} if "score_urn" in columns else None,
            low_memory=False,
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise TableDataError(f"Could not infer TSV columns from {path}: {exc}") from exc

    dtype_map: dict[str, str] = {}
    for column in columns:
        if column == "score_urn":
            dtype_map[column] = "string"
            continue
        if column == "top_model" or column in _TEXT_COLUMNS:
            dtype_map[column] = "object"
            continue
        if column not in sample.columns or not sample[column].notna().any():
            dtype_map[column] = "object"
            continue
        dtype = sample[column].dtype if column in sample.columns else object
        if is_bool_dtype(dtype):
            dtype_map[column] = "boolean"
        elif is_integer_dtype(dtype):
            dtype_map[column] = "Int64" if column in _INTEGER_COLUMNS else "float64"
        elif is_float_dtype(dtype):
            dtype_map[column] = "float64"
        else:
            dtype_map[column] = "object"
    return dtype_map


def _normalise_top_model(frame: pd.DataFrame, path: Path) -> pd.DataFrame:
    if "top_model" not in frame.columns:
        return frame

    normalized = frame["top_model"].astype("string").str.strip().str.lower()
    values = normalized.map(
        {
            "true": True,
            "t": True,
            "1": True,
            "yes": True,
            "false": False,
            "f": False,
            "0": False,
            "no": False,
        }
    )
    invalid = normalized.notna() & values.isna()
    if invalid.any():
        examples = sorted(normalized[invalid].dropna().unique().tolist())[:5]
        raise TableDataError(f"Invalid top_model boolean value in {path}: {examples}")

    frame["top_model"] = values.astype("boolean")
    return frame


def _read_tsv(
    path: Path,
    dtype_map: dict[str, str],
) -> pd.DataFrame:
    try:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=dtype_map,
            keep_default_na=False,
            na_values=_READR_NA_VALUES,
            low_memory=False,
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise TableDataError(f"Could not parse TSV table {path}: {exc}") from exc


def _iter_tsv_chunks(
    path: Path,
    dtype_map: dict[str, str],
    *,
    usecols: list[str] | None = None,
) -> Iterator[pd.DataFrame]:
    try:
        chunks = pd.read_csv(
            path,
            sep="\t",
            dtype=dtype_map,
            usecols=usecols,
            keep_default_na=False,
            na_values=_READR_NA_VALUES,
            chunksize=_STREAM_CHUNK_ROWS,
        )
        yield from chunks
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise TableDataError(f"Could not parse TSV table {path}: {exc}") from exc


def _empty_table(path: Path, dtype_map: dict[str, str]) -> pd.DataFrame:
    try:
        frame = pd.read_csv(
            path,
            sep="\t",
            nrows=0,
            dtype=dtype_map,
            keep_default_na=False,
            na_values=_READR_NA_VALUES,
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise TableDataError(f"Could not read empty TSV table {path}: {exc}") from exc
    return _normalise_top_model(frame, path)


def _filtered_table(
    path: Path, dtype_map: dict[str, str], urns: set[str]
) -> pd.DataFrame:
    if not urns:
        return _empty_table(path, dtype_map)

    parts: list[pd.DataFrame] = []
    for chunk in _iter_tsv_chunks(path, dtype_map):
        chunk = _normalise_top_model(chunk, path)
        selected = chunk.loc[chunk["score_urn"].isin(urns)]
        if not selected.empty:
            parts.append(selected)

    if not parts:
        return _empty_table(path, dtype_map)
    return pd.concat(parts, ignore_index=True, copy=False)


def _nonmissing_values(series: pd.Series) -> set[str]:
    return {str(value) for value in series.dropna().tolist()}


def _ordered_unique_urns(path: Path) -> list[str]:
    seen: set[str] = set()
    urns: list[str] = []
    for chunk in _iter_tsv_chunks(path, {"score_urn": "string"}, usecols=["score_urn"]):
        for value in chunk["score_urn"].dropna().tolist():
            urn = str(value)
            if urn not in seen:
                seen.add(urn)
                urns.append(urn)
    return urns


def load_tables(
    data_root: Path,
    model: str,
    index: int = 0,
    batch_size: int = 100,
) -> tuple[dict[str, pd.DataFrame], list[str]]:
    """Load one selected URN batch and the tables needed to render its plots.

    Scores retain all rows for the selected score URNs and for any score sets
    sharing their target name or UniProt identifier, as required for the
    correlation panel. Other GMM tables retain only rows for selected URNs.
    """
    if model not in {"modes", "flexible"}:
        raise ValueError("model must be one of: modes, flexible")
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError("index must be a nonnegative integer")
    if (
        isinstance(batch_size, bool)
        or not isinstance(batch_size, int)
        or batch_size < 1
    ):
        raise ValueError("batch_size must be a positive integer")

    root = Path(data_root)
    paths = _table_paths(root, model)
    headers = {name: _read_header(path) for name, path in paths.items()}

    required_columns = {
        "meta": {
            "score_urn",
            "target_name",
            "identifier_uniprot",
            "organism",
            "n",
            "n_positions",
            "multiple_mutations",
        },
        "scores": {"score_urn", "hgvs_pro", "s", "w"},
        "singles": {"score_urn", "position", "wt", "mut", "s", "w"},
        "density": {"score_urn", "name", "mode", "x", "density"},
        "posterior": {"score_urn", "iteration", "chain", "draw"},
        "fit_summary": {"score_urn", "top_model", "variable", "rhat"},
        "params": {"score_urn", "top_model", "median", "mad"},
        "loo": {"score_urn", "model", "elpd_loo", "se_elpd_loo"},
    }
    for name, required in required_columns.items():
        missing = sorted(required - set(headers[name]))
        if missing:
            raise TableDataError(
                f"Input table {paths[name]} is missing required column(s): "
                f"{', '.join(missing)}"
            )

    schemas = {name: _infer_dtype_map(paths[name], headers[name]) for name in paths}
    meta = _read_tsv(paths["meta"], schemas["meta"])
    assert isinstance(meta, pd.DataFrame)

    all_urns = _ordered_unique_urns(paths["scores"])
    if index == 0:
        selected_urns = all_urns
    else:
        start = (index - 1) * batch_size
        selected_urns = all_urns[start : start + batch_size]
    selected_set = set(selected_urns)

    related_urns = set(selected_set)
    if selected_set:
        selected_meta = meta.loc[meta["score_urn"].isin(selected_set)]
        target_names = _nonmissing_values(selected_meta["target_name"])
        uniprot_ids = _nonmissing_values(selected_meta["identifier_uniprot"])
        related_meta_mask = meta["target_name"].isin(target_names) | meta[
            "identifier_uniprot"
        ].isin(uniprot_ids)
        related_urns.update(
            _nonmissing_values(meta.loc[related_meta_mask, "score_urn"])
        )

    tables: dict[str, pd.DataFrame] = {"meta": meta}
    tables["scores"] = _filtered_table(paths["scores"], schemas["scores"], related_urns)
    for name in TABLE_NAMES[1:]:
        tables[name] = _filtered_table(paths[name], schemas[name], selected_set)

    return tables, selected_urns


def _native_file_paths(gmm_dir: Path) -> tuple[list[int], dict[str, dict[int, Path]]]:
    if not gmm_dir.is_dir():
        raise FileNotFoundError(f"GMM output directory not found: {gmm_dir}")

    def validate_paths(
        paths: dict[str, dict[int, Path]], indexes: list[int]
    ) -> tuple[list[int], dict[str, dict[int, Path]]]:
        expected = set(indexes)
        for table_name, table_paths in paths.items():
            actual = set(table_paths)
            if actual != expected:
                missing = sorted(expected - actual)
                extra = sorted(actual - expected)
                details = []
                if missing:
                    details.append(f"missing indices {missing}")
                if extra:
                    details.append(f"unexpected indices {extra}")
                raise TableDataError(
                    f"Indexed GMM table family {table_name!r} in {gmm_dir} has "
                    + "; ".join(details)
                )
        return indexes, paths

    paths: dict[str, dict[int, Path]] = {name: {} for name in _NATIVE_PREFIXES}
    batch_root = gmm_dir / "batches"
    if not batch_root.is_dir() and any(
        child.is_dir() and re.fullmatch(r"batch\d+", child.name)
        for child in gmm_dir.iterdir()
    ):
        batch_root = gmm_dir

    if batch_root.is_dir():
        batch_directories = []
        for child in batch_root.iterdir():
            match = re.fullmatch(r"batch(\d+)", child.name)
            if match and child.is_dir():
                batch_directories.append((int(match.group(1)), child))
        for index, batch_directory in batch_directories:
            for table_name, prefix in _NATIVE_PREFIXES.items():
                path = batch_directory / f"{prefix}.tsv"
                if path.is_file():
                    paths[table_name][index] = path
        indexes = sorted(paths["scores"])
        if indexes:
            return validate_paths(paths, indexes)

    for merged_root in (gmm_dir, gmm_dir / "merged"):
        merged_paths = {
            name: merged_root / f"{prefix}.tsv"
            for name, prefix in _NATIVE_PREFIXES.items()
        }
        if all(path.is_file() for path in merged_paths.values()):
            return [1], {name: {1: path} for name, path in merged_paths.items()}
        if any(path.is_file() for path in merged_paths.values()):
            missing = sorted(
                path.name for path in merged_paths.values() if not path.is_file()
            )
            raise TableDataError(
                f"Merged GMM table family in {merged_root} is missing: "
                f"{', '.join(missing)}"
            )

    for table_name, prefix in _NATIVE_PREFIXES.items():
        pattern = re.compile(rf"^{re.escape(prefix)}(\d+)\.tsv$")
        for path in gmm_dir.glob(f"{prefix}*.tsv"):
            match = pattern.fullmatch(path.name)
            if match:
                paths[table_name][int(match.group(1))] = path

    indexes = sorted(paths["scores"])
    if not indexes:
        raise FileNotFoundError(
            f"No GMM batches, merged tables, or indexed score<N>.tsv exports "
            f"found in {gmm_dir}"
        )
    return validate_paths(paths, indexes)


def _load_canonical_meta(meta_path: Path) -> pd.DataFrame:
    if not meta_path.is_file():
        raise FileNotFoundError(f"Canonical metadata file not found: {meta_path}")
    from data_sources.mavedb.analysis._mave_qc import (
        _global_text,
        _is_multi_json,
        _scalar_json_text,
    )

    try:
        frame = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as exc:
        raise TableDataError(
            f"Could not read canonical metadata {meta_path}: {exc}"
        ) from exc

    if "dataset_id" not in frame.columns:
        raise TableDataError(
            f"Canonical metadata {meta_path} is missing required column 'dataset_id'"
        )

    if "primary_analysis" in frame.columns:
        primary = frame["primary_analysis"].map(_global_text).str.lower()
        frame = frame.loc[primary.eq("maves")].copy()

    if "target_name" in frame.columns:
        target_is_multi = frame["target_name"].map(_is_multi_json)
    else:
        target_is_multi = pd.Series(False, index=frame.index)

    dataset_ids = frame["dataset_id"].map(_global_text)
    duplicated = dataset_ids.notna() & dataset_ids.duplicated(keep=False)
    if duplicated.any():
        examples = sorted(dataset_ids[duplicated].dropna().unique().tolist())[:5]
        raise TableDataError(
            f"Canonical metadata {meta_path} contains duplicate dataset_id values: "
            + ", ".join(examples)
        )

    frame["score_urn"] = dataset_ids
    for column in ("target_name", "identifier_uniprot", "organism"):
        if column not in frame.columns:
            frame[column] = pd.NA
        frame[column] = (
            frame[column]
            .map(
                lambda value: _scalar_json_text(
                    value, preserve_literal_na=(column == "target_name")
                )
            )
            .astype("string")
        )

    if "n_analyzed_variants" in frame.columns:
        frame["n"] = pd.to_numeric(frame["n_analyzed_variants"], errors="coerce")
    else:
        frame["n"] = pd.NA
    if "n_positions" in frame.columns:
        frame["n_positions"] = pd.to_numeric(frame["n_positions"], errors="coerce")
    else:
        frame["n_positions"] = pd.NA
    for column in ("multiple_mutations", "single_mutations"):
        if column not in frame.columns:
            frame[column] = pd.NA
        else:
            values = frame[column].astype("string").str.strip().str.lower()
            frame[column] = values.map(
                {
                    "true": True,
                    "t": True,
                    "1": True,
                    "yes": True,
                    "false": False,
                    "f": False,
                    "0": False,
                    "no": False,
                }
            ).astype("boolean")

    if "offset_uniprot" in frame.columns:
        decoded_offset = frame["offset_uniprot"].map(_scalar_json_text)
        frame["offset_uniprot"] = pd.to_numeric(decoded_offset, errors="coerce")
    else:
        frame["offset_uniprot"] = float("nan")

    frame["target_is_multi"] = target_is_multi.astype(bool)
    ambiguous = frame["target_is_multi"]
    frame.loc[ambiguous, ["target_name", "identifier_uniprot"]] = pd.NA
    frame.loc[ambiguous, "offset_uniprot"] = float("nan")
    return frame


def _native_ordered_urns(
    indexes: list[int], paths: dict[str, dict[int, Path]]
) -> list[str]:
    ordered: list[str] = []
    owners: dict[str, int] = {}
    for batch_index in indexes:
        scores = _ordered_unique_urns(paths["scores"][batch_index])
        params = _ordered_unique_urns(paths["params"][batch_index])
        batch_urns = set(scores) | set(params)
        duplicates = sorted(batch_urns & owners.keys())
        if duplicates:
            examples = ", ".join(duplicates[:5])
            prior = owners[duplicates[0]]
            raise TableDataError(
                f"Study score_urn(s) {examples} appear in multiple GMM batch files "
                f"(indices {prior} and {batch_index})"
            )
        for urn in batch_urns:
            owners[urn] = batch_index
        ordered.extend(scores)
        ordered.extend(urn for urn in params if urn not in set(scores))
    return ordered


def _native_top_model_rows(
    frame: pd.DataFrame, path: Path, *, keep_failed_scores: bool
) -> pd.DataFrame:
    frame = _normalise_top_model(frame, path)
    top = frame["top_model"].fillna(False)
    if keep_failed_scores:
        failed_raw = frame["top_model"].isna() & frame["model"].isna()
        return frame.loc[top | failed_raw]
    return frame.loc[top]


def _filtered_native_table(
    file_paths: list[Path],
    dtype_maps: list[dict[str, str]],
    urns: set[str],
    *,
    top_model_only: bool,
    keep_failed_scores: bool = False,
    model: str,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    valid_models = {"flexible"} if model == "flexible" else {"ko", "del", "neut", "pos"}
    for path, dtype_map in zip(file_paths, dtype_maps, strict=True):
        if not urns:
            continue
        for chunk in _iter_tsv_chunks(path, dtype_map):
            if top_model_only:
                chunk = _native_top_model_rows(
                    chunk, path, keep_failed_scores=keep_failed_scores
                )
            selected = chunk.loc[chunk["score_urn"].isin(urns)]
            if "model" in selected.columns:
                found_models = set(selected["model"].dropna().astype(str))
                unexpected = sorted(found_models - valid_models)
                if unexpected:
                    raise TableDataError(
                        f"GMM table {path} contains candidate model(s) {unexpected} "
                        f"incompatible with --model {model!r}"
                    )
            if not selected.empty:
                parts.append(selected)

    if parts:
        return pd.concat(parts, ignore_index=True, sort=False, copy=False)
    return _empty_table(file_paths[0], dtype_maps[0])


def _derive_native_singles(
    scores: pd.DataFrame,
    meta: pd.DataFrame,
    selected_urns: set[str],
) -> pd.DataFrame:
    columns = ["score_urn", "position", "wt", "mut", "n", "s", "w"]
    if not selected_urns or scores.empty:
        return pd.DataFrame(columns=columns)

    from data_sources.mavedb.scripts.process_mavedb import AA_3_TO_1, SINGLE_PATTERN

    selected_scores = scores.loc[scores["score_urn"].isin(selected_urns)]
    offsets = meta.set_index("score_urn")["offset_uniprot"].to_dict()
    parsed_rows: list[pd.DataFrame] = []
    for score_urn, study_scores in selected_scores.groupby("score_urn", sort=False):
        offset = offsets.get(score_urn)
        if pd.isna(offset):
            # Keep the source HGVS coordinate when no unambiguous UniProt offset exists.
            offset = 0
        hgvs = study_scores["hgvs_pro"].astype("string")
        parsed = hgvs.str.extract(SINGLE_PATTERN)
        parsed.columns = ["wt", "position", "mut"]
        parsed["mut"] = parsed["mut"].mask(parsed["mut"].eq("="), parsed["wt"])
        parsed["wt"] = parsed["wt"].str.title().map(AA_3_TO_1)
        normalized_mut = parsed["mut"].str.title()
        normalized_mut = normalized_mut.mask(parsed["mut"].str.lower().eq("del"), "del")
        parsed["mut"] = normalized_mut.map(AA_3_TO_1)
        parsed["position"] = pd.to_numeric(parsed["position"], errors="coerce") + float(
            offset
        )
        parsed["score_urn"] = score_urn
        parsed["s"] = pd.to_numeric(study_scores["s"], errors="coerce")
        parsed["w"] = pd.to_numeric(study_scores["w"], errors="coerce")
        parsed = parsed.dropna(subset=["wt", "mut", "position"])
        if not parsed.empty:
            parsed_rows.append(parsed)

    if not parsed_rows:
        return pd.DataFrame(columns=columns)
    variants = pd.concat(parsed_rows, ignore_index=True)
    return (
        variants.groupby(
            ["score_urn", "wt", "position", "mut"],
            as_index=False,
            sort=False,
        )
        .agg(n=("s", "size"), s=("s", "mean"), w=("w", "mean"))
        .reindex(columns=columns)
    )


def load_native_tables(
    gmm_dir: Path,
    meta_path: Path,
    model: str,
    index: int = 0,
    batch_size: int = 100,
) -> tuple[dict[str, pd.DataFrame], list[str]]:
    """Load native GMM outputs for one report batch.

    ``index`` selects a non-overlapping report batch of unique studies across
    every fitter export in ``gmm_dir``. It does not select one fitter batch file.
    The directory may contain collected batch folders, merged tables, or legacy
    indexed fitter outputs.
    """
    if model not in {"modes", "flexible"}:
        raise ValueError("model must be one of: modes, flexible")
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError("index must be a nonnegative integer")
    if (
        isinstance(batch_size, bool)
        or not isinstance(batch_size, int)
        or batch_size < 1
    ):
        raise ValueError("batch_size must be a positive integer")

    gmm_dir = Path(gmm_dir)
    meta_path = Path(meta_path)
    indexes, paths_by_table = _native_file_paths(gmm_dir)
    meta = _load_canonical_meta(meta_path)

    dtype_maps: dict[str, dict[int, dict[str, str]]] = {}
    for table_name, indexed_paths in paths_by_table.items():
        dtype_maps[table_name] = {}
        for batch_index in indexes:
            path = indexed_paths[batch_index]
            header = _read_header(path)
            missing = sorted(_NATIVE_REQUIRED_COLUMNS[table_name] - set(header))
            if missing:
                raise TableDataError(
                    f"Input table {path} is missing required column(s): "
                    + ", ".join(missing)
                )
            dtype_maps[table_name][batch_index] = _infer_dtype_map(path, header)

    all_urns = _native_ordered_urns(indexes, paths_by_table)
    if index == 0:
        selected_urns = all_urns
    else:
        start = (index - 1) * batch_size
        selected_urns = all_urns[start : start + batch_size]
    selected_set = set(selected_urns)

    related_urns = set(selected_set)
    if selected_set:
        selected_meta = meta.loc[meta["score_urn"].isin(selected_set)]
        target_names = _nonmissing_values(selected_meta["target_name"])
        uniprot_ids = _nonmissing_values(selected_meta["identifier_uniprot"])
        related_meta_mask = meta["target_name"].isin(target_names) | meta[
            "identifier_uniprot"
        ].isin(uniprot_ids)
        related_urns.update(
            _nonmissing_values(meta.loc[related_meta_mask, "score_urn"])
        )

    batch_paths = lambda table: [paths_by_table[table][item] for item in indexes]
    batch_schemas = lambda table: [dtype_maps[table][item] for item in indexes]
    tables: dict[str, pd.DataFrame] = {"meta": meta}
    tables["scores"] = _filtered_native_table(
        batch_paths("scores"),
        batch_schemas("scores"),
        related_urns,
        top_model_only=True,
        keep_failed_scores=True,
        model=model,
    )
    tables["density"] = _filtered_native_table(
        batch_paths("density"),
        batch_schemas("density"),
        selected_set,
        top_model_only=True,
        model=model,
    )
    tables["posterior"] = _filtered_native_table(
        batch_paths("posterior"),
        batch_schemas("posterior"),
        selected_set,
        top_model_only=True,
        model=model,
    )
    tables["fit_summary"] = _filtered_native_table(
        batch_paths("fit_summary"),
        batch_schemas("fit_summary"),
        selected_set,
        top_model_only=True,
        model=model,
    )
    tables["params"] = _filtered_native_table(
        batch_paths("params"),
        batch_schemas("params"),
        selected_set,
        top_model_only=True,
        model=model,
    )
    tables["loo"] = _filtered_native_table(
        batch_paths("loo"),
        batch_schemas("loo"),
        selected_set,
        top_model_only=False,
        model=model,
    )
    tables["singles"] = _derive_native_singles(tables["scores"], meta, selected_set)
    return tables, selected_urns
