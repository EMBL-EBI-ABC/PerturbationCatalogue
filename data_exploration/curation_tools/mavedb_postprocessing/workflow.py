"""Composable file and dataframe steps for MaveDB post-processing."""

from __future__ import annotations

import json
import os
from contextlib import redirect_stdout
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as parquet

from curation_tools.curation_tools import concatenate_parquet_files
from curation_tools.mavedb_curation_tools import process_mavedb
from curation_tools.mavedb_postprocessing.standardization import (
    DEFAULT_MAPPING_PATH,
    DEFAULT_OVERRIDE_PATH,
    filter_curated_metadata,
    standardize_metadata,
    validate_with_unique_values,
)
from curation_tools.perturbseq_anndata_schema import ObsSchema

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MAVEDB_CSV_DIR = (
    REPO_ROOT / "data_exploration/MaveDB/Dump/mavedb-dump.20250612164404/csv"
)
EXCLUDED_INPUT_COLUMNS = ("perturbed_target_number", "__source_files", "__source_urns")


@dataclass(frozen=True)
class FinalArtifacts:
    """Paths to the checked metadata and data products from one curation run."""

    metadata_path: Path
    data_path: Path
    prevalidation_metadata_path: Path
    standardized_metadata_path: Path
    combined_data_path: Path
    joint_data_metadata_path: Path | None = None


def _validate_source(
    source: pd.DataFrame,
    source_name: str,
    *,
    allow_duplicate_dataset_ids: bool = False,
) -> pd.DataFrame:
    frame = source.copy()
    if "dataset_id" not in frame.columns:
        raise ValueError(f"{source_name} is missing the dataset_id column")
    frame = frame.drop(columns=list(EXCLUDED_INPUT_COLUMNS), errors="ignore")
    frame = frame.replace(r"^\s*$", pd.NA, regex=True)
    frame["dataset_id"] = frame["dataset_id"].astype("string").str.strip()
    missing_ids = frame["dataset_id"].isna() | frame["dataset_id"].eq("").fillna(False)
    if missing_ids.any():
        raise ValueError(f"{source_name} contains rows without dataset_id")
    duplicates = frame.loc[frame["dataset_id"].duplicated(), "dataset_id"].unique()
    if len(duplicates) and not allow_duplicate_dataset_ids:
        raise ValueError(
            f"{source_name} contains duplicate dataset_id values: {duplicates[:10].tolist()}"
        )
    return frame


def load_and_merge_metadata(
    llm_metadata_path: str | Path,
    manual_metadata_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load sources and preserve manual rows, with non-empty human values taking precedence."""
    llm = pd.read_csv(llm_metadata_path, dtype=str)
    llm["is_curated"] = "Yes"
    llm = _validate_source(llm, "LLM metadata").set_index("dataset_id", drop=True)

    if manual_metadata_path is None:
        empty_conflicts = pd.DataFrame(
            columns=[
                "dataset_id",
                "manual_row_number",
                "field",
                "manual_value",
                "llm_value",
            ]
        )
        return llm.reset_index(), empty_conflicts

    manual = pd.read_excel(manual_metadata_path, dtype=str)
    manual = _validate_source(
        manual, "Manual metadata", allow_duplicate_dataset_ids=True
    )
    columns = ["dataset_id"] + [
        column
        for column in list(manual.columns) + list(llm.columns)
        if column != "dataset_id"
    ]
    columns = list(dict.fromkeys(columns))

    manual_aligned = manual.reindex(columns=columns).reset_index(drop=True)
    llm_for_manual = llm.reindex(index=manual["dataset_id"])
    llm_aligned = llm_for_manual.reset_index(drop=True).reindex(columns=columns)

    conflicts: list[dict[str, str]] = []
    shared_columns = manual.columns.intersection(llm.columns)
    for field in shared_columns:
        if field == "dataset_id":
            continue
        manual_values = manual_aligned[field].astype("string")
        llm_values = llm_aligned[field].astype("string")
        both_present = manual_values.notna() & llm_values.notna()
        differs = manual_values.ne(llm_values).fillna(False)
        for row_position in (both_present & differs).to_numpy().nonzero()[0]:
            conflicts.append(
                {
                    "dataset_id": str(manual_aligned.at[row_position, "dataset_id"]),
                    # Excel row 1 contains the column headers.
                    "manual_row_number": str(row_position + 2),
                    "field": str(field),
                    "manual_value": str(manual_aligned.at[row_position, field]),
                    "llm_value": str(llm_aligned.at[row_position, field]),
                }
            )

    merged = manual_aligned.combine_first(llm_aligned)
    merged = merged.reindex(columns=columns)
    manual_ids = set(manual["dataset_id"].tolist())
    llm_only = llm.loc[~llm.index.isin(manual_ids)].reset_index()
    if not llm_only.empty:
        merged = pd.concat(
            [merged, llm_only.reindex(columns=columns)], ignore_index=True
        )
    conflict_table = pd.DataFrame(
        conflicts,
        columns=[
            "dataset_id",
            "manual_row_number",
            "field",
            "manual_value",
            "llm_value",
        ],
    )
    return merged.reset_index(drop=True), conflict_table


def create_run_directory(output_dir: str | Path) -> Path:
    """Create an empty run directory and refuse to mix runs in an existing folder."""
    run_dir = Path(output_dir).expanduser().resolve()
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(
            f"Output directory is not empty: {run_dir}. Choose a new run directory."
        )
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def process_mavedb_datasets(
    metadata: pd.DataFrame,
    mavedb_csv_dir: str | Path,
    output_dir: str | Path,
    progress_interval: int = 25,
) -> Path:
    """Create per-dataset AnnData and parquet outputs; write verbose logs to disk."""
    if "dataset_id" not in metadata.columns:
        raise ValueError("Metadata is missing the dataset_id column")
    if progress_interval < 1:
        raise ValueError("progress_interval must be positive")

    run_dir = Path(output_dir).expanduser().resolve()
    non_curated_h5ad_dir = run_dir / "non_curated" / "h5ad"
    non_curated_h5ad_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "logs" / "per_dataset_processing.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    dataset_ids = metadata["dataset_id"].dropna().astype(str).drop_duplicates().tolist()
    if not dataset_ids:
        raise ValueError("No dataset IDs were selected for processing")

    with log_path.open("w", encoding="utf-8") as log_handle:
        for index, dataset_id in enumerate(dataset_ids, start=1):
            with redirect_stdout(log_handle):
                process_mavedb(
                    mavedb_dataset_id=dataset_id,
                    mavedb_csv_dir=str(mavedb_csv_dir),
                    curated_metadata_df=metadata,
                    non_curated_h5ad_dir=str(non_curated_h5ad_dir),
                    overwrite=False,
                )
            if index % progress_interval == 0 or index == len(dataset_ids):
                print(f"Processed {index:,} of {len(dataset_ids):,} datasets")

    return log_path


def combine_dataset_parquets(output_dir: str | Path) -> tuple[Path, Path]:
    """Stream the per-dataset parquet files into combined metadata and data files."""
    run_dir = Path(output_dir).expanduser().resolve()
    parquet_dir = run_dir / "curated" / "parquet"
    combined_dir = run_dir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = combined_dir / "mavedb_all_curated_metadata.parquet"
    data_path = combined_dir / "mavedb_all_curated_data.parquet"
    concatenate_parquet_files(
        parquet_dir=str(parquet_dir),
        output_path=str(metadata_path),
        pattern="*_curated_metadata.parquet",
        verbose=False,
    )
    concatenate_parquet_files(
        parquet_dir=str(parquet_dir),
        output_path=str(data_path),
        pattern="*_curated_data.parquet",
        verbose=False,
    )
    return metadata_path, data_path


def finalize_outputs(
    standardized_metadata_path: str | Path,
    combined_data_path: str | Path,
    output_dir: str | Path,
    save_joint_artifact: bool = False,
) -> FinalArtifacts:
    """Filter, validate, and save final metadata and matching score data.

    When ``save_joint_artifact`` is true, also save score data joined to its
    validated metadata on ``dataset_id`` and, where needed, target symbol.
    """
    final_dir = Path(output_dir).expanduser().resolve() / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_parquet(standardized_metadata_path)
    keep_metadata = metadata["curation_agent_type"].notna()
    metadata = metadata.loc[keep_metadata].copy()

    prevalidation_path = final_dir / "mavedb_all_curated_metadata_prevalidation.parquet"
    metadata.to_parquet(prevalidation_path, index=False)
    validated_metadata = validate_with_unique_values(
        ObsSchema, metadata, n_failure_cases=3
    )
    metadata_path = final_dir / "mavedb_all_curated_metadata_postvalidation.parquet"
    validated_metadata.to_parquet(metadata_path, index=False)

    data = pd.read_parquet(combined_data_path)
    if "dataset_id" not in data.columns:
        raise ValueError("Combined MaveDB data is missing the dataset_id column")
    data = data.loc[data["dataset_id"].isin(validated_metadata["dataset_id"])].copy()
    data_path = final_dir / "mavedb_all_curated_data_postfilter.parquet"
    data.to_parquet(data_path, index=False)

    joint_data_metadata_path = None
    if save_joint_artifact:
        join_columns = ["dataset_id", "sample_id"]
        metadata_columns = [
            column
            for column in validated_metadata.columns
            if column in join_columns or column not in data.columns
        ]
        joint = data.merge(
            validated_metadata[metadata_columns],
            on=join_columns,
            how="left",
            validate="many_to_one",
        )
        joint_data_metadata_path = (
            final_dir / "mavedb_all_curated_data_with_metadata_postfilter.parquet"
        )
        joint.to_parquet(joint_data_metadata_path, index=False)

    return FinalArtifacts(
        metadata_path=metadata_path,
        data_path=data_path,
        prevalidation_metadata_path=prevalidation_path,
        standardized_metadata_path=Path(standardized_metadata_path),
        combined_data_path=Path(combined_data_path),
        joint_data_metadata_path=joint_data_metadata_path,
    )


def run_pipeline(
    llm_metadata_path: str | Path,
    output_dir: str | Path,
    manual_metadata_path: str | Path | None = None,
    mavedb_csv_dir: str | Path = DEFAULT_MAVEDB_CSV_DIR,
    mapping_path: str | Path = DEFAULT_MAPPING_PATH,
    override_path: str | Path = DEFAULT_OVERRIDE_PATH,
    save_joint_artifact: bool = False,
) -> dict[str, Any]:
    """Run the local MaveDB metadata and score-data workflow in a fresh directory."""
    run_dir = create_run_directory(output_dir)
    source_dir = run_dir / "sources"
    normalized_dir = run_dir / "normalized"
    qc_dir = run_dir / "qc"
    for directory in (source_dir, normalized_dir, qc_dir):
        directory.mkdir(parents=True, exist_ok=True)

    merged, conflicts = load_and_merge_metadata(llm_metadata_path, manual_metadata_path)
    merged.to_csv(source_dir / "merged_metadata.csv", index=False)
    conflicts.to_csv(qc_dir / "source_conflicts.csv", index=False)

    curated = filter_curated_metadata(merged)
    curated.to_csv(source_dir / "metadata_for_processing.csv", index=False)

    processing_log = process_mavedb_datasets(
        curated, mavedb_csv_dir=mavedb_csv_dir, output_dir=run_dir
    )
    combined_metadata_path, combined_data_path = combine_dataset_parquets(run_dir)

    combined_metadata = pd.read_parquet(combined_metadata_path)
    standardized, metadata_qc, override_audit = standardize_metadata(
        combined_metadata,
        mapping_path=mapping_path,
        override_path=override_path,
        inplace=True,
    )
    standardized_metadata_path = (
        normalized_dir / "mavedb_all_curated_metadata_standardized.parquet"
    )
    standardized.to_parquet(standardized_metadata_path, index=False)
    missing_ids = metadata_qc.loc[metadata_qc["issue_type"].eq("missing_expected_id")]
    unmapped = metadata_qc.loc[~metadata_qc["issue_type"].eq("missing_expected_id")]
    unmapped.to_csv(qc_dir / "unmapped_terms.csv", index=False)
    missing_ids.to_csv(qc_dir / "missing_expected_ids.csv", index=False)
    override_audit.to_csv(qc_dir / "override_audit.csv", index=False)

    artifacts = finalize_outputs(
        standardized_metadata_path=standardized_metadata_path,
        combined_data_path=combined_data_path,
        output_dir=run_dir,
        save_joint_artifact=save_joint_artifact,
    )

    manifest = {
        "llm_metadata_path": str(Path(llm_metadata_path).expanduser().resolve()),
        "manual_metadata_path": (
            str(Path(manual_metadata_path).expanduser().resolve())
            if manual_metadata_path is not None
            else None
        ),
        "mavedb_csv_dir": str(Path(mavedb_csv_dir).expanduser().resolve()),
        "mapping_path": str(Path(mapping_path).expanduser().resolve()),
        "override_path": str(Path(override_path).expanduser().resolve()),
        "dataset_count_merged": int(merged["dataset_id"].nunique()),
        "dataset_count_processed": int(curated["dataset_id"].nunique()),
        "source_conflict_count": len(conflicts),
        "unmapped_term_rows": len(unmapped),
        "missing_expected_id_rows": len(missing_ids),
        "raw_combined_metadata_path": str(combined_metadata_path),
        "standardized_metadata_path": str(standardized_metadata_path),
        "processing_log": str(processing_log),
        "artifacts": {
            key: str(value) if value is not None else None
            for key, value in asdict(artifacts).items()
        },
    }
    (run_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def publish_to_bigquery(
    metadata_path: str | Path,
    data_path: str | Path,
    project_id: str | None = None,
) -> None:
    """Explicitly merge the final parquet files into the existing MaveDB tables."""
    from google.cloud import bigquery

    from curation_tools.curation_tools import _upload_parquet_to_bq

    project = project_id or os.environ.get("BQ_PROJECT")
    if not project:
        raise ValueError(
            "project_id must be provided or set in the BQ_PROJECT environment variable"
        )
    client = bigquery.Client(project=project)
    table_inputs = (
        (metadata_path, "metadata"),
        (data_path, "data"),
    )
    for parquet_path, table_name in table_inputs:
        source_columns = {
            column.lower()
            for column in parquet.ParquetFile(parquet_path).schema_arrow.names
        }
        table = client.get_table(f"{project}.mavedb.{table_name}")
        target_columns = {field.name.lower() for field in table.schema}
        missing_from_table = source_columns - target_columns - {"ingested_at"}
        missing_from_file = target_columns - source_columns - {"row_id", "ingested_at"}
        if missing_from_table or missing_from_file:
            raise ValueError(
                f"BigQuery schema mismatch for {table_name}: "
                f"columns missing from table={sorted(missing_from_table)}, "
                f"columns missing from parquet={sorted(missing_from_file)}. "
                "Review the MaveDB metadata schema migration before publishing."
            )

    _upload_parquet_to_bq(
        parquet_path=str(metadata_path),
        project_id=project,
        bq_dataset_id="mavedb",
        bq_table_name="metadata",
        key_columns=["dataset_id", "sample_id"],
    )
    _upload_parquet_to_bq(
        parquet_path=str(data_path),
        project_id=project,
        bq_dataset_id="mavedb",
        bq_table_name="data",
        key_columns=["dataset_id", "sample_id", "score_name"],
    )
