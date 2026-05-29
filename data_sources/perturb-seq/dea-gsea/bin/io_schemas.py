from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


DEA_SCHEMA = pa.schema(
    [
        pa.field("dataset_id", pa.string(), nullable=False),
        pa.field("perturbed_target_symbol", pa.string(), nullable=False),
        pa.field("gene", pa.string(), nullable=False),
        pa.field("padj", pa.float64(), nullable=False),
        pa.field("log2FoldChange", pa.float64(), nullable=False),
        pa.field("score_name", pa.string(), nullable=False),
        pa.field("score_value", pa.float64(), nullable=False),
        pa.field("cell_type", pa.string(), nullable=True),
        pa.field("ingested_at", pa.timestamp("us", tz="UTC"), nullable=True),
    ]
)

GSEA_SCHEMA = pa.schema(
    [
        pa.field("dataset_id", pa.string(), nullable=False),
        pa.field("term", pa.string(), nullable=False),
        pa.field("perturbed_target_symbol", pa.string(), nullable=False),
        pa.field("es", pa.float64(), nullable=False),
        pa.field("nes", pa.float64(), nullable=False),
        pa.field("pval", pa.float64(), nullable=False),
        pa.field("sidak", pa.float64(), nullable=False),
        pa.field("fdr", pa.float64(), nullable=False),
        pa.field("geneset_size", pa.int64(), nullable=False),
        pa.field("leading_edge", pa.list_(pa.string()), nullable=True),
        pa.field("cell_type", pa.string(), nullable=True),
        pa.field("ingested_at", pa.timestamp("us", tz="UTC"), nullable=True),
    ]
)


def _coerce_timestamp_values(values: pd.Series) -> list[pd.Timestamp | None]:
    out: list[pd.Timestamp | None] = []
    for value in values.tolist():
        if value is None or pd.isna(value):
            out.append(None)
        else:
            timestamp = pd.Timestamp(value)
            if timestamp.tzinfo is None:
                timestamp = timestamp.tz_localize("UTC")
            else:
                timestamp = timestamp.tz_convert("UTC")
            out.append(timestamp)
    return out


def dataframe_to_table(df: pd.DataFrame, schema: pa.Schema) -> pa.Table:
    arrays = []
    for field in schema:
        if field.name in df.columns:
            values = df[field.name]
        elif field.nullable:
            values = pd.Series([None] * len(df), dtype=object)
        else:
            raise ValueError(f"Missing required output column: {field.name}")

        if pa.types.is_timestamp(field.type):
            array = pa.array(_coerce_timestamp_values(values), type=field.type)
        elif pa.types.is_list(field.type):
            list_values = []
            for value in values.tolist():
                if value is None or (isinstance(value, float) and pd.isna(value)):
                    list_values.append([])
                elif hasattr(value, "tolist") and not isinstance(value, str):
                    list_values.append(
                        [str(item) for item in value.tolist() if str(item)]
                    )
                elif isinstance(value, str):
                    list_values.append([value] if value else [])
                else:
                    list_values.append([str(item) for item in value if str(item)])
            array = pa.array(list_values, type=field.type)
        else:
            array = pa.array(values.tolist(), type=field.type)
        arrays.append(array)
    return pa.Table.from_arrays(arrays, schema=schema)


def write_parquet(df: pd.DataFrame, path: str | Path, schema: pa.Schema) -> None:
    table = dataframe_to_table(df, schema)
    pq.write_table(table, path, compression="zstd")
