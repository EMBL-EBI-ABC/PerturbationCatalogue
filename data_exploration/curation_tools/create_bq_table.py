import ibis
import re
import argparse
import polars as pl

def parquet_to_bq_type(parquet_dtype):
    """Map parquet datatypes to BigQuery SQL types."""
    # This mapping can be extended based on your schema
    parquet_str = str(parquet_dtype).lower()
    if 'int' in parquet_str:
        return 'INT64'
    if 'float' in parquet_str or 'double' in parquet_str or 'decimal' in parquet_str:
        return 'FLOAT64'
    if 'string' in parquet_str or 'text' in parquet_str:
        return 'STRING'
    if 'boolean' in parquet_str or 'bool' in parquet_str:
        return 'BOOL'
    if 'timestamp' in parquet_str:
        return 'TIMESTAMP'
    if 'date' in parquet_str:
        return 'DATE'
    if 'time' in parquet_str:
        return 'TIME'
    # fallback
    return 'STRING'

def generate_create_table_sql(
    table_name: str,
    schema: ibis.expr.schema.Schema,
    dataset_name: str = None,
    partition_column: str = None,
    partition_range_start: int = None,
    partition_range_end: int = None,
    partition_range_interval: int = None,
    cluster_columns: list = None,
):
    # Compose full table name with dataset if provided
    full_table_name = f"{dataset_name}.{table_name}" if dataset_name else table_name
    # Generate column definitions
    columns_sql = []
    for col_name, col_type in schema.items():
        bq_type = parquet_to_bq_type(col_type)
        # BigQuery reserved keywords or spaces require backticks
        safe_col_name = f"`{col_name}`" if re.match(r'\W', col_name) else col_name
        columns_sql.append(f"{safe_col_name} {bq_type}")
    columns_def = ",\n  ".join(columns_sql)

    # Prepare partition clause
    partition_clause = ""
    if partition_column:
        # Validate partition range parameters
        if (partition_range_start is None or partition_range_end is None or partition_range_interval is None):
            raise ValueError("For integer range partitioning, you must specify start, end, and interval.")
        # BigQuery requires partition column identifier to be backticked if needed
        partition_col_safe = f"`{partition_column}`" if re.match(r'\W', partition_column) else partition_column
        partition_clause = (
            f"\nPARTITION BY RANGE_BUCKET({partition_col_safe}, GENERATE_ARRAY("
            f"{partition_range_start}, {partition_range_end}, {partition_range_interval}))"
        )

    # Prepare clustering clause
    cluster_clause = ""
    if cluster_columns:
        cluster_cols_safe = []
        for ccol in cluster_columns:
            cluster_cols_safe.append(f"`{ccol}`" if re.match(r'\W', ccol) else ccol)
        cluster_clause = f"\nCLUSTER BY {', '.join(cluster_cols_safe)}"

    create_table_sql = (
        f"CREATE TABLE IF NOT EXISTS {full_table_name} (\n"
        f"  {columns_def}\n"
        f"){partition_clause}{cluster_clause};"
    )
    return create_table_sql

def create_bq_table(project_id, dataset_name, ddl_sql):
    """Create a BigQuery table using the provided DDL SQL."""
    client = ibis.bigquery.connect(project_id=project_id, dataset_id=dataset_name, location='europe-west2')
    try:
        client.raw_sql(ddl_sql)
        print(f"Table {dataset_name} created successfully in {project_id}.{dataset_name}.")
    except Exception as e:
        print(f"Error creating table: {e}")

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description="Generate BigQuery CREATE TABLE SQL from parquet schema.")
    parser.add_argument('--parquet_file', type=str, required=True, help="Path to the parquet file.")
    parser.add_argument('--bq_project_id', type=str, required=True, help="BigQuery project ID.")
    parser.add_argument('--bq_dataset_name', type=str, required=True, help="Name of the BigQuery dataset.")
    parser.add_argument('--bq_table_name', type=str, required=True, help="Name of the BigQuery table to create.")
    parser.add_argument('--partition_column', type=str, help="Column to partition the table by.")
    parser.add_argument('--partition_range_start', type=int, help="Start of the partition range.")
    parser.add_argument('--partition_range_end', type=int, help="End of the partition range.")
    parser.add_argument('--partition_range_interval', type=int, help="Interval for the partition range.")
    parser.add_argument('--cluster_columns', type=str, nargs='*', help="List of columns to cluster the table by.")
    args = parser.parse_args()

    # Load parquet schema
    parquet_file_path = args.parquet_file
    bq_project_id = args.bq_project_id
    bq_dataset_name = args.bq_dataset_name
    bq_table_name = args.bq_table_name
    
    partition_column = args.partition_column if args.partition_column else None
    partition_range_start = args.partition_range_start if args.partition_range_start else None
    partition_range_end = args.partition_range_end if args.partition_range_end else None
    partition_range_interval = args.partition_range_interval if args.partition_range_interval else None
    cluster_columns = args.cluster_columns if args.cluster_columns else None

    # Load parquet schema
    df = pl.scan_parquet(parquet_file_path)
    schema = df.schema

    # Generate DDL SQL
    ddl_sql = generate_create_table_sql(
        dataset_name=bq_dataset_name,
        table_name=bq_table_name,
        schema=schema,
        partition_column=partition_column,
        partition_range_start=partition_range_start,
        partition_range_end=partition_range_end,
        partition_range_interval=partition_range_interval,
        cluster_columns=cluster_columns
    )
    
    # Print the generated SQL
    print("Generated CREATE TABLE SQL:")    
    print(ddl_sql)

    # Create the BigQuery table
    create_bq_table(bq_project_id, bq_dataset_name, ddl_sql)

    final_print = f"Table {bq_dataset_name}.{bq_table_name} created successfully"
    if partition_column:
        final_print += f", partitioned on {partition_column}"
    if cluster_columns:
        final_print += f", clustered on {', '.join(cluster_columns)}"

    print(final_print)
    
if __name__ == "__main__":
    main()