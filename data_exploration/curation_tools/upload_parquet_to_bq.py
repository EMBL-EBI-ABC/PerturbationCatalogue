from google.cloud import bigquery

def upload_parquet_to_bq(parquet_path, dataset_id, table_name, key_columns, verbose=True):
    client = bigquery.Client()
    target_table_id = f"{dataset_id}.{table_name}"
    staging_table_id = target_table_id+"_staging"
    if verbose:
        print(f"Loading Parquet file {parquet_path} to {staging_table_id}...")

    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
        )

    with open(parquet_path, "rb") as parquet_file:
        load_job = client.load_table_from_file(
            file_obj=parquet_file,
            destination=staging_table_id,
            job_config=job_config,
            rewind=True
        )
    load_job.result()
    if verbose:
        dest_table = client.get_table(staging_table_id)
        print(f"Loaded {dest_table.num_rows} rows to {staging_table_id}")

    # merge staging to target
    key_columns = [key.lower() for key in key_columns]
    update_columns = [col for col in dest_table.schema if col.name not in key_columns]
    update_columns = [col.name for col in update_columns if col.name != "row_id"]
    merge_staging_to_target(client, staging_table_id, target_table_id, key_columns, update_columns)
    
    

def merge_staging_to_target(client, staging_table_id, target_table_id, key_columns, update_columns):
    """
    Merge staging table into target table using BigQuery MERGE statement.

    Args:
        client: BigQuery client
        staging_table_id (str): fully qualified staging table id (project.dataset.table)
        target_table_id (str): fully qualified target table id
        key_columns (list[str]): list of column names used as unique keys
        update_columns (list[str]): list of columns to update on match
    """
    join_condition = " AND ".join([f"T.{col} = S.{col}" for col in key_columns])
    update_set = ", ".join([f"T.{col} = S.{col}" for col in update_columns])
    insert_columns = ", ".join(key_columns + update_columns)
    insert_values = ", ".join([f"S.{col}" for col in key_columns + update_columns])

    merge_sql = f"""
    MERGE `{target_table_id}` T
    USING `{staging_table_id}` S
    ON {join_condition}
    WHEN MATCHED THEN
      UPDATE SET {update_set}
    WHEN NOT MATCHED THEN
      INSERT ({insert_columns}) VALUES ({insert_values})
    """

    query_job = client.query(merge_sql)
    query_job.result()  # Wait for completion
    print(f"Merge completed: staging data merged into {target_table_id}")