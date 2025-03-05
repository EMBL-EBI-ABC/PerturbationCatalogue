import os
import subprocess
import sys
import json
import argparse


def main():
    parser = argparse.ArgumentParser(description="Load data into Elasticsearch.")
    parser.add_argument(
        "--elastic-endpoint", required=True, help="Elasticsearch endpoint"
    )
    parser.add_argument(
        "--elastic-index", required=True, help="Elasticsearch index name"
    )
    parser.add_argument("--jsonl-data", required=True, help="Path to JSONL data file")
    parser.add_argument(
        "--field-properties", required=True, help="Field properties in JSON format"
    )
    args = parser.parse_args()

    ELASTIC_ENDPOINT = args.elastic_endpoint
    ELASTIC_INDEX = args.elastic_index
    JSONL_DATA = args.jsonl_data
    FIELD_PROPERTIES = args.field_properties

    print("Removing the index...")
    subprocess.run(
        [
            "curl",
            "-X",
            "DELETE",
            f"{ELASTIC_ENDPOINT}/{ELASTIC_INDEX}",
            "-H",
            "Content-Type: application/json",
        ],
        check=True,
    )

    print("\nCreating the index...")
    index_settings = {
        "settings": {"index": {"number_of_shards": 1, "number_of_replicas": 1}},
        "mappings": {"properties": json.loads(FIELD_PROPERTIES)},
    }
    subprocess.run(
        [
            "curl",
            "-X",
            "PUT",
            f"{ELASTIC_ENDPOINT}/{ELASTIC_INDEX}",
            "-H",
            "Content-Type: application/json",
            "-d",
            json.dumps(index_settings),
        ],
        check=True,
    )

    print("\nPreparing the data for loading...")
    subprocess.run(
        ["gsutil", "-q", "cp", JSONL_DATA, "/tmp/metadata.jsonl"], check=True
    )

    with open("/tmp/metadata.jsonl", "r") as infile, open(
        "/tmp/formatted_metadata.jsonl", "w"
    ) as outfile:
        for i, line in enumerate(infile, 1):
            outfile.write(
                json.dumps({"index": {"_index": ELASTIC_INDEX, "_id": str(i)}}) + "\n"
            )
            outfile.write(line)

    print("Loading data into index...")
    subprocess.run(
        [
            "curl",
            "-X",
            "POST",
            f"{ELASTIC_ENDPOINT}/{ELASTIC_INDEX}/_bulk",
            "-H",
            "Content-Type: application/x-ndjson",
            "--data-binary",
            "@/tmp/formatted_metadata.jsonl",
        ],
        check=True,
    )

    print("\nAll done.")


if __name__ == "__main__":
    main()
