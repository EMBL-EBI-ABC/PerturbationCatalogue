import argparse
import json
import subprocess

import requests


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
    parser.add_argument(
        "--id-field", required=True, help="Field to use as the document ID"
    )
    args = parser.parse_args()

    print("Removing the index...")
    response = requests.delete(f"{args.elastic_endpoint}/{args.elastic_index}")

    print("Creating the index...")
    index_settings = {
        "settings": {"index": {"number_of_shards": 1, "number_of_replicas": 1}},
        "mappings": {"properties": json.loads(args.field_properties)},
    }
    response = requests.put(
        f"{args.elastic_endpoint}/{args.elastic_index}",
        headers={"Content-Type": "application/json"},
        data=json.dumps(index_settings),
    )
    response.raise_for_status()

    print("Fetching the processed data...")
    subprocess.run(
        ["gsutil", "-q", "cp", args.jsonl_data, "/tmp/metadata.jsonl"], check=True
    )

    print("Preparing the data for loading...")
    with open("/tmp/metadata.jsonl", "r") as infile, open(
        "/tmp/formatted_metadata.jsonl", "w"
    ) as outfile:
        for line in infile:
            record = json.loads(line)
            _id = record.get(args.id_field)
            outfile.write(
                json.dumps({"index": {"_index": args.elastic_index, "_id": _id}}) + "\n"
            )
            outfile.write(line)

    print("Loading data into index...")
    with open("/tmp/formatted_metadata.jsonl", "rb") as data_file:
        response = requests.post(
            f"{args.elastic_endpoint}/{args.elastic_index}/_bulk",
            headers={"Content-Type": "application/x-ndjson"},
            data=data_file,
        )
    response.raise_for_status()

    print("All done.")


if __name__ == "__main__":
    main()
