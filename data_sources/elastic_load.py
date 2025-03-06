import os
import subprocess
import sys
import json


def main():
    if len(sys.argv) != 5:
        print(
            "Usage: python elastic_load.py <ELASTIC_ENDPOINT> <ELASTIC_INDEX> <JSONL_DATA> <FIELD_PROPERTIES>"
        )
        sys.exit(1)

    ELASTIC_ENDPOINT = sys.argv[1]
    ELASTIC_INDEX = sys.argv[2]
    JSONL_DATA = sys.argv[3]
    FIELD_PROPERTIES = sys.argv[4]

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
