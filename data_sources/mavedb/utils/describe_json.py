"""Print a key/type tree of a MaveDB main.json export.

Usage:
    python utils/describe_json.py [--input data/v4/main.json] > data/v4/schema/json_summary.txt
"""
import argparse
import json


def describe_json(data, indent=0):
    pad = '  ' * indent
    if isinstance(data, dict):
        for key, value in data.items():
            print(f"{pad}{key}: {type(value).__name__}")
            describe_json(value, indent + 1)
    elif isinstance(data, list):
        print(f"{pad}List[{len(data)}] of {type(data[0]).__name__}" if data else f"{pad}Empty List")
        if data:
            describe_json(data[0], indent + 1)
    else:
        print(f"{pad}{type(data).__name__}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print a key/type tree of a MaveDB JSON export.")
    parser.add_argument("--input", default="data/v4/main.json", help="Path to main.json (default: data/v4/main.json)")
    args = parser.parse_args()

    with open(args.input, 'r', encoding='utf-8') as f:
        describe_json(json.load(f))
