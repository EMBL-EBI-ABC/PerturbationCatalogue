"""Infer a JSON Schema from a MaveDB main.json export (via genson).

Usage:
    python utils/infer_schema.py [--input data/v4/main.json] > data/v4/schema/json_schema.json
"""
import argparse
import json

from genson import SchemaBuilder

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Infer a JSON Schema from a MaveDB JSON export.")
    parser.add_argument("--input", default="data/v4/main.json", help="Path to main.json (default: data/v4/main.json)")
    args = parser.parse_args()

    with open(args.input, encoding='utf-8') as f:
        builder = SchemaBuilder()
        builder.add_object(json.load(f))
        print(json.dumps(builder.to_schema(), indent=2))
