import argparse
import json
import pandas as pd

def process_json_to_csv(input_filename, output_filename):
    # Load JSON data from file.
    with open(input_filename, "r") as f:
        data = json.load(f)
    
    # Prepare a list to store the rows.
    rows = []

    # Navigate through the JSON structure and extract required fields.
    for experiment_set in data.get("experimentSets", []):
        for experiment in experiment_set.get("experiments", []):
            for score_set in experiment.get("scoreSets", []):
                # Append row data with selected fields.
                rows.append({
                    "title": score_set.get("title", ""),
                    "methodText": score_set.get("methodText", ""),
                    "abstractText": score_set.get("abstractText", ""),
                    "shortDescription": score_set.get("shortDescription", ""),
                    "urn": score_set.get("urn", ""),
                    "numVariants": score_set.get("numVariants", 0)
                })

    # Convert list of rows to a DataFrame.
    df = pd.DataFrame(rows)
    
    # Save DataFrame to JSONL.
    with open(output_filename, 'w') as f:
        for row in rows:
            f.write(json.dumps(row) + '\n')

def main():
    # Set up argument parser.
    parser = argparse.ArgumentParser(description="Process JSON to CSV")
    parser.add_argument("--input-filename", type=str, required=True, help="Input JSON filename")
    parser.add_argument("--output-filename", type=str, required=True, help="Output CSV filename")
    args = parser.parse_args()
    
    # Process JSON and save to CSV.
    process_json_to_csv(args.input_filename, args.output_filename)

if __name__ == "__main__":
    main()
