import argparse
from collections import Counter
import json
import pandas as pd


def process_json_to_csv(input_filename, output_filename):
    # Load JSON data from file.
    with open(input_filename, "r") as f:
        data = json.load(f)

    # Prepare a list to store the rows.
    rows = []

    # Keep track of the various target gene combinations.
    target_gene_cnt = Counter()

    # Navigate through the JSON structure and extract required fields.
    for experiment_set in data.get("experimentSets", []):
        for experiment in experiment_set.get("experiments", []):
            for score_set in experiment.get("scoreSets", []):
                # Analyse target genes.
                target_genes = score_set.get("targetGenes", [])
                match len(target_genes):
                    case 0:
                        target_gene_category = "No target gene"
                    case 1:
                        match target_genes[0]["targetSequence"]["taxonomy"]["taxId"]:
                            case 9606:
                                target_gene_category = "Single target gene, human"
                            case _:
                                target_gene_category = "Single target gene, non-human"
                    case _:
                        target_gene_category = "Multiple target genes"
                # Track the counts of different target gene cases.
                target_gene_cnt[target_gene_category] += 1
                # Skip all cases except when a single human gene is targeted.
                if target_gene_category != "Single target gene, human":
                    continue
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

    # Output target gene information.
    for key, value in target_gene_cnt.most_common():
        print(key, value)

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
