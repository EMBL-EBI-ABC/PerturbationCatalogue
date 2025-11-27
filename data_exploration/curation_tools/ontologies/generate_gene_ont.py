# Generate gene ontology parquet file
# This script generates a parquet file containing gene ontology information
# and saves it to the specified path.
# Note: this requires internet access to fetch data from external databases and approximately 30 gb of RAM.

from argparse import ArgumentParser
from curation_tools.curation_tools import generate_gene_ont


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--save_parquet_path",
        default="data_exploration/curation_tools/ontologies/gene_ont.parquet",
    )
    parser.add_argument(
        "--json_url",
        default="https://ftp.ensembl.org/pub/current/json/homo_sapiens/homo_sapiens.json",
    )
    parser.add_argument(
        "--json_path",
        default="data_exploration/curation_tools/ontologies/homo_sapiens.json",
    )
    
    args = parser.parse_args()

    generate_gene_ont(
        json_url=args.json_url,
        json_path=args.json_path,
        save_parquet_path=args.save_parquet_path,
    )

