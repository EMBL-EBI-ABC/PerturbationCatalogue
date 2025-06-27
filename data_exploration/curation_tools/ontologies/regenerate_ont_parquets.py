# regenerate ontologies parquets

from biomart import BiomartServer
import bionty as bt
import pandas as pd
from datetime import datetime
import os

# this is needed for lamindb to connect
os.environ["PGGSSENCMODE"] = "disable"


def regenerate_gene_ont() -> pd.DataFrame:

    server = BiomartServer("http://www.ensembl.org/biomart")
    server.verbose = True
    mart = server.datasets["hsapiens_gene_ensembl"]

    response = mart.search(
        {
            "attributes": [
                "ensembl_gene_id",
                "hgnc_symbol",
                "gene_biotype",
                "external_synonym",
                "description",
                "chromosome_name",
                "start_position",
                "end_position",
                "strand",
                "transcription_start_site"
            ]
        }
    )

    biomart_lines = []
    for line in response.iter_lines():
        line = line.decode("utf-8")
        biomart_lines.append(line.split("\t"))

    all_genes = pd.DataFrame(
        biomart_lines,
        columns=[
            "ensembl_gene_id",
            "hgnc_symbol",
            "gene_biotype",
            "external_synonym",
            "description",
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
            "transcription_start_site"
        ],
    )

    # dataframe is exploded on synonyms
    # collapse the synonyms on "|" to create a single row for each ENSG
    all_genes = (
        all_genes.groupby(
            [
                "ensembl_gene_id",
                "gene_biotype",
                "description",
                "chromosome_name",
                "start_position",
                "end_position",
                "strand",
                "transcription_start_site",
            ]
        )
        .agg(
            {
                "external_synonym": lambda x: "|".join(set(x)),
                "hgnc_symbol": lambda x: "|".join(set(x)),
            }
        ).reset_index()
    )
    
    all_genes['external_synonym'] = all_genes['external_synonym'].str.strip("|")
    all_genes['hgnc_symbol'] = all_genes['hgnc_symbol'].str.strip("|")

    # create a new column for condensed gene coordinates
    all_genes["gene_coord"] = all_genes.apply(
        lambda x: f"chr{x['chromosome_name']}:{x['start_position']}-{x['end_position']};{x['strand']}",
        axis=1,
    )

    # rename columns
    all_genes = all_genes.rename(
        columns={
            "hgnc_symbol": "symbol",
            "gene_biotype": "biotype",
            "external_synonym": "synonyms",
        }
    )
    # reorder columns
    all_genes = all_genes[
        [
            "ensembl_gene_id",
            "symbol",
            "synonyms",
            "biotype",
            "description",
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
            "gene_coord",
            "transcription_start_site",
        ]
    ]

    # add control gene as a new row
    all_genes = pd.concat(
        [
            all_genes,
            pd.DataFrame(
                {
                    "ensembl_gene_id": ["control"],
                    "symbol": ["control"],
                    "synonyms": [None],
                    "biotype": ["control"],
                    "description": ["control"],
                }
            ),
        ],
        ignore_index=True,
    )

    # replace empty strings with None
    all_genes = all_genes.replace("", None)

    return all_genes


# import ontologies from bionty
bt.CellType.import_source()
bt.CellLine.import_source()
bt.Tissue.import_source()
bt.Disease.import_source()

# get all ontologies
all_genes = regenerate_gene_ont()
all_cell_types = bt.CellType.filter().df()
all_cell_lines = bt.CellLine.filter().df()
all_diseses = bt.Disease.filter().df()
all_tissues = bt.Tissue.filter().df()

# subset the columns to keep
all_cell_types = all_cell_types[
    ["name", "ontology_id", "synonyms", "description"]
].reset_index(drop=True)
all_cell_lines = all_cell_lines[
    ["name", "ontology_id", "synonyms", "description"]
].reset_index(drop=True)
all_diseses = all_diseses[
    ["name", "ontology_id", "synonyms", "description"]
].reset_index(drop=True)
all_tissues = all_tissues[
    ["name", "ontology_id", "synonyms", "description"]
].reset_index(drop=True)


# add timestamp to all dataframes
all_genes["created_at"] = datetime.now()
all_cell_types["created_at"] = datetime.now()
all_cell_lines["created_at"] = datetime.now()
all_diseses["created_at"] = datetime.now()
all_tissues["created_at"] = datetime.now()

# save to parquet
save_dir = "."
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

all_genes.to_parquet(f"{save_dir}/genes.parquet")
all_cell_types.to_parquet(f"{save_dir}/cell_types.parquet")
all_cell_lines.to_parquet(f"{save_dir}/cell_lines.parquet")
all_diseses.to_parquet(f"{save_dir}/diseases.parquet")
all_tissues.to_parquet(f"{save_dir}/tissues.parquet")
