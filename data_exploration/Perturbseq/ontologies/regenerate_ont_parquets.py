# regenerate ontologies parquets

from biomart import BiomartServer
import bionty as bt
import pandas as pd
from datetime import datetime
import os


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
        ],
    )

    # dataframe is exploded on the synonyms
    # collapse the synonyms on "|" to create a single row for each gene
    all_genes = (
        all_genes.groupby(["ensembl_gene_id", "gene_biotype", "description"])
        .agg(
            {
                "external_synonym": lambda x: "|".join(set(x)),
                "hgnc_symbol": lambda x: "|".join(set(x)),
            }
        )
        .reset_index()
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
        ["ensembl_gene_id", "symbol", "synonyms", "biotype", "description"]
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
                    "biotype": [None],
                    "description": ["control"],
                }
            ),
        ],
        ignore_index=True,
    )

    # replace empty strings with None
    all_genes = all_genes.replace("", None)

    # handle cases where there are two symbols mapping to the same ENSG
    # get the length of the symbol column
    sym_len = (
        all_genes["symbol"]
        .str.split("|")
        .apply(lambda x: len(x) if isinstance(x, list) else 0)
    )
    # get the second symbol for the genes with more than one symbol
    # this will be moved to the synonyms column
    second_sym = all_genes.iloc[sym_len[sym_len > 1].index]["symbol"].str.split(
        "|", expand=True
    )[1]

    # add the second symbol to the synonyms column
    all_genes.loc[sym_len[sym_len > 1].index, "synonyms"] = (
        all_genes.loc[sym_len[sym_len > 1].index, "synonyms"]
        + "|"
        + all_genes.loc[sym_len[sym_len > 1].index, "symbol"].str.split("|", expand=True)[1]
    )

    # replace the symbol with the first symbol
    all_genes.loc[sym_len[sym_len > 1].index, "symbol"] = all_genes.loc[
        sym_len[sym_len > 1].index, "symbol"
    ].str.split("|", expand=True)[0]

    
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
save_dir = "Ontologies"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

all_genes.to_parquet(f"{save_dir}/genes.parquet")
all_cell_types.to_parquet(f"{save_dir}/cell_types.parquet")
all_cell_lines.to_parquet(f"{save_dir}/cell_lines.parquet")
all_diseses.to_parquet(f"{save_dir}/diseases.parquet")
all_tissues.to_parquet(f"{save_dir}/tissues.parquet")
