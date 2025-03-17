import os

import dash
from dash import html

from .elastic_table import ElasticTable, Column


# Common parameters.

api_base_url = os.getenv("PERTURBATION_CATALOGUE_BE")


# DepMap.


def high_dependency_genes(data):
    """Dynamic layout for the list of high dependency genes."""
    print(data)
    return " ".join([g["name"] for g in data])


depmap_table = ElasticTable(
    id="depmap",
    api_endpoint=f"{api_base_url}/depmap/search",
    columns=[
        # Special column, displayed nowhere, used only for cross-ref filtering.
        Column(
            field_name="xref",
            display_name="Cross-references",
            filterable=True,
            display_table=False,
            display_details=False,
        ),
        # Special columns: title and subtitle, displayed in both table and detail views.
        Column(
            field_name="ModelID",
            display_name="Model ID",
            display_details="title",
            display_table=lambda model_id: html.A(
                model_id,
                href=f"/data-portal/depmap/{model_id}",
                className="text-decoration-none text-nowrap",
            ),
        ),
        Column(
            field_name="OncotreePrimaryDisease",
            display_name="Oncotree Primary Disease",
            sortable=True,
            default_sort="asc",
            display_details="subtitle",
        ),
        # Data columns.
        Column(
            field_name="CellLineName",
            display_name="Cell Line Name",
            display_details="text",
        ),
        Column(
            field_name="OncotreeLineage",
            display_name="Oncotree Lineage",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="OncotreeSubtype",
            display_name="Oncotree Subtype",
            display_details="text",
        ),
        Column(
            field_name="Age",
            display_name="Age",
            sortable=True,
            display_details="text",
        ),
        Column(
            field_name="AgeCategory",
            display_name="Age Category",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="Sex",
            display_name="Sex",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="PrimaryOrMetastasis",
            display_name="Primary or Metastasis",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="SampleCollectionSite",
            display_name="Sample Collection Site",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="CatalogNumber",
            display_name="Catalog Number",
            display_details="text",
        ),
        Column(
            field_name="ModelType",
            display_name="Model Type",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="high_dependency_genes",
            display_name="High Dependency Genes",
            display_table=high_dependency_genes,
            display_details="text",
        ),
    ],
    details_button_name="View on DepMap",
    details_button_link=lambda ModelID: f"https://depmap.org/portal/model/{ModelID}/",
    title="DepMap",
)

dash.register_page(
    "data-portal-details-depmap",
    path_template="/data-portal/depmap/<record_id>",
    layout=depmap_table.details_layout,
)


# MaveDB.

mavedb_table = ElasticTable(
    id="mavedb",
    api_endpoint=f"{api_base_url}/mavedb/search",
    columns=[
        # Special columns: title and subtitle, for details view only.
        Column(
            field_name="title",
            display_table=False,
            display_details="title",
        ),
        Column(
            field_name="shortDescription",
            display_table=False,
            display_details="subtitle",
        ),
        # Special column: URN, to serve as the primary key and link construction.
        Column(
            field_name="urn",
            display_name="URN",
            display_details="text",
            display_table=lambda urn: html.A(
                urn,
                href=f"/data-portal/mavedb/{urn}",
                className="text-decoration-none text-nowrap",
            ),
        ),
        # Data columns.
        Column(
            field_name="sequenceType",
            display_name="Sequence Type",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="normalisedGeneName",
            display_name="Gene",
            display_details="text",
        ),
        Column(
            field_name="geneName",
            display_name="Domain",
            display_details="text",
        ),
        Column(
            field_name="geneCategory",
            display_name="Gene Category",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="publicationUrl",
            display_name="Publication URL",
            display_table=False,
            display_details="link",
        ),
        Column(
            field_name="publicationYear",
            display_name="Publication Year",
            display_details="text",
            filterable=True,
            sortable=True,
            default_sort="desc",
        ),
        Column(
            field_name="numVariants",
            display_name="Number of Variants",
            display_details="text",
            sortable=True,
        ),
    ],
    details_button_name="View on MaveDB",
    details_button_link=lambda urn: f"https://www.mavedb.org/score-sets/{urn}/",
    title="MaveDB",
)

dash.register_page(
    "data-portal-details-mavedb",
    path_template="/data-portal/mavedb/<record_id>",
    layout=mavedb_table.details_layout,
)


# Main data portal page.


def complete_layout():
    return html.Div(
        [datasource.table_layout() for datasource in (depmap_table, mavedb_table)]
    )


dash.register_page(
    __name__,
    path="/data-portal",
    name="Data Portal",
    button="Open Data Portal",
    description="The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities.",
    icon="bi-table",
    layout=complete_layout,
)


# Callbacks.


def register_callbacks(app):
    mavedb_table.register_callbacks(app)
    depmap_table.register_callbacks(app)
