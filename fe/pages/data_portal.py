from collections import namedtuple

import dash

from .elastic_table import ElasticTable, ColumnDefinition

# MaveDB.

columns = [
    ColumnDefinition(
        field_name="urn",
        display_name="URN",
        filterable=False,
        sortable=False,
        default_sort=False,
    ),
    ColumnDefinition(
        field_name="sequenceType",
        display_name="Sequence Type",
        filterable=True,
        sortable=False,
        default_sort=False,
    ),
    ColumnDefinition(
        field_name="geneName",
        display_name="Gene Name",
        filterable=False,
        sortable=False,
        default_sort=False,
    ),
    ColumnDefinition(
        field_name="geneCategory",
        display_name="Gene Category",
        filterable=True,
        sortable=False,
        default_sort=False,
    ),
    ColumnDefinition(
        field_name="publicationYear",
        display_name="Publication Year",
        filterable=True,
        sortable=True,
        default_sort=True,
    ),
    ColumnDefinition(
        field_name="numVariants",
        display_name="Number of Variants",
        filterable=False,
        sortable=True,
        default_sort=False,
    ),
]

mavedb_table = ElasticTable(
    api_endpoint="https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search",
    columns=columns,
)

dash.register_page(
    __name__,
    path="/data-portal",
    name="Data Portal",
    button="Open Data Portal",
    description="The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities.",
    icon="bi-table",
    layout=mavedb_table.table_layout,
)

dash.register_page(
    "data-portal-details",
    path_template="/data-portal/<urn>",
    layout=mavedb_table.details_layout,
)

# Callbacks.


def register_callbacks(app):
    mavedb_table.register_callbacks(app)
