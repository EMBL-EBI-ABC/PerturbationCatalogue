from collections import namedtuple

import dash
from dash import html, Input, Output, State, dcc
import dash_bootstrap_components as dbc
import requests

from .elastic_table import details, table
from .elastic_table.elastic_table import ElasticTable


FilterField = namedtuple("FilterField", ["id", "title"])
filter_fields = [
    FilterField(id="sequenceType", title="Sequence Type"),
    FilterField(id="geneCategory", title="Gene Category"),
    FilterField(id="publicationYear", title="Publication Year"),
]

sortable_columns = {
    "publicationYear": "Publication Year",
    "numVariants": "Number of Variants",
}

columns = [
    ("URN", "urn"),
    ("Sequence Type", "sequenceType"),
    ("Gene Name", "geneName"),
    ("Gene Category", "geneCategory"),
    ("Publication Year", "publicationYear"),
    ("Number of Variants", "numVariants"),
]

elastic_table = ElasticTable(
    api_endpoint="https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search",
    columns=columns,
    filter_fields=filter_fields,
    sortable_columns=sortable_columns,
    default_sort_field="publicationYear",
    default_sort_order="desc",
)


dash.register_page(
    __name__,
    path="/data-portal",
    name="Data Portal",
    button="Open Data Portal",
    description="The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities.",
    icon="bi-table",
    layout=table.layout(filter_fields),
)
dash.register_page(
    "data-portal-details", path_template="/data-portal/<urn>", layout=details.layout
)


def create_table_header(column_name, field_name, current_sort):
    return elastic_table.create_table_header(column_name, field_name, current_sort)


def create_table(data, sort_data):
    return elastic_table.create_table(data, sort_data)


def register_callbacks(app):
    elastic_table.register_callbacks(app)
