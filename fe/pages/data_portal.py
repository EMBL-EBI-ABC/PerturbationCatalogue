from collections import namedtuple

import dash

from .elastic_table import ElasticTable


FilterField = namedtuple("FilterField", ["id", "title"])


# MaveDB.

mavedb_table = ElasticTable(
    api_endpoint="https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search",
    columns=[
        ("URN", "urn"),
        ("Sequence Type", "sequenceType"),
        ("Gene Name", "geneName"),
        ("Gene Category", "geneCategory"),
        ("Publication Year", "publicationYear"),
        ("Number of Variants", "numVariants"),
    ],
    filter_fields=[
        FilterField(id="sequenceType", title="Sequence Type"),
        FilterField(id="geneCategory", title="Gene Category"),
        FilterField(id="publicationYear", title="Publication Year"),
    ],
    sortable_columns={
        "publicationYear": "Publication Year",
        "numVariants": "Number of Variants",
    },
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
