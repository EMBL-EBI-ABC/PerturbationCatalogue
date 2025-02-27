import dash

from .elastic_table import ElasticTable, Column

# MaveDB.

mavedb_table = ElasticTable(
    api_endpoint="https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search",
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
        ),
        # Data columns.
        Column(
            field_name="sequenceType",
            display_name="Sequence Type",
            display_details="text",
            filterable=True,
        ),
        Column(
            field_name="geneName",
            display_name="Gene Name",
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
    details_page_url=lambda urn: f"/data-portal/{urn}",
    details_button_name="View on MaveDB",
    details_button_link=lambda urn: f"https://www.mavedb.org/score-sets/{urn}/",
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
