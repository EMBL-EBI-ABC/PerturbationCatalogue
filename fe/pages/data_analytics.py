import dash

from . import _iframe

dash.register_page(
    __name__,
    path="/data-analytics",
    name="Data Analytics",
    button="Data Analytics",
    description="The Data Analytics tab allows users to search the data using the Data Warehouse",
    icon="bi-graph-up",
)

layout = _iframe.layout(
    "https://lookerstudio.google.com/embed/reporting/8e98079d-144b-4583-a108-844b2ed3adf7/page/6eWZE"
)
