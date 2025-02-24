import dash

from . import _iframe

dash.register_page(
    __name__,
    path="/dashboards",
    name="Dashboards",
    button="Explore Dashboards",
    description="The Dashboards tab provides a visual overview of the existing data.",
    icon="bi-bar-chart-line",
)

layout = _iframe.layout(
    "https://lookerstudio.google.com/embed/reporting/86ab32f9-151b-4f91-87eb-060a22f2f890/page/QJFZE"
)
