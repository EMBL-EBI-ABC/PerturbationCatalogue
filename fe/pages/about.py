import dash

from . import _iframe

dash.register_page(__name__)

layout = _iframe.layout(
    "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/redoc"
)
