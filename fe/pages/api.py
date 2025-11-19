import dash
import os
from pages import _iframe

dash.register_page(
    __name__,
    path="/api",
    name="API documentation",
    button="API documentation",
    description="Here you can find API documentation",
    icon="bi-terminal",
)

# Use the same environment variable as the rest of the app
api_base_url = os.getenv(
    "PERTURBATION_CATALOGUE_BE",
    "https://perturbation-catalogue-be-november-prototype-959149465821.europe-west2.run.app",
)
if api_base_url:
    layout = _iframe.layout(f"{api_base_url}/docs")
else:
    # Fallback layout if backend URL is not set
    from dash import html

    layout = html.Div(
        "API documentation is not available. Please set PERTURBATION_CATALOGUE_BE environment variable."
    )
