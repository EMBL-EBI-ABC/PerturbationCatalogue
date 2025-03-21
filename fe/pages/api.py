import dash
import os
from . import _iframe

dash.register_page(
    __name__,
    path="/api",
    name="API documentation",
    button="API documentation",
    description="Here you can find API documentation",
    icon="bi-terminal",
)

api_base_url = os.getenv("PERTURBATION_CATALOGUE_BE")
layout = _iframe.layout(f"{api_base_url}/redoc")
