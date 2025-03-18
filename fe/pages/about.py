import dash
import os
from . import _iframe

dash.register_page(
    __name__,
    path="/about",
    name="About",
    button="Contact Us",
    description="Here you can find more information about how to get help or propose new features.",
    icon="bi-question-circle",
)

api_base_url = os.getenv("PERTURBATION_CATALOGUE_BE")
layout = _iframe.layout(f"{api_base_url}/redoc")
