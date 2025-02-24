import dash

from . import _iframe

dash.register_page(
    __name__,
    path="/about",
    name="About",
    button="Contact Us",
    description="Here you can find more information about how to get help or propose new features.",
    icon="bi-question-circle",
)

layout = _iframe.layout(
    "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/redoc"
)
