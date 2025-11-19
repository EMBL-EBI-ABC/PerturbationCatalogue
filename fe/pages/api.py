import dash
from dash import html
import os

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
    # Redirect to the backend API documentation
    docs_url = f"{api_base_url}/docs"
    layout = html.Div(
        [
            html.Meta(httpEquiv="refresh", content=f"0; url={docs_url}"),
            html.Div(
                [
                    html.P("Redirecting to API documentation..."),
                    html.P(
                        [
                            "If you are not redirected automatically, ",
                            html.A("click here", href=docs_url, target="_blank"),
                            ".",
                        ]
                    ),
                ],
                style={"textAlign": "center", "padding": "2em"},
            ),
        ]
    )
else:
    # Fallback layout if backend URL is not set
    layout = html.Div(
        html.Div(
            "API documentation is not available. Please set PERTURBATION_CATALOGUE_BE environment variable.",
            style={"textAlign": "center", "padding": "2em"},
        )
    )
