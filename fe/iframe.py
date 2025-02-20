from dash import html

# URLs for iframes.
DASHBOARDS_URL = "https://lookerstudio.google.com/embed/reporting/86ab32f9-151b-4f91-87eb-060a22f2f890/page/QJFZE"
DATA_ANALYTICS_URL = "https://lookerstudio.google.com/embed/reporting/8e98079d-144b-4583-a108-844b2ed3adf7/page/6eWZE"
ABOUT_URL = "https://perturbation-catalogue-be-959149465821.europe-west2.run.app/redoc"


# Layout for iframe-based pages.
def layout(url):
    return html.Div(
        [
            html.Iframe(
                src=url,
                style={
                    "width": "100%",
                    "height": "100%",
                    "border": "none",
                },
            )
        ],
        style={
            "width": "100%",
            "height": "100%",
            "border": "none",
            "display": "flex",
            "flexDirection": "column",
        },
        className="mt-0",
    )
