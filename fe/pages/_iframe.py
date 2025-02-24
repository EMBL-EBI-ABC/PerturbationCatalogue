from dash import html


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
            "flex": "1",
        },
        className="mt-0",
    )
