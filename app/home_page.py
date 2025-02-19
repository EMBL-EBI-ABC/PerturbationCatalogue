from dash import html
import dash_bootstrap_components as dbc


def layout(pages):
    return html.Div(
        [
            html.Div(
                style={"position": "relative", "text-align": "center"},
                children=[
                    html.Img(
                        src="https://acxngcvroo.cloudimg.io/v7/https://www.embl.org/files/wp-content/uploads/roundels.png",
                        style={
                            "width": "100%",
                            "color": "rgb(24, 29, 25)",
                            "background-color": "rgb(0, 112, 73, 0.8)",
                        },
                    ),
                    html.H1(
                        "Perturbation Catalogue",
                        style={
                            "font-size": "300%",
                            "position": "absolute",
                            "top": "50%",
                            "left": "50%",
                            "transform": "translate(-50%, -50%)",
                            "color": "white",
                        },
                    ),
                ],
            ),
            html.H1("", className="text-center my-4"),
            dbc.Container(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(page.name),
                                                html.P(page.description),
                                                dbc.Button(
                                                    page.button,
                                                    href=f"/{page.selector}",
                                                    color="success",
                                                ),
                                            ]
                                        )
                                    ]
                                ),
                                md=4,
                                className="mb-4",
                            )
                            for page in pages
                            if page.button is not None
                        ],
                    ),
                ]
            ),
        ]
    )
