import dash
from dash import html
import dash_bootstrap_components as dbc

from ._order import get_pages

dash.register_page(
    __name__,
    path="/",
    name="Home",
    button=None,
    description=None,
    icon="bi-house",
)


def layout():
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
                            "object-fit": "cover",
                            "height": "clamp(200px, 25vw, 400px)",
                        },
                    ),
                    html.H1(
                        "Perturbation Catalogue",
                        style={
                            "font-size": "clamp(200%, 8vw, 300%)",
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
                                                html.H5(page["name"]),
                                                html.P(page["description"]),
                                                dbc.Button(
                                                    page["button"],
                                                    href=page["supplied_path"],
                                                    color="success",
                                                ),
                                            ]
                                        )
                                    ]
                                ),
                                md=4,
                                className="mb-4",
                            )
                            for page in get_pages(
                                include_home=False, require_button=True
                            )
                        ]
                    )
                ]
            ),
        ],
        style={"width": "100%", "height": "100%"},
    )
