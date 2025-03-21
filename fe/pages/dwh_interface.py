import dash
from dash import html, callback, Output, Input
import dash_bootstrap_components as dbc

dash.register_page(
    __name__,
    path="/dwh_interface",
    name="DWH Interface",
    button="DWH Interface",
    description="Data warehousing advance query interfaces and presentations",
    icon="bi-database-fill-gear",
)

tab1_content = dbc.Card(
    dbc.CardBody(
        [
            html.P("Write SQL query in a box below", className="card-text"),
            dbc.Textarea(
                id="sql_query",
                size="lg",
                className="mb-3",
                rows=15,
                debounce=True,
                style={
                    "backgroundColor": "#000",
                    "border": "1px solid #000",
                    "color": "#00ff00",
                    "padding": "8px",
                    "fontFamily": "courier new",
                },
            ),
            dbc.Button("Validate SQL", color="info"),
            dbc.Button("Run Query", color="success", style={"marginLeft": "5px"}),
        ]
    ),
    className="mt-3",
)

tab2_content = dbc.Card(
    dbc.CardBody(
        [
            html.P("This is tab 2!", className="card-text"),
            dbc.Button("Don't click here", color="danger"),
        ]
    ),
    className="mt-3",
)

layout = dbc.Container(
    dbc.Row(
        dbc.Col(
            dbc.Tabs(
                [
                    dbc.Tab(tab1_content, label="SQL Interface"),
                    dbc.Tab(tab2_content, label="Advance Query Builder"),
                    dbc.Tab(
                        "This tab's content is never seen",
                        label="Natural Language Query Interface",
                    ),
                ]
            ),
            md=12,
        ),
        style={"marginTop": "2em", "marginBottom": "2em"},
    )
)


@callback(
    Output("sql_query", "valid"),
    Input("sql_query", "value"),
)
def validate_query(query):
    print(query)
    return True
