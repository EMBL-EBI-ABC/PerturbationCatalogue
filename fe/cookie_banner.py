import dash
from dash import html, dcc

store = dcc.Store(id="cookie-consent-store", storage_type="local")

layout = html.Div(
    id="cookie-banner",
    children=[
        html.Div(
            [
                html.P(
                    [
                        "This website requires cookies, and the limited processing of your personal data in order to function. By using the site you are agreeing to this as outlined in our ",
                        html.A(
                            "Privacy Notice",
                            href="https://www.ebi.ac.uk/data-protection/privacy-notice/perturbation-catalogue",
                            target="_blank",
                        ),
                        " and our ",
                        html.A(
                            "Terms of Use",
                            href="https://www.ebi.ac.uk/about/terms-of-use",
                            target="_blank",
                        ),
                        ".",
                    ],
                    className="mb-0 me-3",
                ),
                html.Button(
                    "Accept",
                    id="accept-cookies",
                    className="btn btn-success me-2",
                ),
                html.Button("Reject", id="reject-cookies", className="btn btn-danger"),
            ],
            className="alert alert-dark d-flex justify-content-center align-items-center gap-2 px-3 py-2",
            style={
                "position": "fixed",
                "bottom": "0",
                "left": "0",
                "width": "100%",
                "zIndex": "1050",
                "borderRadius": "0",
                "margin": "0",  # Remove any potential margin
                "padding": "0",  # Ensure no padding is added by the container
                "min-height": "70px",
            },
        )
    ],
    hidden=True,
)


def register_callbacks(app):
    @app.callback(
        [
            dash.Output("cookie-banner", "hidden"),
            dash.Output("cookie-consent-store", "data"),
        ],
        [
            dash.Input("accept-cookies", "n_clicks"),
            dash.Input("reject-cookies", "n_clicks"),
        ],
        [dash.State("cookie-consent-store", "data")],
        prevent_initial_call=False,
    )
    def handle_cookie_consent(accept_clicks, reject_clicks, store_data):
        ctx = dash.callback_context
        if store_data:
            return True, dash.no_update  # Hide banner if choice was already made.

        if ctx.triggered:
            button_id = ctx.triggered[0]["prop_id"].split(".")[0]
            if button_id in ["accept-cookies", "reject-cookies"]:
                return True, {"cookie_consent": button_id}

        return False, dash.no_update  # Show banner if no choice stored.
