import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc

from utils import BACKEND_URL

dash.register_page(
    __name__,
    path="/login",
    name="Login",
    description="Sign in to access the AI Explorer.",
)

layout = html.Div(
    [
        # Hidden div to pass backend URL to clientside callback
        html.Div(
            id="login-backend-url",
            **{"data-url": BACKEND_URL},
            style={"display": "none"},
        ),
        # Hidden store for redirect trigger
        dcc.Store(id="login-redirect"),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.I(
                                    className="bi bi-shield-lock",
                                    style={
                                        "fontSize": "1.5rem",
                                        "color": "#007B53",
                                    },
                                ),
                                html.Span(
                                    "AI Explorer",
                                    style={
                                        "fontWeight": "700",
                                        "fontSize": "1.1rem",
                                        "color": "#1a1d21",
                                    },
                                ),
                            ],
                            style={
                                "display": "flex",
                                "alignItems": "center",
                                "gap": "0.5rem",
                                "justifyContent": "center",
                                "marginBottom": "0.25rem",
                            },
                        ),
                        html.P(
                            "Sign in to continue",
                            style={
                                "color": "#54585A",
                                "fontSize": "0.9rem",
                                "textAlign": "center",
                                "margin": "0",
                            },
                        ),
                    ],
                    style={"marginBottom": "1.5rem"},
                ),
                # Error message
                html.Div(
                    id="login-error",
                    style={
                        "color": "#A6093D",
                        "fontSize": "0.85rem",
                        "textAlign": "center",
                        "marginBottom": "1rem",
                        "minHeight": "1.2em",
                    },
                ),
                # Email field
                dbc.Label("Email", html_for="login-email", className="mb-1"),
                dbc.Input(
                    id="login-email",
                    type="email",
                    placeholder="you@example.com",
                    className="mb-3",
                    autofocus=True,
                ),
                # Password field
                dbc.Label("Password", html_for="login-password", className="mb-1"),
                dbc.Input(
                    id="login-password",
                    type="password",
                    placeholder="Password",
                    className="mb-4",
                ),
                # Submit button
                dbc.Button(
                    "Sign in",
                    id="login-btn",
                    color="success",
                    className="w-100",
                    style={
                        "backgroundColor": "#007B53",
                        "borderColor": "#007B53",
                        "fontWeight": "600",
                    },
                ),
            ],
            style={
                "maxWidth": "380px",
                "margin": "80px auto",
                "padding": "2rem",
                "backgroundColor": "#fff",
                "borderRadius": "12px",
                "border": "1px solid #e9ecef",
                "boxShadow": "0 2px 12px rgba(0,0,0,0.06)",
            },
        ),
    ]
)
