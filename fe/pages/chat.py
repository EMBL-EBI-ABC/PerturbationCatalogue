import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
import os

from utils import BACKEND_URL

dash.register_page(
    __name__,
    path="/chat",
    relative_path="/chat",
    name="AI Explorer (Beta)",
    description="Ask questions about perturbation data using natural language.",
)

layout = dbc.Container(
    [
        # Hidden Graph to force Dash to load Plotly.js
        dcc.Graph(style={"display": "none"}, id="chat-plotly-loader"),
        # Hidden div to pass backend URL to JavaScript
        html.Div(
            id="chat-backend-url",
            **{"data-url": BACKEND_URL},
            style={"display": "none"},
        ),
        # Page title
        html.Div(
            [
                html.H3("AI Explorer (Beta)", className="mb-1", style={"fontWeight": "700"}),
                html.P(
                    "Ask questions about perturbation data in natural language",
                    className="text-muted mb-0",
                ),
            ],
            className="mt-4 mb-3",
        ),
        # Chat panel
        html.Div(
            [
                # Messages area
                html.Div(id="chat-messages", className="chat-messages"),
                # Spinner
                html.Div(
                    [
                        html.I(className="bi bi-flower1 spinner-icon"),
                        html.Span("", className="spinner-status"),
                    ],
                    id="chat-spinner",
                    className="chat-spinner",
                    style={"display": "none"},
                ),
                # Input area
                html.Div(
                    [
                        dcc.Input(
                            id="chat-input",
                            type="text",
                            placeholder="Ask about perturbation data...",
                            className="chat-input",
                            autoComplete="off",
                        ),
                        html.Button(
                            [html.I(className="bi bi-send-fill")],
                            id="chat-send-btn",
                            className="chat-send-btn",
                        ),
                    ],
                    className="chat-input-area",
                ),
                # Suggestion buttons
                html.Div(
                    [
                        html.Span(
                            "Try: ",
                            className="text-muted me-2",
                            style={"fontSize": "0.85rem"},
                        ),
                        html.Button(
                            "What data modalities are available?",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "What data modalities are available and how many datasets are there for each?"
                            },
                        ),
                        html.Button(
                            "BRCA2 perturbation data",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "What perturbation data is available for BRCA2?"
                            },
                        ),
                        html.Button(
                            "Show CRISPR screen results for TP53",
                            className="chat-suggestion-btn",
                            **{"data-query": "Show me CRISPR screen results for TP53"},
                        ),
                    ],
                    id="chat-suggestions",
                    className="chat-suggestions",
                ),
            ],
            className="chat-panel",
        ),
        # Data portal
        html.Div(
            [
                html.Div(
                    [
                        html.H5(
                            "Data Portal", className="mb-0", style={"fontWeight": "600"}
                        ),
                        html.Button(
                            "Clear all",
                            id="chat-clear-portal-btn",
                            className="btn btn-sm btn-outline-secondary",
                            style={"display": "none"},
                        ),
                    ],
                    className="d-flex justify-content-between align-items-center mb-3",
                ),
                html.Div(id="chat-data-portal", className="data-portal-content"),
                html.Div(
                    html.P(
                        "Visualizations will appear here as you ask questions",
                        className="text-muted text-center py-4",
                    ),
                    id="chat-portal-placeholder",
                ),
            ],
            className="data-portal-panel",
        ),
    ],
    className="content-container pb-4",
)
