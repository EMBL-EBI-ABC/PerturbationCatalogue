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
                            "What perturbation data exists for TP53?",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "What perturbation data is available for TP53 across all modalities?"
                            },
                        ),
                        html.Button(
                            "Tell me about the BRCA2 protein",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "Look up the BRCA2 protein and show me a gene card with its function, domains, and disease associations"
                            },
                        ),
                        html.Button(
                            "Show 3D structure of KRAS",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "Show me the predicted 3D protein structure for KRAS"
                            },
                        ),
                        html.Button(
                            "Is EGFR druggable?",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "Use Pharos to check if EGFR is druggable. What is its target development level?"
                            },
                        ),
                        html.Button(
                            "Papers on CRISPR screens in cancer",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "Find recent papers about CRISPR screens in cancer"
                            },
                        ),
                        html.Button(
                            "Known variants of BRCA1",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "What are the known protein variants and mutagenesis data for BRCA1?"
                            },
                        ),
                        html.Button(
                            "Disease associations for MYC",
                            className="chat-suggestion-btn",
                            **{
                                "data-query": "What diseases are associated with MYC according to Open Targets?"
                            },
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
