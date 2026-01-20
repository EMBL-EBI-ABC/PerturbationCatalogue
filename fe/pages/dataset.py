"""Dataset Details page showing all dataset metadata fields."""

from __future__ import annotations

import json
from typing import Any, Dict, Optional

import dash
from dash import Input, Output, callback, dcc, html
import dash_bootstrap_components as dbc

from utils import COLORS, fetch_dataset


dash.register_page(
    __name__,
    path_template="/dataset/<dataset_id>",
    name="Dataset Details",
    title="Dataset Details",
)

DATASET_STORE = "dataset-store"


def _format_dataset_id(dataset_id: str) -> str:
    """Format dataset ID by replacing underscores with spaces and capitalizing."""
    return dataset_id.replace("_", " ").title()


def _format_field_name(field_name: str) -> str:
    """Convert field name like 'cell_type_labels' to 'Cell Type'."""
    # Remove _labels or _ids suffix
    base_name = field_name.replace("_labels", "").replace("_ids", "")
    # Convert snake_case to Title Case
    return " ".join(word.capitalize() for word in base_name.split("_"))


def _ontology_id_to_url(ontology_id: str) -> str:
    """Convert ontology ID (e.g., 'UBERON:0002113') to a URL."""
    if not ontology_id or ":" not in ontology_id:
        return ontology_id
    
    # Split into prefix and ID
    parts = ontology_id.split(":", 1)
    if len(parts) != 2:
        return ontology_id
    
    prefix, term_id = parts
    # Convert to OBO format (UBERON:0002113 -> UBERON_0002113)
    obo_id = f"{prefix}_{term_id}"
    # Create OLS URL
    return f"https://www.ebi.ac.uk/ols4/ontologies/{prefix.lower()}/terms?iri=http://purl.obolibrary.org/obo/{obo_id}"


def _format_value(value: Any) -> str:
    """Format a value for display."""
    if isinstance(value, list):
        if not value:
            return "N/A"
        return ", ".join(str(v) for v in value)
    elif value is None:
        return "N/A"
    else:
        return str(value)


def _create_associated_datasets_display(associated_datasets: Any) -> html.Div:
    """Create display for associated datasets field."""
    if not associated_datasets:
        return _create_field_display("associated_datasets", None)
    
    # Parse JSON strings if needed
    items = []
    if isinstance(associated_datasets, list):
        for item in associated_datasets:
            if isinstance(item, str):
                try:
                    parsed = json.loads(item)
                    if isinstance(parsed, list):
                        items.extend(parsed)
                    else:
                        items.append(parsed)
                except json.JSONDecodeError:
                    continue
            elif isinstance(item, dict):
                items.append(item)
    
    if not items:
        return _create_field_display("associated_datasets", None)
    
    # Create display for each item
    display_items = []
    for item in items:
        description = item.get("dataset_description", "N/A")
        file_name = item.get("dataset_file_name", "")
        uri = item.get("dataset_uri", "")
        
        if uri and file_name:
            link = html.A(
                file_name,
                href=uri,
                target="_blank",
                rel="noopener noreferrer",
                className="text-decoration-none",
                style={"color": COLORS["primary"]},
            )
            display_items.append(
                html.Div(
                    [
                        html.Span(f"{description}: ", className="text-muted"),
                        link,
                    ],
                    className="mb-2",
                )
            )
        else:
            display_items.append(
                html.Div(
                    [
                        html.Span(f"{description}: ", className="text-muted"),
                        html.Span(file_name or "N/A", className="text-muted"),
                    ],
                    className="mb-2",
                )
            )
    
    return html.Div(
        [
            html.Dt(
                _format_field_name("associated_datasets") + ":",
                className="col-sm-4 fw-semibold",
            ),
            html.Dd(
                html.Div(display_items),
                className="col-sm-8",
            ),
        ],
        className="row mb-3",
    )


def _create_field_display(
    field_name: str, value: Any, ontology_ids: Optional[list] = None
) -> html.Div:
    """Create a display element for a field, with links if ontology IDs are available."""
    if isinstance(value, list) and not value:
        value = None

    formatted_value = _format_value(value)

    # If we have ontology IDs and values, create links
    if ontology_ids and isinstance(value, list) and len(value) == len(ontology_ids):
        links = []
        for label, ontology_id in zip(value, ontology_ids):
            if ontology_id:
                ontology_url = _ontology_id_to_url(str(ontology_id))
                links.append(
                    html.A(
                        str(label),
                        href=ontology_url,
                        target="_blank",
                        rel="noopener noreferrer",
                        className="text-decoration-none",
                        style={"color": COLORS["primary"]},
                    )
                )
            else:
                links.append(html.Span(str(label)))
            # Add comma separator except for last item
            if label != value[-1]:
                links.append(html.Span(", "))

        display_value = html.Span(links)
    elif ontology_ids and isinstance(value, list) and ontology_ids:
        # If we have some IDs but not matching length, try to match by index
        links = []
        for i, label in enumerate(value):
            if i < len(ontology_ids) and ontology_ids[i]:
                ontology_url = _ontology_id_to_url(str(ontology_ids[i]))
                links.append(
                    html.A(
                        str(label),
                        href=ontology_url,
                        target="_blank",
                        rel="noopener noreferrer",
                        className="text-decoration-none",
                        style={"color": COLORS["primary"]},
                    )
                )
            else:
                links.append(html.Span(str(label)))
            if i < len(value) - 1:
                links.append(html.Span(", "))
        display_value = html.Span(links)
    else:
        display_value = html.Span(formatted_value, className="text-muted")

    return html.Div(
        [
            html.Dt(
                _format_field_name(field_name) + ":",
                className="col-sm-4 fw-semibold",
            ),
            html.Dd(display_value, className="col-sm-8"),
        ],
        className="row mb-3",
    )


def layout(dataset_id: Optional[str] = None, **kwargs):
    """Page layout for dataset details."""
    return dbc.Container(
        [
            dcc.Store(id=DATASET_STORE, data={"dataset_id": dataset_id}),
            dcc.Loading(
                html.Div(id="dataset-content"),
                type="default",
                color=COLORS["primary"],
            ),
        ],
        className="py-4",
    )


@callback(
    Output("dataset-content", "children"),
    Input(DATASET_STORE, "data"),
)
def render_dataset(data: Optional[Dict[str, Any]]):
    """Render dataset details."""
    if not data or not data.get("dataset_id"):
        return html.Div(
            "No dataset ID provided.",
            className="text-muted text-center",
        )

    dataset_id = data["dataset_id"]
    formatted_dataset_id = _format_dataset_id(dataset_id)
    dataset_data, error = fetch_dataset(dataset_id)

    if error:
        return html.Div(
            [
                html.H1(
                    f"Dataset: {formatted_dataset_id}",
                    className="display-5 fw-bold mb-2 text-center",
                    style={"color": COLORS["primary"]},
                ),
                html.Div(
                    f"Error loading dataset: {error}",
                    className="alert alert-danger mt-4",
                ),
            ]
        )

    if not dataset_data:
        return html.Div(
            [
                html.H1(
                    f"Dataset: {formatted_dataset_id}",
                    className="display-5 fw-bold mb-2 text-center",
                    style={"color": COLORS["primary"]},
                ),
                html.Div(
                    "Dataset not found.",
                    className="text-muted text-center mt-4",
                ),
            ]
        )

    # Group fields: collect _labels and _ids pairs
    field_groups: Dict[str, Dict[str, Any]] = {}
    standalone_fields: Dict[str, Any] = {}

    for key, value in dataset_data.items():
        if key.endswith("_labels"):
            base_name = key.replace("_labels", "")
            if base_name not in field_groups:
                field_groups[base_name] = {}
            field_groups[base_name]["labels"] = value
            # Check if corresponding _ids exists
            ids_key = f"{base_name}_ids"
            if ids_key in dataset_data:
                field_groups[base_name]["ids"] = dataset_data[ids_key]
        elif key.endswith("_ids"):
            # Only process if we haven't already seen the _labels version
            base_name = key.replace("_ids", "")
            if base_name not in field_groups:
                field_groups[base_name] = {"ids": value}
        else:
            # Standalone field
            standalone_fields[key] = value

    # Build the field list
    field_elements = []

    # First, add standalone fields (excluding dataset_id which we show in the header)
    for key, value in sorted(standalone_fields.items()):
        if key != "dataset_id":
            # Special handling for associated_datasets
            if key == "associated_datasets":
                field_elements.append(_create_associated_datasets_display(value))
            else:
                field_elements.append(_create_field_display(key, value))

    # Then, add grouped fields (those with _labels and potentially _ids)
    for base_name in sorted(field_groups.keys()):
        group = field_groups[base_name]
        labels = group.get("labels")
        ids = group.get("ids")
        # Use _labels field name if labels exist, otherwise use _ids
        if labels is not None:
            field_elements.append(
                _create_field_display(f"{base_name}_labels", labels, ids)
            )
        elif ids is not None:
            # If only IDs exist, display them as standalone values
            field_elements.append(
                _create_field_display(f"{base_name}_ids", ids, None)
            )

    return html.Div(
        [
            html.H1(
                f"Dataset: {formatted_dataset_id}",
                className="display-5 fw-bold mb-4 text-center",
                style={"color": COLORS["primary"]},
            ),
            html.Dl(
                field_elements,
                className="mt-4",
                style={"maxWidth": "900px", "margin": "0 auto"},
            ),
        ]
    )
