"""Dataset Details page showing all dataset metadata fields."""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional
from urllib.parse import quote, urlencode

import isodate
import dash
from dash import ALL, Input, Output, callback, dcc, html
import dash_bootstrap_components as dbc

from components.target_data_table import crispr_table, mave_heatmap, perturb_seq_table
from utils import (
    BACKEND_URL,
    COLORS,
    fetch_dataset,
    fetch_dataset_rows,
    reprocessed_badge,
)


dash.register_page(
    __name__,
    path_template="/dataset/<dataset_id>",
    name="Dataset Details",
    title="Dataset Details",
)

DATASET_STORE = "dataset-store"
DATASET_DATA_STORE = "dataset-data-store"
ROWS_PER_PAGE = 15
MAVE_POSITION_RANGE_SIZE = 20

# Search input IDs for Perturb-seq
PERTURB_SEARCH_PERTURBATION_GENE = "dataset-perturb-search-perturbation-gene"
PERTURB_SEARCH_EFFECT_GENE = "dataset-perturb-search-effect-gene"

# Search input IDs for CRISPR screen
CRISPR_SEARCH_PERTURBATION_GENE = "dataset-crispr-search-perturbation-gene"

# Download link IDs
DATASET_DOWNLOAD_LINK = "dataset-download-link"
DATASET_METADATA_DOWNLOAD_LINK = "dataset-metadata-download-link"
DATASET_PARQUET_DOWNLOAD_LINK = "dataset-parquet-download-link"
DATASET_FILTERED_DOWNLOAD_LINK = "dataset-filtered-download-link"

# Map display names to API endpoint modality names
MODALITY_DISPLAY_TO_API = {
    "Perturb-seq": "perturb-seq",
    "CRISPR screen": "crispr-screen",
    "MAVE": "mave",
}


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
    """Convert ontology ID (e.g., 'UBERON:0002113') to an OLS4 URL.

    OLS4 requires:
    - URL format: /ontologies/{ontology}/classes/{double-encoded-IRI}
    - Different ontologies use different IRI bases:
      - EFO: http://www.ebi.ac.uk/efo/EFO_XXXXXXX
      - OBO Foundry (UBERON, CL, etc.): http://purl.obolibrary.org/obo/PREFIX_ID
    """
    if not ontology_id or ":" not in ontology_id:
        return ontology_id

    parts = ontology_id.split(":", 1)
    if len(parts) != 2:
        return ontology_id

    prefix, term_id = parts
    prefix_upper = prefix.upper()
    prefix_lower = prefix.lower()

    # Ontologies that use their own IRI scheme (not OBO Foundry)
    # Format: prefix -> (iri_base, include_prefix_in_id)
    # If include_prefix_in_id is True, ID will be PREFIX_termid, else just termid
    non_obo_iri_bases = {
        "EFO": "http://www.ebi.ac.uk/efo/",
        "SWO": "http://www.ebi.ac.uk/swo/license/",
    }

    # Build the full IRI based on the ontology
    if prefix_upper in non_obo_iri_bases:
        iri = f"{non_obo_iri_bases[prefix_upper]}{prefix_upper}_{term_id}"
    else:
        # Default to OBO Foundry format
        iri = f"http://purl.obolibrary.org/obo/{prefix_upper}_{term_id}"

    # Double-encode the IRI for OLS4 URL path
    # First encode, then encode again
    encoded_iri = quote(quote(iri, safe=""), safe="")

    return f"https://www.ebi.ac.uk/ols4/ontologies/{prefix_lower}/classes/{encoded_iri}"


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


def _format_iso_duration(duration_str: str) -> str:
    """Convert ISO 8601 duration string to human-readable format.

    Examples:
        "P7DT0H0M0S" -> "Day 7"
        "P0DT12H0M0S" -> "Hour 12"
        "P1DT6H0M0S" -> "Day 1, Hour 6"
    """
    try:
        duration = isodate.parse_duration(duration_str)
        parts = []

        # Extract days and hours from timedelta
        days = duration.days
        hours = duration.seconds // 3600
        minutes = (duration.seconds % 3600) // 60

        if days > 0:
            parts.append(f"Day {days}")
        if hours > 0:
            parts.append(f"Hour {hours}")
        if minutes > 0:
            parts.append(f"Minute {minutes}")

        if not parts:
            return "Day 0"

        return ", ".join(parts)
    except (ValueError, isodate.ISO8601Error):
        return duration_str


def _format_timepoints(timepoints: List[str]) -> List[str]:
    """Format a list of ISO 8601 duration strings to human-readable format."""
    return [_format_iso_duration(tp) for tp in timepoints]


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
            dcc.Store(id=DATASET_DATA_STORE, data=None),
            dcc.Loading(
                html.Div(id="dataset-content"),
                type="circle",
                color=COLORS["primary"],
                target_components={"dataset-content": "children"},
            ),
        ],
        className="py-4",
    )


def _get_modality_from_dataset(dataset_data: Dict[str, Any]) -> Optional[str]:
    """Extract API modality name from dataset metadata."""
    data_modalities = dataset_data.get("data_modalities") or []
    for display_name in data_modalities:
        if display_name in MODALITY_DISPLAY_TO_API:
            return MODALITY_DISPLAY_TO_API[display_name]
    return None


def _dataset_download_url(
    dataset_id: str,
    modality: str,
    download_format: str = "csv.gz",
    filters: Optional[Dict[str, str]] = None,
) -> str:
    url = f"{BACKEND_URL}/v1/{modality}/{dataset_id}/download"
    params = {"format": download_format}
    params.update({key: value for key, value in (filters or {}).items() if value})
    return f"{url}?{urlencode(params)}"


def _dataset_download_link(
    label: str,
    url: str,
    link_id: str,
) -> html.A:
    return html.A(
        dbc.Button(
            [html.I(className="bi bi-download me-2"), label],
            color="primary",
            size="sm",
            style={
                "backgroundColor": COLORS["primary"],
                "borderColor": COLORS["primary"],
                "borderRadius": "6px",
            },
        ),
        id=link_id,
        href=url,
        target="_blank",
        rel="noopener noreferrer",
        className="text-decoration-none me-3",
    )


def _dataset_download_controls(dataset_id: str, modality: str) -> html.Div:
    return html.Div(
        [
            _dataset_download_link(
                "Download full data (Parquet)",
                _dataset_download_url(dataset_id, modality, "parquet"),
                DATASET_PARQUET_DOWNLOAD_LINK,
            ),
            _dataset_download_link(
                "Download full data (CSV.gz)",
                _dataset_download_url(dataset_id, modality),
                DATASET_DOWNLOAD_LINK,
            ),
        ],
        className="d-flex align-items-center justify-content-center flex-wrap gap-2 mb-3",
    )


def _perturb_seq_header_filters(
    perturbation_value: str = "", effect_value: str = ""
) -> Dict[str, Any]:
    return {
        "perturbation": dbc.Input(
            id=PERTURB_SEARCH_PERTURBATION_GENE,
            type="text",
            placeholder="Filter by name or ENSG",
            value=perturbation_value,
            size="sm",
            debounce=True,
            style={"width": "200px"},
        ),
        "effect": dbc.Input(
            id=PERTURB_SEARCH_EFFECT_GENE,
            type="text",
            placeholder="Filter by name or ENSG",
            value=effect_value,
            size="sm",
            debounce=True,
            style={"width": "200px"},
        ),
    }


@callback(
    [
        Output("dataset-content", "children"),
        Output(DATASET_DATA_STORE, "data"),
    ],
    Input(DATASET_STORE, "data"),
)
def render_dataset(data: Optional[Dict[str, Any]]):
    """Render dataset details and initialize data store."""
    if not data or not data.get("dataset_id"):
        return (
            html.Div(
                "No dataset ID provided.",
                className="text-muted text-center",
            ),
            None,
        )

    dataset_id = data["dataset_id"]
    formatted_dataset_id = _format_dataset_id(dataset_id)
    dataset_data, error = fetch_dataset(dataset_id)

    if error:
        return (
            html.Div(
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
            ),
            None,
        )

    if not dataset_data:
        return (
            html.Div(
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
            ),
            None,
        )

    # Detect modality from dataset metadata
    modality = _get_modality_from_dataset(dataset_data)

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

    # Define primary fields shown by default, in display order.
    # study_uri is consumed by study_title link and hidden separately.
    PRIMARY_STANDALONE_ORDER = [
        "study_title",
        "first_author",
        "last_author",
        "study_year",
        "experiment_title",
        "experiment_summary",
        "data_modalities",
    ]
    PRIMARY_GROUPED = {"perturbation_type"}
    SKIP_STANDALONE = {
        "dataset_id",
        "max_ingested_at",
        "study_uri",
        # Rendered as a badge next to the title instead of in the field list.
        "perturb_seq_reprocessed",
    }

    study_uri = standalone_fields.get("study_uri")

    # Build primary fields in the specified order
    primary_elements = []
    for key in PRIMARY_STANDALONE_ORDER:
        value = standalone_fields.get(key)
        if value is None:
            continue
        if key == "study_title":
            # Render as clickable link using study_uri
            display_value = (
                html.A(
                    _format_value(value),
                    href=study_uri,
                    target="_blank",
                    rel="noopener noreferrer",
                    className="text-decoration-none",
                    style={"color": COLORS["primary"]},
                )
                if study_uri
                else html.Span(_format_value(value), className="text-muted")
            )
            primary_elements.append(
                html.Div(
                    [
                        html.Dt("Study Title:", className="col-sm-4 fw-semibold"),
                        html.Dd(display_value, className="col-sm-8"),
                    ],
                    className="row mb-3",
                )
            )
        elif key == "data_modalities":
            primary_elements.append(
                html.Div(
                    [
                        html.Dt("Data Modality:", className="col-sm-4 fw-semibold"),
                        html.Dd(
                            html.Span(_format_value(value), className="text-muted"),
                            className="col-sm-8",
                        ),
                    ],
                    className="row mb-3",
                )
            )
        else:
            primary_elements.append(_create_field_display(key, value))

    # Add primary grouped fields (perturbation_type)
    for base_name in PRIMARY_GROUPED:
        group = field_groups.get(base_name)
        if not group:
            continue
        labels = group.get("labels")
        ids = group.get("ids")
        if labels is not None:
            primary_elements.append(
                _create_field_display(f"{base_name}_labels", labels, ids)
            )
        elif ids is not None:
            primary_elements.append(
                _create_field_display(f"{base_name}_ids", ids, None)
            )

    # Build secondary fields (everything else)
    primary_standalone_set = set(PRIMARY_STANDALONE_ORDER)
    secondary_elements = []

    for key, value in sorted(standalone_fields.items()):
        if key in SKIP_STANDALONE or key in primary_standalone_set:
            continue
        if key == "associated_datasets":
            secondary_elements.append(_create_associated_datasets_display(value))
        elif key == "timepoints" and isinstance(value, list):
            formatted_timepoints = _format_timepoints(value)
            secondary_elements.append(_create_field_display(key, formatted_timepoints))
        else:
            secondary_elements.append(_create_field_display(key, value))

    for base_name in sorted(field_groups.keys()):
        if base_name in PRIMARY_GROUPED:
            continue
        group = field_groups[base_name]
        labels = group.get("labels")
        ids = group.get("ids")
        if labels is not None:
            secondary_elements.append(
                _create_field_display(f"{base_name}_labels", labels, ids)
            )
        elif ids is not None:
            secondary_elements.append(
                _create_field_display(f"{base_name}_ids", ids, None)
            )

    # Build metadata section: primary fields + collapsible secondary
    metadata_children = [html.Dl(primary_elements)]
    if secondary_elements:
        metadata_children.append(
            html.Details(
                [
                    html.Summary(
                        "Other Metadata Fields",
                        className="fw-semibold",
                        style={
                            "cursor": "pointer",
                            "color": COLORS["primary"],
                            "fontSize": "1.1rem",
                            "padding": "0.5rem 0",
                        },
                    ),
                    html.Dl(secondary_elements, className="mt-3"),
                ],
                className="mt-2",
            )
        )

    metadata_section = html.Div(
        metadata_children,
        className="mt-4",
        style={"maxWidth": "900px", "margin": "0 auto"},
    )
    metadata_download_control = html.Div(
        _dataset_download_link(
            "Download metadata (JSON)",
            _dataset_download_url(dataset_id, modality, "metadata"),
            DATASET_METADATA_DOWNLOAD_LINK,
        ),
        className="mt-3 d-flex justify-content-center",
        style={"maxWidth": "900px", "margin": "0 auto"},
    )

    # Build search controls based on modality
    is_perturb_seq = modality == "perturb-seq"
    is_crispr = modality == "crispr-screen"
    data_download_controls = (
        _dataset_download_controls(dataset_id, modality) if modality else html.Div()
    )

    if is_perturb_seq:
        search_controls = html.Div(
            [
                # Hidden CRISPR input for callback compatibility
                dbc.Input(
                    id=CRISPR_SEARCH_PERTURBATION_GENE,
                    type="text",
                    value="",
                    style={"display": "none"},
                ),
            ],
            className="d-flex align-items-center flex-wrap gap-2 mb-3",
        )
    elif is_crispr:
        search_controls = html.Div(
            [
                dbc.Input(
                    id=CRISPR_SEARCH_PERTURBATION_GENE,
                    type="text",
                    placeholder="Filter by name or ENSG",
                    value="",
                    size="sm",
                    debounce=True,
                    style={"width": "200px"},
                ),
                # Hidden Perturb-seq inputs for callback compatibility
                dbc.Input(
                    id=PERTURB_SEARCH_PERTURBATION_GENE,
                    type="text",
                    value="",
                    style={"display": "none"},
                ),
                dbc.Input(
                    id=PERTURB_SEARCH_EFFECT_GENE,
                    type="text",
                    value="",
                    style={"display": "none"},
                ),
            ],
            className="d-flex align-items-center flex-wrap gap-2 mb-3",
        )
    else:
        # MAVE has no gene search, but retains the common hidden inputs.
        search_controls = html.Div(
            [
                dbc.Input(
                    id=PERTURB_SEARCH_PERTURBATION_GENE,
                    type="text",
                    value="",
                ),
                dbc.Input(
                    id=PERTURB_SEARCH_EFFECT_GENE,
                    type="text",
                    value="",
                ),
                dbc.Input(
                    id=CRISPR_SEARCH_PERTURBATION_GENE,
                    type="text",
                    value="",
                ),
            ],
            style={"display": "none"},
        )

    # Build score interpretation banner for CRISPR modality
    score_interpretation_banner = None
    if is_crispr:
        score_interpretation = dataset_data.get("score_interpretation")
        if score_interpretation:
            score_interpretation_banner = html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className="bi bi-info-circle-fill me-2",
                                style={
                                    "fontSize": "1.1rem",
                                    "color": COLORS["primary"],
                                },
                            ),
                            html.Strong("Score Interpretation: ", className="me-1"),
                            html.Span(score_interpretation),
                        ],
                        className="d-flex align-items-center",
                    ),
                ],
                className="alert mb-3",
                style={
                    "backgroundColor": "#e8f5e9",
                    "borderColor": COLORS["primary"],
                    "borderLeft": f"4px solid {COLORS['primary']}",
                    "color": "#1b5e20",
                    "borderRadius": "6px",
                },
            )

    # Build data visualization section
    initial_data_content = (
        perturb_seq_table(
            [],
            section_id="dataset_perturb_seq",
            header_filters=_perturb_seq_header_filters(),
        )
        if is_perturb_seq
        else None
    )
    data_section = html.Div(
        [
            html.Hr(className="my-4"),
            html.H2(
                "Data",
                className="fw-bold mb-3 text-center",
                style={"color": COLORS["primary"]},
            ),
            score_interpretation_banner,
            data_download_controls,
            search_controls,
            dcc.Loading(
                html.Div(initial_data_content, id="dataset-data-content"),
                type="circle",
                color=COLORS["primary"],
            ),
        ],
        style={"maxWidth": "1200px", "margin": "0 auto"},
    )

    # Initialize data store for fetching
    initial_data_store = None
    if modality:
        initial_data_store = {
            "dataset_id": dataset_id,
            "modality": modality,
            "results": None,
            "current_page": 1,
            "current_offset": 0,
            "current_position_range": "1_20" if modality == "mave" else None,
            "has_more": False,
            "total_rows_count": None,
            "error": None,
            "needs_fetch": True,
            # Perturb-seq search filters
            "perturbation_gene_search": "",
            "effect_gene_search": "",
            # CRISPR search filters
            "crispr_perturbation_gene_search": "",
        }

    provenance_badge = reprocessed_badge(dataset_data.get("perturb_seq_reprocessed"))
    title_children = [
        html.H1(
            f"Dataset: {formatted_dataset_id}",
            className="display-5 fw-bold mb-0",
            style={"color": COLORS["primary"]},
        )
    ]
    if provenance_badge is not None:
        title_children.append(provenance_badge)
    header_section = html.Div(
        title_children,
        className="d-flex align-items-center justify-content-center gap-2 mb-4 flex-wrap",
    )

    # Hidden inputs for no-modality case (needed for callback)
    hidden_search_inputs = html.Div(
        [
            dbc.Input(
                id=PERTURB_SEARCH_PERTURBATION_GENE,
                type="text",
                value="",
            ),
            dbc.Input(
                id=PERTURB_SEARCH_EFFECT_GENE,
                type="text",
                value="",
            ),
            dbc.Input(
                id=CRISPR_SEARCH_PERTURBATION_GENE,
                type="text",
                value="",
            ),
        ],
        style={"display": "none"},
    )

    content = html.Div(
        [
            header_section,
            metadata_download_control,
            metadata_section,
            (
                data_section
                if modality
                else html.Div(
                    [
                        html.Hr(className="my-4"),
                        # Include hidden inputs when no modality
                        hidden_search_inputs,
                        html.Div(
                            "No data visualization available for this dataset.",
                            className="text-muted text-center py-4",
                        ),
                    ],
                    style={"maxWidth": "900px", "margin": "0 auto"},
                )
            ),
        ]
    )

    return content, initial_data_store


def _fetch_data_rows(store_data: Dict[str, Any]) -> Dict[str, Any]:
    """Fetch data rows from the API."""
    dataset_id = store_data["dataset_id"]
    modality = store_data["modality"]

    # Build filters
    filters: Dict[str, Any] = {}

    if modality == "mave":
        # MAVE uses position range filter
        position_range = store_data.get("current_position_range", "1_20")
        filters["perturbation_position"] = position_range
        filters["effect_score_name"] = "score"

        response = fetch_dataset_rows(
            modality,
            dataset_id,
            filters=filters,
            limit=None,
            offset=None,
        )
    else:
        # CRISPR and Perturb-seq use limit/offset
        offset = store_data.get("current_offset", 0)

        # Add Perturb-seq search filters
        if modality == "perturb-seq":
            perturbation_gene = store_data.get("perturbation_gene_search", "").strip()
            effect_gene = store_data.get("effect_gene_search", "").strip()
            if perturbation_gene:
                filters["perturbation_gene_name"] = perturbation_gene
            if effect_gene:
                filters["effect_gene_name"] = effect_gene

        # Add CRISPR search filters
        if modality == "crispr-screen":
            crispr_perturbation_gene = store_data.get(
                "crispr_perturbation_gene_search", ""
            ).strip()
            if crispr_perturbation_gene:
                filters["perturbation_gene_name"] = crispr_perturbation_gene

        response = fetch_dataset_rows(
            modality,
            dataset_id,
            filters=filters,
            limit=ROWS_PER_PAGE,
            offset=offset,
        )

    # Update store with results
    results = response.get("results") or []
    total_count = response.get("total_rows_count")

    store_data["results"] = results
    store_data["error"] = response.get("error")
    store_data["total_rows_count"] = total_count

    # Calculate has_more
    if modality == "mave":
        store_data["has_more"] = len(results) > 0
    else:
        if total_count is not None:
            current_offset = store_data.get("current_offset", 0)
            store_data["has_more"] = current_offset + len(results) < total_count
        else:
            store_data["has_more"] = len(results) >= ROWS_PER_PAGE

    return store_data


def _paginate_data(store_data: Dict[str, Any], direction: str) -> Dict[str, Any]:
    """Update pagination state based on direction."""
    modality = store_data.get("modality")
    current_page = store_data.get("current_page", 1)

    if modality == "mave":
        # MAVE pagination by position range
        current_range = store_data.get("current_position_range", "1_20")
        try:
            start, end = map(int, current_range.split("_"))
        except (ValueError, AttributeError):
            start, end = 1, 20

        if direction == "next":
            new_start = end + 1
            new_end = new_start + MAVE_POSITION_RANGE_SIZE - 1
            store_data["current_page"] = current_page + 1
        elif direction == "previous":
            new_end = start - 1
            new_start = max(1, new_end - MAVE_POSITION_RANGE_SIZE + 1)
            store_data["current_page"] = max(1, current_page - 1)
        else:
            return store_data

        store_data["current_position_range"] = f"{new_start}_{new_end}"
    else:
        # CRISPR/Perturb-seq pagination by offset
        current_offset = store_data.get("current_offset", 0)

        if direction == "next":
            store_data["current_offset"] = current_offset + ROWS_PER_PAGE
            store_data["current_page"] = current_page + 1
        elif direction == "previous":
            store_data["current_offset"] = max(0, current_offset - ROWS_PER_PAGE)
            store_data["current_page"] = max(1, current_page - 1)

    return store_data


def _build_pagination_controls(store_data: Dict[str, Any]) -> html.Div:
    """Build pagination controls for data rows."""
    current_page = store_data.get("current_page", 1)
    has_more = store_data.get("has_more", False)
    has_previous = current_page > 1
    total_count = store_data.get("total_rows_count")
    modality = store_data.get("modality")

    if not has_previous and not has_more:
        return html.Div()  # No pagination needed

    buttons = []

    if has_previous:
        buttons.append(
            dbc.Button(
                "← Previous",
                id={"type": "dataset-paginate", "direction": "previous"},
                color="secondary",
                outline=True,
                size="sm",
                className="me-2",
            )
        )

    # Page info
    if total_count and modality != "mave":
        total_pages = (total_count + ROWS_PER_PAGE - 1) // ROWS_PER_PAGE
        page_text = f"Page {current_page} of {total_pages}"
    else:
        page_text = f"Page {current_page}"

    buttons.append(
        html.Span(
            page_text,
            className="align-self-center mx-2",
            style={"fontSize": "0.875rem"},
        )
    )

    if has_more:
        buttons.append(
            dbc.Button(
                "Next →",
                id={"type": "dataset-paginate", "direction": "next"},
                color="secondary",
                outline=True,
                size="sm",
            )
        )

    return html.Div(
        buttons,
        className="d-flex align-items-center justify-content-end mt-3 w-100",
    )


def _render_data_visualization(store_data: Dict[str, Any]) -> html.Div:
    """Render the appropriate visualization based on modality."""
    modality = store_data.get("modality")
    results = store_data.get("results") or []
    error = store_data.get("error")
    dataset_id = store_data.get("dataset_id")

    if error:
        return html.Div(
            f"Error loading data: {error}",
            className="alert alert-danger",
        )

    if not results and modality != "perturb-seq":
        return html.Div(
            "No data rows found for this dataset.",
            className="text-muted fst-italic text-center py-4",
        )

    # Render based on modality
    if modality == "mave":
        table_component = mave_heatmap(results, download_url=None)
    elif modality == "perturb-seq":
        # Download controls are in the main layout for perturb-seq.
        filters = {
            "perturbation_gene_name": (
                store_data.get("perturbation_gene_search") or ""
            ).strip(),
            "effect_gene_name": (store_data.get("effect_gene_search") or "").strip(),
        }
        header_download = None
        if any(filters.values()):
            header_download = _dataset_download_link(
                "Download filtered data (CSV)",
                _dataset_download_url(dataset_id, modality, "csv.gz", filters),
                DATASET_FILTERED_DOWNLOAD_LINK,
            )
        table_component = perturb_seq_table(
            results,
            section_id="dataset_perturb_seq",
            download_url=None,
            header_filters=_perturb_seq_header_filters(
                store_data.get("perturbation_gene_search", ""),
                store_data.get("effect_gene_search", ""),
            ),
            header_download=header_download,
        )
    elif modality == "crispr-screen":
        # Download button is in the main layout for crispr-screen
        table_component = crispr_table(results, download_url=None)
    else:
        table_component = html.Div(
            f"Unknown modality: {modality}",
            className="text-muted fst-italic",
        )

    # Build pagination controls
    pagination = _build_pagination_controls(store_data)

    return html.Div([table_component, pagination])


@callback(
    [
        Output("dataset-data-content", "children"),
        Output(DATASET_DATA_STORE, "data", allow_duplicate=True),
    ],
    [
        Input(DATASET_DATA_STORE, "data"),
        Input({"type": "dataset-paginate", "direction": ALL}, "n_clicks"),
        Input(PERTURB_SEARCH_PERTURBATION_GENE, "value"),
        Input(PERTURB_SEARCH_EFFECT_GENE, "value"),
        Input(CRISPR_SEARCH_PERTURBATION_GENE, "value"),
    ],
    prevent_initial_call=True,
)
def render_dataset_data(
    store_data: Optional[Dict[str, Any]],
    _paginate_clicks,
    perturbation_gene_search: Optional[str],
    effect_gene_search: Optional[str],
    crispr_perturbation_gene_search: Optional[str],
):
    """Fetch data rows and render visualization."""
    if not store_data:
        return (
            html.Div(
                "No modality data available for this dataset.",
                className="text-muted fst-italic text-center py-4",
            ),
            store_data,
        )

    dataset_id = store_data.get("dataset_id")
    modality = store_data.get("modality")

    if not dataset_id or not modality:
        return (
            html.Div(
                "Unable to determine data modality.",
                className="text-muted fst-italic text-center py-4",
            ),
            store_data,
        )

    # Determine if we need to fetch data
    ctx = dash.callback_context
    triggered_id = ctx.triggered_id if hasattr(ctx, "triggered_id") else None
    needs_fetch = store_data.get("needs_fetch", False)

    # Handle pagination
    if (
        isinstance(triggered_id, dict)
        and triggered_id.get("type") == "dataset-paginate"
    ):
        direction = triggered_id.get("direction")
        if direction:
            store_data = _paginate_data(store_data, direction)
            needs_fetch = True

    # Handle Perturb-seq search inputs
    if modality == "perturb-seq":
        old_perturbation_search = store_data.get("perturbation_gene_search", "")
        old_effect_search = store_data.get("effect_gene_search", "")
        new_perturbation_search = perturbation_gene_search or ""
        new_effect_search = effect_gene_search or ""

        # Check if search terms changed
        if (
            new_perturbation_search != old_perturbation_search
            or new_effect_search != old_effect_search
        ):
            store_data["perturbation_gene_search"] = new_perturbation_search
            store_data["effect_gene_search"] = new_effect_search
            # Reset pagination when search changes
            store_data["current_page"] = 1
            store_data["current_offset"] = 0
            needs_fetch = True

    # Handle CRISPR search inputs
    if modality == "crispr-screen":
        old_crispr_search = store_data.get("crispr_perturbation_gene_search", "")
        new_crispr_search = crispr_perturbation_gene_search or ""

        # Check if search term changed
        if new_crispr_search != old_crispr_search:
            store_data["crispr_perturbation_gene_search"] = new_crispr_search
            # Reset pagination when search changes
            store_data["current_page"] = 1
            store_data["current_offset"] = 0
            needs_fetch = True

    # Fetch data if needed
    if needs_fetch or store_data.get("results") is None:
        store_data = _fetch_data_rows(store_data)
        store_data["needs_fetch"] = False

    # Render the visualization
    content = _render_data_visualization(store_data)

    return content, store_data


@callback(
    [
        Output(DATASET_METADATA_DOWNLOAD_LINK, "href"),
        Output(DATASET_PARQUET_DOWNLOAD_LINK, "href"),
        Output(DATASET_DOWNLOAD_LINK, "href"),
    ],
    Input(DATASET_DATA_STORE, "data"),
)
def update_dataset_download_link(store_data: Optional[Dict[str, Any]]):
    """Point full downloads at release artifacts."""
    if not store_data:
        return "#", "#", "#"

    dataset_id = store_data.get("dataset_id")
    modality = store_data.get("modality")
    if not dataset_id or not modality:
        return "#", "#", "#"

    return (
        _dataset_download_url(dataset_id, modality, "metadata"),
        _dataset_download_url(dataset_id, modality, "parquet"),
        _dataset_download_url(dataset_id, modality, "csv.gz"),
    )
