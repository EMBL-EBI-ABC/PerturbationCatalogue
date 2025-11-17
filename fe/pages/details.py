"""New Target Details page implementing the reusable data table component."""

from __future__ import annotations

import copy
import json
from typing import Any, Dict, List, Optional

import dash
from dash import ALL, MATCH, Input, Output, State, callback, dcc, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from components.target_data_table import TargetDataTable
from utils import (
    COLORS,
    fetch_dataset_rows,
    fetch_modality_datasets,
)


dash.register_page(
    __name__,
    path_template="/target/<target_name>",
    name="Target Details (Grid)",
    title="Target Details",
)

TARGET_NAME_STORE = "target-details-v2-name"
SECTION_CONFIGS = [
    {
        "id": "crispr",
        "title": "CRISPR screen data",
        "modality": "crispr-screen",
        "filter_field": "perturbation_gene_name",
    },
    {
        "id": "mave",
        "title": "MAVE data",
        "modality": "mave",
        "filter_field": "perturbation_gene_name",
    },
    {
        "id": "perturb_seq_perturbed",
        "title": "Perturb-Seq (Perturbed)",
        "modality": "perturb-seq",
        "filter_field": "perturbation_gene_name",
    },
    {
        "id": "perturb_seq_affected",
        "title": "Perturb-Seq (Affected)",
        "modality": "perturb-seq",
        "filter_field": "effect_gene_name",
    },
]
SECTION_LOOKUP = {config["id"]: config for config in SECTION_CONFIGS}
MODALITY_BUTTON_GROUP = "target-modality-buttons"
SCROLL_VISIBILITY_STORE = "target-scroll-visible"
SCROLL_TOP_STYLE = {
    "position": "fixed",
    "bottom": "1.5rem",
    "right": "1.5rem",
    "zIndex": 1050,
    "display": "none",
}

DATASET_LIMIT = 4
ROWS_PER_DATASET_LIMIT = 5
DATASET_LOAD_MORE_SIZE = 5


def layout(target_name: Optional[str] = None, **kwargs):
    """Page layout with data stores and section shells."""
    stores = [
        dcc.Store(id=TARGET_NAME_STORE, data=target_name),
        dcc.Store(id=SCROLL_VISIBILITY_STORE, data=False),
    ]
    stores.extend(
        dcc.Store(
            id={"type": "target-section-store", "section": config["id"]},
            data={"section": config["id"], "dataset_search": ""},
        )
        for config in SECTION_CONFIGS
    )
    stores.extend(
        dcc.Store(
            id={"type": "section-collapse-store", "section": config["id"]},
            data=False,  # False = collapsed, True = expanded
        )
        for config in SECTION_CONFIGS
    )

    heading = target_name or "Target"

    sections = []
    for config in SECTION_CONFIGS:
        sections.append(
            html.Section(
                [
                    html.Div(
                        [
                            html.H3(
                                config["title"],
                                className="fw-semibold",
                                style={"margin": 0, "cursor": "pointer"},
                            ),
                            html.Span(
                                "▶",
                                id={"type": "section-arrow", "section": config["id"]},
                                className="ms-2",
                                style={
                                    "cursor": "pointer",
                                    "fontSize": "0.8em",
                                    "transition": "transform 0.2s",
                                    "display": "inline-block",
                                },
                            ),
                        ],
                        id={"type": "section-header", "section": config["id"]},
                        className="d-flex align-items-center mb-3",
                    ),
                    html.Div(
                        [
                            dcc.Input(
                                id={"type": "dataset-search", "section": config["id"]},
                                type="text",
                                placeholder="Search datasets…",
                                debounce=True,
                                className="form-control form-control-sm",
                                style={
                                    "maxWidth": "240px",
                                    "position": "absolute",
                                    "top": "2.5rem",
                                    "left": "0.75rem",
                                    "zIndex": 5,
                                    "marginLeft": "15rem",
                                },
                            ),
                            dcc.Loading(
                                html.Div(
                                    id={
                                        "type": "target-table-container",
                                        "section": config["id"],
                                    }
                                ),
                                type="default",
                                color=COLORS["primary"],
                            ),
                        ],
                        id={"type": "section-content", "section": config["id"]},
                        className="position-relative",
                        style={"paddingTop": "2.5rem", "display": "none"},
                    ),
                ],
                id=f"section-{config['id']}",
                className="mb-5",
            )
        )

    return dbc.Container(
        stores
        + [
            html.Div(id="page-top"),
            html.Div(
                [
                    html.H1(
                        f"Target: {heading}",
                        className="display-5 fw-bold mb-2 text-center",
                        style={"color": COLORS["primary"]},
                    ),
                    html.P(
                        "Explore perturbation datasets across modalities.",
                        className="text-muted text-center mb-4",
                    ),
                    html.Div(
                        id=MODALITY_BUTTON_GROUP,
                        className="d-flex flex-wrap gap-3 justify-content-center mb-4",
                    ),
                ]
            ),
            *sections,
            dbc.Button(
                "↑ Back to top",
                id="target-scroll-top",
                href="#page-top",
                color="primary",
                className="position-fixed shadow",
                style=SCROLL_TOP_STYLE,
            ),
        ],
        fluid=True,
        className="py-4 target-details-v2",
    )


def _get_triggered_id():
    """Return the parsed triggered ID for Dash versions with/without triggered_id."""
    ctx = dash.callback_context
    triggered = getattr(ctx, "triggered_id", None)
    if triggered is not None:
        return triggered
    if not ctx.triggered:
        return None
    prop_id = ctx.triggered[0]["prop_id"].split(".")[0]
    try:
        return json.loads(prop_id)
    except json.JSONDecodeError:
        return prop_id


def _empty_store(config: Dict[str, Any], error: Optional[str] = None) -> Dict[str, Any]:
    return {
        "section": config["id"],
        "modality": config["modality"],
        "title": config["title"],
        "datasets": [],
        "total_datasets_count": 0,
        "filters": {},
        "rows_per_dataset_limit": ROWS_PER_DATASET_LIMIT,
        "error": error,
        "target_name": None,
        "dataset_search": "",
    }


@callback(
    Output({"type": "target-table-container", "section": MATCH}, "children"),
    Input({"type": "target-section-store", "section": MATCH}, "data"),
)
def render_section(store_data: Optional[Dict[str, Any]]):
    """Render each section's summary plus data table."""
    triggered = _get_triggered_id()
    section_id = triggered.get("section") if isinstance(triggered, dict) else None

    if not store_data:
        return html.Div("No data available.", className="text-muted fst-italic")

    datasets = store_data.get("datasets") or []
    total_datasets = store_data.get("total_datasets_count", len(datasets))
    rows_limit = store_data.get("rows_per_dataset_limit", ROWS_PER_DATASET_LIMIT)

    summary = html.Div(
        [
            html.Span(
                f"{len(datasets)} dataset(s) shown (of {total_datasets})",
                className="fw-semibold text-muted",
            ),
            html.Span(
                f"Top {rows_limit} rows per dataset shown by default",
                className="text-muted small",
            ),
        ],
        className="d-flex flex-column flex-md-row gap-2 justify-content-between mb-3",
    )

    def control_factory(dataset_id: str, dataset_entry: Dict[str, Any]):
        if not dataset_id or not section_id:
            return None

        current_page = dataset_entry.get("current_page", 1)
        total_count = dataset_entry.get("total_rows_count")
        has_next = dataset_entry.get("has_more", False)
        has_previous = current_page > 1

        # Create page info text
        if total_count is not None and total_count > 0:
            total_pages = (
                total_count + DATASET_LOAD_MORE_SIZE - 1
            ) // DATASET_LOAD_MORE_SIZE
            page_text = f"Page {current_page} of {total_pages}"
        else:
            page_text = f"Page {current_page}"

        buttons = []

        # Previous button
        if has_previous:
            buttons.append(
                dbc.Button(
                    "← Previous",
                    id={
                        "type": "target-table-paginate",
                        "section": section_id,
                        "dataset": dataset_id,
                        "direction": "previous",
                    },
                    color="secondary",
                    outline=True,
                    size="sm",
                    className="me-2",
                )
            )

        # Page info
        buttons.append(
            html.Span(
                page_text,
                className="align-self-center me-2",
                style={"fontSize": "0.875rem"},
            )
        )

        # Next button
        if has_next:
            buttons.append(
                dbc.Button(
                    "Next →",
                    id={
                        "type": "target-table-paginate",
                        "section": section_id,
                        "dataset": dataset_id,
                        "direction": "next",
                    },
                    color="secondary",
                    outline=True,
                    size="sm",
                )
            )

        if not buttons:
            return None

        return html.Div(
            buttons,
            className="d-flex align-items-center",
        )

    section_id = store_data.get("section")

    table = TargetDataTable(
        data=datasets,
        modality=store_data.get("modality", ""),
        rows_per_dataset_limit=rows_limit,
        error_message=store_data.get("error"),
        dataset_control_factory=control_factory,
        effect_gene_source=(
            "perturbation" if section_id == "perturb_seq_affected" else "effect"
        ),
        section_id=section_id,
    )

    return html.Div([summary, table])


@callback(
    Output(MODALITY_BUTTON_GROUP, "children"),
    [
        Input({"type": "target-section-store", "section": config["id"]}, "data")
        for config in SECTION_CONFIGS
    ],
)
def render_modality_buttons(*stores):
    """Display per-modality dataset counts as anchor buttons."""
    buttons = []
    color_map = {
        "crispr-screen": "secondary",
        "mave": "danger",
        "perturb-seq": "success",
    }

    for store, config in zip(stores, SECTION_CONFIGS):
        total = store.get("total_datasets_count") if store else None
        count_text = "…" if total is None else total
        buttons.append(
            html.A(
                dbc.Button(
                    f"{config['title']} ({count_text})",
                    id={"type": "modality-nav-button", "section": config["id"]},
                    color=color_map.get(config["modality"], "primary"),
                    outline=True,
                    className="px-3 py-2 fw-semibold",
                    style={"borderWidth": "2px"},
                ),
                href=f"#section-{config['id']}",
                className="text-decoration-none",
            )
        )

    if not buttons:
        return html.Div("Loading modality counts...", className="text-muted fst-italic")

    return buttons


@callback(
    Output(SCROLL_VISIBILITY_STORE, "data"),
    [
        Input({"type": "modality-nav-button", "section": config["id"]}, "n_clicks")
        for config in SECTION_CONFIGS
    ],
    State(SCROLL_VISIBILITY_STORE, "data"),
    prevent_initial_call=True,
)
def mark_scroll_visible(*args):
    *clicks, current = args
    if any(click and click > 0 for click in clicks):
        return True
    return current


@callback(
    [
        Output(
            {"type": "section-collapse-store", "section": config["id"]},
            "data",
            allow_duplicate=True,
        )
        for config in SECTION_CONFIGS
    ],
    [
        Input({"type": "modality-nav-button", "section": config["id"]}, "n_clicks")
        for config in SECTION_CONFIGS
    ],
    [
        State({"type": "section-collapse-store", "section": config["id"]}, "data")
        for config in SECTION_CONFIGS
    ],
    prevent_initial_call=True,
)
def expand_section_from_button(*args):
    """Expand the section when its modality button is clicked."""
    *clicks_and_states, _ = args
    num_sections = len(SECTION_CONFIGS)
    clicks = clicks_and_states[:num_sections]
    states = clicks_and_states[num_sections:]

    triggered = _get_triggered_id()

    # Determine which section was clicked
    clicked_section_id = None
    if isinstance(triggered, dict) and triggered.get("type") == "modality-nav-button":
        clicked_section_id = triggered.get("section")

    # Update all section states
    results = []
    for i, config in enumerate(SECTION_CONFIGS):
        if config["id"] == clicked_section_id and clicks[i] and clicks[i] > 0:
            # Expand the clicked section
            results.append(True)
        else:
            # Keep current state
            results.append(dash.no_update)

    return results


@callback(
    Output("target-scroll-top", "style"),
    Input(SCROLL_VISIBILITY_STORE, "data"),
)
def toggle_scroll_button(visible: bool):
    style = SCROLL_TOP_STYLE.copy()
    if visible:
        style["display"] = "block"
    return style


@callback(
    Output({"type": "section-collapse-store", "section": MATCH}, "data"),
    Input({"type": "section-header", "section": MATCH}, "n_clicks"),
    Input({"type": "section-arrow", "section": MATCH}, "n_clicks"),
    State({"type": "section-collapse-store", "section": MATCH}, "data"),
    prevent_initial_call=True,
)
def toggle_section_collapse(header_clicks, arrow_clicks, is_expanded):
    """Toggle section collapse state when header or arrow is clicked."""
    if header_clicks or arrow_clicks:
        return not is_expanded
    return is_expanded


@callback(
    [
        Output({"type": "section-content", "section": MATCH}, "style"),
        Output({"type": "section-arrow", "section": MATCH}, "style"),
    ],
    Input({"type": "section-collapse-store", "section": MATCH}, "data"),
)
def update_section_visibility(is_expanded):
    """Update section content visibility and arrow rotation based on collapse state."""
    if is_expanded:
        content_style = {"paddingTop": "2.5rem", "display": "block"}
        arrow_style = {
            "cursor": "pointer",
            "fontSize": "0.8em",
            "transition": "transform 0.2s",
            "transform": "rotate(90deg)",
            "display": "inline-block",
        }
    else:
        content_style = {"paddingTop": "2.5rem", "display": "none"}
        arrow_style = {
            "cursor": "pointer",
            "fontSize": "0.8em",
            "transition": "transform 0.2s",
            "transform": "rotate(0deg)",
            "display": "inline-block",
        }
    return content_style, arrow_style


@callback(
    Output({"type": "target-section-store", "section": MATCH}, "data"),
    Input(TARGET_NAME_STORE, "data"),
    Input(
        {
            "type": "target-table-paginate",
            "section": MATCH,
            "dataset": ALL,
            "direction": ALL,
        },
        "n_clicks",
    ),
    Input({"type": "dataset-search", "section": MATCH}, "value"),
    State({"type": "target-section-store", "section": MATCH}, "data"),
    State({"type": "target-section-store", "section": MATCH}, "id"),
)
def update_section_store(
    target_name: Optional[str],
    _paginate_clicks,
    search_value: Optional[str],
    store_data: Optional[Dict[str, Any]],
    store_id: Optional[Dict[str, Any]],
):
    """Hydrate or paginate rows for a section store."""
    if not store_id:
        raise PreventUpdate

    section_id = store_id.get("section")
    config = SECTION_LOOKUP.get(section_id)
    if not config:
        raise PreventUpdate

    if not target_name:
        return _empty_store(config, "No target specified.")

    triggered = _get_triggered_id()

    search_changed = False
    previous_query = ""
    if store_data is not None:
        previous_query = store_data.get("dataset_search", "")
        search_changed = (search_value or "") != (previous_query or "")

    should_refresh = (
        triggered == TARGET_NAME_STORE
        or not store_data
        or store_data.get("target_name") != target_name
    )

    if isinstance(triggered, dict) and triggered.get("type") == "dataset-search":
        should_refresh = True

    if should_refresh or search_changed:
        current_search = (search_value or "") if search_changed else previous_query
        return _fetch_section_payload(config, target_name, current_search)

    if isinstance(triggered, dict) and triggered.get("type") == "target-table-paginate":
        dataset_id = triggered.get("dataset")
        direction = triggered.get("direction")
        if not dataset_id or not direction:
            raise PreventUpdate
        return _paginate_dataset_rows(store_data, dataset_id, direction)

    return store_data


def _fetch_section_payload(
    config: Dict[str, Any],
    target_name: str,
    dataset_search: Optional[str] = None,
) -> Dict[str, Any]:
    filters = {config["filter_field"]: target_name}
    if dataset_search:
        cleaned = dataset_search.strip()
        if cleaned:
            filters["dataset_metadata"] = cleaned
    response = fetch_modality_datasets(
        config["modality"],
        filters=filters,
        dataset_limit=DATASET_LIMIT,
        rows_per_dataset_limit=ROWS_PER_DATASET_LIMIT,
    )
    datasets = response.get("datasets") or []
    for dataset in datasets:
        results = dataset.get("results") or []
        dataset["rows_per_dataset_limit"] = ROWS_PER_DATASET_LIMIT
        dataset["truncated"] = len(results) == ROWS_PER_DATASET_LIMIT
        # Initialize pagination state
        dataset["current_page"] = 1
        dataset["current_offset"] = 0
        # Get total count from API response, or infer from results
        total_count = dataset.get("total_rows_count")
        if total_count is None:
            # If not provided, assume there are more if we got the full limit
            total_count = (
                len(results) if len(results) < ROWS_PER_DATASET_LIMIT else None
            )
        dataset["total_rows_count"] = total_count
        # Determine if there are more rows: either total > current, or we got full limit and total is unknown
        if total_count is not None:
            dataset["has_more"] = len(results) < total_count
        else:
            dataset["has_more"] = len(results) >= ROWS_PER_DATASET_LIMIT

    return {
        "section": config["id"],
        "modality": config["modality"],
        "title": config["title"],
        "datasets": datasets,
        "total_datasets_count": response.get("total_datasets_count", 0),
        "filters": filters,
        "rows_per_dataset_limit": ROWS_PER_DATASET_LIMIT,
        "error": response.get("error"),
        "target_name": target_name,
        "dataset_search": dataset_search or "",
    }


def _paginate_dataset_rows(
    store_data: Dict[str, Any], dataset_id: str, direction: str
) -> Dict[str, Any]:
    """Paginate rows for a specific dataset by calling the API with pagination."""
    updated_store = copy.deepcopy(store_data)
    datasets = updated_store.get("datasets") or []

    # Find the target dataset - check both dataset_id and id (fallback) in nested dataset object
    target_dataset = None
    for entry in datasets:
        dataset_meta = entry.get("dataset") or {}
        entry_dataset_id = dataset_meta.get("dataset_id") or dataset_meta.get("id")
        if entry_dataset_id == dataset_id:
            target_dataset = entry
            break

    if not target_dataset:
        # Debug: log if dataset not found
        print(
            f"Warning: Dataset {dataset_id} not found in datasets. Available dataset IDs: {[((e.get('dataset') or {}).get('dataset_id') or (e.get('dataset') or {}).get('id')) for e in datasets]}"
        )
        return updated_store

    # Get current pagination state
    current_page = target_dataset.get("current_page", 1)
    current_offset = target_dataset.get("current_offset", 0)

    # Calculate new offset based on direction
    if direction == "next":
        new_offset = current_offset + DATASET_LOAD_MORE_SIZE
        new_page = current_page + 1
    elif direction == "previous":
        new_offset = max(0, current_offset - DATASET_LOAD_MORE_SIZE)
        new_page = max(1, current_page - 1)
    else:
        return updated_store

    # Get the target symbol (perturbed target) for the API call
    target_symbol = updated_store.get("target_name")

    # Get the correct filter field based on the section
    section_id = updated_store.get("section")
    config = SECTION_LOOKUP.get(section_id) if section_id else None
    filter_field = (
        config.get("filter_field", "perturbation_gene_name")
        if config
        else "perturbation_gene_name"
    )

    # Call API: /v1/{modality}/{dataset_id}/search?{filter_field}={target_symbol}&limit=5&offset=X
    response = fetch_dataset_rows(
        updated_store.get("modality", ""),
        dataset_id,
        filters={filter_field: target_symbol} if target_symbol else {},
        offset=new_offset,
        limit=DATASET_LOAD_MORE_SIZE,  # Always 5 rows per page
    )

    # Replace results instead of appending
    new_results = response.get("results") or []

    # Update total count if provided by API
    api_total_count = response.get("total_rows_count")
    if api_total_count is not None:
        target_dataset["total_rows_count"] = api_total_count
        total_count = api_total_count
    else:
        # Keep existing total_count or None if not set
        total_count = target_dataset.get("total_rows_count")

    # Update dataset state
    target_dataset["results"] = new_results
    target_dataset["current_page"] = new_page
    target_dataset["current_offset"] = new_offset
    target_dataset["truncated"] = False

    # Determine if there are more rows available
    if total_count is not None:
        target_dataset["has_more"] = new_offset + len(new_results) < total_count
    else:
        # If total is unknown, assume there are more if we got a full page
        target_dataset["has_more"] = len(new_results) >= DATASET_LOAD_MORE_SIZE

    if response.get("error"):
        updated_store["error"] = response["error"]

    return updated_store
