"""Shared search table components extracted from home.py."""

import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from urllib.parse import quote
from utils import (
    COLORS,
    DATA_MODALITIES_COLOURS,
    FACET_FIELDS,
    results_store,
    format_value,
    reprocessed_badge,
)

SEARCH_RESULTS_PAGE_SIZE = 15


def format_count(value):
    """Format numeric counts for display, falling back to '0' when missing."""
    if value in (None, "", "N/A"):
        return "0"
    if isinstance(value, (int, float)):
        return f"{int(value):,}" if isinstance(value, int) else format_value(value)
    return format_value(value)


def render_targets_table(results):
    """Render the search results table for target entries."""
    results_store.clear()

    if not results:
        return dbc.Alert(
            [
                html.I(className="bi bi-info-circle me-2"),
                "No results found. Try another target name.",
            ],
            color="info",
            className="shadow-sm",
            style={"borderRadius": "10px"},
        )

    header = html.Thead(
        html.Tr(
            [
                html.Th(
                    [
                        "Target name ",
                        html.Span(
                            html.I(className="bi bi-question-circle"),
                            id="target-info-icon",
                            style={"cursor": "pointer", "color": "#6c757d"},
                        ),
                        dbc.Popover(
                            [
                                dbc.PopoverHeader("Target Gene Aggregation"),
                                dbc.PopoverBody(
                                    "Results are aggregated by target gene name. Each row represents all datasets and modalities available for a specific target, combining data from Perturb-seq, CRISPR screens, and MAVE experiments across different tissues, cell types, and experimental conditions."
                                ),
                            ],
                            target="target-info-icon",
                            trigger="click",
                            placement="bottom",
                        ),
                    ],
                    className="fw-semibold",
                ),
                html.Th("Perturb-seq datasets", className="fw-semibold text-center"),
                html.Th(
                    "CRISPR-screen datasets",
                    className="fw-semibold text-center",
                ),
                html.Th(
                    "MAVE datasets",
                    className="fw-semibold text-center",
                ),
                html.Th(
                    [
                        "Top GSEA Terms ",
                        html.Span(
                            html.I(className="bi bi-question-circle"),
                            id="gsea-info-icon",
                            style={"cursor": "pointer", "color": "#6c757d"},
                        ),
                        dbc.Popover(
                            [
                                dbc.PopoverHeader("Pathway enrichment (GSEA)"),
                                dbc.PopoverBody(
                                    "Shows biological pathways whose genes are collectively up- or down-regulated after a genetic perturbation, based on single-cell Perturb-seq data and MSigDB Hallmark gene sets."
                                ),
                            ],
                            target="gsea-info-icon",
                            trigger="click",
                            placement="bottom",
                        ),
                    ],
                    className="fw-semibold",
                ),
                html.Th("Data Modalities", className="fw-semibold"),
            ],
            style={"backgroundColor": "#f1f3f5"},
        )
    )

    rows = []
    for record in results:
        target_id = record.get("ensembl_gene_id") or "N/A"
        approved_symbol = record.get("approved_symbol") or target_id
        approved_name = record.get("approved_name")
        results_store[target_id] = record

        n_sc_perturb_seq = format_count(record.get("n_perturb_seq"))
        n_sc_perturb_seq_up = format_count(record.get("n_sig_perturb_pairs_up"))
        n_sc_perturb_seq_down = format_count(record.get("n_sig_perturb_pairs_down"))
        n_crispr = format_count(record.get("n_crispr"))
        n_sig_crispr = format_count(record.get("n_sig_crispr"))
        n_mave = format_count(record.get("n_mave"))
        top_gsea_terms = record.get("top_gsea_terms") or []
        data_modalities = record.get("data_modalities") or []

        gsea_badges = []
        for term in top_gsea_terms[:5]:
            display_term = term if len(term) <= 25 else term[:22] + "..."
            gsea_badges.append(
                html.Span(
                    dbc.Badge(
                        display_term,
                        color="light",
                        className="me-1 mb-1",
                        style={
                            "fontSize": "0.65rem",
                            "border": "1px solid #6c757d",
                            "color": "#495057",
                            "backgroundColor": "#f8f9fa",
                            "whiteSpace": "nowrap",
                            "overflow": "hidden",
                            "textOverflow": "ellipsis",
                            "maxWidth": "150px",
                            "display": "inline-block",
                        },
                    ),
                    title=term,
                )
            )
        if not gsea_badges:
            gsea_badges = [
                html.Span(
                    "No significant hits were found",
                    className="text-muted",
                    style={"fontSize": "0.8rem"},
                )
            ]

        modalities_badges = []
        for modality in data_modalities:
            badge_color = DATA_MODALITIES_COLOURS.get(modality, COLORS["primary"])
            modalities_badges.append(
                dbc.Badge(
                    modality,
                    color="light",
                    className="me-1 mb-1",
                    style={
                        "fontSize": "0.75rem",
                        "border": f"1px solid {badge_color}",
                        "color": badge_color,
                        "backgroundColor": f"{badge_color}1A",
                    },
                )
            )
        if not modalities_badges:
            modalities_badges = [html.Span("N/A", className="text-muted")]

        rows.append(
            html.Tr(
                [
                    html.Td(
                        dcc.Link(
                            html.Div(
                                [
                                    html.Div(
                                        approved_symbol,
                                        className="fw-bold",
                                        style={"color": COLORS["primary"]},
                                    ),
                                    html.Div(
                                        target_id,
                                        className="text-muted",
                                        style={"fontSize": "0.75rem"},
                                    ),
                                ],
                                title=approved_name or approved_symbol,
                            ),
                            href=f"/perturbation-catalogue/target/{quote(target_id, safe='')}",
                            className="text-decoration-none",
                        )
                    ),
                    html.Td(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        f"Datasets: {n_sc_perturb_seq} ",
                                        style={
                                            "marginRight": "5px",
                                        },
                                    ),
                                    html.Span(
                                        [
                                            html.Span("("),
                                            html.Span(
                                                "Sig. genes up \u2191: ",
                                                style={
                                                    "color": "#1f9d55",
                                                    "marginRight": "0.25rem",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(
                                                n_sc_perturb_seq_up,
                                                style={
                                                    "color": "#1f9d55",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(";"),
                                        ],
                                        className="d-inline-flex align-items-center me-3",
                                    ),
                                    html.Span(
                                        [
                                            html.Span(
                                                "Sig. genes down \u2193: ",
                                                style={
                                                    "color": "#c53030",
                                                    "marginRight": "0.25rem",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(
                                                n_sc_perturb_seq_down,
                                                style={
                                                    "color": "#c53030",
                                                    "fontWeight": "600",
                                                },
                                            ),
                                            html.Span(")"),
                                        ],
                                        className="d-inline-flex align-items-center",
                                    ),
                                ],
                                className="d-flex justify-content-center flex-wrap",
                            )
                        ],
                        className="text-center",
                    ),
                    html.Td(
                        html.Div(
                            [
                                html.Span(
                                    f"Datasets: {n_crispr}",
                                    style={
                                        "marginRight": "5px",
                                    },
                                ),
                                html.Span(
                                    f"Sig. hits: {n_sig_crispr}",
                                    style={
                                        "color": "#1f9d55",
                                        "marginRight": "0.25rem",
                                        "fontWeight": "600",
                                    },
                                ),
                            ],
                            className="d-flex justify-content-center flex-wrap",
                        ),
                        className="text-center",
                    ),
                    html.Td(n_mave, className="text-center"),
                    html.Td(
                        html.Div(gsea_badges, className="d-flex flex-wrap"),
                        style={
                            "maxWidth": "180px",
                            "minWidth": "120px",
                            "overflow": "hidden",
                        },
                    ),
                    html.Td(modalities_badges),
                ]
            )
        )

    table = dbc.Table(
        [header, html.Tbody(rows)],
        bordered=False,
        hover=True,
        responsive=True,
        striped=False,
        className="align-middle shadow-sm",
        style={"borderRadius": "12px", "overflow": "hidden"},
    )

    return table


def render_datasets_table(results):
    """Render the search results table for dataset search results."""
    if not results:
        return dbc.Alert(
            [
                html.I(className="bi bi-info-circle me-2"),
                "No datasets found. Try another search term.",
            ],
            color="info",
            className="shadow-sm",
            style={"borderRadius": "10px"},
        )

    header = html.Thead(
        html.Tr(
            [
                html.Th("Dataset ID", className="fw-semibold"),
                html.Th("Study Title", className="fw-semibold"),
                html.Th("Study Year", className="fw-semibold text-center"),
                html.Th("Data Modality", className="fw-semibold"),
                html.Th("Provenance", className="fw-semibold"),
            ],
            style={"backgroundColor": "#f1f3f5"},
        )
    )

    rows = []
    for record in results:
        dataset_id = record.get("dataset_id", "N/A")
        study_title = (
            record.get("study_title") or record.get("experiment_title") or "N/A"
        )
        display_title = (
            study_title if len(study_title) <= 60 else study_title[:57] + "..."
        )
        modality = record.get("data_modalities") or []
        year = record.get("study_year") or "N/A"

        modality_badges = []
        for mod in modality:
            badge_color = DATA_MODALITIES_COLOURS.get(mod, COLORS["primary"])
            modality_badges.append(
                dbc.Badge(
                    mod,
                    color="light",
                    className="me-1 mb-1",
                    style={
                        "fontSize": "0.75rem",
                        "border": f"1px solid {badge_color}",
                        "color": badge_color,
                        "backgroundColor": f"{badge_color}1A",
                    },
                )
            )
        if not modality_badges:
            modality_badges = [html.Span("N/A", className="text-muted")]

        provenance = reprocessed_badge(
            record.get("perturb_seq_reprocessed"), class_name=""
        )

        rows.append(
            html.Tr(
                [
                    html.Td(
                        dcc.Link(
                            dataset_id,
                            href=f"/perturbation-catalogue/dataset/{dataset_id}",
                            className="text-decoration-none fw-semibold",
                            style={"color": COLORS["primary"]},
                        )
                    ),
                    html.Td(
                        html.Span(display_title, title=study_title),
                        style={"maxWidth": "300px"},
                    ),
                    html.Td(str(year), className="text-center"),
                    html.Td(html.Div(modality_badges, className="d-flex flex-wrap")),
                    html.Td(provenance if provenance is not None else ""),
                ]
            )
        )

    table = dbc.Table(
        [header, html.Tbody(rows)],
        bordered=False,
        hover=True,
        responsive=True,
        striped=False,
        className="align-middle shadow-sm",
        style={"borderRadius": "12px", "overflow": "hidden"},
    )

    return table


def filter_placeholder(message="Search to enable filters."):
    """Placeholder message when filters are unavailable."""
    return dbc.Alert(
        [html.I(className="bi bi-funnel me-2"), message],
        color="light",
        className="shadow-sm",
        style={"borderRadius": "10px"},
    )


def build_filter_controls(
    facets,
    selected_filters=None,
    facet_fields=None,
    search_mode="targets",
    id_prefix="",
):
    """Build filter controls for facet fields.

    Args:
        facets: Facet data from the search API.
        selected_filters: Currently selected filter values.
        facet_fields: List of facet field names to display.
        search_mode: 'targets' or 'datasets'.
        id_prefix: Prefix for component IDs to avoid collisions across pages.
    """
    if not facets:
        return filter_placeholder()

    if selected_filters is None:
        selected_filters = {}

    if facet_fields is None:
        facet_fields = FACET_FIELDS

    filter_type = f"{id_prefix}-facet-filter" if id_prefix else "facet-filter"

    # Facet values that should be relabelled for display (value sent to the API is
    # unchanged; only the visible label differs).
    boolean_value_labels = {
        "perturb_seq_reprocessed": {"true": "Yes", "false": "No"},
    }

    field_icons = {
        "perturb_seq_reprocessed": "bi-arrow-repeat",
        "license": "bi-award-fill",
        "data_modalities": "bi-database",
        "tissues_tested": "bi-universal-access-circle",
        "cell_types_tested": "bi-puzzle",
        "cell_lines_tested": "bi-puzzle-fill",
        "diseases_tested": "bi-virus2",
        "sex_tested": "bi-gender-ambiguous",
        "developmental_stages_tested": "bi-graph-up-arrow",
        "license_labels": "bi-award-fill",
        "library_perturbation_type_labels": "bi-database",
        "tissue_labels": "bi-universal-access-circle",
        "cell_type_labels": "bi-puzzle",
        "cell_line_labels": "bi-puzzle-fill",
        "disease_labels": "bi-virus2",
        "sex_labels": "bi-gender-ambiguous",
        "developmental_stage_labels": "bi-graph-up-arrow",
    }

    field_explanations = {
        "license": {
            "header": "License Filter",
            "body": "Filter targets by the license types of their underlying datasets. Selecting a license shows all targets that have at least one dataset released under that license. Note: A single target may have datasets with different licenses.",
        },
        "data_modalities": {
            "header": "Data Modalities Filter",
            "body": "Filter targets by experimental approach (Perturb-seq, CRISPR screens, MAVE). Selecting a modality shows all targets that have data from that experimental type. A target may appear in multiple modalities.",
        },
        "tissues_tested": {
            "header": "Tissues Filter",
            "body": "Filter targets by the tissues used in experiments. Selecting a tissue shows all targets that have been studied in at least one dataset using that tissue. The aggregated results for a target may include data from multiple tissues.",
        },
        "cell_types_tested": {
            "header": "Cell Types Filter",
            "body": "Filter targets by cell types used in experiments. Selecting a cell type shows all targets studied in at least one dataset using that cell type. A single target's aggregated data may span multiple cell types.",
        },
        "cell_lines_tested": {
            "header": "Cell Lines Filter",
            "body": "Filter targets by cell lines used in experiments. Selecting a cell line shows all targets studied in at least one dataset using that cell line. A target's aggregated results may include data from multiple cell lines.",
        },
        "diseases_tested": {
            "header": "Diseases Filter",
            "body": "Filter targets by disease context of experiments. Selecting a disease shows all targets studied in at least one dataset related to that disease. A target may have been studied across multiple disease contexts.",
        },
        "sex_tested": {
            "header": "Sex Filter",
            "body": "Filter targets by the biological sex of samples used in experiments. Selecting a sex shows all targets studied in at least one dataset using samples of that sex. A target's aggregated data may include both sexes.",
        },
        "developmental_stages_tested": {
            "header": "Developmental Stages Filter",
            "body": "Filter targets by developmental stage of samples. Selecting a stage shows all targets studied in at least one dataset at that developmental stage. A target may have data across multiple developmental stages.",
        },
        "license_labels": {
            "header": "License Filter",
            "body": "Filter datasets by the license types under which they are released. Selecting a license shows all datasets released under that license.",
        },
        "library_perturbation_type_labels": {
            "header": "Data Modalities Filter",
            "body": "Filter datasets by experimental approach (Perturb-seq, CRISPR screens, MAVE). Selecting a modality shows all datasets that have data from that experimental type.",
        },
        "tissue_labels": {
            "header": "Tissues Filter",
            "body": "Filter datasets by the tissues used in experiments. Selecting a tissue shows all datasets that used that tissue.",
        },
        "cell_type_labels": {
            "header": "Cell Types Filter",
            "body": "Filter datasets by cell types used in experiments. Selecting a cell type shows all datasets that used that cell type.",
        },
        "cell_line_labels": {
            "header": "Cell Lines Filter",
            "body": "Filter datasets by cell lines used in experiments. Selecting a cell line shows all datasets that used that cell line.",
        },
        "disease_labels": {
            "header": "Diseases Filter",
            "body": "Filter datasets by disease context of experiments. Selecting a disease shows all datasets related to that disease.",
        },
        "sex_labels": {
            "header": "Sex Filter",
            "body": "Filter datasets by the biological sex of samples used in experiments. Selecting a sex shows all datasets using samples of that sex.",
        },
        "developmental_stage_labels": {
            "header": "Developmental Stages Filter",
            "body": "Filter datasets by developmental stage of samples. Selecting a stage shows all datasets at that developmental stage.",
        },
    }

    controls = []
    total_fields = len(facet_fields)
    for idx, field in enumerate(facet_fields):
        values = facets.get(field, [])
        if not values:
            continue

        display_name_map = {
            "perturb_seq_reprocessed": "Reprocessed",
            "license": "License",
            "data_modalities": "Data Modalities",
            "tissues_tested": "Tissues",
            "cell_types_tested": "Cell Types",
            "cell_lines_tested": "Cell Lines",
            "sex_tested": "Sex",
            "developmental_stages_tested": "Developmental Stages",
            "diseases_tested": "Diseases",
            "license_labels": "License",
            "library_perturbation_type_labels": "Data Modalities",
            "tissue_labels": "Tissues",
            "cell_type_labels": "Cell Types",
            "cell_line_labels": "Cell Lines",
            "disease_labels": "Diseases",
            "sex_labels": "Sex",
            "developmental_stage_labels": "Developmental Stages",
        }
        if field in display_name_map:
            display_name = display_name_map[field]
        else:
            display_name = field.replace("_", " ").title()
        options = []
        option_value_map = {}
        field_selected = [
            str(v).strip().lower()
            for v in selected_filters.get(field, [])
            if v is not None
        ]
        for item in values:
            raw_value = item.get("value")
            count = item.get("count", 0)
            if raw_value is None:
                continue
            value = str(raw_value).strip()
            if not value:
                continue
            if count <= 0 and value.lower() not in field_selected:
                continue
            value_label_map = boolean_value_labels.get(field)
            label_text = (
                value_label_map.get(value.lower(), value) if value_label_map else value
            )
            options.append({"label": f"{label_text} ({count})", "value": value})
            option_value_map[value.lower()] = value

        if not options:
            continue

        selected_values = []
        for item in selected_filters.get(field, []):
            if item is None:
                continue
            cleaned_value = str(item).strip()
            if not cleaned_value:
                continue
            mapped_value = option_value_map.get(cleaned_value.lower())
            if mapped_value and mapped_value not in selected_values:
                selected_values.append(mapped_value)

        if len(options) <= 10:
            control = dcc.Checklist(
                id={"type": filter_type, "field": field},
                options=options,
                value=selected_values,
                inputStyle={"marginRight": "0.5rem"},
                labelStyle={
                    "display": "block",
                    "marginBottom": "0.35rem",
                    "fontSize": "0.9rem",
                },
            )
        else:
            control = dcc.Dropdown(
                id={"type": filter_type, "field": field},
                options=options,
                value=selected_values,
                multi=True,
                placeholder=f"Filter by {display_name}",
                className="facet-dropdown",
                style={
                    "fontSize": "0.9rem",
                    "zIndex": 2000,
                    "position": "relative",
                },
            )

        icon_class = field_icons.get(field)
        explanation = field_explanations.get(field)
        info_icon_id = (
            f"{id_prefix}-facet-info-{field}" if id_prefix else f"facet-info-{field}"
        )

        header_children = []
        if icon_class:
            header_children.append(html.I(className=f"bi {icon_class} me-2"))

        header_children.append(f"{display_name} ")

        if explanation and search_mode != "datasets":
            header_children.append(
                html.Span(
                    html.I(className="bi bi-question-circle me-2"),
                    id=info_icon_id,
                    style={"cursor": "pointer", "marginLeft": "2px"},
                )
            )
            header_children.append(
                dbc.Popover(
                    [
                        dbc.PopoverHeader(explanation["header"]),
                        dbc.PopoverBody(explanation["body"]),
                    ],
                    target=info_icon_id,
                    trigger="click",
                    placement="right",
                )
            )

        controls.append(
            dbc.Card(
                [
                    dbc.CardHeader(
                        header_children,
                        className="fw-semibold",
                        style={"backgroundColor": "#f8f9fa"},
                    ),
                    dbc.CardBody(control, style={"padding": "0.75rem"}),
                ],
                className="mb-3 shadow-sm facet-filter-card",
                style={
                    "borderRadius": "12px",
                    "overflow": "visible",
                    "position": "relative",
                    "zIndex": total_fields - idx,
                },
            )
        )

    if not controls:
        return filter_placeholder("No filters available for these results.")

    return html.Div(controls, className="facet-controls")
