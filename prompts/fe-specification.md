# Front-end specification

## Introduction

Study general instructions for LLM agents in @prompts/README.md.

Study the data model and API contract in @prompts/data-model.md.

Then, read the BE implementation which is based on the above data model in @be/main.py.

Be aware that FE is configured with the $PERTURBATION_CATALOGUE_BE environment variable, which has a BE which was implemented according to this data model. This is what you should use as your endpoint.

Study how relevant FE parts are organised in in @fe/app.py, @fe/pages/home.py, @fe/pages/_order.py. Do not peek into other FE components, as they are irrelevant and can be quite large when loaded into context.

## Pages update

Add a new page named "Perturbations", flag it with a conspicous yellow "New!" label on the front page, and make its URL /perturbations.

For the new page icon in the navbar, use the "tropical-storm" Bootstrap icon.

## New page specification

### Overall layout

The page layout, top to bottom, is: title, subtitle, grouping controls, table header, search fields, and a data grid obtained from BE according to filters.

The page content should span the entire width of the container using `fluid=True`. The main content should have horizontal padding (e.g., `px-5`).

Make sure to use Bootstrap components and styles as much as possible. Make sure that everything uses the same font (the default), is visually pleasing and coherent.

### 1. Title

A very big title, centered: "Perturbation Catalogue". Use an `h1` tag with classes `"text-center display-4 mt-5 mb-3"` and an inline style for `font-size: "60px"`.

### 2. Subtitle

A subtitle below the title, centered: "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.". Use a `p` tag with classes `"text-center lead mb-5"` and an inline style for `font-size: "25px"`.

### 3. Grouping controls

Below the subtitle, add controls for data grouping. This should consist of a `html.Span` with the text "Group by: " followed by a `dbc.ButtonGroup`.

The button group will contain two buttons: "Perturbation" and "Phenotype". Their width should be determined by the text, not stretched. The default selection is "Perturbation". The active button should have `color="primary"` and the inactive one should have `color="light"`.

### 4. Table header

The four column headers are:
* (According to / Dataset)
* (Introducing / Perturbation)
* (Causes / Change)
* (In / Phenotype)

This should be implemented as a two-line header where the top line is a `h4` with classes `"fw-normal fst-italic mb-0"` and the bottom line is an `h3`.

### 5. Search and filter fields

Immediately below the headers, provide search and filter fields for each column:

- **Dataset:** `dbc.Input` with placeholder "Filter by dataset metadata".
- **Perturbation:** `dbc.Input` with placeholder "Filter by perturbed gene".
- **Change:** A `dbc.ButtonGroup` with three buttons: "Up", "Down", and "Both". "Both" is the default. The active button has `color="primary"`, inactive is `color="light"`.
- **Phenotype:** `dbc.Input` with placeholder "Filter by phenotype gene".

When any of these fields are modified, the data table should update immediately.

### 6. Data Grid Layout

To ensure precise and consistent column alignment, the layout for sections 4-6 (headers, filters, data grid) must use a **CSS Grid**. This replaces nested `dbc.Row`/`dbc.Col` structures.

- The entire table structure is a single `html.Div` with `display: 'grid'` and `grid-template-columns: 4fr 3fr 2fr 3fr`.
- Use `row-gap: '0.5rem'` for spacing between grid rows and `column-gap: '1rem'` for spacing between columns.
- There should be **no** border between the grouping controls and the headers.
- A `border-top` must be applied **only** to the dataset cells to create a visual separation between different dataset groups. No other borders should be present in the data grid.
- All cells must be top-aligned using `align-self: 'start'`.

### 7. Data Column Rendering

This section specifies the exact rendering for the content within each data column. This rendering must be applied consistently, regardless of whether the item is being grouped on or is part of a nested data row.

#### Dataset Column
- The `dataset_id` is displayed first in a bold `h5` tag.
- Below the ID, all other metadata fields (Tissue, Cell type, etc.) are listed. Each field must be on a **new line** and formatted as "Label: Value" (e.g., "Tissue: Blood").

#### Perturbation Column
- The `perturbation_gene_name` is displayed in a bold `h5` tag.
- On a new line: "Affects **N** phenotypes", where N is the `n_total` value in bold.
- On a third line: "▲ **N** up | ▼ **N** down", where N is `n_up` and `n_down` respectively, in bold.

#### Change Column
- All content must be on a **single line**.
- Display `log2fc` first, with the label "log2fc:" and the value formatted to **two decimal places**.
- Next, display `padj`, with the label "padj:" and the value formatted in **scientific notation** (e.g., 1.23e-05).
- Finally, display an arrow: ▲ for "increased" or ▼ for "decreased". The arrow must be colored **green** for increased and **red** for decreased.

#### Phenotype Column
- The `phenotype_gene_name` is displayed in a bold `h5` tag.
- On a new line: "Base mean: **N**", where N is the `base_mean` value, formatted as an **integer** and in bold.
- On a new line: "Affected by **N** perturbations", where N is `n_total` in bold.
- On a final line: "▲ by **N** up | ▼ by **N** down", where N is `n_up` and `n_down` respectively, in bold.