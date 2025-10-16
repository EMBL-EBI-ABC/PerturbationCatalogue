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

Below the subtitle, add controls for data grouping. This should consist of a `html.Span` with the text "Group by " followed by a `dbc.ButtonGroup`.

The button group will contain two buttons: "Perturbation" and "Phenotype". Their width should be determined by the text, not stretched. The default selection is "Perturbation". The active button should have `color="primary"` and the inactive one should have `color="light"`.

### 4. Table header

The four column headers are:
* (According to / Dataset)
* (Introducing / Perturbation)
* (Leads to / Change)
* (Affecting / Phenotype)

This should be implemented as a two-line header where the top line is a `h4` with classes `"fw-normal fst-italic mb-0"` and the bottom line is an `h3`. The actual / must not be displayed.

### 5. Search and filter fields

Immediately below the headers, provide search and filter fields for each column:

- **Dataset:** `dbc.Input` with placeholder "Filter by dataset metadata".
- **Perturbation:** `dbc.Input` with placeholder "Filter by perturbed gene".
- **Change:** A `dbc.ButtonGroup` with three buttons: "Up", "Down", and "Both". "Both" is the default. The active button has `color="primary"`, inactive is `color="light"`.
- **Phenotype:** `dbc.Input` with placeholder "Filter by phenotype gene".

When any of these fields are modified, the data table should update immediately.

### 6. Data Grid Layout

To ensure precise and consistent column alignment, the layout for sections 4-6 (headers, filters, and the data itself) must be implemented as a single **CSS Grid**.

- The container for the headers, filters, and data should be an `html.Div` with `display: 'grid'`.
- Define the column widths with `grid-template-columns: 4fr 3fr 2fr 3fr`.
- Set `row-gap: '0.5rem'` and `column-gap: '1rem'`.
- To create visual separation, a spacer row with a height of `1.5rem` must be added between the filter controls and the data grid itself.
- To handle aggregated data correctly, where one cell in an outer group (like Dataset or a grouped Perturbation/Phenotype) corresponds to multiple rows of child data, the implementation must use CSS Grid's `grid-row: span N` property.
- The callback that generates the data grid must first calculate the number of rows each aggregated cell needs to span.
- The callback must then return a **flat list** of cell components. It must be constructed carefully so that for rows where a column is spanned by a cell from a previous row, that item is **omitted** from the list for the current row, maintaining the grid structure.
- Any message, such as "No results" or an error, should also be a child of the grid and span all columns (`style={'grid-column': '1 / -1'}`).
- To create a clear visual separation between datasets, a **full-width horizontal divider** must be placed *between* each top-level dataset group. This should be an `html.Div` that spans all four columns (`grid-column: '1 / -1'`) and has a `border-top`. A small amount of `padding` should be added to the divider to create vertical space.
- There should be **no other borders or dividers** within the data grid (e.g., no borders on individual cells).
- All cells must be strictly top-aligned using `align-self: 'start'`. Do not add any `padding-top` or `margin-top` to individual cell rendering functions.

### 7. Data Column Rendering

This section specifies the exact rendering for the content within each data column.

#### Universal Styling for Labels and Values
A consistent styling pattern must be used for all informational text across all four columns.
- **Labels** (e.g., "Tissue", "Affects", "logfc"): Rendered using a light font weight (e.g., `fw-light`).
- **Values** (e.g., "Blood", "502", "1.23"): Rendered in a `html.Span` with `className="fw-semibold"`.
- There should be **no colon** between a label and its value.
- For negative numeric values (`log2fc`, `padj` exponents), the standard hyphen-minus (-) must be replaced with the Unicode minus sign (`\u2212`).

#### Layout Principles
- For the Perturbation, Change, and Phenotype columns, a two-component layout should be used. This should be implemented using a flexbox (`display: 'flex', align-items: 'start'`) to ensure the components are top-aligned.
- The secondary informational text (e.g., "Affects X phenotypes", "logfc 1.23") should be made slightly smaller than the default font size (e.g., by wrapping it in a `html.Small` tag).

#### Vertical Spacing
A 3-tier vertical spacing model must be implemented:
1.  **Sub-item Spacing:** The smallest space, between individual rows of data (e.g., between multiple Change/Phenotype rows for a single Perturbation), is handled by the grid's `row-gap` (`0.5rem`).
2.  **Group Spacing:** A larger vertical space (`1.5rem`) must be created between aggregated groups within the same dataset. This is achieved by adding a `margin-top` to the primary cell of each group (excluding the first one).
3.  **Dataset Spacing:** The largest separation is between datasets. This is achieved with a full-width horizontal line that has `2rem` of `padding-top`.

#### Dataset Column
- The `dataset_id` is displayed first in a bold `h5` tag.
- Below the ID, list all other metadata fields. Each field must be on a **new line**. Apply the universal label/value styling.
- The **first letter** of each metadata *value* must be **capitalized**.

#### Perturbation Column
- A two-component, top-aligned layout.
- **Left Component:** The `perturbation_gene_name` displayed in a large, semi-bold font (e.g., `html.H3` with `fw-semibold`).
- **Right Component:** A block of smaller text containing two lines:
    1. "Affects **N** phenotypes" (where N is `n_total`).
    2. "**X** up, **Y** down" (where X/Y are `n_up`/`n_down`).
- The universal label/value styling must be applied. Do not use arrow symbols here.

#### Change Column
- A two-component, top-aligned layout, with the order reversed.
- **Left Component:** A block of smaller text containing two lines:
    1. "logfc **X**" (formatted to two decimal places).
    2. "padj **Y**" (formatted in scientific notation).
- **Right Component:** A large, equilateral triangle icon (e.g., `bi-triangle-fill` from Bootstrap Icons) with a `font-size` of `2.5rem`. The icon for "decreased" must be rotated 180 degrees. The icon must be colored with distinct, soft shades: a pleasant green (`#34d399`) for increased and a soft red (`#f87171`) for decreased.
- The universal label/value styling must be applied.

#### Phenotype Column
- A two-component, top-aligned layout, identical in style to the Perturbation column.
- **Left Component:** The `phenotype_gene_name` displayed in a large, semi-bold font (e.g., `html.H3` with `fw-semibold`).
- **Right Component:** A block of smaller text containing three lines:
    1. "Affected by **N** perturbations" (where N is `n_total`).
    2. "**X** up, **Y** down" (where X/Y are `n_up`/`n_down`).
    3. "Base expression **Z**" (where Z is `base_mean`, formatted as an integer with thousands separators).
- The universal label/value styling must be applied. Do not use arrow symbols here.