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
- Define the column widths with `grid-template-columns: 4fr 3fr 3fr 3fr`.
- Set `row-gap: '0.5rem'` and `column-gap: '1rem'`.
- To create visual separation, a spacer row with a height of `1.5rem` must be added between the filter controls and the data grid itself.
- To handle aggregated data correctly, where one cell in an outer group (like Dataset) corresponds to multiple rows of child data, the implementation must use CSS Grid's `grid-row: span N` property.
- The callback that generates the data grid must first calculate the number of rows each aggregated cell needs to span.
- The callback must then return a **flat list** of cell components. It must be constructed carefully so that for rows where a column is spanned by a cell from a previous row, that item is **omitted** from the list for the current row, maintaining the grid structure.
- Any message, such as "No results" or an error, should also be a child of the grid and span all columns (`style={'grid-column': '1 / -1'}`).

### 7. Data Column Rendering

This section specifies the exact rendering for the content within each data column.

#### Font Sizes
Only two font sizes should be used within the data grid cells:
- **Large:** A font size equivalent to `H3`. This is used for the `dataset_id` and the gene name of the primary, grouped-on item (e.g., the main Perturbation in a group).
- **Regular:** The default browser font size. This is used for all other text, including dataset metadata, all informational text (e.g., "Affects phenotypes", "up, down", "log2fc"), and direction arrows.

#### Universal Styling for Labels and Values
A consistent styling pattern must be used for all informational text.
- **Labels** (e.g., "Tissue", "Affects", "log2fc"): Rendered using a light font weight (`fw-light`).
- **Values** (e.g., "Blood", "502", "1.23"): Rendered in a `html.Span` with `className="fw-semibold"`.
- For negative numeric values, the standard hyphen-minus (-) must be replaced with the Unicode minus sign (`\u2212`).

#### Vertical Spacing
A 3-tier vertical spacing model must be implemented:
1.  **Sub-item Spacing:** The smallest space, between individual rows of data, is handled by the grid's `row-gap` (`0.5rem`).
2.  **Group Spacing:** A larger vertical space (`1.5rem`) with a horizontal line must be created between aggregated groups within the same dataset. This is achieved by applying a `border-top` and `padding-top` to the first row of cells belonging to a new group (spanning columns 2, 3, and 4).
3.  **Dataset Spacing:** The largest separation is between datasets. This is achieved with a full-width horizontal line that has `2rem` of `padding-top`.

#### Dataset Column
- The `dataset_id` is displayed in a **Large** font size (`H3`) and bold weight.
- Below the ID, list all other metadata fields on new lines, using the **Regular** font size. Apply the universal label/value styling.
- The **first letter** of each metadata *value* must be **capitalized**.

#### Perturbation Column
This cell has two rendering modes:
- **When Grouped On:** Displays the `perturbation_gene_name` in a **Large**, semi-bold font (`H3` with `fw-semibold`). Below this, on new lines, it shows informational text in the **Regular** font size:
    1. "Affects **N** phenotypes"
    2. "**X** up, **Y** down"
- **When in a Sub-Row:** Displays only the `perturbation_gene_name` in the **Regular** font size.

#### Change Column
- All content must be on a **single line** and in the **Regular** font size.
- Displays `log2fc` and its value, then `padj` and its value, with extra margin (`me-3`) between them.
- The `logfc` label must be corrected to `log2fc`.
- Finally, display a Unicode arrow symbol (`🡱` for up, `🡳` for down).
- The arrows must be colored with bright, specific shades: green (`#2acc06`) for up and red (`#ff4824`) for down.

#### Phenotype Column
This cell has two rendering modes:
- **When Grouped On:** Displays the `phenotype_gene_name` in a **Large**, semi-bold font (`H3` with `fw-semibold`). Below this, on new lines, it shows informational text in the **Regular** font size:
    1. "Base expression **Z**"
    2. "Affected by **N** perturbations"
    3. "**X** up, **Y** down"
- **When in a Sub-Row:** Displays only the `phenotype_gene_name` in the **Regular** font size.