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

Below the subtitle and above the main table headers, add a `dbc.ButtonGroup` to control data grouping. This component should be centered.

The button group will contain two buttons: "Group by Perturbation" and "Group by Phenotype".

By default, the view should be grouped by "Perturbation" when the page first loads. The `group_by` parameter in the API call should be set to `perturbation_gene_name` by default.

### 4. Table header

The main concept is four columns, left to right: Dataset > Perturbation > Change > Phenotype.

The column headers need to be conspicuously labeled in a large font, with smaller explanation on top. So each header is two-line: first line in a smaller font provides a grammatical connectivity and context, while the second line is in a much larger font to serve as the column header.

The four titles are:
* (In / Dataset)
* (Introducing / Perturbation)
* (Leads to / Change)
* (In / Phenotype)

This is what it should look like on the page:
[in small font] In                   Introducing      Leads to     In
[in large font] Dataset              Perturbation     Change       Phenotype

The top element is an `h4` with classes `"fw-normal fst-italic mb-0"`.

The bottom element is an `h3` with default styling.

### 5. Search and filter fields

Immediately below the headers, there are corresponding search and filter fields.

- **Dataset:** `dbc.Input` with placeholder "Filter by dataset metadata".
- **Perturbation:** `dbc.Input` with placeholder "Filter by perturbed gene".
- **Change:** A `dbc.ButtonGroup` with three buttons: "Up", "Down", and "Both". "Both" should be the default selection. This controls the `change_direction` API parameter.
- **Phenotype:** `dbc.Input` with placeholder "Filter by phenotype gene".

When any of these fields are modified, the data table should update immediately.

### 6. Data grid

Below the search fields, data is displayed, as provided by the BE API according to the structure in @prompts/data-model.md. The layout is detailed in the "Layout and Alignment" section below.

Dataset ID, perturbed gene name and phenotype gene name should be most visible (for example in bold) because these are the main properties of each of the columns.

#### Dataset column representation
First, dataset ID is listed in prominent bold. Do not use a label. For example, instead of "Dataset: some_dataset_123", display "some_dataset_123" only.

Then, dataset metadata is listed, one on new line. Do not use colons, as they are visually distracting. This is bad: "Tissue type: Blood".

Instead, display labels (such as "Tissue type") in light, thin italics (`html.I`), and values (such as "Blood") in semibold text (`className="fw-semibold"`). Do not use a colon between them.

#### Perturbation, Change, and Phenotype column representation

Gene names in the "Perturbation" and "Phenotype" columns should be displayed in bold.

The "Change" column should display the direction of change with an arrow (▲ for "increased", ▼ for "decreased"), followed by the `log2fc` and `padj` values.

Similarly to dataset information, do not use colons for labels (e.g., "log2fc: 0.63"). Instead, display the label in thin, pale italics, and the value in a semibold font. Ensure consistent styling across all columns. Use good spacing between label-value pairs (e.g. `me-3`).

#### Value formatting

This applies to all columns.

For any values that are displayed, make sure their first letter only is capitalised, *except* for the dataset ID.

When displaying float values, make sure that - (hyphen) is replaced with a minus sign (`−`) for both negative values such as -0.13, and negative exponents of padj such as 1.03e-05.

### Layout, Alignment, and Grouping

To ensure precise and consistent column alignment that is robust to content and grouping changes, the layout for sections 4-6 (table header, search fields, data grid) must be implemented using a **CSS Grid**. This approach replaces the previous nested `dbc.Row`/`dbc.Col` structure with a single, flat grid container, providing direct control over column widths and cell placement.

#### Grid Definition
- The entire table-like structure (headers, filters, and data rows) should be contained within a single `html.Div` acting as the grid container.
- This container's style must include `display: 'grid'` and `grid-template-columns: 4fr 3fr 2fr 3fr`. This sets up the four main columns with the required `4:3:2:3` ratio, making them flexible and proportional.
- Use `gap` property (e.g. `row-gap: '0.5rem'`) for spacing between rows, but no column gap. Padding should be handled inside the cells. A `border-top` on the grid container can separate it from the filters.

#### Grid Content: Headers and Filters
- Headers and filter fields are straightforward: they will occupy the first four columns of the first two rows of the grid.

#### Grid Content: Data and Vertical Merging
The key to this layout is dynamically calculating the vertical span of cells that are "merged".

- **Dataset Cells:** A dataset cell will always be in the first column and will span multiple rows. You must calculate how many data rows belong to that dataset and set the `grid-row: 'span <number_of_rows>'` style property on the dataset cell's `html.Div`. The cell should also have a `border-top` to visually separate dataset groups.

- **Grouped Cells:** When grouping by "Perturbation" or "Phenotype", the grouped cell acts as a merged cell. It will span multiple rows corresponding to its children. This also requires calculating the number of child rows and applying the `grid-row: 'span <number_of_rows>'` style. The content of the grouped cell must be updated to show summary information:
    - **When grouping by perturbation:**
        * The cell is in the "Perturbation" column.
        * Display the `perturbation_gene_name` in a larger font.
        * Below, write in normal font: "Affects XXX phenotypes", where XXX is `n_total` (in bold).
        * Below: "▲ XXX up" (`n_up` in bold).
        * Below: "▼ XXX down" (`n_down` in bold).
    - **When grouping by phenotype:**
        * The cell is in the "Phenotype" column.
        * Display the `phenotype_gene_name` in a larger font.
        * Below, write in normal font: "Affected by XXX perturbations", where XXX is `n_total` (in bold).
        * Below: "▲ by XXX up" (`n_up` in bold).
        * Below: "▼ by XXX down" (`n_down` in bold).
        * Below: "Base mean expression is XX.XX" (`base_mean` in bold).

- **Standard Cells:** Regular, non-merged cells (like "Change" and the non-grouped "Phenotype"/"Perturbation") simply occupy a single cell in the grid. The remaining columns in the row will display the corresponding data for each entry in the group.

- **Alignment:** All cells should be top-aligned. This can be achieved by styling the grid items with `align-self: 'start'`. For content within cells, use flexbox (`display: 'flex', flexDirection: 'column'`) to control layout and ensure content is top-aligned and wraps correctly.

This CSS Grid-based structure ensures all columns across all rows are perfectly aligned to the master grid definition, solving the width consistency problem regardless of grouping.