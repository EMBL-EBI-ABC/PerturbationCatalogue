# Front-end specification

## Introduction

Study general instructions for LLM agents in @prompts/README.md.

Study the data model and API contract in @prompts/data-model.md.

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

### Layout and Alignment

The layout for sections 4-6 (table header, search fields, data grid) must be strictly aligned using a 12-column Bootstrap grid system (`dash-bootstrap-components`).

#### Column Ratios
The four main columns (Dataset, Perturbation, Change, Phenotype) should consistently use a `4:3:2:3` width ratio. This applies to the headers, search fields, and the data grid.

#### Data Grid Structure
The data grid requires a specific nested structure to achieve the desired layout where the dataset information acts as a vertically merged cell.

- Each dataset entry should be a single `dbc.Row` with `align="start"` to ensure its child columns are top-aligned and have independent heights. A `border-top` should separate these main dataset rows.
- This main row is split into two columns:
    1. A `dbc.Col` with `width=4` for the "Dataset" information.
    2. A `dbc.Col` with `width=8` that will contain all the data rows for that dataset.
- Inside the `width=8` column, each data entry is rendered as its own nested `dbc.Row`.
- These nested rows should have `align="start"`. No borders should be used between them.
- Each nested row is split into three columns that maintain the `3:2:3` ratio relative to their parent. Since the parent column has a width of 8, the nested columns for "Perturbation", "Change", and "Phenotype" should have widths of `3`, `2`, and `3` respectively.
- All data content within columns should be top-aligned. For columns containing multiple horizontal elements, use `d-flex align-items-start flex-wrap` to ensure content is top-aligned and wraps correctly.

# Note on grouping

Grouping is controlled by the `dbc.ButtonGroup` as described in section 3. The BE provides aggregated data based on the `group_by` parameter.

When data is grouped, the grouped column ("Perturbation" or "Phenotype") becomes a "vertically-merged" cell, similar to the "Dataset" column, and displays additional summary statistics.

## UI specifics for grouped views
When grouping by perturbation:
* Display the `perturbation_gene_name` in a larger font.
* Below, write in normal font: "Affects XXX phenotypes", where XXX is `n_total` (in bold).
* Below: "▲ XXX up" (`n_up` in bold).
* Below: "▼ XXX down" (`n_down` in bold).

When grouping by phenotype:
* Display the `phenotype_gene_name` in a larger font.
* Below, write in normal font: "Affected by XXX perturbations", where XXX is `n_total` (in bold).
* Below: "▲ by XXX up" (`n_up` in bold).
* Below: "▼ by XXX down" (`n_down` in bold).
* Below: "Base mean expression is XX.XX" (`base_mean` in bold).

The remaining columns in the row will display the corresponding data for each entry in the group.