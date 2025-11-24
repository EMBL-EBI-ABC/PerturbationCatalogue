# Front-end specification

## Introduction

Study general instructions for LLM agents in @prompts/README.md.

Study the data model and API contract in @prompts/data-model.md.

Then, read the BE implementation which is based on the above data model in @be/data_query.py.

Be aware that FE is configured with the $PERTURBATION_CATALOGUE_BE environment variable, which has a BE which was implemented according to this data model. This is what you should use as your endpoint.

This document specifies a "Target Details" page and a reusable data table component it uses.

## Target Details Page

### URL

The page will be accessible at the URL `/target/<TARGET_NAME>`, where `<TARGET_NAME>` is the name of the gene to be displayed.

### Overall Layout

The page will display information for the target gene specified in the URL. It will not have a banner, title, or subtitle. The main content of the page will consist of four sections, each displaying a table of perturbation data relevant to the target gene.

The page should have a main title, e.g. `html.H1(f"Target: {target_name}")`.

Each of the four sections should have a subtitle, e.g. `html.H3(...)`.

The four sections are:
1.  **CRISPR screen data**: A table showing CRISPR screen data where the target gene was perturbed.
2.  **MAVE data**: A table showing MAVE data where the target gene was perturbed.
3.  **Perturb-Seq (Perturbed)**: A table showing Perturb-Seq data where the target gene was the one being perturbed (`perturbation_gene_name`).
4.  **Perturb-Seq (Affected)**: A table showing Perturb-Seq data where the target gene was one of the affected genes (`effect_gene_name`).

Each section will consist of a title and an instance of the reusable `DataTable` component.

## Reusable Data Table Component

To avoid code duplication, the data table must be defined as a separate, customisable component in its own Python file. The target details page will then import this component and create four instances of it, each with different data and configuration.

### Component Inputs (Properties)

The component should accept the following properties:
- `data`: The list of `DatasetResult` objects from the API.
- `modality`: The modality of the data being displayed (`crispr-screen`, `mave`, or `perturb-seq`).

### Data Fetching

- **Initial Population**: The table should initially be populated using the `/v1/{MODALITY}/search` API endpoint, filtered by the target gene from the URL.
- **Fetching More Rows**: To fetch additional rows for a specific dataset within the table (e.g., for pagination or infinite scrolling), the `/v1/{MODALITY}/{DATASET_ID}/search` API endpoint should be used.

### Layout

The component will render the data in a grid layout.

#### 1. Table Header

The table will have a three-column header:
*   (According to / **Dataset**)
*   (Introducing / **Perturbation**)
*   (Resulting in / **Effect**)

This should be implemented as a two-line header where the top line is a `h4` with classes `"fw-normal fst-italic mb-0"` and the bottom line is an `h3`. The slashes (`/`) must not be displayed.

#### 2. Data Grid Layout

- The layout for the header and the data must be implemented as a single **CSS Grid**.
- The container should be an `html.Div` with `display: 'grid'`.
- Define column widths with `grid-template-columns: 4fr 4fr 4fr`.
- Use a `row-gap` of `0.25rem` and a `column-gap` of `1.5rem`.
- Set `grid-auto-rows: 'min-content'` to prevent rows from stretching.
- To handle aggregated data where one dataset cell corresponds to multiple result rows, the implementation must use CSS Grid's `grid-row: span N` property for the Dataset column.
- The callback that generates the data grid must first calculate the number of rows each dataset cell needs to span.
- The callback must return a **flat list** of cell components. For rows where a column is spanned by a cell from a previous row, that item must be **omitted** from the list for the current row.
- Any message, such as "No results" or an error, should also be a child of the grid and span all columns (`style={'grid-column': '1 / -1'}`).

#### 3. Data Column Rendering

This section specifies the rendering for the content within each data column.

##### Font Styles
- **Large & Bold:** A font size equivalent to `H3` with `fw-bold`. Used for `dataset_id` and primary gene names.
- **Regular & Semi-bold:** The default browser font size with `fw-semibold`. Used for secondary gene names.
- **Regular:** The default browser font size. Used for all other informational text.

##### Universal Styling for Labels and Values
- **Labels** (e.g., "Tissue", "padj"): Rendered using a light font weight (`fw-light`).
- **Values** (e.g., "Blood", "1.23e-5"): Rendered in a `html.Span` with `className="fw-semibold"`.
- For negative numeric values, use the Unicode minus sign (`−`) instead of a hyphen-minus.

##### Vertical Spacing
- **Sub-item Spacing:** Handled by the grid's `row-gap` (`0.25rem`).
- **Dataset Spacing:** A thick, full-width horizontal line (`2px`) with significant vertical margin should separate each dataset group.

##### Dataset Column
- The `dataset_id` is displayed in a **Large & Bold** font. It should be formatted to replace underscores with spaces and capitalize the first letter (e.g., `adamson_2016_pilot` becomes `Adamson 2016 pilot`).
- Below the ID, list other metadata fields (`tissue`, `cell_type`, etc.) on new lines using the **Regular** font size and the label/value styling.
- The first letter of each metadata *value* must be capitalized.

##### Perturbation Column
- Displays the `perturbation_gene_name` in **Large & Bold** font.
- For `perturb-seq`, display the summary stats (`n_total`, `n_up`, `n_down`) below the gene name, using the label/value styling. The up/down counts and icons (`△`/`▽`) should be colored green (`#2acc06`) and red (`#ff4824`).
- For `mave`, if `perturbation_name` exists, display it below the gene name.

##### Effect Column
This column's content depends on the modality.

- **For `perturb-seq`:**
    - Display the `effect_gene_name` in **Regular & Semi-bold** font.
    - Display `log2fc`, `padj`, and `base_mean` using the label/value style.
    - Display the `effect_direction` with a colored triangle (`△` for "increased", `▽` for "decreased").
    - Display the summary stats (`n_total`, `n_up`, `n_down`) for the effect gene.

- **For `crispr-screen` and `mave`:**
    - Display the `effect_score_name` and `effect_score_value` using the label/value style.

#### 4. Truncation Notice

If the number of results within a dataset is limited by the API (`rows_per_dataset_limit`), a truncation notice should be displayed at the end of that dataset's rows.
- The notice should have a light grey background (`#f8f9fa`), padding, and rounded corners.
- The text should be center-aligned, italicized, and read "Displaying top X results".
- This notice should span the Perturbation and Effect columns.