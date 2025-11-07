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

The page is structured into a full-width banner at the top, followed by a two-column main content area.

- The page content should span the entire width of the container using `fluid=True`.
- The main content area below the banner should be a `dbc.Row` with horizontal padding (e.g., `px-5`) and a bottom margin (`mb-5`).
- The left column (`dbc.Col` with `width=3`) will contain the facet filters.
- The right column (`dbc.Col` with `width=9`) will contain all other controls and the main data grid (grouping/limit controls, results summary, headers, search fields, and the data itself).

Make sure to use Bootstrap components and styles as much as possible. Make sure that everything uses the same font (the default), is visually pleasing and coherent.

### Banner

A full-width banner should be at the top of the page.
- It should use the background image from `https://new-pc-fe-959149465821.europe-west2.run.app/assets/banner.png`.
- The background image should cover the banner area and be centered.
- The banner should contain the title and subtitle.
- The text inside the banner should be white with a text shadow to ensure readability against the background image.

### 1. Title

A very big title, centered, inside the banner: "Perturbation Catalogue". Use an `h1` tag with classes `"text-center display-4 mb-3 fw-bold"` and an inline style for `font-size: "60px"`.

### 2. Subtitle

A subtitle below the title, centered, inside the banner: "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.". Use a `p` tag with class `"text-center lead"` and inline styles for `font-size: "25px"`.

### 3. Facet Filters

A panel on the left side of the page for filtering data based on dataset metadata facets.

- The filter panel should be populated from the `facet_counts` object in the API response.
- For each key in `facet_counts` (e.g., "tissue", "cell_type"), a filter group should be created.
- Each group will have a header (e.g., `html.H5` with the capitalized and underscore-replaced facet name like "Cell Type").
- Below the header, a list of `dbc.Checkbox` components will be displayed for each available filter value. The label for each checkbox should show the value and the count (e.g., "blood (50)").
- When a checkbox is toggled, the page must trigger a new search with the corresponding filter applied.
- The list of available filter options in the facet panel must be updated after each query to reflect the new results.
- **Important:** The code for generating these filter groups must be procedural. Define a list or dictionary of the facet fields and their properties, and then iterate over it to create the components. Do not duplicate the layout code for each facet type.

### 4. Grouping and Limit Controls

Above the main data table (in the right-hand column), provide controls for grouping and limiting the results. This section should be a flex container to align items neatly.

- **Grouping:** This should consist of a `html.Span` with the text "Group by " followed by a `dbc.ButtonGroup`.
    - The button group will contain two buttons: "Perturbation" and "Phenotype".
    - The default selection is "Perturbation". The active button should have `color="primary"` and the inactive one should have `color="light"`.
    - The `dbc.ButtonGroup` should have `size="lg"`. The "Group by" text should have a font size of `1.25rem` to match the buttons.
- **Result Limits:** Near the grouping controls, add three `dbc.Input` fields for controlling the number of results returned. Each should have a label and be of `type="number"`.
    - **Max Datasets:** Label "Max datasets per modality", for the `max_datasets_per_modality` parameter. Default to 10.
    - **Max Groups:** Label "Max top-level groups", for the `max_top_level` parameter. Default to 10.
    - **Max Rows:** Label "Max rows per group", for the `max_rows` parameter. Default to 10.

This entire container should have a larger bottom margin (`mb-5`) to create space before the table headers.

### 5. Results Summary

Below the grouping/limit controls and before the table headers, display a summary of the total number of datasets found for each modality.

- This should be generated by iterating through the `modalities` list in the API response.
- For each modality, display a text element like: "Found **123** datasets for **perturb-seq**".
- The counts and modality names should be dynamically inserted from the `total_datasets_count` and `modality` fields.

### 6. Table header

The four column headers are:
* (According to / Dataset)
* (Introducing / Perturbation)
* (Leads to / Change)
* (Affecting / Phenotype)

This should be implemented as a two-line header where the top line is a `h4` with classes `"fw-normal fst-italic mb-0"` and the bottom line is an `h3`. The actual / must not be displayed.

### 7. Search and filter fields

Immediately below the headers, provide search and filter fields for each column:

- **Dataset:** `dbc.Input` with placeholder "Filter by dataset metadata".
- **Perturbation:** `dbc.Input` with placeholder "Filter by perturbed gene".
- **Change:** A `dbc.ButtonGroup` with three buttons: "△ Up", "▽ Down", and "△▽ Both". "Both" is the default. The active button has `color="secondary"`, inactive is `color="light"`.
- **Phenotype:** `dbc.Input` with placeholder "Filter by phenotype gene".

When any of these fields are modified, the data table should update immediately.

### 8. Data Grid Layout

To ensure precise and consistent column alignment, the layout for sections 6-8 (headers, filters, and the data itself) must be implemented as a single **CSS Grid** within the right-hand column.

- The container for the headers, filters, and data should be an `html.Div` with `display: 'grid'`.
- Define the column widths with `grid-template-columns: 4fr 3fr 3fr 3fr`.
- To make nested rows more compact, set `row-gap: '0.25rem'`.
- Set `column-gap: '1.5rem'` for more space between columns.
- To prevent rows from stretching vertically to fill space, set `grid-auto-rows: 'min-content'`.
- To create a clear visual break, a prominent, full-width separator must be placed after the filter controls and before the data grid begins. This should be a `div` with a thick top border (`2px`) and significant vertical margin.
- To handle aggregated data correctly, where one cell in an outer group (like Dataset) corresponds to multiple rows of child data, the implementation must use CSS Grid's `grid-row: span N` property.
- The callback that generates the data grid must first calculate the number of rows each aggregated cell needs to span, including rows for separators.
- The callback must then return a **flat list** of cell components. It must be constructed carefully so that for rows where a column is spanned by a cell from a previous row, that item is **omitted** from the list for the current row, maintaining the grid structure.
- Any message, such as "No results" or an error, should also be a child of the grid and span all columns (`style={'grid-column': '1 / -1'}`).

### 9. Data Column Rendering

This section specifies the exact rendering for the content within each data column.

#### Font Styles
- **Large & Bold:** A font size equivalent to `H3` with `fw-bold`. This is used for the `dataset_id` and the gene name of the primary, grouped-on item (e.g., the main Perturbation in a group).
- **Regular & Semi-bold:** The default browser font size with `fw-semibold`. This is used for sub-row gene names (the non-grouped items).
- **Regular:** The default browser font size. This is used for all other text, including dataset metadata and all informational text (e.g., "Affects phenotypes", "log2fc").

#### Universal Styling for Labels and Values
A consistent styling pattern must be used for all informational text.
- **Labels** (e.g., "Tissue", "Affects", "padj"): Rendered using a light font weight (`fw-light`).
- **Values** (e.g., "Blood", "502", "1.23e-5"): Rendered in a `html.Span` with `className="fw-semibold"`.
- For negative numeric values, the standard hyphen-minus (-) must be replaced with the Unicode minus sign (`−`).

#### Vertical Spacing
A 3-tier vertical spacing model must be implemented:
1.  **Sub-item Spacing:** The smallest space, between individual rows of data, is handled by the grid's `row-gap` (`0.25rem`).
2.  **Group Spacing:** A continuous horizontal line (`1px`) must separate aggregated groups within the same dataset. This should be implemented as a dedicated `html.Div` spanning the three rightmost columns (`grid-column: 2 / 5`), with vertical margins (e.g., `margin-top: 1rem`, `margin-bottom: 1rem`) to create even spacing.
3.  **Dataset Spacing:** The largest separation is between datasets.. This is achieved with a thicker, full-width horizontal line (`2px`) with significant, even vertical spacing above and below it (e.g., `margin-top: 2rem`, `padding-top: 2rem`).

#### Dataset Column
- The `dataset_id` is displayed in a **Large & Bold** font. It must be formatted to replace underscores with spaces and capitalize the first letter (e.g., `adamson_2016_pilot` becomes `Adamson 2016 pilot`).
- Below the ID, list all other metadata fields on new lines, using the **Regular** font size. Apply the universal label/value styling.
- The **first letter** of each metadata *value* must be **capitalized**.

#### Perturbation Column
This cell has two rendering modes:
- **When Grouped On:** Displays the `perturbation_gene_name` in a **Large & Bold** font with reduced bottom margin (`mb-1`). Below this, on a new line, it shows a conditional label: "Multiple gene knockout" if the gene name contains a `|`, otherwise "Single gene knockout". This is followed by the informational text in the **Regular** font size:
        1. "Affects **N** phenotypes"
        2. "△ **X** ▽ **Y**" (e.g., "△ 225 ▽ 277"). The up-regulated count and icon must be green (`#2acc06`), and the down-regulated count and icon must be red (`#ff4824`).- **When in a Sub-Row:** Displays only the `perturbation_gene_name` in **Regular & Semi-bold**.

#### Change Column
- To ensure vertical alignment of values across rows, the content must be structured using a flexbox container (`display: flex`).
- Each component (`padj` info, `log2fc` info, direction arrow) should be in its own `html.Div` with a `min-width` to create stable columns.
- Displays `padj` and its value, then `log2fc` and its value.
- Finally, display a Unicode triangle symbol (`△` for up, `▽` for down).
- The triangles must be colored with bright, specific shades: green (`#2acc06`) for up and red (`#ff4824`) for down.

#### Phenotype Column
This cell has two rendering modes:
- **When Grouped On:** Displays the `phenotype_gene_name` in a **Large & Bold** font with reduced bottom margin (`mb-1`). Below this, on new lines, it shows informational text in the **Regular** font size:
    1. "Base expression **Z**"
    2. "Affected by **N** perturbations"
    3. "△ **X** ▽ **Y**". The up-regulated count and icon must be green (`#2acc06`), and the down-regulated count and icon must be red (`#ff4824`).
- **When in a Sub-Row:** Displays only the `phenotype_gene_name` in **Regular & Semi-bold**.

### 10. Truncation Notice

When the number of returned elements for any level (datasets per modality, top-level groups, or rows per group) is less than the total available elements (as indicated by the API response, e.g., `total_datasets_count` for modalities, or the `n_total` fields within perturbation/phenotype objects), a truncation notice must be displayed.

- The notice should have a light grey background (`#f8f9fa`), padding, and rounded corners.
- The text should be center-aligned, italicized, and read "Displaying top X elements", where X is the number of elements actually displayed.
- This notice should be placed appropriately within the data grid, spanning the relevant columns. For example, if datasets are truncated, the notice should appear at the end of the dataset list. If top-level groups are truncated within a dataset, it should appear at the end of that group's list.