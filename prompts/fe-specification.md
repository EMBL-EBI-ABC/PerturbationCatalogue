# Front-end specification

## Introduction

Study general instructions for LLM agents in @prompts/README.md.

Study how relevant BE parts are organised in @be/main.py. Be aware that FE is configured with the $PERTURBATION_CATALOGUE_BE env variable - this is what you should use as your endpoint. The endpoint to be used is /v1/search.

Study how relevant FE parts are organised in in @fe/app.py, @fe/pages/home.py, @fe/pages/_order.py. Do not peek into other FE components, as they are irrelevant and can be quite large when loaded into context.

## Pages update

Add a new page named "Perturbations", flag it with a conspicous yellow "New!" label on the front page, and make its URL /perturbations.

For the new page icon in the navbar, use the "tropical-storm" Bootstrap icon.

## New page specification

### Overall layout

The page layout, top to bottom, is: title, subtitle, table header, search fields, and a data grid obtained from BE according to filters.

The page content should span the entire width of the container using `fluid=True`. The main content should have horizontal padding (e.g., `px-5`).

Make sure to use Bootstrap components and styles as much as possible. Make sure that everything uses the same font (the default), is visually pleasing and coherent.

### 1. Title

A very big title, centered: "Perturbation Catalogue". Use an `h1` tag with classes `"text-center display-4 mt-5 mb-3"` and an inline style for `font-size: "60px"`.

### 2. Subtitle

A subtitle below the title, centered: "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.". Use a `p` tag with classes `"text-center lead mb-5"` and an inline style for `font-size: "25px"`.

### 3. Table header

The main concept is three columns, left to right: dataset > perturbation > effect.

The column headers need to be conspicuously labeled in a large font, with smaller explanation on top. So each header is two-line: first line in a smaller font provides a grammatical connectivity and context, while the second line is in a much larger font to serve as the column header.

The three titles are:
* (According to this / Dataset)
* (Introducing this / Perturbation)
* (Leads to the following) / (Effect)

This is what it should look like on the page:
[in small font] According to this    Introducing this    Leads to the following
[in large font] Dataset              Perturbation        Effect

The top element is an `h4` with classes `"fw-normal fst-italic mb-0"`.

The bottom element is an `h3` with default styling.

### 4. Search fields

Immediately below the headers, there are three corresponding search fields. Make sure these are nice Bootstrap search fields (`dbc.Input`) and not the bare, default ones.

Their descriptions (placeholders which are inside the field, in a pale colour) shoud read: "Filter by dataset metadata", "Filter by perturbed gene", "Filter by affected gene".

When either of the search fields is affected, the data table is immediately updated. Hence, no need for a separate Search button.

### 5. Data grid

Below the search fields, data is displayed, as provided by the BE API. The layout is detailed in the "Layout and Alignment" section below.

Note that in the BE response, the effect + perturbation values are currently all in a single object. Example of such an object:

```
{"perturbation_gene_name":"MARS1","effect_gene_name":"BRCA1","log2foldchange":0.6373134697742381,"padj":8.753715655002655e-15,"base_mean":69.076813168962}
```

Here, `perturbation_gene_name` should be in the "Perturbation" section, and the rest of the fields should be in the "Effect" section.

Dataset ID, perturbed gene name and affected gene name should be most visible (for example in bold) because these are the main properties of each of the columns.

#### Dataset column representation
First, dataset ID is listed in prominent bold. Do not use a label. For example, instead of "Dataset: some_dataset_123", display "some_dataset_123" only.

Then, dataset metadata is listed, one on new line. Do not use colons, as they are visually distracting. This is bad: "Tissue type: Blood".

Instead, display labels (such as "Tissue type") in light, thin italics (`html.I`), and values (such as "Blood") in semibold text (`className="fw-semibold"`). Do not use a colon between them.

#### Perturbation and effect column represetntation

Label should not be used for the gene name in both columns, it should simply be displayed as the very leftmost value in bold.

Importantly, all contents for the "Effect" section should be in one row to make the table visually compact. No newline breaks unless the content wraps. In comparison, the "Datasource" information mentioned above can and should be in multiple lines.

Depending on whether log2fc is positive (increased) or negative (decreased), display an arrow up or arrow down next to the effect gene name. Use not regular arrows, but ▲▼ arrows.

Similarly to dataset information, do not use colons such as "base mean: 0.288", instead again display label (base mean) in thin, pale italics, and value in semibold font. Make sure the labels and values styling is consistent between all three main columns. Good spacing should be between label-value pairs for the Effect column (e.g. `me-3`).

#### Value formatting

This applies to all columns.

For any values that are displayed, make sure their first letter only is capitalised, *except* for datasource ID.

When displaying float values, make sure that - (hyphen) is replaced with a minus sign (`−`) for both negative values such as -0.13, and negative exponents of padj such as 1.03e-05.

### Layout and Alignment

The layout for sections 3-5 (table header, search fields, data grid) must be strictly aligned using a 12-column Bootstrap grid system (`dash-bootstrap-components`).

#### Column Ratios
The three main columns (Dataset, Perturbation, Effect) should consistently use a `4:2:6` width ratio. This applies to the headers, search fields, and the data grid.

#### Data Grid Structure
The data grid requires a specific nested structure to achieve the desired layout where the dataset information acts as a vertically merged cell.

- Each dataset entry should be a single `dbc.Row` with `align="start"` to ensure its child columns are top-aligned and have independent heights. A `border-top` should separate these main dataset rows.
- This main row is split into two columns:
    1. A `dbc.Col` with `width=4` for the "Dataset" information.
    2. A `dbc.Col` with `width=8` that will contain all the perturbation rows for that dataset.
- Inside the `width=8` column, each perturbation is rendered as its own nested `dbc.Row`.
- These nested rows should have `align="start"`. No borders should be used between these nested rows; the only visual separator should be the main `border-top` between datasets. This applies to all views (ungrouped and grouped).
- Each nested row is split into two columns that maintain the `2:6` ratio relative to their parent. Since the parent column has a width of 8, the nested columns for "Perturbation" and "Effect" should have widths of `3` and `9` respectively (as `2/8 = 3/12` and `6/8 = 9/12`).
- All data content within columns should be top-aligned. For columns containing multiple horizontal elements (like the Effect column), use `d-flex align-items-start flex-wrap` to ensure content is top-aligned and wraps correctly.

# Note on grouping

Another feature of the BE which has not yet been discussed is grouping. Data can be grouped optionally by either Perturbation or Effect in the back-end. This means that a given UI side (Perturbation or Effect) becomes "vertically-merged", similarly to a Dataset.

Also, when grouping is enabled, additional metrics are displayed. In all cases, it's n_total, n_up, and n_down: how many records there are for sigificant total change of expression (total, upregulation, downregulation). Display this information in the grouped column, too.

To highlight that the columns can be grouped by, add a dotted-line, rounded-corner rectangle outline around the "Perturbation" and "Effect" `H3` headers. By default, the view should be grouped by Perturbation when the page first loads.

The grouping selector should be clickable. Clicking an unselected element selects it and unselects the other. Clicking an already-selected element unselects it, returning to the ungrouped view.

The styling of the selector is critical. The `H3` element itself should be styled with `display: 'inline-block'` to ensure the border wraps the text tightly. The default style for the selector should be `border: '1px dotted black'`, `border-radius: '8px'`, `padding: '2px'`, and `cursor: 'pointer'`. When an element is selected for grouping, it should receive an additional `background-color: '#FFF3CD'`.

## UI specifics
When grouping by perturbation:
* Display the perturbed gene name in a larger font
* Below, write in normal font: "Affects XXX genes", where number is in bold (in all following examples number also bold)
* Below: "(green arrow up) XXX genes", 
* Below: "(red arrow down) XXX genes"

When grouping by effect:
* Display the effect gene name in a larger font.
* Below, write in normal font: "Affected by XXX genes"
* Below: "(green arrow up) by XXX genes",
* Below: "(red arrow down) by XXX genes"
* Below: "Base mean expression is XX.XX", number also in bold

Then, display the perturbation rows. Note that perturbed gene should be displayed in the perturbation column, and log2fc, padj must be displayed in the Effect column. These should be well aligned.
