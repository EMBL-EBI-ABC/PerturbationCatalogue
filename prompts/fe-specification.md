# Front-end specification

## Introduction

Study general instructions for LLM agents in @prompts/README.md.

Study how relevant BE parts are organised in @be/search.py. Be aware that API endpoint is configured with the $PERTURBATION_CATALOGUE_BE env variable.

Study how relevant FE parts are organised in in @fe/app.py and @fe/pages/home.py.

## Pages update

Rename existing Data Portal page to Data Portal (legacy). Only the name of the button & menu needs to change.

Add another page named "Perturbations", flag it with a conspicous yellow "New!" label on the front page, and make its URL /perturbations.

For the new page icon, use the "Whirlpool" icon.

## New page specification

### Overall layout

The page layout, top to bottom, is: title, subtitle, table header, search fields, and a data table obtained from BE according to filters.

The page content should span the entire width of the container. If the page is too wide, it might be reduced to avoid overstretching the table width (use this dynamically using Boostrap styles). In this case, the page contents should be centered, so that margins are to the left and right and are equal. 

Make sure to use Bootstrap components and styles as much as possible. Make sure that everything uses the same font (the default), is visually pleasing and coherent.

### 1. Title

A very big title, centered: "Perturbation Catalogue". h1, font-size 60px, mt-8, mb-3.

### 2. Subtitle

"A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results.", font-size 25px, mb-7.

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

The top element is h4, regular coursive, not muted, mb-0.

The bottom element is h3.

### 4. Search fields

Immediately below the headers, there are three corresponding search fields. Make sure these are nice Bootstrap search fields and not the bare, default ones.

Their descriptions (placeholders which are inside the field, in a pale colour) shoud read: "Filter by dataset metadata", "Filter by perturbed gene", "Filter by affected gene".

When either of the search fields is affected, the data table is immediately updated. Hence, no need for a separate Search button.

### 5. Data table

Below the search fields, data is displayed, as provided by the BE API.

Essentially, it's a table where 1 row = 1 perturbation from the API. The dataset metadata is a vertically-merged cell which spans all perturbations + effects for this given dataset.

Across the entire table, not add *any* vertical or horizontal cell dividers. Instead, use padding to separate contents from each other. Rows of a single dataset should be closer together, while datasets should be reasonably apart from each other.

Note also that in the BE response, the effect + perturbation values are currently all in a single object. Example of such an object:

```
{"perturbation_gene_name":"MARS1","effect_gene_name":"BRCA1","log2foldchange":0.6373134697742381,"padj":8.753715655002655e-15,"base_mean":69.076813168962}
```

Here, perturbation_gene_name should be in the "Perturbation" section, and the rest of the fields should be in the "Effect" section.

Dataset ID, perturbed gene name and affected gene name should be most visible (for example in bold) because these are the main properties of each of the columns. 

#### Dataset column representation
First, dataset ID is listed in prominent bold. Do not use a label. For example, instead of "Dataset: some_dataset_123", display "some_dataset_123" only.

Then, dataset metadata is listed, one on new line. Do not use colons, as they are visually distracting. This is bad: "Tissue type: Blood".

Instead, display labels (such as "Tissue type") in light, thin italics, and values (such as "Blood") in semibold text. Do not use a colon between them.

#### Perturbation and effect column represetntation

Label should not be used for the gene name in both columns, it should simply be displayed as the very leftmost value.

Importantly, all contents for the "Perturbation" and "Effect" sections should be all in one row to make the table visually compact. No newline breaks in those columns. In comparison, the "Datasource" information mentioned above can and should be in multiple lines, as it spans multiple perturbation + effect rows.

Depending on whether log2fc is positive (increased) or negative (decreased), display an arrow up or arrow down to the left of the value. Use not regular arrows, but ▲▼ arrows.

Similarly to dataset information, do not use colons such as "base mean: 0.288", instead again display label (base mean) in thin, pale coursive, and value in semibold font. Make sure the labels and values styling is consistent between all three main columns.

#### Value formatting

This applies to all columns.

For any values that are displayed, make sure their first letter only is capitalised, *except* for datasource ID, which 

When displaying float values, make sure that - (hyphen) is replaced with a minus sign for both negative values such as -0.13, and negative exponents of padj such as 1.03e-05.

The label comes first (in thin italics as described above), then value in semibold. Good, but not excessive spacing should be between label-value pairs for the Effect column.

### Special notes on layout

It is essential to have *no* divider lines anywhere in the table part, nor horizontal or vertical. Make sure none are present.

Because the Perturbation column only has gene name, while both Dataset and Effect columns have lots of information, this column should be less wide.

For sections 3-5 (table header, search fields, data table) you need to direct extreme attention to make sure they are propertly aligned horizontally and vertically. For example, in the Perturbation column, all elements (header, search field, data) must be all aligned horizontally. For any given dataset, the data in all three columns needs to be properly aligned vertically. The best way to make it work is to wrap *everything* into a single table.
