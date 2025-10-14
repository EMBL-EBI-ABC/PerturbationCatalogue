# Front-end specification

## Introduction

Study general instructions for LLM agents in @prompts/README.md.

Study how relevant BE parts are organised in @be/search.py. Be aware that API endpoint is configured with the $PERTURBATION_CATALOGUE_BE env variable.

Study how relevant FE parts are organised in in @fe/app.py and @fe/pages/home.py.

## Pages update

Rename existing Data Portal page to Data Portal (legacy). Only the name of the button & menu needs to change.

Add another page named "Perturbations", flag it with a conspicous yellow "New!" label on the front page, and make its URL /perturbations.

## New page specification

### Overall layout

The page layout, top to bottom, is: title, subtitle, table header, search fields, and a data table obtained from BE according to filters.

The page content should span the entire width of the container. If the page is too wide, it might be reduced to avoid overstretching the table width (use this dynamically using Boostrap styles). In this case, the page contents should be centered, so that margins are to the left and right and are equal. 

Make sure to use Bootstrap components and styles as much as possible. Make sure that everything uses the same font (the default), is visually pleasing and coherent.

### 1. Title

A very big title, centered: "Perturbation Catalogue". h1, font-size 48px, mt-8, mb-3.

### 2. Subtitle

"A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results". h2, font-size 25px, mb-5.

### 3. Table header

The main concept is three columns, left to right: dataset > perturbation > effect.

The column headers need to be conspicuously labeled in a large font, with smaller explanation on top. So each header is two-line: first line in a smaller font provides a grammatical connectivity and context, while the second line is in a much larger font to serve as the column header.

The three titles are:
* (According to this / Dataset)
* (Introducing this / Perturbation)
* (Leads to the following) / (Effect)

This is what it should look like on the page:
[in small font] According to this    Introducing this    Leads to the following
[in large font] Dataset..............Perturbation........Effect

Note the dotted lines between "Dataset", "Perturbation" and "Effect". They should serve as a visual connecting element, as well as serving a dual purpose of separating table header and data.

The top element is regular coursive, regular size.

The bottom element is h3.

### 4. Search fields

Immediately below the headers, there are three corresponding search fields. Make sure these are nice Bootstrap search fields and not the bare, default ones.

Their descriptions (placeholders which are inside the field, in a pale colour) shoud read: "Filter by dataset metadata", "Filter by perturbed gene", "Filter by affected gene".

When either of the search fields is affected, the data table is immediately updated. Hence, no need for a separate Search button.

### 5. Data table

Below the search fields, data is displayed, as provided by the BE API.

Essentially, it's a table where 1 row = 1 perturbation from the API. The dataset metadata is a vertically-merged cell which spans all perturbations + effects for this given dataset.

Across the entire table, not add *any* vertical or horizontal cell dividers. Instead, use padding to separate contents from each other. Rows of a single dataset should be closer together, while datasets should be reasonably apart from each other.

Even though it's a table, it might be easier to use <div>s instead <table> because their configuration is more flexible.

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

Importantly, all contents for the "Perturbation" and "Effect" sections should be all in one row to make the table visually compact. In comparison, the "Datasource" information mentioned above can and should be in multiple lines, as it spans multiple perturbation + effect rows.

Depending on whether log2fc is positive (increased) or negative (decreased), display an arrow up or arrow down to the left of the value. Use not regular arrows, but a Unicode wide up/down arrows.

Similarly to dataset information, do not use colons such as "base mean: 0.288", instead again display label (base mean) in thin, pale coursive, and value in semibold font. Make sure the labels and values styling is consistent between all three main columns.

#### Value formatting

This applies to all columns.

For any values that are displayed, make sure they are capitalised, *except* for datasource ID.

When displaying floata values, make sure that - (hyphen) is replaced with a minus sign for both negative values such as -0.13, and negative exponents of padj such as 1.03e-05.

## Final details and particulars to keep in mind
* It is essential to have *no* divider lines anywhere in the table part, nor horizontal or vertical. Make sure none are present.
* Make sure column layout is fully aligned: explanation ("According to this"), header ("Dataset"), the corresponding search field, and then the data fields should be perfectly aligned horizontally. Make sure also that the start of every dataset is aligned in all three columns by the top.
* Do not forget a break between explanation and header: "According to this" should be above "Dataset", they should not be next to each other.
* Use Bootstrap themes as much as possible, including for large title and subtitle. 
* Overall layout must be strict, exact, and professional looking.
