# Front-end specification

## Introduction

First, study general instructions for LLM agents in @prompts/README.md.

Then, study basic fe organisation in @fe/app.py and @fe/pages/home.py (you may look at other files in the fe folder as necessary).

Rename existing Data Portal page to Data Portal (legacy). Only the name of the button & menu needs to change.

Add another page named "Perturbations", flag it with a conspicous "New!" label on the front page, and make its URL /perturbations.

This page should first contain a very big title, centered: "Perturbation Catalogue"

Then, a subtitle: "A unified engine to search and filter CRISPR, MAVE and Perturb-Seq perturbation results"

The page content is a searchable + filterable data representation.

Data comes from the back-end: consult @be/search.py for the relevant API.

API endpoint is configured with the $PERTURBATION_CATALOGUE_BE env variable.

## Page header

The main concept is three columns, left to right: dataset > perturbation > effect.

The columns need to be conspicuously labeled in a large font, with explanation on top. So each header is two-line: first line in a smaller font provides a grammatical coherence, while the second line is in a much larger font to serve as the column header.

The three titles are:
* (According to this / Dataset)
* (Introducing this / Perturbation)
* (Causes the following) / (Effect)

This is how it should look like together:
According to this   Introducing this  Causes the following
DATASET.............PERTURBATION......EFFECT
(Note that the last line should not be in caps, but just in a large conspicuous font)
(Also note that a dotted lined between the Dataset, Perturbation, and Effect should be included to serve both for visually linking them into the dataset > perturbation > effect flow, and to serve as a visual table header separator)

## Search fields

Immediately below the headers, there are three corresponding search fields. Their descriptions (placeholders which are inside the field in light grey) shoud read: "Filter by dataset metadata", "Filter by perturbed gene", "Filter by affected gene".

When either of the search fields is affected, the data table is immediately updated. Hence, no need for a separate Search button.

## Data section

Below the search fields, data is displayed, as provided by the BE API.

Essentially, it's a table where 1 row = 1 perturbation from the API. The dataset metadata is a vertically-merged cell which spans all perturbations + effects for this given dataset.

Do not add either vertical or horizontal cell dividers: use padding to separate contents from each other. Rows of a single dataset should be closer together, while datasets should be reasonably apart from each other.

Even though it's a table, it might be easier to use <div>s instead <table> because their configuration is more flexible.

Note also that in the BE response, the effect + perturbation values are currently all in a single object. Example of such an object:

```
{"perturbation_gene_name":"MARS1","effect_gene_name":"BRCA1","log2foldchange":0.6373134697742381,"padj":8.753715655002655e-15,"base_mean":69.076813168962}
```

Here, perturbation_gene_name and base_mean should be in the "Perturbation" section, and the rest of the fields should be in the "Effect" section.

## Representation of data in columns

For now, there is no detailed specification of how contents of each column should look like. Use your best judgement for usability and for things to be nice visually. Some suggestions:
* Dataset ID, perturbed gene name and affected gene name should be most visible (for example in bold) because these are the main properties of each of the columns
* Depending on whether log2fc is positive (increased) or negative (decreased), display an arrow up or arrow down to the left of the value
