#!/usr/bin/env Rscript

# 1. Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# 2. Check for the correct number of arguments
if (length(args) != 2) {
  stop("Usage: Rscript convert_to_h5ad.R <input_rds_path> <output_h5ad_path>\n", call. = FALSE)
}

input_rds <- args[1]
output_h5ad <- args[2]

# scCustomize needs the directory and filename split up
out_dir <- dirname(output_h5ad)
out_file <- basename(output_h5ad)

# 3. Load Packages quietly to keep logs clean
cat("[1/4] Loading required packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(scCustomize)
  library(reticulate)
})

# 4. Initialize Python
cat("[2/4] Linking Python environment...\n")
use_python("/homes/zakirov/.r-venv/bin/python", required = TRUE)

if (!py_module_available("anndata")) {
  stop("Error: The 'anndata' module was not found in the virtual environment.", call. = FALSE)
}

# 5. Load the Seurat Object
cat(sprintf("[3/4] Loading Seurat object from: %s\n", input_rds))
sobj <- readRDS(input_rds)

# 6. Convert and Save
cat(sprintf("[4/4] Converting to AnnData and saving to: %s\n", output_h5ad))
as.anndata(
  x = sobj, 
  file_path = out_dir, 
  file_name = out_file
)

cat("Success! Conversion complete.\n")
