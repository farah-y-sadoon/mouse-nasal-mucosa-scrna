# Differential expression analysis of mouse nasal mucosa after influenza infection with single cell RNA-seq
# This analysis is based on the following tutorials of single-cell RNA-seq data analysis in R: 
# BINF 6110 Lecture 18 Tutorial

# LOAD LIBRARIES ----
library(tidyverse)
library(patchwork)
# install.packages("Seurat")
library(Seurat)
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
library(presto)
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# LOAD DATA ----

# DIFFERENTIAL EXPRESSION ANALYSIS ----
cat("Cells with missing mouse_id:", sum(seurat_obj$mouse_id == ""))

# Remove cells with missing mouse_ids
# Get cells with valid mouse_id
valid_cells <- colnames(seurat_obj)[seurat_obj$mouse_id != ""]
seurat_de <- seurat_obj[, valid_cells]

cat(sprintf("Cells after mouse_id filter: %d\n", ncol(seurat_de))