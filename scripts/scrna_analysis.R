# Analysis of mouse nasal mucosa after influenza infection with single cell RNA-seq
# This analysis is based on a tutorial of single-cell RNA-seq data analysis in R by Zhisong He and Barbara Treutlein: https://github.com/quadbio/scRNAseq_analysis_vignette/blob/master/Tutorial.md#preparation 
# and BINF 6110 Lecture 17 Tutorial 

# LOAD LIBRARIES ----
library(tidyverse)
library(patchwork)
# install.packages("Seurat")
library(Seurat)
# if (!requireNamespace("glmGamPoi", quietly = TRUE))
#   BiocManager::install("glmGamPoi")
library(glmGamPoi)


# LOAD DATA ----
seurat_obj <- readRDS("../data/seurat_ass4.rds")
cat(sprintf("Cells after loading: %d\n", ncol(seurat_obj)))
# QUALITY CONTROL ----
# Check for missing values
colSums(is.na(seurat_obj@meta.data) | seurat_obj@meta.data == "")

# # Remove records with missing mouse id (look for mice ids that are required for downstream batch effect corrections)
# seurat_rm <- subset(seurat_rm, subset = mouse_id %in% c("m1_RT_D5", "m2_RT_D5", "m3_RT_D5",
#                                                         "m1_RT_Naive", "m2_RT_Naive", "m3_RT_Naive"))
# cat(sprintf("Cells remaining after mouse_id filter: %d\n", ncol(seurat_rm)))

# Add percent.mt column to investigate percentage of mitochondrial DNA
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Visualize QC metrics as a violin plot
qc1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Filter data for degraded cells and potential doublets
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)
cat(sprintf("Cells after QC filtering: %d\n", ncol(seurat_obj)))

qc2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

qc_violin_plots <- (qc1 / qc2) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1))

ggsave("../figs/02_qc_violin.png", plot = qc_violin_plots, width = 12, height = 8, dpi = 300)

# NORMALIZATION ----
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
