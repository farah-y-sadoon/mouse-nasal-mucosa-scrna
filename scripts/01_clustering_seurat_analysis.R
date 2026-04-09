# Clustering analysis of mouse nasal mucosa after influenza infection with single cell RNA-seq
# This analysis is based on the following tutorials of single-cell RNA-seq data analysis in R: 
# Zhisong He and Barbara Treutlein: https://github.com/quadbio/scRNAseq_analysis_vignette/blob/master/Tutorial.md#preparation 
# BINF 6110 Lecture 17 Tutorial 
# https://swbioinf.github.io/scRNAseqInR_Doco/clustermarkers.html
# https://swbioinf.github.io/scRNAseqInR_Doco/singler.html

# SET WORKING DIRECTORY ----
# Assumes user is in the top level of the directory
current_dir <- basename(getwd())

if (current_dir != "scripts") {
  setwd("./scripts/")
}

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
seurat_obj <- readRDS("../data/seurat_ass4.rds")
cat(sprintf("Cells after loading: %d\n", ncol(seurat_obj)))
# QUALITY CONTROL ----
# Add percent.mt column to investigate percentage of mitochondrial DNA
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Visualize QC metrics as a violin plot
qc1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Filter data for potentially degraded cells
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
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Visualize top features
top_features <- head(VariableFeatures(seurat_obj), 20)
p1 <- VariableFeaturePlot(seurat_obj)
p2 <- LabelPoints(plot = p1, points = top_features, repel = TRUE)
ggsave("../figs/03_variable_features.png", plot = p2, width = 12, height = 6, dpi = 300)

# PCA ----
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave("../figs/04_elbow_plot.png", plot = elbow_plot, width = 8, height = 6, dpi = 300)

# CLUSTERING ----
# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
head(seurat_obj$seurat_clusters)

# Create UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) + 
  ggtitle("")
umap_tissue <- DimPlot(seurat_obj, reduction = "umap", group.by = "organ_custom", label = TRUE, raster = FALSE) + 
  ggtitle("")
umap_time <- DimPlot(seurat_obj, reduction = "umap", group.by = "time", raster = FALSE) + 
  ggtitle("")
ggsave("../figs/umap_tissue.png", plot = umap_tissue, width = 10, height = 8, dpi = 300)
ggsave("../figs/umap_time.png", plot = umap_time, width = 10, height = 8, dpi = 300)
ggsave("../figs/05_umap_clusters.png", plot = umap_clusters, width = 10, height = 8, dpi = 300)

# Find cluster markers
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

write.csv(all_markers, "../output/all_markers.csv", row.names = FALSE)

# ANNOTATION ----
# Loading the mouse immune gene reference library 
mouse_imm_gen_ref <- ImmGenData()

# Create single cell experiment object
sce <- as.SingleCellExperiment(seurat_obj)

# Run SingleR to annotate with reference library
singler_results <- SingleR(test = sce, ref = mouse_imm_gen_ref, labels = mouse_imm_gen_ref$label.main)
table(singler_results$labels)

# Add SingleR labels to the Seurat object
seurat_obj$SingleR.labels <- singler_results$labels

# Look at cell types across time points and tissue types
cat("Frequency of cell types across time:")
table(seurat_obj$SingleR.labels, seurat_obj$time)

cat("Frequency of cell types across tissue types:")
table(seurat_obj$SingleR.labels, seurat_obj$organ_custom)

# Look at macrophage counts across clusters
cat("Frequency of Macrophages across clusters:")
table(seurat_obj$seurat_clusters[seurat_obj$SingleR.labels == "Macrophages"])

# UMAP with annotations
umap_annotated <- DimPlot(seurat_obj, reduction = "umap", 
                          group.by = "SingleR.labels", 
                          label = TRUE, repel = TRUE, raster = FALSE) + 
  ggtitle("")
ggsave("../figs/06_umap_annotated.png", plot = umap_annotated, width = 12, height = 8, dpi = 300)

# Look at top 20 genes (with statistical significance) in cluster 1 to see if markers line up with Macrophages
print("Top 20 genes in cluster 1 (by log2 fold change):")
all_markers %>% 
  filter(cluster == 1, p_val_adj < 0.05) %>% 
  slice_max(order_by = avg_log2FC, n = 20) %>%
  print()

# Feature plot to show cluster 1 markers (Macrophage markers)
feature_plot <- FeaturePlot(seurat_obj, features = c("Fcrls", "Trem2", "C1qa"), ncol = 3)
ggsave("../figs/07_feature_plot_macrophages.png", 
       plot = feature_plot, 
       width = 24, height = 8, dpi = 300)

# Save Processed Seurat Object
saveRDS(seurat_obj, "../data/seurat_processed.rds")
