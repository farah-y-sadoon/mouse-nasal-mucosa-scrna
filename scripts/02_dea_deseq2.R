# Differential expression and gene set enrichment analysis of mouse nasal mucosa after influenza infection with single cell RNA-seq
# This analysis is based on the following tutorials of single-cell RNA-seq data analysis in R: 
# BINF 6110 Lecture 18 Tutorial
# https://support.bioconductor.org/p/9150793/ 
# BINF 6110 Assignment 2 Differential Expression Analysis Script 

# SET WORKING DIRECTORY ----
# Assumes user is in the top level of the directory
current_dir <- basename(getwd())

if (current_dir != "scripts") {
  setwd("./scripts/")
}

# LOAD LIBRARIES ----
library(tidyverse)
library(patchwork)
library(Seurat)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

# LOAD DATA ----
seurat_obj <- readRDS("../data/seurat_processed.rds")
cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))

# PREPARE DATA ----
# Subset to only include macrophage cells for GSEA
mac_obj <- subset(seurat_obj, subset = seurat_clusters == 1)
cat(sprintf("Macrophage cells: %d\n", ncol(mac_obj)))

# remove Seurat obj to free memory ()
rm(seurat_obj)
gc()

# Look at missing mouse ids
cat("Macrophages with missing mouse_id:", sum(mac_obj$mouse_id == ""))

# Remove cells with missing mouse_ids
mac_obj <- subset(mac_obj, subset = mouse_id != "")
cat("Macrophages left after filtering out missing mouse_ids:", ncol(mac_obj))

# DIFFERENTIAL EXPRESSION ANALYSIS ----

# Pseudobulk samples 
pseudo_mac <- AggregateExpression(mac_obj, assays = "RNA", return.seurat = T, group.by = c("time", "mouse_id", "organ_custom"))

# Differential expression with DESeq2
# Get counts data
counts <- GetAssayData(pseudo_mac, assay = "RNA", layer = "counts")

# Extract metadata 
meta <- pseudo_mac@meta.data

# Set reference levels - Naive as baseline
meta$time <- factor(meta$time, levels = c("Naive", "D02", "D05", "D08", "D14"))
meta$organ_custom <- factor(meta$organ_custom, levels = c("RM", "OM", "LNG"))
meta$mouse_id <- factor(meta$mouse_id)

# Create DESeq2 object and run differential analysis for each time point (reference = Naive)
dds <- DESeqDataSetFromMatrix(counts,
                              colData = meta,
                              design = ~ organ_custom + time)
dds <- DESeq(dds)
res_D02 <- lfcShrink(dds, coef = "time_D02_vs_Naive", type = "apeglm")
res_D05 <- lfcShrink(dds, coef = "time_D05_vs_Naive", type = "apeglm")
res_D08 <- lfcShrink(dds, coef = "time_D08_vs_Naive", type = "apeglm")
res_D14 <- lfcShrink(dds, coef = "time_D14_vs_Naive", type = "apeglm")

# Save the results for each time period
write.csv(as.data.frame(res_D02), "../results/de_results_D02.csv")
write.csv(as.data.frame(res_D05), "../results/de_results_D05.csv")
write.csv(as.data.frame(res_D08), "../results/de_results_D08.csv")
write.csv(as.data.frame(res_D14), "../results/de_results_D14.csv")

# Visualize with heat map
# Select genes with highest log2-fold change between naive and d14
res_clean_D14 <- na.omit(as.data.frame(res_D14))

# Filter significant genes first, then rank by fold change
top_genes <- res_clean_D14 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)

gene_names <- rownames(top_genes)

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)

# Store counts in a matrix for the heat map
mat_vsd <- assay(vsd)[gene_names, ]

# Create an annotation data frame 
annotation_df <- data.frame(Time = meta$time, Tissue = meta$organ_custom, row.names = colnames(mat_vsd))

# Define column order with Naive first
col_order <- c(grep("Naive", colnames(mat_vsd), value = TRUE), 
               grep("D02", colnames(mat_vsd), value = TRUE), 
               grep("D05", colnames(mat_vsd), value = TRUE),
               grep("D08", colnames(mat_vsd), value = TRUE),
               grep("D14", colnames(mat_vsd), value = TRUE))

# Reorder matrix and annotation together
mat_vsd_ordered <- mat_vsd[, col_order]
annotation_df <- annotation_df[col_order, , drop = FALSE]

# Define where each time point ends to create gaps in heat map
gaps <- c(9, 18, 27, 36)

# Plot heat map
png("../figs/08_heatmap_top20.png", width = 10, height = 6, units = "in", res = 600)
pheatmap(mat_vsd_ordered, 
         scale = "row",
         cluster_cols = FALSE,
         annotation_col = annotation_df,
         show_colnames = FALSE,
         gaps_col = gaps,
         annotation_names_col = FALSE, 
         color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(100))
dev.off()

# GENE SET ENRICHMENT ANALYSIS (GSEA) ----

