setwd("./Desktop/seurat-clustering/data")

# Libraries
library(Seurat)
library(patchwork)
library(tidyverse)

# Load dataset
lung_cancer.sparse = Read10X_h5(filename = "./data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")

# Check the modality of the dataset
str(lung_cancer.sparse)

# Save gene exp matrix to var counts
counts <- lung_cancer.sparse$`Gene Expression`

# Keep all the features that are expressed in at least 3 cells, keep all genes that have at least 200 genes
lung_cancer.seurat.obj <- CreateSeuratObject(counts=counts, project="NSCLC", min.cells=3, min.features=200)
str(lung_cancer.seurat.obj)

# Quality Control -- filter out low-quality cells
# A poor quality cell has a low number of genes or molecules detected.
lung_cancer.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(lung_cancer.seurat.obj, pattern="^MT-")
View(lung_cancer.seurat.obj@meta.data)

VlnPlot(lung_cancer.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(lung_cancer.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# Filtering
lung_cancer.seurat.obj <- subset(lung_cancer.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# Normalization
lung_cancer.seurat.obj <- NormalizeData(lung_cancer.seurat.obj)
str(lung_cancer.seurat.obj)

# Identify highly variable features
lung_cancer.seurat.obj <- FindVariableFeatures(lung_cancer.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(lung_cancer.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(lung_cancer.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Scaling
all.genes <- rownames(lung_cancer.seurat.obj)
lung_cancer.seurat.obj <- ScaleData(lung_cancer.seurat.obj, features = all.genes)

str(lung_cancer.seurat.obj)

# Linear dimensionality reduction
lung_cancer.seurat.obj <- RunPCA(lung_cancer.seurat.obj, features = VariableFeatures(object = lung_cancer.seurat.obj))

# Identify some PCs
print(lung_cancer.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualizing cells and features
VizDimLoadings(lung_cancer.seurat.obj, dims = 1:2, reduction = "pca")
DimPlot(lung_cancer.seurat.obj, reduction = "pca")

# DimHeatmap visualizes sources of heterogeneity; helpful for deciding PCs
DimHeatmap(lung_cancer.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(lung_cancer.seurat.obj, dims = 1:10, cells = 500, balanced = TRUE)

# Determine dimensionality -- two ways to do this:

# 1) JackStraw procedure: random permutation of a subset of the data
# and rerun PCA; construct a null-distribution of feature scores;
# take PCs with low p-value features

lung_cancer.seurat.obj <- JackStraw(lung_cancer.seurat.obj, num.replicate = 100)
lung_cancer.seurat.obj <- ScoreJackStraw(lung_cancer.seurat.obj, dims = 1:20)

# 2) Elbow plot: heuristic method
ElbowPlot(lung_cancer.seurat.obj)

# Based on this visualization, the "elbow" can be observed 
# from PC 6-10, so we can choose a number of dimensions in this range for
# clustering analysis.

# Clustering: based on the elbow plot heuristic the first 10 PCs were selected,
# erring on the higher side.
lung_cancer.seurat.obj <- FindNeighbors(lung_cancer.seurat.obj, dims = 1:10)
lung_cancer.seurat.obj <- FindClusters(lung_cancer.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(lung_cancer.seurat.obj@meta.data)

DimPlot(lung_cancer.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

# setting identity of clusters
Idents(lung_cancer.seurat.obj)
Idents(lung_cancer.seurat.obj) <- "RNA_snn_res.0.1"
Idents(lung_cancer.seurat.obj)

# Non-linear dimensionality reduction: uses the same number of PCs as 
# chosen for clustering input, determined by PCA and visualization.

# UMAP
lung_cancer.seurat.obj <- RunUMAP(lung_cancer.seurat.obj, dims = 1:10)
DimPlot(lung_cancer.seurat.obj, reduction = "umap")

# tSNE
lung_cancer.seurat.obj <- RunTSNE(lung_cancer.seurat.obj, dims = 1:10)
DimPlot(lung_cancer.seurat.obj, reduction = "tsne")

# Differential Gene Expression (DGE) Analysis

# Groups to perform DGE on. Groups 0-7 displayed here which correspond
# to our 8 clusters.
levels(lung_cancer.seurat.obj)

lung_cancer.seurat.obj.markers <- FindAllMarkers(lung_cancer.seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung_cancer.seurat.obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# FeaturePlot(lung_cancer.seurat.obj, features = c("APOE", "SNHG9", "CCL5", "IL7R", "MS4A1", "CD37", "S100A9"))
FeaturePlot(lung_cancer.seurat.obj, features = c("APOE", "SNHG9"))
            