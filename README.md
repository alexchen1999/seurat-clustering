# Non-Small Cell Lung Cancer scRNA Data Visualization and Clustering Analysis

This project implements the Seurat standard workflow following the tutorial on the [Satija lab website](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). 

The dataset used was the CellRanger 6.1.2 [20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1 (without intronic reads)](https://www.10xgenomics.com/resources/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0).

The dataset was loaded into a SeuratObject, and low quality cells and genes were filtered out. PCA was run and 10 principal components were selected for UMAP/tSNE visualization.
