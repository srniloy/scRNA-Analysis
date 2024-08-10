
library(conflicted)
library(dplyr)
library(Seurat)
library(tidyverse)

nsclc.sparse.m <- Read10X_h5(filename = "Datasets/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5")
str(nsclc.sparse.m)
#cts <- nsclc.sparse.m$`Gene Expression`


# Initialize the seurat obj with the raw (non-normalized)

nsclc.seurat.obj <- CreateSeuratObject(counts = nsclc.sparse.m, min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj

# 1. QC
# % mt reads

nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)


VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

# filtering -------------

nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                             percent.mt < 8)

FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

# Normalizing--------------------
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)

# identify highly variable genes -----------

nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 1000)

top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features--------

plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# scaling -----------

all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)


# perform linear dimensionality reduction -----

nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))


# Visualize the PCA results


print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data


ElbowPlot(nsclc.seurat.obj)

# clustering --------

nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution

nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)


# setting identity of clusters

Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"

# non-linear dimensionality reduction ----------

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

# label can be set using 'label = TRUE' or by using 'LabelClusters' function

DimPlot(nsclc.seurat.obj, reduction = "umap", label = TRUE)







