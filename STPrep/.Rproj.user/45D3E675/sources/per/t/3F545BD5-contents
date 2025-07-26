library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(rhdf5)
library(data.table)
source("STPrep.R")

h5ls("GSE208253/rawData/S4/filtered_feature_bc_matrix.h5")

# getting data and putting it into Seurat
data.dir <- "../GSE208253/S12/raw_data"
seu <- Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  #bin.size = NULL,
  #filter.matrix = TRUE,
  #to.upper = FALSE,
  #image = NULL,
)
img <- Read10X_Image(
  image.dir = "../GSE208253/S12/raw_data/spatial",
  filter.matrix = TRUE
)

# 
DefaultAssay(seu) <- "Spatial"
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- ScaleData(seu , verbose = FALSE)
seu <- FindVariableFeatures(seu)
seu <- RunPCA(seu , npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
seu <- FindClusters(seu, verbose = TRUE)

DimPlot(seu)
#SpatialFeaturePlot(seu, features = c("SDF4"))

exp <- getExpression(seu)

markers = getExpGenes(seu)
filtered_markers <- markers %>% filter(p_val_adj < 0.05)

# Get list of gene names from filtered markers
marker_genes <- unique(filtered_markers$gene)
marker_expr <- exp[rownames(exp) %in% marker_genes, ,drop = FALSE]

marker_expr_df <- as.data.frame(as.matrix(marker_expr))
marker_expr_df$gene <- rownames(marker_expr_df)
marker_expr_df <- marker_expr_df[, c("gene", setdiff(names(marker_expr_df), "gene"))]

fwrite(marker_expr_df, "GSE208253/info/S4/expression_filtered.csv", row.names = TRUE)
#spatial
c <- getSpatial(seu)
fwrite(c, "../GSE208253/S12/info/coordinates.csv", row.names = TRUE)


#pathologist annotations from Seurat objects
obj <- load("../GSE208253/SeuratObjects/sample_6.Robj")
fwrite(sample_6@meta.data, "../GSE208253/S6/info/pathologist_annotations.csv", row.names = TRUE)

