library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(rhdf5)
library(data.table)
source("STPrep.R")
library(aricode)

# getting data and putting it into Seurat
data.dir <- "../GSE208253/S1/raw_data"
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
  image.dir = "../GSE208253/S1/raw_data/spatial",
  filter.matrix = TRUE
)

load("../GSE208253/SeuratObjects/sample_1.Robj")
anno <- sample_1@meta.data
colnames(anno)
anno$sample_id.x=NULL
anno$cluster_annotations=NULL

cl.anno <- sample_1@meta.data
cl.anno$sample_id.x=NULL
cl.anno$pathologist_anno.x=NULL


seu <- AddMetaData(object = seu,
                          metadata = anno,
                          col.name = "pathologist_annotation")

seu <- AddMetaData(object = seu,
                   metadata = cl.anno,
                   col.name = "cluster_annotation")

# QC
seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^MT-")
VlnPlot(seu_raw, features = "percent.mt", pt.size = 0.1)

plot1 <- VlnPlot(seu_raw, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seu_raw, features = "nCount_Spatial", image.scale = "hires") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

VlnPlot(seu, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size = 0.1)

seu <- subset(seu_raw, subset = 
                nFeature_Spatial > 200 &
                nFeature_Spatial < 9000 &
                percent.mt < 6
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

#ARI
clusters <- seu$seurat_clusters
annotations <- seu$pathologist_annotation
valid_idx <- !is.na(clusters) & !is.na(annotations)

ari_score <- ARI(clusters[valid_idx], annotations[valid_idx])
print(ari_score)

SpatialDimPlot(seu,group.by = "cluster_annotation", image.scale = "hires")
DimPlot(seu,group.by = "pathologist_annotation")

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

