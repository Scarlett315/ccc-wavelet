library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(rhdf5)
library(data.table)
source("STPrep.R")
rm(list = ls()) 


# getting data and putting it into Seurat
data.dir <- "../Data/human_breast_cancer/raw_data"
seu_raw <- Load10X_Spatial(
  data.dir,
  filename = "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  bin.size = NULL,
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

# Quality Control
seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^MT-")
VlnPlot(seu_raw, features = "percent.mt", pt.size = 0.1)

plot1 <- VlnPlot(seu_raw, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seu_raw, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

VlnPlot(seu_raw, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size = 0.1)

seu <- subset(seu_raw, subset = 
                nFeature_Spatial > 200 &
                nFeature_Spatial < 9000 &
                percent.mt < 8
)

#compare
p1 <- VlnPlot(seu_raw, features = "nFeature_Spatial") + ggtitle("before QC")
p2 <- VlnPlot(seu, features = "nFeature_Spatial")+ ggtitle("after QC")
p1 + p2



# getting info
DefaultAssay(seu) <- "Spatial"

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- ScaleData(seu , verbose = FALSE)
seu <- FindVariableFeatures(seu)
seu <- preprocessNormalized(seu)


exp <- getExpression(seu)
head(exp)

markers = getExpGenes(seu, min.pct=0)

filtered_markers <- markers %>% filter(p_val_adj < 0.05)

# Get list of gene names from filtered markers
marker_genes <- unique(filtered_markers$gene)
marker_expr <- exp[rownames(exp) %in% marker_genes, ,drop = FALSE]

marker_expr_df <- as.data.frame(as.matrix(marker_expr))
marker_expr_df$gene <- rownames(marker_expr_df)
marker_expr_df <- marker_expr_df[, c("gene", setdiff(names(marker_expr_df), "gene"))]

fwrite(marker_expr_df, "../Data/human_breast_cancer/info/expression_filtered_less.csv")

#unfiltered export
exp <- as.data.frame(as.matrix(exp))
exp$gene <- rownames(exp)
exp <- exp[, c("gene", setdiff(names(exp), "gene"))]
fwrite(exp, "../Data/human_breast_cancer/info/expression_unfiltered.csv")


#spatial
c <- getSpatial(seu)
fwrite(c, "../Data/human_breast_cancer/info/coordinates.csv", row.names=TRUE)

Layers(seu)
summary(as.vector(GetAssayData(seu, assay = "Spatial", layer = "counts")))
summary(as.vector(GetAssayData(seu, assay = "Spatial", layer = "data")))
summary(marker_expr_df[-1])
str(marker_expr_df)


VlnPlot(seu, features = "LXA4")
