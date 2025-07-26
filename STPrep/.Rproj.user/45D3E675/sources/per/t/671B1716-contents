library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(data.table)

#data.dir <- "SampleSTData"
#heart <- Load10X_Spatial(
  #data.dir,
  #filename = "V1_Human_Heart_filtered_feature_bc_matrix.h5",
  #assay = "Spatial",
  #slice = "slice1",
  #bin.size = NULL,
  #filter.matrix = TRUE,
  #to.upper = FALSE,
  #image = NULL,
#)

preprocess <- function(seuratObj){
  #sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
  seuratObj <- SCTransform(seuratObj, assay = "Spatial", verbose = TRUE)
  
  #preprocessing
  seuratObj <- RunPCA(seuratObj, assay = "SCT", verbose = FALSE)
  seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:30)
  seuratObj <- FindClusters(seuratObj, verbose = FALSE)
  seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:30)
  print(SpatialFeaturePlot(seuratObj, features = c("SDF4")))
  return(seuratObj)
}

preprocessNormalized <- function(seuratObj){
  #preprocessing
  seuratObj <- RunPCA(seuratObj, assay = "Spatial", verbose = TRUE)
  seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:30)
  seuratObj <- FindClusters(seuratObj, verbose = TRUE)
  seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:30)
  print(SpatialFeaturePlot(seuratObj, features = c("SDF4")))
  return(seuratObj)
}

getExpression <- function(seuratObj){
  #getting expression matrix and exporting to csv
  expression_matrix <- GetAssayData(object = seuratObj, assay = "Spatial", layer = "data")
  dense_expr <- as.matrix(expression_matrix)
  return(dense_expr)
}

getSpatial <- function(seuratObj){
  #getting spatial data and exporting
  coords <- GetTissueCoordinates(seuratObj, image = "slice1")
  return(coords)
  
}

# get most expressed genes (for testing)
getExpGenes <- function(seuratObj){
  cluster_markers <- FindAllMarkers(seuratObj, 
                                    ident.use = "orig.ident", # or "orig.ident" for samples
                                    only.pos = TRUE, # Only find up-regulated genes
                                    #min.pct = 0.5, # Minimum percentage of cells expressing the gene
                                    logfc.threshold = 1, # Minimum log fold-change in expression

  )
  head(cluster_markers)
  return(cluster_markers)
}


