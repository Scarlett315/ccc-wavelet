library(sceasy)
library(reticulate)
library(Seurat)
library(SeuratDisk)
library(zellkonverter)
library(SingleCellExperiment)
source("STPrep.R")

use_condaenv('/opt/anaconda3/envs/sceasy_env' , required = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

loompy <- reticulate::import('loompy')
convertFormat("../Data/human_lung/human_lung.h5ad", from="anndata", to="seurat", outFile="human_lung.h5seurat")

Sys.setenv(HDF5_PLUGIN_PATH = "")
file.exists("../Data/human_lung/human_lung.h5seurat")


sce <- readH5AD("../Data/human_lung/human_lung.h5ad")

assayNames(sce)  # Check all available assays
assay(sce, "counts") <- assay(sce, "X")
assay(sce, "logcounts") <- assay(sce, "X")

lung <- as.Seurat(sce)
head(lung@meta.data)
samples <- SplitObject(lung, split.by = "Donor")
sample1 <- samples[["A32"]]
head(sample1)
