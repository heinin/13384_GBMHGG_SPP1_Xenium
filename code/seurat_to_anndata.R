#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 10/11/2025
# Description: Convert Seurat to AnnData with SeuratDisk
#==============================================================================#

#==============================================================================
# Import packages
#==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(anndata)})

#==============================================================================
# Environment variables and helper functions
#==============================================================================

setwd("/home/hnatri/13384_GBMHGG_SPP1_Xenium/code/RSC_latest_EDM_2025-08-06/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

#==============================================================================
# Import data and convert
#==============================================================================

seurat_data <- readRDS("/scratch/hnatri/13384_GBMHGG_SPP1_Xenium/cell_immune_subset.rds")

# Seurat v5 assay causes an error
#seurat_data[["RNA"]] <- as(object = seurat_data[["RNA"]], Class = "Assay")
#seurat_data <- FindVariableFeatures(seurat_data)

assay_v3 <- CreateAssayObject(
  counts = seurat_data[["RNA"]]$counts,
)

seurat_data[["RNA"]] <- assay_v3

head(seurat_data@meta.data)

SaveH5Seurat(seurat_data, filename = "seurat_data_immune.h5Seurat", overwrite = TRUE)
Convert("seurat_data_immune.h5Seurat", dest = "h5ad", overwrite = TRUE)

