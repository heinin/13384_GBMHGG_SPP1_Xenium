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

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/spatial_filtered_splitsamples.rds")

# Seurat v5 assay causes an error
seurat_data[["RNA"]] <- as(object = seurat_data[["RNA"]], Class = "Assay")
seurat_data <- FindVariableFeatures(seurat_data)

head(seurat_data@meta.data)

SaveH5Seurat(seurat_data, filename = "seurat_data.h5Seurat", overwrite = TRUE)
Convert("seurat_data.h5Seurat", dest = "h5ad", overwrite = TRUE)

