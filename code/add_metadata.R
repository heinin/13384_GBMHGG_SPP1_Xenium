#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 01/09/2026
# Description: Adding metadata
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

suppressPackageStartupMessages({
  library(workflowr)
  library(Seurat)
  library(plyr)})

#==============================================================================#
# Import libraries
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/Seurat/immune_nonimmune_annotated_merged_Seurat_niches.rds")

#==============================================================================#
# Add metadata
#==============================================================================#

unique(seurat_data$Sample)

seurat_data$UPN <- as.numeric(gsub("UPN", "", sapply(strsplit(seurat_data$Sample,"-"), `[`, 1)))

SPP1hi <- c(146, 248, 234)

seurat_data$SPP1_high <- ifelse(seurat_data$UPN %in% SPP1hi, "SPP1high", "SPP1low")

# Saving
saveRDS(seurat_data, "/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/Seurat/immune_nonimmune_annotated_merged_Seurat_niches.rds")
