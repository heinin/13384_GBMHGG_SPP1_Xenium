#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/05/2025
# Description: Converting AnnData to Seurat
#==============================================================================#

#==============================================================================
# Libraries
#==============================================================================

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(ggplot2)

#==============================================================================
# Convert data
#==============================================================================

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/spatial_filtered_splitsamples.rds")
#cell_seurat_immune <- readRDS("/scratch/hnatri/PIPAC/cell_immune_subset.rds")
#cell_seurat_nonimmune <- readRDS("/scratch/hnatri/PIPAC/cell_nonimmune_subset.rds")

#seurat_data <- cell_seurat_nonimmune

# Importing each element and building the seurat object
raw_counts <- read.csv("/home/hnatri/13384_GBMHGG_SPP1_Xenium/code/RSC_latest_EDM_2025-08-06/raw_counts.csv", header = F)
log1p_counts <- read.csv("/home/hnatri/13384_GBMHGG_SPP1_Xenium/code/RSC_latest_EDM_2025-08-06/log1p_counts.csv", header = F)
obs <- read.csv("/home/hnatri/13384_GBMHGG_SPP1_Xenium/code/RSC_latest_EDM_2025-08-06/obs.csv")
pcs <- as.matrix(read.csv("/home/hnatri/13384_GBMHGG_SPP1_Xenium/code/RSC_latest_EDM_2025-08-06/pcs.csv", header = F))
umap <- as.matrix(read.csv("/home/hnatri/13384_GBMHGG_SPP1_Xenium/code/RSC_latest_EDM_2025-08-06/umap.csv", header = F))

rownames(obs) <- obs$cell_id
obs <- obs[,-which(names(obs) %in% c("X"))]

colnames(raw_counts) <- rownames(seurat_data)
rownames(raw_counts) <- obs$cell_id
colnames(log1p_counts) <- rownames(seurat_data)
rownames(log1p_counts) <- obs$cell_id
rownames(pcs) <- colnames(seurat_data)
rownames(umap) <- colnames(seurat_data)

seurat <- CreateSeuratObject(counts = t(as.matrix(raw_counts)), 
                             meta.data = obs)

# Add PCA reduction
colnames(pcs) <- paste0("PCA_", 1:100)
seurat@reductions[["pca"]] <- CreateDimReducObject(embeddings = pcs, key = "PCA_", assay = DefaultAssay(seurat))

# Add spatial coordinates
position_xy <- cbind((seurat$x_centroid)*-1,
                     (seurat$y_centroid))
rownames(position_xy) <- rownames(seurat@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
seurat[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(seurat))

DimPlot(seurat,
        group.by = "TMA",
        reduction = "sp",
        raster = T,
        label = F) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  NoLegend() +
  RotatedAxis()

# Add UMAP reduction
#umap <- data$obsm$X_umap
colnames(umap) <- paste0("UMAP_", 1:2)
#rownames(umap) <- colnames(seurat)
seurat@reductions[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(seurat))

DimPlot(seurat,
        reduction = "umap",
        group.by = "Sample",
        raster = T) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  NoLegend() +
  ggtitle("")

# Saving as Seurat
saveRDS(seurat, "/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/Seurat/spatial_clustered_NN30_PC50_Seurat.rds")
