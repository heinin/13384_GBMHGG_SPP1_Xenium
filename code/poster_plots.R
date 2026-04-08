#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 11/30/2023
# Description: Plotting Xenium data
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(googlesheets4)
library(gprofiler2)
library(scCustomize)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import data
#==============================================================================#

# The final object
seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/Seurat/clustered_obj_w_scRNA_labels_03_04_2026.rds")

seurat_data$binary_response <- ifelse(seurat_data$MRPplus.Overall.Best.Response == "Progression Disease (PD)",
                                      "PD", "CR_SD")

FeaturePlot(seurat_data, features = c("CD3D", "CD3E"), blend = T)
FeaturePlot(seurat_data, features = c("CD44", "COL1A1"), blend = T)

corrs <- LayerData(seurat_data,
                   assay = "RNA",
                   layer = "data")

plot(as.numeric(corrs["CD3E",]), as.numeric(corrs["CD3D",]))

# Saving annotations and coordinates
#coords <- seurat_data@meta.data[,c("predicted.id", "x_centroid", "y_centroid")]
#colnames(coords) <- c("Selection","X","Y")
#
#write.csv(coords, "/tgen_labs/banovich/personal_space/hnatri/Xenium_annots.csv",
#          quote = F, row.names = F)

head(seurat_data@meta.data)

ct_sp <- DimPlot(seurat_data,
                 group.by = "predicted.id",
                 cols = c(immune_fibro_celltype_col, tumor_celltype_col),
                 reduction = "sp",
                 raster = T,
                 ncol = 1,
                 raster.dpi = c(2048, 2048),
                 pt.size = 2) +
           ggtitle("") +
           theme_classic() +
           manuscript_theme + 
           NoLegend() +
           NoAxes() &
           coord_fixed()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/Xenium_celltype_sp.pdf"
pdf(file = filename,
    width = 4,
    height = 4)

ct_sp

dev.off()


FeaturePlot(seurat_data,
            reduction = "sp",
            features = c("nCount_Spatial", "nFeature_Spatial"),
            split.by = "TMA") &
  theme_classic() &
  coord_fixed()

# GBM genes
gbm_genes <- c("MGMT", "IDH1", "IDH2", "EGFR", "PTEN", "TP53", "TERT")

FeaturePlot(seurat_data,
            reduction = "sp",
            features = gbm_genes,
            ncol = 1,
            raster = T,
            raster.dpi = c(1024, 1024)) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  theme_classic() &
  coord_fixed()

FeaturePlot_scCustom(seurat_object = seurat_data,
                     features =  c("nCount_Spatial", "nFeature_Spatial"),
                     reduction = "sp",
                     num_columns = 1,
                     raster = F,
                     #raster.dpi = c(1024, 1024),
                     slot = "data",
                     pt.size = 0.1,
                     #na_cutoff = 1,
                     order = T)

FeaturePlot_scCustom(seurat_object = seurat_data,
                     features = c("SPP1"),
                     reduction = "sp",
                     num_columns = 1,
                     raster = T,
                     raster.dpi = c(2048, 2048),
                     slot = "data",
                     pt.size = 1,
                     #na_cutoff = 20,
                     order = T)

resp_cols <- c("PD" = "deeppink3",
               "CR_SD" = "aquamarine3")

VlnPlot(seurat_data,
        features = c("SPP1", "FN1", "CD44"),
        group.by = "binary_response",
        pt.size = 0,
        cols = resp_cols,
        ncol = 4) &
  theme_classic() &
  xlab("") &
  NoLegend()

barp <- create_barplot(seurat_data,
               group_var = "binary_response",
               plot_var = "predicted.id",
               plot_levels = sort(unique(seurat_data$predicted.id)),
               group_levels = sort(unique(seurat_data$binary_response)),
               plot_colors = c(immune_fibro_celltype_col, tumor_celltype_col),
               var_names =  c("Frequency (%)", ""),
               legend_title = "") + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/Xenium_ctprop_barplot_presentation.pdf"
pdf(file = filename,
    width = 3,
    height = 4)

barp

dev.off()

# Using scProportionTest
prop_test <- sc_utils(seurat_data)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.id",
  sample_1 = "CR_SD", sample_2 = "PD",
  sample_identity = "binary_response")

perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "gray89")) + NoLegend()

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/Xenium_scProptest.pdf"
pdf(file = filename,
    width = 3,
    height = 5)

perm_plot + scale_colour_manual(values = c("tomato", "gray89")) + NoLegend() + ylab("") + xlab("")

dev.off()

# TME celltype marker dotplot
plot_features <- c("PTMA", "PFN1", "CFL1", "TMSB4X", "TPT1", "TMSB10", "MIF",
                   "PDPN", "NLRP3", "IL1B", "CCL4", "S100A8", "S100A9",
                   "S100A10", "TYROBP", "CD68", "ICAM1", "C1QA", "C1QB", "C1QC",
                   "CD74", "AREG", "CD4", "APOE", "FABP5", "SPP1", "CD274",
                   "CD96", "PTPRC", "CEMIP2", "KLRD1", "CD8A", "NKG7", "IL32",
                   "CD3D", "BTG1", "IFITM2", "ITM2A", "SELL", "GZMB", "CD79A",
                   "ACTA2", "PDGFRB", "COL1A1", "CD163", "MRC1", "ITGAM", "CD14",
                   "CD279", "PDCD1", "TREM2", "TMEM119", "P2RY12", "CX3CR1", "CD44")

# CPVL, FLT3, CLEC9A, C1orf54, CSF1R, CD1C, FCER1A, LTB, LAG3, BTLA, TIGIT,
# HAVCR2, ICOS, IL10RA, CD80, CD86, FOXP3, TMSB4X, SELL, LY6C2 (?), CD160, GZMB,
# ROCR (?), AREG, KLRG1, CDH1 (E-Cadherin), TACR1 (NK1.1, NK1R1), NKX1-2 (NK1.2),
# NCAM (CD56), KLRB1 (CD161), CD27 (TNFRSF7), TNFSF7 (?), CD74, MIF, IL1-B
additional_features <- c("CPVL", "FLT3", "CLEC9A", "C1orf54", "CSF1R", "CD1C",
                         "FCER1A", "LTB", "LAG3", "BTLA", "TIGIT", "HAVCR2",
                         "ICOS", "IL10RA", "CD80", "CD86", "FOXP3", "TMSB4X",
                         "SELL", "CD160", "GZMB", "AREG", "KLRG1", "CDH1",
                         "TACR1", "KLRB1", "CD27", "CD74", "MIF", "IL1B")

# CD15 (FUT4),CD45 (PTPRC),CD66B (not in our data), CD62L (SELL), Arg1 and PDL-1 (CD274)
# Dotplot
immune_fibro <- subset(seurat_data, subset = predicted.id %in% names(immune_fibro_celltype_col))
immune_fibro_dotplot <- DotPlot(immune_fibro,
                                features = plot_features,
                                group.by = "predicted.id",
                                cols = c("gray89", "tomato3"))

immune_fibro_dotplot <- create_dotplot_heatmap_horizontal(seurat_object = immune_fibro,
                                                          plot_features = sort(plot_features),
                                                          group_var = "predicted.id",
                                                          group_colors = immune_fibro_celltype_col,
                                                          column_title = "",
                                                          #row_km = 5,
                                                          col.order = sort(plot_features))

filename <- "/home/hnatri/Xenium_TME_dotplot.pdf"
pdf(file = filename,
    width = 10,
    height = 4)

immune_fibro_dotplot

dev.off()

# Quick niche analysis
coords <- seurat_data@meta.data[, c("adj_x_centroid", "adj_y_centroid")]
seurat_data[["fov"]]  <- CreateFOV(coords,
                                   assay = "sp_adj",
                                   type = "centroids")

seurat_data <- BuildNicheAssay(object = seurat_data,
                               fov = "fov",
                               group.by = "predicted.id",
                               neighbors.k = 20)

DimPlot(seurat_data,
        group.by = "predicted.id",
        cols = c(immune_fibro_celltype_col, tumor_celltype_col),
        reduction = "sp",
        raster = T,
        ncol = 1,
        raster.dpi = c(2048, 2048),
        pt.size = 2) +
  ggtitle("") +
  theme_classic() +
  manuscript_theme + 
  NoLegend() +
  NoAxes() &
  coord_fixed()


