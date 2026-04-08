#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 11/30/2023
# Description: Plotting proximity results
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

suppressPackageStartupMessages({
  library(workflowr)
  library(dplyr)
  library(plyr)
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(googlesheets4)
  library(workflowr)
  library(patchwork)
  library(spatstat)
  library(missMDA)
  library(FactoMineR)})

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

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/13384_GBMHGG_Xenium/Seurat/clustered_obj_w_scRNA_labels_03_04_2026.rds")

seurat_data$binary_response <- ifelse(seurat_data$MRPplus.Overall.Best.Response == "Progression Disease (PD)",
                                      "PD", "CR_SD")

enrichment_scores_compiled_tumor <- read.table("/home/hnatri/13384_GBMHGG_SPP1_Xenium/enrichment_tumor_scores_compiled.tsv",
                                               header = T)
proximity_compiled <- read.table("/home/hnatri/13384_GBMHGG_SPP1_Xenium/enrichment_tumor_scores_compiled.tsv",
                                 header = T)

#==============================================================================#
# Heatmap
#==============================================================================#

enrichment_scores_compiled_tumor <- enrichment_scores_compiled_tumor[complete.cases(enrichment_scores_compiled_tumor),]

pval <- reshape2::dcast(prox_to ~ prox_ct, value.var = 'pvale_adjust', data = enrichment_scores_compiled_tumor) 
or <- reshape2::dcast(prox_to ~ prox_ct, value.var = 'log_or', data = enrichment_scores_compiled_tumor)

rownames(pval) <- pval$prox_to
pval$prox_to <- NULL

rownames(or) <- or$prox_to
or$prox_to <- NULL

ct_pval<- pval
ct_pval[is.na(ct_pval)] = 1

ct_or<- or
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

ct_or[ct_pval > 0.05] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

identical(colnames(ct_or), rownames(ct_or))

hp <- Heatmap(t(ct_or),
              name = 'log(OR)',
              row_title = 'Anchor',
              row_title_side = 'right',
              row_names_side = 'left',
              row_dend_side = 'right',
              row_title_gp = gpar(fontsize = 15),
              column_title = 'Nearest Neighbor',
              column_title_side = "bottom",
              column_names_side = 'top',
              column_names_rot = 45,
              column_dend_side = 'bottom',
              column_title_gp = gpar(fontsize = 15),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              width = unit(calc_heatmap_w(ct_or), "in"),
              height = unit(calc_heatmap_h(ct_or), "in"),
              na_col = 'white',
              col = col_fun,
              row_dend_width = unit(0.5, 'cm'),
              column_dend_height = unit(0.5, 'cm'))

draw(hp)

filename <- "/home/hnatri/13384_CART/13384_Tumors/Plots/Xenium_proximity_heatmap.pdf"
pdf(file = filename,
    width = 10,
    height = 10)

draw(hp)

dev.off()

