#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/05/2025
# Description: 13384 GBM/HGG Xenium colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)

#==============================================================================
# Colors and themes
#==============================================================================

# ggplot theme
plot_theme <- theme(text = element_text(size = 10),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 10),  
                    axis.title.x = element_text(size = 10),
                    axis.title.y = element_text(size = 10))

# ggplot theme for manuscript plots
manuscript_theme <- theme(text = element_text(size = 6),
                          axis.text.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),  
                          axis.title.x = element_text(size = 6),
                          axis.title.y = element_text(size = 6))

# Colors for plotting
# Define colors for each level of categorical variables

# Cluster colors
main_clusters <- as.factor(c(0, seq(1:20)))
main_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(main_clusters))
names(main_cluster_col) <- levels(main_clusters)

immune_clusters <- as.factor(c(0, seq(1:33)))
immune_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(immune_clusters))
names(immune_cluster_col) <- levels(immune_clusters)

nonimmune_clusters <- as.factor(c(0, seq(1:29)))
nonimmune_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(nonimmune_clusters))
names(nonimmune_cluster_col) <- levels(nonimmune_clusters)

# Cell type colors
#gs4_deauth()
#metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1sXXwOreLxjMSUoPt79c6jmaQpluWkaxA5P5HfDsed3I/edit?usp=sharing")
#pipac_celltypes <- read_sheet(metadata, sheet = "Cell type colors")
#
#pipac_celltype_col <- pipac_celltypes$Color
#names(pipac_celltype_col) <- pipac_celltypes$Annotation

# Tissue
tissue_col <- c("Tumor" = "brown3",
                "Normal" = "cyan4")

