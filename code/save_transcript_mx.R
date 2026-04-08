
### Packages and environment variables
suppressPackageStartupMessages({
  library(workflowr)
  library(arrow)
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(workflowr)})

### Environment variables and helper functions
setwd("/home/hnatri/13384_GBMHGG_SPP1_Xenium/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

### Import data
data_dir <- "/tgen_labs/banovich/xenium_run_folders/20251003__191226__COMBO_NB-TGenTMA_CB-hBrainTMA1"
tma_list <- c(
  TMA1_0069295 = "output-XETG00048__0069295__NB_TGenTMA1__20251003__191241",
  TMA1_0069754 = "output-XETG00048__0069754__CB_hBrainTMA1__20251003__191242")

# Get subdirectory names for obtaining file paths
subdir <- "output-XETG00048__0069754__CB_hBrainTMA1__20251003__191242"

# Get transcript counts and metadata
all_files <- list.files(file.path(data_dir, subdir), full.names = TRUE)
h5_files <- all_files[grep(".h5", all_files)]
transcript_files <- all_files[grep("transcripts.parquet", all_files)]
meta_files <- all_files[grep("cells.csv.gz", all_files)]

# Read in files
counts <- lapply(h5_files, Read10X_h5)

transcripts <- lapply(transcript_files, function(XX) {
  read_parquet(XX) })

metadata <- lapply(meta_files, function(XX) {
  tmp_meta <- read.delim(XX, sep = ",", colClasses = c(cell_id = "character"))
  rownames(tmp_meta) <- tmp_meta$cell_id
  tmp_meta })

# Rename files in lists
tma_ids <- names(tma_list)
names(counts) <- tma_ids
names(transcripts) <- tma_ids
names(metadata) <- tma_ids

all_transcripts <- list()
nuc_transcripts <- list()
updated_metadata <- list()
nuc_updated_metadata <- list()

# Filter out low quality transcripts
tma <- 1
all_transcripts[[tma]] <- transcripts[[tma]][transcripts[[tma]]$qv > 20, ]

# Create cell x gene dataframe
all_transcripts[[tma]] <- as.data.frame(table(all_transcripts[[tma]]$cell_id, 
                                              all_transcripts[[tma]]$feature_name))
names(all_transcripts[[tma]]) <- c("cell_id", "feature_name", "Count")
all_transcripts[[tma]] <- all_transcripts[[tma]] %>% 
  pivot_wider(names_from = "feature_name", values_from = "Count")

# Get blanks count per cell
blank_ids <- all_transcripts[[tma]]$cell_id
blank_mat <- all_transcripts[[tma]][, grep("BLANK", 
                                           colnames(all_transcripts[[tma]]))]
blank_counts <- as.data.frame(rowSums(blank_mat))
blank_counts$cell_id <- blank_ids

# Remove negative controls and convert to cell x gene matrix
all_transcripts[[tma]][, grep("NegControl", 
                              colnames(all_transcripts[[tma]]), 
                              invert = FALSE)]
all_transcripts[[tma]][, grep("BLANK", 
                              colnames(all_transcripts[[tma]]), 
                              invert = FALSE)]

all_transcripts[[tma]] <- all_transcripts[[tma]][, grep("NegControl", 
                                                        colnames(all_transcripts[[tma]]), 
                                                        invert = TRUE)]
all_transcripts[[tma]] <- all_transcripts[[tma]][, grep("BLANK", 
                                                        colnames(all_transcripts[[tma]]), 
                                                        invert = TRUE)]
keep_cells <- all_transcripts[[tma]]$cell_id
all_transcripts[[tma]] <- as.data.frame(all_transcripts[[tma]])
rownames(all_transcripts[[tma]]) <- keep_cells
all_transcripts[[tma]] <- all_transcripts[[tma]][, -1]
all_transcripts[[tma]] <- as.matrix(t(all_transcripts[[tma]]))

# Subset nuclear metadata to "cells" with transcripts that overlap nuclei
updated_metadata[[tma]] <- metadata[[tma]][metadata[[tma]]$cell_id %in% keep_cells, ]

# Add blank counts to metadata
updated_metadata[[tma]] <- full_join(updated_metadata[[tma]], blank_counts,
                                     by = "cell_id")
updated_metadata[[tma]] <- updated_metadata[[tma]] %>%
  dplyr::rename(num.blank = `rowSums(blank_mat)`)
rownames(updated_metadata[[tma]]) <- updated_metadata[[tma]]$cell_id
