# --------------------------------------------------------------------------------------------------------
# Script: 00_setup_environment.R
# Project: Metaproteomics MM Data Analasys - Enviroment set-up.
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
#
# Purpose:
#   Standard setup for reproducible analysis in this project
#   Load packages, sets paths, sources utility functions, and define global settings.
# --------------------------------------------------------------------------------------------------------

### 1. Environment setup ---------------------------------------------------------------------------------

# Clean environment
rm(list = ls())

# Reproducibility: fix a random seed
set.seed(1234)

### 2. Load required packages ----------------------------------------------------------------------------
## Libraries
pacman::p_load(
  # Loading and data manipulation.
  readr, openxlsx, janitor, naniar, skimr, reshape2, dplyr, tidyverse, forcats, vroom, here,
  
  # Omics and microbiomics data analysis.
  phyloseq, microbiome, MSnSet.utils, MSnbase, microViz, 
  
  # Statistical analysis and dimensional reduction.
  car, PCAtools, factoextra, DESeq2, limma, DEqMS, MSstats, vegan, hagis, ranacapa, ALDEx2, Maaslin2, rstatix, 
  MicrobiomeProfiler, clusterProfiler, KEGGREST,
  
  # Phylogeny and evolutive analysis,
  ape,  
  
  # Data visualization
  ggplot2, ggpubr, ggrepel, ggalt, ggdist, gghalves, gplots, ComplexHeatmap,  
  paletteer, RColorBrewer, scales, gridExtra, VennDiagram, shape, ggsci, patchwork, 
  ggstatsplot, ggVennDiagram,
  
  # Reports
  gtsummary
  )

### 3. Paths (relative to project root) ------------------------------------------------------------------
# Data
path_data_raw <- here("data", "raw")
path_data_processed <- here("data", "processed")

# Metadata
# TODO: Turn "run" into parameter
path_metadata <- here("metadata", "run_2025")

# Project-specific results
path_results <- here("results")

# Utility functions
utils_files <- list.files(here("code", "utils"), pattern = "\\.R$", full.names = TRUE)

### 4. Source custom functions ----------------------------------------------------------------------------
if (length(utils_files) > 0) {
  for (f in utils_files) {
    message("Loading utils: ", basename(f))
    source(f)
    }
  } else {
  warning("No utils scripts found: project functions unavailable.")
}

### 5. Global plotting options
theme_set(
  theme_bw() +
    theme(
      title = element_text(size = 24),
      legend.title = element_blank(),
      legend.text = element_text(size = 22 ),
      axis.text = element_text(size = 22, colour = "black" ),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title.x = element_text(size = 22, colour = "black", face = "plain"),
      axis.title.y = element_text(size = 22, colour = "black", face = "plain"),
    ) 
)

# Define a default palette for consistency
palette_mm_status <- c(
  "NDMM" = "red",
  "RRMM" = "blueviolet"
 )

# Proteomics metrics
metrics_by_level <- list(
  peptide = c(
    intens = "intensity",
    max_intens = "max_lfq_intensity",
    spc = "spectral_count"
  ),
  protein = c(
    intens = "intensity",
    max_intens = "max_lfq_intensity",
    spc = "spectral_count",
    uniq_spc = "unique_spectral_count",
    total_spc = "total_spectral_count"
  )
)

metrics_to_use <- c("intens", "spc")

prev_threshold <- 0.5

# Current run
# TODO: Evalute if create a config file
current_run <- "run_2025"


                  