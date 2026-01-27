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
# Clean environment
rm(list = ls())

# Reproducibility: fix a random seed
set.seed(1234)

# -----------------------------
# 1. Load libraries
# -----------------------------
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

# -----------------------------
# 2. Experiment/run info
# -----------------------------
current_run <- "run_2025"
subrun <- NULL    # opcional, si hay subrun

# -----------------------------
# 3. Paths
# -----------------------------
path_data_raw       <- here("data", "raw", current_run)
path_data_processed <- here("data", "processed", current_run)
path_metadata       <- here("metadata", current_run)
path_metadata_processed       <- here("results", current_run)
metadata_file       <- "metaproteomics_MM_metadata.csv"
path_results        <- here("results", current_run)


# -----------------------------
# 4. Metadata filtering rules
# -----------------------------
# Column to use for filtering active samples
metadata_include_col <- "include"

# Allowed values for downstream analysis
metadata_include_yes <- c("YES")
metadata_condition_allowed <- c("NDMM", "RRMM")  


# Utility functions
utils_files <- list.files(here("code", "utils"), pattern = "\\.R$", full.names = TRUE)


# ----------------------------------
# 5. Loading source custom functions
# ----------------------------------
if (length(utils_files) > 0) {
  for (f in utils_files) {
    message("Loading utils: ", basename(f))
    source(f)
    }
  } else {
  warning("No utils scripts found: project functions unavailable.")
}

# ----------------------------------
# 6. Metric setup
# ----------------------------------

# Allowed proteomics metrics
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
  ),
  functional = c(
    intens = "intensity",
    max_intens = "max_lfq_intensity"
  ),
  taxonomy = c(
    intens = "intensity"
  )
)

# Metrics to use
metrics_to_use_pep <- c("intens", "spc")
metrics_to_use_prot <- c("intens", "spc")
metrics_to_use_func <- c("intens")
metrics_to_use_tax <- c("intens")

# Prevalence threshold
prev_threshold <- 0.5


# ----------------------------------
# 7. Visualization set up
# ----------------------------------
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



                  