# --------------------------------------------------------------------------------------------------------
# Script: 01_peptide_preprocess.R
# Description: Metaproteomics data preprocessing at peptide level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------
library(here)

# -----------------------------
# Objetive
# -----------------------------




# -----------------------------
# Load Setup and Utilities
# -----------------------------
source(here::here("code", "00_setup_environment.R"))

# -----------------------------
# Load Peptide Data
# -----------------------------
peptides <- read_metap_data("raw", "run_2025", "peptide")

# -----------------------------
# Load Project Metadata
# -----------------------------
metadata <- read.csv(here(path_metadata, "metaproteomics_MM_metadata.csv"))







