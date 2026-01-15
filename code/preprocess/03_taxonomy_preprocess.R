# --------------------------------------------------------------------------------------------------------
# Script: 03_taxonomy_preprocess.R
# Description: Metaproteomics data preprocessing at taxonomy level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------

# -----------------------------
# Objetive
# -----------------------------


# -----------------------------
# 1. Load Taxonomical Data
# -----------------------------
taxonomy <- read_metap_data("raw", current_run, "taxonomy")

# -----------------------------
# 2. Clean names
# -----------------------------
# TODO: This should be an external step where names should be properly stablished because this patters could change
# between runs.

colnames(taxonomy) <- gsub(
  intensity_pattern, 
  intensity_replacement, 
  janitor::make_clean_names(colnames(taxonomy))
)


# -----------------------------
# 3. Filter samples
# -----------------------------

# Keep only samples retained after metadata-level QC
taxonomy_processed <- filter_samples(taxonomy, samples_to_include = samples, verbose = T)


# -----------------------------
# 4. Validate metrics. 
# -----------------------------

# Ensure functional metrics are consistent with protein-level definitions
metrics <- validate_metrics(
  df = taxonomy_processed,
  data_level = "taxonomy",
  metrics_requested = metrics_to_use_tax,
  metrics_by_level = metrics_by_level
)







