# --------------------------------------------------------------------------------------------------------
# Script: 05_metadata_preprocess.R
# Description: Metaproteomics data preprocessing at taxonomy level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# Purpose:
#     Exclude user-specified samples and generate specific metadata subset for 
#     downstream analysis
# --------------------------------------------------------------------------------------------------------

# -----------------------------
# 1. Load Metadata
# -----------------------------
metadata <- read.csv(file.path(path_metadata, metadata_file))

# -----------------------------
# 2. Sample filtering
# -----------------------------
metadata_clean <- metadata %>%
  
  dplyr::mutate(
    
    select_samples = if_else(
      !!sym(metadata_include_col) %in% metadata_include_yes &
        condition %in% metadata_condition_allowed,
      TRUE, FALSE
    )
  ) %>%
  dplyr::filter(select_samples) %>%
  dplyr::select(sample_id, everything(), -select_samples)

# -----------------------------
# 3. Export cleaned file
# -----------------------------
write.csv(metadata_clean, file.path(path_metadata_processed, paste0("metadata_cleaned_", current_run, ".csv")))
