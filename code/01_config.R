# --------------------------------------------------------------------------------------------------------
# Script: 01_config.R
# Project: Metaproteomics MM Data Analysis - Config file
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
#
# Purpose:
#   - Load and filter sample metadata for downstream analysis
#   - Define included samples and conditions
#   - Generate condition-based sample groupings
#
# Outputs:
#   - samples        : character vector of included sample IDs
#   - metadata_filt  : filtered metadata dataframe
#   - cond_list      : named list of sample IDs by condition
# --------------------------------------------------------------------------------------------------------


# -----------------------------
# Metadata
# -----------------------------
metadata <- read.csv(
  here(path_metadata, "metaproteomics_MM_metadata.csv"), 
  header = T,
  stringsAsFactors = FALSE
  )


stopifnot(
  all(c("sample_id", "condition", "include") %in% colnames(metadata))
)

message("Metadata loaded: ", nrow(metadata), " samples")

# -----------------------------
# Filtering strategy
# -----------------------------
# Filter only selected samples for downstream analysis
samples <- metadata %>% 
  filter(include == "YES") %>% 
  pull(sample_id)

metadata_filt <- metadata %>% filter(sample_id %in% samples)

excluded_samples <- setdiff(metadata$sample_id, samples)

if (length(excluded_samples) > 0) {
  message("Excluded samples (include == NO): ", paste(excluded_samples, collapse = ", "))
  } else {
  message("No samples excluded based on 'include' flag")
}

# 2. Filter by condition.
# TODO:
# This is intentionally kept minimal.
# In the future this can be moved to a YAML or CLI argument.

exclude_condition <- "SHIELD" 

# Conditions previous filtering
conditions_before_filt <- unique(metadata_filt$condition)

# Filtering step
metadata_filt <- metadata_filt %>% filter(condition != exclude_condition)

# Condition after filtering
conditions_after_filt <- unique(metadata_filt$condition)

# Removed conditions
removed_conditions <- setdiff(conditions_before_filt, conditions_after_filt)

# Message
if (length(removed_conditions) > 0) {
  message("Excluded condition: ", exclude_condition)
  } else {
  message("No condition-level filtering applied")
}

samples <- metadata_filt$sample_id

# 3. List with samples by condition.
cond_list <- metadata_filt %>%
  group_by(condition) %>%
  summarise(sample_ids = list(sample_id)) %>%
  deframe()

# -----------------------------
# Outputs senior-friendly
# -----------------------------
# 1. samples -> vector de todas las muestras incluidas
# 2. cond_list -> lista nombrada de condiciones con sus IDs
# 3. metadata_filt -> dataframe de metadata filtrado



