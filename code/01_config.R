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
# TODO: Path to specific metaproteomic metadata file.
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

# Summary
message("Configuration summary")
message("Included samples (n = ", length(samples), ")")
message("Samples by condition:")
for (cond in names(cond_list)) {
  message(
    "  - ", cond, " (n = ", length(cond_list[[cond]]), "): ",
    paste(cond_list[[cond]], collapse = ", ")
  )
}


# 4. Filtering sets criteria

# Peptides
peptide_sets_defs <- list(
  peptides_core = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity",
    "keep_min_prev"
  ),
  peptides_core_relaxed = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity"
  ),
  peptides_NDMM_only = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity",
    "keep_all_nothing_NDMM"
  ),
  peptides_RRMM_only = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity",
    "keep_all_nothing_RRMM"
  )
)

# Proteins
protein_sets_defs <- list(
  proteins_core = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity",
    "keep_comb_uniq_spc",
    "keep_comb_total_peptides",
    "keep_n_unique_peptides",
    "keep_min_prev"
  ),
  proteins_NDMM_only = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity",
    "keep_comb_uniq_spc",
    "keep_comb_total_peptides",
    "keep_n_unique_peptides",
    "keep_all_nothing_NDMM"
  ),
  proteins_RRMM_only = c(
    "keep_non_human",
    "keep_unique",
    "keep_intensity",
    "keep_comb_uniq_spc",
    "keep_comb_total_peptides",
    "keep_n_unique_peptides",
    "keep_all_nothing_RRMM"
  )
)

# Functional
protein_flags_for_functions <- c(
  "keep_non_human",
  "keep_unique",
  "keep_intensity",
  "keep_comb_uniq_spc",
  "keep_comb_total_peptides",
  "keep_n_unique_peptides",
  "keep_min_prev",
  "keep_all_nothing_NDMM",
  "keep_all_nothing_RRMM"
)


# Taxonomy flags
# Proteins
taxonomy_sets_defs <- list(
  taxa_prev_50 = c(
    "keep_min_prev"
  ),
  taxa_core = c(
    "keep_intensity",
    "keep_prokaryota",
    "keep_classified_any",
    "keep_min_prev"
  ),
  taxa_NDMM_only = c(
    "keep_intensity",
    "keep_prokaryota",
    "keep_classified_any",
    "keep_all_nothing_NDMM"
  ),
  taxa_RRMM_only = c(
    "keep_intensity",
    "keep_prokaryota",
    "keep_classified_any",
    "keep_all_nothing_RRMM"
  )
)


# Functional constants
intensity_pattern <- "intensity_(\\d+\\w*)(?:_(max_lfq))?"
intensity_replacement <- "ID_\\1\\2_intensity"

# Taxonomy information and prefixes
tax_levels <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
tax_prefixes   <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
tax_miss_pref <- "Unclassified"

# Pipeline steps
# Preprocessing
run_peptide <- TRUE
run_protein <- TRUE
run_taxonomy <- TRUE
run_functional <- TRUE
run_metadata <- TRUE


# -----------------------------
# Outputs senior-friendly
# -----------------------------
# 1. samples -> vector de todas las muestras incluidas
# 2. cond_list -> lista nombrada de condiciones con sus IDs
# 3. metadata_filt -> dataframe de metadata filtrado

