# --------------------------------------------------------------------------------------------------------
# Script: 01_peptide_preprocess.R
# Description: Metaproteomics data preprocessing at peptide level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------

# -----------------------------
# Objetive
# -----------------------------

# -----------------------------
# 1. Load Peptide Data
# -----------------------------
peptides <- read_metap_data("raw", current_run, "peptide")

# -----------------------------
# 2. Clean names
# -----------------------------

# TODO: This should be an external step where names should be properly stablished because this patters could change
# between runs.
# Ex. Samples should follow the same patter from begging: ID_203
colnames(peptides) <-  gsub("^x", "ID_", janitor::make_clean_names(colnames(peptides)))

# -----------------------------
# 3. Filter samples
# -----------------------------

# Sample-metrics filtering according to sample filtering in metadata.
peptides_processed <- filter_samples(peptides, samples_to_include = samples, verbose = T)


# -----------------------------
# 4. Validate metrics. 
# -----------------------------

# Metrics to use in preprocessing
metrics <- validate_metrics(df = peptides_processed,
                 data_level = "peptide",
                 metrics_requested = metrics_to_use_pep,
                 metrics_by_level = metrics_by_level)


# -----------------------------
# 5.Feature-level annotations
# -----------------------------
peptides_processed <- peptides_processed %>%
  mutate(
    human_protein  = is_human_feature(protein, "HUMAN"),
    unique_peptide = is_unique_feature(mapped_proteins)
  ) %>%
  dplyr::rename(feature_id = peptide_sequence)

# -----------------------------
# 6. Intensity-derived metrics
# -----------------------------
peptides_processed <- peptides_processed %>%
  
  mutate(
     # Sum intensity per each peptide. 
    intensity_sum = rowSums(across(matches("^ID_\\d+_intensity")), na.rm = T), 
    
    # Intensity equal to 0.
    intensity_sum_0 = intensity_sum == 0
    ) 

# -----------------------------
# 7. Prevalence-derived annotations
# -----------------------------
peptides_processed <- peptides_processed %>%
  
  # Peptide prevalence by condition.
  mutate(!!!count_presence_by_condition(
      df = .,
      cond_list,
      metrics = metrics
    )) %>%
  
  # Peptide prevalence filtering
  mutate(
    !!!filter_by_min_prevalence(
      df = .,
      cond_list = cond_list,
      metric = metrics_to_use_pep,
      min_prop = prev_threshold
    )
  ) %>% 
  
  # Peptide all/nothing selection.
  mutate(
    !!!filter_all_nothing(
      df = .,
      cond_list = cond_list,
      metric = metrics_to_use_pep,
      min_prop = prev_threshold
    )
  )

# -----------------------------
# 8. Atomic flags.
# -----------------------------
peptides_processed <- peptides_processed %>%
  mutate(
    keep_non_human = !human_protein,
    keep_unique = unique_peptide,
    keep_intensity = intensity_sum > 0,
    keep_min_prev = feature_min_prev_intens,
    keep_all_nothing_NDMM = NDMM_all_nothing_intens,
    keep_all_nothing_RRMM = RRMM_all_nothing_intens,
  )

# -----------------------------
# 9. Build final biological sets.
# -----------------------------
peptides_processed <- build_sets(peptides_processed, peptide_sets_defs)

# -------------------------------------
# 10. Export preprocessed peptide file
# -------------------------------------
write_metap_data(peptides_processed, state = "processed", run = current_run, level = "peptide")
