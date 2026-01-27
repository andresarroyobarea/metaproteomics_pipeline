# --------------------------------------------------------------------------------------------------------
# Script: 02_protein_preprocess.R
# Description: Metaproteomics data preprocessing at protein level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------

# -----------------------------
# Objetive
# -----------------------------

# -----------------------------
# 1. Load Protein Data
# -----------------------------
proteins <- read_metap_data("raw", current_run, "protein")


# -----------------------------
# 2. Clean names
# -----------------------------

# TODO: This should be an external step where names should be properly stablished because this patters could change
# between runs.
# Ex. Samples should follow the same patter from begging: ID_203
colnames(proteins) <- gsub("^x", "ID_", janitor::make_clean_names(colnames(proteins)))


# -----------------------------
# 3. Filter samples
# -----------------------------

# Sample-metrics filtering according to sample filtering in metadata.
proteins_processed <- filter_samples(proteins, samples_to_include = samples, verbose = T)


# -----------------------------
# 4. Validate metrics. 
# -----------------------------

# Metrics to use in protein preprocessing
metrics <- validate_metrics(df = proteins_processed,
                            data_level = "protein",
                            metrics_requested = metrics_to_use_prot,
                            metrics_by_level = metrics_by_level)


# -----------------------------
# 5.Feature-level annotations
# -----------------------------
proteins_processed <- proteins_processed %>%
  mutate(
    human_protein  = is_human_feature(protein, "HUMAN"),
    unique_protein = is_unique_feature(indistinguishable_proteins)
  )


# -----------------------------
# 6. Intensity-derived metrics
# -----------------------------
proteins_processed <- proteins_processed %>%
  
  mutate(
    # Sum intensity per each peptide. 
    intensity_sum = rowSums(across(matches("^ID_\\d+_intensity")), na.rm = T), 
    
    # Intensity equal to 0.
    intensity_sum_0 = intensity_sum == 0
  ) 


# -----------------------------
# 7. Prevalence-derived annotations
# -----------------------------
proteins_processed <- proteins_processed %>%
  
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
      metric = metrics_to_use_prot,
      min_prop = prev_threshold
    )
  ) %>% 
  
  # Peptide all/nothing selection.
  mutate(
    !!!filter_all_nothing(
      df = .,
      cond_list = cond_list,
      metric = metrics_to_use_prot,
      min_prop = prev_threshold
    )
  )

# ---------------------------------------
# 8. Peptide-derived protein metrics
# ---------------------------------------

# Number of unique peptides (relaxed criteria) by protein
pep_to_prot_relaxed <- count_peptides_per_proteins(
  peptides_processed,
  peptide_flag = "peptides_core_relaxed",
  out_var = "n_unique_peptides_relaxed"
)

# Number of unique peptides (strict criteria with prev) by protein
pep_to_prot_strict <- count_peptides_per_proteins(
  peptides_processed,
  peptide_flag = "peptides_core",
  out_var = "n_unique_peptides_strict"
)

# Add peptide-level information to protein processed data.
proteins_processed <- proteins_processed %>%
  left_join(pep_to_prot_relaxed, by = "protein") %>%
  left_join(pep_to_prot_strict, by = "protein") %>%
  mutate(
    n_unique_peptides_relaxed = replace_na(n_unique_peptides_relaxed, 0),
    n_unique_peptides_strict = replace_na(n_unique_peptides_strict, 0)
  )


# -----------------------------
# 9. Atomic flags.
# -----------------------------
proteins_processed <- proteins_processed %>%
  # MAYBE ANY OF THESE FILTERING CRITERIA SHOULD BE IN CONFIG? RELAXED OR STRICT
  # PEPTIDE SELECTION?
  mutate(
    keep_non_human = !human_protein,
    keep_unique = unique_protein,
    keep_intensity = intensity_sum > 0,
    keep_min_prev = feature_min_prev_intens,
    keep_comb_uniq_spc = combined_unique_spectral_count > 0,
    keep_comb_total_peptides = combined_total_peptides >= 2,
    keep_n_unique_peptides = n_unique_peptides_relaxed >= 1,
    keep_all_nothing_NDMM = NDMM_all_nothing_intens,
    keep_all_nothing_RRMM = RRMM_all_nothing_intens,
  ) %>%
  dplyr::rename(feature_id = protein)

# -----------------------------
# 10. Build final biological sets.
# -----------------------------
proteins_processed <- build_sets(proteins_processed, protein_sets_defs)

# -------------------------------------
# 11. Export preprocessed proteins file
# -------------------------------------
write_metap_data(proteins_processed, state = "processed", run = current_run, level = "protein")

