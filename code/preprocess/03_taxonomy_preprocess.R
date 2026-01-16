# --------------------------------------------------------------------------------------------------------
# Script: 03_taxonomy_preprocess.R
# Description: Metaproteomics data preprocessing at taxonomy level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------

# -----------------------------
# Objective
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


# -----------------------------
# 5. Add taxa IDs.
# -----------------------------
taxonomy_processed <- taxonomy_processed %>% 
  mutate(
    taxa_id = paste0("taxa_", row_number()),
    .before = name
) 

# -----------------------------
# 6. Intensity-derived metrics
# -----------------------------
taxonomy_processed <- taxonomy_processed %>%
  
  mutate(
    # Sum intensity per each peptide. 
    intensity_sum = rowSums(across(matches("^ID_\\d+_intensity")), na.rm = T), 
    
    # Intensity equal to 0.
    intensity_sum_0 = intensity_sum == 0
  ) 


# -----------------------------
# 7. Prevalence-derived annotations
# -----------------------------
taxonomy_processed <- taxonomy_processed %>%
  
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
      metric = metrics_to_use_func,
      min_prop = prev_threshold
    )
  ) %>% 
  
  # Peptide all/nothing selection.
  mutate(
    !!!filter_all_nothing(
      df = .,
      cond_list = cond_list,
      metric = metrics_to_use_func,
      min_prop = prev_threshold
    )
  )

# -----------------------------
# 8. Taxonomic annotations (atomic biological flags)
# -----------------------------
taxonomy_processed <- taxonomy_processed %>%
  mutate(
    is_bacteria  = superkingdom == "Bacteria",
    is_eukaryota = superkingdom == "Eukaryota",
    is_virus     = superkingdom == "Viruses",
    is_archea   = superkingdom == "Archea",
  )

# -----------------------------
# 9. Atomic flags.
# -----------------------------
taxonomy_processed <- taxonomy_processed %>%
  # MAYBE ANY OF THESE FILTERING CRITERIA SHOULD BE IN CONFIG? RELAXED OR STRICT
  # PEPTIDE SELECTION?
  mutate(
    keep_intensity = intensity_sum > 0,
    keep_min_prev = feature_min_prev_intens,
    keep_all_nothing_NDMM = NDMM_all_nothing_intens,
    keep_all_nothing_RRMM = RRMM_all_nothing_intens,
    keep_bacteria = is_bacteria,
    keep_eukaryota = is_eukaryota,
    keep_virus = is_virus,
    keep_archea = is_archea
  )

# -----------------------------
# 10. Build final biological sets.
# -----------------------------
taxonomy_processed <- build_sets(taxonomy_processed, taxonomy_sets_defs)


# -----------------------------
# 11. Create specific objects
# -----------------------------

# Taxonomy sets.
taxonomy_sets <- taxonomy_processed %>% select(taxa_id, taxa_core:taxa_RRMM_only)

# Taxonomy table
taxonomy_tax_table <- taxonomy_processed %>% 
  # TODO: SELECT VARIABLES TO RETAIN
  select(-starts_with("keep")) %>%
  mutate(
    across(
      .cols = tax_levels, 
      .fns = ~ if_else(
        is.na(.x),
        paste0(tax_prefixes[which(tax_levels == cur_column())], tax_miss_pref),
        paste0(tax_prefixes[which(tax_levels == cur_column())], .x)
      )
    )
  ) %>%
  tidyr::unite(
    "taxa_name", 
    all_of(tax_levels), 
    sep = ";", 
    remove = F)


# ASV table
taxonomy_asv_table <- taxonomy_processed %>% 
  select(taxa_id, matches(samples)) %>%
  # TODO: Update in the future if more metrics are added at taxonomical level.
  rename_with(~ gsub("_intensity$", "", .x), ends_with("_intensity"))


  

