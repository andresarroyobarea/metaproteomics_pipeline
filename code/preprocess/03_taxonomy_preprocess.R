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
# 5. Add taxa IDs and Turn NA into "Unclassified"
# -----------------------------
taxonomy_processed <- taxonomy_processed %>% 
  mutate(
    feature_id = paste0("taxa_", row_number()),
    .before = name
) %>%
  mutate(across(superkingdom:species, ~replace_na(., tax_miss_pref)))

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
    
    # Domain
    is_bacteria  = superkingdom == "Bacteria",
    is_eukaryota = superkingdom == "Eukaryota",
    is_virus     = superkingdom == "Viruses",
    is_archaea   = superkingdom == "Archaea",
    is_prokaryota  = superkingdom %in% c("Bacteria", "Archaea"),
    
    # Taxonomic level classification flags
    # Flags por nivel taxonómico
    is_phylum_classified  = phylum  != "Unclassified",
    is_class_classified   = class   != "Unclassified",
    is_order_classified   = order   != "Unclassified",
    is_family_classified  = family  != "Unclassified",
    is_genus_classified   = genus   != "Unclassified",
    is_species_classified = species != "Unclassified",
    
    # General classification flag: From phylum to species
    is_any_tax_level_classified = phylum  != "Unclassified" |
      class   != "Unclassified" |
      order   != "Unclassified" |
      family  != "Unclassified" |
      genus   != "Unclassified" |
      species != "Unclassified"
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
    keep_archaea = is_archaea,
    keep_prokaryota = is_prokaryota,
    keep_phylum_classified = is_phylum_classified,
    keep_class_classified = is_class_classified,
    keep_order_classified = is_order_classified,
    keep_family_classified = is_family_classified,
    keep_genus_classified = is_genus_classified,
    keep_species_classified = is_species_classified,
    keep_classified_any = is_any_tax_level_classified
  )

# -----------------------------
# 10. Build final biological sets.
# -----------------------------
taxonomy_processed <- build_sets(taxonomy_processed, taxonomy_sets_defs)

# -----------------------------
# 11. Create specific objects
# -----------------------------

# TODO: Think about exporting two separte objects for taxonomy and ASVs instead of
# exporting a complete one
taxonomy_processed <- taxonomy_processed %>% 
  # TODO: SELECT VARIABLES TO RETAIN
  select(-starts_with("is")) %>%
  mutate(
    across(
      .cols = tax_levels, 
      .fns = ~paste0(tax_prefixes[which(tax_levels == cur_column())], .x)
      )
    )%>%
  tidyr::unite(
    "taxa_name", 
    all_of(tax_levels), 
    sep = ";", 
    remove = F) %>% 
  # TODO: Update in the future if more metrics are added at taxonomical level.
  rename_with(~ gsub("_intensity$", "", .x), starts_with("ID_")) 


# -------------------------------------
# 12. Export preprocessed proteins file
# -------------------------------------
write_metap_data(taxonomy_processed, state = "processed", run = current_run, level = "taxonomy")