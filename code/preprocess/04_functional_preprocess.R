# --------------------------------------------------------------------------------------------------------
# Script: 04_functional_preprocess.R
# Description: Metaproteomics data preprocessing at functional level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------

# -----------------------------
# Objetive
# -----------------------------

# The functional information is a dataframe with one row per protein with funcional information (KEGG, COG, NOG...)
# added in different variable. Protein id from functional data can be mapped with functional IDs from protein data
# so we can recycle conditions and flags from the protein information to create the same variables.


# -----------------------------
# 1. Load Functional Data
# -----------------------------
functional <- read_metap_data("raw", current_run, "functional")

# -----------------------------
# 2. Clean names
# -----------------------------
# TODO: This should be an external step where names should be properly stablished because this patters could change
# between runs.

colnames(functional) <- gsub(
  intensity_pattern, 
  intensity_replacement, 
  janitor::make_clean_names(colnames(functional))
)

# -----------------------------
# 3. Filter samples
# -----------------------------

# Keep only samples retained after metadata-level QC
functional_processed <- filter_samples(functional, samples_to_include = samples, verbose = T)


# -----------------------------
# 4. Validate metrics. 
# -----------------------------

# Ensure functional metrics are consistent with protein-level definitions
metrics <- validate_metrics(
  df = proteins_processed,
  data_level = "functional",
  metrics_requested = metrics_to_use_func,
  metrics_by_level = metrics_by_level
  )


# -----------------------------
# 5. Variables generation and flags according to proteins informatino.
# -----------------------------
# Project protein-level flags onto the functional table.
# Each functional protein inherits the flag if it belongs to the
# corresponding protein-level set.
functional_processed <- functional_processed %>%
  {fp <- . 
    # Add a variable per flag.
    bind_cols(
      
      fp,
      
      map_dfc(protein_flags_for_functions, function(flag){
        
        # Proteins satisfying the given protein-level condition
        prot_flag <- proteins_processed %>%
          filter(.data[[flag]]) %>%
          pull(protein)
        
        # TRUE if this functional protein belongs to that set  
        tibble(!!flag := fp$protein_id %in% prot_flag)
    })
    )
}


# -----------------------------
# 6. Build final biological sets.
# -----------------------------
# Functional sets mirror protein-level biological definitions
functional_processed <- build_sets(functional_processed, protein_sets_defs)

# -------------------------------------
# 7.Export preprocessed functions file
# -------------------------------------
write_metap_data(functional_processed, state = "processed", run = current_run, level = "functional")
