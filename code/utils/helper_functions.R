# --------------------------------------------------------------------------------------------------------
# Script: helper_functions.R
# Project: Metaproteomics MM Data Analasys - Helper functions
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
#
# Purpose:
#   General functions which could be used in the analysis
# --------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------
#' Filter samples input files
#' 
#' Filter samples from biological datasets according to use specifications.
#'
#' @param df data.frame. 
#'    Biological dataset (eg. peptide, protein, taxonomy or functional).
#'    
#' @param samples_to_include 
#'    Character. Samples selected by user. They come from 01_config.R.
#'    
#' @param sample_pattern Character. 
#'    Initial identifier to identify sample information in Fragpipe/Imetalab files.
#'    
#' @param verbose logical, default TRUE.
#'   Whether to print informative messages about excluded or missing metrics.
#'   
#' @return tibble 
#'    Filtered requested dataset.
#'    
#' @example
#' read_metap_data(df = peptides, samples_to_include = samples, sample_pattern = "^ID_[^_]_", verbose = TRUE)

filter_samples <- function(
    df,
    samples_to_include,
    sample_pattern = "^ID_[^_]+_",
    verbose = T
) {
  
  # 1. Identify sample variables.
  sample_cols <- grep(sample_pattern, colnames(df), value = TRUE)
  
  if (length(sample_cols) == 0) {
    stop("No sample variables detected using pattern: ", sample_pattern)
  }
  
  # 2. Extract sample_id from colnames.
  # TODO: This patterns does not capture complex id (eg. 203_sh). Maybe the easiest
  # way to solve it is to create an specific IDs system with only numbers.
  sample_ids <- gsub("^(ID_[^_]+)_.*", "\\1", sample_cols)
  
  # 3. Keep variables according to selected samples.
  keep_sample_cols <- sample_cols[sample_ids %in% samples]
  
  # 4. Structural variables (peptides/protein/taxonomy/functional information)
  structural_cols <- setdiff(colnames(peptides), sample_cols)
  
  # 5. Informative checks
  missing_samples <- setdiff(samples, unique(sample_ids))
  dropped_samples <- setdiff(unique(sample_ids), samples)
  
  if (verbose) {
    if (length(missing_samples) > 0) {
      warning(
        "Samples expected but not found in data: ",
        paste(missing_samples, collapse = ",")
      )
    }
    if (length(dropped_samples) > 0) {
      message(
        "Samples removed from dataset: ",
        paste(dropped_samples, collapse = ",")
      )
    }
    
  }
  
  df_filtered <- df %>% dplyr::select(all_of(c(structural_cols, keep_sample_cols)))
  
  message("Sample-related columns filtering complete")
  
  return(df_filtered)
}

# --------------------------------------------------------------------------------------------------------
#' Function to define is a protein has a human origin.
#' 
#' #' This function checks if the protein name matches a given human pattern and
#' returns a factor with levels YES, NO, UNKNOWN.
#'
#' @param protein_var_name character. 
#'    Protein name variable retrieved from biological datasets (peptide, protein, taxonomy, functional)
#'    
#' @param human_pattern character. 
#'    Specific pattern to identify human proteins. It could change between different reference databases.
#'    
#' @return factor.
#'    YES/NO factor if satisfy the condition.
#'
#' @example 
#' is_human_protein(c("HUMAN_P12345", "NLM010_GL0045818", NA), "HUMAN")

is_human_protein <- function(protein_var_name, human_pattern) {
  
  # Ensure input is character
  protein_var_name <- as.character(protein_var_name)
  
  factor(
    case_when(
      is.na(protein_var_name) ~ "UNKNOWN",
      grepl(human_pattern, protein_var_name, ignore.case = TRUE) ~ "YES",
      TRUE ~ "NO"
    ), levels = c("YES", "NO", "UNKNOWN")
  )
  }

  
  
# --------------------------------------------------------------------------------------------------------
#' Function to define if a peptide map to only one protein
#' 
#'
#' @param mapped_proteins_var_name character. Additional proteins mapped by peptides.
#' @return YES/NO factor


is_unique_peptide <- function(mapped_protein_var_name) {
  
  case_when(
    # TODO: Initially this varible has "" instead of NA and maybe it should be modified.
    mapped_protein_var_name == "" ~ "YES",
    mapped_protein_var_name != "" ~ "NO",
    TRUE ~ NA
    
  ) %>% factor()
}

# --------------------------------------------------------------------------------------------------------
#' Function to define if a peptide maps to only one protein
#'
#' Checks the "mapped_proteins" variable. Returns YES if no additional proteins are mapped,
#' NO if the peptide maps to multiple proteins, NA if unknown/missing.
#'
#' @param mapped_protein_var_name character. Additional proteins mapped by peptides.
#' 
#' @return Factor with levels "YES", "NO".
#' 
#' @examples
#' is_unique_peptide(c("", "P12345;P67890", NA))
# --------------------------------------------------------------------------------------------------------
is_unique_peptide <- function(mapped_protein_var_name) {
  
  # Ensure input is character
  mapped_protein_var_name <- as.character(mapped_protein_var_name)
  
  # Compute uniqueness
  factor(
    case_when(
      is.na(mapped_protein_var_name) | mapped_protein_var_name == "" ~ "YES",
      TRUE ~ "NO"
    ),
    levels = c("YES", "NO")
  )
}

# --------------------------------------------------------------------------------------------------------
#' Validate requested proteomics metrics against dataset and data level
#'
#' This function checks whether the metrics requested by the user:
#'   1) Are defined for the given biological data level (e.g. peptide, protein),
#'   2) Actually exist as columns in the provided dataset.
#'
#' Only metrics that satisfy both conditions are returned. Optionally, informative
#' messages are printed to report missing or invalid metrics.
#' 
#' @param df data.frame. 
#'    Biological dataset (eg. peptide, protein, taxonomy or functional)
#'    
#' @param data_level character. 
#'    Biological information level. Must match a name in `metrics_by_level`. 
#'    
#' @param metrics_requested character vector. 
#'    Short metrics keys requested by the user (e.g. c("intens", "spc")).
#'    
#' @param metrics_by_level named list. 
#'    Registry of available metrics per data level. Each element should be a named character
#'    vector mapping short metric keys to column name suffixes.
#'    
#' @param verbose logical, default TRUE
#'    Whether to print informative messages about excluded or missing metrics.
#'    
#' @return Named character vector. 
#'    A subset of `metrics_by_level[[data_level]]` containing only metrics that
#'    were requested and are present in the dataset.
#'    
#' @examples
#' metrics_by_level <- list(
#'   peptide = c(intens = "intensity", spc = "spectral_count"),
#'   protein = c(intens = "intensity", spc = "spectral_count",
#'               uniq_spc = "unique_spectral_count")
#' )
#' 
#' validate_metrics(
#'   df = peptide_df, 
#'   data_level = "peptide", 
#'   metrics_requested = c("intens", "spc"),
#'   metrics_by_level = metrics_by_level
#' @export
#' --------------------------------------------------------------------------------------------------------

validate_metrics <- function(
    df, 
    data_level,
    metrics_requested,
    metrics_by_level,
    verbose = TRUE
) {
  
  # 1. Available metrics defined in config file.
  available_metrics <- metrics_by_level[[data_level]]
  
  if (is.null(available_metrics)) {
    stop("Unknown data level: ", data_level)
  }
  
  # 2. Valid metrics requested for an specific level
  metrics_valid <- intersect(metrics_requested, names(available_metrics))
  
  if (length(metrics_valid) == 0) {
    stop("None of the requested metrics are valid for level: ", data_level)
  }
  
  # 3. Real checking in the dataset
  metrics_present <- metrics_valid[
    sapply(
      available_metrics[metrics_valid],
      function(m)
        any(grepl(paste0("_", m, "$"), colnames(df)))
    )
  ]
  
  # 4. Identify metrics not present in the dataset
  if (verbose && length(setdiff(metrics_valid, metrics_present)) > 0) {
    message(
      "Metrics ignored (not present in dataset): ",
      paste(setdiff(metrics_valid, metrics_present), collapse = ",")
    )
  }
  
  return(available_metrics[metrics_present])
  
}

# --------------------------------------------------------------------------------------------------------
#' Count sample-level presence by condition and metric
#'
#' This function computes, for each biological entity (row), the number of samples
#' in which it is present within each condition group, separately for one or more
#' proteomics metrics.
#' 
#' #' Presence is defined as values strictly greater than a given threshold.
#' New columns are added to the dataset following the pattern:
#' 
#'  <condition>_sum_<metric>
#'  (eg. NDMM_sum_intensity, RRMM_sum_spectral_count)
#'  
#' @param df data.frame. 
#'    Biological dataset in wide format (e.g. peptide or protein table).
#'    
#' @param cond_list named list.
#'    Each element corresponds to a biological condition and contins the sample
#'    IDs associated with that condition (eg. list(NDMM = c("ID_100", "ID_203)))
#'    
#' @param metrics named character vector. 
#'    Proteomics metrics to use, mapping short names to column suffixes
#'    (e.g. c(intens = "intensity", spc = "spectral_count")).
#'
#' @param threshold numeric, default 0.
#'    Minimum value to consider a feature as present in a sample.

#' @return data.frame
#'    The input dataset with additional variables containing presence counts per 
#'    condition and metric.
#'    
#' @examples
#' metrics <- c(intens = "intensity", spc = "spectral_count")
#' 
#' peptides <- count_presence_by_condition(
#'   df = peptides,
#'   cond_list = cond_list,
#'   metrics = metrics,
#'   threshold = 0
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------

count_presence_by_condition <- function(
  df,
  cond_list,
  metrics,
  threshold = 0
) {
  
  # 1. Iterate by condition.
  cond_metric_sums <- lapply(names(cond_list), function(cond) {
    
    # Extraction condition ids.
    ids <- cond_list[[cond]]
    
    # Extract metrics variables for previous condition ids.
    lapply(names(metrics), function(metric_key) {
      
      # Dataset metric name
      metric_name <- metrics[[metric_key]]
      
      # Specific variables for metric and condition
      pattern <- paste0(
        "^(", paste(ids, collapse = "|"), ")_", metric_name, "$"
      )
      
      cols <- grep(pattern, colnames(df), value = TRUE)
      
      # Check if any variable was retrieved
      if (length(cols) == 0) {
        return(rep(NA_integer_, nrow(df)))
      }
      
      # Sum presence of each metric value in each condition.
      rowSums(df[, cols, drop = FALSE] > threshold, na.rm = T)
      
    }) %>% setNames(paste0(cond, "_feature_prev_", names(metrics)))
    
  }) %>% unlist(recursive = FALSE)
  
  return(cond_metric_sums)

}


# --------------------------------------------------------------------------------------------------------
#' Filter features by prevalence.
#'
#' Description
#' 
#' @param df data.frame. 
#'    Biological dataset in wide format (e.g. peptide or protein table).
#'    
#' @param cond_list named list.
#'    Each element corresponds to a biological condition and contins the sample
#'    IDs associated with that condition (eg. list(NDMM = c("ID_100", "ID_203)))
#'    
#' @param metric character vector. 
#'    Proteomics metric to use in the filtering process.
#'
#' @param min_prop numeric, default 0.5.
#'    Minimum prevalence value to meet in each group to retain the feature.
#' @return data.frame
#'    The input dataset with additional variables for filtering criteria.
#'    
#' @examples
#' 
#' peptides <- filter_by_min_prevalence(
#'   df = peptides,
#'   cond_list = cond_list,
#'   metrics = "intens",
#'   min.prop = 0.5
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------

filter_by_min_prevalence <- function(
  df,
  cond_list,
  metric = metrics,
  min_prop = 0.5
) {
  
  # Check parameter entities
  stopifnot(is.data.frame(df))
  stopifnot(is.list(cond_list))
  stopifnot(is.character(metric))
  stopifnot(min_prop > 0 && min_prop < 1)
  
  # Evaluate the prevalence criteria per condition and create a flag to identify
  # feature which meets the condition.
  res <- lapply(metric, function(metric) {
  
    presence_suffix <- paste0("_feature_prev_", metric)
    
    value <- if_else(
        rowSums(
          sapply(names(cond_list), function(cond) {
            
            col <- paste0(cond, presence_suffix)
            
            if (!col %in% colnames(df)) {
              stop("Column not found: ", col)
            }
            
            df[[col]] >= ceiling(length(cond_list[[cond]]) * min_prop)
          }) 
        ) == length(cond_list), 
        "YES",
        "NO")
    
    setNames(list(value),  paste0("feature_min_prev_", metric))
    
  }) %>% unlist(recursive = FALSE)

}

# --------------------------------------------------------------------------------------------------------
#' Filter all/nothing features
#'
#' Description
#' 
#' @param df data.frame. 
#'    Biological dataset in wide format (e.g. peptide or protein table).
#'    
#' @param cond_list named list.
#'    Each element corresponds to a biological condition and contins the sample
#'    IDs associated with that condition (eg. list(NDMM = c("ID_100", "ID_203)))
#'    
#' @param metric character vector. 
#'    Proteomics metric to use in the filtering process.
#'
#' @param min_prop numeric, default 0.5.
#'    Minimum prevalence value to meet in each group to retain the feature.

#' @return data.frame
#'    The input dataset with additional variables for filtering criteria.
#'    
#' @examples
#' 
#' peptides <- filter_all_nothing(
#'   df = peptides,
#'   cond_list = cond_list,
#'   metrics = "intens",
#'   min.prop = 0.5
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------


filter_all_nothing <- function(
    df, 
    cond_list, 
    metric, 
    min_prop = 0.5
    ) {
  
  # Check parameter entities
  stopifnot(is.data.frame(df))
  stopifnot(is.list(cond_list))
  stopifnot(is.character(metric))
  stopifnot(min_prop > 0 && min_prop < 1)
  
  # Evaluate the prevalence criteria per condition and create a flag to identify
  # feature which meets the condition.
  
  # Iterate per metric
  all_noth_res <- lapply(metric, function(met) {
    
    # Iterate per condition
    sapply(names(cond_list), function(cond) {
      
      # Current condition variable
      cond_col <- paste0(cond, "_feature_prev_", met)
      if (!cond_col %in% colnames(df)) 
        stop("Column not found: ", cond_col)
      
      # Remaining condition/s variable/s
      other_cols <- setdiff(paste0(names(cond_list), "_feature_prev_", met), cond_col)
      
      # "all/nothing" checks.
      if_else(
        df[[cond_col]] >= ceiling(length(cond_list[[cond]]) * min_prop) &
          rowSums(df[other_cols] != 0) == 0,
        "YES",
        "NO"
      )
      
    }, simplify = FALSE) %>% 
      
      # Nombramos cada vector con la convención deseada
      setNames(paste0(names(cond_list), "_all_nothing_", met))
    
  }) %>% 
    unlist(recursive = FALSE)
  
  return(all_noth_res)
}

# --------------------------------------------------------------------------------------------------------
#' Build biological feature sets based on logical combinations of criteria.
#'
#' @description 
#' This function generates new categorical variables ("sets") by combiniing
#' multiple precomputed logical/boolean criteria using an AND logic.
#' A feature is assigned to a set if and only if it satisfied **all**
#' the conditions defined for that set.
#' 
#' The criteria themselves must already exist as logical (TRUE/FALSE)
#' columns provided in the input dataframe. The definition of each set is
#' provided externally via a named list, allowing flexible and scalable filtering
#' strategies without hardcoding logic into the analysis pipeline.
#' 
#' @param df data.frame. 
#'   Biological dataset in wide format (e.g. peptide- or protein-level table)
#'   containing logical columns representing filtering criteria.
#'    
#' @param sets named list.
#'   A named list where each element defines a feature set.
#'   The name of each element is the name of the new variable to be created,
#'   and its value is a character vector with the column names (criteria)
#'   that must all be TRUE for a feature to belong to that set.
#'@return
#'   A data.frame identical to the input but augmented with one new factor
#'   column per defined set.
#' 
#' @example
#' feature_processed <- build_feature_sets(
#'   df = peptides,
#'   sets = peptide_sets,
#' )
#' 
#' @examples
#' peptide_sets <- list(
#'   peptides_filtered = c(
#'     "keep_human",
#'     "keep_unique",
#'     "keep_intensity",
#'     "feature_min_prev_intens"
#'   ),
#'   NDMM_all_nothing = c(
#'     "keep_human",
#'     "keep_unique",
#'     "keep_intensity",
#'     "NDMM_all_nothing_intens"
#'   )
#' )
#'
#' peptides_processed <- build_sets(
#'   df = peptides_processed,
#'   set_defs = peptide_sets
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------

build_sets <- function(df, set_defs) {
  
  # Iterate over each set definition (one output variable per set)
  for (set_name in names(set_defs)) {
    
    # Vector of column names defining the criteria for this set
    cols <- set_defs[[set_name]]
    
    # A feature belongs to the set only if all criteria are TRUE
    df[[set_name]] <- factor(
      if_else(Reduce(`&`, lapply(cols, function(c) df[[c]])), "YES", "NO")
    )
  }
  df
}






