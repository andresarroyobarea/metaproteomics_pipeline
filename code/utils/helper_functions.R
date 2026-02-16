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
    log_error(paste0("No sample variables detected using pattern: ", sample_pattern))
  }
  
  # 2. Extract sample_id from colnames.
  # TODO: This patterns does not capture complex id (eg. 203_sh). Maybe the easiest
  # way to solve it is to create an specific IDs system with only numbers.
  sample_ids <- gsub("^(ID_[^_]+)_.*", "\\1", sample_cols)
  
  # 3. Keep variables according to selected samples.
  keep_sample_cols <- sample_cols[sample_ids %in% samples]
  
  # 4. Structural variables (peptides/protein/taxonomy/functional information)
  structural_cols <- setdiff(colnames(df), sample_cols)
  
  # 5. Informative checks
  missing_samples <- setdiff(samples, unique(sample_ids))
  dropped_samples <- setdiff(unique(sample_ids), samples)
  
  if (verbose) {
    if (length(missing_samples) > 0) {
      log_warn(paste0(
        "Samples expected but not found in data: ",
        paste(missing_samples, collapse = ",")
      )
      )
    }
    if (length(dropped_samples) > 0) {
      log_info(paste0(
        "Samples removed from dataset: ",
        paste(dropped_samples, collapse = ",")
      )
      )
    }
    
  }
  
  df_filtered <- df %>% dplyr::select(all_of(c(structural_cols, keep_sample_cols)))
  
  log_info("Sample-related columns filtering complete")
  
  return(df_filtered)
}

# --------------------------------------------------------------------------------------------------------
#' Function to determine if a feature (peptide or protein) is human.
#' 
#' #' This function checks if the protein name matches a given human pattern and
#' returns a factor with levels YES, NO, UNKNOWN.
#'
#' @param feature_var character. 
#'    Protein name variable retrieved from biological datasets (peptide, protein, taxonomy, functional)
#'    
#' @param human_pattern character. 
#'    Regex pattern to identify human proteins
#'    
#' @return Logical vector: TRUE if human, FALSE if non-human, NA if unknown
#'
#' @example 
#' is_human_feature(c("HUMAN_P12345", "NLM010_GL0045818", NA), "HUMAN")

is_human_feature <- function(feature_var, human_pattern) {
  
  # Ensure input is character
  feature_var <- as.character(feature_var)
  
  # Capture human features
  result <- if_else(is.na(feature_var), NA,
                    grepl(human_pattern, feature_var, ignore.case = TRUE))
    
  return(result)
  }



# --------------------------------------------------------------------------------------------------------
#' Determine if a feature (peptide or protein) is unique
#'
#' @param mapped_feature_var Character vector of mapped proteins
#' 
#' @return Logical vector: TRUE if unique, FALSE if maps to multiple, NA if unknown
#' 
#' @examples
#' is_unique_peptide(c("", "P12345;P67890", NA))
# --------------------------------------------------------------------------------------------------------
is_unique_feature <- function(mapped_feature_var) {
  
  # Ensure input is character
  mapped_feature_var <- as.character(mapped_feature_var)
  
  # Unique if empty or NA
  result <- is.na(mapped_feature_var) | mapped_feature_var == ""
  
  return(result)
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
    log_error(paste0("Unknown data level: ", data_level))
  }
  
  # 2. Valid metrics requested for an specific level
  metrics_valid <- intersect(metrics_requested, names(available_metrics))
  
  if (length(metrics_valid) == 0) {
    log_error(paste0("None of the requested metrics are valid for level: ", data_level))
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
    log_info(paste0(
      "Metrics ignored (not present in dataset): ",
      paste(setdiff(metrics_valid, metrics_present), collapse = ",")
    )
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
    
    value <- rowSums(
          sapply(names(cond_list), function(cond) {
            
            col <- paste0(cond, presence_suffix)
            
            if (!col %in% colnames(df)) {
              log_error(paste0("Column not found: ", col))
            }
            
            df[[col]] >= ceiling(length(cond_list[[cond]]) * min_prop)
          }) 
        ) == length(cond_list)
    
    setNames(list(value),  paste0("feature_min_prev_", metric))
    
  }) %>% unlist(recursive = FALSE)
  
  return(res)

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

#' @return logical.  
#'    TRUE if the feature is present in this condition according to min_prop
#     AND absent in all other conditions ("all/nothing")
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
        log_error(paste0("Column not found: ", cond_col))
      
      # Remaining condition/s variable/s
      other_cols <- setdiff(paste0(names(cond_list), "_feature_prev_", met), cond_col)
      
      # "all/nothing" checks.
      value <-  df[[cond_col]] >= ceiling(length(cond_list[[cond]]) * min_prop) &
          rowSums(df[other_cols] != 0) == 0
      
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
    df[[set_name]] <- Reduce(
      `&`, lapply(cols, function(c) df[[c]])
      )
  }
  df
}

# --------------------------------------------------------------------------------------------------------
#' Count peptides per protein based on a peptide-level flag
#'
#' @description 
#' This function aggregates peptide-level information to protein level by counting
#' the number of peptides that satisfy a given logical condition (`peptide_flag`)
#' for each protein.
#' 
#' The filtering criterion must already exist as a logical (TRUE/FALSE) column
#' in the peptide dataset. This design decouples biological definitions
#' (e.g. core peptides, all/nothing peptides, relaxed filters) from the aggregation
#' logic, allowing flexible and scalable reuse across analyses
#' 
#' @param peptides_df data.frame.
#'   Peptide-level dataset (typically processed peptides).
#'   
#' @param peptide_flag character.
#'   Name of a logical column used to select peptides for aggregation.
#'
#' @param protein_var character, default = "protein"
#'   Name of the column identifying proteins in the peptide dataset.
#'
#' @param out_var character, default = "n_peptides"
#'   Name of the output column containing the peptide count per protein.
#'
#' @return data.frame
#'   A data.frame with one row per protein and a single column containing
#'   the number of peptides satisfying the specified condition.
#'
#' @examples
#' count_peptides_per_protein(
#'   peptides_processed,
#'   peptide_flag = "peptides_core",
#'   protein_var = "protein",
#'   out_var = "n_unique_peptides_core"
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
count_peptides_per_proteins <- function(
    peptides_df,
    peptide_flag,
    protein_var = "protein",
    out_var = "n_unique_peptides") {
  
  # Summarize the information.
  peptides_df_agg <- peptides_df %>%
    dplyr::filter(.data[[peptide_flag]]) %>%
    dplyr::count(.data[[protein_var]], name = out_var)
  
  return(peptides_df_agg)  
}


# --------------------------------------------------------------------------------------------------------
#' Save a ggplot object to disk optionally.
#' 
#' @description
#' This helper function wraps the process of saving a ggplot object to a specified
#' directory, creating it recursively if it does not exist. It supports a flag 
#' `save_plots` to enable/disable saving, making it ideal for previewing plots
#' during development or batch runs without writing files.
#' 
#' @param plot_obj ggplot
#'  The ggplot object to save
#'  
#' @param plot_name character
#'  The filename (without extension) to save the plot as.
#'
#' @param type_dir character
#'  The directory path where the plot should be saved. Will be created recursively
#'  if it does not exist.
#' 
#' @param save_plots logical
#'  Flag to control whether the plot is actually saved. Useful for previewing 
#'  plots in pipelines writing files. 
#'
#' @return ggplot
#'  Return the same ggplot2 invisibly, allowing it to be used in list or further
#'  manipulations.
#'
#' @details
#' The function ensures reproducibility in pipelines by:
#' - Creating output directories if missing (`recursive = TRUE`)
#' - Standardizing file naming (`.png`)
#' - Providing width and height defaults suitable for publication-quality figures
#'
#' @examples
#' p <- ggplot(mtcars, aes(mpg, hp)) + geom_point()
#' plot_and_save(p, "mtcars_scatter", "/tmp/plots", save_plots = TRUE)
#'
#' @export
# --------------------------------------------------------------------------------------------------------
plot_and_save <- function(plot_obj, plot_name, type_dir, save_plots) {
  
  if (save_plots) {
    dir.create(type_dir, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(
      filename = file.path(type_dir, paste0(plot_name, ".png")),
      plot = plot_obj,
      width = 12,
      height = 10
    )
  }
  
  return(plot_obj)
}


# --------------------------------------------------------------------------------------------------------
#' Create a grid of parameters for feature-centric EDA
#' 
#' @description
#'  This function generates a comprehensive grid of parameters to systematically
#'  perform Exploratory Data Analysis (EDA) across multiple biological levels,
#'  metrics, subsets, and transformation modes. It is designed to automate and
#'  standardize EDA at feature level.
#'  Each row of the resulting tibble represents a unique combination of:
#'    - Biological level (e.g., peptide, protein, functional, taxonomy)
#'    - Metric (e.g., intensity, spectral count)
#'    - Subset of features (e.g., core proteins)
#'    - Transformation mode (e.g., log2, log10, median_norm)
#' 
#'  This design decouples the definitions of biological levels, metrics, and
#'  preprocessing steps from the execution of EDA, enabling flexible, reproducible,
#'  and scalable analyses.
#'  
#' @param feature_data_list list
#'   A named list of feature-level datasets, where names correspond to biological levels.
#' 
#' @param eda_metrics_user list
#'   Named list defining the metrics to analyze per biological level. For example:
#'   ```r
#'   list(
#'     peptide = c(intensity = "intensity", spc = "spectral_count"),
#'     protein = c(intensity = "intensity")
#'   )
#'   ````
#' @param eda_transform_modes list
#'   Named list of transformation modes to apply during EDA (e.g., log2, log10),
#'   where each entry is a list of parameters defining the mode.
#'   
#' @param eda_levels character
#'   Character vector defining the "level" of EDA. Options: feature, sample or condition.
#'   
#' @param eda_feature_subset list
#'   Named list defining subsets of features to include per biological level. For exmample:
#'   ```r
#'   list(
#'     peptide = c(intensity = "intensity", spc = "spectral_count"),
#'     protein = c(intensity = "intensity")
#'   )
#'   ````
#'   
#' @param eda_modes_filter character or NULL
#'   Optional vector specifying a subset of transformation mode names to keep.
#'   If NULL (default), all modes in `eda_transform_modes` are used.  
#'   
#' @return tibble
#'   A tibble where each row corresponds to a unique combination of:
#'   `eda_level`, `biological_level`, `metric`, `subset`, `eda_transform_mode_name`,
#'   and `eda_transform_modes` (list-column containing the parameters for each mode).
#' 
#' @examples
#' create_eda_grid(
#'   feature_data_list = metap_feature_data,
#'   eda_metrics_user = list(protein = c(intensity = "intensity")),
#'   eda_transform_modes = list(log2 = list(log_base = 2)),
#'   eda_feature_subset = list(protein = "proteins_core")
#' )
#' 
#' @export
# --------------------------------------------------------------------------------------------------------
create_eda_grid <- function(
    feature_data_list,
    eda_metrics_user,
    eda_transform_modes,
    eda_levels = c("feature"),
    eda_feature_subset,      # ahora viene desde config
    eda_modes_filter = NULL
) {
  
  # filtrar modos si se especifica
  if (!is.null(eda_modes_filter)) {
    eda_transform_modes <- eda_transform_modes[eda_modes_filter]
  }
  
  purrr::imap_dfr(eda_metrics_user, function(metrics, level) {
    subsets <- eda_feature_subset[[level]]
    
    tidyr::expand_grid(
      eda_level        = eda_levels,
      biological_level = level,
      metric           = unname(metrics),
      subset           = subsets,
      eda_transform_mode_name    = names(eda_transform_modes)
    )
  }) %>%
    dplyr::mutate(
      eda_transform_modes = purrr::map(eda_transform_mode_name, ~ eda_transform_modes[[.x]])
    )
}

# --------------------------------------------------------------------------------------------------------
#' Aggregate proteomics metrics according to covariates
#' 
#' 
#' 
#' 
get_agg_fun <- function(aggregation = c("mean", "median")) {
  aggregation <- match.arg(aggregation)
  switch (
    aggregation,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE)
  )
}

# --------------------------------------------------------------------------------------------------------
#' Aggregate feature abundance data
#' 
#' 
#' 
#' 
aggregate_feature_abundance <- function(
    data_long, 
    group_level = c("feature", "sample", "condition"),
    agg_fun = mean,
    feature_col = "feature_id",
    sample_col = "sample_id",
    condition_col = "condition",
    cumm_abundance = TRUE) {
  
  
  group_level <- match.arg(group_level)
  
  # 1. Define groupping variables
  group_vars <- switch(group_level,
                       feature = feature_col,
                       sample = c(feature_col, sample_col),
                       condition = c(feature_col, condition_col)
                       )
    
  # 2. Determine if apply aggregation
  use_agg <- group_level %in% c("feature", "condition") 
  
  # 3. Data aggregation
  agg_data <- data_long %>%
    dplyr::group_by(dplyr::across(all_of(group_vars))) %>%
    dplyr::summarise(abundance = if (use_agg) agg_fun(value) else value, .groups = "drop")
  
  # 4. Order and ranking data
  agg_data <- agg_data %>%
    dplyr::group_by(dplyr::across(all_of(setdiff(group_vars, feature_col)))) %>%
    dplyr::arrange(desc(abundance)) %>%
    dplyr::mutate(rank = row_number()) %>%
    dplyr::ungroup()
  
  # 5. Calculate cummulative abundance if needed.
  if (cumm_abundance ) {
    agg_data <- agg_data %>%
      dplyr::group_by(dplyr::across(all_of(setdiff(group_vars, feature_col)))) %>%
      dplyr::mutate(cum_abundance = cumsum(abundance) / sum(abundance) * 100) %>%
      dplyr::ungroup()
  }
    
    return(agg_data)
}






