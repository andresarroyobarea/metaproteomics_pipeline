# --------------------------------------------------------------------------------------------------------
# Script: eda_qc_checks.R
# Project: Metaproteomics MM - EDA QC Sanity Checks
# Author: Andrés Arroyo Barea
# Date: 2026-02-06
#
# Purpose:
#   Perform quick sanity checks on EDA-prepared datasets to detect
#   issues in normalization, feature selection, zero handling, and log transforms.
# --------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------
#' edad_qc_check
#'
#'  QC sanity checks for EDA-prepared datasets.
#'   
#' @param df data.frame
#'   A biological-level feature table (peptide, protein, taxonomy or functional)
#'   containing sample-level proteomics metrics.
#'    
#' @param metric_name character.
#'   Proteomics metric key to retrieve (e.g. \code{"intensity"}, \code{"spectral_count"}).
#'   The metric must be defined for the selected biological level in
#'   \code{eda_metrics}.
#'   
#' @param eda_verbose logical, default TRUE
#'   Whether to print informative messages describing applied filters,
#'   transformations and potential issues.
#' 
#' @return data.frame
#'   A stats-by-sample data frame containing \code{sample} and the statistics
#'   summary columns.
#'   
#' @examples
#' # Example 1: QC stats for proteins intensity
#' protein_data <- get_eda_data(
#'   feature_data_list = list(protein = proteins_processed),
#'   biological_level = "protein",
#'   metric_name = "intensity",
#'   filter_subset = "proteins_core",
#'   eda_mode = list(transform_mode = "quantitative",
#'                   zero_strategy = "na",
#'                   log_transform = "log2",
#'                   normalization = "median",
#'                   post_zero_strategy = "keep")
#'                   )
#'
# qc_results <- eda_qc_checks(protein_data, metric = "intensity")


eda_qc_check <- function(df, 
                          metric_name = "intensity", 
                          verbose = TRUE) {
  
  # 1. Check if dataframe exits
  if (is.null(df) || nrow(df) == 0) log_error("Empty data frame supplied.")
  
  # 2. Identify sample columns (feature_id is first)
  sample_cols <- setdiff(colnames(df), "feature_id")
  
  # 3. Statistics summary
  stats_per_sample <- df %>%
    dplyr::select(-feature_id) %>%
    dplyr::summarise(across(everything(), list(
      min = ~ min(., na.rm = TRUE),
      mean = ~ mean(., na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      max = ~ max(., na.rm = TRUE),
      sd = ~ sd(., na.rm = TRUE),
      iqr = ~ IQR(., na.rm = TRUE),
      zeros = ~ sum(. == 0, na.rm = T),
      zeros_perc = ~ round(sum(. == 0, na.rm = T)/nrow(df)* 100, 2),
      NAs = ~ sum(is.na(.)),
      NAs_per = ~ round(sum(is.na(.))/nrow(df)* 100, 2)
    ))) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("sample", "stat"),
      names_pattern = "^(.*)_(min|max|mean|median|sd|iqr|cv|zeros|zeros_perc|NAs|NAs_per)$"
    ) %>% 
    pivot_wider(
      names_from = stat,
      values_from = value
    ) %>%
    extract(
      col = sample,
      into = c("sample_id", "metric"),
      regex = "(ID_\\d+)_(.+)"
    )
  
  # 2. Verbose reporting.
  if (verbose) {
    log_info("=== QC summary ===")
    log_info(paste0("Metric: ", metric_name))
    log_info("\n# Number of features in the dataset:")
    print(nrow(df))
    log_info("\n# Sample statistics:")
    print(stats_per_sample)
  }
  
  if (any(stats_per_sample$NAs_per > 50) | any(stats_per_sample$zeros_perc > 50)) {
    log_warn("More than 50% features in some samples.")
  }
}
