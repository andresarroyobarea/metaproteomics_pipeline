# --------------------------------------------------------------------------------------------------------
# Script: 01_io.R
# Project: Metaproteomics MM Data Analasys - IO functions.
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
#
# Purpose:
#   Standardized and reproducible IO functions for reading metaproteomics data.
#   Supports raw and processed data, multiple runs, and multiple biological levels.
# --------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------
#' Read metaproteomics input files
#' 
#' Reads raw or processed data at a given biological level for a specified run.
#' This function centralizes file paths and enforces reproducibility.
#'
#' @param state Character. Either "raw" or "processed".
#' @param run Character. Experimental run, e.g. "run_2024", "run_2025", "run_human".
#' @param level Character. Biological level: "peptide", "protein", "functional", or "taxonomy".
#' @param filename Optional character. If NULL, defaults to "level.tsv" (e.g., peptide.tsv).
#' @return tibble with the requested dataset.
#' @examples
#' read_metap_data(state = "raw", run = "run_2025", level = "peptide")
#' read_metap_data(state = "processed", run = "run_2024", level = "protein", filename = "custom_proteins.tsv")

read_metap_data <- function(state, 
                            run, 
                            level, 
                            filename = NULL
                            ) {
  
  # --- Validate inputs --- #
  valid_states <- c("raw", "processed")
  if (!state %in% valid_states){
    stop(paste0("'state' muste be one of: ", paste(valid_states, collapse = ",")))
  }
  
  valid_levels <- c("peptide", "protein", "taxonomy", "functional")
  if (!level %in% valid_levels) stop(
    paste0("'level' muste be: ", paste(valid_levels, collapse = ","))
  )
  
  
  # --- Defalut filename --- #
  if (is.null(filename)) {
    # Default filenames
    filename <- switch(level,
                       "peptide" = "combined_peptide.tsv",
                       "protein" = "combined_protein.tsv",
                       "taxonomy" = "taxonomy.tsv",
                       "functional" = "functional.tsv"
    )
  }
  # --- Build path --- #
  path <- here("data", state, run, level, filename)
  
  # Check if file exists
  if (!file.exists(path)) stop(paste0("Error: File does not exist: ", path))
  
  # --- Read tsv files --- #
  df <- vroom(path, 
              delim = "\t", 
              col_names = TRUE)
  
  return(df)
}






