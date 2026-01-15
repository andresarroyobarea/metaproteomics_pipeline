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
#' @param state Character. 
#'    Either "raw" or "processed".
#'    
#' @param run Character. 
#'    Experimental run, e.g. "run_2024", "run_2025", "run_human".
#'    
#' @param level Character. 
#'    Biological level: "peptide", "protein", "functional", or "taxonomy".
#'    
#' @param filename Optional character. 
#'    If NULL, defaults to "level.tsv" (e.g., peptide.tsv).
#'
#' @param verbose logical, default FALSE
#'   Whether to print informative messages data importing.
#'    
#' @return tibble with the requested dataset.
#' 
#' @examples
#' read_metap_data(state = "raw", run = "run_2025", level = "peptide")
#' read_metap_data(state = "processed", run = "run_2024", level = "protein", filename = "custom_proteins.tsv")

read_metap_data <- function(state, 
                            run, 
                            level, 
                            filename = NULL,
                            verbose = FALSE
                            ) {
  
  # 1. Validate inputs
  valid_states <- c("raw", "processed")
  if (!state %in% valid_states){
    stop(paste0("'state' muste be one of: ", paste(valid_states, collapse = ",")))
  }
  
  valid_levels <- c("peptide", "protein", "taxonomy", "functional")
  if (!level %in% valid_levels) stop(
    paste0("'level' muste be: ", paste(valid_levels, collapse = ","))
  )
  
  
  # 2. Defalut filename
  if (is.null(filename)) {
    # Default filenames
    filename <- switch(level,
                       "peptide" = "combined_peptide.tsv",
                       "protein" = "combined_protein.tsv",
                       "taxonomy" = "taxonomy.tsv",
                       "functional" = "functions.tsv"
    )
  }
  # --- Build path --- #
  path <- here("data", state, run, level, filename)
  
  # Check if file exists
  if (!file.exists(path)) stop(paste0("Error: File does not exist: ", path))
  
  # --- Read tsv files --- #
  if (verbose) {
    df <- vroom(path, 
                delim = "\t", 
                col_names = TRUE)
  } else {
    df <- vroom(path, 
                delim = "\t", 
                col_names = TRUE,
                progress = FALSE,
                show_col_types = FALSE)
  }

  return(df)
}


# --------------------------------------------------------------------------------------------------------
#' Write metaproteomics processed files
#' 
#' Write processed data at a given biological level for a specified run.
#' This function centralizes file paths and enforces reproducibility.
#' 
#' @param df data.frame. 
#'    Biological dataset to export (eg. peptide, protein, taxonomy or functional).
#'
#' @param state Character. 
#'    It should be "processed".
#'    
#' @param run Character. 
#'    Experimental run, e.g. "run_2024", "run_2025", "run_human".
#'    
#' @param level Character. 
#'    Biological level: "peptide", "protein", "functional", or "taxonomy".
#'
#' @param subrun Character. 
#'    Experimental subrun, e.g. "run_2024_no_ID1"
#'    
#' @param filename Optional character. 
#'    If NULL, defaults to "level.tsv" (e.g., peptide.tsv).
#'    
#' @param overwrite Optional logical. 
#'    Prevent or allow file overwriting.
#'
#' @return tibble with the requested dataset.
#' 
#' @examples
#' write_metap_data(state = "raw", run = "run_2025", level = "peptide")
#' write_metap_data(state = "processed", run = "run_2024", level = "protein", filename = "custom_proteins.tsv")

write_metap_data <- function(df,
                             state = "processed", 
                             run, 
                             level,
                             subrun = NULL,
                             filename = NULL,
                             overwrite = TRUE
) {
  
  # 1. Validate inputs
  stopifnot(is.data.frame(df))
  
  valid_states <- c("processed")
  if (!state %in% valid_states){
    stop(paste0("'state' must be processed"))
  }
  
  valid_levels <- c("peptide", "protein", "taxonomy", "functional")
  if (!level %in% valid_levels) stop(
    paste0("'level' muste be: ", paste(valid_levels, collapse = ","))
  )
  
  
  # 2. Default filename 
  if (is.null(filename)) {
    # Default filenames
    filename <- switch(level,
                       "peptide" = "combined_peptide.tsv",
                       "protein" = "combined_protein.tsv",
                       "taxonomy" = "taxonomy.tsv",
                       "functional" = "functional.tsv"
    )
  }
  
  # 3. Build path
  base_path <- here("data", state, run)
  
  if (!is.null(subrun)) {
    base_path <- here(base_path, subrun)
  }
  
  dir_path <- do.call(here::here, as.list(c(base_path, level)))
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  path <- file.path(dir_path, filename)
  
  if(file.exists(path) && !overwrite) {
    stop("File already exists and overwrite = FALSE: ", path)
  }
  
  # --- Read tsv files --- #
  readr::write_tsv(df, path)
  
  message(level, " processed file saved in: ", path)
  
  invisible(path)
}





