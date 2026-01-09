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
#' @param df data.frame. Biological dataset (eg. peptide, protein, taxonomy or functional)
#' @param samples_to_include Character. Samples selected by user. They come from 01_config.R.
#' @param sample_pattern Character. Initial identifier to identify sample information in Fragpipe/Imetalab files.
#' @param verbose Optional logical. Display messages.
#' @return tibble with the requested dataset.
#' @examples
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
#' @param protein_var_name character. Protein name variable retrieved from biological datasets (peptide, protein, taxonomy, functional)
#' @param human_pattern character. Specific pattern to identify human proteins. It could change between different reference databases.
#' @return YES/NO factor
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
#' @return Factor with levels "YES", "NO".
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

