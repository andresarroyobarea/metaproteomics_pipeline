# --------------------------------------------------------------------------------------------------------
# Script: run_preprocessing.R
# Description: Metaproteomics data preprocessing Pipeline
# Author: Andrés Arroyo Barea
# Date: 2026-01-13
#       Purpose: Generate a pipeline to orchestrate the whole metaprotemomics 
#       prepocessing process.
# --------------------------------------------------------------------------------------------------------
log_info("=== Starting preprocessing pipeline ===")

if (run_peptide) {
  log_info(">>> Peptide preprocessing starting...")
  source("code/preprocess/01_peptide_preprocess.R")
  log_info("Peptide preprocessing completed!")
}

if (run_protein) {
  log_info(">>> Protein preprocessing starting...")
  source("code/preprocess/02_protein_preprocess.R")
  log_info("Protein preprocessing completed!")
}

if (run_taxonomy) {
  log_info(">>> Taxonomy preprocessing starting...")
  source("code/preprocess/03_taxonomy_preprocess.R")
  log_info("Taxonomy preprocessing completed!")
}

if (run_functional) {
  log_info(">>> Functional preprocessing starting...")
  source("code/preprocess/04_functional_preprocess.R")
  log_info("Functional preprocessing completed!")
}

if (run_metadata) {
  log_info(">>> Metadata preprocessing starting...")
  source("code/preprocess/05_metadata_preprocess.R")
  log_info("Metadata preprocessing completed!")
}

log_info(">>> Building feature data list...")

metap_feature_data <- list()

if (run_peptide) {
  metap_feature_data$peptide <- peptides_processed
}

if (run_protein) {
  metap_feature_data$protein <- proteins_processed
}

if (run_taxonomy) {
  metap_feature_data$taxonomy <- taxonomy_processed
}

if (run_functional) {
  metap_feature_data$functional <- functional_processed
}

if (length(metap_feature_data) == 0) {
  log_error("No feature table processed at any level")
}

log_info("Feature data list built.")
log_info(sprintf("Biological veles available: %s",
                paste(names(metap_feature_data), collapse = ",")))

log_info("Preprocessing finished.")
