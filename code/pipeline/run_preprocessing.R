# --------------------------------------------------------------------------------------------------------
# Script: run_preprocessing.R
# Description: Metaproteomics data preprocessing Pipeline
# Author: Andrés Arroyo Barea
# Date: 2026-01-13
#       Purpose: Generate a pipeline to orchestrate the whole metaprotemomics 
#       prepocessing process.
# --------------------------------------------------------------------------------------------------------
message("=== Starting preprocessing pipeline ===")

if (run_peptide) {
  message(">>> Peptide preprocessing starting...")
  source("code/preprocess/01_peptide_preprocess.R")
  message("[INFO]: Peptide preprocessing completed!")
}

if (run_protein) {
  message(">>> Protein preprocessing starting...")
  source("code/preprocess/02_protein_preprocess.R")
  message("[INFO]: Protein preprocessing completed!")
}

if (run_taxonomy) {
  message(">>> Taxonomy preprocessing starting...")
  source("code/preprocess/03_taxonomy_preprocess.R")
  message("[INFO]: Taxonomy preprocessing completed!")
}

if (run_functional) {
  message(">>> Functional preprocessing starting...")
  source("code/preprocess/04_functional_preprocess.R")
  message("[INFO]: Functional preprocessing completed!")
}

if (run_metadata) {
  message(">>> Metadata preprocessing starting...")
  source("code/preprocess/05_metadata_preprocess.R")
  message("[INFO]: Metadata preprocessing completed!")
}

message(">>> Building feature data list...")

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
  stop("[ERROR]: No feature table processed at any level")
}

message("[INFO]: Feature data list built.")
message(sprintf("[INFO]: Biological veles available: %s",
                paste(names(metap_feature_data), collapse = ",")))

message("[INFO]: Preprocessing finished.")
