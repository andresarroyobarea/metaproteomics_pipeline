# --------------------------------------------------------------------------------------------------------
# Script: run_preprocessing.R
# Description: Metaproteomics data preprocessing Pipeline
# Author: Andrés Arroyo Barea
# Date: 2026-01-13
#       Purpose: Generate a pipeline to orchestrate the whole metaprotemomics 
#       prepocessing process.
# --------------------------------------------------------------------------------------------------------
source("code/00_setup_environment.R")
source("code/01_config.R")

message("=== Starting preprocessing pipeline ===")

if (run_peptide) {
  message(">>> Peptide preprocessing starting...")
  source("code/preprocess/01_peptide_preprocess.R")
  message("[OK]: Peptide preprocessing completed!")
}

if (run_protein) {
  message(">>> Protein preprocessing starting...")
  source("code/preprocess/02_protein_preprocess.R")
  message("[OK]: Protein preprocessing completed!")
}

if (run_taxonomy) {
  message(">>> Taxonomy preprocessing starting...")
  source("code/preprocess/03_taxonomy_preprocess.R")
  message("[OK]: Taxonomy preprocessing completed!")
}

if (run_functional) {
  message(">>> Functional preprocessing starting...")
  source("code/preprocess/04_functional_preprocess.R")
  message("[OK]: Functional preprocessing completed!")
}