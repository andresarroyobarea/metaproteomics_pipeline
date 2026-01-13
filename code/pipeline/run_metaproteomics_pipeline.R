# --------------------------------------------------------------------------------------------------------
# Script: run_metaproteomics_pipeline.R
# Description: Complete Metaproteomics Data Analysis Pipeline
# Author: Andrés Arroyo Barea
# Date: 2026-01-13
#       Purpose: Generate a pipeline to orchestrate the whole metaprotemomics 
#       data analysis
# --------------------------------------------------------------------------------------------------------
source("code/00_setup_environment.R")
source("code/01_config.R")

# --- Preprocessing --- 
source("code/pipeline/run_preprocessing.R")

