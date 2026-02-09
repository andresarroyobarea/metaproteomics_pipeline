# --------------------------------------------------------------------------------------------------------
# Script: run_metaproteomics_pipeline.R
# Description: Complete Metaproteomics Data Analysis Pipeline
# Author: Andrés Arroyo Barea
# Date: 2026-01-13
#       Purpose: Generate a pipeline to orchestrate the whole metaprotemomics 
#       data analysis
# --------------------------------------------------------------------------------------------------------
# 1. Setup environment (packages, paths)
message("=== Metaproteomics environment setup ===")
source("code/00_setup_environment.R")

# 2. Load global config (project-wide settings)
log_info("=== Metaproteomics global config loading ===")
source("code/01_config.R")

# 3. Load EDA config
log_info("=== EDA config loading ===")
source("code/eda/config/eda_config.R")
#source("code/eda/core/get_eda_data.R")

# 4. Load EDA functions
log_info("=== EDA functions loading ===")
eda_core_files <- list.files("code/eda/core", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(eda_core_files, source))

# 5. Preprocessing
source("code/pipeline/run_preprocessing.R")