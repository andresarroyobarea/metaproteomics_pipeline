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
message("=== Metaproteomics global config loading ===")
source("code/01_config.R")

# 3. Load EDA config
message("=== EDA config loading ===")
source("code/eda/config/eda_config.R")
#source("code/eda/core/get_eda_data.R")

# 4. Load EDA functions
message("=== EDA functions loading ===")
eda_core_files <- list.files("code/eda/core", pattern = "\\.R$", full.names = TRUE)
lapply(eda_core_files, source)

# 5. Preprocessing
source("code/pipeline/run_preprocessing.R")



metap_list <- list(
  "peptide" = peptides_processed,
  "protein" = proteins_processed,
  "taxonomy" = taxonomy_processed,
  "functional" = functional_processed
)



s <- get_eda_data(feature_data_list = metap_list,
             biological_level = "protein",
             metric_name = "intensity",
             filter_subset = "proteins_core",
             eda_mode = list(
               transform_mode = "quantitative",
               zero_strategy   = "na",
               pseudocount = 0.1,
               log_transform = "log2",
               normalization = "none",
               post_zero_strategy = "keep"
             ),
             eda_verbose = TRUE)


eda_qc_check(s)




s %>%
  dplyr::select(-feature_id) %>%
  dplyr::summarise(across(everything(), list(
    min = ~ min(., na.rm = TRUE),
    max = ~ max(., na.rm = TRUE),
    mean = ~ mean(., na.rm = TRUE),
    median = ~ median(., na.rm = TRUE),
    sd = ~ sd(., na.rm = TRUE),
    zeros = ~ sum(. == 0, na.rm = T),
    zeros_perc = ~ round(sum(. == 0, na.rm = T)/nrow(s)* 100, 2),
    NAs = ~ sum(is.na(.)),
    NAs_per = ~ round(sum(is.na(.))/nrow(s)* 100, 2)
    ))) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("sample", "stat"),
    names_pattern = "^(.*)_(min|max|mean|median|sd|zeros|zeros_perc|NAs|NAs_per)$"
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

sapply(s, function(x) sum(x == 0, na.rm = T))
colnames(proteins_processed)

help(separate)

apply_log_transform(metap_list$peptide %>% dplyr::select(matches("^ID_\\d+_intensity")), log_transform = "log10")
