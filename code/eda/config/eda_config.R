# --------------------------------------------------------------------------------------------------------
# Script: eda_config.R
# Project: Metaproteomics MM - EDA Configuration
# Author: Andrés Arroyo Barea
# Date: 2026-01-27
#
# Purpose:
#   Central configuration for Exploratory Data Analysis (EDA) in metaproteomics.
#   Define levels, metrics, filters, normalization options, and plotting preferences.
#
# Usage:
#   Source this script in your EDA scripts to load settings:
#     source("code/eda/config/eda_config.R")
# --------------------------------------------------------------------------------------------------------

# -------------------------------------------
# 1. Levels for EDA
# -------------------------------------------
# Which biological levels could be included in the EDA.
eda_levels <- c("peptide", "protein", "taxonomy", "functional")

# -------------------------------------------
# 2. Metrics per level
# -------------------------------------------
# Keep quantitative measurements available for EDA
# This is a subset of metrics from metrics_by_level, tailored for exploratory plots.
# Changing these does not affect the preprocessing or downstream analysis pipeline.
eda_metrics <- list(
  peptide = c(
    intens = "intensity",
    max_intens = "max_lfq_intensity",
    spc = "spectral_count"
  ),
  protein = c(
    intens = "intensity",
    max_intens = "max_lfq_intensity",
    spc = "spectral_count"
  ),
  functional = c(
    intens = "intensity",
    max_intens = "max_lfq_intensity"
  ),
  taxonomy = c(
    intens = "intensity"
  )
)

# -------------------------------------------
# 3. Filter sets
# -------------------------------------------
# These are predefined filters per level (from preprocessing)
# Define which subsets of features to explore in EDA.
# TODO: Incorporate specific filter not included in this sets: Eg. HUMAN proteins.
eda_filter_sets <- list(
  peptide = names(peptide_sets_defs),
  protein = names(protein_sets_defs),
  functional = names(protein_sets_defs),
  taxonomy = names(taxonomy_sets_defs)
)

# -------------------------------------------
# 4. Normalization / transformation options
# -------------------------------------------
# Define trasnformations that can be applied before plotting.
# EDA scripts can check this config to decide wheter to log-transform, scale or
# normalize the data.
# Typical proteomics normalization approaches are included in this set.
eda_normalization_options <- c(
  "none",
  "log2",
  "log10",
  "log2_median_norm",
  "log10_median_norm"
)

# Default normalization for each level
eda_default_normalization <- list(
  peptide = "log2",
  protein = "log2",
  functional = "log2",
  taxonomy = "log2"
)

# EDA zero handling
eda_zero_handling <- list(
  log_method = "na",      # "na" | "pseudocount"
  pseudocount = 1,
  after_log = "keep_na"   # "keep_na" | "back_to_zero"
)


# -------------------------------------------
# 5. Subset options
# -------------------------------------------
# EDA could be done only in a subset of samples or conditions.
# Default are all conditions included in cleaned metadata.
eda_subsets <- list(
  samples = NULL,
  condition = NULL
)

# -------------------------------------------
# 6. Biology specific plots
# -------------------------------------------
# This option include plots that are only allowed at specific 
# biological levels.
biology_specific_plots <- list(
  peptide = c("plot_peptide_length", "plot_peptide_charge"),
  protein = c("plot_protein_length", "plot_protein_unique_peptides"),
  taxonomy = c("plot_taxonomy_rank_distribution", "plot_taxonomy_relative_abundance"),
  functional = c()  
)

# -----------------------------
# 7. Plot settings (optional defaults)
# -----------------------------
# These are defaults that can be overridden in plotting functions
eda_plot_defaults <- list(
  theme = theme_minimal(base_size = 22),
  point_size = 3,
  line_size  = 1.5,
  colors    = c("NDMM" = "red","RRMM" = "blueviolet", "SHIELD" = "darkgreen"),
  font_size = 22,
  save_width  = 8,
  save_height = 6,
  save_units  = "in",
  save_dpi    = 300
)

# -----------------------------
# 8. Default EDA plots per dimension
# -----------------------------
# This can serve as a blueprint for scripts
eda_plot_plan <- list(
  feature_level   = c("abundance_distribution", "prevalence_distribution", "prevalence_vs_abundance", "filtering_impact", "biology_specific"),
  sample_level    = c("total_metric_per_sample", "feature_counts_per_sample", "missingness_heatmap", "rarefaction", "dominance_plot", "PCA", "clustering"),
  condition_level = c("metric_by_condition", "feature_counts_by_condition", "overlap_venn", "PCA_condition", "beta_diversity", "effect_size_plots")
)

# -----------------------------
# 9. Logging and messaging
# -----------------------------
eda_verbose <- TRUE   # if TRUE, scripts will print messages about subsets, metrics, normalization, filters

# --------------------------------------------------------------------------------------------------------
# End of eda_config.R
# --------------------------------------------------------------------------------------------------------

