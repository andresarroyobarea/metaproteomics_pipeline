# --------------------------------------------------------------------------------------------------------
# Script: run_eda_feature_centric.R
# Description: EDA - Feature-centric approach.
# Author: Andrés Arroyo Barea
# Date: 2026-02-10
#       Purpose: Perform exploratory data analysis of metaproteomics data at 
#       feature level.
# --------------------------------------------------------------------------------------------------------

# GRID FILTRADO PARA PRUEBA
# HAY QUE DEFINIR COMO INTEGRAR ESTE GRID FILTRADO, DONDE UBICARLO Y SI SE PUEDE
# INCLUIR SIEMPRE DESDE EL PRINCIPIO O SE DEBE FILTRAR INTERNAMENTE SEGUN EL TIPO DE GRÁFICO
eda_grid_prueba <- create_eda_grid(
  feature_data_list = metap_feature_data,
  eda_metrics_user = eda_metrics,
  eda_transform_modes = eda_transform_modes,
  eda_levels = c("features"),
  eda_feature_subset = eda_feature_subset,
  eda_modes_filter = "log2") 



run_eda_feature_centric <- function(
    feature_data_list,
    sample_metadata,
    feature_metadata = NULL,
    eda_grid,
    enable_facet_sample = TRUE,
    enable_facet_condition = TRUE,
    condition_vars = condition_vars, 
    save_plots = TRUE,
    path_results_eda,
    verbose = eda_verbose) {
  
  log_info(">>>> Running Feature-Centric EDA")
  
  # 0. # Initialize list to save results.
  eda_feature_centric_res <- list()
  
  # 1. I think first we should create all transformations only one time and save
  # in cache to avoid repetitions.
  #eda_cache <- list()
  
  
  
  # 2. pwalk
  purrr::pwalk(
    eda_grid,
    function(
    eda_level,
    biological_level,
    metric,
    subset,
    eda_transform_mode_name,
    eda_transform_modes
    ) {
      
      if (verbose) {
        log_info(glue::glue(
          "Running feature-centric EDA: {biological_level} - {metric} - {subset} - {eda_transform_mode_name}",
          .envir = environment()
        ))
        
      }
      
      # 1. Data extraction with get_eda_data
      data_wide <- get_eda_data(
        feature_data_list = metap_feature_data,
        biological_level = biological_level,
        metric_name = metric,
        filter_subset = subset,
        eda_mode = eda_transform_modes,
        eda_verbose = verbose
      )
      
      # 2. Transform wide data to long + metadata.
      data_long <- eda_prepare_data(
        data_wide = data_wide,
        sample_metadata = metadata_filt,
        verbose = verbose,
        metric = metric  
      )
      
      # 3. Generate plots
      
      plots <- list()
      
      # Plot abundance distribution
      plots$abundance_distribution <- plot_abundance_distribution(
        data_long = data_long,
        biological_level = biological_level,
        metric = metric,
        log_base = eda_transform_mode_name,
        group_by_sample = enable_facet_sample,
        group_by_condition = enable_facet_condition,
        condition_var = condition_vars, 
        plot_feature_means = T,
        path_results_eda = path_results_eda,
        save_plots = save_plots,
        verbose = verbose
      )
      
      # ADD PLOT FEATURE DOMINANCE.
      
      
      eda_feature_centric_res[[biological_level]][[metric]][[subset]][[eda_transform_mode_name]] <<- plots 
    }
    )
  return(eda_feature_centric_res)
}


# ---- METRICAS CLAVES QUE TENEMOS QUE CONTROLAR -----
# Niveles biológicos que el usuario haya querido añadir.
# Metricas definidas para cada nivel biológico
# Subset de features donde hacerlo
# Tipo de procesamiento para cada momento.

# TODO: Think about tables only one time and then iterate over them. 
# TODO: Test different conditions.

# PRUEBA:
res <- run_eda_feature_centric(
    metap_feature_data,
    metadata_filt,
    feature_metadata = NULL,
    eda_grid_prueba,
    enable_facet_sample = T,
    enable_facet_condition = T,
    condition_vars = condition_vars, # TODO: ADD groupping by different conditions variables. It should be defined in the config file.
    save_plots = T,
    path_results_eda = path_results_eda,
    verbose = eda_verbose)

