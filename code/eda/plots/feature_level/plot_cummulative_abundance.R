# --------------------------------------------------------------------------------------------------------
# Script: plot_cummulative_abundance.R
# Description: Feature dominance.
# Author: Andrés Arroyo Barea
# Date: 2026-02-09
#       Purpose: EDA graphs.
# --------------------------------------------------------------------------------------------------------
plot_cummulative_abundance <- function(
    data_long,
    biological_level = biological_level,
    metric = metric,
    log_base = 2,
    group_by_sample = F,
    group_by_condition = F,
    aggregation = "mean",
    condition_vars = condition_vars,
    path_results_eda = path_results_eda,
    save_plots = TRUE,
    verbose = TRUE
) {
  
  # Aggregation function.
  agg_fun <- get_agg_fun(aggregation)

  plots <- list()
  
  # -------------------------------
  # 1. Aggregate abundance per feature
  # -------------------------------
  if (verbose){
    log_info(paste0("Generating cummulative overall feature rank abudance ", metric, " plot..."))
  }
  
  feature_rank_abund_cumm <- aggregate_feature_abundance(
    data_long,
    group_level = "feature",
    agg_fun = mean,
    feature_col = "feature_id",
    sample_col = "sample_id",
    condition_col = "condition",
    cumm_abundance = T
  )
  
  # -------------------------------
  # 2. Cummulative Rank-abundance plot at feature level
  # -------------------------------
  p_feature_cumm_rank <- ggplot(feature_rank_abund_cumm, aes(x = rank, y = cum_abundance)) +
    geom_line(color = "#1f77b4", size = 1) +
    labs(
      title = "Cumulative Abundance Curve",
      subtitle = paste0("Aggregation: ", aggregation, " | Features: ", nrow(feature_rank_abund_cumm)),
      x = "Feature rank (descending abundance)",
      y = "Cumulative abundance (%)"
    )
  
  # -------------------------------
  # 3. Save plot
  # -------------------------------
  plots$overall <- plot_and_save(p_feature_cumm_rank,
                                 plot_name = paste0("cumm_rank_feature_abundance_", aggregation),
                                 type_dir = file.path(path_results_eda,"feature_centric", unique(data_long$biological_level)),
                                 save_plots = save_plots
  )
  
  # ---------------------------------------------------
  # 4. Feature abundance distribution by sample
  # ---------------------------------------------------
  if (group_by_sample) {
    
    # Agregate by feature and sample
    sample_rank_abund_cumm <- aggregate_feature_abundance(
      data_long = data_long,
      group_level = "sample",
      agg_fun = agg_fun,
      feature_col = "feature_id",
      sample_col = "sample_id",
      condition_col = "condition",
      cumm_abundance = T
    )
    
    if (verbose){
      log_info(paste0("Generating cummulative rank abudance ", metric, " plot at sample level..."))
    }
    
    p_feature_sample_cumm_rank <- ggplot(sample_rank_abund_cumm, aes(x = rank, y = cum_abundance, colour = sample_id)) +
      geom_line(na.rm = TRUE, alpha = 0.8) +
      labs(
        title = "Rank–Cumm-Abundance Curve - Sample level",
        subtitle = paste0(
          "Aggregation: ", aggregation,
          " | Features (", biological_level,  "): ", nrow(feature_rank_abund_cumm)
        ),
        x = "Feature rank (descending abundance)",
        y = paste0(log_base, "(", metric, ")")
      ) +
      scale_color_paletteer_d("rcartocolor::Antique")
    
    plots$per_sample <- plot_and_save(p_feature_sample_cumm_rank,
                                      plot_name = paste0(metric, "cumm_rank_feature_abundance_sample_level_",aggregation),
                                      type_dir = file.path(path_results_eda,"feature_centric",unique(data_long$biological_level)),
                                      save_plots = save_plots
    )
  }
  
  # ---------------------------------------------------
  # 4. Feature abundance distribution by condition
  # ---------------------------------------------------
  if (group_by_condition) {
    for (cond in condition_vars) {
      
      if (verbose) {
        log_info(paste0("Generating cummulative overall rank-abudance ", metric, " plot for ", cond, " level..."))
      }
      # Agregate by feature and sample
      condition_rank_abund_cumm <- aggregate_feature_abundance(
        data_long = data_long,
        group_level = "condition",
        agg_fun = mean,
        feature_col = "feature_id",
        sample_col = "sample_id",
        condition_col = cond,
        cumm_abundance = T
      )
      
      p_feature_condition_rank_cumm <- ggplot(condition_rank_abund_cumm, aes(x = rank, y = cum_abundance, colour = .data[[cond]])) +
        geom_line(na.rm = TRUE, alpha = 0.8) +
        #facet_grid(rows = vars(.data[[cond]])) +
        labs(
          title = paste0("Rank–Cumm-Abundance Curve - ", cond),
          subtitle = paste0(
            "Aggregation: ", aggregation,
            " | Features (", biological_level,  "): ", nrow(feature_rank_abund_cumm)
          ),
          x = "Feature rank (descending abundance)",
          y = paste0(log_base, "(", metric, ")")
        ) +
        scale_color_paletteer_d("rcartocolor::Antique")
      
      
      plots[[paste0("per_condition_", cond)]] <- plot_and_save(p_feature_condition_rank_cumm,
                                                               plot_name = paste0(metric,
                                                                                  "cumm_rank_feature_abundance_condition_level_",
                                                                                  aggregation),
                                                               type_dir = file.path(path_results_eda,"feature_centric",unique(data_long$biological_level)),
                                                               save_plots = save_plots
      ) 
    }
  }
  return(plots)
}