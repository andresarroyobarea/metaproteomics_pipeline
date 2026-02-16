# --------------------------------------------------------------------------------------------------------
# Script: plot_abundance_distribution.R
# Description: Ditribution of global features after log transformation.
# Author: Andrés Arroyo Barea
# Date: 2026-02-09
#       Purpose: EDA graphs.
# --------------------------------------------------------------------------------------------------------

plot_abundance_distribution <- function(
    data_long,
    biological_level = biological_level,
    metric = metric,
    log_base = 2,
    group_by_sample = F,
    group_by_condition = F,
    condition_vars = c("condition"), # Importante porque podríamos tener varias comparaciones.
    plot_feature_means = F,
    path_results_eda = path_results_eda,
    save_plots = TRUE,
    verbose = TRUE
) {
  
  plots <- list()
  # -------------------------------
  # 4.1 Feature abundance global
  # -------------------------------
  if (verbose) log_info(paste0("4.1. Generating overall abudance, ",  metric, "plot..."))
  p_feature <- ggplot(data_long, aes(x = value)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    labs(
      title = "Feature Abundance Distribution - Overall",
      x = paste0(log_base, " (", metric, ")"),
      y = "Number of features"
    )
  
  plots$overall <- plot_and_save(
    p_feature,
    paste0(metric, "_distribution_overall"),
    file.path(path_results_eda, "feature_centric", biological_level),
    save_plots
  )
  
  # 4.2 --- Feature abundance distribution by sample --- 
  if (group_by_sample) {
    log_info(paste0("4.2. Generating overall abudance ", metric, " plot at sample level..."))
    p_sample <- ggplot(data_long, aes(x = value, color = sample_id)) +
      geom_density(na.rm = TRUE, alpha = 0.3, width = 1) +
      labs(
        title = "Feature Abundance Distribution - Sample level",
        x = paste0(log_base, " (", metric, ")"),
        y = "Density"
      ) +
      scale_color_paletteer_d("rcartocolor::Antique")
    
    plots$per_sample <- plot_and_save(
      p_sample,
      paste0(metric, "_distribution_per_sample"),
      file.path(path_results_eda, "sample_centric", biological_level),
      save_plots
    )
  }
  
  # 4.3 --- Feature abundance distribution by condition ---
  if (group_by_condition) {
    for (cond in condition_vars) {
      
      log_info(paste0("4.3. Generating overall abudance ", metric, " plot for ", cond, " level..."))
      
      p_condition <- ggplot(
        data_long, 
        aes(x = value, color = .data[[cond]])
        ) +
        geom_density(na.rm = TRUE, alpha = 0.3, width = 1) +  
        labs(
          title = paste0("Feature Abundance Distribution - ", cond),
          x = paste0(log_base, " (", metric, ")"),
          y = "Density"
        ) +
        scale_color_paletteer_d("rcartocolor::Antique")
      
      plots[[paste0("per_condition_", cond)]] <- plot_and_save(
        p_condition,
        paste0(metric, "_distribution_per_", cond),
        file.path(path_results_eda, "condition_centric", biological_level),
        save_plots
      ) 
    }
  
  }
  
  # 4.4 Feature means plot.
  if (plot_feature_means) {
    log_info("4.4. Generating feature-mean abudance plot...")
    feature_summary <- data_long %>% 
      group_by(feature_id) %>%
      dplyr::summarise(mean_abundance = mean(value, na.rm = T)) 
    
    p_feature_means <- ggplot(feature_summary, aes(x = mean_abundance)) +
      geom_histogram(bins = 50, fill = "skyblue", color = "black") +
      labs(
        title = "Feature Abundance Distribution",
        x = paste0("Feature-mean (",  log_base, "-", metric, ")"),
        y = "Number of features"
      ) 
    
    plots$feature_means <- plot_and_save(
      p_feature_means,
      paste0("metric", "_feature_mean_distribution"),
      file.path(path_results_eda, "feature_centric", biological_level),
      save_plots
    )
  }
  return(plots)
}