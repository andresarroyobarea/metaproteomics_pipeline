# --------------------------------------------------------------------------------------------------------
# Script: eda_prepare_long_dataset.R
# Description: Prepare proteomics data in long format for EDA.
# Author: Andrés Arroyo Barea
# Date: 2026-02-09
#       Purpose: Extract EDA data in long format for EDA.
# --------------------------------------------------------------------------------------------------------
eda_prepare_data <- function(
    data_wide,
    sample_metadata = NULL,
    feature_metadata = NULL,
    verbose = TRUE,
    metric = "spectral_count" # TODO: Quitar parámetro cuando se integre con funciones superiores. 
) {
  
  # Turn data into long format and merge with metadata
  if (verbose) 
    log_info("Pivotting data to long format and metadata aggregation.")
    
    data_long <-  data_wide %>%
      pivot_longer(
        cols = -feature_id,    
        names_to = "sample_id",
        values_to = "value"
      ) %>%
      # TODO: Think about if we should remove pattern here or only for graphs.
      dplyr::mutate(
        sample_id = gsub(paste0("_", metric, "$"), "", sample_id)
      ) %>% 
      dplyr::filter(!is.na(value))
    
    if (!is.null(sample_metadata)) {
      data_long <- dplyr::left_join(data_long, sample_metadata, by = "sample_id")
    }
    
    if (!is.null(feature_metadata)) {
      data_long <- dplyr::left_join(data_long, feature_metadata, by = "feature_id")
    }
  
  # Data long validation.
  required_cols <- c("feature_id","sample_id","value")
  missing <- setdiff(required_cols, colnames(data_long))
  
  if (length(missing) > 0) 
    log_stop("Some relevant variables are missing to plot abundance metrics.")
  
  return(data_long)
}






