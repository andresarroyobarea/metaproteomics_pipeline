# --------------------------------------------------------------------------------------------------------
#' get_eda_data
#'
#' Description
#' 
#' @param feature_data_list 
#'    
#' @param biological_level 
#'    
#' @param metric_name
#'
#' @param normalization 
#' 
#' @param verbose
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
get_eda_data <- function(
  feature_data_list,
  biological_level,
  metric_name,
  filter_subset,
  normalization = "none",
  eda_verbose = TRUE
) {
  
  # Select dataframe (preguntar esto)
  df <- feature_data_list[[biological_level]]
  stopifnot(!is.null(df))
  
  if (eda_verbose) message("[INFO] Starting EDA: Level = ", biological_level, "; Metric = ", metric_name)
  
  # 2. Apply feature filter subset if provided
  if (!is.null(filter_subset)) {
    # TODO: Translate filter subset into filtering criteria
    df <- apply_filter_subset(df, filter_subset, biological_level)
    if (eda_verbose) message("EDA: Applied filter susbet: ", filter_subset)
  }
  
  # 3. Select metric of interest
  df <- df %>% dplyr::select(feature_id, matches(paste0("ID_\\d+_", metric_name, "$")))
 
  # 4. Apply normalization.
  df <- normalize_metric(df, metric_col = metric_name, method = normalization)
  
  # Filter sample in case it was selected in eda config.
  # al mismo tiempo que se extraen las muestras?
  
  return(df)
}

log10(0)
eda_metrics

# --------------------------------------------------------------------------------------------------------
#' apply_filter_subset
#'
#' Description
#' 
#' @param feature_data_list 
#'    
#' @param filter_subset 
#'    
#' @param biological_level
#' 
#' @param verbose
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
apply_feature_filtering <- function(
    feature_data_list,
    biological_level,
    filter_subset,
    verbose = TRUE
) {
  
  # 1. Select dataframe by biological level.
  df <- feature_data_list[[biological_level]]
  
  if(is.null(df)) {
    stop("[ERROR]: Biological level not found", biological_level)
  }
  
  # 2. If no filter requested, return data as-is.
  if (is.null(filter_subset)) {
    if (verbose) message("[INFO] No feature filtering applied (", biological_level, ").")
    return(df)
    }
  
  # 3. Check if filter_subset exits as a colum
  if (!filter_subset %in% colnames(df)) {
    stop("[ERROR] Filter subset not found in data: ", filter_subset)
    }
  
  if (verbose) {
    message("[INFO] Applying filter subset: ", biological_level, " -> ", filter_subset)
  }
  
  # 4. Apply filter
  df_filtered <- df %>% dplyr::filter(.data[[filter_subset]])
  
  # Filter sample in case it was selected in eda config.
  # al mismo tiempo que se extraen las muestras?
  if (verbose) {
    message("[INFO] Filtering completed: ", 
            nrow(df_filtered), " / ", nrow(df), " features retained.")
  } 
  
  return(df_filtered)
}



# --------------------------------------------------------------------------------------------------------
#' normalize_prot_metric
#'
#' Description
#' 
#' @param feature_data_list 
#'    
#' @param filter_subset 
#'    
#' @param biological_level
#' 
#' @param verbose
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
normalize_prot_metric <- function(
    feature_data_list,
    biological_level,
    normalization,
    verbose = TRUE
) {
  # 1. Select dataframe by biological level.
  df <- feature_data_list[[biological_level]]
  
  if(is.null(df)) {
    stop("[ERROR]: Biological level not found", biological_level)
  }
  
  # 2. If no filter requested, return data as-is.
  if (is.null(filter_subset) | filter_subset == "none") {
    if (verbose) message("[INFO] No feature normalization applied (", biological_level, ").")
    return(df)
  }
  

  
}



