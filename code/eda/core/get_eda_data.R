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
  eda_mode = list(
    biological_mode = "raw",        # raw | abundance
    zero_strategy   = "keep",       # keep | na | pseudo
    transform       = "none",       # none | log2 | log10
    normalization   = "none",       # none | median
    post_zero       = "keep"        # keep | zero
  ),
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
  df <- apply_data_transform(df, metric_col = metric_name, method = normalization)
  
  # Filter sample in case it was selected in eda config.
  # al mismo tiempo que se extraen las muestras?
  
  return(df)
}



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
apply_data_transform <- function(
    df,
    biological_level,
    transform_mode = "process",
    verbose = TRUE
) {
  
  if(is.null(df)) {
    stop("[ERROR]: Biological level not found", biological_level)
  }
  
  # 2. If no filter requested, return data as-is.
  if (is.null(filter_subset) | filter_subset == "none") {
    if (verbose) message("[INFO] No feature normalization applied (", biological_level, ").")
    return(df)
  }
  
  # 3. Not further processing
  if (transform_mode == "raw") {
    if (verbose) message("[INFO]: No transformation applied.")
    return(df)
  }
  
  # 4. Absence/presence
  if (transform_mode == "abundance") {
    df_abs_pres <- df %>% mutate_all(~ ifelse(. > 0, 1, 0))
    if (verbose) message("[INFO]: Data was converted into absence/presence.")
    return(df_abs_pres)
  }
  
  # 5. Processing : Log2 or Log2 normalization.
  if (transform_mode == "process") {
    
    
    
    return(df)
  }
  
  
  

}



# --------------------------------------------------------------------------------------------------------
#' apply_zero_handling
#'
#' Description
#' 
#' @param df 
#'    
#' @param strategy 
#'    
#' @param pseudocount
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------

apply_zero_handling <- function(df, strategy = "none", pseudocount = 1) {
  
  # Not change
  if (strategy == "none") {
    message("[INFO]: Zeros were ketp as zeros.")
    return(df)
  }
  
  # Change 0 to NA
  if (strategy == "na") {
    df <- df %>% mutate(across(where(is.numeric), ~ ifelse(. == 0, NA, .)))
    message("[INFO]: Zeros were turned into NA previous log transformation.")
    return(df)
  }
  
  # Add pseudocount if needed
  if (strategy == "pseudocount") {
  # TODO: Add different zero count strategy for intensities and spectral counts.
    message(paste0("[INFO]: A pseudocount of", pseudocount,  "was added previous log transformation."))
    df <- df %>% mutate(across(where(is.numeric), ~ . + pseudocount))
    return(df)
  }
  
  stop("[ERROR]: Unknown strategy: ", strategy)
}

# --------------------------------------------------------------------------------------------------------
#' apply_log_transform
#'
#' Description
#' 
#' @param df 
#'    
#' @param transform 
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------


apply_log_transform <- function(df, transform = "none") {
  
  # No transformation
  if (transform == "none") {
    message("[INFO]: NO log transformation was applied.")
    return(df)
  }
  
  # Log transformations
  if ((transform == "log2" | transform == "log10") & any(df == 0,na.rm = T)) {
    message("[WARNING] Log transformation will be applied over 0: Inf values were returned.")
  }
  
  # Log2 transformation  
  if (transform == "log2") {
    df <- df %>% mutate(across(where(is.numeric), ~ log2(.)))
    message("[INFO]: Log2 transformation was applied.")
    return(df)
    }
    
  # Log10 transformation
  if (transform == "log10") {
    df <- df %>% dplyr::mutate(across(where(is.numeric), ~ log10(.)))
    message("[INFO]: Log10 transformation was applied.")
    return(df)
  }
  
  stop("Unknown transform: ", transform)
}

# --------------------------------------------------------------------------------------------------------
#' apply_normalization
#'
#' Description
#' 
#' @param df 
#'    
#' @param norm_mode 
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
apply_normalization <- function(df, norm_mode = "median") {
  
  # TODO: Add quantile normalization with preprocessCore R package.
  
  # No Normalization
  if (norm_mode == "none"){
    message("[INFO]: NO normalization was applied.")
    return(df)
  }
  
  if (any(df == -Inf | df == Inf, na.rm = T)) {
    stop("[ERROR]: numeric columns contain Inf or -Inf. This usually happens if you log-transformed zeros. Replace zeros with a small value or NA before normalizing.")
  }
  
  # Median normalization.
  if (norm_mode == "median"){
    df_norm <- df %>% dplyr::mutate(across(where(is.numeric), ~ . - median(., na.rm = TRUE)))
    message("[INFO]: Median normalization was applied.")
    return(df_norm)
  }
  
  stop("Unknown normalization mode: ", norm_mode)
  
}

# --------------------------------------------------------------------------------------------------------
#' apply_post_zero_handling
#'
#' Description
#' 
#' @param df 
#'    
#' @param strategy 
#' 
#' @return data.frame
#'    
#' @examples
#' 
#'
#' @export
#'
#'--------------------------------------------------------------------------------------------------------
apply_post_zero_handling <- function(df, strategy = "none"){
  
  # Sometimes zeros or NA were needed for different purposes after normalization.
  
  # Zeros retained
  if (strategy == "none") {
    message("[INFO]: NO post-zero handling after normalization.")
    return(df) 
  }
  
  # NA -> ZERO
  if (strategy == "zero" & any(is.na(df))) {
    df <- df %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))
    message("[INFO]: NAs were turned into zeros")
    return(df)
  }
  
  stop("Unknown post-zero handling mode: ", norm_mode)
  
}






