# --------------------------------------------------------------------------------------------------------
#' get_eda_data
#'
#' Prepare quantitative or detection data for exploratory data analysis (EDA)
#' and downstream statistical analyses in metaproteomics.
#' 
#' This function acts as a unified data access and transformation layer.
#' Starting from preprocessed feature-level data (peptides, proteins, taxonomy
#' or functional annotations), it applies:
#'   - predefined feature filtering subsets,
#'   - metric selection (e.g. intensity, spectral counts),
#'   - optional transformation into detection data (presence/absence),
#'   - controlled zero handling strategies,
#'   - log-transformation and normalization when required.
#'  
#' 
#' @param feature_data_list list.
#'    Named list containing preprocessed biological datasets.
#'    Expected elements include: \code{peptide}, \code{protein},
#'   \code{taxonomy} and \code{functional}.
#'    
#' @param biological_level character.
#'    Biological level to extract data from.
#'    One of: \code{"peptide"}, \code{"protein"}, \code{"taxonomy"},
#'   \code{"functional"}.
#'   
#' @param metric_name character.
#'   Proteomics metric key to retrieve (e.g. \code{"intensity"}, \code{"spectral_count"}).
#'   The metric must be defined for the selected biological level in
#'   \code{eda_metrics}.
#'    
#'TODO: Check if this parameter should be the key or the values of metrics
#'
#' @param filter_subset character.
#'    Name of a logical column generated during preprocessing used to subset
#'    features (e.g. \code{"proteins_core"}, \code{"peptides_core"}).
#'    
#' @param eda_mode list.
#'   List defining how the data should be transformed.
#'   Supported options:
#'    \itemize{
#'     \item \code{transform_mode}: \code{"raw"}, \code{"detection"} (presence/absence)
#'           or \code{"quantitative"}.
#'     \item \code{zero_strategy}: \code{"keep"}, \code{"na"} or \code{"pseudocount"}.
#'     \item \code{pseudocount}: Numeric pseudocount added when
#'           \code{zero_strategy = "pseudocount"}. Default: 0.1.
#'     \item \code{log_transform}: \code{"none"}, \code{"log2"} or \code{"log10"}.
#'     \item \code{normalization}: \code{"none"} or \code{"median"}.
#'     \item \code{post_zero_strategy}: \code{"keep"} or \code{"zero"}, applied
#'           after normalization.
#'   }
#'      
#' @param eda_verbose logical, default TRUE
#'   Whether to print informative messages describing applied filters,
#'   transformations and potential issues.
#' 
#' @return data.frame
#'   A feature-by-sample data frame containing \code{feature_id} and the selected
#'   metric columns, transformed according to \code{eda_mode}.
#'    
#' @examples
#' # Example 1: Protein-level intensity data for core proteins,
#' # log2-transformed and median normalized (typical PCA-ready input)
#'
#' protein_eda <- get_eda_data(
#'   feature_data_list = feature_data_list,
#'   biological_level  = "protein",
#'   metric_name       = "intensity",
#'   filter_subset     = "proteins_core",
#'   eda_mode = list(
#'     transform_mode      = "quantitative",
#'     zero_strategy       = "na",
#'     log_transform       = "log2",
#'     normalization       = "median",
#'     post_zero_strategy  = "keep"
#'   )
#' )
#'
#' # Example 2: Peptide detection matrix (presence/absence)
#' peptide_detection <- get_eda_data(
#'   feature_data_list = feature_data_list,
#'   biological_level  = "peptide",
#'   metric_name       = "spectral_count",
#'   filter_subset     = "peptides_core",
#'   eda_mode = list(
#'     transform_mode = "detection"
#'   )
#' )
#' @export
#'--------------------------------------------------------------------------------------------------------
get_eda_data <- function(
  feature_data_list,
  biological_level,
  metric_name,
  filter_subset,
  eda_mode = list(
    transform_mode = "raw",
    zero_strategy   = "keep",
    pseudocount = 0.1,
    log_transform = "none",
    normalization = "none",
    post_zero_strategy = "keep"
  ),
  eda_verbose = TRUE
) {
  
  # 1. Select biological level dataframe
  df <- feature_data_list[[biological_level]]
  if(is.null(df)) {
    stop("[ERROR]: Biological level not found '", biological_level, "'")
  }
  
  if (eda_verbose) message("[INFO]: Starting EDA data preparation: Level = ", biological_level, "; Metric = ", metric_name)
  
  # 2. Apply feature filter subset if provided
  df <- apply_feature_filtering(df, filter_subset = filter_subset, verbose = eda_verbose)
  
  # 3. Apply metric selection
  df <- apply_metric_selection(df = df, biological_level = biological_level, metric_key = metric_name, 
                               eda_metrics = eda_metrics, verbose = eda_verbose)
  
  
  # 4. Apply data transform
  df <- apply_data_transform(df = df, 
                             #biological_level = biological_level, 
                             transform_mode = eda_mode$transform_mode,
                            zero_strategy = eda_mode$zero_strategy, pseudocount = eda_mode$pseudocount, 
                            log_transform = eda_mode$log_transform, normalization = eda_mode$normalization,
                            post_zero_strategy = eda_mode$post_zero_strategy, verbose = eda_verbose)

  # 5. Final checks
  if (nrow(df) == 0) {
    stop("[ERROR] No features were retrieved.")
  }
  
  if (ncol(df) == 0) {
    stop("[ERROR] Samples/metrics were not properly retrieved")
  }
  
  return(df)
}

# --------------------------------------------------------------------------------------------------------
#' apply_filter_subset
#'
#' Subset a biological feature table using logical filter variables generated
#' during the preprocessing step (e.g. \code{proteins_core}, \code{peptides_core},
#' \code{taxa_core}).
#' 
#' This function is a lightweight wrapper around feature-level filtering logic.
#' It assumes that filtering criteria have already been defined and evaluated
#' upstream during preprocessing, and simply applies the selected logical column
#' to subset the dataset.
#' 
#' @param df data.frame.
#'   A biological-level feature table (peptide, protein, taxonomy or functional)
#'   containing logical filter columns.
#'    
#' @param filter_subset 
#'   Name of the logical column used to subset the data. This column must exist
#'   in \code{df} and evaluate to \code{TRUE/FALSE} at the feature level.
#' 
#' @param verbose
#'   Whether to print informative messages describing the applied filtering and
#'   the number of retained features.
#'    
#' @return data.frame.
#'   A filtered version of \code{df} containing only features that satisfy the
#'   selected filter subset.
#'
#' @examples
#' # Apply core protein filtering
#' protein_filtered <- apply_feature_filtering(
#'   df = proteins_processed,
#'   filter_subset = "proteins_core",
#'   verbose = TRUE
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
apply_feature_filtering <- function(
    df,
    filter_subset,
    verbose = TRUE
) {
  
  # 1. If no filter requested, return data as-is.
  if (is.null(filter_subset)) {
    if (verbose) message("[INFO]: No feature filtering applied")
    return(df)
    }
  
  # 2. Check if filter_subset exits as a colum
  if (!filter_subset %in% colnames(df)) {
    stop("[ERROR]: Filter subset not found in data: ", filter_subset)
    }
  
  if (verbose) {
    message("[INFO]: Applying filter subset: ", filter_subset, " selected")
  }
  
  # 4. Apply filter
  df_filtered <- df %>% dplyr::filter(.data[[filter_subset]])
  
  # Filter sample in case it was selected in eda config.
  # al mismo tiempo que se extraen las muestras?
  if (verbose) {
    message("[INFO]: Filtering completed: ", 
            nrow(df_filtered), " / ", nrow(df), " features retained.")
  } 
  
  return(df_filtered)
}


# --------------------------------------------------------------------------------------------------------
#' apply_metric_selection
#'
#' Select a single proteomics metric from a biological-level feature table.
#' 
#' This function extracts sample-level quantitative measurements corresponding
#' to a specific proteomics metric (e.g. intensity, spectral count).
#' 
#' Metric selection is performed by matching column name patterns that encode
#' sample identifiers and metric names (e.g. \code{ID_123_intensity}).
#' 
#' @param df data.frame
#'   A biological-level feature table (peptide, protein, taxonomy or functional)
#'   containing sample-level proteomics metrics.
#'    
#' @param biological_level character.
#'    Biological level from which the metric is extracted.
#'    One of: \code{"peptide"}, \code{"protein"}, \code{"taxonomy"},
#'   \code{"functional"}.
#'   
#' @param metric_key character.
#'   Proteomics metric key to retrieve (e.g. \code{"intensity"},
#'   \code{"spectral_count"}). The metric must be defined for the selected
#'   biological level in \code{eda_metrics}.
#' 
#' @param eda_metrics list.
#'   Named list defining the allowed proteomics metrics for each biological level.
#' 
#' @param verbose logical, default \code{TRUE}
#'   Whether to print informative messages describing the metric selection process.
#'  
#' @return data.frame
#'   A reduced version of \code{df} containing \code{feature_id} and the selected
#'   sample-level metric columns.
#'   
#' @examples
#' # Extract protein-level intensity values
#' protein_intensity <- apply_metric_selection(
#'   df = proteins_processed,
#'   biological_level = "protein",
#'   metric_key = "intensity",
#'   eda_metrics = eda_metrics,
#'   verbose = TRUE
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
apply_metric_selection <- function(
    df,
    biological_level,
    metric_key,
    eda_metrics,
    verbose = TRUE
) {
  
  # 1. Checking
  
  # 1.1 A proteomic metric must be supplied.
  if (is.null(metric_key)) {
    if (verbose) stop("[ERROR]: Any proteomic metric must be supplied.")
  }
  
  # 1.2. Check available metrics at this biological level.
  available_metrics <- eda_metrics[[biological_level]]
  
  if (!metric_key %in% eda_metrics[[biological_level]]) {
    stop(
      "[ERROR]: Metric '", metric_key,
      "' not available for level '", biological_level,
      "'. Available: ", paste(available_metrics, collapse = ", "))
  }
  
  if (verbose) {
    message(
      "[INFO]: Selecting metrics: ", metric_key
    )
  }
  
  # 2. Apply metric selection
  metric_pattern <- paste0("ID_\\d+_", metric_key, "$")
  
  df_metric <- df %>% 
    dplyr::select(feature_id, matches(metric_pattern)) 
  
  # 3. Final checking
  if (ncol(df_metric) == 1) {
    stop("[ERROR]: No columns found for metric: ", metric_key)
  }
  
  if (verbose) {
    message("[INFO]: Metric selection completed (", 
            ncol(df_metric) - 1, " samples)")
  } 
  
  return(df_metric)
}


# --------------------------------------------------------------------------------------------------------
#' apply_data_transform
#'
#' Apply transformation and normalization strategies to prepare proteomics data
#' for exploratory data analysis (EDA) and downstream statistical analyses.
#'
#' This function implements a modular transformation pipeline including
#' presence/absence conversion, zero handling, log-transformation and
#' normalization.
#' 
#' @param df data.frame
#'   A biological-level feature table (peptide, protein, taxonomy or functional)
#'   containing sample-level proteomics metrics.
#'    
#' @param transform_mode character.
#'   Overall transformation strategy:
#'   \itemize{
#'     \item \code{"raw"}: Return data without any transformation.
#'     \item \code{"detection"}: Convert quantitative values to presence/absence (0/1).
#'     \item \code{"quantitative"}: Apply zero handling, log-transformation and
#'     normalization steps.
#'   }
#'    
#' @param zero_strategy character.
#'  Strategy to handle zeros prior to log-transformation:
#'   \itemize{
#'     \item \code{"keep"}: Keep zeros unchanged.
#'     \item \code{"na"}: Convert zeros to \code{NA}.
#'     \item \code{"pseudocount"}: Add a pseudocount to all values.
#'   }
#' 
#' @param pseudocount numeric. Default \code{0.1}.
#'   Pseudocount value added when \code{"pseudocount"} strategy is selected.
#'    
#' @param log_transform character.
#'   Logarithmic transformation applied after zero handling:
#'   \itemize{
#'     \item \code{"none"}: No log-transformation.
#'     \item \code{"log2"}: Log2 transformation.
#'     \item \code{"log10"}: Log10 transformation.
#'   }
#' 
#' @param normalization character.
#'   Normalization strategy applied after log-transformation:
#'   \itemize{
#'     \item \code{"none"}: No normalization applied.
#'     \item \code{"median"}: Median centering across samples.
#'   }
#'   
#' @param post_zero_strategy character.
#'   Zero-strategy applied after normalization:
#'   \itemize{
#'     \item \code{"keep"}: Retain \code{NA} values.
#'     \item \code{"zero"}: Convert \code{NA} values back to zero.
#'   }
#'  
#' @param verbose logical, default \code{TRUE}
#'   Whether to print informative messages describing each transformation step.
#'      
#' @return data.frame
#'   A data frame containing transformed sample-level metric values suitable
#'   for EDA and downstream analyses (e.g. PCA, clustering, differential expression).
#'    
#' @examples
#' # Protein-level intensity data, log2-transformed and median-normalized
#' protein_median_norm <- apply_data_transform(
#'   df = protein_processed,
#'   transform_mode = "quantitative",
#'   zero_strategy = "na",
#'   log_transform = "log2",
#'   normalization = "median",
#'   post_zero_strategy = "keep",
#'   verbose = TRUE
#' )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
apply_data_transform <- function(
    df,
    #biological_level,
    transform_mode = "quantitative",
    zero_strategy = "keep",
    pseudocount = 0.1,
    log_transform = "none",
    normalization = "none",
    post_zero_strategy = "keep",
    verbose = TRUE
) {
  
  # 1. Check if a data frame exits
  if (is.null(df)) {
    if (verbose) stop("[ERROR]: Data frame not found.")
    return(df)
  }
  
  # 2. Check if transformation option is available.
  if (!transform_mode %in% eda_data_options$transform_mode) {
    stop("[ERROR]: Transformation ", transform_mode, " not found. Avaliable options are ",  
         paste(eda_data_options$transform_mode, collapse = ", "))
  }

  # 3. No transformation.
  if (transform_mode == "raw") {
    if (verbose) message("[INFO]: No transformation applied.")
    return(df)
  }
  
  # 4. Apply absence/presence matrix
  if (transform_mode == "detection") {
    df <- df %>% mutate(across(-feature_id, ~ ifelse(. > 0, 1, 0)))
    if (verbose) message("[INFO]: Data was converted into absence/presence.")
    return(df)
  }
  
  # 5. Apply processing.
  if (transform_mode == "quantitative") {
    
    # 5.1 Apply zero handling
    df <- apply_zero_handling(df, strategy = zero_strategy, pseudocount = pseudocount)
    
    # 5.2 Apply log transform
    df <- apply_log_transform(df, log_transform = log_transform)
    
    # 5.3 Apply normalization
    df <- apply_normalization(df, norm_mode = normalization)
    
    # 5.4 Apply post-zero handling
    df <- apply_post_zero_handling(df, strategy = post_zero_strategy)
  }
  
  return(df)
  
}

# --------------------------------------------------------------------------------------------------------
#' apply_zero_handling
#'
#' Apply an initial zero-handling strategy prior to log-transformation.
#' 
#' @param df data.frame
#'   A biological-level feature table (peptide, protein, taxonomy or functional)
#'   containing sample-level proteomics metrics.
#'    
#' @param strategy character. Default \code{"keep"}
#'  Strategy to handle zeros prior to log-transformation:
#'   \itemize{
#'     \item \code{"keep"}: Keep zeros unchanged.
#'     \item \code{"na"}: Convert zeros to \code{NA}.
#'     \item \code{"pseudocount"}: Add a pseudocount to all values.
#'   }
#'    
#' @param pseudocount numeric. Default \code{0.1}.
#'   Pseudocount value added when \code{"pseudocount"} strategy is selected.
#' 
#' @return data.frame
#'   A data frame with zero values handled according to the selected strategy.
#'    
#' @examples
#' Convert protein-level zeros into NA prior to log-transformation
#' protein_zeros_na <- apply_zero_handling(
#'   df = proteins_processed,
#'   strategy = "na"
#' )
#'
#' # Add a pseudocount before log-transformation
#' protein_pseudo <- apply_zero_handling(
#'   df = proteins_processed,
#'   strategy = "pseudocount",
#'   pseudocount = 0.1
#' )
#'
#' @export
#'--------------------------------------------------------------------------------------------------------

apply_zero_handling <- function(df, strategy = "keep", pseudocount = 0.1) {
  
  # 1. Check if initial zero handling strategy option is available.
  if (!strategy %in% eda_data_options$zero_strategy) {
    stop("[ERROR]: Initial zero handling strategy ", strategy, " not found. Available options are ", 
         paste(eda_data_options$zero_strategy, collapse = ", "))
  }

  # 2. Initial zero handling options
  # 2.1 Not change
  if (strategy == "keep") {
    message("[INFO]: Zeros were ketp as zeros.")
    return(df)
  }
  
  # 2.2 Change 0 to NA
  if (strategy == "na") {
    df <- df %>% mutate(across(matches("^ID_\\d+_"), ~ ifelse(. == 0, NA, .)))
    message("[INFO]: Zeros were turned into NA previous log transformation.")
    return(df)
  }
  
  # 2.3 Add pseudocount if needed
  if (strategy == "pseudocount") {
    if (pseudocount < 0) {
      stop("[ERROR]: A negative pseudocount is not possible.")
    }
    if (!pseudocount %in% eda_data_options$pseudocount) {
      warning("[WARNING]: Common pseudocount optios are ", 
              paste0(eda_data_options$pseudocount, collapse = ","))
    }
  # TODO: Add different zero count strategy for intensities and spectral counts.
    message(paste0("[INFO]: A pseudocount of ", pseudocount,  "was added previous log transformation."))
    df <- df %>% mutate(across(matches("^ID_\\d+_"), ~ . + pseudocount))
    return(df)
  }
  
  }

# --------------------------------------------------------------------------------------------------------
#' apply_log_transform
#'
#' Apply logarithmic transformation to quantitative proteomics data.
#' 
#' This function applies a logarithmic transformation to sample-level
#' proteomics metrics after zero handling. Log-transformation is a standard
#' step in proteomics data analysis to stabilize variance and reduce the
#' influence of extreme values.
#'
#' @param df data.frame.
#'  A biological-level feature table (peptide, protein, taxonomy or functional)
#'   containing sample-level proteomics metrics.
#'   
#' @param log_transform character.
#'  Logarithmic transformation applied after zero handling:
#'   \itemize{
#'     \item \code{"none"}: No log-transformation.
#'     \item \code{"log2"}: Log2 transformation.
#'     \item \code{"log10"}: Log10 transformation.
#'   }
#' 
#' @return data.frame
#'   A data frame with log-transformed sample-level values according to the
#'   selected transformation.
#'  
#' @examples
#' Apply log2 transformation to protein-level data.
#' protein_log2 <- apply_log_transform(
#'   df = protein_processed,
#'   log_transform = "log2"
#'   )
#'
#' @export
#' # --------------------------------------------------------------------------------------------------------
apply_log_transform <- function(df, log_transform = "none") {
  
  # 1. Check if log transform option is available.
  if (!log_transform %in% eda_data_options$log_transform) {
    stop("[ERROR]: Log transformation ", log_transform, " not found. Available options are ", 
         paste(eda_data_options$log_transform, collapse = ", "))
  }
  
  # 2. Log transformations
  # 2.1 No transformation
  if (log_transform == "none") {
    message("[INFO]: NO log transformation was applied.")
    return(df)
  }
  
  if ((log_transform == "log2" | log_transform == "log10") & any(df == 0,na.rm = T)) {
    message("[WARNING]: Log transformation will be applied over 0: Inf values were returned.")
  }
  
  # 2.2 Log2 transformation  
  if (log_transform == "log2") {
    df <- df %>% mutate(across(matches("^ID_\\d+_"), ~ log2(.)))
    message("[INFO]: Log2 transformation was applied.")
    return(df)
    }
    
  # 2.3 Log10 transformation
  if (log_transform == "log10") {
    df <- df %>% dplyr::mutate(across(matches("^ID_\\d+_"), ~ log10(.)))
    message("[INFO]: Log10 transformation was applied.")
    return(df)
  }
  
}

# --------------------------------------------------------------------------------------------------------
#' apply_normalization
#'
#' Apply a normalization strategy to proteomics data.
#'
#' This function applies normalization to sample-level proteomics metrics
#' after log-transformation. Normalization helps to remove global intensity
#' differences between samples and is a standard step before multivariate
#' analyses (e.g., PCA, clustering).
#'
#' TODO: Add additional strategies such as quantile normalization.
#' 
#' @param df data.frame.
#'  A biological-level feature table (peptide, protein, taxonomy or functional)
#'  containing sample-level proteomics metrics.
#'    
#' @param norm_mode character.
#'   Normalization strategy applied after log-transformation:
#'   \itemize{
#'     \item \code{"none"}: No normalization applied.
#'     \item \code{"median"}: Median centering across samples.
#'   }
#'   
#' @return data.frame
#'   A data frame with normalized sample-level values according to the selected strategy.
#'   
#' @examples
#' # Median normalization of protein-level log2 intensities
#' protein_median_norm <- apply_normalization(
#'   df = proteins_processed,
#'   norm_mode = "median"
#' )
#'
#' @export
#' --------------------------------------------------------------------------------------------------------
apply_normalization <- function(df, norm_mode = "median") {
  
  # TODO: Add quantile normalization with preprocessCore R package.
  
  # 1. Check if normalization option is available.
  if (!norm_mode %in% eda_data_options$normalization) {
    stop("[ERROR]: Normalization  ", norm_mode, " not found. Available options are ", 
         paste(eda_data_options$normalization, collapse = ", "))
  }
  
  # 2. Normalization options.
  # 2.1 No Normalization
  if (norm_mode == "none"){
    message("[INFO]: NO normalization was applied.")
    return(df)
  }
  
  if (any(df == -Inf | df == Inf, na.rm = T)) {
    stop("[ERROR]: numeric columns contain Inf or -Inf. This usually happens if you log-transformed zeros. Replace zeros with a small value or NA before normalizing.")
  }
  
  # 2.2 Column-wise median centering normalization.
  if (norm_mode == "median"){
    df_norm <- df %>% dplyr::mutate(across(matches("^ID_\\d+_"), ~ . - median(., na.rm = TRUE)))
    message("[INFO]: Median normalization was applied.")
    return(df_norm)
  }
}

# --------------------------------------------------------------------------------------------------------
#' apply_post_zero_handling
#'
#' Apply a post-normalization zero-handling strategy.
#' 
#' This function modifies NA values after normalization according to the analysis needs.
#' Some downstream methods (e.g., PCA, clustering) may require zeros instead of NA, 
#' while others can handle NA values directly. This step ensures consistency with the 
#' chosen EDA or downstream workflow.
#' 
#' @param df data.frame.
#'   A normalized biological-level feature table (peptide, protein, taxonomy, or functional)
#'   containing sample-level proteomics metrics.
#'  
#' @param strategy 
#'  Zero-strategy applied after normalization:
#'   \itemize{
#'     \item \code{"keep"}: Retain \code{NA} values.
#'     \item \code{"zero"}: Convert \code{NA} values back to zero.
#'   }
#'  
#' @return data.frame
#'   A normalized data frame with sample-level values adjusted according to the
#'   selected post-zero strategy.
#'   
#' @examples
#' # Convert NA values into zeros after median normalization
#' protein_norm_zero <- apply_post_zero_handling(
#'   df = protein_processed_norm,
#'   strategy = "zero"
#' )
#'
#' @export
#'--------------------------------------------------------------------------------------------------------
apply_post_zero_handling <- function(df, strategy = "keep"){
  
  # Sometimes zeros or NA were needed for different purposes after normalization.
  
  # 1. Check if post zero handling strategy is available.
  if (!strategy %in% eda_data_options$post_zero_strategy) {
    stop("[ERROR]: Post-zero handling strategy  ", strategy, " not found. Available options are ", 
         paste(eda_data_options$post_zero_strategy, collapse = ", "))
  }
  
  # 2. Post-zero handling strategy options
  # 2.1 Zeros retained
  if (strategy == "keep") {
    message("[INFO]: NO post-zero handling after normalization.")
    return(df) 
  }
  
  # 2.2 NA -> ZERO
  if (strategy == "zero" & any(is.na(df))) {
    df <- df %>% mutate(across(matches("^ID_\\d+_"), ~ ifelse(is.na(.), 0, .)))
    message("[INFO]: NAs were turned into zeros")
    return(df)
  }

}