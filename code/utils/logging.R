# --------------------------------------------------------------------------------------------------------
# Script: loggin.R
# Description: Metaproteomics data preprocessing Pipeline - Logging
# Author: Andrés Arroyo Barea
# Date: 2026-01-13
#       Purpose: Generate a global logging reporting system.
# --------------------------------------------------------------------------------------------------------
log_info <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg_eval <- glue::glue(msg)
  message(sprintf("[%s] [INFO] %s", timestamp, msg_eval))
}

log_warn <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg_eval <- glue::glue(msg)
  warning(sprintf("[%s] [WARN]: %s", timestamp, msg_eval), call. = FALSE)
}

log_error <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg_eval <- glue::glue(msg) <- glue::glue(msg)
  stop(sprintf("[%s] [ERROR]: %s", timestamp, msg_eval), call. = FALSE)
}