# ------------------------------------------------------------------------------
# 0. Libraries (all in one place; installs are gated so re-running is cheap)
# ------------------------------------------------------------------------------
required_cran_packages <- c("metafor", "plyr", "dplyr")
for (pkg in required_cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("multtest", quietly = TRUE)) BiocManager::install("multtest")

library(metafor)
library(plyr)
library(dplyr)
library(multtest)
library(gemma.R)
# NOTE: this script also assumes GEOquery / GEOmetadb-style functions
# `searchGEO`, `getGEO`, `Meta` are already available in your session
# (load whichever package supplies them, e.g. `library(GEOquery)`).
# ------------------------------------------------------------------------------
# 4. Download differential expression results for the chosen contrasts
# ------------------------------------------------------------------------------

# Returns a list instead of writing several objects into the global
# environment via <<-, so the caller decides what to name things.
download_de_results <- function(result_set_contrasts) {
  
  unique_result_set_ids <- unique(result_set_contrasts$ResultSetIDs)
  message("Result sets of interest: ", paste(unique_result_set_ids, collapse = ", "))
  
  differentials <- lapply(unique_result_set_ids, function(x) {
    get_differential_expression_values(resultSet = x)[[1]]
  })
  
  missing_contrasts <- vapply(differentials, function(d) nrow(d) == 0, logical(1))
  
  differentials          <- differentials[!missing_contrasts]
  unique_result_set_ids  <- unique_result_set_ids[!missing_contrasts]
  
  message("Result sets with differential expression results: ",
          paste(unique_result_set_ids, collapse = ", "))
  
  contrasts_log2fc <- paste0("contrast_", result_set_contrasts$ContrastIDs, "_log2fc")
  contrasts_tstat  <- paste0("contrast_", result_set_contrasts$ContrastIDs, "_tstat")
  
  list(
    differentials         = differentials,
    unique_result_set_ids = unique_result_set_ids,
    contrasts_log2fc      = contrasts_log2fc,
    contrasts_tstat       = contrasts_tstat
  )
}

# --- Example usage for section 4 ---
# de_download <- download_de_results(result_set_contrasts)
# differentials          <- de_download$differentials
# unique_result_set_ids  <- de_download$unique_result_set_ids
# Contrasts_Log2FC       <- de_download$contrasts_log2fc
# Contrasts_Tstat        <- de_download$contrasts_tstat
