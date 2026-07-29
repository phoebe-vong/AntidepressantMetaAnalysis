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
# 7. Identify contrast IDs / experimental factors for a result set
# ------------------------------------------------------------------------------

get_contrast_ids_for_result_set <- function(fold_change_column_names) {
  broken_up <- strsplit(fold_change_column_names, "_")
  matrix_broken_up <- matrix(unlist(broken_up), ncol = 3, byrow = TRUE)
  matrix_broken_up[, 2]  # the contrast id piece
}

# Returns a list of everything the next steps need, instead of five separate
# <<- globals.
extract_de_results_for_contrasts <- function(de_results_good_annotation,
                                             contrasts_log2fc, contrasts_tstat,
                                             result_set_contrasts) {
  
  names_fc_cols    <- colnames(de_results_good_annotation)[colnames(de_results_good_annotation) %in% contrasts_log2fc]
  names_tstat_cols <- colnames(de_results_good_annotation)[colnames(de_results_good_annotation) %in% contrasts_tstat]
  
  contrast_ids_in_df <- get_contrast_ids_for_result_set(names_fc_cols)
  
  datasets_in_df <- result_set_contrasts$ExperimentID[result_set_contrasts$ContrastIDs %in% contrast_ids_in_df]
  gse_id <- datasets_in_df[1]
  
  factors_in_df <- result_set_contrasts$ExperimentalFactors[result_set_contrasts$ContrastIDs %in% contrast_ids_in_df]
  
  comparisons_of_interest <- paste(datasets_in_df, factors_in_df, sep = "_")
  
  list(
    names_fc_cols            = names_fc_cols,
    names_tstat_cols         = names_tstat_cols,
    gse_id                   = gse_id,
    comparisons_of_interest  = comparisons_of_interest
  )
}

# --- Example usage for section 7 ---
# extraction <- extract_de_results_for_contrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, result_set_contrasts)
# NamesOfFoldChangeColumns <- extraction$names_fc_cols
# NamesOfTstatColumns      <- extraction$names_tstat_cols
# GSE_ID                   <- extraction$gse_id
# ComparisonsOfInterest    <- extraction$comparisons_of_interest
#
# # If the auto-derived names are unwieldy, override manually here, e.g.:
# # ComparisonsOfInterest <- c("GSE270831_social isolation")
