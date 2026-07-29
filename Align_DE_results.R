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
# 9. Align per-dataset results into one big data.frame per species
# ------------------------------------------------------------------------------

# Both rat/mouse aligners follow the same shape; factored the shared logic
# out so there's only one place to fix bugs in the future.
.align_species_datasets <- function(list_of_de_results, gene_id_column_name) {
  
  fc_dfs <- lapply(list_of_de_results, function(x) {
    df <- data.frame(row.names(x[[1]]), x[[1]], stringsAsFactors = FALSE)
    names(df)[1] <- gene_id_column_name
    df
  })
  fold_changes <- join_all(fc_dfs, by = gene_id_column_name, type = "full")
  
  sv_dfs <- lapply(list_of_de_results, function(x) {
    df <- data.frame(row.names(x[[4]]), x[[4]], stringsAsFactors = FALSE)
    names(df)[1] <- gene_id_column_name
    df
  })
  sv <- join_all(sv_dfs, by = gene_id_column_name, type = "full")
  
  list(fold_changes = fold_changes, sv = sv)
}

align_rat_datasets <- function(list_of_rat_de_results) {
  .align_species_datasets(list_of_rat_de_results, "Rat_EntrezGene.ID")
}

align_mouse_datasets <- function(list_of_mouse_de_results) {
  .align_species_datasets(list_of_mouse_de_results, "Mouse_EntrezGene.ID")
}

# --- Example usage for section 9 ---
# rat_aligned <- align_rat_datasets(ListOfRatDEResults)
# Rat_MetaAnalysis_FoldChanges <- rat_aligned$fold_changes
# Rat_MetaAnalysis_SV          <- rat_aligned$sv
#
# mouse_aligned <- align_mouse_datasets(ListOfMouseDEResults)
# Mouse_MetaAnalysis_FoldChanges <- mouse_aligned$fold_changes
# Mouse_MetaAnalysis_SV          <- mouse_aligned$sv
