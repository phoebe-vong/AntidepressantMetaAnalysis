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
# 8. Collapse to one Log2FC / Tstat / SE / SV value per gene
# ------------------------------------------------------------------------------

# No more setwd()/dir.create() side effects: pass an explicit output directory
# and use file.path() throughout, so a mid-function error can't leave your R
# session in the wrong working directory.
collapse_de_results_one_per_gene <- function(gse_id, de_results_good_annotation,
                                             comparisons_of_interest,
                                             names_fc_cols, names_tstat_cols,
                                             output_dir = ".") {
  
  gse_dir <- file.path(output_dir, gse_id)
  if (!dir.exists(gse_dir)) dir.create(gse_dir, recursive = TRUE)
  
  fc_avg_list    <- vector("list", length(names_fc_cols))
  tstat_avg_list <- vector("list", length(names_fc_cols))
  se_avg_list    <- vector("list", length(names_fc_cols))
  
  for (i in seq_along(names_fc_cols)) {
    
    fc_col    <- de_results_good_annotation[[names_fc_cols[i]]]
    tstat_col <- de_results_good_annotation[[names_tstat_cols[i]]]
    se_col    <- fc_col / tstat_col
    
    fc_avg_list[[i]]    <- tapply(fc_col, de_results_good_annotation$NCBIid, mean)
    tstat_avg_list[[i]] <- tapply(tstat_col, de_results_good_annotation$NCBIid, mean)
    se_avg_list[[i]]    <- tapply(se_col, de_results_good_annotation$NCBIid, mean)
  }
  
  fc_by_gene    <- do.call(cbind, fc_avg_list)
  tstat_by_gene <- do.call(cbind, tstat_avg_list)
  se_by_gene    <- do.call(cbind, se_avg_list)
  
  colnames(fc_by_gene)    <- comparisons_of_interest
  colnames(tstat_by_gene) <- comparisons_of_interest
  colnames(se_by_gene)    <- comparisons_of_interest
  
  sv_by_gene <- se_by_gene ^ 2
  
  write.csv(fc_by_gene,    file.path(gse_dir, "DE_Results_GoodAnnotation_FoldChange_AveragedByGene.csv"))
  write.csv(tstat_by_gene, file.path(gse_dir, "DE_Results_GoodAnnotation_Tstat_AveragedByGene.csv"))
  write.csv(se_by_gene,    file.path(gse_dir, "DE_Results_GoodAnnotation_SE_AveragedByGene.csv"))
  write.csv(sv_by_gene,    file.path(gse_dir, "DE_Results_GoodAnnotation_SV.csv"))
  
  list(Log2FC = fc_by_gene, Tstat = tstat_by_gene, SE = se_by_gene, SV = sv_by_gene)
}

# --- Example usage for section 8 ---
# DEResults_thisGSE <- collapse_de_results_one_per_gene(
#   GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest,
#   NamesOfFoldChangeColumns, NamesOfTstatColumns
# )
# assign(paste0("DEResults_", GSE_ID), DEResults_thisGSE)  # if you want a per-GSE named object
