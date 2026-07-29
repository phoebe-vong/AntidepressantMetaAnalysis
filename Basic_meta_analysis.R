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
# 10. Run the random-effects meta-analysis, gene by gene
# ------------------------------------------------------------------------------

# Fixed bug: removed the dead code after the first return(); now returns a
# single named list instead of five separate <<- globals plus unreachable
# print() calls.
#
# Also: `Columns_DE` and `Column_GeneName` can now be found by name instead of
# hardcoded position (safer if the join order/column count ever changes) —
# see the example usage below.
run_basic_meta_analysis <- function(number_of_comparisons, cutoff_for_nas,
                                    meta_analysis_fold_changes, meta_analysis_sv,
                                    columns_de, column_gene_name) {
  
  n_stats <- 21  # keep in sync with the colnames() vector below if you add stats
  
  na_per_row <- apply(meta_analysis_fold_changes[, columns_de], 1, function(y) sum(is.na(y)))
  message("Table of # of NAs per row (gene):")
  print(table(na_per_row))
  
  fc_for_meta <- meta_analysis_fold_changes[na_per_row < cutoff_for_nas, ]
  sv_for_meta <- meta_analysis_sv[na_per_row < cutoff_for_nas, ]
  
  meta_output <- matrix(NA, nrow(fc_for_meta), n_stats)
  
  influence_dfbs  <- matrix(NA, nrow(fc_for_meta), ncol(fc_for_meta[columns_de]))
  influence_cookd <- matrix(NA, nrow(fc_for_meta), ncol(fc_for_meta[columns_de]))
  influence_tf    <- matrix(NA, nrow(fc_for_meta), ncol(fc_for_meta[columns_de]))
  
  colnames(influence_dfbs)  <- colnames(fc_for_meta)[columns_de]
  colnames(influence_cookd) <- colnames(fc_for_meta)[columns_de]
  colnames(influence_tf)    <- colnames(fc_for_meta)[columns_de]
  
  row.names(influence_dfbs)  <- fc_for_meta[, column_gene_name]
  row.names(influence_cookd) <- fc_for_meta[, column_gene_name]
  row.names(influence_tf)    <- fc_for_meta[, column_gene_name]
  
  for (i in seq_len(nrow(fc_for_meta))) {
    
    effect <- as.numeric(fc_for_meta[i, columns_de])
    var    <- as.numeric(sv_for_meta[i, columns_de])
    
    temp_meta <- tryCatch(rma(effect, var), error = function(e) NULL)
    
    if (!is.null(temp_meta)) {
      
      meta_output[i, 1]  <- temp_meta$b
      meta_output[i, 2]  <- temp_meta$se
      meta_output[i, 3]  <- temp_meta$pval
      meta_output[i, 4]  <- temp_meta$ci.lb
      meta_output[i, 5]  <- temp_meta$ci.ub
      meta_output[i, 6]  <- number_of_comparisons - sum(is.na(effect))
      meta_output[i, 7]  <- temp_meta$k
      meta_output[i, 8]  <- temp_meta$p
      meta_output[i, 9]  <- temp_meta$tau2
      meta_output[i, 10] <- temp_meta$se.tau2
      meta_output[i, 11] <- temp_meta$QE
      meta_output[i, 12] <- temp_meta$QEp
      meta_output[i, 13] <- temp_meta$I2
      meta_output[i, 14] <- temp_meta$H2
      
      pub_bias <- tryCatch(regtest(temp_meta), error = function(e) NULL)
      if (!is.null(pub_bias)) {
        meta_output[i, 15] <- pub_bias$zval
        meta_output[i, 16] <- pub_bias$pval
        meta_output[i, 17] <- pub_bias$dfs
      }
      
      robustness <- tryCatch(leave1out(temp_meta), error = function(e) NULL)
      if (!is.null(robustness)) {
        meta_output[i, 18] <- min(robustness$estimate, na.rm = TRUE)
        meta_output[i, 19] <- max(robustness$estimate, na.rm = TRUE)
        meta_output[i, 20] <- min(robustness$pval, na.rm = TRUE)
        meta_output[i, 21] <- max(robustness$pval, na.rm = TRUE)
      }
      
      influence_stats <- tryCatch(influence(temp_meta), error = function(e) NULL)
      if (!is.null(influence_stats)) {
        influence_dfbs[i, influence_stats$ids]  <- influence_stats$dfbs$intrcpt
        influence_cookd[i, influence_stats$ids] <- influence_stats$inf$cook.d
        influence_tf[i, influence_stats$ids]    <- influence_stats$is.infl
      }
    }
  }
  
  colnames(meta_output) <- c(
    "Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub",
    "Number_Of_Comparisons", "Number_of_Contrasts", "Number_of_Coefficients",
    "tau2_ResidualHeterogeneity", "SE_tau2_ResidualHeterogeneity",
    "QE_CochransQ_Teststat", "QEp_CochransQ_pval",
    "I2_PercentVar_TrueHeterogeneity", "H2_Ratio_EffectHetero_overSamplVar",
    "PubBias_Egger_Zstat", "PubBias_Egger_pval", "PubBias_Egger_DF",
    "Leave1Out_Min_Log2FC", "Leave1Out_Max_Log2FC",
    "Leave1Out_Min_Pval", "Leave1Out_Max_Pval"
  )
  row.names(meta_output) <- fc_for_meta[, column_gene_name]
  
  list(
    meta_output           = meta_output,
    meta_analysis_annotation = fc_for_meta[, -columns_de],
    influence_dfbs        = influence_dfbs,
    influence_cookd       = influence_cookd,
    influence_tf          = influence_tf,
    fc_for_meta           = fc_for_meta,
    sv_for_meta           = sv_for_meta
  )
}

# --- Example usage for section 10 ---
# colnames(MetaAnalysis_FoldChanges)  # inspect to confirm names below match your data
#
# # Prefer name-based selection over hardcoded positions where possible, e.g.:
# # Columns_DE <- grep("_log2fc$|_tstat$", colnames(MetaAnalysis_FoldChanges))  # adjust pattern to your join
# # Column_GeneName <- which(colnames(MetaAnalysis_FoldChanges) == "MouseRat_GeneSymbol")
# #
# # If your join doesn't leave a clean name to match on, hardcoded positions are
# # fine -- just re-derive them from colnames() every time you regenerate
# # MetaAnalysis_FoldChanges, rather than trusting old numbers.
# Columns_DE <- c(29:34)
# Column_GeneName <- 35
#
# meta_result <- run_basic_meta_analysis(
#   number_of_comparisons = 6,
#   cutoff_for_nas = 2,
#   meta_analysis_fold_changes = MetaAnalysis_FoldChanges,
#   meta_analysis_sv = MetaAnalysis_SV,
#   columns_de = Columns_DE,
#   column_gene_name = Column_GeneName
# )
#
# metaOutput <- meta_result$meta_output
# MetaAnalysis_Annotation <- meta_result$meta_analysis_annotation
#
# write.csv(metaOutput, "metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
# write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
# write.csv(meta_result$influence_dfbs, "influence_dfbs.csv")
# write.csv(meta_result$influence_cookd, "influence_cookd.csv")
# write.csv(meta_result$influence_tf, "influence_TF.csv")
# write.csv(meta_result$fc_for_meta, "MetaAnalysis_FoldChanges_ForMeta.csv")
# write.csv(meta_result$sv_for_meta, "MetaAnalysis_SV_ForMeta.csv")
#
# # FDR correction / forest plots / volcano plot: no changes 
# # source("Function_FalseDiscoveryCorrection.R")
# # FalseDiscoveryCorrection(metaOutput, MetaAnalysis_Annotation)