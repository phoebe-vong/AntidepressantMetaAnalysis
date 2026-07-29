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
# 2. Check which experiments are available in Gemma + basic expression stats
# ------------------------------------------------------------------------------

# `experiment_ids`: character vector of GEO accessions (e.g. from
# read.csv("YourGEOAccessionIDs.csv", header = FALSE, stringsAsFactors = FALSE)[, 1])
#
# `histogram_dir`: folder to write the per-experiment histogram PDFs into
# (created if it doesn't exist). Passed explicitly instead of writing to
# whatever the current working directory happens to be.
check_gemma_expression <- function(experiment_ids, histogram_dir = "histograms") {
  
  if (!dir.exists(histogram_dir)) dir.create(histogram_dir, recursive = TRUE)
  
  n <- length(experiment_ids)
  
  gemma_expression_info <- data.frame(
    experiment_ids   = experiment_ids,
    in_gemma          = character(n),
    min_expression    = numeric(n),
    median_expression = numeric(n),
    max_expression    = numeric(n),
    stringsAsFactors  = FALSE
  )
  
  for (i in seq_along(experiment_ids)) {
    
    message("Processing ", i, "/", n, ": ", experiment_ids[i])
    
    expression <- try(get_dataset_processed_expression(experiment_ids[i]), silent = TRUE)
    
    if (inherits(expression, "try-error")) {
      
      message("Not in Gemma")
      gemma_expression_info$in_gemma[i] <- "N"
      
    } else {
      
      gemma_expression_info$in_gemma[i] <- "Y"
      
      expression_matrix <- as.matrix(expression[, -c(1:4)])
      
      pdf(file.path(histogram_dir, paste0(experiment_ids[i], "_Histogram.pdf")))
      hist(
        expression_matrix,
        main = "Histogram", xlab = "Log2 Expression",
        col = "green", cex.axis = 1.3, cex.lab = 1.3
      )
      dev.off()
      
      min_expr    <- min(expression_matrix, na.rm = TRUE)
      median_expr <- median(expression_matrix, na.rm = TRUE)
      max_expr    <- max(expression_matrix, na.rm = TRUE)
      
      message("Min expression: ", min_expr)
      message("Median expression: ", median_expr)
      message("Max expression: ", max_expr)
      
      gemma_expression_info$min_expression[i]    <- min_expr
      gemma_expression_info$median_expression[i] <- median_expr
      gemma_expression_info$max_expression[i]    <- max_expr
    }
  }
  
  gemma_expression_info
}

# --- Example usage for section 2 ---
# experiment_ids <- read.csv("YourGEOAccessionIDs.csv", header = FALSE, stringsAsFactors = FALSE)[, 1]
# gemma_expression_info <- check_gemma_expression(experiment_ids)
# write.csv(gemma_expression_info, "GemmaExpressionInfo.csv", row.names = FALSE)
#
# --> Manual step: add a "Keep" column (Y/N) to GemmaExpressionInfo.csv,
#     save as e.g. "GemmaExpressionInfo_Keep.csv"
