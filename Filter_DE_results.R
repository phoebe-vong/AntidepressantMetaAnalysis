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
# 6. Filter DE results down to rows with usable gene annotation
# ------------------------------------------------------------------------------

# Returns the filtered data.frame instead of assigning via <<-.
filter_de_results_good_annotation <- function(de_results, write_csv = TRUE,
                                              output_path = "DE_Results_GoodAnnotation.csv") {
  
  message("# of rows in results: ", nrow(de_results))
  message("# of rows with missing NCBI annotation: ",
          sum(de_results$NCBIid == "" | de_results$NCBIid == "null"))
  message("# of rows with NA NCBI annotation: ", sum(is.na(de_results$NCBIid)))
  message("# of rows with missing Gene Symbol annotation: ",
          sum(de_results$GeneSymbol == "" | de_results$GeneSymbol == "null"))
  message("# of rows mapped to multiple NCBI IDs: ", length(grep("\\|", de_results$NCBIid)))
  message("# of rows mapped to multiple Gene Symbols: ", length(grep("\\|", de_results$GeneSymbol)))
  
  no_na <- de_results[
    !(de_results$NCBIid == "" | de_results$NCBIid == "null") & !is.na(de_results$NCBIid),
  ]
  
  multi_mapped <- grep("\\|", no_na$NCBIid)
  good_annotation <- if (length(multi_mapped) == 0) no_na else no_na[-multi_mapped, ]
  
  message("# of rows with good annotation: ", nrow(good_annotation))
  
  if (write_csv) write.csv(good_annotation, output_path, row.names = FALSE)
  
  good_annotation
}

# --- Example usage for section 6 ---
# DE_Results <- differentials[[1]]
# DE_Results_GoodAnnotation <- filter_de_results_good_annotation(DE_Results)
