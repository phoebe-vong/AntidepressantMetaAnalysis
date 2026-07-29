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
# 1. Build and run the GEO search query
# ------------------------------------------------------------------------------

build_geo_query <- function() {
  paste(
    "(",
    "(environ*[All Fields] AND enrich*[All Fields])",
    "OR (enrich*[All Fields] AND housing[All Fields])",
    "OR (enrich*[All Fields] AND housed[All Fields])",
    "OR (social*[All Fields] AND housing[All Fields])",
    "OR (social*[All Fields] AND housed[All Fields])",
    "OR 'run'[All Fields]",
    "OR running[All Fields]",
    "OR exercis*[All Fields]",
    "OR wheel*[All Fields]",
    "OR toy*[All Fields]",
    "OR welfare[All Fields]",
    "OR (social[All Fields] AND enrich*[All Fields])",
    "OR (sensor*[All Fields] AND enrich*[All Fields])",
    "OR (motor[All Fields] AND enrich*[All Fields])",
    "OR (cognitiv*[All Fields] AND enrich*[All Fields])",
    "OR (behav*[All Fields] AND enrich*[All Fields])",
    "OR (experienc*[All Fields] AND novel*[All Fields])",
    "OR (environmen*[All Fields] AND novel*[All Fields])",
    "OR (stimulat*[All Fields] AND novel*[All Fields])",
    "OR (stimulat*[All Fields] AND environmen*[All Fields])",
    "OR (stimulat*[All Fields] AND social*[All Fields])",
    "OR (stimulat*[All Fields] AND cognitiv*[All Fields])",
    "OR (stimulat*[All Fields] AND motor*[All Fields])",
    "OR lifestyle[All Fields]) AND (hippocamp*[All Fields]",
    "OR 'dentate gyrus'[All Fields]",
    "OR CA1[All Fields]",
    "OR CA2[All Fields]",
    "OR CA3[All Fields]",
    "OR 'cornu ammonis'[All Fields])",
    "AND ('Musmusculus'[ORGN] OR 'Rattus norvegicus'[ORGN])",
    "AND ('Expression profiling by high throughput sequencing'[DataSet Type]",
    "OR 'Expression profiling by array'[DataSet Type]) AND 'gse'[Filter]",
    ")",
    sep = ""
  )
}

# Enriches a searchGEO() results data.frame with citation/PMID/contributor/
# date/abstract metadata pulled per-record from getGEO(). Returns a NEW
# data.frame rather than mutating in place, so it's obvious at the call site
# what you're working with.
enrich_geo_metadata <- function(query_results) {
  
  metadata_columns <- c("Citation", "PMID", "Contributor", "Date", "Abstract")
  query_results[metadata_columns] <- ""
  
  for (i in seq_len(nrow(query_results))) {
    
    gse_raw <- getGEO(query_results$`Series Accession`[i], GSEMatrix = FALSE)
    
    query_results$Citation[i]    <- paste(Meta(gse_raw)$citation, collapse = " ")
    query_results$PMID[i]        <- paste(Meta(gse_raw)$pubmed_id, collapse = " ")
    query_results$Contributor[i] <- paste(Meta(gse_raw)$contributor, collapse = " ")
    query_results$Date[i]        <- paste(Meta(gse_raw)$submission_date, collapse = " ")
    query_results$Abstract[i]    <- paste(Meta(gse_raw)$summary, collapse = " ")  # now collapsed, like the others
  }
  
  query_results
}
add_manual_review_columns <- function(query_results) {
  review_columns <- c(
    "Tissue",
    "DevelopmentalStage",
    "ManipulatedVariables",
    "Notes",
    "ManipulationUnrelatedToTopic",
    "WrongTissue",
    "NotBulkDissection_ParticularCellTypeOrSubRegion",
    "IncorrectDevelopmentalStage",
    "NotFullTranscriptome",
    "MetadataIssues_MissingInfo_Retracted_Duplicated",
    "Excluded",
    "WhyExcluded"
  )
  query_results[review_columns] <- ""
  query_results
}

# --- Example usage for section 1 ---
# MyQueryTerms <- build_geo_query()
# QueryResults <- searchGEO(MyQueryTerms)
# str(QueryResults)
# QueryResults <- enrich_geo_metadata(QueryResults)
# QueryResults <- add_manual_review_columns(QueryResults)
# write.csv(QueryResults, "QueryResults.csv", row.names = FALSE)
#
# --> Manual step: open QueryResults.csv, fill in the review columns by hand,
#     and save your final list of GEO accessions as a single-column, no-header
#     CSV, e.g. "YourGEOAccessionIDs.csv"
