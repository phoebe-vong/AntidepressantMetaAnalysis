# Grabbing information from Gemma about
# the distribution of the preprocessed data on Gemma
# the statistical contrasts performed during differential expression analysis by Gemma
# Stephanie Maciejewski
# July 27, 2026

# Load the Gemma.R package
library(gemma.R)

# ---------------------------------------------------------
# Read GEO accession IDs
# ---------------------------------------------------------
# Provide a vector containing your GEO Accession IDs
# ExperimentIDs<-c("GSE135306", "GSE81672", "GSE180465", "GSE179667", "GSE146358", "GSE128255", "GSE187418")

# Or you can copy them into a spreadsheet as a single column 
# and save them as a .csv file (preferably in the same order 
# as your inclusion/exclusion spreadsheet, without a column 
# name (header)) and read them in:
experiment_ids <- read.csv("YourGEOAccessionIDs.csv", 
                           header=FALSE, 
                           stringsAsFactors = FALSE)

str(experiment_ids)

#This is a dataframe, but the following code needs a vector, 
# so let's grab the first (and only) column and make it a vector.

experiment_ids <- experiment_ids[,1]
str(experiment_ids)

#This should be a vector now

# ---------------------------------------------------------
# Initialize results table
# ---------------------------------------------------------
n <- length(experiment_ids)

gemma_expression_info <- data.frame(experiment_ids=experiment_ids, 
                                    in_gemma=character(n), 
                                    min_expression=numeric(n), 
                                    median_expression=numeric(n), 
                                    max_expression=numeric(n)
)

str(gemma_expression_info)

# ---------------------------------------------------------
# Retrieve processed expression from Gemma
# ---------------------------------------------------------
#Let's loop over each of the GEO Accession IDs:

for(i in seq_along(experiment_ids)) {
  
  #To track what we are doing, let's print which GEO Accession ID we are on:
  message(
    "Processing",
    i,
    "/",
    n,
    ": ",
    experiment_ids[i]
    )
  
  #Let's test whether that GEO Accession ID is in the Gemma database. If it isn't, let's just print an error and skip to the end of the loop:
  
  expression <- try(
    get_dataset_processed_expression(experiment_ids[i]),
    silent = TRUE
  )
  
  if (inherits(expression, "try-error")) {
    
    message("Not in Gemma")
    gemma_expression_info$in_gemma[i] <- "N"
    
  } else {
    
    gemma_expression_info$in_gemma[i] <- "Y"
    
    expression_matrix <- as.matrix(expression[, -c(1:4)])
    
    pdf(paste0(experiment_ids[i], "_Histogram.pdf"))
    
    hist(expression_matrix, main="Histogram", xlab="Log2 Expression", col="green", cex.axis=1.3, cex.lab=1.3)
    
    dev.off()
    
    min_expr <- min(expression_matrix, na.rm = TRUE)
    median_expr <- median(expression_matrix, na.rm = TRUE)
    max_expr <- max(expression_matrix, na.rm = TRUE)
    
    message("Min expression", min_expr)
    message("Median expression", median_expr)
    message("Max expression", max_expr)
    
    gemma_expression_info$min_expression[i] <- min_expr
    gemma_expression_info$median_expression[i] <- median_expr
    gemma_expression_info$max_expression[i] <- max_expr
  }
  
}

# ---------------------------------------------------------
# Save summary results
# ---------------------------------------------------------
str(gemma_expression_info)

write.csv(gemma_expression_info, "GemmaExpressionInfo.csv")

# The histograms should be outputted in your working directory
# the summary stats will be printed in your console,
# you should save them somewhere in your working directory.