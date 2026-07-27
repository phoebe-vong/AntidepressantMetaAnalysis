#This is R code for a practice search in GEO
#2026-07-27
# Stephanie Maciejewski

# ============================================================
# Install and load packages
# ============================================================

if (!require(BiocManager, quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

# ============================================================
# Search GEO
# ============================================================

# Uncomment to see function definition
# searchGEO

#Example search code:

MyQueryTerms <- paste(
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
  "AND ('Musmusculus"[ORGN] OR "Rattus norvegicus"[ORGN])", 
  "AND ("Expression profiling by high throughput sequencing"[DataSet Type]", 
  "OR "Expression profiling by array"[DataSet Type]) AND "gse"[Filter]"
  ")",
  sep = ""
)

QueryResults <- searchGEO(MyQueryTerms)

# This will show you an overview of the identified GEO Records:

str(QueryResults)

#Adding columns to hold the additional metadata to our Query Results object:

new_columns <- c(
  "Citation",
  "PMID",
  "Conributor",
  "Date",
  "Abstract"
)

QueryResults[new_columns] <- ""

#Looping over each of the identified GEO records and extracting the desired metadata:

for(i in seq_len(nrow(QueryResults))){
  
  gse_raw <- getGEO(QueryResults$`Series Accession`[i], GSEMatrix=FALSE)
  
  QueryResults$Citation[i] <- paste(Meta(gse_raw)$citation, collapse=" ")
  
  QueryResults$PMID[i] <- paste(Meta(gse_raw)$pubmed_id, collapse=" ")
  
  QueryResults$Contributor[i] <- paste(Meta(gse_raw)$contributor, collapse = " ")
  
  QueryResults$Date[i] <- paste(Meta(gse_raw)$submission_date, collapse= " ")
  
  QueryResults$Abstract[i] <- Meta(gse_raw)$summary

}

#Getting an overview of our Query Result object with its new additions:

str(QueryResults)

#Adding empty columns to hold additional information that we will find while reviewing the dataset records:

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

QueryResults[review_columns] <- ""

#Phrases that indicate the variable manipulated:
# Subjects were treated with ___ and rna-sequencing performed after testing behaviorally for...
# Subjects were divided into groups and one group experienced...
# Subjects received one of two interventions...

#Output the query results as a comma-separated variable file:

write.csv(QueryResults, "QueryResults.csv")

# ============================================================
# Retrieve metadata
# ============================================================

#How to extract the sample metadata in data frame format for a single GEO series:
#Use GEOQuery to pull down the full record (sample metadata and expression data)

gse_raw <- getGEO("GSE237890", GSEMatrix=TRUE)


#You can see all of the goodies stashed in this object using the structure function:
str(gse_raw)

#Grab the "expression set" component of the object (item #1)
eset <- gse_raw[[1]]

#Grab the "phenoData" (sample metadata) for the expression set:
metadata_df <- pData(eset)

#You can see all of the components in the phenoData using str:
str(metadata_df)

#View the column names for the phenoData:
colnames(metadata_df)

#Note: A lot of these column names have been auto named "characteristics_ch1..."
#To find out what those actually are, you can click on the data frame in your global environment (upper right)
#the actual variable name is stashed in the individual cells for the columns followed by a hyphen.

#View the first few rows of the phenoData using head:
head(metadata_df)

#Or you can write out the metadata as a .csv file and peruse it in a spreadsheet program:
getwd()

write.csv(metadata_df, "GSE237890_MetaData.csv")