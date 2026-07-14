#Stephanie Maciejewski

pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, GEOQuery, BiocManager) 

########################

if (!require(BiocManager, quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

#protocols.io: quotation marks in code are a problem

########################

searchGEO

#Query terms:
MyQueryTerms<-'((environ*[All Fields] AND enrich*[All Fields]) OR (enrich*[All Fields] AND housing[All Fields]) OR (enrich*[All Fields] AND housed[All Fields]) OR (social*[All Fields] AND housing[All Fields]) OR (social*[All Fields] AND housed[All Fields]) OR "run"[All Fields] OR running[All Fields] OR exercis*[All Fields] OR wheel*[All Fields] OR toy*[All Fields] OR welfare[All Fields] OR (social[All Fields] AND enrich*[All Fields]) OR (sensor*[All Fields] AND enrich*[All Fields]) OR (motor[All Fields] AND enrich*[All Fields]) OR (cognitiv*[All Fields] AND enrich*[All Fields]) OR (behav*[All Fields] AND enrich*[All Fields]) OR (experienc*[All Fields] AND novel*[All Fields]) OR (environmen*[All Fields] AND novel*[All Fields]) OR (stimulat*[All Fields] AND novel*[All Fields]) OR (stimulat*[All Fields] AND environmen*[All Fields]) OR (stimulat*[All Fields] AND social*[All Fields]) OR (stimulat*[All Fields] AND cognitiv*[All Fields]) OR (stimulat*[All Fields] AND motor*[All Fields]) OR lifestyle[All Fields]) AND 
(hippocamp*[All Fields] OR "dentate gyrus"[All Fields] OR CA1[All Fields] OR CA2[All Fields] OR CA3[All Fields] OR CA4[All Fields] OR "CA field"[All Fields] OR subiculum[All Fields] OR fimbria[All Fields] OR "cornu ammonis"[All Fields]) AND 
("Mus musculus"[ORGN] OR "Rattus norvegicus"[ORGN]) AND 
("Expression profiling by high throughput sequencing"[DataSet Type] OR "Expression profiling by array"[DataSet Type]) AND 
"gse"[Filter]' 

QueryResults <- searchGEO(MyQueryTerms)

str(QueryResults)

# Identified GEO Records:

str(QueryResults)

#Adding columns to hold the additional metadata to the Query Results object:

QueryResults$Citation<-character(length=nrow(QueryResults))
QueryResults$PMID<-character(length=nrow(QueryResults))
QueryResults$Contributor<-character(length=nrow(QueryResults))
QueryResults$Date<-character(length=nrow(QueryResults))
QueryResults$Abstract<-character(length=nrow(QueryResults))

#Looping over each of the identified GEO records and extracting the desired metadata:

for(i in c(1:nrow(QueryResults))){
  
  gse_raw <- getGEO(QueryResults$`Series Accession`[i], GSEMatrix=FALSE)
  
  QueryResults$Citation[i] <- paste(Meta(gse_raw)$citation, collapse=" ")
  
  QueryResults$PMID[i] <- paste(Meta(gse_raw)$pubmed_id, collapse=" ")
  
  QueryResults$Contributor[i] <- paste(Meta(gse_raw)$contributor, collapse = " ")
  
  QueryResults$Date[i] <- paste(Meta(gse_raw)$submission_date, collapse= " ")
  
  QueryResults$Abstract[i] <- Meta(gse_raw)$summary
  
  rm(gse_raw)
}

#Getting an overview of the Query Result object with its new additions:

str(QueryResults)

#Adding empty columns to hold additional information that we will find while reviewing the dataset records:

QueryResults$Tissue<-character(length=nrow(QueryResults))
QueryResults$DevelopmentalStage<-character(length=nrow(QueryResults))
QueryResults$ManipulatedVariables<-character(length=nrow(QueryResults))
QueryResults$Notes<-character(length=nrow(QueryResults))

#Phrases that indicate the variable manipulated:
# Subjects were treated with ___ and rna-sequencing performed after testing behaviorally for...
# Subjects were divided into groups and one group experienced...
# Subjects received one of two interventions...

# Add some empty columns for taking inclusion/exclusion notes:

QueryResults$ManipulationUnrelatedToTopic<-character(length=nrow(QueryResults))
QueryResults$WrongTissue<-character(length=nrow(QueryResults))
QueryResults$NotBulkDissection_ParticularCellTypeOrSubRegion<-character(length=nrow(QueryResults))
QueryResults$IncorrectDevelopmentalStage<-character(length=nrow(QueryResults))
QueryResults$NotFullTranscriptome<-character(length=nrow(QueryResults))
QueryResults$MetadataIssues_MissingInfo_Retracted_Duplicated<-character(length=nrow(QueryResults))

QueryResults$Excluded<-character(length=nrow(QueryResults))
QueryResults$WhyExcluded<-character(length=nrow(QueryResults))

#Output the query results as a comma-separated variable file:

write.csv(QueryResults, "QueryResults2.csv")

# Check file output in working directory.
getwd()
