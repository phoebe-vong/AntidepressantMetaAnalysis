#Gemma Results Extractions - Team Antidepressant (EE/PA on hippocampus)
#Phoebe Vong
#2026-07-23
#wd = "G:/My Drive/BDA Code/GemmaResults"
#workspace = save.image("G:/My Drive/BDA Code/ProcessingGemmaResults/Github source files/AntidepressantTopic_GemmaProcessing_new.RData")

#####################################
setwd("G:/My Drive/BDA Code/GemmaResults")

library(gemma.R)

list.files()
ResultSets_toScreen<-read.csv("Antidepressant_ResultSets_toScreen - ResultSets_toScreen.csv", header=TRUE, stringsAsFactors = FALSE)

str(ResultSets_toScreen)

ResultSet_contrasts<-ResultSets_toScreen[ResultSets_toScreen$Include =="Y", ]

str(ResultSet_contrasts)

############################

DownloadingDEResults<-function(ResultSet_contrasts){
  
    UniqueResultSetIDs<-unique(ResultSet_contrasts$ResultSetIDs)
  
  print("These are the result sets that you identified as being of interest")
  print(UniqueResultSetIDs)

    
  differentials <- UniqueResultSetIDs %>% lapply(function(x){
   
    get_differential_expression_values(resultSet = x)[[1]]
  })
  
  str(differentials)

  
  
  missing_contrasts <- differentials %>% sapply(nrow) %>% {.==0}

  differentials <<- differentials[!missing_contrasts]
  UniqueResultSetIDs<<-UniqueResultSetIDs[!missing_contrasts]
  
  print("These are the result sets that had differential expression results:")
  print(UniqueResultSetIDs)
  
  print("Your differential expression results for each of your result sets are stored in the object named differentials. This object is structured as a list of data frames. Each element in the list represetns a result set, with the data frame containing the differential expression results")
  
  #We already identified which statistical contrasts we wanted during our screening: object with the specific contrast ids that we want:
  #ResultSet_contrasts$ContrastIDs
  
  #Which will be these columns within the listed dataframes of differential expression results:
  
  print("These are the columns for the effect sizes for our statistical contrasts of interest (Log(2) Fold Changes")
  Contrasts_Log2FC<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_log2fc", sep="")
  
  print(Contrasts_Log2FC)
  
  
  print("these are the columns for the T-statistics for our statistical contrasts of interest - we will use that information to derive the sampling variances")
  
  Contrasts_Tstat<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_tstat", sep="")
  
  print(Contrasts_Tstat)
  
}


library(gemma.R)
library(tidyr)

DownloadingDEResults(ResultSet_contrasts)


#################


SavingGemmaDEResults_forEachResultSet<-function(differentials, UniqueResultSetIDs, ResultSet_contrasts){
  
  for (i in c(1:length(differentials))){
    
    ThisResultSet<-UniqueResultSetIDs[i]
    
    ThisDataSet<-ResultSet_contrasts$ExperimentID[ResultSet_contrasts$ResultSetIDs==ThisResultSet][1] 
    
    write.csv(differentials[[i]], paste("DEResults", ThisDataSet, ThisResultSet, ".csv", sep="_"))
    
    rm(ThisDataSet, ThisResultSet)
  }
}












##################
###################
#########New version

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!requireNamespace("gemma.R", quietly = TRUE)) {
  BiocManager::install("gemma.R")
}

library("gemma.R", character.only = TRUE)

DownloadingDEResults<-function(ResultSet_contrasts){
  
  #To pull down the statistical results, we'll only want the unique result set ids:
  UniqueResultSetIDs<-unique(ResultSet_contrasts$ResultSetIDs)
  
  print("These are the result sets that you identified as being of interest")
  print(UniqueResultSetIDs)
  
  differentials <- UniqueResultSetIDs %>% lapply(function(x){
    #   # take the first and only element of the output. the function returns a list 
    #   # because single experiments may have multiple resultSets. Here we use the 
    #   # resultSet argument to directly access the results we need
    get_differential_expression_values(resultSet = x)[[1]]
  })
  
  names(differentials)<-UniqueResultSetIDs
  
  str(differentials)

  
  # # some datasets might not have all the advertised differential expression results
  # # calculated due to a variety of factors. here we remove the empty differentials
  missing_contrasts <- differentials %>% sapply(nrow) %>% {.==0}
  differentials <<- differentials[!missing_contrasts]
  UniqueResultSetIDs<<-UniqueResultSetIDs[!missing_contrasts]
  
  print("These are the result sets that had differential expression results:")
  print(UniqueResultSetIDs)
  
  print("Your differential expression results for each of your result sets are stored in the object named differentials. This object is structured as a list of data frames. Each element in the list represetns a result set, with the data frame containing the differential expression results")
  
  #Within any particular Result Set, there are likely to be some contrasts that we want and others that we don't want
  #For example, a result set might contain a variety of stress interventions
  #And maybe we only want the acute stress contrast results
  
  #We already identified which statistical contrasts we wanted during our screening:
  #This is the object with the specific contrast ids that we want:
  #ResultSet_contrasts$ContrastIDs
  
  #Which will be these columns within the listed dataframes of differential expression results:
  
  print("These are the columns for the effect sizes for our statistical contrasts of interest (Log(2) Fold Changes")
  Contrasts_Log2FC<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_log2fc", sep="")
  
  print(Contrasts_Log2FC)
  
  print("these are the columns for the T-statistics for our statistical contrasts of interest - we will use that information to derive the sampling variances")
  
  Contrasts_Tstat<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_tstat", sep="")
  
  print(Contrasts_Tstat)

}

apply


