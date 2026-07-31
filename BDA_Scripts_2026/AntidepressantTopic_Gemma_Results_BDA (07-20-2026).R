#Gemma Results Extractions - Team Antidepressant (EE/PA on hippocampus)
#Phoebe Vong
# 07/20/2026
#wd = "G:/My Drive/BDA Code/GemmaResults"


setwd("G:/My Drive/BDA Code/GemmaResults")

library(gemma.R)


GemmaExpressionInfo_Keep<-read.csv("AntidepressantGemmaExpressionInfo_Keep.csv", header=TRUE, stringsAsFactors = FALSE)

str(GemmaExpressionInfo_Keep)

ExperimentIDs<-GemmaExpressionInfo_Keep$ExperimentIDs[GemmaExpressionInfo_Keep$Keep=="Y"]



######################################################
#setting up function

GettingResultSetInfoForDatasets<-function(ExperimentIDs){
  
  ResultSets_toScreen<-data.frame(ExperimentID="NA",ResultSetIDs="NA", ContrastIDs="NA", ExperimentIDs="NA", FactorCategory="NA", ExperimentalFactors="NA", BaselineFactors="NA", Subsetted=FALSE, SubsetBy="NA")
  
  str(ResultSets_toScreen)
  # 'data.frame':	1 obs. of  9 variables:
  # $ ExperimentID       : chr "NA"
  # $ ResultSetIDs       : chr "NA"
  # $ ContrastIDs        : chr "NA"
  # $ ExperimentIDs      : chr "NA"
  # $ FactorCategory     : chr "NA"
  # $ ExperimentalFactors: chr "NA"
  # $ BaselineFactors    : chr "NA"
  # $ Subsetted          : logi FALSE
  # $ SubsetBy           : chr "NA"
  
  #loop over each of the datasets:
  
  for(i in c(1:length(ExperimentIDs))){
    
    #For each dataset, we will use Gemma's API to access the experimental design info:
    Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])
    
    if(nrow(Design)>0){
      #Next, we'll make some empty vectors to store the experimental factor and baseline factor information for each result id for the dataset:
      ExperimentalFactors<-vector(mode="character", length(Design$result.ID))
      BaselineFactors<-vector(mode="character", length(Design$result.ID))
      
      #We will then loop over each of the result ids for the dataset:
      for(j in c(1:length(Design$result.ID))){
        
        #And grab the vector of experimental factors associated with that result id
        ExperimentalFactorVector<-Design$experimental.factors[[j]]$summary
        #And collapse that info down to a single entry that will fit in our data.frame
        ExperimentalFactors[j]<-paste(ExperimentalFactorVector, collapse="; ")
        
        #And then grab the vector of baseline/control/reference values associated with that result id
        BaselineFactorVector<-Design$baseline.factors[[j]]$summary
        #And collapse that info down to a single entry that will fit in our data.frame
        BaselineFactors[j]<-paste(BaselineFactorVector, collapse="; ")
      }
      
      SubsetBy<-vector(mode="character", length(Design$result.ID))
      
      if(Design$isSubset[1]==TRUE){
        
        for (j in c(1:length(Design$result.ID))){
          
          SubsetByVector<-Design$subsetFactor[[j]]$summary
          
          SubsetBy[j]<-paste(SubsetByVector, collapse="; ")
        }  
        
      }else{
       
        SubsetBy<-rep(NA, length((Design$result.ID)))
      }
      
      ResultSets_ForExperiment<-cbind.data.frame(ExperimentID=rep(ExperimentIDs[i],length(Design$result.ID)),ResultSetIDs=Design$result.ID, ContrastIDs=Design$contrast.ID, ExperimentIDs=Design$experiment.ID, FactorCategory=Design$factor.category, ExperimentalFactors, BaselineFactors, Subsetted=Design$isSubset, SubsetBy)
      
      ResultSets_toScreen<-rbind.data.frame(ResultSets_toScreen, ResultSets_ForExperiment)
      
      rm(ResultSets_ForExperiment, Design, ExperimentalFactors, BaselineFactors, SubsetBy)
      
    }else{
      rm(Design)
    }
    
  }
  
  ResultSets_toScreen<-ResultSets_toScreen[-1,]
  
  Include<-vector(mode="character", length=nrow(ResultSets_toScreen))
  WrongBaseline<-vector(mode="character", length=nrow(ResultSets_toScreen))
  ResultsNotRegionSpecific<-vector(mode="character", length=nrow(ResultSets_toScreen))
  ReAnalyze<-vector(mode="character", length=nrow(ResultSets_toScreen))
           
  ResultSets_toScreen<-cbind.data.frame(ResultSets_toScreen, Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)
  
    write.csv(ResultSets_toScreen, "AntidepressantResultSets_toScreen.csv")
  
  print("The Result Sets for your Datasets have been outputted into AntidepressantResultSets_toScreen.csv")
  print(str(ResultSets_toScreen))

  rm(Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)
}


###################################################
#Applying the function:

GettingResultSetInfoForDatasets(ExperimentIDs)
