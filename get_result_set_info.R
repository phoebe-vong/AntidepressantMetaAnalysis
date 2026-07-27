# July 27, 2026
# Stephanie Maciejewski

# ---------------------------------------------------------
# Load libraries and files 
# ---------------------------------------------------------

# Starting with "GemmaExpressionInfo.csv"
# Added a column in that file called "Keep" and marked off Y or N, then resaved the file as: "GemmaExpressionInfo_Keep.csv"

gemma_expression_info <- read.csv("GemmaExpressionInfo_Keep.csv", header=TRUE, stringsAsFactors = FALSE)

str(gemma_expression_info)

experiment_ids <- gemma_expression_info$experiment_ids[gemma_expression_info$Keep=="Y"]

library(gemma.R)

# ---------------------------------------------------------
# Retrieve experimental design information
# ---------------------------------------------------------

get_result_set_info <- function(experiment_ids){
  
  #Making an empty data.frame to store results:
  
  ResultSets_toScreen<-data.frame(ExperimentID = character(),
                                  ResultSetIDs = character(), 
                                  ContrastIDs = character(), 
                                  ExperimentIDs = character(), 
                                  FactorCategory = character(), 
                                  ExperimentalFactors = character(), 
                                  BaselineFactors = character(), 
                                  Subsetted = logical(), 
                                  SubsetBy = character()
                                  )
  
  str(ResultSets_toScreen)
  
  # ---------------------------------------------------------
  # Extract contrasts and factors
  # ---------------------------------------------------------
  #We will then loop over each of the datasets:
  
  for(i in seq_along(experiment_ids)){
    
    message(
      "Processing ",
      i,
      "/",
      length(experiment_ids),
      ": ",
      experiment_ids[i]
    )
    
    #For each dataset, we will use Gemma's API to access the experimental design info:
    Design <- gemma.R::get_dataset_differential_expression_analyses(experiment_ids[i])
    
    if(nrow(Design)>0){
      #Next, we'll make some empty vectors to store the experimental factor and baseline factor information for each result id for the dataset:
      ExperimentalFactors <- character(length(Design$result.ID))
      BaselineFactors <- character(length(Design$result.ID))
      
      #We will then loop over each of the result ids for the dataset:
      for(j in seq_along(Design$result.ID)){
        
        #And grab the vector of experimental factors associated with that result id
        ExperimentalFactorVector <- Design$experimental.factors[[j]]$summary
        #And collapse that info down to a single entry that will fit in our data.frame
        ExperimentalFactors[j] <- paste(ExperimentalFactorVector, collapse="; ")
        
        #And then grab the vector of baseline/control/reference values associated with that result id
        BaselineFactorVector <- Design$baseline.factors[[j]]$summary
        #And collapse that info down to a single entry that will fit in our data.frame
        BaselineFactors[j] <- paste(BaselineFactorVector, collapse="; ")
      }
      
      #Some of the datasets are subsetted for the differential expression analyses
      #We will make an empty vector to store subset information for each result id
      SubsetBy <- character(length(Design$result.ID))
      
      #Then we will determine whether the dataset is subsetted:
      if(Design$isSubset[1]==TRUE){
        
        #If it is subsetted, we will loop over each result id for the dataset
        for (j in seq_along(Design$result.ID)){
          
          #And grab the vector of subsetting information
          SubsetByVector <- Design$subsetFactor[[j]]$summary
          
          #And then collapse that information down to a single entry that will fit in our dataframe
          SubsetBy[j] <- paste(SubsetByVector, collapse="; ")
        }  
        
        #if the dataset wasn't subsetted for the differential expression analysis:
      }else{
        #We'll just make a vector of NA values to put in the "Subsetted by" column
        SubsetBy <- rep(NA, length((Design$result.ID)))
      }
      
      #Then we combine all of the information for all of the result sets for the dataset into a dataframe
      ResultSets_ForExperiment <- cbind.data.frame(
        ExperimentID=rep(ExperimentIDs[i], length(Design$result.ID)),
        ResultSetIDs=Design$result.ID, 
        ContrastIDs=Design$contrast.ID, 
        ExperimentIDs=Design$experiment.ID, 
        FactorCategory=Design$factor.category, 
        ExperimentalFactors = ExperimentalFactors, 
        BaselineFactors = BaselineFactors, 
        Subsetted=Design$isSubset, 
        SubsetBy = SubsetBy)
      
      #And add that information as rows to our data frame including the result set information for all datasets:
      ResultSets_toScreen <- rbind.data.frame(ResultSets_toScreen, ResultSets_ForExperiment)
      
    }
    
  }
  
  # ---------------------------------------------------------
  # Store results
  # ---------------------------------------------------------
  n <- nrow(ResultSets_toScreen)
  
  Include <- character(n)
  WrongBaseline <- character(n)
  ResultsNotRegionSpecific <- character(n)
  ReAnalyze <- character(n)
  
  #And add them as columns to our dataframe:            
  ResultSets_toScreen <- cbind.data.frame(
    ResultSets_toScreen, 
    Include, 
    WrongBaseline, 
    ResultsNotRegionSpecific, 
    ReAnalyze)
  
  #And then write everything out as a .csv file that we can easily mark up in a spreadsheet program:
  write.csv(ResultSets_toScreen, "ResultSets_toScreen.csv")
  
  message("Wrote", 
          nrow(ResultSets_toScreen),
          " result sets to ResultSets_toScreen"
  )

}


#Applying the function:
get_result_set_info(experiment_ids)