#Gemma Histogram Extractions - Team Antidepressant (EE/PA on hippocampus)
#Phoebe Vong
# 07/16/2026
#wd = "G:/My Drive/BDA Code/Gemma Histograms"

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gemma.R")

library(gemma.R)

setwd("G:/My Drive/BDA Code/Gemma Histograms")
getwd()

list.files()

ExperimentIDs<-c("GSE135306", "GSE81672", "GSE180465", "GSE179667", "GSE146358", "GSE128255", "GSE187418")
ExperimentIDs

#header=TRUE as my file has no column title
ExperimentIDs<-read.csv("TeamAntidepressant_GSE - Sheet1.csv", header=TRUE, stringsAsFactors = FALSE)
str(ExperimentIDs)

#take first column only
ExperimentIDs<-ExperimentIDs[,1]
str(ExperimentIDs)


for(i in c(1:length(ExperimentIDs))){
  
  print(ExperimentIDs[i])
  

  if(inherits(try(get_dataset_processed_expression(ExperimentIDs[i]), silent=TRUE), "try-error")){
    
    print("Error: Not in Gemma")
    
  }else{
    
    Expression<-gemma.R::get_dataset_processed_expression(ExperimentIDs[i])
    
    ExpressionMatrix<-as.matrix(Expression[,-c(1:4)])
    
    pdf(paste(ExperimentIDs[i], "_Histogram.pdf", sep=""), height=4, width=4)
    
    hist(ExpressionMatrix, main="Histogram", xlab="Log2 Expression", col="green", cex.axis=1.3, cex.lab=1.3)
    
    dev.off()
    
    print(min(ExpressionMatrix, na.rm=TRUE))
    
    print(median(ExpressionMatrix, na.rm=TRUE))
    
    print(max(ExpressionMatrix, na.rm=TRUE))
    
    rm(Expression, ExpressionMatrix)
    
  }
  
}




##############################
##############################
##############################
############################## REVISED for data structure excel




#This is example code demonstrating how to pull down information from Gemma about...
## The distribution of the preprocessed data on Gemma
## The statistical contrasts performed during differential expression analysis by Gemma
# Phoebe Vong
# July 17, 2026



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gemma.R")

library(gemma.R)


ExperimentIDs<-read.csv("TeamAntidepressant_GSE - Sheet1.csv", header=TRUE, stringsAsFactors = FALSE)
str(ExperimentIDs)


ExperimentIDs<-ExperimentIDs[,1]
str(ExperimentIDs)
#This should be a vector now

GemmaExpressionInfo<-data.frame(ExperimentIDs=ExperimentIDs, InGemma=character(length=length(ExperimentIDs)), MinExpression=numeric(length=length(ExperimentIDs)), MedianExpression=numeric(length=length(ExperimentIDs)), MaxExpression=numeric(length=length(ExperimentIDs)))

str(GemmaExpressionInfo)


for(i in c(1:length(ExperimentIDs))){
  
  print(ExperimentIDs[i])
  
  if(inherits(try(get_dataset_processed_expression(ExperimentIDs[i]), silent=TRUE), "try-error")){
    
    print("Error: Not in Gemma")
    GemmaExpressionInfo$InGemma[i]<-"N"
    
  }else{
    
    GemmaExpressionInfo$InGemma[i]<-"Y"
    
    Expression<-gemma.R::get_dataset_processed_expression(ExperimentIDs[i])
    
    ExpressionMatrix<-as.matrix(Expression[,-c(1:4)])
    
    pdf(paste(ExperimentIDs[i], "_Histogram.pdf", sep=""), height=4, width=4)
    
    hist(ExpressionMatrix, main="Histogram", xlab="Log2 Expression", col="green", cex.axis=1.3, cex.lab=1.3)
    
    dev.off()
    
    print(min(ExpressionMatrix, na.rm=TRUE))
    
    GemmaExpressionInfo$MinExpression[i]<-min(ExpressionMatrix, na.rm=TRUE)
    
    print(median(ExpressionMatrix, na.rm=TRUE))
    
    GemmaExpressionInfo$MedianExpression[i]<-median(ExpressionMatrix, na.rm=TRUE)
    
    print(max(ExpressionMatrix, na.rm=TRUE))
    
    GemmaExpressionInfo$MaxExpression[i]<-max(ExpressionMatrix, na.rm=TRUE)
    
    rm(Expression, ExpressionMatrix)
    
  }
  
}

str(GemmaExpressionInfo)

write.csv(GemmaExpressionInfo, "AntidepressantGemmaExpressionInfo.csv")

#The histograms should be outputted in your working directory
# the summary stats will be printed in your console - you should save them somewhere in your working directory.





