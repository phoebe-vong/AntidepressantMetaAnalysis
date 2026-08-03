#Basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
#Updated version: March 10, 2026
#Updated again: July 28, 2026 (to make the code more generalizable for the full 2026 cohort)

######################

#Antidepressant Meta-Analysis
#Phoebe Vong
#07-29-2026
#wd = "G:/My Drive/BDA Code/MetaAnalysis"
#workspace = save.image("G:/My Drive/BDA Code/MetaAnalysis/AntidepressantTopic_MetaAnalysis_BDA_Environ (07-29-2026).RData")

#####################

#Grabbing input data and setting the working directory:

######################

#Installing and loading relevant code packages:

if (!require("metafor", quietly = TRUE)){
  install.packages("metafor")
}

library(metafor)

if (!require("plyr", quietly = TRUE)){
  install.packages("plyr")
}

library(plyr)

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("multtest")                                 

library(multtest)

######################

#Download the two necessary functions from our Github repository:

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_RunBasicMetaAnalysis.R

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_FalseDiscoveryCorrection.R

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_MakeForestPlots.R

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_VolcanoPlot.R

######################

#Read in the functions:

source("Function_RunBasicMetaAnalysis.R")

source("Function_FalseDiscoveryCorrection.R")

source("Function_MakeForestPlots.R")

source("Function_VolcanoPlot.R")

######################

#Example Usage:

#Figure out which columns contain differential expression output:
colnames(MetaAnalysis_FoldChanges)

#Then make a variable that is a vector of with the column numbers containing the DE output:
#Example code for making a variable representing the columns 
Columns_DE<-c(29:45)

#Then make a variable that is the shorthand annotation that we will use to label the genes in our output, e.g. MouseRat_GeneSymbol:
Column_GeneName<-46

NumberOfComparisons=16
CutOffForNAs=6


metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV, Columns_DE, Column_GeneName)
#Note: this function can take a while to run, especially if you have a lot of data  

#Take a peek at the output:

str(metaOutput)
head(metaOutput)
tail(metaOutput)

write.csv(metaOutput, "metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")

colnames(metaOutput)

write.csv(influence_dfbs, "influence_dfbs.csv")
write.csv(influence_cookd, "influence_cookd.csv" )
write.csv(influence_TF, "influence_TF.csv")

###############

#FDR correction:

FalseDiscoveryCorrection(metaOutput, MetaAnalysis_Annotation)

###############

#Writing out the meta-analysis input for record-keeping:

write.csv(MetaAnalysis_FoldChanges_ForMeta, "MetaAnalysis_FoldChanges_ForMeta.csv")
write.csv(MetaAnalysis_SV_ForMeta, "MetaAnalysis_SV_ForMeta.csv")





###############
###############
###############
###############
###############
###############




#Make forest plots:

setwd("G:/My Drive/BDA Code/MetaAnalysis/Forest Plots")


#Set the lower and upper limits for the x-axis:
Xaxis_LowerAndUpperBound<-c(-6,6)

#Set the height of the plot in inches (should be preferably <11.5)
HeightInInches<-7

#Making forest plots - example: eg. "Frem3" is gene name
#check top genes by MetaAnalysis Output FDR sorted by pvalue (high to low), check FDR column
#selected top 19 genes

MakeForestPlots(metaOutputFDR_annotated, "Cygb_Cygb", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Atg16l2_Atg16l2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Plekha6_Plekha6", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "H2bc4_NA", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Ociad2_Ociad2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Bbs10_Bbs10", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Akap12_Akap12", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Slc27a2_Slc27a2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Hook2_Hook2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Syt12_Syt12", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Smyd4_Smyd4", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Trpm2_Trpm2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Smyd3_Smyd3", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Gucy1a1_Gucy1a1", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

#error with this one
#in function MakeForestPlots.R: in rma(effect, var) : Processing terminated since k = 0.
MakeForestPlots(metaOutputFDR_annotated, "D16Ertd472e_C11h21orf9", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)


MakeForestPlots(metaOutputFDR_annotated, "Igfbp5_Igfbp5", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Tfrc_Tfrc", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Fam107a_Fam107a", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

MakeForestPlots(metaOutputFDR_annotated, "Fkbp1a_Fkbp1a", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

###############
###############
###############


#Make a volcano plot (unannotated version)

#Object containing differential expression or meta-analysis results:
DE_Results<-metaOutputFDR_OrderbyPval

#Categorical Variable of Interest:
VariableOfInterest<-"Lifestyle Interventions"

#Name of column containing Log2 Fold Change for the Variable of Interest:
CoefficientCol<-"Log2FC_estimate"

#Name of column containing p-value for the Variable of Interest:
PvalueCol<-"pval"

#Name of column containing FDR for the Variable of Interest:
FDRCol<-"FDR"

VolcanoPlot(DE_Results, VariableOfInterest, CoefficientCol, PvalueCol, FDRCol)


