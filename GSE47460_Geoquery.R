library(GEOquery)
library(affycoretools)
library(dplyr)
library(ggplot2)
library(Biobase)
library(limma)
library(tidyverse)
library(stringr)

myGeo1 <- getGEO("gse47460")[[1]] # pull out first half of expressionset
dim(myGeo1)# 429 samples, 15261 features

myGeo2 <- getGEO("gse47460")[[2]] # pull out second half of expressionset
dim(myGeo2)# 153 samples, 15261 features

# combine two set into one complete set for expression set, phenotype data and feature data
expressionset <- cbind(exprs(myGeo1), exprs(myGeo2))
dim(expressionset) # 15261 582
all(colnames(fData(myGeo1)) == colnames(fData(myGeo2))) # feature data column names match
PhenoData <- rbind(pData(myGeo1), pData(myGeo2))
dim(PhenoData) # 582 55
all(fData(myGeo1) == fData(myGeo2))  # has the same feature data 

# creat a new expression set including the two sets above
totalExpressionSet <- ExpressionSet(assayData = expressionset,
                                    phenoData = AnnotatedDataFrame(PhenoData),
                                    featureData = AnnotatedDataFrame(fData(myGeo1)))
dim(totalExpressionSet) # 15261 582
totalAssayData <- exprs(totalExpressionSet)
dim(totalAssayData)
class(totalAssayData)
dim(totalAssayData)
totalPhenoData <- pData(totalExpressionSet)
dim(totalPhenoData) # 582 55
sum(str_detect(totalPhenoData$title,"CTRL")==TRUE) # 108
sum(str_detect(totalPhenoData$title,"COPD")==TRUE) # 219
sum(str_detect(totalPhenoData$title,"ILD")==TRUE) # 255

# get the numerical expression values
expressionData <- exprs(totalExpressionSet)# pull out expression set as matrix
dim(expressionData)
colnames(expressionData)
rownames(expressionData) # probe ID
class(expressionData) # matrix
head(expressionData) 
expressionData[1:4, 1:5] # subset

# Get feature-level data for ProbIDs and Gene IDs
probInfo <- fData(totalExpressionSet) # genename and sequence
class(probInfo)# data frame
colnames(probInfo)
rownames(probInfo) # probe ID
geneSymbol <- probInfo[,7]# gene symbol for each probe
class(geneSymbol)
geneName <- probInfo[,8]# gene name for eaach probe

# Get the sample-level data for Age and Gender
sampleInfo <- pData(totalExpressionSet)
class(sampleInfo)# data frame
dim(sampleInfo)
colnames(sampleInfo)
patientCategory <- sampleInfo[,1]
patientGender <-  sampleInfo[,"characteristics_ch1.1" ]
patientAge <- sampleInfo[, "characteristics_ch1.2" ]
class(patientGender) # factor
class(patientCategory) # factor
class(patientAge) # factor
rownames(sampleInfo)
head(sampleInfo)

# test
rownames(expressionData) == rownames(probInfo)
rownames(sampleInfo) == colnames(expressionData)


# others from myGEO
sampleList <- c(sampleNames(totalExpressionSet))
sampleList
probeList <- c(featureNames(totalExpressionSet))
probeList
varLabels(totalExpressionSet) # columns of patients information table
annotation(totalExpressionSet)

# transpose expressionData to add in r
expressionData <- as.data.frame(expressionData)
expressionDataTranspose <- t(expressionData)
dim(expressionDataTranspose)

# change column names to gene symbol and row names to patient title
colnames(expressionDataTranspose) <- geneSymbol
rownames(expressionDataTranspose) <- patientCategory
dim(expressionDataTranspose)

# add age and gender 
expressionDataTranspose <- cbind(expressionDataTranspose,patientGender,patientAge)
dim(expressionDataTranspose)

# Cleaning up sample-level data for sex, age and adding condition column
expressionDataTranspose <- as.data.frame(expressionDataTranspose)
expressionDataTranspose$Condition <- gsub(".*_","",row.names(expressionDataTranspose))
expressionDataTranspose$patientGender <- gsub(".*-","",expressionDataTranspose$patientGender)
expressionDataTranspose$patientAge <- as.numeric(gsub(".* ","",expressionDataTranspose$patientAge))

## Analyzing differential expression with LIMMA
condition <- as.factor(expressionDataTranspose$Condition)
design <- model.matrix(~0+condition)
colnames(design) <- c("COPD","CTRL","ILD")

## Which comparisons will be made
cm <- makeContrasts (CTRLvILD = ILD - CTRL,
                     levels = design)
fit1 <- lmFit(totalExpressionSet, design = design)
fit2 <- contrasts.fit(fit1, contrasts = cm)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)
results_all <- topTable(fit2,adjust = "fdr", n=Inf)


#### Subsetting patients with IPF from ILD group based off external .csv file ####
IPF_Sample_Names <- read.csv("path-to-Patient_IDs_IPF.csv")
IPF_Sample_Names <- IPF_Sample_Names$Symbol
# Resetting sampleNames in totalExpressionSet to match formating in .csv file
sampleNames(totalExpressionSet) <- row.names(expressionDataTranspose) 
# Creating a subset of totalExpressionSet only containing IPF and Ctrl samples from .csv file
myGeo_IPF <- totalExpressionSet[, sampleNames(totalExpressionSet) %in% IPF_Sample_Names]
dim(myGeo_IPF) #15261  263 

## Analyzing differential expression with LIMMA
expressionDataTranspose_IPF <- subset(expressionDataTranspose, row.names(expressionDataTranspose) %in% IPF_Sample_Names)
condition <- as.factor(expressionDataTranspose_IPF$Condition)
design <- model.matrix(~0+condition)
colnames(design) <- c("CTRL","IPF")

## Which comparisons will be made
cm <- makeContrasts (CTRLvIPF = IPF - CTRL,
                     levels = design)
fit1 <- lmFit(myGeo_IPF, design = design)
fit2 <- contrasts.fit(fit1, contrasts = cm)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)
results_all_IPF <- topTable(fit2,adjust = "fdr", n=Inf)
