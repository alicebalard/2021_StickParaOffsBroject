## Differential methylation analyses
## A. Balard
## November 2021

library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)
library(ggsignif) ## for significance bars on ggplot
library(lme4) ## for mixed models
library(nlme) ## for mixed models
library(tidyverse)
library(emmeans) ## for post-hoc Tukey tests

## load custom functions
source("customRfunctions.R")

## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknownChr.RData")
## For further analyses: CpG covered in at least 2 then 6 individuals per group
load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov2_woSexAndUnknownChr.RData")
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_woSexAndUnknowChr.RData")

## Load samples metadata
## Metylation metadata merged with Kostas files
fullMetadata <- read.csv("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/fullMetadata137_Alice.csv")
## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))
## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)
#####################################################



################# PCA
## PCA analysis on our samples: plot a scree plot for importance of components
p=PCASamples(uniteCovALL_N137_final, screeplot=TRUE, obj.return = T) # first axis very important
s=summary(p)
#create scree plot
library(ggplot2)
qplot(c(1:137), s$importance[2,]) + 
  geom_line() + 
  xlab("Principal Component") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L), name = "Variance Explained")+
  ggtitle("Scree Plot") +
  theme_bw()

## PCA: validate the absence of strong sex influence on the methylation
uniteCovALL_N137_final_SEX = uniteCovALL_N137_final
uniteCovALL_N137_final_SEX@treatment = as.numeric(as.factor(metadata$Sex))
PCASamples(uniteCovALL_N137_final_SEX); title(sub="colored by sex")
PCASamples(uniteCovALL_N137_final_SEX, comp = c(3,4)); title(sub="colored by sex")

## PCA: check the Family influence on the methylation pattern
uniteCovALL_N137_final_FAM = uniteCovALL_N137_final
uniteCovALL_N137_final_FAM@treatment = as.numeric(as.factor(metadata$Family))
PCASamples(uniteCovALL_N137_final_FAM); title(sub="colored by family")
PCASamples(uniteCovALL_N137_final_FAM, comp = c(3,4)); title(sub="colored by family")

## PCA: check the Pattern treatment influence on the methylation pattern
uniteCovALL_N137_final_PAT = uniteCovALL_N137_final
metadata$PAT="Exposed father group"
metadata$PAT[metadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"
uniteCovALL_N137_final_PAT@treatment = as.numeric(as.factor(metadata$PAT))
PCASamples(uniteCovALL_N137_final_PAT); title(sub="colored by paternal exposure group")
PCASamples(uniteCovALL_N137_final_PAT, comp=c(3,4)); title(sub="colored by paternal exposure group")



