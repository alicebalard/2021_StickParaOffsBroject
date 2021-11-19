## Alice Balard
## November 2021

library(ggplot2)
library(ggrepel)
library(readxl)
library(tidyverse)
library(lme4) # linear mixed model
library(nlme) # linear mixed model with p-values
library(emmeans) # post-hoc Tukey test
library(ggsignif) # plot significance

#########################
## Data loading & prep ##
#########################

# Kostas previous results to compare and have the data
Kostas <- readxl::read_xlsx("../../data/Kostas_G2_info.xlsx")
Kostas$SampleID <- Kostas$ID

rerun = FALSE

if (rerun==TRUE){
  ## Raw reads quality check
  file <- read.csv("../../data/multiqc_general_stats_rawReads.txt") 
  file$SampleID = gsub("_L00.*","", gsub(".*-L1_","", file$Sample.Name))
  names(file) <- c("Sample.Name", "percent_duplicates_rawReads", "percent_GC_rawReads", "M.Seqs_rawReads", "SampleID")
  
  ## Trimmed reads quality check
  fileT <- read.table("../../data/multiqc_general_stats_trimmedReadsCutadapt.txt", header = T) 
  names(fileT) <- gsub("FastQC_mqc.generalstats.fastqc.", "", names(fileT))
  fileT$SampleID = gsub("_L00.*","", gsub(".*-L1_","", fileT$Sample))
  names(fileT) <- c("Sample", "percent_duplicates_trimmedReads", "percent_gc_trimmedReads", "avg_sequence_length_trimmedReads", "percent_fails_trimmedReads", "total_sequences_trimmedReads", "SampleID") 
  
  ## Mapping efficiency after Bismark
  mappBismarck <- read.table(file = "../../data/Report_mapping_efficiency_Bismark.txt")
  mappDatBIS <- data.frame(SampleID = gsub("_L00*.","", gsub(".*-L1_","", gsub("_R1_001_trimmed_cutadapt_bismark_bt2_SE_report.txt", "", mappBismarck$V1))),
                           MappingEfficiency = as.numeric(gsub("%", "", mappBismarck$V4)))
  
  ## Mapping efficiency after BSBolt
  mappBSBolt <- read.table(file = "../../data/Report_mapping_efficiency_BSBolt.txt")
  mappDatBSB <- data.frame(SampleID = gsub("_L00*.","", gsub(".*-L1_","", gsub("_R1_001_trimmed_cutadapt.fastq.gz", "", mappBSBolt$V1))),
                           MappingEfficiency = as.numeric(mappBSBolt$V3))
  
  ## Nbr of methylated sites:
  methylBSdf <- read.delim("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_stats/nbrMethylatedCpGperSample.txt", header = FALSE)
  names(methylBSdf) <- c("Sample_Name", "NbrMethylatedCpG")
  methylBSdf$SampleID <- gsub("_L00*.","", gsub(".*-L1_","", gsub("_R1_001", "", methylBSdf$Sample_Name)))
  
  ## Merge metadata:
  fullMetadata <- merge(merge(merge(merge(methylBSdf, mappDatBSB), file), fileT),Kostas)
  fullMetadata <- fullMetadata[!names(fullMetadata) %in% c("Sample.Name", "Sample")]
  names(fullMetadata)[names(fullMetadata) %in% "MappingEfficiency"] <- "MappingEfficiency%BSBoldvsGynogen"
  
  ## After exploration of raw data, we decide to remove 7 samples from (1) fam 12 (N=4, only in parents) and (2) with bad quality (S12, S118, S142)
  fullMetadata <- fullMetadata[!fullMetadata$Family %in% "Fam12",]
  fullMetadata <- fullMetadata[!fullMetadata$ID %in% c("S12", "S118", "S142"),]
  
  nrow(fullMetadata)
  
  ## Export summary table:
  write.csv(fullMetadata, "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/fullMetadata137_Alice.csv", row.names=FALSE, quote=FALSE)
  
}

#######################################
## Comparison between the 2 aligners ##
#######################################

DO_COMPA = FALSE

if (DO_COMPA==TRUE){
  # Plot mapping efficiency by reads number:
  hist(mappDatBIS$MappingEfficiency, breaks = 100)
  hist(mappDatBSB$MappingEfficiency, breaks = 100)
  
  # Mean and 95% confidence interval: 71.9 +/-0.59 for Bismark
  mean(mappDatBIS$MappingEfficiency)
  qnorm(0.975)*sd(mappDatBIS$MappingEfficiency)/sqrt(nrow(mappDatBIS))
  
  ## 85.4 +/- 0.51 for BSBolt
  mean(mappDatBSB$MappingEfficiency)
  qnorm(0.975)*sd(mappDatBSB$MappingEfficiency)/sqrt(nrow(mappDatBSB))
  
  ## plot:
  AllDFBIS <- merge(file, mappDatBIS)
  AllDFBSB <- merge(file, mappDatBSB)
  
  AllDFBIS$MappingEfficiencyBIS <- AllDFBIS$MappingEfficiency
  AllDFBSB$MappingEfficiencyBSB <- AllDFBSB$MappingEfficiency
  AllDFBSB <- AllDFBSB[,-6]
  AllDFBIS <- AllDFBIS[,-6]
  
  AllDF <- merge(AllDFBIS, AllDFBSB)
  
  plotComparisonAligners <- ggplot(AllDF, aes(x=MappingEfficiencyBIS, y=MappingEfficiencyBSB, label = SampleID))+
    geom_point() +
    geom_label_repel() +
    theme_bw()+
    theme(legend.position = "none") +
    geom_abline(slope=1, intercept=0, col = "red") +
    coord_cartesian(xlim = c(48,90), ylim = c(48,90)) +
    xlab("Mapping efficiency Bismark (% aligned reads)")+
    ylab("Mapping efficiency BSBolt (% aligned reads)")
  pdf("plot_comparison_aligners.pdf")
  plotComparisonAligners
  dev.off()
  
  png("plot_comparison_aligners.png")
  plotComparisonAligners
  dev.off()
  
  ggplot(AllDFBSB, aes(x=M.Seqs, y=MappingEfficiencyBSB, label = SampleID))+
    geom_point() +
    geom_label_repel() +
    theme_bw()+
    theme(legend.position = "none")
  
  # Not good sample: S12 (baaad fastqc too)
  AllDFBSB[AllDFBSB$SampleName %in% "S12",] # 68.96% mappability
  
  AllDFBSB$X..Dups <- as.numeric(gsub("%", "", AllDFBSB$X..Dups))
  
  ggplot(AllDFBSB, aes(x=X..Dups, y=MappingEfficiencyBSB, label = SampleID))+
    geom_point() +
    geom_label_repel() +
    theme_bw()+
    theme(legend.position = "none")
  
  # S12 is a NE-control offspring, to be removed for bad quality
  # the 2 other samples with rather low read counts have an ok quality,
  # and belong to 2 different groups (E_control, NE_exposed) so it should
  # not bias the further steps too much
}

#######################################################
## Nbr/Ratio of Methylated Sites in different groups ##
#######################################################
fullMetadata <- read.csv("../../data/fullMetadata137_Alice.csv")
## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))
## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)

## Nbr samples: 137
nrow(fullMetadata)

# Mean nbr of million reads: 11.2
mean(fullMetadata$M.Seqs_rawReads)
# 95% confidence interval: 0.35
qnorm(0.975)*sd(fullMetadata$M.Seqs_rawReads)/sqrt(nrow(fullMetadata))

# Average mapping efficiency +/-SD = 85.5% +/-0.47
mean(fullMetadata$MappingEfficiency.BSBoldvsGynogen)
qnorm(0.975)*sd(fullMetadata$MappingEfficiency.BSBoldvsGynogen)/sqrt(nrow(fullMetadata))

# Is there a correlation between nbr of methylated sites and coverage? YES
cor.test(fullMetadata$NbrMethylatedCpG,fullMetadata$M.Seqs_rawReads, method="pearson")
## t = 14.296, df = 135, p-value < 2.2e-16; cor = 0.78
ggplot(fullMetadata, aes(x=NbrMethylatedCpG, y=M.Seqs_rawReads, fill = trtG1G2))+
  geom_point(pch=21, size = 3, alpha =.8)+
  theme_minimal()

#############################################################
## Does RMS (ratio of methylated sites) differ per trt group? 
## NO, and our new methylation calling + reference genome does not allow to reproduce 
## Kostas results (even in a reduced dataset)
fullMetadata$RMS <- fullMetadata$NbrMethylatedCpG / (fullMetadata$M.Seqs_rawReads*10e6)

# simple lm, compare to null model: not signif
anova(lm(RMS ~ outcome, data = fullMetadata)) 

# with family as random factor: not significant
mod <- lmer(RMS ~ outcome + (1|Family), data = fullMetadata)
mod0 <- lmer(RMS ~ 1 + (1|Family), data = fullMetadata)
anova(mod, mod0)

ggplot(fullMetadata, aes(x = trtG1G2, y = RMS)) +
  geom_boxplot() +
  facet_wrap(.~Family)+
  theme_minimal()

## Infected and control all together by trt group
ggplot(fullMetadata, aes(x = outcome, y = RMS, fill =outcome)) +
  geom_boxplot() +
  scale_fill_manual(values=c("grey","red")) +
  theme_minimal()

#####################################################################
## Compare fitness traits between the different offsprings groups ###
## Follow up of Sagonas 2020 & Ferre Ortega's master's dissertation #
#####################################################################
fullMetadata_offs <- fullMetadata[fullMetadata$Generat %in% "O",]
fullMetadata_offs$trtG1G2 <- droplevels(fullMetadata_offs$trtG1G2)

## Create variable for offsping and parents separated
fullMetadata_offs$offsTrt <- "controlO"
fullMetadata_offs$offsTrt[fullMetadata_offs$Tr %in% c("TT", "CT")] <- "infectedO"
fullMetadata_offs$patTrt <- "controlP"
fullMetadata_offs$patTrt[fullMetadata_offs$Tr %in% c("TC", "TT")] <- "infectedP"
## Sanity check
table(fullMetadata_offs$offsTrt, fullMetadata_offs$trtG1G2)
table(fullMetadata_offs$patTrt, fullMetadata_offs$trtG1G2)

## Kaufmann et al. 2014: Body condition of the G2 fish, an estimate of fish health and a predictor
# of energy reserves and reproductive success, was calculated using there residuals from the 
# regression of body mass on body length (Chellappaet al.1995).
fullMetadata_offs$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|Family), data=fullMetadata_offs))

## Effect of paternal treatment on body condition of offspring:
## Kaufmann et al. 2014:
# To investigate in which way paternal G1 exposure affected offspring tolerance,
# we tested how the relationship between G2 body condition and infection intensity
# was affected by paternal G1 exposure. This was tested in a linear mixed model on 
# G2 body condition with paternal G1 treatment and the interaction between 
# paternal G1 treatment and G2 infection intensity as fixed effects. Maternal 
# half-sibship identity was set as a random effect

## Effect of paternal exposure on tolerance:
modTol <- lme(BCI ~ patTrt + patTrt:No.Worms,random=~1|Family,data=fullMetadata_offs)
anova(modTol)

myBCdf <- fullMetadata_offs %>% group_by(patTrt, No.Worms) %>% 
  summarise(BCI = mean(BCI)) %>% data.frame()
ggplot(fullMetadata_offs, aes(x=No.Worms, y = BCI, group = patTrt, col = patTrt))+
  geom_point() + geom_line(data=myBCdf)+
  geom_point(data=myBCdf, aes(fill = patTrt), col = "black", size = 3, pch = 21)+
  scale_color_manual(values = c("gray", "red"))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_bw()

## Effect of treatment groups of offspring on body condition:
## Kaufmann et al. 2014:
# The linear mixed effect model (nlme function in R) included G2 body condition as dependent variable, 
# sex, G2 treatment (exposed vs. control), paternal G1 treatment (exposed vs. control) 
# and their interactions as fixed effects as well as maternal G2 half-sibship identity as a random effect
mod1 <- lme(BCI ~ offsTrt * patTrt, random=~1|Family,data=fullMetadata_offs)
anova(mod1) # strong significant effect of both offspring trt & paternal + interactions

mod1.2 <- lme(BCI ~  trtG1G2, random=~1|Family,data=fullMetadata_offs)
## pairwise posthoc test
emmeans(mod1.2, list(pairwise ~ trtG1G2), adjust = "tukey")
## Control father - treatment offspring has a strongly significantly lower BC than 
## every other group, same as Kaufmann et al. 2014

ggplot(fullMetadata_offs, aes(x=trtG1G2, y = BCI))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("NE_control", "NE_exposed")),
              map_signif_level=TRUE, annotations="***",
              y_position = 150, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons = list(c("NE_exposed", "E_control")),
              map_signif_level=TRUE, annotations="***",
              y_position = 200, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons = list(c("NE_exposed", "E_exposed")),
              map_signif_level=TRUE, annotations="***",
              y_position = 250, tip_length = 0, vjust=0.4) +
  theme_bw()


