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

# Raw data including Kostas previous results to compare and have the data
rawData <- read.csv("../../data/cleanedRawData144fishG1G2.csv")
rawData$SampleID <- rawData$ID

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
    fullMetadata <- merge(merge(merge(merge(methylBSdf, mappDatBSB), file), fileT), rawData)
  fullMetadata <- fullMetadata[!names(fullMetadata) %in% c("Sample.Name", "Sample")]
  names(fullMetadata)[names(fullMetadata) %in% "MappingEfficiency"] <- "MappingEfficiency%BSBoldvsGynogen"
  
  ## After exploration of raw data, we decide to remove 7 samples from (1) fam 12 (N=4, only in parents) and (2) with bad quality ("S12", "S22", "S110", "S118", "S142")
  fullMetadata <- fullMetadata[!fullMetadata$Family %in% "Fam12",]
  fullMetadata <- fullMetadata[!fullMetadata$ID %in% c("S12", "S22", "S110", "S118", "S142"),]
  
  nrow(fullMetadata)
  
  ## Export summary table:
  write.csv(fullMetadata, "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/fullMetadata135_Alice.csv", row.names=FALSE, quote=FALSE)
  
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
  pdf("../../data/fig/plot_comparison_aligners.pdf")
  plotComparisonAligners
  dev.off()
  
  png("../../data/fig/plot_comparison_aligners.png")
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
