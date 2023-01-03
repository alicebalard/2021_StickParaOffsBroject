## Alice Balard
## November 2021 (updated Dec22)

# Each script sources the previous script of the pipeline if needed
source("R00_rawDataCleaning.R")

#########################
## Data loading & prep ##
#########################
# Raw data including Kostas previous results to compare and have the data
rawData <- read.csv("../../data/cleanedRawData132fishG1G2.csv")
rawData$SampleID <- rawData$ID

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
methylBSdf <- read.delim("../../data/BSBolt_Methcall_Methstats_nbrMethylatedCpGperSample.txt", header = FALSE)
names(methylBSdf) <- c("Sample_Name", "NbrMethylatedCpG")
methylBSdf$SampleID <- gsub("_L00*.","", gsub(".*-L1_","", gsub("_R1_001", "", methylBSdf$Sample_Name)))

## Merge metadata:
fullMetadata <- merge(merge(merge(merge(methylBSdf, mappDatBSB), file), fileT), rawData)
fullMetadata <- fullMetadata[!names(fullMetadata) %in% c("Sample.Name", "Sample")]
names(fullMetadata)[names(fullMetadata) %in% "MappingEfficiency"] <- "MappingEfficiency%BSBoldvsGynogen"

# Mapping stats:
## Number of raw reads:
mean(fullMetadata$M.Seqs_rawReads) # 11.02M reads
qnorm(0.975)*sd(fullMetadata$M.Seqs_rawReads)/sqrt(nrow(fullMetadata)) # +/-0.38

## Mapping efficiency with BSBolt:
mean(fullMetadata$MappingEfficiency) # 85.6
qnorm(0.975)*sd(fullMetadata$MappingEfficiency)/sqrt(nrow(fullMetadata)) # +/- 0.51

## After exploration of raw data (multiQC files), we decide to remove 4 samples
# with bad quality and/or less than 6M reads after trimming ("S12", "S22", "S110", "S118", "S142")
fullMetadata <- fullMetadata[!fullMetadata$ID %in% c("S12", "S22", "S110", "S118", "S142"),]

sum(table(fullMetadata$Generat))
table(fullMetadata$Generat)
# N=127: 111 offspring + 16 parents

##############################
## Prepare sub metadatasets ##
##############################

## Change one term: NbrMethylatedCpG is the one calculated by BSBolt IN TOTAL
names(fullMetadata)[names(fullMetadata) %in% "NbrMethylatedCpG"] <- "NbrMethylatedCpG_global_BSBolt"

# give a numerical value to treatment, for Methylkit
fullMetadata$trtG1G2_NUM <- as.numeric(as.factor(fullMetadata$trtG1G2))

## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))

## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)

# paternal exposure
fullMetadata$PAT="Exposed father group"
fullMetadata$PAT[fullMetadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"

## Add brother pairs
fullMetadata$brotherPairID <- sapply(str_split(fullMetadata$clutch.ID, "_"), "[", 2 )
## avoid confusion with numeric
fullMetadata$brotherPairID <- paste0("BP",fullMetadata$brotherPairID)

########################
## Parents only metadata
fullMetadata_PAR <- fullMetadata[fullMetadata$Generat %in% "P",]

##########################
## Offspring only metadata
fullMetadata_OFFS <- fullMetadata[fullMetadata$Generat %in% "O",]
fullMetadata_OFFS$trtG1G2 <- droplevels(fullMetadata_OFFS$trtG1G2)

## Create variable for offsping and parents separated
fullMetadata_OFFS$offsTrt <- "controlO"
fullMetadata_OFFS$offsTrt[fullMetadata_OFFS$Tr %in% c("TT", "CT")] <- "infectedO"
fullMetadata_OFFS$patTrt <- "controlP"
fullMetadata_OFFS$patTrt[fullMetadata_OFFS$Tr %in% c("TC", "TT")] <- "infectedP"

## Sanity check
table(fullMetadata_OFFS$offsTrt, fullMetadata_OFFS$trtG1G2)
table(fullMetadata_OFFS$patTrt, fullMetadata_OFFS$trtG1G2)

## REORDER metadata by sample ID
fullMetadata = fullMetadata[order(as.numeric(gsub("S", "", fullMetadata$SampleID))),]
fullMetadata_PAR = fullMetadata_PAR[order(as.numeric(gsub("S", "", fullMetadata_PAR$SampleID))),]
fullMetadata_OFFS = fullMetadata_OFFS[order(as.numeric(gsub("S", "", fullMetadata_OFFS$SampleID))),]

## Export summary table:
write.csv(fullMetadata, "../../data/fullMetadata127_Alice.csv", row.names=FALSE, quote=FALSE)

#######################################
## Comparison between the 2 aligners ##
#######################################
rerun = FALSE
if (rerun == TRUE){
  
  # Plot mapping efficiency by reads number:
  hist(mappDatBIS$MappingEfficiency, breaks = 100)
  hist(mappDatBSB$MappingEfficiency, breaks = 100)
  
  # Mean and 95% confidence interval: 71.8 +/-0.59 for Bismark
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
  
  ggplot(AllDFBSB, aes(x=M.Seqs_rawReads, y=MappingEfficiencyBSB, label = SampleID))+
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

## clean workspace
rm(file, fileT, mappDatBIS, mappBismarck, mappBSBolt, mappDatBSB, methylBSdf, rawData)
