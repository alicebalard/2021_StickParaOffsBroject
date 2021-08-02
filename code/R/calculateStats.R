library(ggplot2)
library(ggrepel)

############ Data loading ############
# Kostas previous results to compare and have the data
Kostas <- readxl::read_xlsx("../../data/Kostas_G2_info.xlsx")

# Raw reads quality check
file <- read.csv("../../data/multiqc_general_stats_rawReads.txt") 

# Trimmed reads quality check
fileT <- read.csv("../../data/multiqc_general_stats_trimmedReadsCutadapt.txt") 

# Mapping efficiency after Bismark
mappBismarck <- read.table(file = "../../data/Report_mapping_efficiency_Bismark.txt")

# Mapping efficiency after BSBolt
mappBSBolt <- read.table(file = "../../data/Report_mapping_efficiency_BSBolt.txt")

############## Data formatting #######################
file$SampleName = gsub("_L00.*","",
                       gsub(".*-L1_","", file$Sample.Name))

mappDatBIS <- data.frame(SampleName = gsub("_L00*.","",
                                           gsub(".*-L1_","", 
                                                gsub("_R1_001_trimmed_cutadapt_bismark_bt2_SE_report.txt", "", mappBismarck$V1))),
                         MappingEfficiency = as.numeric(gsub("%", "", mappBismarck$V4)))

mappDatBSB <- data.frame(SampleName = gsub("_L00*.","",
                                           gsub(".*-L1_","",
                                                gsub("_R1_001_trimmed_cutadapt.fastq.gz", "", mappBSBolt$V1))),
                         MappingEfficiency = as.numeric(mappBSBolt$V3))

################# Data analyses #####################

# Nbr samples: 144
nrow(file)

# Mean nbr of million reads: 11.1
mean(file$M.Seqs)

# 95% confidence interval: 0.36
qnorm(0.975)*sd(file$M.Seqs)/sqrt(nrow(file))

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

plotComparisonAligners <- ggplot(AllDF, aes(x=MappingEfficiencyBIS, y=MappingEfficiencyBSB, label = SampleName))+
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

ggplot(AllDFBSB, aes(x=M.Seqs, y=MappingEfficiencyBSB, label = SampleName))+
  geom_point() +
  geom_label_repel() +
  theme_bw()+
  theme(legend.position = "none")

# Not good sample: S12 (baaad fastqc too)
AllDFBSB[AllDFBSB$SampleName %in% "S12",] # 68.96% mappability

AllDFBSB$X..Dups <- as.numeric(gsub("%", "", AllDFBSB$X..Dups))

ggplot(AllDFBSB, aes(x=X..Dups, y=MappingEfficiencyBSB, label = SampleName))+
  geom_point() +
  geom_label_repel() +
  theme_bw()+
  theme(legend.position = "none")

# S12 is a NE-control offspring, to be removed for bad quality
# the 2 other samples with rather low read counts have an ok quality,
# and belong to 2 different groups (E_control, NE_exposed) so it should
# not bias the further steps too much