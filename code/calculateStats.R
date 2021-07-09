library(ggplot2)
library(ggrepel)

############ Data loading ############
# Kostas previous results to compare and have the data
Kostas <- readxl::read_xlsx("../data/Kostas_G2_info.xlsx")

# Raw reads quality check
file <- read.csv("../data/multiqc_report_rawReads_generalstats.csv") 

# Mapping efficiency after Bismark
mapp <- read.table(file = "../data/Report_mapping_efficiency.txt")
#####################################

mean(file$M.Seqs)

# 95% confidence interval
qnorm(0.975)*sd(file$M.Seqs)/sqrt(nrow(file))

file$SampleName = gsub("_R1_001", "", file$Sample.Name)



mappDat <- data.frame(SampleName = gsub("_R1_001_trimmed_cutadapt_bismark_bt2_SE_report.txt", "", mapp$V1),
                      MappingEfficiency = as.numeric(gsub("%", "", mapp$V4)))

hist(mappDat$MappingEfficiency, breaks = 100)

table(mappDat$MappingEfficiency < 60)

# Mean and 95% confidence interval
mean(mappDat$MappingEfficiency)
qnorm(0.975)*sd(mappDat$MappingEfficiency)/sqrt(nrow(mappDat))

## Both
AllDF <- merge(file, mappDat)
AllDF$SampleID <- sub('.*_', '', sub('_L00.*', '',  AllDF$SampleName))

PA <- ggplot(AllDF, aes(x=M.Seqs, y=MappingEfficiency, label = SampleID))+
  geom_point() +
  geom_label_repel() +
  theme_bw()+
  theme(legend.position = "none")

PA

# Not ideal samples would be:
AllDF[AllDF$MappingEfficiency < 60 | AllDF$M.Seqs < 5,]
notIdeal <- AllDF[AllDF$MappingEfficiency < 60 | AllDF$M.Seqs < 5,"SampleID"]

Kostas[Kostas$ID %in% notIdeal,]
# S12 is a NE-control offspring, to be removed for bad quality
# the 2 other samples with rather low read counts have an ok quality,
# and belong to 2 different groups (E_control, NE_exposed) so it should
# not bias the further steps too much