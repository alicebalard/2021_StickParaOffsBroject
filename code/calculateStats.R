library(ggplot2)
library(ggrepel)

file <- read.csv("../data/multiqc_report_rawReads_generalstats.csv") 
mean(file$M.Seqs)

# 95% confidence interval
qnorm(0.975)*sd(file$M.Seqs)/sqrt(nrow(file))

file$SampleName = gsub("_R1_001", "", file$Sample.Name)

# Mapping efficiency after Bismark
mapp <- read.table(file = "../data/Report_mapping_efficiency.txt")

mappDat <- data.frame(SampleName = gsub("_R1_001_trimmed_bismark_bt2_SE_report.txt", "", mapp$V1),
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

# Kostas previous results to compare
Kostas <- read.csv("G2_info_Kostas.csv")
## he marked in the xls sheet for either low reads or bad mapping: 
Kostas$tagged <- "white"
Kostas$tagged[Kostas$ID %in% c("S110", "S142", "S118", "S22", "S12", "S31", "S41", "S73", "S36")] <- "red"

PK <- ggplot(Kostas, aes(x=MReads, y=Mapping, label = ID))+
  geom_point() +
  geom_label_repel(aes(fill=tagged)) +
  scale_fill_manual(values = c("red", "white")) +
  theme_bw() +
  theme(legend.position = "none")

library(cowplot)
cowplot::plot_grid(PA, PK, labels = c("Alice", "Kostas"), nrow = 2)

a <- AllDF[AllDF$SampleID %in% Kostas$ID[Kostas$tagged %in% "red"], c("SampleID","M.Seqs", "MappingEfficiency")]
a <-a[order(as.numeric(gsub("S","", a$SampleID))),]

b <- Kostas[Kostas$tagged %in% "red", c("ID", "MReads", "Mapping")]
names(b) <- names(a)
b <-b[order(as.numeric(gsub("S","", b$SampleID))),]

a$obs <- "Alice"
b$obs <- "Kostas"
cbind(a,b)

a$SampleID

