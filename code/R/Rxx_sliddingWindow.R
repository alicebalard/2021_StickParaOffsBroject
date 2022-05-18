## A. Balard
## May 2022
## Slidding window analysis: detecting peaks of methylation 
###########################################################
machine="mythinkpad" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
loadannot = FALSE # load genome annotations
source("R02.3_DATALOAD.R")

#################################
# Starting data set: methylations at CpG covered in all fish (N=55530CpG)
uniteObj = uniteCovALL_woSexAndUnknowChr

# calculate % methylation: each col is a sample, each row a position
perc_uniteObj = percMethylation(uniteObj, rowids = TRUE)

#################################
## 1. With average methylation ##
#################################
# calculate average methylation per treatment group at each position
perc_uniteObj_ave = calcAveMeth(perc_uniteObj)

# simplify names of CpG windows
rownames(perc_uniteObj_ave) <- sub('^([^.]+.[^.]+).*', '\\1', rownames(perc_uniteObj_ave))

# order positions by chromosomes & position
ordered_perc_uniteObj <- apply(
  perc_uniteObj_ave, 2, reorderByChrom)

# Check order of chromosomes
unique(sapply(strsplit(row.names(ordered_perc_uniteObj),"\\."), `[`, 1))

# Slide over 100bp complete windows, shift of 20bp
SWresults <- apply(ordered_perc_uniteObj, 2, 
                   function(x) unlist(slide(x, mean, .after = 100, .step = 20, .complete = TRUE)))

nrow(SWresults)# 2772 slidding windows

# Define annotation colors
mycolors <- c("black", "red", colOffs)
# names(mycolors) <- levels(metadata$Treatment)
names(mycolors) <- gsub("ave","", colnames(SWresults))
mycolors <- list(Treatment = mycolors)

# keep only names of the first chromosome, to display on heatmap
SWresults_hm = SWresults
rownames(SWresults_hm) <- data.frame(x = rownames(SWresults_hm)) %>%
  separate(x, c("A", "chr", "C")) %>%
  group_by(chr) %>%
  dplyr::mutate(numbering = row_number()) %>%
  dplyr::mutate(printRow = if_else(numbering == 1, chr, "")) %>% 
  .$printRow

# print heatmap
pheatmap(t(SWresults_hm), cluster_cols=F)

# print plot
# SWresults_mp=melt(SWresults)
# ggplot(SWresults_mp, aes(x=Var1, y=value))+
#   geom_line(aes(group=Var2, col=Var2)) +
#   scale_color_manual(values=c("black", "red", colOffs))

####################################
## 2. With individual methylation ##
####################################
# Starting data set: methylations at CpG covered in all fish (N=55530CpG)
uniteObj

# Tile for DMR calculation
tiles = tileMethylCounts(uniteObj, win.size=100,step.size=20,cov.bases = 10)

# Select only parents to get parDMR:
tilesG1 = reorganize(methylObj = tiles,
                     treatment = tiles@treatment[
                       tiles@sample.ids %in% fullMetadata$SampleID[fullMetadata$Generat %in% "P"]],
                     sample.ids = tiles@sample.ids[
                       tiles@sample.ids %in% fullMetadata$SampleID[fullMetadata$Generat %in% "P"]])

DMRG1 = getDiffMeth(tilesG1, fullMetadata[fullMetadata$Generat %in% "P",], mccores=3, mydif = 5)
  
# Select these tiles for G1 & G2
tiles_atDMRG1 = methylKit::select(tiles, 
                                  which(paste(tiles$chr, tiles$start) %in%
                                          paste(DMRG1$chr, DMRG1$start)))

# calculate % methylation: each col is a sample, each row a position
perc_uniteObj2 = percMethylation(tiles_atDMRG1, rowids = TRUE)

# calculate average methylation per treatment group at each parDMS
perc_uniteObj_ave_parDMS = calcAveMeth(perc_uniteObj2)

# Reorder by chromosome
ordered_perc_uniteObj2_ave <- apply(
  perc_uniteObj2_ave, 2, reorderByChrom)

# print heatmap
pheatmap(t(ordered_perc_uniteObj2_ave), cluster_cols=F)


ordered_perc_uniteObj2_ave




  
  ## bonus
  # if individuals are considered instead of grouping per treatment:
  # # Prepare metadata to compare the observed structure with an expected one: add colors for trt on the side of heatmap
  # metadata = data.frame(Treatment = fullMetadata$trtG1G2[match(colnames(perc_uniteObj), fullMetadata$SampleID)],
  #                       BP = fullMetadata$brotherPairID[match(colnames(perc_uniteObj), fullMetadata$SampleID)])
  # rownames(metadata) = fullMetadata$SampleID[match(colnames(perc_uniteObj), fullMetadata$SampleID)]
  
  # # Slidding window analysis on CpG covered in all fish which belong to parental DMS (N=123)
  # coveredDMSpar = methylKit::select(uniteCovALL_woSexAndUnknowChr, 
  #                                   which(paste(uniteCovALL_woSexAndUnknowChr$chr, uniteCovALL_woSexAndUnknowChr$start) %in% 
  #                                           paste(DMS15pc_G1_half$chr, DMS15pc_G1_half$start)))
  


