# Each script sources the previous script of the pipeline if needed
source("R06_GlobalMethylationProfile.R")

message("R07 starting...\n")
# Generates:
## figure S1 in dataOut/fig/FigS1_featuresCpGandDMS.pdf

# Load the calculated DMS
## G1 diff
load("../../gitignore/bigdata/PQLseq/fitPQLseqG2_fit_G2_CC.TC.RData")
load("../../gitignore/bigdata/PQLseq/fitPQLseqG2_fit_G2_CT.TT.RData")
## G2 diff
load("../../gitignore/bigdata/PQLseq/fitPQLseqG2_fit_G2_CC.CT.RData")
load("../../gitignore/bigdata/PQLseq/fitPQLseqG2_fit_G2_TC.TT.RData")

## Statistical setup
# PARENTAL effect: DMS found in either CC-TC or CT-TT comparisons
# OFFSPRING effect: DMS found in either CC-CT or TC-TT comparisons
# INTERACTION effects: DMS found in CC-CT which show a differential methylation (not necessarily significant)
# in the opposite direction in TC-TT, or inversely

################################
## Define the threshold for pass
myqval = 0.05
mydif = 5

message(paste0("We define a threshold for DMS: q-value > ", myqval, 
               " , difference of methylation > ", mydif, "%"))

makeVolcano = function(fit){
  ggplot()+
    geom_point(data = fit[fit$qvalue < 0.1,],
               aes(x = aveDiffMeth_ab, y = qvalue), alpha = .3) +
    scale_y_continuous(trans = "log1p") + # log(x+1)
    geom_rect(aes(xmin=-Inf, xmax=-5, ymin=0, ymax=0.05), fill = "red", alpha = .3)+
    geom_rect(aes(xmin=5, xmax=Inf, ymin=0, ymax=0.05), fill = "red", alpha = .3)
}

volcanoPlotThreshold <- cowplot::plot_grid(makeVolcano(fit_G2_CC.CT), makeVolcano(fit_G2_CC.TC),
                                           makeVolcano(fit_G2_CT.TT), makeVolcano(fit_G2_TC.TT), 
                                           nrow = 2, labels = c("CC.CT", "CC.TC", "CT.TT", "TC.TT"))
volcanoPlotThreshold

## Select based on threshold
CC.TC_pass = fit_G2_CC.TC[fit_G2_CC.TC$converged==T & fit_G2_CC.TC$qvalue <= myqval & 
                            (fit_G2_CC.TC$aveDiffMeth_ab>mydif | fit_G2_CC.TC$aveDiffMeth_ab< -mydif),]

CT.TT_pass = fit_G2_CT.TT[fit_G2_CT.TT$converged==T & fit_G2_CT.TT$qvalue <= myqval & 
                            (fit_G2_CT.TT$aveDiffMeth_ab>mydif | fit_G2_CT.TT$aveDiffMeth_ab< -mydif),]

CC.CT_pass = fit_G2_CC.CT[fit_G2_CC.CT$converged==T & fit_G2_CC.CT$qvalue <= myqval & 
                            (fit_G2_CC.CT$aveDiffMeth_ab>mydif | fit_G2_CC.CT$aveDiffMeth_ab< -mydif),]

TC.TT_pass = fit_G2_TC.TT[fit_G2_TC.TT$converged==T & fit_G2_TC.TT$qvalue <= myqval & 
                            (fit_G2_TC.TT$aveDiffMeth_ab>mydif | fit_G2_TC.TT$aveDiffMeth_ab< -mydif),]

PATERNAL_DMS = rbind(CC.TC_pass, CT.TT_pass) %>% dplyr::select(chrom,start,end, pos)
row.names(PATERNAL_DMS) = NULL
PATERNAL_DMS[duplicated(PATERNAL_DMS)] ## no duplicates
message(paste0("Nbr of paternal DMS: ", nrow(PATERNAL_DMS))) #314

OFFSPRING_DMS = rbind(CC.CT_pass, TC.TT_pass) %>% dplyr::select(chrom,start,end, pos)
row.names(OFFSPRING_DMS) = NULL
OFFSPRING_DMS[duplicated(OFFSPRING_DMS)] ## no duplicates
message(paste0("Nbr of offspring DMS: ", nrow(OFFSPRING_DMS))) #400

PATERNAL_DMS$effect = "PATERNAL"
OFFSPRING_DMS$effect = "OFFSPRING"

intersect(PATERNAL_DMS$pos, OFFSPRING_DMS$pos) # 7 intersection

###########################
## Create our 4 categories:
outersect <- function(x, y) {
  sort(c(setdiff(x, y), setdiff(y, x)))
}

## 1. INFECTION_INDUCED: methylation change only due to offspring treatment
INFECTION_INDUCED = OFFSPRING_DMS[OFFSPRING_DMS$pos %in% outersect(OFFSPRING_DMS$pos, PATERNAL_DMS$pos),]
INFECTION_INDUCED$effect = "INFECTION_INDUCED"

## 2. INTERGENERATIONAL: methylation change only due to paternal treatment
INTERGENERATIONAL = PATERNAL_DMS[PATERNAL_DMS$pos %in% outersect(OFFSPRING_DMS$pos, PATERNAL_DMS$pos),]
INTERGENERATIONAL$effect = "INTERGENERATIONAL"

## OVERLAP=DMS found in PATERNAL_DMS & OFFSPRING_DMS
OVERLAP = PATERNAL_DMS[PATERNAL_DMS$pos %in% intersect(PATERNAL_DMS$pos, OFFSPRING_DMS$pos),]
OVERLAP$effect = "OVERLAP"

## 3.INTERACTION: interaction effects in the overlap, i.e. DMS found in CC-CT which show a differential methylation
# in the opposite direction in TC-TT, or inversely (reaction norm are inverse)

## 4.ADDITIVE: additive effect in the overlap, no slope inversion

# INTERACTION = DMS found in comparison 1) “G1control-G2control (CC) vs G1control-G2infected (TC)” 
# and showed a mean differential methylation in the opposite direction in comparison 
# 2) “G1infected-G2control (CT) vs G1infected-G2infected (TT)”, or inversely

A=fit_G2_CC.TC[fit_G2_CC.TC$pos %in% OVERLAP$pos,]
B=fit_G2_CT.TT[fit_G2_CT.TT$pos %in% OVERLAP$pos,]

INTERACTION = OVERLAP[sign(A$aveDiffMeth_ab) != sign(B$aveDiffMeth_ab),]
INTERACTION$effect = "INTERACTION"

ADDITIVE = OVERLAP[sign(A$aveDiffMeth_ab) == sign(B$aveDiffMeth_ab),]
ADDITIVE$effect = "ADDITIVE"

rm(A,B)

## DF of all our effects
EffectsDF <- rbind(INTERGENERATIONAL, INFECTION_INDUCED, INTERACTION, ADDITIVE)
table(duplicated(EffectsDF$pos)) # all should be unique
## order effects factor
EffectsDF$effect <- factor(EffectsDF$effect, levels = c("INFECTION_INDUCED", "INTERGENERATIONAL", "ADDITIVE", "INTERACTION" ))

table(EffectsDF$effect)
# INFECTION_INDUCED INTERGENERATIONAL          ADDITIVE       INTERACTION 
# 393               307                 1                 6 

message(paste0("We identified ", sum(table(EffectsDF$effect)), " DMS out of ",
               nrow(fit_G2_CC.CT), " CpG sequenced (",
               round(sum(table(EffectsDF$effect))/nrow(fit_G2_CC.CT)*100,2), "%) including ",
               table(EffectsDF$effect)[1], " infection induced (", 
               round(table(EffectsDF$effect)[1]/sum(table(EffectsDF$effect))*100,2), "%), ",
               table(EffectsDF$effect)[2], " intergenerational (", 
               round(table(EffectsDF$effect)[2]/sum(table(EffectsDF$effect))*100,2), "%), ",
               table(EffectsDF$effect)[3], " additive (", 
               round(table(EffectsDF$effect)[3]/sum(table(EffectsDF$effect))*100,2), "%), ",
               table(EffectsDF$effect)[4], " interaction (", 
               round(table(EffectsDF$effect)[4]/sum(table(EffectsDF$effect))*100,2), "%)"))
               
## Add chromosome length and relative position for future plots
EffectsDF <- merge(EffectsDF, GYgynogff[c("chrom", "length", "gstart", "gend")])

###################################
## On which features are these DMS?

getFeature <- function(DMSvec){
  # Change the DMS vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=paste(sapply(strsplit(DMSvec, "_"), `[`, 1), sapply(strsplit(DMSvec, "_"), `[`, 2), sep = "_"), 
                                                  start=sapply(strsplit(DMSvec, "_"), `[`, 3),
                                                  end=sapply(strsplit(DMSvec, "_"), `[`, 3),
                                                  DMS=DMSvec), keep.extra.columns = T)
  annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = annotBed12)
}

print(paste0("Positions of the ", length(EffectsDF$pos)," DMS:"))
print(getFeature(DMSvec = EffectsDF$pos))

print(paste0("Positions of the ", length(fit_G2_CC.CT$pos)," CpG sequenced:"))
allCpG <- fit_G2_CC.CT$pos
getFeatallCpG <- getFeature(DMSvec = allCpG)
print(getFeatallCpG)

## Are the positions of DMS on features random? 
## Comparison with sequenced CpGs which are not DMS
A=getFeature(DMSvec = EffectsDF$pos)
AallCpG=getFeatallCpG

ChiTable1 = merge((A@members %>% data.frame() %>% mutate(feature=ifelse(prom==1, "promoter", 
                                                                        ifelse(exon==1, "exon",
                                                                               ifelse(intron==1, "intron", "intergenic")))) %>% 
                     dplyr::select(feature) %>% 
                     table %>% melt %>% dplyr::rename(DMS=value)),
                  (AallCpG@members %>% data.frame() %>% mutate(feature=ifelse(prom==1, "promoter", 
                                                                              ifelse(exon==1, "exon",
                                                                                     ifelse(intron==1, "intron", "intergenic")))) %>% 
                     dplyr::select(feature) %>% 
                     table %>% melt %>% dplyr::rename(allCpG=value)))
row.names(ChiTable1) <- ChiTable1$feature
ChiTable1 <-ChiTable1[!names(ChiTable1) %in% "feature"]

chisq.test(ChiTable1)
## Corrected when annotation finished
# Pearson's Chi-squared test
# data:  ChiTable1
# X-squared = 53.747, df = 3, p-value = 1.271e-11

##run a post hoc analysis for Pearson’s Chi-squared Test for Count Data
## ref T. Mark Beasley & Randall E. Schumacker (1995) Multiple Regression Approach to Analyzing Contingency Tables: Post Hoc and Planned Comparison Procedures, The Journal of Experimental Education, 64:1, 79-93, DOI: 10.1080/00220973.1995.9943797

ChiTable1 %>% mutate(DMSprop=DMS/sum(DMS),allCpGprop=allCpG/sum(allCpG))
chisq.posthoc.test(ChiTable1, method = "bonferroni")
## Exon and introns are the same, but proportionally less in promoters and more in intergenic regions

#######################
## Are the positions of DMS on chromosomes random? Comparison with sequenced CpGs which are not DMS
a=table(paste(sapply(strsplit(EffectsDF$pos, "_"), `[`, 1), sapply(strsplit(EffectsDF$pos, "_"), `[`, 2), sep = "_"))
b=table(paste(sapply(strsplit(allCpG, "_"), `[`, 1), sapply(strsplit(allCpG, "_"), `[`, 2), sep = "_"))

a=as.data.frame(a)
names(a)=c("chr", "nbrDMS")
b=as.data.frame(b)
names(b)=c("chr", "nbrCpG")

df=merge(a,b)

# Number of positions sequenced on each chromosome
df2 = table(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr) %>% data.frame() %>% 
  mutate(chromosome=str_remove(Var1,"Gy_chr"), nCpG=Freq) 
names(df2)[1]="chr"

df=merge(df, df2) %>% dplyr::select(c(chromosome, nbrDMS, nbrCpG))%>% 
  mutate(percent=nbrDMS/(nbrDMS+nbrCpG))

chisq.test(df[c("nbrDMS", "nbrCpG")])
# Pearson's Chi-squared test
# data:  df[c("nbrDMS", "nbrCpG")]
# X-squared = 38.427, df = 19, p-value = 0.005235

df %>% arrange(percent) %>% mutate(perc = round(percent*100,2))
# from 0.05% of the CpG sequenced (chromosome XIV) to 0.11% (chromosome V)

#######################
## Same question but split by effect: 
## Are the positions of DMS for both IG & INT on chromosomes random? Comparison with sequenced CpGs which are not DMS
a1=table(paste(sapply(strsplit(EffectsDF$pos[EffectsDF$effect %in% "INFECTION_INDUCED"], "_"), `[`, 1),
              sapply(strsplit(EffectsDF$pos[EffectsDF$effect %in% "INFECTION_INDUCED"], "_"), `[`, 2), 
              sep = "_"))

a2=table(paste(sapply(strsplit(EffectsDF$pos[EffectsDF$effect %in% "INTERGENERATIONAL"], "_"), `[`, 1),
               sapply(strsplit(EffectsDF$pos[EffectsDF$effect %in% "INTERGENERATIONAL"], "_"), `[`, 2), 
               sep = "_"))

b=table(paste(sapply(strsplit(allCpG, "_"), `[`, 1), sapply(strsplit(allCpG, "_"), `[`, 2), sep = "_"))

a1=as.data.frame(a1)
names(a1)=c("chr", "nbrDMS_infind")

a2=as.data.frame(a2)
names(a2)=c("chr", "nbrDMS_interg")

b=as.data.frame(b)
names(b)=c("chr", "nbrCpG")

df=merge(merge(a1,a2), b)

## distribution of DMS along chromosomes?
M <- as.table(cbind(df$nbrDMS_infind, df$nbrDMS_interg, df$nbrCpG))
dimnames(M) <- list(chr = df$chr,
                    cat = c("DMS_infind", "DMS_interg", "CpG"))
chisq.test(M)
# Pearson's Chi-squared test
# 
# data:  M
# X-squared = 74.756, df = 38, p-value = 0.0003433
chisq.posthoc.test(M, method = "bonferroni")

## significantly higher number of infection-induced DMS at chromosome XVIII (P=0.018),
## and significantly higher number of intergenerational DMS at chromosome XX (P=8e-04*)

# We compared the distribution of DMS on the chromosomes using a chisquare test with number 
# of infection-induced DMS, number of intergenerational DMS and total number of CpG 
# sequenced as columns and chromosome as rows, followed by ebbertd/chisq.posthoc.test
# with Bonferroni correction.

df = df %>% mutate(percDMS_interg=nbrDMS_interg/nbrCpG, percDMS_infind=nbrDMS_infind/nbrCpG)

df2 <- df %>%
  pivot_longer(
    cols = c(percDMS_interg, percDMS_infind),
    names_to = "category",
    values_to = "value"
  ) %>% data.frame()

df2$chromosome <- gsub("Gy_chr", "", df2$chr)

###################
#### Plot sup fig 1
donutDF = ChiTable1 %>% dplyr::mutate(feature = row.names(ChiTable1), percDMS=DMS/sum(DMS)*100, percCpG=allCpG/sum(allCpG)*100) %>%
  dplyr::select("feature", "percDMS", "percCpG") %>% melt

P1 <- ggplot(donutDF, aes(x = variable, y = value, fill = feature)) +
  geom_col() +  scale_x_discrete(limits = c(" ", "percCpG","percDMS")) +
  coord_polar("y")+
  theme_void()+
  scale_fill_manual(values = c("#1F449C", "#EEBAB4", "#A8B6CC","#E57A77"))+
  annotate("text", x = 3.8, y = 0, label = "DMS")+
  annotate("text", x = 1, y = 0, label = "all CpGs")

col = colEffects[c("infection induced", "intergenerational")]
names(col) = NULL
P2 <- ggplot(df2, aes(x=chromosome, y=value, group = category, fill = category))+
  geom_bar(stat = "identity", col= "black", position = "dodge")+
  ylab("Percentage of DMS among CpG sequenced")+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values = col)

message("Save figure S1 in dataOut/fig/FigS1_featuresCpGandDMS.pdf")
pdf(file = "../../dataOut/fig/FigS1_featuresCpGandDMS.pdf", width = 8, height = 6)
cowplot::plot_grid(P1, P2, nrow = 2, labels = c("A", "B"))
dev.off()

message("R07 done.\n")
