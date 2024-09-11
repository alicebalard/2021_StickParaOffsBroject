# Each script sources the previous script of the pipeline if needed
source("R06_GlobalMethylationProfile.R")
# Generates figure S1

## Load genome annotation
## NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
## The default option is to take -1000,+1000bp around the TSS and you can change that. 
## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream
annotBed12=readTranscriptFeatures("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12",
                                  remove.unusual = FALSE, up.flank = 1500, down.flank = 500)
## Change recursively the gene names to keep only ID
getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
for (i in 1:length(annotBed12)){
  annotBed12[[i]]$name <- getName(annotBed12[[i]]$name)
}

# Load the calculated DMS
load("../../dataOut/fitPQLseqG2_fit_G2_CC.TC.RData")
load("../../dataOut/fitPQLseqG2_fit_G2_CT.TT.RData")

fit_G2_CC.TC = fit_G2_CC.TC$fit ## to rm after clean raw output, rm fitpass
fit_G2_CT.TT = fit_G2_CT.TT$fit
##
fit_G2_CC.CT = head(fit_G2_CT.TT, 10000)## TO TEST ONLY RMMMM LATER
fit_G2_TC.TT = head(fit_G2_CC.TC, 10000) ## TO TEST ONLY RMMMM LATER

## Statistical setup
# PARENTAL effect: DMS found in either CC-TC or CT-TT comparisons
# OFFSPRING effect: DMS found in either CC-CT or TC-TT comparisons
# INTERACTION effects: DMS found in CC-CT which show a differential methylation (not necessarily significant)
# in the opposite direction in TC-TT, or inversely

## Define the threshold for pass
mypval = 0.05
mydif = 10

CC.TC_pass = fit_G2_CC.TC[fit_G2_CC.TC$converged==T & fit_G2_CC.TC$pvalue < mypval & 
                            (fit_G2_CC.TC$aveDiffMeth_ab>mydif | fit_G2_CC.TC$aveDiffMeth_ab< -mydif),]

CT.TT_pass = fit_G2_CT.TT[fit_G2_CT.TT$converged==T & fit_G2_CT.TT$pvalue < mypval & 
                            (fit_G2_CT.TT$aveDiffMeth_ab>mydif | fit_G2_CT.TT$aveDiffMeth_ab< -mydif),]

CC.CT_pass = fit_G2_CC.CT[fit_G2_CC.CT$converged==T & fit_G2_CC.CT$pvalue < mypval & 
                            (fit_G2_CC.CT$aveDiffMeth_ab>mydif | fit_G2_CC.CT$aveDiffMeth_ab< -mydif),]

TC.TT_pass = fit_G2_TC.TT[fit_G2_TC.TT$converged==T & fit_G2_TC.TT$pvalue < mypval & 
                            (fit_G2_TC.TT$aveDiffMeth_ab>mydif | fit_G2_TC.TT$aveDiffMeth_ab< -mydif),]

PATERNAL_DMS = rbind(CC.TC_pass, CT.TT_pass)
PATERNAL_DMS[duplicated(PATERNAL_DMS)] ## no duplicates

OFFSPRING_DMS = rbind(CC.CT_pass, TC.TT_pass)
OFFSPRING_DMS[duplicated(OFFSPRING_DMS)] ## tbc

INFECTION_INDUCED = row.names(dplyr::anti_join(OFFSPRING_DMS, PATERNAL_DMS))
INTERGENERATIONAL = row.names(dplyr::anti_join(PATERNAL_DMS,OFFSPRING_DMS))

## OVERLAP=DMS found in PATERNAL_DMS & OFFSPRING_DMS
# case 1 -> INTERACTION effects DMS found in CC-CT which show a differential methylation
# in the opposite direction in TC-TT, or inversely (reaction norm are inversed)
# case 2 -> ADDITIVE effect: no slope inversion

OVERLAP = dplyr::intersect(PATERNAL_DMS, OFFSPRING_DMS)

# INTERACTION = DMS found in comparison 1) “G1control-G2control (CC) vs G1control-G2infected (TC)” 
# and showed a mean differential methylation in the opposite direction in comparison 
# 2) “G1infected-G2control (CT) vs G1infected-G2infected (TT)”, or inversely

A=fit_G2_CC.TC[row.names(fit_G2_CC.TC) %in% row.names(OVERLAP),"aveDiffMeth_ab"]
B=fit_G2_CT.TT[row.names(fit_G2_CT.TT) %in% row.names(OVERLAP),"aveDiffMeth_ab"]
C=row.names(fit_G2_CC.TC[row.names(fit_G2_CC.TC) %in% row.names(OVERLAP), ])

INTERACTION = C[sign(A) != sign(B)]
ADDITIVE = C[sign(A) == sign(B)]

rm(A,B,C)

dataFinal = data.frame(NbrDMS = c(length(INTERGENERATIONAL),
                                  length(INFECTION_INDUCED),
                                  length(ADDITIVE),
                                  length(INTERACTION)),
                       DMSgroup = as.factor(c("Intergenerational", "Infection-induced",  "Additive", "Interaction")))

# # Plot a Venn diagram
# ggVennDiagram(list("Paternal effect" = row.names(PATERNAL_DMS), "Offspring effect" = row.names(OFFSPRING_DMS),
#                    "InteractionEffects" = INTERACTION),
#               label_alpha = 0) + scale_color_manual(values = c(1,1,1))+
#   scale_fill_gradient(low="white",high = "yellow") + theme(legend.position = "none")
# 
# # Save:
# # pdf(file = "../../dataOut/DMS3groupsVenn.pdf", width = 7, height = 6)
# ggVennDiagram(list("Paternal effect" = row.names(PATERNAL_DMS), 
#                    "Offspring effect" = row.names(OFFSPRING_DMS),
#                    "INTER" = INTERACTION,
#                    "ADD" = ADDITIVE),
#               label_alpha = 0) + scale_color_manual(values = c(1,1,1))+
#   scale_fill_gradient(low="white",high = "yellow") + theme(legend.position = "none")
# # dev.off()

#######################
## Where are these DMS?
DMSvec=unique(c(INTERGENERATIONAL, INFECTION_INDUCED, ADDITIVE, INTERACTION))

getFeature <- function(DMSvec){
  # Change the DMS vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=paste(sapply(strsplit(DMSvec, "_"), `[`, 1), sapply(strsplit(DMSvec, "_"), `[`, 2), sep = "_"), 
                                                  start=sapply(strsplit(DMSvec, "_"), `[`, 3),
                                                  end=sapply(strsplit(DMSvec, "_"), `[`, 3),
                                                  DMS=DMSvec), keep.extra.columns = T)
  annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = annotBed12)
}

A=getFeature(DMSvec = DMSvec)
print(paste0("Positions of the ", length(DMSvec)," DMS:"))
print(A)

print(paste0("Positions of the ", length(INTERGENERATIONAL)," intergenerational DMS:"))
print(getFeature(DMSvec = INTERGENERATIONAL)@precedence)

print(paste0("Positions of the ", length(INFECTION_INDUCED)," infection-induced DMS:"))
print(getFeature(DMSvec = INFECTION_INDUCED)@precedence)

print(paste0("Positions of the ", length(ADDITIVE)," additive DMS:"))
print(getFeature(DMSvec = ADDITIVE)@precedence)

print(paste0("Positions of the ", length(INTERACTION)," interaction DMS:"))
print(getFeature(DMSvec = INTERACTION)@precedence)

#######################
## Are the positions of DMS on features random? Comparison with sequenced CpGs which are not DMS
allCpG = paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, 
               uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$start, sep = "_")

A=getFeature(DMSvec = DMSvec)
AallCpG= getFeature(DMSvec = allCpG)

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

chisq.test(ChiTable1[c("DMS", "allCpG")])
## TBCorrected when all run

#######################
## Are the positions of DMS on chromosomes random? Comparison with sequenced CpGs which are not DMS
a=table(paste(sapply(strsplit(DMSvec, "_"), `[`, 1), sapply(strsplit(DMSvec, "_"), `[`, 2), sep = "_"))
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
## TBC

df %>% arrange(percent) 
# range from 0.13% (XV) to 0.39% of CpG beign DMS (XVIII) TBC++++

#### Plot sup fig 1
donutDF = ChiTable1 %>% dplyr::mutate(percDMS=DMS/sum(DMS)*100, percCpG=allCpG/sum(allCpG)*100) %>%
  dplyr::select("feature", "percDMS", "percCpG") %>% melt

P1 <- ggplot(donutDF, aes(x = variable, y = value, fill = feature)) +
  geom_col() +  scale_x_discrete(limits = c(" ", "percCpG","percDMS")) +
  coord_polar("y")+
  theme_void()+
  scale_fill_manual(values = c("#1F449C", "#EEBAB4", "#A8B6CC","#E57A77"))+
  annotate("text", x = 3.8, y = 0, label = "DMS")+
  annotate("text", x = 1, y = 0, label = "all CpGs")

P2 <- ggplot(df, aes(x=chromosome, y=percent))+
  geom_bar(stat = "identity", col= "black", fill="black")+
  ylab("Percentage of DMS among CpG sequenced")+
  scale_y_continuous(labels=scales::percent)

pdf(file = "../../dataOut/suppl1.pdf", width = 8, height = 6)
cowplot::plot_grid(P1, P2, nrow = 2, labels = c("A", "B"))
dev.off()
