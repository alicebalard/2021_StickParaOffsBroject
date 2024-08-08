# # Each script sources the previous script of the pipeline if needed
# source("R09_CorrelationPCAmethPhenotype.R")
# 
# # Differential methylations "Top down": DM in parents -\> how do they look in offspring?
# 
# NB: was finally not exploited, code has change, if we go back to that they would be some cleaning to do
#
# ### TREATMENT EFFECT: Differential methylation between control and infected, within each parental treatment
# 
# ## Features Annotation (use package genomation v1.24.0)
# ## NB Promoters are defined by options at genomation::readTranscriptFeatures function.
# ## The default option is to take -1000,+1000bp around the TSS and you can change that.
# ## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream
# 
# ## Parents comparison:
# diffAnn_PAR = annotateWithGeneParts(as(DMSlist$DMS_15pc_G1_C_T,"GRanges"),annotBed12)
# ## Offspring from control parents comparison:
# diffAnn_G2_controlG1 = annotateWithGeneParts(as(DMSlist$DMS_15pc_CC_CT,"GRanges"),annotBed12)
# ## Offspring from infected parents comparison:
# diffAnn_G2_infectedG1 = annotateWithGeneParts(as(DMSlist$DMS_15pc_TC_TT,"GRanges"),annotBed12)
# ```
# 
# ```{r, fig.height = 6, fig.width = 10}
# par(mfrow=c(1,3))
# par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
# ## Parents comparison:
# diffAnn_PAR = annotateWithGeneParts(as(DMSlist$DMS_15pc_G1_C_T,"GRanges"),annotBed12)
# diffAnn_PAR
# genomation::plotTargetAnnotation(diffAnn_PAR,precedence=TRUE, main="DMS G1", cex.legend = 1, border="white")
# 
# ## Offspring from control parents comparison:
# diffAnn_G2_controlG1 = annotateWithGeneParts(as(DMSlist$DMS_15pc_CC_CT,"GRanges"),annotBed12)
# diffAnn_G2_controlG1
# genomation::plotTargetAnnotation(diffAnn_G2_controlG1,precedence=TRUE, main="DMS G2-G1c", cex.legend = 1, border="white")
# ## Offspring from infected parents comparison:
# diffAnn_G2_infectedG1 = annotateWithGeneParts(as(DMSlist$DMS_15pc_TC_TT,"GRanges"),annotBed12)
# diffAnn_G2_infectedG1
# genomation::plotTargetAnnotation(diffAnn_G2_infectedG1,precedence=TRUE, main="DMS G2-G1i", cex.legend = 1, border="white")
# par(mfrow=c(1,1))
# ```
# 
# ```{r}
# ## Run the function to get DMS info
# DMS_info_G1 <- myDMSinfo(DMSlist$DMS_15pc_G1_C_T, uniteCov6_G1_woSexAndUnknowChrOVERLAP)
# DMS_info_G2_G1c_final <- myDMSinfo(DMSlist$DMS_15pc_CC_CT, uniteCov14_G2_woSexAndUnknowChrOVERLAP)
# DMS_info_G2_G1i_final <- myDMSinfo(DMSlist$DMS_15pc_TC_TT,uniteCov14_G2_woSexAndUnknowChrOVERLAP)
# ```
# 
# NB Kostas' results: "We found a total of 1,973 CpG sites out of
# 1,172,887 CpGs (0.17%) across the genome that showed at least 15%
# differential fractional methylation (differentially methylated site
# [DMS]; q \< 0.01) between infected and uninfected fish"
# 
# Here: we obtain out a total of
# `r nrow(uniteCov14_G2_woSexAndUnknowChrOVERLAP)` CpG sites:
# `r length(DMS_info_G1$DMS)` (`r round(DMS_info_G1$percentDMS, 2)`%)
# showed at least 15% differential fractional methylation (differentially
# methylated site [DMS]; q \< 0.01) between infected and uninfected fish
# in the parental group; `r length(DMS_info_G2_G1c_final$DMS)`
# (`r round(DMS_info_G2_G1c_final$percentDMS, 2)`%) in the offspring from
# **control parents** comparison; `r length(DMS_info_G2_G1i_final$DMS)`
# (`r round(DMS_info_G2_G1i_final$percentDMS, 2)`%) in the offspring from
# **infected parents** comparison.
# 
# #### Intersections of DMS: Venn diagrams
# 
# ```{r}
# ## Chi2 test: are the number of DMS from G2-G1C and G2-G1I overlapping with DMSpar statistically different?
# 
# A=length(intersect(DMS_info_G1$DMS,DMS_info_G2_G1c_final$DMS))
# B=length(DMS_info_G2_G1c_final$DMS)
# C=length(intersect(DMS_info_G1$DMS,DMS_info_G2_G1i_final$DMS))
# D=length(DMS_info_G2_G1i_final$DMS)
# 
# Observed=matrix(c(A, B-A, C, D-C),nrow=2)
# Observed
# 
# chisq.test(Observed)
# ## not statistically different
# ```
# 
# ```{r, fig.width = 12, fig.height = 15}
# ## output Venn diagrams
# allVenn <- ggVennDiagram(list("DMS G1" = DMS_info_G1$DMS, "DMS G2-c" = DMS_info_G2_G1c_final$DMS, "DMS G2-i" = DMS_info_G2_G1i_final$DMS), label_alpha = 0) +
#   scale_fill_gradient(low="white",high = "red")
# hypoVenn <- ggVennDiagram(list("DMS G1\nhypo" = DMS_info_G1$DMS[DMS_info_G1$direction %in% "hypo"],
#                                "DMS G2-c\nhypo" = DMS_info_G2_G1c_final$DMS[DMS_info_G2_G1c_final$direction %in% "hypo"],
#                                "DMS G2-i\nhypo" = DMS_info_G2_G1i_final$DMS[DMS_info_G2_G1i_final$direction %in% "hypo"]), label_alpha = 0) +
#   scale_fill_gradient(low="white",high = "red")
# hyperVenn <- ggVennDiagram(list("DMS G1\nhyper" = DMS_info_G1$DMS[DMS_info_G1$direction %in% "hyper"],
#                                 "DMS G2-c\nhyper" = DMS_info_G2_G1c_final$DMS[DMS_info_G2_G1c_final$direction %in% "hyper"],
#                                 "DMS G2-i\nhyper" = DMS_info_G2_G1i_final$DMS[DMS_info_G2_G1i_final$direction %in% "hyper"]), label_alpha = 0) +
#   scale_fill_gradient(low="white",high = "red")
# 
# ggarrange(allVenn,
#           ggarrange(hypoVenn, hyperVenn, ncol = 2, legend = "none"),
#           nrow = 2, widths = c(.5,1))
# ```
# 
# ##### Venn with annotated features
# 
# ```{r}
# runHyperHypoAnnot <- function(){
#   par(mfrow=c(2,3))
#   par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
#   ####### HYPO
#   ## Parents comparison:
#   A = annotateWithGeneParts(
#     as(DMSlist$DMS_15pc_G1_C_T[DMS_info_G1$direction %in% "hypo",],"GRanges"),annotBed12)
#   genomation::plotTargetAnnotation(A,precedence=TRUE, main="DMS G1\nhypo",
#                                    cex.legend = .4, border="white")
#   ## Offspring from control parents comparison:
#   B = annotateWithGeneParts(
#     as(DMSlist$DMS_15pc_CC_CT[DMS_info_G2_G1c_final$direction %in% "hypo",],"GRanges"),annotBed12)
#   genomation::plotTargetAnnotation(B,precedence=TRUE, main="DMS G2-G1c\nhypo",
#                                    cex.legend = .4, border="white")
#   ## Offspring from infected parents comparison:
#   C = annotateWithGeneParts(
#     as(DMSlist$DMS_15pc_TC_TT[DMS_info_G2_G1i_final$direction %in% "hypo",],"GRanges"),annotBed12)
#   genomation::plotTargetAnnotation(C,precedence=TRUE, main="DMS G2-G1i\nhypo",
#                                    cex.legend = .4, border="white")
#   ####### HYPER
#   ## Parents comparison:
#   D = annotateWithGeneParts(
#     as(DMSlist$DMS_15pc_G1_C_T[DMS_info_G1$direction %in% "hyper",],"GRanges"),annotBed12)
#   genomation::plotTargetAnnotation(D,precedence=TRUE, main="DMS G1\nhyper",
#                                    cex.legend = .4, border="white")
#   ## Offspring from control parents comparison:
#   E = annotateWithGeneParts(
#     as(DMSlist$DMS_15pc_CC_CT[DMS_info_G2_G1c_final$direction %in% "hyper",],"GRanges"),annotBed12)
#   genomation::plotTargetAnnotation(E,precedence=TRUE, main="DMS G2-G1c\nhyper",
#                                    cex.legend = .4, border="white")
#   ## Offspring from infected parents comparison:
#   f = annotateWithGeneParts(
#     as(DMSlist$DMS_15pc_TC_TT[DMS_info_G2_G1i_final$direction %in% "hyper",],"GRanges"),annotBed12)
#   genomation::plotTargetAnnotation(f,precedence=TRUE, main="DMS G2-G1i\nhyper",
#                                    cex.legend = .4, border="white")
#   par(mfrow=c(1,1))
#   return(list(G1hypo=A, G2G1chypo=B, G2G1ihypo=C, G1hyper=D, G2G1chyper=E, G2G1ihyper=f))
# }
# 
# myannot=runHyperHypoAnnot()
# ```
# 
# ```{r}
# ############################################################
# ## Venn diagram of overlapping features by their annotation:
# table(rowSums(as.data.frame(myannot$G1hypo@members))) # NB: some positions are labelled with several features!
# ## as in MBE 2021: "giving precedence to the following order promoters, exons,
# ## introns, and intergenic regions when features overlapped"
# 
# myAnnotateDMS <- function(DMS, annot){
#   ## sanity check
#   if (nrow(DMS) != nrow(annot)){"STOP error in arguments"}
#   DMS$pos <- paste(DMS$chr, DMS$start, DMS$end)
#   ## NB as in MBE 2021: "giving precedence to the following order promoters, exons,
#   ## introns, and intergenic regions when features overlapped"
#   DMS$feature <- NA
#   ## 1. promoters
#   DMS$feature[which(annot$prom == 1)] = "promoter"
#   ## 2. exons
#   DMS$feature[which(annot$exon == 1 & annot$prom ==0)] = "exon"
#   ## 3. intron
#   DMS$feature[which(annot$intro == 1 & annot$exon == 0 & annot$prom ==0)] = "intron"
#   ## 4. intergenic regions
#   DMS$feature[which(annot$intro == 0 & annot$exon == 0 & annot$prom ==0)] = "intergenic"
#   return(DMS)
# }
# 
# DMSlist$DMS_15pc_G1_C_T = myAnnotateDMS(DMSlist$DMS_15pc_G1_C_T, as.data.frame(diffAnn_PAR@members))
# DMSlist$DMS_15pc_G1_C_T_HYPO = myAnnotateDMS(DMSlist$DMS_15pc_G1_C_T[DMS_info_G1$direction %in% "hypo",],
#                                              as.data.frame(myannot$G1hypo@members))
# DMSlist$DMS_15pc_G1_C_T_HYPER = myAnnotateDMS(DMSlist$DMS_15pc_G1_C_T[DMS_info_G1$direction %in% "hyper",],
#                                               as.data.frame(myannot$G1hyper@members))
# 
# DMSlist$DMS_15pc_CC_CT = myAnnotateDMS(DMSlist$DMS_15pc_CC_CT, as.data.frame(diffAnn_G2_controlG1@members))
# DMSlist$DMS_15pc_CC_CT_HYPO = myAnnotateDMS(DMSlist$DMS_15pc_CC_CT[DMS_info_G2_G1c_final$direction %in% "hypo",],
#                                             as.data.frame(myannot$G2G1chypo@members))
# DMSlist$DMS_15pc_CC_CT_HYPER = myAnnotateDMS(DMSlist$DMS_15pc_CC_CT[DMS_info_G2_G1c_final$direction %in% "hyper",],
#                                              as.data.frame(myannot$G2G1chyper@members))
# 
# DMSlist$DMS_15pc_TC_TT = myAnnotateDMS(DMSlist$DMS_15pc_TC_TT, as.data.frame(diffAnn_G2_infectedG1@members))
# DMSlist$DMS_15pc_TC_TT_HYPO = myAnnotateDMS(DMSlist$DMS_15pc_TC_TT[DMS_info_G2_G1i_final$direction %in% "hypo",],
#                                             as.data.frame(myannot$G2G1ihypo@members))
# DMSlist$DMS_15pc_TC_TT_HYPER = myAnnotateDMS(DMSlist$DMS_15pc_TC_TT[DMS_info_G2_G1i_final$direction %in% "hyper",],
#                                              as.data.frame(myannot$G2G1ihyper@members))
# 
# ## Make Venn diagram for each feature
# getFeatureDFHYPO <- function(myfeat){
#   a = DMSlist$DMS_15pc_G1_C_T_HYPO$pos[DMSlist$DMS_15pc_G1_C_T_HYPO$feature %in% myfeat]
#   b = DMSlist$DMS_15pc_CC_CT_HYPO$pos[DMSlist$DMS_15pc_CC_CT_HYPO$feature %in% myfeat]
#   c = DMSlist$DMS_15pc_TC_TT_HYPO$pos[DMSlist$DMS_15pc_TC_TT_HYPO$feature %in% myfeat]
#   return(list(a=a,b=b,c=c))
# }
# 
# getFeatureDFHYPER <- function(myfeat){
#   a = DMSlist$DMS_15pc_G1_C_T_HYPER$pos[DMSlist$DMS_15pc_G1_C_T_HYPER$feature %in% myfeat]
#   b = DMSlist$DMS_15pc_CC_CT_HYPER$pos[DMSlist$DMS_15pc_CC_CT_HYPER$feature %in% myfeat]
#   c = DMSlist$DMS_15pc_TC_TT_HYPER$pos[DMSlist$DMS_15pc_TC_TT_HYPER$feature %in% myfeat]
#   return(list(a=a,b=b,c=c))
# }
# 
# getVenn <- function(feat, direction){
#   if (direction == "hypo"){
#     ggVennDiagram(list(A = getFeatureDFHYPO(feat)[["a"]],
#                        B = getFeatureDFHYPO(feat)[["b"]],
#                        C = getFeatureDFHYPO(feat)[["c"]]), label_alpha = 0,
#                   category.names = c(paste0("DMS G1\nhypo\n", feat), paste0("DMS G2-c\nhypo\n", feat), paste0("DMS G2-i\nhypo\n", feat))) +
#       scale_fill_gradient(low="white",high = "red")
#   } else if (direction == "hyper"){
#     ggVennDiagram(list(A = getFeatureDFHYPER(feat)[["a"]],
#                        B = getFeatureDFHYPER(feat)[["b"]],
#                        C = getFeatureDFHYPER(feat)[["c"]]), label_alpha = 0,
#                   category.names = c(paste0("DMS G1\nhyper\n", feat), paste0("DMS G2-c\nhyper\n", feat), paste0("DMS G2-i\nhyper\n", feat))) +
#       scale_fill_gradient(low="white",high = "red")
#   }
# }
# ```
# 
# ```{r, fig.width = 12, fig.height = 15}
# ggarrange(getVenn("promoter", "hypo"), getVenn("exon", "hypo"),
#           getVenn("intron", "hypo"), getVenn("intergenic", "hypo"),
#           nrow = 2, ncol = 2)
# ```
# 
# ```{r, fig.width = 12, fig.height = 15}
# ggarrange(getVenn("promoter", "hyper"), getVenn("exon", "hyper"),
#           getVenn("intron", "hyper"), getVenn("intergenic", "hyper"),
#           nrow = 2, ncol = 2)
# ```
# 
# #### Manhattan plot of DMS
# 
# ```{r, fig.width= 10, fig.height=3}
# ## Parents trt-ctrl
# # load annotation
# annot_PAR <- as.data.frame(diffAnn_PAR@members)
# makeManhattanPlots(DMSfile = DMSlist$DMS_15pc_G1_C_T, annotFile = annot_PAR, GYgynogff = GYgynogff,
#                    mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")
# 
# ## G2-G1c trt-ctrl
# # load annotation
# annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)
# makeManhattanPlots(DMSfile = DMSlist$DMS_15pc_CC_CT, annotFile = annot_G2_G1c, GYgynogff = GYgynogff,
#                    mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")
# 
# ## G2-G1i trt-ctrl
# # load annotation
# annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)
# makeManhattanPlots(DMSfile = DMSlist$DMS_15pc_TC_TT, annotFile = annot_G2_G1i, GYgynogff = GYgynogff,
#                    mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")
# 
# ## Outliers in Manhattan plot: 15% diff + 2SD
# outliers_G1_final <- which(abs(DMSlist$DMS_15pc_G1_C_T$meth.diff) > 15 + 2*sd(abs(DMSlist$DMS_15pc_G1_C_T$meth.diff)))
# outliers_annot_G1 <- as.data.frame(diffAnn_PAR@members)[outliers_G1_final,]
# makeManhattanPlots(DMSfile = DMSlist$DMS_15pc_G1_C_T[outliers_G1_final, ],
#                    annotFile = outliers_annot_G1, GYgynogff = GYgynogff,
#                    mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")
# 
# outliers_G2_G1c_final <- which(abs(DMSlist$DMS_15pc_CC_CT$meth.diff) > 15 + 2*sd(abs(DMSlist$DMS_15pc_CC_CT$meth.diff)))
# outliers_annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)[outliers_G2_G1c_final,]
# makeManhattanPlots(DMSfile = DMSlist$DMS_15pc_CC_CT[outliers_G2_G1c_final, ],
#                    annotFile = outliers_annot_G2_G1c, GYgynogff = GYgynogff,
#                    mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")
# 
# outliers_G2_G1i_final <- which(abs(DMSlist$DMS_15pc_TC_TT$meth.diff) > 15 + 2*sd(abs(DMSlist$DMS_15pc_TC_TT$meth.diff)))
# outliers_annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)[outliers_G2_G1i_final,]
# makeManhattanPlots(DMSfile = DMSlist$DMS_15pc_TC_TT[outliers_G2_G1i_final, ],
#                    annotFile = outliers_annot_G2_G1i, GYgynogff = GYgynogff,
#                    mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")
# ```
# 
# ## I. Focus on CpG positions at parental (Ctrl-Inf) DMS
# 
# Question: how are the beta values in the different groups at the
# parental DMS?
# 
# ```{r}
# ##############
# ## Prepare dataset
# PM_G1 <- getPMdataset(uniteCov = uniteCov6_G1_woSexAndUnknowChrOVERLAP, MD = fullMetadata_PAR, gener="parents")
# PM_G2 <- getPMdataset(uniteCov = uniteCov14_G2_woSexAndUnknowChrOVERLAP, MD = fullMetadata_OFFS, gener="offspring")
# ```
# 
# ### What is the relative contribution of methylation to brother pair & paternal treatment?
# 
# Test of VCA: variance component analysis
# <https://cran.r-project.org/web/packages/VCA/vignettes/VCA_package_vignette.html>
# 
# #### Hypo methylation
# 
# ```{r, fig.width=10}
# PM_G2_mean_hypo <- PM_G2[PM_G2$hypohyper %in% "hypo", ] %>%
#   group_by(brotherPairID, G1_trt, G2_trt, ID) %>%
#   dplyr::summarize(MeanBetaValue = mean(BetaValue, na.rm=TRUE)) %>% data.frame()
# 
# varPlot(form = MeanBetaValue~(G1_trt* G2_trt*brotherPairID), Data = PM_G2_mean_hypo,
#         MeanLine=list(var=c("G1_trt", "G2_trt"),
#                       col=c("white", "blue"), lwd=c(2,2)),
#         BG=list(var="G2_trt", col=paste0("gray", c(80, 90))),
#         YLabel=list(cex = .8, text="Mean beta value at parDMS \n hypomethylated upon infection"))
# ```
# 
# ```{r}
# myfitVCA_hypo <- fitVCA(form = MeanBetaValue~(G1_trt* G2_trt*brotherPairID), Data = PM_G2_mean_hypo)
# 
# ### Real values
# trtEffect <- sum(myfitVCA_hypo$aov.tab[2:4, 5])
# genEffect <- sum(myfitVCA_hypo$aov.tab[5:8, 5])
# error <- sum(myfitVCA_hypo$aov.tab[9, 5])
# realValHypoVCA <- data.frame(trtEffect=trtEffect, genEffect=genEffect,error=error)
# 
# ### Randomisation
# myrandomVCA <- function(df=PM_G2_mean_hypo){
#   randomDF = df
#   randomDF$G1_trt = sample(PM_G2_mean_hypo$G1_trt, replace = F)
#   randomDF$G2_trt = sample(PM_G2_mean_hypo$G2_trt, replace = F)
#   randomDF$brotherPairID = sample(PM_G2_mean_hypo$brotherPairID, replace = F)
#   myfitVCA <- fitVCA(form = MeanBetaValue~(G1_trt* G2_trt*brotherPairID), Data = randomDF)
#   trtEffect <- sum(myfitVCA$aov.tab[2:4, 5])
#   genEffect <- sum(myfitVCA$aov.tab[5:8, 5])
#   error <- sum(myfitVCA$aov.tab[9, 5])
#   return(data.frame(trtEffect=trtEffect, genEffect=genEffect,error=error))
# }
# 
# # randomHypoVCA = do.call(rbind, lapply(1:1000, function(x) {
# #   df=myrandomVCA(PM_G2_mean_hypo)
# #   df$rep=x
# #   return(df)}))
# 
# # randomHypoVCA = melt(randomHypoVCA, id.vars = "rep")
# # saveRDS(randomHypoVCA, file = "Rdata/randomHypoVCA.RDS")
# randomHypoVCA <- readRDS(file = "Rdata/randomHypoVCA.RDS")
# df2=reshape2::melt(realValHypoVCA)
# 
# sumDF <- randomHypoVCA %>%
#   group_by(variable) %>%
#   dplyr::summarize(value = mean(value)) %>% data.frame()
# 
# ggplot(randomHypoVCA, aes(x=variable, y=value))+
#   geom_boxplot()+
#   geom_jitter(width=.1, alpha=.2)+
#   geom_point(data = df2, col = "red", size = 6)+
#   geom_text(data=sumDF, aes(label=round(value)), col="white")+
#   geom_text(data = df2, aes(label=round(value)), col="white")+
#   theme_cleveland()+
#   ggtitle("VCA with bootstrap N=1000 at hypo-parDMS", subtitle = "red: observed values")
# 
# # estimate 95% confidence intervals, request CI for
# # all variance components via 'VarVC=TRUE'
# VCAinference(myfitVCA_hypo, VarVC=TRUE)
# ```
# 
# #### Hyper methylation
# 
# ```{r, fig.width=10}
# PM_G2_mean_hyper <- PM_G2[PM_G2$hypohyper %in% "hyper", ] %>%
#   group_by(brotherPairID, G1_trt, G2_trt, ID) %>%
#   dplyr::summarize(MeanBetaValue = mean(BetaValue, na.rm=TRUE)) %>% data.frame()
# 
# varPlot(form = MeanBetaValue~(G1_trt* G2_trt*brotherPairID), Data = PM_G2_mean_hyper,
#         MeanLine=list(var=c("G1_trt", "G2_trt"),
#                       col=c("white", "blue"), lwd=c(2,2)),
#         BG=list(var="G2_trt", col=paste0("gray", c(80, 90))),
#         YLabel=list(cex = .8, text="Mean beta value at parDMS \n hypermethylated upon infection"))
# ```
# 
# ```{r}
# myfitVCA_hyper <- fitVCA(form = MeanBetaValue~(G1_trt* G2_trt*brotherPairID), Data = PM_G2_mean_hyper)
# 
# ### Real values
# trtEffect <- sum(myfitVCA_hyper$aov.tab[2:4, 5])
# genEffect <- sum(myfitVCA_hyper$aov.tab[5:8, 5])
# error <- sum(myfitVCA_hyper$aov.tab[9, 5])
# realValHyperVCA <- data.frame(trtEffect=trtEffect, genEffect=genEffect,error=error)
# 
# ### Randomisation
# # randomHyperVCA = do.call(rbind, lapply(1:1000, function(x) {
# #   df=myrandomVCA(PM_G2_mean_hyper)
# #   df$rep=x
# #   return(df)}))
# #
# # randomHyperVCA = melt(randomHyperVCA, id.vars = "rep")
# # saveRDS(randomHyperVCA, file = "Rdata/randomHyperVCA.RDS")
# randomHyperVCA <- readRDS(file = "Rdata/randomHyperVCA.RDS")
# df2=reshape2::melt(realValHyperVCA)
# 
# sumDF <- randomHyperVCA %>%
#   group_by(variable) %>%
#   dplyr::summarize(value = mean(value)) %>% data.frame()
# 
# ggplot(randomHyperVCA, aes(x=variable, y=value))+
#   geom_boxplot()+
#   geom_jitter(width=.1, alpha=.2)+
#   geom_point(data = df2, col = "red", size = 6)+
#   geom_text(data=sumDF, aes(label=round(value)), col="white")+
#   geom_text(data = df2, aes(label=round(value)), col="white")+
#   theme_cleveland()+
#   ggtitle("VCA with bootstrap N=1000 at hyper-parDMS", subtitle = "red: observed values")
# 
# # estimate 95% confidence intervals, request CI for
# # all variance components via 'VarVC=TRUE'
# VCAinference(myfitVCA_hypo, VarVC=TRUE)
# ```
# 
# ## DMS values in parents
# 
# ```{r}
# parmod <- lmer(data = PM_G1, BetaValue ~ meth.diff.parentals : Treatment + (1|CpGSite) + (1|brotherPairID))
# 
# ## check normality of residuals assumption
# qqnorm(resid(parmod))
# qqline(resid(parmod))
# 
# pred <- ggpredict(parmod, terms = c("meth.diff.parentals", "Treatment"))
# plot(pred, add.data = T)+
#   scale_color_manual(values = c("black", "red"))+
#   scale_y_continuous(name = "Beta values")+
#   scale_x_continuous(name = "Methylation difference between infected and control parents in percentage")+
#   ggtitle("Predicted methylation ratio (Beta) values in parents\n as a function of differential methylation between exposed and control groups")+
#   theme_bw()
# ```
# 
# ### Linear model: does the beta value of offspring at DMS depends on treatment Parent x Offspring?
# 
# ```{r}
# modFull <- lmer(BetaValue ~ (G1_trt * G2_trt):hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID),data = PM_G2, REML = F) # REML =F for model comparison
# mod_noG1trt <- lmer(BetaValue ~ G2_trt:hypohyper + (1|CpGSite)+ (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
# mod_noG2trt <-lmer(BetaValue ~ G1_trt:hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
# mod_noInteractions <- lmer(BetaValue ~ (G1_trt + G2_trt):hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
# mod_noHypoHyper <- lmer(BetaValue ~ (G1_trt * G2_trt) + (1|CpGSite) + (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
# 
# ## check normality of residuals assumption
# qqnorm(resid(modFull))
# qqline(resid(modFull))
# 
# ## Likelihood ratio tests for all variables:
# lrtest(modFull, mod_noG1trt) # G1 trt is VERY VERY significant (LRT: χ² (4) = 1163.6, p < 0.001)
# lrtest(modFull, mod_noG2trt) # G2 trt is VERY VERY significant (LRT: χ² (4) = 30.02, p < 0.001) NB that changed when brotherpair is used instead of family!
# lrtest(modFull, mod_noInteractions) # interactions are significant (LRT: χ² (2) = 9.21, p < 0.01)
# lrtest(modFull, mod_noHypoHyper) # hypo/hyper VERY VERY significant (LRT: χ² (4) = 1140, p < 0.001)
# 
# ## Post-hoc tests between treatments
# modFull <- lmer(BetaValue ~ (G1_trt * G2_trt):hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID),data = PM_G2)
# modFull_emmeans <- emmeans(modFull, list(pairwise ~ (G1_trt:G2_trt):hypohyper), adjust = "tukey")
# modFull_emmeans
# 
# P1 <- plot(modFull_emmeans, by = "hypohyper", comparisons = TRUE) +
#   # coord_flip()+
#   theme_bw() +
#   ggtitle("Estimated marginal means of methylation ratio (beta)\n of offspring at parental DMS")+
#   theme(legend.position = "none", axis.title.x = element_blank()) +
#   scale_x_continuous("Beta value (methylation ratio)", limits = c(47,69.5))
# 
# ## NB: Comparison arrows: https://cran.r-project.org/web/packages/emmeans/vignettes/xplanations.html
# ## two estimated marginal means (EMMs) differ significantly if, and only if, their respective comparison arrows do not overlap
# ## These comparison arrows are decidedly not the same as confidence intervals.
# ## Confidence intervals for EMMs are based on the statistical properties of the individual EMMs, whereas comparison arrows
# ## are based on the statistical properties of differences of EMMs.
# 
# ## Add the PARENTAL DMS value
# ## Same test on ALL, G1 and G2 fish
# modFullG1 <- lmer(BetaValue ~ G1_trt:hypohyper + (1|CpGSite) + (1|brotherPairID), data = PM_G1)
# 
# modFullG1_emmeans <- emmeans(modFullG1, list(pairwise ~ G1_trt:hypohyper), adjust = "tukey")
# modFullG1_emmeans
# 
# P2 <- plot(modFullG1_emmeans, by = "hypohyper", comparisons = TRUE) +
#   theme_bw() +
#   ggtitle("Estimated marginal means of methylation ratio (beta)\n of parents at DMS")+
#   theme(legend.position = "none", axis.title.x = element_blank()) +
#   scale_x_continuous("Beta value (methylation ratio)", limits = c(47,69.5))
# 
# ggarrange(P2, P1, labels = c("A", "B"), ncol = 1, nrow = 2)
# ```
# 
# ### Are the nbr of residuals methylation AT PARENTAL DMS different in the 4 G2 trt? (for hypo vs hypermeth)?
# 
# ```{r}
# length(unique(PM_G1$CpGSite))# 3648 positions
# PM_G1 %>% dplyr::count(ID)## NB: not all covered in all samples
# length(unique(PM_G2$CpGSite[PM_G2$hypohyper %in% "hypo"]))# 1176 positions hypomethylated upon parental inf
# length(unique(PM_G2$CpGSite[PM_G2$hypohyper %in% "hyper"]))# 2472 positions hypermethylated upon parental inf
# 
# myfun <- function(X){
#   ## Calculate nbr of CpG hypo/hypermethylated per individual, and nbr of covered CpG:
#   X <- X %>% group_by(ID, Treatment, brotherPairID, clutch.ID, Sex) %>%
#     dplyr::summarise(ncov = n(),
#                      hypoMeth = sum(BetaValue < 0.3),
#                      hyperMeth = sum(BetaValue > 0.7)) %>% data.frame()
#   # Calculate residuals of nbr of methCpG by nbr of covered CpG
#   X$res_Nbr_methCpG_Nbr_coveredCpG_HYPO = residuals(lm(X$hypoMeth ~ X$ncov))
#   X$res_Nbr_methCpG_Nbr_coveredCpG_HYPER = residuals(lm(X$hyperMeth ~ X$ncov))
#   
#   ## Statistical model: is it different in each offspring trt group?
#   mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPO ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
#                data = X, REML = F)
#   mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPO ~ 1 + (1|brotherPairID/clutch.ID) + (1|Sex),
#                data = X, REML = F)
#   print(lrtest(mod1, mod0))
#   
#   ## Post-hoc test:
#   modhypo <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPO ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
#                   data = X)
#   ## pairwise posthoc test
#   print(emmeans(modhypo, list(pairwise ~ Treatment), adjust = "tukey"))
#   
#   mod3 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPER ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
#                data = X, REML = F)
#   mod4 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPER ~ 1 + (1|brotherPairID/clutch.ID) + (1|Sex),
#                data = X, REML = F)
#   print(lrtest(mod3, mod4))
#   
#   ## Post-hoc test:
#   modhyper <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPER ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
#                    data = X)
#   ## pairwise posthoc test
#   print(emmeans(modhyper, list(pairwise ~ Treatment), adjust = "tukey"))
#   
#   ## Plot
#   phypo <- plot(ggpredict(modhypo, terms = c("Treatment")), add.data = TRUE)+
#     scale_y_continuous("Residuals of number of hypomethylated methylated \ncytosines on number of cytosines covered") +
#     ggtitle("Predicted residuals nbr of hypomethylated CpG")+
#     theme_bw()
#   
#   phyper <- plot(ggpredict(modhyper, terms = c("Treatment")), add.data = TRUE)+
#     scale_y_continuous("Residuals of number of hypermethylated methylated \n cytosines on number of cytosines covered") +
#     ggtitle("Predicted residuals nbr of hypermethylated CpG")+
#     theme_bw()
#   return(list(phypo, phyper))
# }
# 
# listplots <- myfun(X = PM_G2[PM_G2$hypohyper %in% "hypo",])
# ## NOT significant
# annotate_figure(ggarrange(listplots[[1]], listplots[[2]],ncol = 2, nrow = 1),
#                 top = text_grob("Parental DMS hypomethylated upon infection, in offspring"))
# ```
# 
# ```{r}
# listplots <- myfun(X = PM_G2[PM_G2$hypohyper %in% "hyper",])
# ## Treatment SIGNIFICANT in both excess hypo/hyper methylation **
# 
# # $`pairwise differences of Treatment`
# ## HYPO
# # 1                       estimate   SE   df t.ratio p.value
# # NE_control - E_control     23.71 7.28 10.3   3.257  0.0353
# # NE_control - E_exposed     26.88 7.26 10.3   3.701  0.0172
# 
# ## HYPER
# # 1                       estimate   SE   df t.ratio p.value
# # NE_control - E_control    -24.06 7.36 10.3  -3.269  0.0348
# # NE_control - E_exposed    -27.07 7.34 10.3  -3.687  0.0177
# 
# annotate_figure(ggarrange(listplots[[1]], listplots[[2]],ncol = 2, nrow = 1),
#                 top = text_grob("Parental DMS hypermethylated upon infection, in offspring"))
# ```
# 
# -\> The beta values in parentalDMS in offspring follow the parental
# pattern hypo/hyper methylated upon infection
