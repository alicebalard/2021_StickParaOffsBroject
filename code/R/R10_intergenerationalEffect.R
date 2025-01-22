# Each script sources the previous script of the pipeline if needed
source("R09_GeneOntology.R")

message("R10 starting...\n")

## Are the parental DMS found in offspring DMS? In which effect? Likely infection-induced!
load("../../gitignore/bigdata/PQLseq/fitPQLseqG1_fit_G1_C.T.RData")

## Select based on threshold in R07
C.T_pass = fit_G1_C.T[fit_G1_C.T$converged==T & fit_G1_C.T$qvalue <= myqval & 
                        (fit_G1_C.T$aveDiffMeth_ab>mydif | fit_G1_C.T$aveDiffMeth_ab< -mydif),]

hist(fit_G1_C.T$aveDiffMeth_ab)
## no difference above 5%

C.T_pass = fit_G1_C.T[fit_G1_C.T$converged==T & fit_G1_C.T$qvalue <= myqval,]

## Venn diagram of overlaps
ggvenn(data = list(C.T=row.names(C.T_pass),
                   intergen=EffectsDF$pos[EffectsDF$effect %in% "INTERGENERATIONAL"],
                   infind=EffectsDF$pos[EffectsDF$effect %in% "INFECTION_INDUCED"]), 
       fill_color = c("blue", "#EFC000FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4
)


## Select based on threshold
ggvenn(data = 
         list(C.T = fit_G1_C.T[fit_G1_C.T$converged==T & fit_G1_C.T$qvalue <= myqval,"pos"],
              G1 = c(fit_G2_CC.TC[fit_G2_CC.TC$converged==T & fit_G2_CC.TC$qvalue <= myqval,"pos"],
                     fit_G2_CT.TT[fit_G2_CT.TT$converged==T & fit_G2_CT.TT$qvalue <= myqval,"pos"]),
              G2 = c(fit_G2_CC.CT[fit_G2_CC.CT$converged==T & fit_G2_CC.CT$qvalue <= myqval,"pos"],
                     fit_G2_TC.TT[fit_G2_TC.TT$converged==T & fit_G2_TC.TT$qvalue <= myqval,"pos"])))

## Which 8 positions are differentially methylated in all 3 cases?
full =c(unique(fit_G1_C.T[fit_G1_C.T$converged==T & fit_G1_C.T$qvalue <= myqval,"pos"]),
        unique(c(fit_G2_CC.TC[fit_G2_CC.TC$converged==T & fit_G2_CC.TC$qvalue <= myqval,"pos"],
                 fit_G2_CT.TT[fit_G2_CT.TT$converged==T & fit_G2_CT.TT$qvalue <= myqval,"pos"])),
        unique(c(fit_G2_CC.CT[fit_G2_CC.CT$converged==T & fit_G2_CC.CT$qvalue <= myqval,"pos"],
                 fit_G2_TC.TT[fit_G2_TC.TT$converged==T & fit_G2_TC.TT$qvalue <= myqval,"pos"])))

DMS8AllComp <- myHomebrewDMSannotation(
  DMSvec = names(table(full))[table(full) == 3], myannotBed12 = annotBed12, myannotGff3 = annotGff3, rerun = FALSE)

#######################################################
## Hyp: PATERNAL EFFECT ONLY CpG persist despite G2 trt
# Test: is the correlation between G1 and G2 methylation at these CpG stronger
## than for the infection-induced effect?

nrow(EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED",]) # 393
nrow(EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL",]) # 307

## Get AB: methylation values for G1 and G2 for the DMS associated to a given effect
getG1G2methByEffect <- function(DMSeffect){
  # G1
  A=methylKit::select(
    uniteCovHALF_G1_woSexAndUnknowChrOVERLAP,
    which(paste(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr,
                uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$start, sep = "_") %in%
            DMSeffect))
  # G2
  B=methylKit::select(
    uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
    which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr,
                uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$start, sep = "_") %in%
            DMSeffect))
  
  # Methylation values:
  DMS=paste(A$chr, A$start, sep = "_")
  
  A=percMethylation(A) %>% data.frame
  A$DMS = DMS
  B=percMethylation(B) %>% data.frame
  B$DMS = DMS
  
  # Prep G1 df
  A=melt(A)
  names(A)[names(A)%in% "variable"]="SampleID"
  names(A)[names(A)%in% "value"]="G1methylation"
  A=merge(A, fullMetadata_PAR[c("SampleID", "brotherPairID", "trtG1G2")])
  names(A)[names(A)%in% "trtG1G2"]="patTrt"
  A=A[!names(A)%in% "SampleID"] # rm sample ID for fathers, as only one per BP
  
  # Prep G2 df
  B=melt(B)
  names(B)[names(B)%in% "variable"]="SampleID"
  names(B)[names(B)%in% "value"]="G2methylation"
  B=merge(B, fullMetadata_OFFS[c("SampleID", "brotherPairID", "trtG1G2", "patTrt")])
  B$patTrt[B$patTrt %in% "controlP"]="Control"
  B$patTrt[B$patTrt %in% "infectedP"]="Exposed"
  
  AB=merge(A, B)
  
  # Rm NA
  AB=AB[!is.na(AB$G1methylation) & !is.na(AB$G2methylation), ]
  return(AB)
}

AB_infectionInduced = getG1G2methByEffect(
  DMSeffect = EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED","DMS"]
)
AB_infectionInduced$effect = "infection-induced"

AB_intergenerational = getG1G2methByEffect(
  DMSeffect = EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL","DMS"]
)
AB_intergenerational$effect = "intergenerational"

AB=rbind(AB_intergenerational, AB_infectionInduced)

## important order for the plots
AB$effect = factor(AB$effect, levels = c("infection-induced", "intergenerational"))

# Test: is the correlation between G1 and G2 methylation at intergenerational CpG stronger
## than for the infection-induced effect?

###################################################
# When the variables are not continuous but could be ranked then we do not use pearson correlation
# coefficient to find the linear relationship, in this case spearman correlation coefficient comes
# into the scene. Since the spearman correlation coefficient considers the rank of values, the
# correlation test ignores the same ranks to find the p-values as a result we get the warning
# “Cannot compute exact p-value with ties”. This can be avoided by using exact = FALSE inside the cor.test function.

getCorG1G2methByEffect <- function(AB){
  list(cor.test(x=AB$G1methylation, y=AB$G2methylation, method="spearman", exact=F),
       spearmanRho(x=AB$G1methylation, y=AB$G2methylation, method="spearman", R=1000,ci=T))
}

# NB: Takes a few minutes to run
# getCorG1G2methByEffect(AB_infectionInduced)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# rho = 0.624 [0.617-0.632] 95%CI

# getCorG1G2methByEffect(AB_intergenerational)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# 0.621[0.611-0.63 95%CI]

mod1=lme4::lmer(G2methylation~G1methylation * effect + (1|brotherPairID/SampleID), data=AB)
mod2=lme4::lmer(G2methylation~G1methylation + effect + (1|brotherPairID/SampleID), data=AB)
mod3=lme4::lmer(G2methylation~G1methylation + (1|brotherPairID/SampleID), data=AB)
mod4=lme4::lmer(G2methylation~effect + (1|brotherPairID/SampleID), data=AB)

lmtest::lrtest(mod1, mod2)
# no effect of interaction
lmtest::lrtest(mod1, mod3)
# effect of "effect" on G2 methylation
lmtest::lrtest(mod1, mod4)
# effect of "G1methylation" on G2 methylation

plot_model(mod2, type = "pred", terms = c("G1methylation", "effect"))+
  geom_abline(slope = 1, linetype = 2)

###############################
## Check some genes of interest
ggplot(AB, aes(x=trtG1G2, y=G1methylation, col=trtG1G2))+
  geom_violin()+
  facet_wrap(.~effect)+
  scale_color_manual(values = colOffs) 

plotG1byG2pergene <- function(gene){
  df = AB[AB$DMS %in% EffectsDF_ANNOT[EffectsDF_ANNOT$GeneSymbol %in% gene,"DMS"],]
  
  # Find the convex hull of the points being plotted 
  hull <- df %>% group_by(trtG1G2, DMS) %>% 
    slice(chull(G1methylation, G2methylation))
  
  ## plot
  ggplot(df,
         aes(x=G1methylation, y = G2methylation, col = trtG1G2, group = trtG1G2, fill=trtG1G2))+
    scale_color_manual(values = colOffs) +
    scale_fill_manual(values = colOffs) +
    facet_wrap(DMS~., ncol = 2)+
    # geom_smooth(method = "lm", alpha=.1) +
    geom_abline(slope = 1, linetype = 2)+
    # geom_rect(aes(xmin=0, xmax=30, ymin=0, ymax=30), col = "lightgrey", fill = "lightgrey")+
    # geom_rect(aes(xmin=0, xmax=30, ymin=70, ymax=100), col = "lightgrey", fill = "lightgrey")+
    # geom_rect(aes(xmin=70, xmax=100, ymin=0, ymax=30), col = "lightgrey", fill = "lightgrey")+
    # geom_rect(aes(xmin=70, xmax=100, ymin=70, ymax=100), col = "lightgrey", fill = "lightgrey")+
    geom_polygon(data = hull, alpha = 0.2, 
                 aes(fill = trtG1G2,colour = trtG1G2)) +
    geom_point(size = 3, alpha =.5)
}

## Intergenerational
plotG1byG2pergene("CD4") 
plotG1byG2pergene("Chd5") 
plotG1byG2pergene("Mosmo")
plotG1byG2pergene("B3GAT1")
plotG1byG2pergene("Stk24")

## infection induced
plotG1byG2pergene("bmp2") 
plotG1byG2pergene("LRFN2")
plotG1byG2pergene("Gzf1")
plotG1byG2pergene("PCDH7")
plotG1byG2pergene("Rab11fip5")
plotG1byG2pergene("TRIM16")
plotG1byG2pergene("ST6GALNAC3")
plotG1byG2pergene("NDFIP2")

message("R10 done. \n")
