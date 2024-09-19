# Each script sources the previous script of the pipeline if needed
source("R09_GeneOntology.R")

message("R10 starting...\n")

## produces dataOut/fig/FigS2_plotmodelG1methG2meth

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

## See if the correlation is different for both effects
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
getCorG1G2methByEffect(AB_infectionInduced)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# rho = 0.624 [0.617-0.632] 95%CI

getCorG1G2methByEffect(AB_intergenerational)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# 0.621[0.611-0.63 95%CI]

## And as a regression? We put only this one in the paper, it's the most sophisticated version of Spearman
AB=rbind(AB_intergenerational, AB_infectionInduced)

## important order for the plots
AB$effect = factor(AB$effect, levels = c("infection-induced", "intergenerational"))

mod1=lme4::lmer(G2methylation~G1methylation * effect + (1|brotherPairID/SampleID), data=AB)
mod2=lme4::lmer(G2methylation~G1methylation + effect + (1|brotherPairID/SampleID), data=AB)
mod3=lme4::lmer(G2methylation~G1methylation + (1|brotherPairID/SampleID), data=AB)
mod4=lme4::lmer(G2methylation~effect + (1|brotherPairID/SampleID), data=AB)

lmtest::lrtest(mod1, mod2)
# no effect of interaction (p=0.452)

lmtest::lrtest(mod1, mod3)
# effect of "effect" on G2 methylation  chi2 = 42.033, df=2. p<0.001

lmtest::lrtest(mod1, mod4)
# effect of "G1methylation" on G2 methylation chi2 = 27561, df=2. p<0.001

### We select mod2

library(emmeans)
emmp <- emmeans(mod2, pairwise ~ effect, adjust = "tukey")
summary(emmp, infer=TRUE)$contrast
# contrast                                estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
# (infection-induced) - intergenerational    -1.05 0.157 Inf     -1.35    -0.739  -6.670  <.0001

dfTukey = summary(emmp, infer=TRUE)$contrast
dfTukey = dfTukey %>% 
  dplyr::mutate_if(is.numeric, round, digits = 4)

#                               contrast estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
# 1 (infection-induced) - intergenerational  -1.0469 0.157 Inf   -1.3545   -0.7392 -6.6698       0

# Set up modelplot
modelPlot <- plot_model(mod2, type = "pred", terms = c("G1methylation", "effect"))+
  geom_abline(slope = 1, linetype = 2)+
  theme_bw()+
  guides(col = "none") +
  xlim(0,100) + ylim(0,100)+
  # xlim(30,70)+
  # ylim(30,70)+
  xlab("Methylation value of father (%)")+
  ylab("Methylation value of offspring (%)")+
  ggtitle(NULL)+
  scale_color_manual(values = colEffects)+
  scale_fill_manual(values = colEffects)+
  theme(plot.margin = margin())

# Define marginal histogram
marginal_distribution2 <- function(x, var, group) {
  ggplot(x, aes_string(x = var, fill = group)) +
    geom_density(adjust=1.5, alpha=.3) +
    guides(fill = "none") +
    scale_fill_manual(values = colEffects) +
    theme(axis.text = element_text(size = 7))+ # reduce size text in axis
    theme(axis.ticks = element_blank())+
    theme(plot.margin = margin())
}

marginal_distribution3 <- function(x, var, group) {
  ggplot(x, aes_string(x = var, fill = group)) +
    geom_histogram(col="black")+
    guides(fill = "none") +
    scale_fill_manual(values = colEffects) +
    theme(axis.text = element_text(size = 7))+ # reduce size text in axis
    theme(axis.ticks = element_blank())+
    theme(plot.margin = margin())
}

# Set up marginal histograms & density
x_dens <- marginal_distribution2(AB, "G1methylation", "effect")+ xlab("")+
  theme(axis.text.x = element_blank())
x_hist <- marginal_distribution3(AB, "G1methylation", "effect")+ xlab("")+
  theme(axis.text.x = element_blank())
y_dens <- marginal_distribution2(AB, "G2methylation", "effect")+ 
  coord_flip()+xlab("")+
  theme(axis.text.y = element_blank()) 
y_hist <- marginal_distribution3(AB, "G2methylation", "effect")+ 
  coord_flip() + xlab("")+
  theme(axis.text.y = element_blank())


# Align histograms with scatterplot
aligned_x_hist <- align_plots(x_hist, modelPlot, align = "v")[[1]]
aligned_y_hist <- align_plots(y_hist, modelPlot, align = "h")[[1]]

# Align density with scatterplot
aligned_x_dens <- align_plots(x_dens, modelPlot, align = "v")[[1]]
aligned_y_dens <- align_plots(y_dens, modelPlot, align = "h")[[1]]

pdf(file = "../../dataOut/fig/FigS2_plotmodelG1methG2meth.pdf", width = 12, height = 6)
cowplot::plot_grid(
  aligned_x_hist, NULL, NULL,
  aligned_x_dens, NULL, NULL,
  modelPlot, aligned_y_dens, aligned_y_hist,
  ncol = 3, nrow = 3, 
  rel_heights = c(0.2, 0.2, 0.6), rel_widths = c(0.6, 0.2, 0.2)
)
dev.off()

message("R10 done. \n")