# Each script sources the previous script of the pipeline if needed
source("R08_DMSannotation.R")


########### TBC here
### Gene Ontology analysis OVERALL (all genes with G1, G2 or both effects):
# create gene universe
gene_universe <- data.frame(
  subsetByOverlaps(GRanges(annotGff3), GRanges(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP))) %>% # subselect covered CpGs
  filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
  filter(type %in% "gene")  %>% # keep all the 7404 genes with GO terms
  dplyr::select(c("Name", "Ontology_term")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, it's from Interproscan after Maker. Sticklebacks (biomart) have 82701 IEA and 63 ISS.
  relocate("Ontology_term","go_linkage_type","Name") %>%
  unnest(Ontology_term) %>% # one GO per line (was a list before in this column)
  data.frame()

gene_universe$Name %>% unique %>% length #7404 genes

# Create gene set collection
goFrame <- GOFrame(gene_universe, organism="Gasterosteus aculeatus")
goAllFrame <- GOAllFrame(goFrame)
gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())

# **IMPORTANT NOTE from Mel: why conditional hypergeometric test?** The GO ontology is set up as a directed acyclic graph, where a parent term is comprised of all its child terms. If you do a standard hypergeometric, you might e.g., find 'positive regulation of kinase activity' to be significant. If you then test 'positive regulation of catalytic activity', which is a parent term, then it might be significant as well, but only because of the terms coming from positive regulation of kinase activity.
# The conditional hypergeometric takes this into account, and only uses those terms that were not already significant when testing a higher order (parent) term.

# select genes which contains DMS only associated with ONE effect:
GO_G1only = makedfGO(DMS_G1onlyEffect_4BPmin_ANNOT %>%
                       distinct(feature.name,.keep_all = TRUE), gene_universe, 
                     effect = paste0(length(unique(DMS_G1onlyEffect_4BPmin_ANNOT$feature.name)),
                                     " genes with DMS associated with G1 effect only"))

GO_G2only = makedfGO(DMS_G2onlyEffect_4BPmin_ANNOT %>%
                       distinct(feature.name,.keep_all = TRUE) , gene_universe, 
                     effect = paste0(length(unique(DMS_G2onlyEffect_4BPmin_ANNOT$feature.name)),
                                     " genes with DMS associated with G2 effect only"))

GO_G1G2addit = makedfGO(DMS_G1G2additiveEffect_4BPmin_ANNOT %>%
                          distinct(feature.name,.keep_all = TRUE), gene_universe, 
                        effect = paste0(length(unique(DMS_G1G2additiveEffect_4BPmin_ANNOT$feature.name)),
                                        " genes with DMS associated with additive effect"))

GO_G1G2inter = makedfGO(DMS_G1G2interactionEffect_4BPmin_ANNOT %>%
                          distinct(feature.name,.keep_all = TRUE), gene_universe, 
                        effect = paste0(length(unique(DMS_G1G2interactionEffect_4BPmin_ANNOT$feature.name)),
                                        " genes with DMS associated with interaction effect"))

dfGO = rbind(GO_G1only, GO_G2only, GO_G1G2addit, GO_G1G2inter)

### GO plot
GOplot <- dfGO %>% ggplot(aes(x=Effect, y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name="adjusted\np-value", low = "red", high = "blue") +
  scale_size_continuous(name = "% of genes")+
  theme_bw() + ylab("") + xlab("Treatments comparison") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), legend.position="left") + # grey box for legend
  # scale_y_discrete(labels = function(x) str_wrap(x, width = 10))+ # split too long GO names in half
  facet_grid(fct_inorder(GO.category)~., scales="free",space = "free")+
  scale_y_discrete(limits=rev, # revers axis to have alphabetical order
                   labels = function(x) str_wrap(x, width = 40)) # split too long GO names in half

pdf(GOplot, file = "../../dataOut/Supplementary_GOplot_figS3VennGenes.pdf", width = 6, height = 18)
GOplot
dev.off()

## Hyp: PATERNAL EFFECT ONLY CpG persist despite G2 trt
# Test: is the correlation between G1 and G2 methylation at these CpG stronger than for the other effects?

length(DMS_G1onlyEffect_4BPmin) # 1640 positions different ONLY following paternal treatments
length(DMS_G2onlyEffect_4BPmin) # 309 positions different ONLY following offspring treatments
length(DMS_G1G2additiveEffect_4BPmin) # 173 positions
length(DMS_G1G2interactionEffect_4BPmin) # 151 positions 

getCorG1G2methByEffect <- function(myeffect){
  # G1
  A=methylKit::select(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP, 
                      which(paste(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$start) %in% 
                              myeffect))
  # G2
  B=methylKit::select(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 
                      which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$start) %in% 
                              myeffect))
  
  # Methylation values:
  DMS=paste(A$chr, A$start)
  
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
  AB=na.omit(AB[c("G1methylation", "G2methylation")])
  # When the variables are not continuous but could be ranked then we do not use pearson correlation 
  # coefficient to find the linear relationship, in this case spearman correlation coefficient comes 
  # into the scene. Since the spearman correlation coefficient considers the rank of values, the 
  # correlation test ignores the same ranks to find the p-values as a result we get the warning
  # “Cannot compute exact p-value with ties”. This can be avoided by using exact = FALSE inside the cor.test function.
  
  return(list(cor.test(x=AB$G1methylation, y=AB$G2methylation, method="spearman"),
              spearmanRho(x=AB$G1methylation, y=AB$G2methylation, method="spearman", R=1000,ci=T)))
}

# NB: Takes 4 minutes to run the BS each time
# getCorG1G2methByEffect(myeffect = DMS_G1onlyEffect_4BPmin)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# 0.481[0.477-0.484 95%CI]

# getCorG1G2methByEffect(myeffect = DMS_G2onlyEffect_4BPmin)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# 0.407[0.398-0.417 95%CI]

# getCorG1G2methByEffect(myeffect = DMS_G1G2additiveEffect_4BPmin)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# 0.367[0.354-0.379 95%CI]

# getCorG1G2methByEffect(myeffect = DMS_G1G2interactionEffect_4BPmin)
# Spearman's rank correlation rho $ 95% CI (BC 1000)
# 0.391[0.378-0.404 95%CI]

## And as a regression? We put only this one in the paper, it's the most sophisticated version of Spearman
DMSalleffectsDF = data.frame(DMS=c(DMS_G1onlyEffect_4BPmin, DMS_G2onlyEffect_4BPmin, DMS_G1G2additiveEffect_4BPmin, DMS_G1G2interactionEffect_4BPmin),
                             effect=c(rep("G1", length(DMS_G1onlyEffect_4BPmin)),
                                      rep("G2", length(DMS_G2onlyEffect_4BPmin)),
                                      rep("addit", length(DMS_G1G2additiveEffect_4BPmin)),
                                      rep("inter", length(DMS_G1G2interactionEffect_4BPmin))))

A=methylKit::select(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP, 
                    which(paste(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$start) %in% 
                            DMSalleffectsDF$DMS))
# G2
B=methylKit::select(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 
                    which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$start) %in% 
                            DMSalleffectsDF$DMS))
# Methylation values:
DMS=paste(A$chr, A$start)

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
AB=merge(AB, DMSalleffectsDF)
AB$SampleID=as.character(AB$SampleID)

# level effects
AB$effect <- factor(AB$effect,levels = c("G1", "G2", "addit", "inter"))

AB=na.omit(AB[c("G1methylation", "G2methylation", "effect", "brotherPairID", "SampleID")])

mod1=lme4::lmer(G2methylation~G1methylation * effect + (1|brotherPairID/SampleID), data=AB)
mod2=lme4::lmer(G2methylation~G1methylation + effect + (1|brotherPairID/SampleID), data=AB)
mod3=lme4::lmer(G2methylation~G1methylation + (1|brotherPairID/SampleID), data=AB)
mod4=lme4::lmer(G2methylation~effect + (1|brotherPairID/SampleID), data=AB)

lmtest::lrtest(mod1, mod2)
# Likelihood ratio test
# 
# Model 1: G2methylation ~ G1methylation * effect + (1 | brotherPairID/SampleID)
# Model 2: G2methylation ~ G1methylation + effect + (1 | brotherPairID/SampleID)
# #Df   LogLik Df  Chisq Pr(>Chisq)    
# 1  11 -1123068                         
# 2   8 -1123235 -3 335.32  < 2.2e-16 ***

lmtest::lrtest(mod1, mod3)
# Likelihood ratio test
# 
# Model 1: G2methylation ~ G1methylation * effect + (1 | brotherPairID/SampleID)
# Model 2: G2methylation ~ G1methylation + (1 | brotherPairID/SampleID)
# #Df   LogLik Df  Chisq Pr(>Chisq)    
# 1  11 -1123068                         
# 2   5 -1123334 -6 532.44  < 2.2e-16 ***

lmtest::lrtest(mod1, mod4)
# Likelihood ratio test
# 
# Model 1: G2methylation ~ G1methylation * effect + (1 | brotherPairID/SampleID)
# Model 2: G2methylation ~ effect + (1 | brotherPairID/SampleID)
# #Df   LogLik Df Chisq Pr(>Chisq)    
# 1  11 -1123068                        
# 2   7 -1152337 -4 58539  < 2.2e-16 ***

mod1

library(emmeans)
emmp <- emmeans(mod1, pairwise ~ effect, adjust = "tukey")
summary(emmp, infer=TRUE)$contrast
# contrast      estimate    SE  df asymp.LCL asymp.UCL z.ratio p.value
# G1 - G2          1.885 0.174 Inf     1.438     2.332  10.830  <.0001
# G1 - addit       1.397 0.218 Inf     0.837     1.957   6.410  <.0001
# G1 - inter       2.262 0.233 Inf     1.665     2.860   9.723  <.0001
# G2 - addit      -0.488 0.262 Inf    -1.159     0.184  -1.865  0.2434
# G2 - inter       0.378 0.274 Inf    -0.326     1.081   1.379  0.5127
# addit - inter    0.865 0.304 Inf     0.085     1.646   2.849  0.0228

dfTukey = summary(emmp, infer=TRUE)$contrast
dfTukey = dfTukey %>% 
  dplyr::mutate_if(is.numeric, round, digits = 4)

write.csv(dfTukey, file = "../../dataOut/SuppTableTukeyG1G2.csv", row.names = F)

# Set up modelplot
modelPlot <- plot_model(mod1, type = "pred", terms = c("G1methylation", "effect"))+
  geom_abline(y=1)+
  theme_bw()+
  xlim(1,100)+
  ylim(0,100)+
  xlab("Methylation value of father (%)")+
  ylab("Methylation value of offspring (%)")+
  ggtitle(NULL)+
  guides(color = "none") +
  # ggtitle("Predicted values of offspring CpG methylation\n associated with each effect versus their father's")+
  scale_color_manual(values = c("#e69f00", "#56b4e9", "#009e73", "#cc79a7"))+
  scale_fill_manual(values = c("#e69f00", "#56b4e9", "#009e73", "#cc79a7"))+
  theme(plot.margin = margin()) + theme_bw() +
  theme(legend.position = "none")

# Define marginal histogram
marginal_distribution2 <- function(x, var, group) {
  ggplot(x, aes_string(x = var, fill = group)) +
    geom_density(adjust=1.5, alpha=.3) +
    # geom_density(adjust=1.5, alpha=.7, aes(y=..count..))+
    # geom_histogram(col="black")+
    guides(fill = "none") +
    scale_fill_manual(values = c("#e69f00", "#56b4e9", "#009e73", "#cc79a7")) +
    # theme_void() +
    theme(plot.margin = margin())
}

marginal_distribution3 <- function(x, var, group) {
  ggplot(x, aes_string(x = var, fill = group)) +
    geom_histogram(col="black")+
    guides(fill = "none") +
    scale_fill_manual(values = c("#e69f00", "#56b4e9", "#009e73", "#cc79a7")) +
    # theme_void() +
    theme(plot.margin = margin())
}

# Set up marginal histograms & density
x_dens <- marginal_distribution2(AB, "G1methylation", "effect")
x_hist <- marginal_distribution3(AB, "G1methylation", "effect")
y_dens <- marginal_distribution2(AB, "G2methylation", "effect")+coord_flip()
y_hist <- marginal_distribution3(AB, "G2methylation", "effect")+coord_flip()

# Align histograms with scatterplot
aligned_x_hist <- align_plots(x_hist, modelPlot, align = "v")[[1]]
aligned_y_hist <- align_plots(y_hist, modelPlot, align = "h")[[1]]

# Align density with scatterplot
aligned_x_dens <- align_plots(x_dens, modelPlot, align = "v")[[1]]
aligned_y_dens <- align_plots(y_dens, modelPlot, align = "h")[[1]]

pdf(file = "../../dataOut/plotmodelG1methG2meth.pdf", width = 12, height = 6)
# Arrange plots
cowplot::plot_grid(cowplot::plot_grid(
  aligned_x_hist, NULL, modelPlot, aligned_y_hist, ncol = 2, nrow = 2, rel_heights = c(0.3, 1), rel_widths = c(1, 0.3)),
  cowplot::plot_grid(
    aligned_x_dens, NULL, modelPlot, aligned_y_dens, ncol = 2, nrow = 2, rel_heights = c(0.3, 1), rel_widths = c(1, 0.3)
  ), ncol=2)
dev.off()


