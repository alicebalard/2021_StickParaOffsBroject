# Each script sources the previous script of the pipeline if needed
source("R07_DMSeffects.R")

source("homebrewDMSannotation.R")# needed for annotation, slight modification of genomation
## Load file containing length of each gynogen chromosomes
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")

## Bed12 file loaded in script R06

## Load curated gff file
annotGff3 <- rtracklayer::readGFF("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff")
#########################

### Annotate the different DMS groups
DMS_G1onlyEffect_4BPmin_ANNOT = myHomebrewDMSannotation(DMSvec = DMS_G1onlyEffect_4BPmin,
                                                        myannotBed12 = annotBed12, myannotGff3 = annotGff3)
DMS_G1onlyEffect_4BPmin_ANNOT=DMS_G1onlyEffect_4BPmin_ANNOT %>% mutate(effect = "G1")
# "check that these features are identical:"
# "gasAcul16628-RA" "gasAcul16627-RA" -> overlapping: Protein of unknown function & Tp63
# "gasAcul15294-RA" "gasAcul15295-RA" -> overlapping: Cdh23 & Vsir (immune!)
# "gasAcul19985-RA" "gasAcul19984-RA" -> overlapping: ST3GAL1 & ST3GAL1, all good for this one
DMS_G2onlyEffect_4BPmin_ANNOT = myHomebrewDMSannotation(DMSvec = DMS_G2onlyEffect_4BPmin,
                                                        myannotBed12 = annotBed12, myannotGff3 = annotGff3)
DMS_G2onlyEffect_4BPmin_ANNOT=DMS_G2onlyEffect_4BPmin_ANNOT %>% mutate(effect = "G2")

DMS_G1G2additiveEffect_4BPmin_ANNOT = myHomebrewDMSannotation(DMSvec = DMS_G1G2additiveEffect_4BPmin,
                                                              myannotBed12 = annotBed12, myannotGff3 = annotGff3)
DMS_G1G2additiveEffect_4BPmin_ANNOT=DMS_G1G2additiveEffect_4BPmin_ANNOT %>% mutate(effect = "addit")

DMS_G1G2interactionEffect_4BPmin_ANNOT = myHomebrewDMSannotation(DMSvec = DMS_G1G2interactionEffect_4BPmin,
                                                                 myannotBed12 = annotBed12, myannotGff3 = annotGff3)
DMS_G1G2interactionEffect_4BPmin_ANNOT=DMS_G1G2interactionEffect_4BPmin_ANNOT%>% mutate(effect = "inter")
# "check that these features are identical:"
# "gasAcul04256-RA" "gasAcul04255-RA" -> overlapping: Protein if unknown function & Proteolipid protein DM beta

# Plot a Venn diagram to see genes in common
pdf(file = "../../dataOut/DMSgroupsVenn_geneLevel.pdf", width = 7, height = 6)
ggVennDiagram(list("G1" = DMS_G1onlyEffect_4BPmin_ANNOT$feature.name, "G2" = DMS_G2onlyEffect_4BPmin_ANNOT$feature.name,
                   "addit" = DMS_G1G2additiveEffect_4BPmin_ANNOT$feature.name, "inter" = DMS_G1G2interactionEffect_4BPmin_ANNOT$feature.name),
              label_alpha = 0) + scale_color_manual(values = c(1,1,1,1))+
  scale_fill_gradient(low="white",high = "yellow") + theme(legend.position = "none") + ggtitle("Genes in each effect")
dev.off()

## NB some genes have DMS in different effects!
# A gene has DMSs in the 4 effects!
geneAll4 = intersect(intersect(intersect(DMS_G1onlyEffect_4BPmin_ANNOT$feature.name, DMS_G2onlyEffect_4BPmin_ANNOT$feature.name),
                               DMS_G1G2additiveEffect_4BPmin_ANNOT$feature.name), DMS_G1G2interactionEffect_4BPmin_ANNOT$feature.name)
DMS_G1onlyEffect_4BPmin_ANNOT[DMS_G1onlyEffect_4BPmin_ANNOT$feature.name %in% geneAll4,]
## FKBP3 & DPP6

## Make a summary table with genes and to which effect they belong:
allDMSAnnot = rbind(DMS_G1onlyEffect_4BPmin_ANNOT,
                    DMS_G2onlyEffect_4BPmin_ANNOT,
                    DMS_G1G2additiveEffect_4BPmin_ANNOT,
                    DMS_G1G2interactionEffect_4BPmin_ANNOT)

# add summary in select
allDMSAnnot = allDMSAnnot %>% 
  dplyr::select(GeneSymbol, feature.name, chrom, start.gene, end.gene, geneLengthkb, nDMSperGene, effect, summary) %>%
  unique %>% tidyr::spread(key = effect, value = nDMSperGene) %>% 
  dplyr::mutate(nDMS = rowSums(across(c(G1, G2, 'addit', 'inter')), na.rm = T),
                nDMSperGeneKb=round(nDMS/geneLengthkb, 2)) %>%
  rowwise() %>% dplyr::mutate(effect = paste0(c("G1", "G2", "addit", "inter")[
    !is.na(c_across(all_of(c("G1", "G2", "addit", "inter"))))], collapse = ' ')) %>% # add effect
  data.frame() %>%
  arrange(desc(nDMSperGeneKb))#arrange in descending order
# reorder columns with summary at the end
allDMSAnnot = allDMSAnnot %>% 
  dplyr::select(c(names(allDMSAnnot)[!names(allDMSAnnot) %in% "summary"], "summary"))

# Get the first gene for each effect, i.e. the one with the most DMS/kb
allDMSAnnot_top = allDMSAnnot %>% arrange(effect)
allDMSAnnot_top = allDMSAnnot_top[!duplicated(allDMSAnnot_top$effect),]%>% arrange(desc(nDMSperGeneKb))

# Write out
write.csv(allDMSAnnot, file = "../../dataOut/allDMSAnnot_supTabS1.csv", row.names = F)
write.csv(allDMSAnnot_top, file = "../../dataOut/allDMSAnnot_top.csv", row.names = F)

#### Focus on gene FKBP3
plotGeneTarget <- function(myTargetGene, myannotBed12=annotBed12){
  # plotdf
  dfplot = rbind(DMS_G1onlyEffect_4BPmin_ANNOT[DMS_G1onlyEffect_4BPmin_ANNOT$GeneSymbol %in% myTargetGene,],
                 DMS_G2onlyEffect_4BPmin_ANNOT[DMS_G2onlyEffect_4BPmin_ANNOT$GeneSymbol %in% myTargetGene,],
                 DMS_G1G2additiveEffect_4BPmin_ANNOT[DMS_G1G2additiveEffect_4BPmin_ANNOT$GeneSymbol %in% myTargetGene,],
                 DMS_G1G2interactionEffect_4BPmin_ANNOT[DMS_G1G2interactionEffect_4BPmin_ANNOT$GeneSymbol %in% myTargetGene,])
  # Find TSS position of the gene
  dfplot$TSSpos = myannotBed12$TSSes[myannotBed12$TSSes$name %in% dfplot$feature.name]@ranges@start
  # Set TSS as origin
  dfplot$start_distToTSS = dfplot$start - dfplot$TSSpos 
  dfplot$end_distToTSS = dfplot$end - dfplot$TSSpos 
  dfplot$start.gene_distToTSS = dfplot$start.gene - dfplot$TSSpos 
  dfplot$end.gene_distToTSS = dfplot$end.gene - dfplot$TSSpos 
  # Reorder effects factor for legend
  dfplot$effect <- factor(dfplot$effect, levels = c("G1", "G2", "addit", "inter"))
  # Prepare rectangles
  mini=min(dfplot$start.gene_distToTSS, dfplot$start_distToTSS)
  maxi=max(dfplot$end.gene_distToTSS, dfplot$end_distToTSS)+100
  start.gene_distToTSS = unique(dfplot$start.gene_distToTSS)
  end.gene_distToTSS = unique(dfplot$end.gene_distToTSS)
  # Plot
  plotGeneTarget= ggplot(dfplot) +
    geom_rect(xmin=mini, xmax=maxi, ymin=0, ymax =.5, fill="#bfb6b6")+
    geom_rect(aes(xmin=start.gene_distToTSS, xmax=end.gene_distToTSS, ymin=0, ymax =.5), fill = "black")+
    geom_point(aes(x = start_distToTSS, y = .8, col = effect, pch=featureType), size = 5) +
    geom_segment(aes(x = start_distToTSS, xend = start_distToTSS, y=0, yend=.8, col = effect)) +
    geom_segment(aes(x = 0, xend = 0, y=0, yend=.5), col = "red", size = 3) + # TSS
    theme_blank() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    labs(title = paste(unique(dfplot$GeneSymbol), ":", dfplot$description),
         subtitle = str_wrap(dfplot$summary, width = 150))+
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))#cb friendly palette
  ####### And by brother pairs
  dfplot_BP = merge(df_effects_full[df_effects_full$DMS %in% dfplot$DMS,c("BP", "DMS", "effectBPlevel")], dfplot)
  # Rm potision with no effect in this BP
  dfplot_BP=dfplot_BP[!is.na(dfplot_BP$effectBPlevel),]
  # Same order of Father's family as in figure 1 (clusters)
  dfplot_BP$BP = factor(dfplot_BP$BP, levels = c("BP05", "BP31", "BP04", "BP30", "BP16", "BP34", "BP49","BP46"))
  # Plot
  plotGeneTargetBP = ggplot(dfplot_BP) +
    geom_rect(xmin=mini, xmax=maxi, ymin=0, ymax =.5, fill="#bfb6b6")+
    geom_rect(aes(xmin=start.gene_distToTSS, xmax=end.gene_distToTSS, ymin=0, ymax =.5), fill = "black")+
    geom_point(aes(x = start_distToTSS, y = .8, col = effectBPlevel, pch=featureType), size = 5) +
    geom_segment(aes(x = start_distToTSS, xend = start_distToTSS, y=0, yend=.8, col = effectBPlevel)) +
    geom_segment(aes(x = 0, xend = 0, y=0, yend=.5), col = "red", size = 3) + # TSS
    theme_blank() +
    facet_grid(BP~.)+ theme(panel.spacing = unit(1.5, "lines"))+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    scale_y_continuous(expand=expansion(mult=c(0,0.15))) # increase space up
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))#cb friendly palette
  
  return(list(myTargetGene_DMSdf=dfplot, plotGeneTarget=plotGeneTarget, plotGeneTargetBP=plotGeneTargetBP))
}

## All 4 effects:
P = plotGeneTarget(myTargetGene = "FKBP3")

P$myTargetGene_DMSdf

# Zoom in 1 
P1 = P$plotGeneTarget + coord_cartesian(xlim=c(3588,3636)) + theme(legend.position = "none") +  
  ggtitle("") + labs(subtitle = "")+ theme(plot.background = element_rect(colour = "black", fill=NA, size=1))
# Zoom in 2
P2 = P$plotGeneTarget + coord_cartesian(xlim=c(4890,5060)) + theme(legend.position = "none") +  
  ggtitle("") + labs(subtitle = "")+ theme(plot.background = element_rect(colour = "black", fill=NA, size=1))

## For all BP, zoomed
Pbp1 = P$plotGeneTargetBP + theme(legend.position = "none") + coord_cartesian(xlim=c(3588,3636))+
  theme(plot.background = element_rect(colour = "black", fill=NA, size=1))
Pbp2 = P$plotGeneTargetBP + theme(legend.position = "none") + coord_cartesian(xlim=c(4890,5060))+ 
  theme(plot.background = element_rect(colour = "black", fill=NA, size=1))

fullplotFKBP3 = cowplot::plot_grid(P$plotGeneTarget,
                                   cowplot::plot_grid(P1, P2, nrow = 1, ncol =2), 
                                   cowplot::plot_grid(Pbp1, Pbp2, nrow = 1, ncol =2), 
                                   nrow = 3, rel_heights = c(1,1,3), labels = c("A", "B", "C"))
fullplotFKBP3

# Plot gene and zooms
pdf("../../dataOut/FKBP3_DMS.pdf", width = 15, height = 15)
fullplotFKBP3
dev.off()

#### Manhattan plot of the genes in the four main effects
# select genes which contains DMS only associated with ONE effect:
allDMSAnnot_uniqueEffect <- allDMSAnnot[rowSums(is.na(allDMSAnnot[c("G1", "G2", "addit", "inter")])) %in% 3,]

## Prepare annotation
annotFile=allDMSAnnot_uniqueEffect %>% dplyr::select(c("start.gene", "end.gene", "nDMSperGeneKb", "effect", "chrom", "GeneSymbol", "effect"))

## Select top 5 genes per effect:
top5perEffect <- annotFile[!is.na(annotFile$GeneSymbol),] %>% 
  arrange(desc(nDMSperGeneKb)) %>% 
  group_by(effect) %>%
  slice(1:5) %>%
  dplyr::pull(GeneSymbol)
annotFile$isTop = ifelse(test = annotFile$GeneSymbol %in% top5perEffect, "top", NA)

## Prepare genome for Manhattan plots:
genome4Manhattan = GYgynogff %>%
  #genome without chrXIX and unknown re-type:
  filter(chrom!="Gy_chrXIX" & chrom!= "Gy_chrUn")%>%
  mutate(chrom_nr=chrom %>% deroman(),
         chrom_order=factor(chrom_nr) %>% as.numeric()) %>% arrange(chrom_order) %>%
  mutate(gstart=lag(length,default=0) %>% cumsum(),
         gend=gstart+length,
         typeBG=LETTERS[2-(chrom_order%%2)],
         gmid=(gstart+gend)/2)

# Prepare data and change gene position to start at the good chromosome
data4Manhattan = dplyr::left_join(annotFile, genome4Manhattan) %>%
  dplyr::mutate(posInPlot=((end.gene+start.gene)/2)+gstart)

# level effects in order
data4Manhattan$effect <- factor(data4Manhattan$effect,levels = c("G1", "G2", "addit", "inter"))

pdf(file = "../../dataOut/ManhattanFig3.pdf", width = 8, height = 5)
# Manhattan plot
ggplot()+
  # add grey background every second chromosome
  geom_rect(data=genome4Manhattan[genome4Manhattan$typeBG %in% "B",],
            aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf), fill=rgb(.9,.9,.9))+  xlab(NULL)+
  scale_x_continuous(name = "Chromosomes", breaks=genome4Manhattan$gmid,labels=genome4Manhattan$chrom %>% str_remove(.,"Gy_chr"),
                     position = "bottom",expand = c(0,0))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ # add frame
  # scale_y_continuous(breaks = seq(0, 20), expand = expansion(mult = 0.5)) + # increase size under plot for labels
  ylab("Number of differentially methylated CpG per gene kb")+
  scale_y_continuous(breaks = 0:5)+
  geom_point(data = data4Manhattan, aes(x=posInPlot, y = nDMSperGeneKb)) +
  geom_label_repel(data = data4Manhattan[data4Manhattan$isTop %in% "top",], 
                   aes(x=posInPlot, y = nDMSperGeneKb, label = GeneSymbol, fill = factor(effect)),
                   color="white",max.overlaps = Inf, box.padding = 0.6, segment.color = "grey")+
  scale_fill_manual(values=c("#e69f00", "#56b4e9", "#009e73", "#cc79a7"))
dev.off()

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
  
  # When the variables are not continuous but could be ranked then we do not use pearson correlation 
  # coefficient to find the linear relationship, in this case spearman correlation coefficient comes 
  # into the scene. Since the spearman correlation coefficient considers the rank of values, the 
  # correlation test ignores the same ranks to find the p-values as a result we get the warning
  # “Cannot compute exact p-value with ties”. This can be avoided by using exact = FALSE inside the cor.test function.
  return(cor.test(x = AB$G1methylation, AB$G2methylation, method = "spearman", exact = FALSE ))
}

getCorG1G2methByEffect(myeffect = DMS_G1onlyEffect_4BPmin)
# Spearman's rank correlation rho
# S = 4.224e+14, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4808556 

getCorG1G2methByEffect(myeffect = DMS_G2onlyEffect_4BPmin)
# Spearman's rank correlation rho
# S = 3.0929e+12, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4071759 

getCorG1G2methByEffect(myeffect = DMS_G1G2additiveEffect_4BPmin)
# Spearman's rank correlation rho
# S = 6.8324e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3671226 

getCorG1G2methByEffect(myeffect = DMS_G1G2interactionEffect_4BPmin)
# Spearman's rank correlation rho
# S = 4.2648e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3908235 

## And as a regression?
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

mod1=lme4::lmer(G2methylation~G1methylation * effect + (1|brotherPairID/SampleID), data=AB)
mod2=lme4::lmer(G2methylation~G1methylation + effect + (1|brotherPairID/SampleID), data=AB)
mod3=lme4::lmer(G2methylation~G1methylation + (1|brotherPairID/SampleID), data=AB)

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




