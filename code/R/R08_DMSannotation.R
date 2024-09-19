# Each script sources the previous script of the pipeline if needed
source("R07_DMSeffects.R")

## Produces: 
## dataOut/fig/Fig3A_DMSgroupsVenn_geneLevel.pdf
## dataOut/fig/Fig3B_positionsAnnotatedGenes.pdf

## Annotate the different DMS groups:

## NB: set rerun = T when change list to reload (>30min!)
EffectsDF_ANNOT= myHomebrewDMSannotation(
  DMSvec = EffectsDF$pos, myannotBed12 = annotBed12, myannotGff3 = annotGff3, rerun = FALSE)

# NB: Some positions have the same GeneName but different GeneSymbol, if they were annotated from different species
# (namely Rho guanine nucleotide exchange factor 28, Collagen alpha-1, Metabotropic glutamate receptor 7, Titin)

# we have annotations for 647 DMS on 548 unique known genes, and 60 DMS on genes coding for proteins with unknown function (Total= 707 DMS).
length(unique(EffectsDF_ANNOT$GeneName)) # 548 known genes
length(unique(EffectsDF_ANNOT$uniprotID)) # 556 known protein
length(unique(EffectsDF_ANNOT$Function)) # 493 functions known

## Add the effect detected in the previous script by merging
EffectsDF_ANNOT$pos = EffectsDF_ANNOT$DMS
EffectsDF_ANNOT = merge(EffectsDF_ANNOT, EffectsDF, by = c("pos"))

# Plot a Venn diagram to see genes in common
pdf(file = "../../dataOut/fig/Fig3A_DMSgroupsVenn_geneLevel.pdf", width = 7, height = 6)
ggVennDiagram(list(
  "Intergenerational" = EffectsDF_ANNOT$feature.name[EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL"],
  "Infection-induced" = EffectsDF_ANNOT$feature.name[EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED"],
  "Additive" = EffectsDF_ANNOT$feature.name[EffectsDF_ANNOT$effect %in% "ADDITIVE"],
  "Interaction" = EffectsDF_ANNOT$feature.name[EffectsDF_ANNOT$effect %in% "INTERACTION"]),
  label_alpha = 0)+
  scale_color_manual(values = c(1,1,1,1))+
  scale_fill_gradient(low="white",high = "red") + 
  theme(legend.position = "none") + ggtitle("Genes associated with each effect")
dev.off()

## Look for keywords
EffectsDF_ANNOT[grep("immune", EffectsDF_ANNOT$Function, ignore.case = T),
                c("uniprotID")] %>% table # P01730 has 6 DMS
EffectsDF_ANNOT[EffectsDF_ANNOT$uniprotID %in% "P01730",] # CD4 T-cell gene!! effect INTERGENERATIONAL

EffectsDF_ANNOT[grep("transcription", EffectsDF_ANNOT$Function, ignore.case = T),
                c("uniprotID")] %>% table %>% data.frame() %>% arrange(Freq) # D3ZUU2 and O19006 have 4 DMS, A2A8L1 3 DMS
EffectsDF_ANNOT[EffectsDF_ANNOT$uniprotID %in% "D3ZUU2",] # Transcriptional repressor, effect INFECTION_INDUCED
EffectsDF_ANNOT[EffectsDF_ANNOT$uniprotID %in% "O19006",] # bmp2 gene, Growth factor, transcriptional activity, effect INFECTION_INDUCED
EffectsDF_ANNOT[EffectsDF_ANNOT$uniprotID %in% "A2A8L1",] # Chd5 gene, Chromatin-remodeling protein that binds DNA through histones and regulates gene transcription, effect INTERGENERATIONAL

EffectsDF_ANNOT[grep("MAPK", EffectsDF_ANNOT$Function, ignore.case = T),
                c("uniprotID")] %>% table %>% data.frame() %>% arrange(Freq) # P35400 and B0LT89 have 3 DMS
EffectsDF_ANNOT[EffectsDF_ANNOT$uniprotID %in% "B0LT89",] # Chd5 gene, Chromatin-remodeling protein that binds DNA through histones and regulates gene transcription, effect INTERGENERATIONAL

## NB 20 genes have DMS in different effects!
my20genes <- intersect(EffectsDF_ANNOT$feature.name[EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL"],
                       EffectsDF_ANNOT$feature.name[EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED"])

table(EffectsDF_ANNOT[EffectsDF_ANNOT$feature.name %in% my20genes,"uniprotID"])
# A2A8L5 B0LT89 B5DE93 E7F5E1 O02747 O14763 P35400 P43141 Q6DFV8 Q6KEQ9 Q70EL4 
# 3      3      2      2      2      2      2      2      2      2      2 
# Q8BMQ2 Q8BPM0 Q8R5G7 Q99K41 Q9JI18 
# 2      2      2      2      2 

df20=EffectsDF_ANNOT[EffectsDF_ANNOT$feature.name %in% my20genes,
                     c("uniprotID", "GeneSymbol", "Species", "Function")]

df20[grep("immune", df20$Function, ignore.case = T),] # Stk24 x3, AHRx2
df20[grep("transcription", df20$Function, ignore.case = T),] # cdk12 x2, AHRx2
df20[grep("MAPK", df20$Function, ignore.case = T),] # Stk24 x3, GRM7x2
df20[grep("rho", df20$Function, ignore.case = T),] # Stk24 x3, Arap3x2, DAAM1x2
# Stk24 immune regulation
# AHR transcription factor, immunity
# cdk12 key regulator of transcription elongation. Regulates the expression of genes 
# involved in DNA repair and is required for the maintenance of genomic stability.


# add summary in select
EffectsDF_ANNOT = EffectsDF_ANNOT %>% 
  dplyr::mutate(nDMSperGeneKb = EffectsDF_ANNOT$nDMSperGene / EffectsDF_ANNOT$geneLengthkb) %>%
  arrange(desc(nDMSperGeneKb)) %>% #arrange in descending order
  dplyr::select(c(names(EffectsDF_ANNOT)[!names(EffectsDF_ANNOT) %in% "Function"], "Function")) # reorder columns with summary at the end

# Write out
write.csv(EffectsDF_ANNOT, file = "../../dataOut/fig/TableS1_EffectsDF_ANNOT.csv", row.names = F)

#### Focus on gene CD4 & TRIM16
plotGeneTarget <- function(myTargetGene){
  dfplot = EffectsDF_ANNOT[EffectsDF_ANNOT$GeneSymbol %in% myTargetGene,]
  # Find TSS position of the gene
  dfplot$TSSpos = annotBed12$TSSes[grep(unique(dfplot["feature.name"]),
                                        annotBed12$TSSes$name)]@ranges@start
  # Set TSS as origin
  dfplot$start_distToTSS = dfplot$start - dfplot$TSSpos 
  dfplot$end_distToTSS = dfplot$end - dfplot$TSSpos 
  dfplot$start.gene_distToTSS = dfplot$start.gene - dfplot$TSSpos 
  dfplot$end.gene_distToTSS = dfplot$end.gene - dfplot$TSSpos 
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
    geom_segment(data = data.frame(a=0, b=0),
                 aes(x = 0, xend = 0, y=0, yend=.5), col = "red", size = 3) + # TSS
    theme_blank() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    labs(title = paste(unique(dfplot$GeneSymbol), ":", dfplot$description),
         subtitle = str_wrap(dfplot$summary, width = 150))
  plotGeneTarget
}
  
plotGeneTarget("CD4")
plotGeneTarget("TRIM16")
plotGeneTarget("bmp2")

#### Manhattan plot of the genes in the four main effects
## Prepare data and change gene position to start at the good chromosome
data4Manhattan = dplyr::left_join(EffectsDF_ANNOT, GYgynogff) %>%
  dplyr::mutate(posInPlot=((end.gene+start.gene)/2)+gstart)

data4Manhattan =
  data4Manhattan %>% dplyr::select(type,gstart,gend,gmid,chrom,nDMSperGene,
                                   posInPlot,nDMSperGenekb,GeneSymbol,effect) %>%
  unique
# level effects in order
# data4Manhattan$effect <- factor(data4Manhattan$effect,levels = c("G1", "G2", "addit", "inter"))

# Manhattan plot
ggplot()+
  # add grey background every second chromosome
  geom_rect(data=data4Manhattan[data4Manhattan$type %in% "B",],
            aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf), fill=rgb(.9,.9,.9))+
  xlab(NULL)+
  scale_x_continuous(name = "Chromosomes", 
                     breaks=data4Manhattan[!duplicated(data4Manhattan$chrom),"gmid"],
                     labels=data4Manhattan[!duplicated(data4Manhattan$chrom),"chrom"] %>% str_remove(.,"Gy_chr"),
                     position = "bottom",expand = c(0,0))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ # add frame
  ylab("Number of differentially methylated CpG per gene kb")+
  scale_y_continuous(breaks = 0:5)+
  geom_point(data = data4Manhattan, aes(x=posInPlot, y = nDMSperGenekb)) +
  geom_label_repel(data = data4Manhattan[data4Manhattan$nDMSperGenekb >1,], 
                   aes(x=posInPlot, y = nDMSperGenekb, label = GeneSymbol, fill = factor(effect)),
                   color="white",max.overlaps = Inf, box.padding = 0.6, segment.color = "grey")+
  scale_fill_manual(values=colEffects)

#########################################################
## Better than the Manhattan plot: plot along chromosomes
unique=data.frame(chrom=unique(EffectsDF$chrom),
                  length=unique(EffectsDF$length))

EffectsDF_ANNOT = EffectsDF_ANNOT %>% 
  mutate(GeneSymbol = ifelse(!is.na(GeneSymbol), GeneSymbol, "unknown"))

message("save plot of the DMS position along the chromosomes in dataOut/fig/Fig3C_positionsAnnotatedGenes.pdf")
pdf("../../dataOut/fig/Fig3B_positionsAnnotatedGenes.pdf", width = 10, height = 6)

ggplot(EffectsDF_ANNOT) +
  geom_col(data = unique, aes(x=length, y=chrom), width = .7, fill="#e0ebeb")+ # plot full chromosome
  geom_tile(aes(x=start, y=chrom, fill=effect, width = 100000, height = .8))+
  geom_label_repel(data = EffectsDF_ANNOT[!duplicated(EffectsDF_ANNOT$feature.name) &
                                            EffectsDF_ANNOT$nDMSperGene >2,],
                   aes(x=start, y = chrom, 
                       label = paste0(GeneSymbol, " (", nDMSperGene, "DMS)"),
                       fill = factor(effect)),
                   color="white",max.overlaps = Inf, box.padding = 0.6, segment.color = "grey") +
  theme_blank()+
  scale_fill_manual(values = as.vector(palette.colors(palette = "Okabe-Ito")[1:4]),
                    name = "Effect:", labels = c("infection induced", "intergenerational", "additive", "interaction"))+ # rename legend
  scale_x_continuous(breaks = seq(0, 3e+07, by = 0.5e+7),
                     labels = paste0(seq(0, 3e+07, by = 0.5e+7)/1e+6, "Mb"),
                     expand = c(0, 0)) + # Remove space before 0
   ylab(NULL) + xlab(NULL)+
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.85, 0.75))
dev.off()

message("R08 done.\n")

