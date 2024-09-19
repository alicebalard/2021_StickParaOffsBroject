# Each script sources the previous script of the pipeline if needed
source("R08_DMSannotation.R")

message("R09 starting...\n")

#################################################
### Gene Ontology analysis OVERALL then by effect

# create gene universe from all the covered CpGs
gene_universe <- data.frame(
  subsetByOverlaps(GRanges(annotGff3), 
                   GRanges(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP))) %>% # subselect covered CpGs
  filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
  filter(type %in% "gene")  %>% # keep all the 7404 genes with GO terms
  dplyr::select(c("Name", "Ontology_term")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, it's from Interproscan after Maker. Sticklebacks (biomart) have 82701 IEA and 63 ISS.
  relocate("Ontology_term","go_linkage_type","Name") %>%
  unnest(Ontology_term) %>% # one GO per line (was a list before in this column)
  data.frame()

gene_universe$Name %>% unique %>% length #9234 genes (before re-annot: 7404)

# Create gene set collection
goFrame <- GOFrame(gene_universe, organism="Gasterosteus aculeatus")
goAllFrame <- GOAllFrame(goFrame)
gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())

# **IMPORTANT NOTE from Mel: why conditional hypergeometric test?** 
# The GO ontology is set up as a directed acyclic graph, where a parent term is 
# comprised of all its child terms. If you do a standard hypergeometric, you might e.g., 
# find 'positive regulation of kinase activity' to be significant. If you then test 'positive 
# regulation of catalytic activity', which is a parent term, then it might be significant
# as well, but only because of the terms coming from positive regulation of kinase activity.
# The conditional hypergeometric takes this into account, and only uses those terms 
# that were not already significant when testing a higher order (parent) term.

## Run GO on genes (use "unique") level (same if done on CpG level, the method
## considers one occurence of each)
dfGO_allDMS = makedfGO(
  annot = EffectsDF_ANNOT %>% distinct(feature.name,.keep_all = TRUE), 
  gene_universe = gene_universe,
  effect = "allDMS",
  label = paste0(length(unique(EffectsDF_ANNOT$feature.name)), " genes with DMS"))

dfGO_intergen = makedfGO(
  annot = EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL",] %>%
    distinct(feature.name,.keep_all = TRUE), 
  gene_universe = gene_universe, 
  effect = "intergenerational",
  label = paste0(length(unique(
    EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL","feature.name"])),
    " genes with DMS associated with intergenerational effect only"))

dfGO_infectind = makedfGO(
  annot = EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED",] %>%
    distinct(feature.name,.keep_all = TRUE), 
  gene_universe = gene_universe, 
  effect = "infection-induced",
  label = paste0(length(unique(
    EffectsDF_ANNOT[EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED","feature.name"])),
    " genes with DMS associated with infection-induced effect only"))

dfGO = rbind(dfGO_intergen, dfGO_infectind)

# ###############
# ## Check GO slim terms for easy interpretation 9not very meaningful...)
# ## GO subsets (also known as GO slims) are condensed versions of the GO containing a subset of the terms. 
# # dl the GO slim Developed by GO Consortium for the Alliance of Genomes Resources
# # download.file(url = "https://current.geneontology.org/ontology/subsets/goslim_agr.obo",
#               # destfile = "../../data/goslim_agr.obo")
# 
# slim <- GSEABase::getOBOCollection("../../data/goslim_agr.obo")
# 
# GSEABase::goSlim(idSrc = GOCollection(dfGO$GO.term),
#                  slimCollection = slim, 
#                  ontology = "BP") %>%  filter(Count !=0)
# 
# GSEABase::goSlim(idSrc = GOCollection(dfGO$GO.term),
#                                slimCollection = slim, 
#                                ontology = "MF") %>%  filter(Count !=0)
# 
# GSEABase::goSlim(idSrc = GOCollection(dfGO$GO.term),
#                                slimCollection = slim, 
#                                ontology = "CC") %>%  filter(Count !=0)

## Select interesting terms for our study
listTermsSelect <- unique(c(dfGO[grep("transcription", dfGO$GO.name),"GO.name"],
                            dfGO[grep("expression", dfGO$GO.name),"GO.name"],
                            dfGO[grep("immun", dfGO$GO.name),"GO.name"],
                            dfGO[grep("methyl", dfGO$GO.name),"GO.name"],
                            dfGO[grep("RNA", dfGO$GO.name),"GO.name"]))

###########
## GO plots
makeGOplot <- function(dfGO){
  dfGO %>%
  dplyr::filter(p.value.adjusted < 0.05) %>% 
  ggplot(aes(x=Effect, y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(
    name="adjusted\np-value", low = "red", high = "blue", 
    limits = c(0, 0.05), breaks = c(0, 0.02, 0.04), labels =c("0", "0.02", "0.04")) +
  scale_size_continuous(name = "% of genes")+
  theme_bw() + ylab("") + xlab("") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), # grey box for legend
        legend.position="top",
        axis.text.y = element_text(size = 8),  # Decrease y-axis text size
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1)  # Increase x-axis text size and rotate
  )+
  facet_grid(.~fct_inorder(GO.category), scales="free",space = "free")+
  coord_flip() + # flip axes
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+ # split long text
  scale_y_discrete(limits=rev, # revers axis to have alphabetical order
                   labels = function(x) str_wrap(x, width = 30)) # split too long GO names in half
}

pdf(GOplot, file = "../../dataOut/fig/FigS3_GOplot_complete.pdf", width = 30, height = 4)
makeGOplot(dfGO)
dev.off()

pdf(GOplot, file = "../../dataOut/fig/Fig3C_GOplot_subset.pdf", width = 7, height = 4)
makeGOplot(dfGO[dfGO$GO.name %in% listTermsSelect,] )
dev.off()

message("R09 done. \n")
