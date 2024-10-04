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

###########
## GO plots

pdf(file = "../../dataOut/fig/FigS2_GOplot_complete.pdf", width = 30, height = 4)
makeGOplot(dfGO)
dev.off()

## With slim GO terms
pdf(file = "../../dataOut/fig/Fig3C_GOplot_slim.pdf", width = 7, height = 4)
makeGOplotslim(dfGO)$GOplot
dev.off()

dfGOslim=makeGOplotslim(dfGO)$dfGOslim

## Highlight interesting GO terms for our study
listTermsSelect <- unique(c(dfGO[grep("transcription", dfGO$GO.name),"GO.name"],
                            dfGO[grep("expression", dfGO$GO.name),"GO.name"],
                            dfGO[grep("immun", dfGO$GO.name),"GO.name"],
                            dfGO[grep("methyl", dfGO$GO.name),"GO.name"],
                            dfGO[grep("RNA", dfGO$GO.name),"GO.name"]))
makeGOplotslim(dfGO[dfGO$GO.name %in% listTermsSelect,])$GOplot

getGo2Goslim("DNA metabolic process", "BP", "intergenerational")
getGo2Goslim("DNA metabolic process", "BP", "infection-induced")

getGo2Goslim("response to stimulus", "BP", "intergenerational")

getGo2Goslim("RNA metabolic process", "BP", "intergenerational")
getGo2Goslim("RNA metabolic process", "BP", "infection-induced")

getGo2Goslim("catalytic activity", "MF", "infection-induced")
getGo2Goslim("DNA binding", "MF", "intergenerational")



message("R09 done. \n")
