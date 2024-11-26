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

dfGO[dfGO$GO.name %in% "oxidative DNA demethylase activity", ]
EffectsDF_ANNOT[grep("GO:0035516", EffectsDF_ANNOT$Ontology_term),]

## With slim GO terms
dfGOslim = makeGOslim(dfGO = dfGO)

pdf(file = "../../dataOut/fig/Fig3C_GOplot_slim.pdf", width = 7, height = 4)
makeGOslimPlot(dfGOslim)
dev.off()

###########################
## Notes on important terms
getNumGO <- function(list) {
  print("Count infection-induced:")
  print(sum(dfGO[grepl(paste(list, collapse = "|"),
                 dfGO$GO.name) & dfGO$Effect %in% "infection-induced","Count"]))
  print("Size infection-induced:")
  print(sum(dfGO[grepl(paste(list, collapse = "|"),
                 dfGO$GO.name) & dfGO$Effect %in% "infection-induced","Size"]))
  print("Count intergenerational:")
  print(sum(dfGO[grepl(paste(list, collapse = "|"),
                 dfGO$GO.name) & dfGO$Effect %in% "intergenerational","Count"]))
  print("Size intergenerational:")
  print(sum(dfGO[grepl(paste(list, collapse = "|"),
                 dfGO$GO.name) & dfGO$Effect %in% "intergenerational","Size"]))}

dfGOslim[dfGOslim$go_slim_full_name %in% "protein metabolic process", "mapped_go_names"]
getNumGO(c("N-terminal protein myristoylation",
           "N-terminal peptidyl-methionine acetylation",
           "N-terminal protein amino acid modification",
           "positive regulation of protein modification process",
           "positive regulation of cyclin-dependent protein serine/threonine kinase activity",
           "protein ufmylation"))

dfGOslim[dfGOslim$go_slim_full_name %in% "immune system process", "mapped_go_names"]
getNumGO(c("antigen processing and presentation of endogenous antigen",
           "antigen processing and presentation of endogenous peptide antigen via MHC class I",
           "antigen processing and presentation of peptide antigen"))

dfGOslim[dfGOslim$go_slim_full_name %in% "catalytic activity", "mapped_go_names"]
getNumGO(c("galactosylceramide sulfotransferase activity",
"prenylcysteine oxidase activity", "helicase activity", "protein tyrosine phosphatase activity",
"sulfotransferase activity", "procollagen-proline 3-dioxygenase activity",
"ATP-dependent chromatin remodeler activity", "sulfuric ester hydrolase activity",
"protoheme IX farnesyltransferase activity", "N-acyltransferase activity",
"oxidative DNA demethylase activity", "hyaluronan synthase activity", "UFM1 ligase activity",
"fatty-acyl-CoA reductase (alcohol-forming) activity"))

dfGOslim[dfGOslim$go_slim_full_name %in% "signaling receptor activity", "mapped_go_names"]
getNumGO(c("G protein-coupled receptor activity", "parathyroid hormone receptor activity"))

dfGOslim[dfGOslim$go_slim_full_name %in% "protein-containing complex", "mapped_go_names"]
getNumGO(c("transcription factor TFIIIC complex","kinesin complex","EARP complex",
           "NatC complex", "Nem1-Spo7 phosphatase complex", "LUBAC complex"))

dfGOslim[dfGOslim$go_slim_full_name %in% "cellular component organization", "mapped_go_names"]
getNumGO(c("regulation of neuron projection development",
  "regulation of cell projection organization",
  "phospholipid translocation", 
  "regulation of dendrite morphogenesis",
  "positive regulation of axonogenesis",
  "regulation of membrane lipid distribution", 
  "plasma membrane bounded cell projection organization"))

dfGOslim[dfGOslim$go_slim_full_name %in% "developmental process", "mapped_go_names"]
getNumGO(c("positive regulation of cell development",
  "regulation of neuron projection development",
  "neuron differentiation",
  "cell morphogenesis involved in neuron differentiation",
  "regulation of dendrite morphogenesis",
  "regulation of neurogenesis",
  "positive regulation of axonogenesis",
  "positive regulation of nervous system development"))

dfGOslim[dfGOslim$go_slim_full_name %in% "transporter activity", "mapped_go_names"]
getNumGO(c("phosphate ion transmembrane transporter activity",
  "phospholipid transporter activity", "cation transmembrane transporter activity",
  "calcium ion transmembrane transporter activity",
  "antiporter activity","P-type ion transporter activity",
  "large conductance calcium-activated potassium channel activity", "intramembrane lipid transporter activity")
  )


########

dfGOslim[dfGOslim$go_slim_full_name %in% "cell differentiation", "mapped_go_names"]

dfGO[grepl(paste(c("positive regulation of cell development",
                   "regulation of neuron projection development",
                   "neuron differentiation",
                   "cell morphogenesis involved in neuron differentiation",
                   "regulation of dendrite morphogenesis",
                   "regulation of neurogenesis",
                   "positive regulation of axonogenesis"), collapse = "|"),
           dfGO$GO.name),"Size"] %>% sum
dfGO[grepl(paste(c("positive regulation of cell development",
                   "regulation of neuron projection development",
                   "neuron differentiation",
                   "cell morphogenesis involved in neuron differentiation",
                   "regulation of dendrite morphogenesis",
                   "regulation of neurogenesis",
                   "positive regulation of axonogenesis"), collapse = "|"), 
           dfGO$GO.name),"Count"] %>% sum



# infection-induced: antigen processing and presentation of endogenous antigen;antigen processing and presentation of endogenous peptide antigen via MHC class I;antigen processing and presentation of peptide antigen
dfGO[grep("antigen", dfGO$GO.name), ]



dfGOslim[dfGOslim$go_slim_full_name %in% "immune system process", ]
# infection-induced: antigen processing and presentation of endogenous antigen;antigen processing and presentation of endogenous peptide antigen via MHC class I;antigen processing and presentation of peptide antigen
dfGO[grep("antigen", dfGO$GO.name), ]

dfGOslim[dfGOslim$go_slim_full_name %in% "DNA binding", "mapped_go_names"]
dfGO[grep("RNA polymerase II cis-regulatory region sequence-specific DNA binding", dfGO$GO.name), ]
dfGO[grep("ATP-dependent DNA/DNA annealing activity", dfGO$GO.name), ]

dfGOslim[dfGOslim$go_slim_full_name %in% "cell differentiation", "mapped_go_names"]

dfGO[grep("positive regulation of cell development", dfGO$GO.name), "Count"]+
  dfGO[grep("regulation of neuron projection development", dfGO$GO.name), "Count"]+
  sum(dfGO[grep("neuron differentiation", dfGO$GO.name), "Count"])+
  dfGO[grep("regulation of dendrite morphogenesis", dfGO$GO.name), "Count"]+
  dfGO[grep("regulation of neurogenesis", dfGO$GO.name), "Count"]+
  dfGO[grep("positive regulation of axonogenesis", dfGO$GO.name), "Count"]/
  
  (dfGO[grep("positive regulation of cell development", dfGO$GO.name), "Size"]+
     dfGO[grep("regulation of neuron projection development", dfGO$GO.name), "Size"]+
     sum(dfGO[grep("neuron differentiation", dfGO$GO.name), "Size"])+
     dfGO[grep("regulation of dendrite morphogenesis", dfGO$GO.name), "Size"]+
     dfGO[grep("regulation of neurogenesis", dfGO$GO.name), "Size"]+
     dfGO[grep("positive regulation of axonogenesis", dfGO$GO.name), "Size"])
## 13% (14/126)

dfGOslim[grep("transc", dfGOslim$go_slim_full_name),]
dfGOslim[grep("transc", dfGOslim$mapped_go_names),]


dfGOslim[dfGOslim$go_slim_full_name %in% "DNA binding", ]
# intergenerational: "RNA polymerase II cis-regulatory region sequence-specific DNA binding", "ATP-dependent DNA/DNA annealing activity"

dfGOslim[dfGOslim$go_slim_full_name %in% "cell differentiation", ]
# intergenerational: "positive regulation of cell development", "regulation of neuron projection development", "neuron differentiation", "cell morphogenesis involved in neuron differentiation", "regulation of dendrite morphogenesis", "regulation of neurogenesis", "positive regulation of axonogenesis"

dfGOslim[dfGOslim$go_slim_full_name %in% "RNA metabolic process", ]
# infection-induced: transcription initiation at RNA polymerase III promoter;oxidative single-stranded RNA demethylation;RNA repair RNA metabolic process
# intergenerational: transcription initiation at RNA polymerase III promoter;positive regulation of transcription by RNA polymerase II RNA metabolic process

dfGOslim[dfGOslim$go_slim_full_name %in% "DNA metabolic process", ]
# infection-induced: "DNA dealkylation involved in DNA repair", "DNA integration", "oxidative single-stranded DNA demethylation"
# intergenerational: "transcription-coupled nucleotide-excision repair", "replication fork processing"

dfGOslim[dfGOslim$go_slim_full_name %in% "response to stimulus", ]
# infection-induced: "DNA dealkylation involved in DNA repair", "G protein-coupled receptor signaling pathway"
# intergenerational: "transcription-coupled nucleotide-excision repair", "integrin-mediated signaling pathway", "regulation of Ral protein signal transduction", "negative regulation of TORC1 signaling"

dfGOslim[dfGOslim$go_slim_full_name %in% "protein-containing complex", ]
# infection-induced: "transcription factor TFIIIC complex", "NatC complex", "Nem1-Spo7 phosphatase complex", "LUBAC complex protein-containing complex"
# intergenerational: "transcription factor TFIIIC complex", "kinesin complex", "EARP complex"

dfGOslim[dfGOslim$go_slim_full_name %in% "catalytic activity", ]
# infection-induced: "sulfuric ester hydrolase activity", "protoheme IX farnesyltransferase activity", "N-acyltransferase activity", 
# "oxidative DNA demethylase activity", "hyaluronan synthase activity", "UFM1 ligase activity", "fatty-acyl-CoA reductase (alcohol-forming) activity"
# intergenerational: "galactosylceramide sulfotransferase activity", "prenylcysteine oxidase activity", 
# "helicase activity", "protein tyrosine phosphatase activity", "sulfotransferase activity", "procollagen-proline 3-dioxygenase activity",
# "ATP-dependent chromatin remodeler activity"

message("R09 done. \n")
