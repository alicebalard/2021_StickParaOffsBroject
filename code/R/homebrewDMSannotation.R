###########################
## Annotate a vector of DMS

# Usage: 
# myHomebrewDMSannotation(DMSvec = c("Gy_chrI_26116565", "Gy_chrII_2635914"),
#                         myannotBed12 = annotBed12, myannotGff3 = annotGff3)

myHomebrewDMSannotation <- function(DMSvec, myannotBed12, myannotGff3){
  # Change the vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=paste(sapply(strsplit(DMSvec, "_"), `[`, 1), 
                                                            sapply(strsplit(DMSvec, "_"), `[`, 2), sep = "_"),
                                                  start=sapply(strsplit(DMSvec, "_"), `[`, 3),
                                                  end=sapply(strsplit(DMSvec, "_"), `[`, 3),
                                                  DMS=DMSvec), keep.extra.columns = T)
  
  # Annotate with features
  GRangeOBJ_annot = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), 
                                          feature = myannotBed12)

  ## We assign the feature type to GRangeOBJ
  GRangeOBJ$featureType = ifelse(GRangeOBJ_annot@members[,1]==1, "promoters",
                                 ifelse(GRangeOBJ_annot@members[,2]==1, "exons",
                                        ifelse(GRangeOBJ_annot@members[,3]==1, "introns", "intergenic")))
  ## We assign distance to TSS
  GRangeOBJ$dist.to.TSS = GRangeOBJ_annot@dist.to.TSS$dist.to.feature
  
  ## valid only if trusted, to check
  GRangeOBJ$feature.name = as.character(GRangeOBJ_annot@dist.to.TSS$feature.name)
  
  # Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((GRangeOBJ$dist.to.TSS>10000 | 
                     GRangeOBJ$dist.to.TSS < -10000) &
                    GRangeOBJ$featureType %in% "intergenic")
  if (is_empty(rows2rm)){
    GRangeOBJ = GRangeOBJ
  } else {
    GRangeOBJ = GRangeOBJ[-rows2rm,]
  }
  
  # Change recursively the gene names to keep only ID
  getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
  
  for (i in 1:length(GRangeOBJ)){
    GRangeOBJ$feature.name[i] <- getName(GRangeOBJ$feature.name[i])
  }
  
  ## Get annotations for our DMS
  annotDF = data.frame(GRangeOBJ)
  
  annotDF$Note = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Note"] %>% unlist()
  annotDF$Ontology_term = sapply(myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Ontology_term"], 
                                 function(x) paste(x, collapse = " ")) ## keep multiple GO terms together
  annotDF$start = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"start"] 
  annotDF$end = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"end"] 
  annotDF$strand = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"strand"] 
  annotDF$Parent = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Parent"] %>% unlist()
  
  annotDF = annotDF %>% dplyr::rename(start.gene = start, end.gene = end)
  annotDF = annotDF[!names(annotDF) %in% "width"]

  ## How many CpG per gene?
  annotDF = merge(annotDF,
                  data.frame(table(annotDF$feature.name)) %>% 
                    dplyr::rename(feature.name=Var1, nDMSperGene=Freq))
  
  # Add extra info
  annotDF = annotDF  %>%
    mutate(geneLengthkb = (end.gene - start.gene)/1000,
           nDMSperGenekb = nDMSperGene/geneLengthkb,
           GeneSymbol =   str_extract(Note, "(?<=Similar to ).*?(?=[:()])")) # extract after Similar to and before : or (
  
  annotDF$GeneSymbol = annotDF$GeneSymbol %>% toupper # upper case all gene symbols to fit human DB
  
  message("we have ",   
          length(unique(annotDF$feature.name[!annotDF$Note %in% "Protein of unknown function"])),
          " DMS on ", length(unique(annotDF$GeneSymbol[!is.na(annotDF$GeneSymbol)])), " unique known genes",
          " and ", table(is.na(annotDF$GeneSymbol))[2], " DMS on ",  
          length(unique(annotDF$feature.name[annotDF$Note %in% "Protein of unknown function"])), 
          " genes with unknown function")
  
  # Convert the uniprot gene names to entrez ids
  ENTREZIDlist = mapIds(org.Hs.eg.db, keys = unique(na.omit(annotDF$GeneSymbol)),
                        column = "ENTREZID", keytype = "SYMBOL")
  
  ##### MANUAL CURATION AREA
  message("Check NA in ENTREZIDlist and manually add it! To do after annotation!!!")
  message(paste("There are ",
                nrow(annotDF[is.na(annotDF$GeneSymbol) & !annotDF$Note %in% "Protein of unknown function",]),
                " genes to curate manually"))
  
  # change of gene ID
  annotDF[annotDF$GeneSymbol %in% "RNF165", "GeneSymbol"] = "ARK2C" 
  
  ENTREZIDlist = ENTREZIDlist[!is.na(ENTREZIDlist)]

  # Retrieve gene summary & description IN HUMANS (more annotation)
  
  # Throw in an error if too toooo big (needs then a bigger fix e.g. rollaply)
  if (length(ENTREZIDlist) > 1200){print("ERROR nCpG too numerous, please keep it below ~1200 genes")}
  ## Fix if DB too big:
  if (length(ENTREZIDlist) > 300 & length(ENTREZIDlist) < 800){
    SummaENTREZ1 = entrez_summary(db="gene", id=ENTREZIDlist[1:400])
    SummaENTREZ2 = entrez_summary(db="gene", id=ENTREZIDlist[401:length(ENTREZIDlist)])
    SummaENTREZ = c(SummaENTREZ1, SummaENTREZ2)
  } else if (length(ENTREZIDlist) > 300 & length(ENTREZIDlist) > 800){
    SummaENTREZ1 = entrez_summary(db="gene", id=ENTREZIDlist[1:400])
    SummaENTREZ2 = entrez_summary(db="gene", id=ENTREZIDlist[401:800])
    SummaENTREZ3 = entrez_summary(db="gene", id=ENTREZIDlist[801:length(ENTREZIDlist)])
    SummaENTREZ = c(SummaENTREZ1, SummaENTREZ2,SummaENTREZ3)
  } else {
    SummaENTREZ = entrez_summary(db="gene", id=ENTREZIDlist)
  }
  
  ## Select useful
  SummaENTREZ = SummaENTREZ[names(SummaENTREZ) %in% ENTREZIDlist]
  
  SummaDF = unique(data.frame(GeneSymbol=  sapply(SummaENTREZ, function(x) x[["name"]]),
                       ENTREZID = names(SummaENTREZ),
                       description=sapply(SummaENTREZ, function(x) x[["description"]]),
                       summary=sapply(SummaENTREZ, function(x) x[["summary"]])))
    
  # merge with Notes from uniprot (contained in the gff3)
  SummaDF = unique(merge(annotDF, SummaDF, all=T))
  
  # Order by nDMSperGenekb
  SummaDF = SummaDF[order(SummaDF$nDMSperGenekb, decreasing = T),]
  # set rownames NULL
  rownames(SummaDF) = NULL
  # unlist "Note"
  SummaDF$Note=unlist(SummaDF$Note)
  # remove NA and character(0) columns
  SummaDF = SummaDF[apply(SummaDF, 2, function(x) length(unlist(x)))!=0]
  SummaDF = SummaDF[apply(SummaDF, 2, function(x) sum(!is.na(x))!=0)]
  return(SummaDF)
}

                                                      