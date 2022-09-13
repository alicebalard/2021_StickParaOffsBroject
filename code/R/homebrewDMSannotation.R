###########################
## Annotate a vector of DMS

# Usage: 
# myHomebrewDMSannotation(DMSvec = c(paste("Gy_chrXVIII", 14108812), paste("Gy_chrXVIII", 14108815), paste("Gy_chrXVIII", 14008815)),
#                         myannotBed12 = annotBed12, myannotGff3 = annotGff3)

myHomebrewDMSannotation <- function(DMSvec, myannotBed12, myannotGff3){
  
  # Change the vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=sapply(strsplit(DMSvec, " "), `[`, 1), 
                                                  start=sapply(strsplit(DMSvec, " "), `[`, 2),
                                                  end=sapply(strsplit(DMSvec, " "), `[`, 2),
                                                  DMS=DMSvec), keep.extra.columns = T)
  A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = myannotBed12)
  ## We assign the feature type to GRangeOBJ
  GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoter",
                                 ifelse(A@members[,2]==1, "exon",
                                        ifelse(A@members[,3]==1, "intron", "intergenic")))
  # Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((A@dist.to.TSS$dist.to.feature>10000 | A@dist.to.TSS$dist.to.feature< -10000) &
                    rowSums(A@members) %in% 0)
  if (is_empty(rows2rm)){
    GRangeOBJ = GRangeOBJ
  } else {
    GRangeOBJ = GRangeOBJ[-rows2rm,]
  }
  
  ## Case 1: the feature is NOT intergenic: we get the annotation by intersection with bed file
  GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]
  names(myannotBed12) = c("exon", "intron", "promoter", "TSSes")
  
  GRangeOBJ1$feature.name = NA
  for (i in 1:length(GRangeOBJ1$featureType)){
    a=myannotBed12[[GRangeOBJ1$featureType[i]]][
      queryHits(GenomicRanges::findOverlaps(myannotBed12[[GRangeOBJ1$featureType[i]]],
                                            GRangeOBJ1[i]))]
    if(length(unique(a$name))>1){
      print("check that these features are identical:"); print(unique(a$name))}
    GRangeOBJ1[i]$feature.name = a$name[1]
  }
  
  ## Case 2: the feature is intergenic: we get the annotation by proximity to nearest TSS
  GRangeOBJ2=GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]
  a = annotateWithGeneParts(target = as(GRangeOBJ2,"GRanges"), feature = annotBed12)
  # Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((a@dist.to.TSS$dist.to.feature>10000 | a@dist.to.TSS$dist.to.feature< -10000) &
                    rowSums(a@members) %in% 0)
  if (is_empty(rows2rm)){  GRangeOBJ2 = GRangeOBJ2
  } else { GRangeOBJ2 = GRangeOBJ2[-rows2rm,] }
  ## Re-annotate the subsetted object
  b = annotateWithGeneParts(as(GRangeOBJ2,"GRanges"),annotBed12)
  ## Get genes associated with these TSS
  c = getAssociationWithTSS(b)
  GRangeOBJ2$feature.name=c$feature.name
  
  ## Merge back 2 cases
  GRangeOBJ=c(GRangeOBJ1, GRangeOBJ2)
  
  ## Get annotations for these genes
  annotDF = merge(GRangeOBJ %>% data.frame()%>% dplyr::rename(chrom = seqnames),
                  data.frame(subset(annotGff3, Name %in% GRangeOBJ$feature.name)) %>%
                    dplyr::select(c("Name", "Note","Ontology_term", "start","end","strand", "Parent"))  %>%
                    dplyr::rename(feature.name=Name, start.gene = start, end.gene = end), by="feature.name")
  
  ## How many CpG per gene?
  annotDF = merge(annotDF,
                  data.frame(table(annotDF$feature.name)) %>% dplyr::rename(feature.name=Var1, nDMSperGene=Freq))
  
  # Add extra info (nbr CpG per gene length, gene length, chrom name)
  annotDF = annotDF  %>%
    mutate(geneLengthkb = (end.gene - start.gene)/1000, nDMSperGenekb = round(nDMSperGene/geneLengthkb,2))
  
  # Add full genes descriptions whenever possible
  # Extract gene symbol from the "Note" attribute
  annotDF$Note = unlist(annotDF$Note)
  annotDF$GeneSymbol = str_extract(annotDF$Note, "(?<=Similar to )(\\w+)")
  # MANUAL CURATION!! All the genes that are weirdly named, with "-" or so
  #check = annotDF[c("Note", "GeneSymbol")]
  listOfWeirdos = c("Type-4 ice-structuring protein LS-12", "MNCb-2990", "Trypsin-3", "anxa2-b", "unc5b-b", "QtsA-11015", # comp1
                    "Type-4 ice-structuring protein LS-12", "MNCb-2990", "en2-a", "draxin-B", "Trypsin-3", " Protein C1orf43 homolog", "Uncharacterized protein FLJ43738", "tlcd4-b", #comp2
                    "tlcd4-b", # comp3
                    "anxa2-b", "QtsA-11015", "unc5b-b") # comp4
  for (i in 1:length(listOfWeirdos)){
    annotDF[grepl(listOfWeirdos[i],annotDF$Note),"GeneSymbol"] = listOfWeirdos[i]
  }
  annotDF$GeneSymbol = annotDF$GeneSymbol %>% toupper # upper case all gene symbols to fit human DB
  # Convert the uniprot gene names to entrez ids
  ENTREZIDlist = mapIds(org.Hs.eg.db, keys = annotDF$GeneSymbol, column = "ENTREZID", keytype = "SYMBOL")
  
  # Retrieve gene summary & description IN HUMANS (more annotation)
  ## Fix if DB too big:
  if (length(ENTREZIDlist) > 300){
    SummaENTREZ1 = entrez_summary(db="gene", id=ENTREZIDlist[1:400])
    SummaENTREZ2 = entrez_summary(db="gene", id=ENTREZIDlist[401:800])
    SummaENTREZ3 = entrez_summary(db="gene", id=ENTREZIDlist[801:length(ENTREZIDlist)])
    SummaENTREZ = c(SummaENTREZ1, SummaENTREZ2,SummaENTREZ3)
  } else {
    SummaENTREZ = entrez_summary(db="gene", id=ENTREZIDlist)
  }
  SummaDF = data.frame(GeneSymbol=sapply(SummaENTREZ, function(x) x[["name"]]) %>% unlist(),
                       ENTREZID = names(SummaENTREZ),
                       description=sapply(SummaENTREZ, function(x) x[["description"]]) %>% unlist(),
                       summary=sapply(SummaENTREZ, function(x) x[["summary"]]) %>% unlist())
  
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
