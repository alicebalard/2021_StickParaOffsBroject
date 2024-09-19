#####################################################################
## Annotate a vector of DMS using a bed12 and a gff3 annotation files
#####################################################################

# Usage:
# myHomebrewDMSannotation(DMSvec = c("Gy_chrI_26116565", "Gy_chrII_2635914"),
#                         myannotBed12 = annotBed12, myannotGff3 = annotGff3)

## Issue in the original "annotateWithGeneParts" from genomation function:
## involved associating all CpGs by proximity to a TSS first, then assigning
## gene ID/feature type. However, if a CpG is at the end of a large gene, it could
## be closer to the TSS of the next gene along, so would get associated incorrectly.

## Solution: 2 cases, one for the genic part, overlapping directly with the annotation
## via intersection; then the TSS proximity method is only applied for intergenic sites.

## Charley Yen's (https://github.com/eugeniecyen) improvements on my original function:
## - removing a filter by distance <10kb from TSS before splitting into the 2 cases,
## as this might also remove some genic sites on large genes.
## - having to set unique.prom=FALSE for promoters to come up in readTranscriptFeatures
## to load bed12 annotation.
## - for some reason, the gene ID association to DMS was getting jumbled up in the
## final output table. This seemed to be occurring with the old looping code, but
## no longer with the adapted homebrew function

# set rerun = T when change list to reload
myHomebrewDMSannotation <- function(DMSvec, myannotBed12, myannotGff3, rerun=F){
  myGRanges = getGRange_corrected(DMSvec, myannotBed12)
  message(paste0("We have ", length(myGRanges), " positions extracted from the bed12 file."))
  if (rerun == T){
    message("We rerun the annotation, it will take some time...")  
  } else {
    message("annotation is simply reloaded.")
  }
  res = getAnnotDMS(myGRanges,myannotGff3, rerun=rerun)
  return(res)
}

getGRange_corrected <- function(DMSvec, myannotBed12){
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
  
  ## CASE 1: genic --> intersection with BED12 file
  GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]
  
  # Add empty column called geneInfo to fill with gene IDs
  GRangeOBJ1$feature.name <- NA
  
  add_geneInfo_genic <- function(x, GRangeOBJ, annotBed12){
    ov = GenomicRanges::findOverlaps(
      annotBed12[[x]],
      GRangeOBJ[GRangeOBJ$featureType %in% x,])
    ## Add gene annotation to subject GRanges (i.e. left join)
    mcols(GRangeOBJ[GRangeOBJ$featureType %in% x,])[subjectHits(ov), "feature.name"] =
      mcols(annotBed12[[x]])[queryHits(ov), "name"]
    return(GRangeOBJ)
  }
  
  GRangeOBJ_ex = add_geneInfo_genic("exons", GRangeOBJ1, myannotBed12) # Add gene names for sites on exons
  GRangeOBJ_ex_in = add_geneInfo_genic("introns", GRangeOBJ_ex, myannotBed12) # Add gene names for sites on introns
  GRangeOBJ_genic = add_geneInfo_genic("promoters", GRangeOBJ_ex_in, myannotBed12) # Add gene names for sites on promoters
  
  ## CASE 2: intergenic --> original method
  GRangeOBJ2=GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]
  
  ### Extract annotation from the overall annotation object
  GRangeOBJ_annot_2 = GRangeOBJ_annot@dist.to.TSS[GRangeOBJ$featureType %in% "intergenic",]
  
  GRangeOBJ_annot_2$feature.name
  
  ## Filter out the DMS non associated with a gene:
  # Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((GRangeOBJ2$dist.to.TSS>10000 |
                     GRangeOBJ2$dist.to.TSS < -10000) &
                    GRangeOBJ2$featureType %in% "intergenic")
  if (!is_empty(rows2rm)){
    GRangeOBJ2 = GRangeOBJ2[-rows2rm,]
    GRangeOBJ_annot_2 = GRangeOBJ_annot_2[-rows2rm,]
  }
  
  GRangeOBJ_intergenic = GRangeOBJ2
  GRangeOBJ_intergenic$feature.name = GRangeOBJ_annot_2$feature.name
  
  ## Merge both cases
  GRangeOBJ_both = c(GRangeOBJ_genic, GRangeOBJ_intergenic)
  
  # Change recursively the gene names to keep only ID
  getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
  
  for (i in 1:length(GRangeOBJ_both)){
    GRangeOBJ_both$feature.name[i] <- getName(GRangeOBJ_both$feature.name[i])
  }
  
  return(GRangeOBJ_both)
}

# set rerun = T when change list to reload
getAnnotDMS <- function(GRangeOBJ, myannotGff3, rerun){ 
  
  ## Get annotations for our DMS
  annotDF = data.frame(GRangeOBJ)
  annotDF = annotDF[!names(annotDF) %in% "width"]
  
  ### Extract info from gff3 by overlap
  annotDF$GeneName = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Name"]
  annotDF$Note = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Note"] %>% unlist()
  annotDF$Ontology_term = sapply(myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Ontology_term"],
                                 function(x) paste(x, collapse = " ")) ## keep multiple GO terms together
  annotDF$start = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"start"]
  annotDF$end = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"end"]
  annotDF$strand = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"strand"]
  annotDF$Parent = myannotGff3[match(annotDF$feature.name,myannotGff3$Name),"Parent"] %>% unlist()
  
  ### Calculate number of DMS per gene
  annotDF = merge(annotDF, data.frame(table(annotDF$feature.name)) %>%
                    dplyr::rename(feature.name=Var1, nDMSperGene=Freq)) %>%
    dplyr::rename(start.gene = start, end.gene = end) %>%
    dplyr::mutate(geneLengthkb = (end.gene - start.gene)/1000,
                  nDMSperGenekb = nDMSperGene/geneLengthkb)
  
  ### Extract info from note
  annotDF = annotDF  %>%
    mutate(
      GeneName = trimws(stringr::str_extract(Note, "(?<=:).*?(?=\\()")), # extract after Similar to and before : or (
      GeneSymbol = stringr::str_extract(Note, "(?<=Similar to ).*?(?=[:()])"), # extract after Similar to and before : or (
      Species = stringr::str_extract(Note, "(?<=\\()[^()]*(?= OX)"))# extract the organism
  
  temp = annotDF[!is.na(annotDF$GeneName),c("GeneName", "GeneSymbol")] %>%
    arrange(GeneSymbol) %>% unique
  message("NB: Some positions have the same GeneName but different GeneSymbol, if they were annotated from different species
          (namely ", paste(unique(temp$GeneName[duplicated(temp$GeneName)]), collapse = ', '), ")")
  rm(temp)
  
  message(
    "we have annotations for ",
    length(annotDF$GeneName[!annotDF$Note %in% "Protein of unknown function"]),
    " DMS on ",
    length(unique(annotDF$GeneName[!annotDF$Note %in% "Protein of unknown function"])),
    " unique known genes, and ",
    length(annotDF$GeneName[annotDF$Note %in% "Protein of unknown function"]),
    " DMS on genes coding for proteins with unknown function",
    " (Total= ", length(annotDF$GeneName), " DMS).")
  ## we have annotations for 647 DMS on 548 unique known genes,
  ## and 60 DMS on genes coding for proteins with unknown function (Total= 707 DMS).
  
  message(paste0("rerun=", rerun,
                 " : rerun the search of protein function, new list of genes (takes >30min minutes)"))
  
  if (rerun == T){
    ## Get protein functions from Uniprot
    ### Keys for all species in uniprot
    uniprotOrgs = UniProt.ws::availableUniprotSpecies()
    
    ### Function to retrieve the description of the protein function
    findFunctionProt <-  function(GeneSymbol, Species){
      ## we will output the protein function and uniprotID
      res=data.frame(Function=NA, uniprotID=NA)
      
      message(paste0("Gene symbol:", GeneSymbol))
      # return NA for both if unknow protein
      if (!is.na(GeneSymbol)==T){ 
        tid = uniprotOrgs[uniprotOrgs$`Official (scientific) name` %in% 
                            Species,"Taxon Node"]
        ### Create a UniProt.ws object:
        up = UniProt.ws(taxId = tid)
        ### Select the first entry
        dfdesc = UniProt.ws::select(x = up,  keys = GeneSymbol, 
                                    keytype = "Gene_Name",
                                    columns = c("protein_name", "cc_function"))[1,]
        res$Function = dfdesc$Function..CC.
        res$uniprotID = dfdesc$Entry
      }
      return(res)
    }
    
    system.time(
      annotDF_complete <- annotDF %>%
        rowwise() %>%
        dplyr::mutate(Function = findFunctionProt(GeneSymbol, Species)$Function,
                      uniprotID = findFunctionProt(GeneSymbol, Species)$uniprotID) %>%
        data.frame()
    )
    
    ## Clean, manual curation for missing uniprotID
    annotDF_complete$GeneName = trimws(annotDF_complete$GeneName)
    annotDF_complete$GeneSymbol = trimws(annotDF_complete$GeneSymbol)
    
    annotDF_complete = annotDF_complete %>%
      mutate(uniprotID = case_when(
        GeneSymbol == "Galactose-specific lectin nattectin" & Species == "Thalassophryne nattereri" ~ "Q66S03",
        GeneSymbol == "Nuclear factor 7, ovary" & Species == "Xenopus laevis" ~ "Q91431",
        GeneSymbol == "Zinc-binding protein A33" & Species == "Pleurodeles waltl" ~ "Q02084",
        GeneSymbol == "Uncharacterized protein CXorf38 homolog" & Species == "Mus musculus" ~ "Q8C5K5",
        GeneSymbol == "Alpha-1-antitrypsin homolog" & Species == "Cyprinus carpio" ~ "P32759",
        TRUE ~ uniprotID  # Keep other values unchanged
      ))
    
    ## and save
    saveRDS(annotDF_complete, "../../dataOut/annotDF_complete.RDS")
  } else {
    annotDF_complete = readRDS("../../dataOut/annotDF_complete.RDS")
  }
  
  message(paste0("We have ",
                 sum(is.na(annotDF_complete$uniprotID[!annotDF_complete$Note %in% "Protein of unknown function"])),
                 " gene for which we miss a uniprotID"))
  
  # Order by nDMSperGenekb
  annotDF_complete = annotDF_complete[order(annotDF_complete$nDMSperGenekb, decreasing = T),]
  
  # set rownames NULL
  rownames(annotDF_complete) = NULL
  
  return(annotDF_complete)
}
