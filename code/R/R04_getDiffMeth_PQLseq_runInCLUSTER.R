## Data preparation for PQLseq
## A. Balard
## 2nd of August 2024

## load previous scripts
source("R03_prepObjectMethylkit_runInCLUSTER.R")

## set run = TRUE at the last block of code to re-run the comparisons

#############################
## PQLSeq node alarm issue ##
## library(RhpcBLASctl) # To deal with alarming nodes issue (in librarie script)

## Set multi-threading parameters to avoid alarming nodes
Sys.setenv(MKL_NUM_THREADS = 1)
blas_set_num_threads(1)
omp_set_num_threads(1)

## New function made by Charley
source("/data/SBCS-EizaguirreLab/Turtle_WGBS/00_Scripts/Functions/Custom_PQLseq_Function.R")

# Set up environment for custom function: 
# https://stackoverflow.com/questions/24331690/modify-package-function
## assures that the function will be able to call other hidden functions from the package:
environment(custom_pqlseq) <- asNamespace('PQLseq')

## assures that other functions from the package will call your updated version of the function:
assignInNamespace("pqlseq", custom_pqlseq, ns = "PQLseq")

######################################
## Parents file filtered and normalised with methylkit
## Sanity check
uniteCovHALF_G1_woSexAndUnknowChrOVERLAP %>% nrow == 1002565

## Offspring file filtered and normalised with methylkit
## Sanity check
uniteCovHALF_G2_woSexAndUnknowChrOVERLAP %>% nrow == 1002565

##################################
## PQLseq differential analysis ##
##################################

## PQLseq requires four input files containing 
# (1) methylated read counts, (2) total read counts, 
# (3) relatedness matrix, and (4) predictor variable of interest
## Optional: covariate files (we use sex as covariate)

####################
## Prepare relatedness matrix: 1 against oneself, 0.5 between sibling, 0.25 uncle/nephew-niece, 0.125 cousins
makefile.3 <- function(met = fullMetadata){
  relatMat <- matrix(nrow=nrow(met), ncol=nrow(met), data = 0)
  rownames(relatMat)<-met$SampleID
  colnames(relatMat)<-met$SampleID
  ## coef of relatedness with oneself = 1
  diag(relatMat) <- 1
  ## Same clutch: coef of relatedness with your sibling/father/offspring = 0.5 
  for (i in 1:nrow(met)) {
    for (j in 1:nrow(met)) {
      if (relatMat[i, j]==0 && met$clutch.ID[i] == met$clutch.ID[j]) {
        relatMat[i, j] = 0.5
      }}
  }
  ## Same BP, diff clutch (populated step before), gener P: coef of relatedness with your brother in G1 = 0.5
  for (i in 1:nrow(met)) {
    for (j in 1:nrow(met)) {
      if (relatMat[i, j]==0 &&  met$brotherPairID[i] == met$brotherPairID[j] &&
          met$Generat[i] == met$Generat[j] && met$Generat[i] == "P") {
        relatMat[i, j] = 0.5
      }}
  }
  ## Same BP, diff clutch (populated 2 steps before), gener O: coef of relatedness with your cousins in G2 = 0.125
  for (i in 1:nrow(met)) {
    for (j in 1:nrow(met)) {
      if (relatMat[i, j] == 0 && met$brotherPairID[i] == met$brotherPairID[j] &&
          met$Generat[i] == met$Generat[j] && met$Generat[i] == "O"){
        relatMat[i, j] = 0.125
      }}
  }
  ## Same BP, diff clutch, diff gener: coef of relatedness with your uncle/niece/nephew = 0.25
  for (i in 1:nrow(met)) {
    for (j in 1:nrow(met)) {
      if (relatMat[i, j] == 0 && met$brotherPairID[i] == met$brotherPairID[j]&&
          met$Generat[i] != met$Generat[j]) {
        relatMat[i, j] = 0.25
      }}
  }
  return(relatMat)
}

relatMat <- makefile.3()

## Sanity check done in IDE
## load ggcorr function instead of full GGally
## source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R")
## ggcorr(data = NULL, cor_matrix = relatMat, size = 2, color = "grey50")
## Ok, all the relationships within the 8 families are well represented

#############################################
### 5 comparisons to do:
## 0. PARENTS trt-ctrl
## 1. CC-TC = CONTROL fish (parent CvsT)
## 2. CT-TT = TREATMENT fish (parent CvsT)
## 3. CC-CT = fish from CONTROL parents (G2 CvsT)
## 4. TC-TT = fish from TREATMENT parents (G2 CvsT)

## 1&2: paternal effect 
## 3&4: G2 direct treatment effect

################################################################
## Correspondance trtG1G2 and numerical values used by MethylKit
## table(fullMetadata$trtG1G2_NUM, fullMetadata$trtG1G2)
## Control C 1
## Exposed T 4
## NE_control CC 5
## NE_exposed CT 6
## E_control TC 2
## E_exposed TT 3

########################################################################
####################### launch PQLseq analysis #########################
## Calculate DMS/DMR accounting for covariates: brotherPairID and sex ##
########################################################################
## NB: make sure the order of individuals in the kinship.txt, the order of individuals
## in the read counts and total counts, matches the order of individuals in the phenotype file
## pheno.txt.
getDiffMeth_PQLseq <- function(uniteCov, groupa, groupb, subset = FALSE){
  
  sub_uniteCov = methylKit::reorganize(
    methylObj = uniteCov,
    sample.ids = uniteCov@sample.ids[uniteCov@treatment %in% c(groupa, groupb)], 
    treatment = uniteCov@treatment[uniteCov@treatment %in% c(groupa, groupb)])
  
  if (subset == TRUE){ # subset the first 100 sites for debugging
    sub_uniteCov = sub_uniteCov[1:100,]
  }
  
  ## (1) methylated read counts
  ## C represents a methylated cytosine that remained unchanged during bisulfite treatment
  subcounts=methylKit::getData(sub_uniteCov)[ ,grep("numCs", names(sub_uniteCov))]
  colnames(subcounts)=sub_uniteCov@sample.ids
  row.names(subcounts)=paste(sub_uniteCov$chr, sub_uniteCov$start, sep ="_")
  
  ## (2) total read counts
  subtotalcounts=methylKit::getData(sub_uniteCov)[ ,grep("coverage", names(sub_uniteCov))]
  colnames(subtotalcounts)=sub_uniteCov@sample.ids
  row.names(subtotalcounts)=paste(sub_uniteCov$chr, sub_uniteCov$start, sep ="_")
  
  ## (3) relatedness matrix
  subrelatMat = relatMat[rownames(relatMat) %in% sub_uniteCov@sample.ids, 
                         colnames(relatMat) %in% sub_uniteCov@sample.ids]
  
  ## (4) phenotype file
  subpheno = sub_uniteCov@treatment
  ### make binary:
  subpheno[subpheno %in% groupa] = 0
  subpheno[subpheno %in% groupb] = 1
  
  ## (5) covariate file (sex)
  subcovariate = as.matrix(data.frame(1, sex=ifelse(fullMetadata$Sex[
    fullMetadata$SampleID %in% sub_uniteCov@sample.ids] =="M",1,2)))
  
  subcovariate = as.matrix(ifelse(
    fullMetadata$Sex[fullMetadata$SampleID %in% sub_uniteCov@sample.ids] =="M",1,2))
  
  ## Calculate Wald test diff meth 
  message(paste("Comparison between group ", groupa, " and ", groupb))

    ## here the custom pqlseq should be loaded, see start of the script
  fit = pqlseq(RawCountDataSet=subcounts, Phenotypes=subpheno, RelatednessMatrix=subrelatMat,
               LibSize=subtotalcounts, fit.model="BMM", Covariates = subcovariate)
  
  message(paste0(sum(fit$converged), " CpG sites converged out of ", length(fit$converged), " sites"))
  
  ## Calculate q-value from FDR B-H test
  fit$qvalue = qvalue::qvalue(fit$pvalue)$qvalue
  
  ## Add mean differential methylation
  fit$aveDiffMeth_ab = rowSums(subcounts[subpheno %in% 0]/subtotalcounts[subpheno %in% 0], na.rm = T) - 
    rowSums(subcounts[subpheno %in% 1]/subtotalcounts[subpheno %in% 1], na.rm = T)
  
  return(fit=fit)
}

## NE_control CC 5
## NE_exposed CT 6
## E_control TC 2
## E_exposed TT 3

run = F
if (run == T){
    ## 1. CC-TC = CONTROL fish (parent CvsT)
    fit_G2_CC.TC = getDiffMeth_PQLseq(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 5, 2)

    ## 2. CT-TT = TREATMENT fish (parent CvsT)
    fit_G2_CT.TT = getDiffMeth_PQLseq(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 6, 3)

    ## 3. CC-CT = fish from CONTROL parents (G2 CvsT)
    fit_G2_CC.CT = getDiffMeth_PQLseq(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 5, 6)

    ## 4. TC-TT = fish from TREATMENT parents (G2 CvsT)
    fit_G2_TC.TT = getDiffMeth_PQLseq(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 2, 3)

    ## Add the chromosomal position as data columns
    addChrPos <- function(df){
      df$chrom=paste(sapply(strsplit(row.names(df), "_"), `[`, 1),
                   sapply(strsplit(row.names(df), "_"), `[`, 2), sep = "_")
      df$start=as.numeric(sapply(strsplit(row.names(df), "_"), `[`, 3))
      df$end=as.numeric(sapply(strsplit(row.names(df), "_"), `[`, 3))
      df$pos=paste(df$chrom, df$start, sep = "_")
      return(df)
    }
    
    fit_G2_CC.CT = addChrPos(fit_G2_CC.CT)
    fit_G2_CC.TC = addChrPos(fit_G2_CC.TC)
    fit_G2_CT.TT = addChrPos(fit_G2_CT.TT)
    fit_G2_TC.TT = addChrPos(fit_G2_TC.TT)
    
    save(fit_G2_CC.CT = fit_G2_CC.CT, file = "../../dataOut/fitPQLseqG2_fit_G2_CC.CT.RData")
    save(fit_G2_CC.TC = fit_G2_CC.TC, file = "../../dataOut/fitPQLseqG2_fit_G2_CC.TC.RData")
    save(fit_G2_CT.TT = fit_G2_CT.TT, file = "../../dataOut/fitPQLseqG2_fit_G2_CT.TT.RData")
    save(fit_G2_TC.TT = fit_G2_TC.TT, file = "../../dataOut/fitPQLseqG2_fit_G2_TC.TT.RData")
}

