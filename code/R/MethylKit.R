## Data preparation for methylKit
## A. Balard
## 3rd of August 2021

library(methylKit)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html

################ PART I 

##### Load and prepare dataset #####
temp = list.files(path="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted",
                  pattern = "*_trimmed_cutadapt.fastq.gz.bam.sorted.CGmap.gz.CG.map.gz",
                  full.names = T)

length(temp) # check: must be 144

head(temp)

## import files into a list (long)
myfiles = lapply(temp, read.csv, sep="\t", header=F)

## name the list myfiles
names(myfiles) = lapply(temp, function(x) gsub(pattern = "(.*BSBoltAlignments_)(.*)(_trimmed.*)", replacement = "\\2", x))

names(myfiles)

class(myfiles)

## add column names
new.names=c("chrom", "nucleotide", "position", "context", "sub-context", "methylation_value", "methylated_bases", "all_bases")
myfiles=lapply(myfiles, setNames, new.names)

## Transform BSBolt output format into MethylKit input format:
myrenameFUN <- function(BSBDF){
    MKDF=data.frame(chrBase=paste(BSBDF$chrom,BSBDF$position, sep = "."),
                chr=BSBDF$chrom,
                base=BSBDF$position,
                strand=ifelse(BSBDF$nucleotide=="C", yes = "F", no = "R"),
                coverage=BSBDF$all_bases,
                freqC=round(BSBDF$methylation_value*100, 2),
                freqT=round((1-BSBDF$methylation_value)*100,2))
    return(MKDF)
}

## LONG: make compatible DF
myfilesMK=lapply(myfiles, myrenameFUN)       

## Output the transformed files
length(myfilesMK)

for(i in 1:length(myfilesMK)){
    write.table(myfilesMK[[i]],
                file=paste0("/data/scratch/btx915/BSBolt/MethylationCallingClean/formatCG4methylKit/", names(myfilesMK)[[i]],".CG4methylkit.txt"),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep= "\t")
}

################# PART 2 

# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#load-datasets

# /data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject_archive/08CompGenomes_mCextr/04RAnalyses_Methylome/methylbam_peichel...

## Reading the methylation call files and store them as flat file database
## read the files to a methylRawListDB object: myobjDB and save in databases in folder methylDB
temp2 = list.files(path="/data/scratch/btx915/BSBolt/MethylationCallingClean/formatCG4methylKit",
                   pattern = ".CG4methylkit.txt",
                   full.names = T)

temp2

## Metadata (manually added)
treatment=c("O_Ec","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","P_NE","P_E","O_Ec","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","O_NEe","P_NE","P_E","O_Ec","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","P_E","P_NE","O_Ec","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","O_NEe","P_NE","P_E","P_E","P_NE","P_NE","P_E","P_E","P_NE","P_NE","P_E","P_NE","P_E","O_Ec","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","P_NE","P_E","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","O_NEe","P_NE","P_E","O_Ec","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","O_NEe","P_E","P_NE","P_NE","P_E","O_Ec","O_Ec","O_Ec","O_Ee","O_Ee","O_Ee","O_Ee","O_NEc","O_NEc","O_NEc","O_NEc","O_NEe","O_NEe","O_NEe","O_NEe","P_NE","P_E")

### MakeDB
myobjDB=methRead(as.list(temp2),
                 mincov=10,
                 sample.id=as.list(names(myfilesMK)),
                 assembly="Gynogen_pchrom_assembly_all",
                 treatment=treatment,
                 context="CpG",
                 dbtype = "tabix",
                 dbdir = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylDB")

# to use -> myobjDB

print(myobjDB[[1]]@dbpath)

## NB: Most if not all functions in this package will work with methylDB objects the same way as it does with normal methylKit objects. Functions that return methylKit objects, will return a methylDB object if provided, but there are a few exceptions such as the select, the [ and the selectByOverlap functions.

############### End data preparation #########

dbdir = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylDB/G05016-L1_S144_L008_R1_001.bgz.tbi"

myobjDB = readMethylDB(dbdir)




getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

baseDB.obj <- makeMethylDB(methylBase.obj,"my/path")
mydbpath <- getDBPath(baseDB.obj)
rm(baseDB.obj)
readMethylDB(mydbpath)

dbdir = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylDB"

myobjDB = readMethylDB(dbdir)


You can also convert methylDB objects to their in-memory equivalents. Since that requires an additional parameter (the directory where the files will be located), we have a different function, named makeMethylDB to achieve this goal. Below, we convert a methylBase object to methylBaseDB and saving it at “exMethylDB” directory.

data(methylKit)
 
objDB=makeMethylDB(methylBase.obj,"exMethylDB")
