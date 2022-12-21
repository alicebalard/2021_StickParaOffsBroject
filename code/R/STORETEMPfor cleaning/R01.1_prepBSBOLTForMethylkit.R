## Data preparation for methylKit
## A. Balard
## 3rd of August 2021

temp = list.files(path="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted",
                  pattern = "*_trimmed_cutadapt.fastq.gz.bam.sorted.CGmap.gz.CG.map.gz",
                  full.names = T)

length(temp) # check: must be 144

head(temp)

## import files into a list (10 minutes)
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

### Files moved after to:
# /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted/formatCG4methylKit/
