# Adapted from Melanie Heckwolf
# create Manhattan plots over the genome from a DMS file

#load packages:
library(dplyr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

#the stickleback genome is annotated as chrI, chrII, chrIII ... to chrUn (unknown, includes all scaffolds)
#this function takes the roman rumbers and turns them into actual numbers:
deroman <- function(x){ x %>% str_remove(.,"Gy_chr") %>% 
    ifelse(. %in% c("Un","M"),., as.roman(.) %>% 
             as.integer() %>% as.character() %>% str_pad(.,2,"left",0))
}

get_pos <- function(CHROM,genome){
  tibble(CHROM = CHROM, 
         POS = sample(x=1:genome$length[genome$chrom == CHROM],size=1))
}

## Load file containing length of each gynogen chromosomes 
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff <- read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) <- c("chrom","length")

#GA_genome.fa.sizes.txt is a file with chromosome sizes and names
genome <- GYgynogff %>%
  mutate(chrom_nr=chrom %>% deroman(),
         chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
  arrange(chrom_order) %>%
  mutate(gstart=lag(length,default=0) %>% cumsum(),
         gend=gstart+length, 
         type=LETTERS[2-(chrom_order%%2)],
         gmid=(gstart+gend)/2)

#genome without M re-type:
genome2=genome[genome$chrom_nr!="M",] %>%
  mutate(type=rep(c("A","B"),length(length)/2))

# load file with your DMS
DMS_PAR <- readRDS("../../data/myDiff1_15p_parentalDiffMeth.RDS")

dataPAR = tibble(chrom=DMS_PAR$chr,
                 pos=DMS_PAR$start,
                 meth.diff=DMS_PAR$meth.diff,
                 qval=DMS_PAR$qvalue)

table(DMS_PAR$chr)## ERROR chrUN should be removed!!

# join DMS and genomic position
data = left_join(dataPAR, genome2) %>% 
  mutate(gpos=pos+gstart,significance= ifelse(abs(qval>0.0125) | abs(meth.diff)<15,"not significant","significant"))

table(data$significance) # all signif

anno_range <- c(max(genome2$gend),max(genome2$gend)*1.2)

#plot

#only significant DMS:
ggplot()+
  geom_rect(data=genome2,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type))+
  geom_point(data=data[abs(data$meth.diff)>15 & data$significance=="significant",],aes(x=gpos,y=meth.diff))+
  scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none")+
  # scale_color_manual(values = clr,guide="none")+
  scale_x_continuous(breaks=genome2$gmid,
                     labels=genome2$chrom %>% str_remove(.,"chr"),
                     position = "top",
                     expand = c(0,0),
                     limits = c(0,anno_range[2]))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line=element_line(),
        axis.title = element_blank(),
        strip.placement = "outside")

#pdf saved as 8/14




# #only significant DMS, colored in genomic regions (intergenic, intron, exon, promoter)
# ggplot()+
#   geom_rect(data=genome2,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type))+
#   #annotation_custom(stickl_B6,xmin=anno_range[1],xmax=anno_range[2])+
#   geom_hypo_grob(data = anno_tab2,aes(x = .92, y = .6, grob = grob,angle=angle),height=.3)+
#   geom_point(data=data[abs(data$meth.diff)>15 & data$significance=="significant",],aes(x=gpos,y=meth.diff,col=region,shape=region),fill="white")+
#   facet_grid(pop~.)+
#   scale_fill_manual(values=c(A=bgcol,B=NA),guide=F)+
#   scale_color_manual(values = c("grey50","grey50","darkred","darkred"))+
#   scale_shape_manual(values=c(21,19,21,19),guide=F)+
#   labs(color="Genomic region:")+
#   scale_x_continuous(breaks=genome2$gmid,
#                      labels=genome2$chrom %>% str_remove(.,"chr"),
#                      position = "top",
#                      expand = c(0,0),
#                      limits = c(0,anno_range[2]))+
#   guides(color=guide_legend(override.aes = list(shape=c(21,19,21,19))))+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         axis.line=element_line(),
#         axis.title = element_blank(),
#         strip.placement = "outside",
#         legend.position = c(0.93, 0.55))




