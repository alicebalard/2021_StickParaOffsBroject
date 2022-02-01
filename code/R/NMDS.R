library(vegan)
library(ggplot2)
library(RColorBrewer)
library(MASS)
library(plyr)
#load libraries - I haven't checked if all are actually required for what I am sending, I quickly removed a lot of other plots I was creating in the same R script

P1 <- read.csv("~/UCL/NHM_Project/Data/Pres_Phylum.csv", header=T, row.names=1)
metadata <- read.csv("~/UCL/NHM_Project/Data/MTMapping_16S2.csv", header=T)
#Read in data and metadata files

TP1 <- as.data.frame(t(P1))
#transposing data - you may not need this step

P1_Dist  = as.matrix((vegdist(TP1, "bray"))) 
#Creating a bray-curtis distance matrix - ?vegdist() for other distance calculation methods

NMDS <- metaMDS(P1_Dist)
#Create NMDS based on bray-curtis distances - metaMDS finds the most stable NMDS solution by randomly starting from different points in your data

stressplot(NMDS)
#check to see stress of NMDS

MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
#extract plotting coordinates

NMDS2 = data.frame( MDS1 = MDS1, MDS2 = MDS2, Name = metadata$ï..SampleID, Group2=as.factor(metadata$Group2), Group=as.factor(metadata$Group), Depth = as.factor(metadata$Depth), Mat = metadata$Mat_Type,Type = metadata$Type, Conductivity = metadata$Conductivity, NH4=metadata$NH4, NO3=metadata$NO3)
#create new dataframe with plotting coordinates and variables I was testing to see if I can group points by
#You only need one 



################### This step is optional, and only necessary if you need to create your own convex hulls  ############
find_hull <- function(my_data) my_data[chull(my_data[,1], my_data[,2]), ]
#function to create convex hulls, though ggplot has the stat_ellipse() function that can do this automatically. 
hulls <- ddply(NMDS2, "Group2", find_hull)
#genetating convex hulls by "Group2" in my metadata

ggplot(NMDS2, aes(x=MDS1, y=MDS2, size=Conductivity)) +
  geom_point(aes(col=Depth, shape=Type)) +
  #this plots the points, and also changes their size based on my Conductivity measurements - I left it in to provide information on how adaptable ggplot is
  stat_ellipse() +
  #this should automatically draw 95% confidence elipses around your groupings, but you can play around with the settings
  geom_polygon(data = hulls, alpha = 0.3, aes(fill=factor(Group2))) +  
  #use this to draw hulls if you are not using the above function
  scale_fill_manual(values = c("#66FBFF","#6666FF","#8D8D8D",col),label=c("Shallow","Middle","Deep"))+
  #manually adjusting the colours of the points
  geom_text(aes(label=Name, size=11, col=Depth),hjust=1.3, vjust=0, show.legend=F) +
  #adding names next to each point h- and vjust change the horizontal and verical placement of the names, and I turned the legend off for this
  scale_size_continuous(range=c(4,8), breaks=c(0.75,5,10,15,30,60)) +
  #manually changing the size of the points as I'm using size to indicate conductivity 
  scale_shape_manual(values=c(16,17) ,labels = c("Whole Mat","Mat Cross Section")) +
  #changing the shape of certain points to indicate what type of microbial mat they came from
  theme_bw() + 
  #use a more publication friendly theme, another good one is theme_classic()
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+
  #increase size of the text - this may have been for a poster, so 20 could be huge
  guides(shape=guide_legend(override.aes=list(size=5)),col=guide_legend(override.aes=list(size=5)))
#Changing the text size of the legend. 
