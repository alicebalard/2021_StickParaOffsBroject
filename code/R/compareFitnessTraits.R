### Compare fitness traits between the different offsprings groups
### A. Balard nov. 2021
### Follow up of Sagonas 2020 & Ferre Ortega's master's dissertation

library(ggplot2)
library(tidyverse)
library(lme4) # linear mixed model

##############################
## Data import and cleaning ##
##############################
metadata <- readxl::read_xlsx("../../data/Kostas_G2_info.xlsx")
metadata_offsprings <- metadata[metadata$Generat  %in% "O",] # N=116
metadata_offsprings$Tr <- as.factor(metadata_offsprings$Tr)


table(metadata_offsprings$Tr) 
# CC CT TC TT 
# 30 28 29 29

## NB: Carles has a few fish less. Why?

## Kaufmann et al. 2014: Body condition of the G2fish, an estimate of fish health and a predictor
# of energy reserves and reproductive success, was calculated using there residuals from the 
# regression of body mass on body length (Chellappaet al.1995).

###########################
## Check BCI calculation ##
###########################

### to check with what's on Pangea for Kostas paper
metadata2 <- read.delim2("../../data/Exp-parasite-infections_fitness_nocom.tab")
# remove 2 fish with no weigth recorded
metadata2 <- metadata2[!metadata2$G..aculeatus.m..g. %in% "",]
metadata2$residuals <- residuals(lm(data = metadata2, formula = G..aculeatus.m..g. ~ G..aculeatus.TL..mm.))
## Compare previously calculated BCI and residuals
metadata2$BCI
round(metadata2$residuals, 5)
### Weird, not exactly the same. Must have been calculated on the bigger dataset

metadata_offsprings$newBCI <- residuals(lm(data = metadata_offsprings, formula = Wnettofin ~ Tlfin))

## Problem: inverse of Kostas's values
plot(metadata_offsprings$BCgen ~ metadata_offsprings$newBCI)

plot(metadata_offsprings$No.Worms...11, metadata_offsprings$No.Worms...23)
## Which one is correct????

#################################################################################
## Effect of parasite infection and parental background on fish body condition ##
#################################################################################
## ad a variable "father treatment"
metadata_offsprings$fatherTrt <- "Exposed"
metadata_offsprings$fatherTrt[metadata_offsprings$Tr %in% c("CC", "CT")] <- "Control"

# using Carles values (the ones present in the table)
## find the mean then plot
carlesF2DF <- metadata_offsprings %>% group_by(fatherTrt, No.Worms...11) %>% 
  summarise(BCgen = mean(BCgen)) %>% data.frame()
## Let's reproduce Carles figure 2:
ggplot(metadata_offsprings, aes(x=No.Worms...11, y = BCgen, group = fatherTrt, col = fatherTrt))+
  geom_point() +
  geom_line(data=carlesF2DF)+
  geom_point(data=carlesF2DF, aes(fill = fatherTrt), col = "black", size = 3, pch = 21)+
  scale_color_manual(values = c("gray", "red"))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_bw()

# using my values, it's reversed!
myBCdf <- metadata_offsprings %>% group_by(fatherTrt, No.Worms...11) %>% 
  summarise(newBCI = mean(newBCI)) %>% data.frame()
ggplot(metadata_offsprings, aes(x=No.Worms...11, y = newBCI, group = fatherTrt, col = fatherTrt))+
  geom_point() + geom_line(data=myBCdf)+
  geom_point(data=myBCdf, aes(fill = fatherTrt), col = "black", size = 3, pch = 21)+
  scale_color_manual(values = c("gray", "red"))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_bw()

## And separating all 4 treatments:
# using my values:
all4trtBCdf <- metadata_offsprings %>% group_by(Tr, No.Worms...11) %>% 
  summarise(newBCI = mean(newBCI)) %>% data.frame()
ggplot(metadata_offsprings, aes(x=No.Worms...11, y = newBCI, group = Tr, col = Tr))+
  geom_point() + geom_line(data=all4trtBCdf)+
  geom_point(data=all4trtBCdf, aes(fill = Tr, pch = Tr), col = "black", size = 3)+
  scale_color_manual(values = c("gray", "gray", "red","red"))+
  scale_fill_manual(values = c("gray", "gray", "red","red"))+
  scale_shape_manual(values = c(21,22,21,22))+
  theme_bw() 

## Statistical testing:
## Carles measures:
## all models
modCar <- lmer(BCgen ~ No.Worms...11 * fatherTrt + (1|Family), data=metadata_offsprings)
modCar_noparas <- lmer(BCgen ~ fatherTrt + (1|Family), data=metadata_offsprings)
modCar_nofather <- lmer(BCgen ~ No.Worms...11 + (1|Family), data=metadata_offsprings)
modCar_noint <- lmer(BCgen ~ No.Worms...11 + fatherTrt + (1|Family), data=metadata_offsprings)
## anovas
anova(modCar, modCar_noparas) # parasite load significant p=0.004
anova(modCar, modCar_nofather) # parental treatment significant p=<0.001
anova(modCar, modCar_noint) # interaction significant p =0.04

## New measures for BCI:
## all models
mod <- lmer(newBCI ~ No.Worms...11 * fatherTrt + (1|Family), data=metadata_offsprings)
mod_noparas <- lmer(newBCI ~ fatherTrt + (1|Family), data=metadata_offsprings)
mod_nofather <- lmer(newBCI ~ No.Worms...11 + (1|Family), data=metadata_offsprings)
mod_noint <- lmer(newBCI ~ No.Worms...11 + fatherTrt + (1|Family), data=metadata_offsprings)
## anovas
anova(mod, mod_noparas) # parasite load significant p=0.001
anova(mod, mod_nofather) # parental treatment significant p<0.001
anova(mod, mod_noint) # interaction significant p =0.03

## Redo the fig 3 of Kaufmann et al. 2014 with this subgroup
ggplot(metadata_offsprings, aes(x=Tr, y = newBCI))+
  geom_boxplot() +
  facet_grid(.~fatherTrt, scales = "free_x") + theme_bw()
## -> same results than Kaufmann et al.

mod2 <- lmer(newBCI ~ Tr + (1|Family), data = metadata_offsprings)
mod2_notrt <- lmer(newBCI ~ (1|Family), data = metadata_offsprings)
anova(mod2, mod2_notrt) # p < 0.001

summary(mod2)

library(emmeans)
emmeans(mod2, list(pairwise ~ Tr), adjust = "tukey")
# $`pairwise differences of Tr`
# 1       estimate   SE  df t.ratio p.value
# CC - CT    73.08 16.2 109   4.498  0.0001 ***
# CC - TC   -25.07 16.1 109  -1.557  0.4073
# CC - TT    -3.21 16.1 109  -0.199  0.9972
# CT - TC   -98.15 16.4 109  -5.987  <.0001 ***
# CT - TT   -76.29 16.4 109  -4.644  0.0001 ***
# TC - TT    21.86 16.3 109   1.345  0.5365

## Control father - treatment offsprng has a lower BC than eery other group: same as Kaufmann et al. 2014



