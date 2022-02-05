Explanation of the different files:

customRfunctions.R for dendogram plots 
R00_calculateStats.R general stats on previously done alignment, methylation call etc.
R01.1_prepBSBOLTForMethylkit.R format BSBOLT output to feed to MethylKit
R01.2_prepObjectMethylkit.R output MethylKit objects with CpG shared by all, 2 or 6 fish per group
!! R01.3_prepMetadata.R TO BE SOURCED IN FURTHER SCRIPTS format the different metadata
!! R01.4_prepMethyldata.R TO BE SOURCED IN FURTHER SCRIPTS format the methylkit objects (all, offsprings, parents)
R02_methylationAndFitness.R link methylation and BCI
R03.1_globalMethylation.R cluster based on methylation profiles, NMDS
R03.2_differentialMethylation.R calculate DMS between different groups

