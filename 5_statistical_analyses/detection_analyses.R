########################
# Ana Miller-ter Kuile
#Potential prey detection analyses ####
# May 27, 2020
########################

#These code look at how surface sterilization treatment changes our ability to
#detect potential prey items for consumers we offered a potential prey item
#in mesocosms and for consumers which were free to roam in the field and feed
#on potential prey (and come into contact with potential environmental contaminants)
#in a natural way.

##########################
# Required packages ####
library(here) #tidy data
library(tidyverse) #tidy data
library(ggplot2) #visualize data
#library(emmeans) #marginal means of models
#library(rcompanion) #pairwise comparisions
#library(effects) #dotplots and stuff
#library(lattice) #again dotplots
library(DHARMa) #model diagnostics
library(MuMIn) #model diagnostics
library(glmmTMB) #mixed models
library(lme4) #mixed models
#library(ggeffects) #ggpredict function
library(performance) #for binomail model fit
library(see) #for binomial model fit
#library(coin) #wilcoxon test
#library(MASS) #glm.nb
#library(vegan) #rrarefy
library(cowplot) #plot grid at end
#library(picante) #phylogenetic analyses
#library(ape) #phylogenetic analyses
#library(phytools) #phylogenetic analyses
#library(phyloseq) #phylogenetic analyses
#library(hciR) #for as_matrix function
#library(esc) #effect sizes
#library(effsize) #alternative effect size package
###########################



#We will be using our taxonomic assignment dataframe created in the 
#taxonomic_assignment.R script as well as the rarefied community dataframes 
#created in the divide_rarefy.R script.
#we will also need our sample data so we know which samples were surface sterilized
#versus left unsterilized

lab_rare <- read.csv(here("data", "outputs", "lab_comm_rare.csv"))

field_rare <- read.csv(here("data", "outputs", "field_comm_rare.csv"))

#sample metadata
metadata <- read.csv(here("data", "Sample_Metadata.csv"))

#taxonomies
taxonomies <- read.csv(here("data", "outputs", "taxonomic_assignments.csv"))



u3_comm_rare$ASV <- rownames(u3_comm_rare)

u3_comm_rare %>%
  left_join(u3_id_all, by = "ASV") %>%
  filter(ID_ncbi == "Acrididae") %>%
  tally()

#this binds them all to the community matrix for analyses
u3_rare <- u3_comm_rare %>%
  left_join(u3_id_all, by = "ASV") %>% #combine with taxonomy DF
  filter(ID_ncbi == "Acrididae") %>% #select known diet only
  dplyr::select(-ID_bold, -ID_ncbi) %>% #remove taxonomy columns
  group_by(ASV) %>% 
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>% #make long while grouped by ASV
  group_by(sample) %>%
  summarise(reads = sum(reads)) #summarise total prey reads by sample 


all_rare <- u3_rare %>%
  left_join(metadata, by = "sample")