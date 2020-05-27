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

##########################
# Import Community matrices, metadata, and taxonomic data ####
##########################

#the lab rarefied dataframe from the divide_rarefy.R script
lab_rare <- read.csv(here("data", "outputs", "lab_comm_rare.csv"))

#the field rarefied dataframe from the divide_rarefy.R script
field_rare <- read.csv(here("data", "outputs", "field_comm_rare.csv"))

#sample metadata
metadata <- read.csv(here("data", "Sample_Metadata.csv"))

#taxonomies dataframe from taxonomic_assignments.R
taxonomies <- read.csv(here("data", "outputs", "taxonomic_assignments.csv"))

##########################
# Mesocosms: Combine and tidy dataframes for analysis ####
##########################

#bind community matrix to taxonomies, subset prey item,
#gather dataframe so it's ready for a glm analysis, and 
#attach metadata on surface sterilization treatments
lab_rare <- lab_rare %>%
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #combine with taxonomy DF
  filter(ID_ncbi == "Acrididae") %>% #select known diet only
  dplyr::select(-ID_bold, -ID_ncbi) %>% #remove taxonomy columns
  group_by(ASV) %>% 
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>% #make long while grouped by ASV
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>% #summarise total prey reads by sample 
  left_join(metadata, by = "sample") #attach metadata on sterilization treatment
 
##########################
# Field: Combine and tidy dataframes for analysis ####
##########################

field_rare <- field_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "prey") #filter out just the prey taxonomies

#These prey ASVs contain duplicates of the same species, which we want to
#concatenate before analyses

#make the ID columns characters for the ifelse statements below
field_rare$ID_bold <- as.character(field_rare$ID_bold)
field_rare$ID_ncbi <- as.character(field_rare$ID_ncbi)

#ifelse statement that assigns each ASV to a specific taxonomy
#I want the BOLD ids to be the last option because these are to
#species in many cases where NCBI/Genbank blasted to order, family, or genus
#therefore, I'm going to say that if there is no BOLD id, give it the NCBI
#ID, otherwise, give it the BOLD ID.
unique_ID <- ifelse(is.na(field_rare$ID_bold), field_rare$ID_ncbi, 
                 field_rare$ID_bold)

#look at what these are:
unique_ID

#make this a dataframe so we can attach it back to the community
#matrix for analyses at the sample level and later Jaccard dissimiliarity
unique_ID <- as.data.frame(unique_ID)
unique_ID$ASV <- field_rare$ASV #give it an ASV column for attaching later

field_rare <- field_rare %>%
  left_join(unique_ID, by = "ASV")  %>%
  gather(sample, reads, HEV65:HEV100) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  left_join(metadata, by = "sample")
  
#Some of these have zero counts across all samples, so I'd like to delete
#them for later. 

fld_zeros <- field_rare %>%
  ungroup() %>%
  group_by(unique_ID) %>% #group by unique_ID
  summarise(reads = sum(reads)) %>% #sum across all samples
  filter(reads == 0) %>% #find the ones that are equal to zero
  dplyr::select(unique_ID) #create a dataframe of just their names

field_rare <- field_rare %>%
  anti_join(fld_zeros, by = "unique_ID") #anti join to remove the IDs that are zero across all smaples

##########################
#species richness of prey in each pipeline
all_rich <- test %>%
  group_by(sample, unique_ID, Sterilized) %>%
  summarize(reads = sum(reads)) %>% #this summarizes reads by each sp in each sample
  ungroup() %>%
  group_by(sample, Sterilized) %>% #then if we group by just sample, can count number of >0
  summarize(SR = sum(reads > 0, na.rm=TRUE)) %>% #count number of species >0
  mutate(pipeline = "u3") #add pipeline column to bind later

all_abund <- test %>%
  group_by(sample, Sterilized) %>% 
  summarize(abund = sum(reads)) %>% #this counts total abundance of prey reads in a sample
  mutate(pipeline = "u3", presence = ifelse(abund > 0, 1, 0)) #this says wheter prey detected
