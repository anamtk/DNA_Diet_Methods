########################
# Ana Miller-ter Kuile
#Rarefied taxonomically sorted ####
# May 27, 2020
########################

#These code create dataframes of the different "kinds" of DNA in each sample for
#each environment. This includes known prey (for mesocosms), all prey, and predator
#DNA (for both mesocosm and field predators)

##########################
# Required packages ####
library(here) #tidy data
library(tidyverse) #tidy data
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
lab_rare <- read.csv(here("data", "outputs", "divide_rarefy", "lab_comm_rare.csv"))

#the field rarefied dataframe from the divide_rarefy.R script
field_rare <- read.csv(here("data", "outputs", "divide_rarefy", "field_comm_rare.csv"))

#sample metadata
metadata <- read.csv(here("data", "Sample_Metadata.csv"))

#taxonomies dataframe from taxonomic_assignments.R
taxonomies <- read.csv(here("data", "outputs", "taxonomic_assignments", "taxonomic_assignments.csv"))

##########################
# Mesocosms: Known Prey ####
##########################

#bind community matrix to taxonomies, subset prey item,
#gather dataframe so it's ready for a glm analysis, and 
#attach metadata on surface sterilization treatments
lab_known_prey <- lab_rare %>%
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #combine with taxonomy DF
  filter(ID_ncbi == "Acrididae") %>% #select known diet only
  dplyr::select(-ID_bold, -ID_ncbi) %>% #remove taxonomy columns
  group_by(ASV) %>% 
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>% #make long while grouped by ASV
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>% #summarise total prey reads by sample 
  left_join(metadata, by = "sample") #attach metadata on sterilization treatment
 
#write output to a file for later analysis
write.csv(lab_known_prey, here("data", "outputs", "rarefied_taxonomic_sort", "lab_known_prey_rare.csv"))
##########################
# Mesocosms: All Prey ####
##########################
lab_prey <- lab_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "prey") #filter out just the prey taxonomies

#make the ID columns characters for the ifelse statements below
lab_prey$ID_bold <- as.character(lab_prey$ID_bold)
lab_prey$ID_ncbi <- as.character(lab_prey$ID_ncbi)

#ifelse statement that assigns each ASV to a specific taxonomy
#I want the BOLD ids to be the last option because these are to
#species in many cases where NCBI/Genbank blasted to order, family, or genus
#therefore, I'm going to say that if there is no BOLD id, give it the NCBI
#ID, otherwise, give it the BOLD ID.
lab_unique_ID <- ifelse(is.na(lab_prey$ID_bold), lab_prey$ID_ncbi, 
                    lab_prey$ID_bold)

#look at what these are:
lab_unique_ID

#make this a dataframe so we can attach it back to the community
#matrix for analyses at the sample level and later Jaccard dissimiliarity
lab_unique_ID <- as.data.frame(lab_unique_ID)
lab_unique_ID$ASV <- lab_prey$ASV #give it an ASV column for attaching later

lab_unique_ID <- lab_unique_ID %>%
  rename("unique_ID" = "lab_unique_ID")

lab_prey <- lab_prey %>%
  left_join(lab_unique_ID, by = "ASV")  %>%
  gather(sample, reads, HEV07:HEV29) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  left_join(metadata, by = "sample")

#Some of these have zero counts across all samples, so I'd like to delete
#them for later. 

lab_zeros <- lab_prey %>%
  ungroup() %>%
  group_by(unique_ID) %>% #group by unique_ID
  summarise(reads = sum(reads)) %>% #sum across all samples
  filter(reads == 0) %>% #find the ones that are equal to zero
  dplyr::select(unique_ID) #create a dataframe of just their names

lab_prey <- lab_prey %>%
  anti_join(lab_zeros, by = "unique_ID") #anti join to remove the IDs that are zero across all smaples

write.csv(lab_prey, here("data", "outputs", "rarefied_taxonomic_sort", "lab_all_prey_rare.csv"))

##########################
# Mesocosms: Predator ####
##########################

lab_predator <- lab_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "predator") %>% #filter out just the prey taxonomies
  gather(sample, reads, HEV07:HEV29, factor_key =TRUE) %>%
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>%
  left_join(metadata, by = "sample")

write.csv(lab_predator, here("data", "outputs", "rarefied_taxonomic_sort", "lab_predator_rare.csv"))
##########################
# Field: Prey ####
##########################

field_prey <- field_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "prey") #filter out just the prey taxonomies

#These prey ASVs contain duplicates of the same species, which we want to
#concatenate before analyses

#make the ID columns characters for the ifelse statements below
field_prey$ID_bold <- as.character(field_prey$ID_bold)
field_prey$ID_ncbi <- as.character(field_prey$ID_ncbi)

#ifelse statement that assigns each ASV to a specific taxonomy
#I want the BOLD ids to be the last option because these are to
#species in many cases where NCBI/Genbank blasted to order, family, or genus
#therefore, I'm going to say that if there is no BOLD id, give it the NCBI
#ID, otherwise, give it the BOLD ID.
unique_ID <- ifelse(is.na(field_prey$ID_bold), field_prey$ID_ncbi, 
                 field_prey$ID_bold)

#look at what these are:
unique_ID

#make this a dataframe so we can attach it back to the community
#matrix for analyses at the sample level and later Jaccard dissimiliarity
unique_ID <- as.data.frame(unique_ID)
unique_ID$ASV <- field_prey$ASV #give it an ASV column for attaching later

field_prey <- field_prey %>%
  left_join(unique_ID, by = "ASV")  %>%
  gather(sample, reads, HEV65:HEV100) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  left_join(metadata, by = "sample")
  
#Some of these have zero counts across all samples, so I'd like to delete
#them for later. 

fld_zeros <- field_prey %>%
  ungroup() %>%
  group_by(unique_ID) %>% #group by unique_ID
  summarise(reads = sum(reads)) %>% #sum across all samples
  filter(reads == 0) %>% #find the ones that are equal to zero
  dplyr::select(unique_ID) #create a dataframe of just their names

field_prey <- field_prey %>%
  anti_join(fld_zeros, by = "unique_ID") #anti join to remove the IDs that are zero across all smaples

write.csv(field_prey, here("data", "outputs", "rarefied_taxonomic_sort", "field_prey_rare.csv"))

##########################
# Field: Predator ####
##########################

field_predator <- field_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "predator") %>% #filter out just the prey taxonomies
  gather(sample, reads, HEV65:HEV100, factor_key =TRUE) %>%
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>%
  left_join(metadata, by = "sample")

write.csv(field_predator, here("data", "outputs", "rarefied_taxonomic_sort", "field_predator_rare.csv"))

