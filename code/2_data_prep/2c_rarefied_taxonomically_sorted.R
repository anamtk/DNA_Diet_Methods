########################
# Ana Miller-ter Kuile
#Rarefied taxonomically sorted ####
# May 27, 2020
########################

#These code create dataframes of the different "kinds" of DNA in each sample for
#each environment. This includes known prey (for mesocosms), all prey, and predator
#DNA (for both mesocosm and field predators) and is an assessment of the rarefied data

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
lab_rare <- read.csv(here("data", 
                          "outputs", 
                          "2b_divide_rarefy", 
                          "lab_comm_rare.csv"))

#the field rarefied dataframe from the divide_rarefy.R script
field_rare <- read.csv(here("data", 
                            "outputs", 
                            "2b_divide_rarefy", 
                            "field_comm_rare.csv"))

#sample metadata
metadata <- read.csv(here("data", 
                          "data_raw",
                          "2_sample_data",
                          "Sample_Metadata.csv"))

#taxonomies dataframe from taxonomic_assignments.R
taxonomies <- read.csv(here("data", 
                            "outputs", 
                            "2a_taxonomic_assignments", 
                            "taxonomic_assignments.csv"))

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
write.csv(lab_known_prey, here("data", 
                               "outputs", 
                               "2c_rarefied_taxonomic_sort", 
                               "lab_known_prey_rare.csv"))
##########################
# Mesocosms: All Prey ####
##########################
lab_prey <- lab_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "prey") %>%#filter out just the prey taxonomies
  gather(sample, reads, HEV07:HEV29) %>% #make long
  group_by(sample, Family_ncbi, Order_ncbi) %>% #group by sample and family
  summarise(reads = sum(reads)) %>% #summarise at family level
  left_join(metadata, by = "sample") #join to consumer metadata
  
#Some of these have zero counts across all samples, so I'd like to delete
#them for later. 

lab_zeros <- lab_prey %>%
  ungroup() %>%
  group_by(Family_ncbi) %>% #group by family
  summarise(reads = sum(reads)) %>% #sum across all samples
  filter(reads == 0) %>% #find the ones that are equal to zero
  dplyr::select(Family_ncbi) #create a dataframe of just their names

lab_prey <- lab_prey %>%
  anti_join(lab_zeros, by = "Family_ncbi") #anti join to remove the IDs that are zero across all smaples

#export for analysis later
write.csv(lab_prey, here("data", 
                         "outputs", 
                         "2c_rarefied_taxonomic_sort", 
                         "lab_all_prey_rare.csv"))

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

write.csv(lab_predator, here("data", 
                             "outputs", 
                             "2c_rarefied_taxonomic_sort", 
                             "lab_predator_rare.csv"))
##########################
# Field: Prey ####
##########################

field_prey <- field_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "prey") %>%#filter out just the prey taxonomies
  gather(sample, reads, HEV65:HEV100) %>% #make long
  group_by(sample, Family_ncbi, Order_ncbi) %>% #group by sample and family
  summarise(reads = sum(reads)) %>% #summarise at family level
  left_join(metadata, by = "sample") #join to consumer metadata

field_prey_ASVs <- field_rare %>% 
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #join by taxonomy
  filter(taxonomy == "prey") %>%#filter out just the prey taxonomies
  gather(sample, reads, HEV65:HEV100) %>% #make long
  group_by(sample, ASV) %>% #group by sample and family
  summarise(reads = sum(reads)) %>% #summarise at family level
  left_join(metadata, by = "sample") #join to consumer metadata

#Some of these have zero counts across all samples, so I'd like to delete
#them for later. 

field_zeros <- field_prey %>%
  ungroup() %>%
  group_by(Family_ncbi) %>% #group by family
  summarise(reads = sum(reads)) %>% #sum across all samples
  filter(reads == 0) %>% #find the ones that are equal to zero
  dplyr::select(Family_ncbi) #create a dataframe of just their names

field_zeros_ASVs <- field_prey_ASVs %>%
  ungroup() %>%
  group_by(ASV) %>% #group by family
  summarise(reads = sum(reads)) %>% #sum across all samples
  filter(reads == 0) %>% #find the ones that are equal to zero
  dplyr::select(ASV) #create a dataframe of just their names

field_prey <- field_prey %>%
  anti_join(field_zeros, by = "Family_ncbi") #anti join to remove the IDs that are zero across all smaples

field_prey_ASVs <- field_prey_ASVs %>%
  anti_join(field_zeros_ASVs, by = "ASV")
#export for analysis later
write.csv(field_prey, here("data", 
                           "outputs", 
                           "2c_rarefied_taxonomic_sort", 
                           "field_prey_rare.csv"))

write.csv(field_prey_ASVs, 
          here("data", 
               "outputs", 
               "2c_rarefied_taxonomic_sort", 
               "field_prey_ASVs_rare.csv"))

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

write.csv(field_predator, here("data",
                               "outputs", 
                               "2c_rarefied_taxonomic_sort",
                               "field_predator_rare.csv"))

##########################
# Lab: Non-diet ####
##########################
lab_nondiet <- lab_rare %>%
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #combine with taxonomy DF
  filter(taxonomy == "non-diet") %>% #select non-diet only
  dplyr::select(-ID_bold, -ID_ncbi) %>% #remove taxonomy columns
  group_by(ASV) %>% 
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>% #make long while grouped by ASV
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>% #summarise total prey reads by sample 
  left_join(metadata, by = "sample") #attach metadata on sterilization treatment

write.csv(lab_nondiet, here("data", 
                            "outputs", 
                            "2c_rarefied_taxonomic_sort",
                            "lab_nondiet_rare.csv"))

##########################
# Field: Non-diet ####
##########################
fld_nondiet <- field_rare %>%
  rename("ASV" = "X") %>%
  left_join(taxonomies, by = "ASV") %>% #combine with taxonomy DF
  filter(taxonomy == "non-diet") %>% #select non-diet only
  dplyr::select(-ID_bold, -ID_ncbi) %>% #remove taxonomy columns
  group_by(ASV) %>% 
  gather(sample, reads, HEV65:HEV100, factor_key = TRUE) %>% #make long while grouped by ASV
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>% #summarise total prey reads by sample 
  left_join(metadata, by = "sample") #attach metadata on sterilization treatment

write.csv(fld_nondiet, here("data", 
                            "outputs", 
                            "2c_rarefied_taxonomic_sort", 
                            "field_nondiet_rare.csv"))


