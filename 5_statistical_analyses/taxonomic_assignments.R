########################
# Ana Miller-ter Kuile
#Taxonomic Assignments####
# May 27, 2020
########################

# This code takes our matrix of samples and read abundances, assigns taxonomies to each of the 
#ASVs (splitting into pred-potential prey-other IDs) so that we can subset them at the end 
#to just do analyses on some sets of this data (namely, potential prey items, but perhaps predator
#as well)

##########################
# Required packages ####
library(here) #tidy data
library(tidyverse) #tidy data
###########################

###########################
# Import Community Data and taxonomies from BOLD and NCBI
###########################

comm <- read.csv(here("data", "denoised_data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = '\t')
#rename X to ASV in all these tables.
comm <- rename(comm, "ASV" = "X.OTU.ID")

#rename all the sample names across dataframes for consistency
colnames(comm) <- sapply(str_split(colnames(comm), "S"), function(x){return(x[[1]])})
comm <- rename(comm, "ASV" = "A")

#import both NCBI and BOLD taxonomies
ncbi <- read.csv(here("3_taxonomic_assignment", "MEGAN", "unoise_unclean",
                         "unoise_uc_ncbi_taxa.csv"))
bold <- read.csv(here("3_taxonomic_assignment", "BOLD", "usearch_uc", 
                         "unoise_uc_bold_taxa.csv"))
bold <- rename(bold, "ASV" = "Query.ID")

###########################
# Taxonomy categories
###########################

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it

#sort into categories for bold:
bold$taxonomy <- ifelse(
  bold$ID == "Heteropoda venatoria", "predator",
  ifelse(bold$Kingdom == "Fungi" | bold$Class == "Mammalia", "non-diet", 
         ifelse(bold$Phylum == "Arthropoda" | bold$Class == "Reptilia" & bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
ncbi$taxonomy <- ncbi$Category

#thus, the taxonomy variable has 4 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table)
#subset only the variables of interest from these two dataframes before merging
bold_id <- bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

ncbi_id <- ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
id_all <- bold_id %>%  
  full_join(ncbi_id, by = c("ASV", "taxonomy"))

id_all <- id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y)

###########################
# Bind community data and taxonomies
###########################

#this binds them all to the community matrix for analyses
id_tab <- comm %>%
  left_join(id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
id_tab$taxonomy <- replace_na(id_tab$taxonomy, "no hit")

###########################
# Create output for future analyses
###########################

#write taxonomies to a file for import into future analyses
taxonomies <- id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, taxonomy)

write.csv(taxonomies, here("data", "outputs", "taxonomic_assignments.csv"))
