########################
# Ana Miller-ter Kuile
#Rarefaction and subset mesocosm from field ####
# May 27, 2020
########################

# This code will provide subset field and mesocosm read abundances as well as rarefy both of these 
# dataframes for future analyses

##########################
# Required packages ####
library(here) #tidy data
library(tidyverse) #tidy data
library(vegan) #rrarefy function
###########################

###########################
# Import Community Data and metadata for subsetting
###########################

comm <- read.csv(here("data", "denoised_data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = '\t')
#rename X to ASV in all these tables.
comm <- rename(meso, "ASV" = "X.OTU.ID")

#rename all the sample names across dataframes for consistency
colnames(comm) <- sapply(str_split(colnames(comm), "S"), function(x){return(x[[1]])})
comm <- rename(comm, "ASV" = "A")

#sample metadata
metadata <- read.csv(here("data", "Sample_Metadata.csv"))

###########################
# Subset mesocosm from field spiders and create outputs
###########################

#vector of the sample names for lab samples
lab <- as.character(metadata$sample[which(metadata$Source == "LAB")]) 

field <- as.character(metadata$sample[which(metadata$Source == "FIELD")]) 

comm_lab <- comm
rownames(comm_lab) <- comm_lab$ASV
comm_lab <- comm_lab %>%
  dplyr::select(all_of(lab))

write.csv(comm_lab, here("data", "outputs", "lab_comm_raw.csv"))

comm_field <- comm
rownames(comm_field) <- comm_field$ASV
comm_field <- comm_field %>%
  dplyr::select(all_of(field))

write.csv(comm_field, here("data", "outputs", "field_comm_raw.csv"))

###########################
# Rarefy field and lab separately
###########################

#Assessing the variation in sequencing depth
hist(colSums(comm_lab))
hist(colSums(comm_field))

colSums(comm_lab)
colSums(comm_field)

#comparing max and min, huge diff
#max read abundance 
max(colSums(comm_lab))
max(colSums(comm_field))
#min read abundance
min(colSums(comm_lab))
min(colSums(comm_field))

#we will rarefy based on lowest sample read abundance for both lab and field 
lab_rare <- min(colSums(comm_lab))
lab_rare #55205

field_rare <- min(colSums(comm_field))
field_rare #16004

#rarefaction is a random process, so set seed for consistent results
set.seed(1)

lab_comm_rare <- as.data.frame(t(rrarefy(t(comm_lab), sample = lab_rare)))

field_comm_rare <- as.data.frame(t(rrarefy(t(comm_field), sample = field_rare)))

#now we see that the column sums are the same across, so we can be confident in comparing 
#these samples to each other
colSums(comm_lab)
colSums(comm_field)

###########################
# Create output for future analyses
###########################

write.csv(lab_comm_rare, here("data", "outputs", "lab_comm_rare.csv"))
write.csv(field_comm_rare, here("data", "outputs", "field_comm_rare.csv"))