######################
#BBSplit Reference Creation####
#February 3, 2020
#Ana Miller-ter Kuile
######################

#This R code takes ASV output from dada2 on uncleaned sequences and creates reference
#databases to use in BBSplit to remove predator sequences from prey sequences. The goal
#of this is to see if this increases the number of prey ASVs in the final dataset

#This code requires
#1. The ASV file output from the unclean dada2pipeline (written file of "asv_fasta" object)
#2. Taxonomies matched to that ASV file from MEGAN (will have to BLAST the file against the 
#nucleotide database on GenBank first, and export that file to MEGAN, and then export with
#File -> Export -> Text (CSV) format -> readName_to_taxonName with summarized reads). I manually
#added higher taxonomies to this file since MEGAN v6 only outputs the last of this path to
#this file

######################
#Required packages####
######################
library(here)
library(tidyverse)

######################
#Load and tidy data####
######################
#import dada2 ASVs
ASV_dada <- read_delim(here("data",
                            "data_raw",
                            "1_denoised_data",
                            "ASV_lists",
                            "dada2_uc_asvs.fasta"), delim="\r", col_names = F)

#change formatting so that it is two columns - 1 with ASV, one with sequence
ASV_dada <- lapply(ASV_dada[,1], function(x) {
  odd <- seq_along(x) %% 2 == 1
  o <- x[odd]
  e <- x[!odd]
  length(e) <- length(o)
  data.frame(o, e)
})[[1]]; colnames(ASV_dada) <- c('Label', 'Text')
ASV_dada[,1] <- as.character(ASV_dada[,1]); ASV_dada[,2] <- as.character(ASV_dada[,2])

######################
#Load taxonomic assignments####
######################
#the file of IDs that match different taxonomies
#I manually entered higher taxonomies in this file since MEGAN v6 only gives the 
#last part of the taxon path now
ASV_id_NCBI <- read.csv(here("3_taxonomies", 
                             "MEGAN", 
                             "dada2_unclean", 
                             "dada2_uc_ncbi_taxa.csv"))
ASV_id_BOLD <- read.csv(here("3_taxonomies", 
                             "BOLD", 
                             "dada2_uc", 
                             "dada2_uc_taxa_bold.csv"))
ASV_id_BOLD <- dplyr::rename(ASV_id_BOLD, "ASV" = "Query.ID")

bold_id <- ASV_id_BOLD %>%
  dplyr::select(ASV, ID)

ncbi_id <- ASV_id_NCBI %>%
  dplyr::select(ASV, ID)

id_both <- bold_id %>%  
full_join(ncbi_id, by = "ASV") #combine assignments from MEGAN and BOLD

######################
#Sort predator from prey####
######################
#filters out just the predator ASVs from the ID file
predator <- id_both %>%
  filter(ID.x == "Heteropoda venatoria" | ID.y == "Sparassidae") %>%
  select("ASV")


#you may need to change the ASV name to include a > here, I also added
#it in text edit at some point, so two options here. 
#change name so that it matches that in ASV_dada
predator$ASV <- paste(">", predator$ASV, sep="")
predator$ASV <- as.character(predator$ASV)

#filters out just the prey ASVs from the ID file and matches name
prey <- id_both %>%
  filter(ID.x != "Heteropoda venatoria" | ID.y != "Sparassidae") %>%
  select("ASV")

prey$ASV <- paste(">", prey$ASV, sep="")
prey$ASV <- as.character(prey$ASV)

#subsets ASV_dada for just predator
predator_ASV <- ASV_dada %>%
  semi_join(predator, by = c("Label" = "ASV"))

#subsets ASV_dada for just prey
prey_ASV <- ASV_dada %>%
  semi_join(prey, by = c("Label" = "ASV"))

######################
#Write files to export and use in BBSPlit####
######################
#write predator and prey reference files
write.table(prey_ASV, file = here("data", 
                                  "outputs",
                                  "1_bioinformatics",
                                  "bbsplit_references",
                                  "prey_ref.fasta"), 
            row.names=F, quote=F, col.names=F)

write.table(predator_ASV, file = here("data",
                                      "outputs",
                                      "1_bioinformatics",
                                      "bbsplit_references",
                                      "predator_ref.fasta"), 
            row.names=F, quote=F, col.names=F)









