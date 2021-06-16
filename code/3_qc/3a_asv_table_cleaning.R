############################
#Correcting Error using positives (Austen)####
#Ana Miller-ter Kuile
#June 15, 2020
############################

############################
#Load Packages####
############################
library(here)
library(tidyverse)

############################
#Load Data####
############################

#load the dataframe
u3_error_check <- read.csv(here("data", 
                                "data_raw",
                                "1_denoised_data", 
                                "ASV_tables", 
                                "unoise_uc_zotu_tab.txt"), sep = "\t")
#rename columns for simplicity
colnames(u3_error_check) <- sapply(str_split(colnames(u3_error_check), "S"), function(x){return(x[[1]])})
u3_error_check <- rename(u3_error_check, "ASV" = "X.OTU.ID")

############################
#Error Checking####
############################

#The ASVs assigned to the positive controls need to be corrected, but they also
#look really good, so...

#non-zero reads should only be in the positive control samples:
u3 <- u3_error_check %>%
  filter(ASV %in% c("Zotu4", "Zotu2", "Zotu3")) %>%
  gather(sample, reads, CL1:QC1) %>%
  filter(reads > 0) %>%
  mutate(type = "u3")
