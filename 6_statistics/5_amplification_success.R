###########################
# Amplification success####
# Ana Miller-ter Kuile
# June 1, 2020
############################

#just stats on amplification success for the results section of paper
library(here)
library(tidyverse)
amp <- read.csv(here("data", "Amplification_Success.csv"))
metadata <- read.csv(here("data", "Sample_Metadata.csv"))

amp %>%
  left_join(metadata, by = c("sample", "Source", "Sterilized")) %>%
  group_by(Source, Sterilized) %>%
  summarise(total = n(), proportion = sum(Amp_success/n()), success = sum(Amp_success >0))

amp %>%
  summarise(total = n(), proportion = sum(Amp_success/n()), success = sum(Amp_success >0))
