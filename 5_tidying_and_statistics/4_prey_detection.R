###########################
# Detection analyses ####
# Ana Miller-ter Kuile
# May 27, 2020 
############################

# This script is the detection analyses of potential ingestion of 
# the offered prey (mesocosm) and the wild prey (field) predators
# assessing whether surface contamination alters detection of these
#prey items in these environments.

###########################
# Required packages ####
############################

library(here) #tidy data
library(tidyverse) #tidy data
library(ggplot2) #visualize data
library(effects) #dotplots and allEffects plots
library(DHARMa) #model diagnostics
library(MuMIn) #model diagnostics
library(glmmTMB) #mixed models
library(performance) #for binomial model fit binned residuals
library(cowplot) #plot grid at end

###########################
# Load datasets ####
############################

#the lab dataframe is already in the format it needs to be for analysis, with
#the addition of a presence-absence column
lab_detect <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", "lab_known_prey_rare.csv"))

#this field dataframe is currently by-species (so it can be used in PERMANOVA later),
#so we will need to summarise it a bit further to be able to do presence-absence analyses
#of all prey
field_detect <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", "field_prey_rare.csv"))

###########################
# Create presence-absence column in both DFs ####
############################
lab_detect <- lab_detect %>%
  mutate(presence = ifelse(reads > 0, 1, 0))

field_detect <- field_detect %>%
  group_by(sample, Sterilized) %>% #group by sample, and also sterilized to keep in DF
  summarise(reads = sum(reads)) %>% #summarize total prey reads per sample
  mutate(presence = ifelse(reads > 0, 1, 0)) #set a binary presence-absence column
  
###########################
# Mesocosm Known prey detection ####
############################           

lab_detect_mod <- glmmTMB(presence ~ Sterilized,
                        data = lab_detect, 
                        family = binomial)

lab_null_model <- glmmTMB(presence ~ 1,
                         data = lab_detect,
                         family = binomial)

AICc(lab_detect_mod, lab_null_model)

#The sterilized term is marginally significant at p-value 0.07
summary(lab_detect_mod)

#assess model fits, which look pretty good:
simulationOutput_lab <- simulateResiduals(fittedModel = lab_detect_mod) 
fit_lab <- plot(simulationOutput_lab, asFactor=TRUE)
binned_residuals(lab_detect_mod)

#Surface sterilization slightly reduces detection of the offered 
#potential prey item
plot(allEffects(lab_detect_mod))

###########################
# Field prey detection ####
############################           

fld_detect_mod <- glmmTMB(presence ~ Sterilized,
                          data = field_detect, 
                          family = binomial)

fld_null_model <- glmmTMB(presence ~ 1,
                          data = field_detect,
                          family = binomial)

AICc(fld_detect_mod, fld_null_model)

#The sterilized term is clearly non-significant
summary(fld_detect_mod)
summary(fld_null_model)
#assess model fits, which look pretty good:
simulationOutput_fld <- simulateResiduals(fittedModel = fld_null_model) 
fit_fld <- plot(simulationOutput_fld, asFactor=TRUE)
binned_residuals(fld_null_model)

###########################
# Visualizations and summaries for paper ####
############################

#Mesocosms visualations
ASVs_known <- lab_detect %>%
  group_by(Sterilized) %>%
  summarise(total = n(), presence = sum(presence), detection = presence/total)

my_y_title <- expression(paste(italic("O. japonica"), " ASV detection (%)"))

(a <- ggplot(ASVs_known, aes(x = Sterilized, y = detection*100)) +
    geom_bar(stat="identity", position= "dodge", fill = "#7B8D65", color = "black") +theme_bw() +
    labs(y = my_y_title) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()))

#Mesocosm summaries
#by sterilization detection:
ASVs_known

#overall detection:
lab_detect %>%
  summarize(total = n(), present = sum(presence), detection = present/total)

#Field visualization
prey_presence <- field_detect %>%
  group_by(Sterilized) %>%
  summarise(total = n(), presence = sum(presence), detection = presence/total)

(c <- ggplot(prey_presence, aes(x = Sterilized, y = detection*100)) +
    geom_bar(stat="identity", position= "dodge", fill = "#F29979", color = "black") +theme_bw() +
    labs(x = "Surface sterilization treatment", y = "Prey ASV detection (%)") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
    theme(legend.position = "none"))

#Field summaries
#by sterilization
prey_presence

#overall detection:
field_detect %>%
  ungroup() %>%
  summarise(total = n(), presence = sum(presence), detection = presence/total)

###########################
# Plot to export ####
############################

plot_grid(a, c, nrow= 2, align = "hv")
  
