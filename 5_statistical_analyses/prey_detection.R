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
library(lattice) #again dotplots
library(DHARMa) #model diagnostics
library(MuMIn) #model diagnostics
library(glmmTMB) #mixed models
library(lme4) #mixed models
library(performance) #for binomail model fit
library(see) #for binomial model fit
library(MASS) #glm.nb
library(cowplot) #plot grid at end

###########################
# Load datasets ####
############################
lab_detect <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", "lab_known_prey_rare.csv"))

field_detect <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", "field_prey_rare.csv"))

