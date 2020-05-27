###########################
# Richness Analyses ####
# Ana Miller-ter Kuile
# May 27, 2020 
############################

# This script is looks at species richness for field predators, and then
#looks at species composition of this richness using a presence-absence
#based PERMANOVA
# Ending is additional analyses based on abundance for field, and including
#presence-absence and abundance-based PERMANOVAs for the mesocosm spiders as well

###########################
# Required packages ####
############################

library(here) #tidy data
library(tidyverse) #tidy data
library(ggplot2) #visualize data
library(effects) #dotplots and stuff
#library(lattice) #again dotplots
library(DHARMa) #model diagnostics
library(MuMIn) #model diagnostics
library(glmmTMB) #mixed models
#library(lme4) #mixed models
#library(ggeffects) #ggpredict function
library(performance) #for binomail model fit
#library(see) #for binomial model fit
#library(coin) #wilcoxon test
#library(MASS) #glm.nb
#library(vegan) #rrarefy
library(cowplot) #plot grid at end
#library(picante) #phylogenetic analyses
#library(ape) #phylogenetic analyses
#library(phytools) #phylogenetic analyses
#library(phyloseq) #phylogenetic analyses
#library(hciR) #for as_matrix function
library(esc) #effect sizes
#library(effsize) #alternative effect size package

###########################
# Load data ####
############################

field <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", "field_prey_rare.csv"))

###########################
# Richness analysis ####
############################

richness <- field %>%
  group_by(sample, Sterilized) %>%
  summarise(SR = sum(reads > 0, na.rm=TRUE))

rich_mod <- glmmTMB(SR ~ Sterilized,
                    data = richness, 
                    family = poisson,
                    REML = FALSE)

rich_null <- glmmTMB(SR ~ 1,
                    data = richness, 
                    family = poisson, 
                    REML = FALSE)

AICc(rich_mod, rich_null)

rich_null <- glmmTMB(SR ~ 1,
                     data = richness, 
                     family = poisson)

plot(residuals(rich_null))
summary(rich_mod)

###5/27 START FROM HERE####

simulationOutput_rich <- simulateResiduals(fittedModel = rich_ster_null) 
rich_fit <- plot(simulationOutput_rich, asFactor=TRUE)
SRzi <- testZeroInflation(simulationOutput_rich) 
SRod <- testDispersion(simulationOutput_rich) 
