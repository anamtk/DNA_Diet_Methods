###########################
# Abundance Analyses ####
# Ana Miller-ter Kuile
# June 1, 2020
############################

#This script looks at DNA abundance in samples, specifically examining
#whether surface sterilization altered any other types of DNA (e.g.
#consumer or other prey(for mesocosm consumers)), suggesting that 
#surface sterilization was altered samples in ways OTHER than through
#the removal of contaminants (e.g. in ways that might suggest that 
#surface sterilizaiton is in some way degrading the DNA in samples)

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

#these are all by potential prey species or all predator sequences, so I will
#want to summarise them to be just a value per sample of read abundance
#of each type, removing any samples with zero values? I guess? maybe not at 
#first, i'll have to think about it.

lab_prey <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", 
                          "lab_known_prey_rare.csv"))

lab_prey <- lab_prey %>%
  group_by(sample, Sterilized) %>%
  summarise(known = sum(reads))

lab_all_prey <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", 
                          "lab_all_prey_rare.csv"))

lab_all_prey <- lab_all_prey %>%
  group_by(sample, Sterilized) %>%
  summarise(prey = sum(reads))

lab_pred <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", 
                              "lab_predator_rare.csv"))

lab_pred <- lab_pred %>%
  group_by(sample, Sterilized) %>%
  summarise(pred = sum(reads))

lab_all <- lab_prey %>%
  left_join(lab_all_prey, by = c("sample", "Sterilized")) %>%
  left_join(lab_pred, by = c("sample", "Sterilized")) %>%
  mutate(total = 55205) 

fld_prey <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", 
                          "field_prey_rare.csv"))

fld_prey <- fld_prey %>%
  group_by(sample, Sterilized) %>%
  summarise(prey = sum(reads))

fld_pred <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", 
                          "field_predator_rare.csv"))

fld_pred <- fld_pred %>%
  group_by(sample, Sterilized) %>%
  summarise(pred = sum(reads))

fld_all <- fld_prey %>%
  left_join(fld_pred, by = c("sample", "Sterilized")) %>%
  mutate(total = 16004)

###########################
# Lab model ####
############################

#For all models, I'm going to subset just samples that had some prey
#items detected, as what we are asking here is whether some aspect
#of the surface sterilization treatment altered DNA abundance 
#proportions/composition when prey is detected

lab_all_nz <- lab_all %>%
  filter(known > 0)

lab_mod <- glmmTMB(known ~ Sterilized,
                   data = lab_all_nz,
                   offset = log(total),
                   REML = FALSE,
                   family = "genpois")

lab_null <- glmmTMB(known ~ 1,
                   data = lab_all_nz,
                   offset = log(total),
                   REML = FALSE,
                   family = "genpois")

AICc(lab_mod, lab_null)

lab_mod <- glmmTMB(known ~ Sterilized,
                   data = lab_all_nz,
                   offset = log(total),
                   family = "genpois")

lab_null <- glmmTMB(known ~ 1,
                    data = lab_all_nz,
                    offset = log(total),
                    family = "genpois")

summary(lab_mod)

plot(allEffects(lab_mod))

plot(residuals(lab_null))
simulationOutput <- simulateResiduals(fittedModel = lab_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput)

###########################
# Field model ####
############################

fld_all_nz <- fld_all %>%
  filter(prey > 0)

fld_mod <- glmmTMB(prey ~ Sterilized,
                   data = fld_all_nz,
                   offset = log(total),
                   REML = FALSE,
                   family = "genpois")

fld_null <- glmmTMB(prey ~ 1,
                   data = fld_all_nz,
                   offset = log(total),
                   REML = FALSE,
                   family = "genpois")

AIC(fld_mod, fld_null)

fld_mod <- glmmTMB(prey ~ Sterilized,
                   data = fld_all_nz,
                   offset = log(total),
                   family = "genpois")

fld_null <- glmmTMB(prey ~ 1,
                    data = fld_all_nz,
                    offset = log(total),
                    family = "genpois")

summary(fld_mod)
plot(allEffects(fld_mod))

plot(residuals(fld_null))
simulationOutput <- simulateResiduals(fittedModel = fld_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput)

###########################
# Lab Predator model ####
############################

#These models are not nice...

lab_pred_mod <- glmmTMB(log10(pred) ~ Sterilized,
                   data = lab_all,
                   REML = FALSE,
                   family = "gaussian")

lab_pred_null <- glmmTMB(log10(pred) ~ 1,
                    data = lab_all,
                    REML = FALSE,
                    family = "gaussian")

AICc(lab_pred_mod, lab_pred_null)

lab_pred_mod <- glmmTMB(pred ~ Sterilized,
                        data = lab_all,
                        family = "gaussian")

lab_pred_null <- glmmTMB(pred ~ 1,
                         data = lab_all,
                         family = "gaussian")

summary(lab_pred_mod)
plot(allEffects(lab_pred_mod))

plot(residuals(lab_pred_null))
simulationOutput <- simulateResiduals(fittedModel = lab_pred_null) 
fit <- plot(simulationOutput, asFactor=TRUE)


###########################
# Field Predator model ####
############################
fld_pred_mod <- glmmTMB(pred ~ Sterilized,
                        data = fld_all,
                        offset = log(total),
                        REML = FALSE,
                        family = "genpois")

fld_pred_null <- glmmTMB(pred ~ 1,
                         data = fld_all,
                         offset = log(total),
                         REML = FALSE,
                         family = "genpois")

AICc(fld_pred_mod, fld_pred_null)

fld_pred_mod <- glmmTMB(pred ~ Sterilized,
                        data = fld_all,
                        offset = log(total),
                        family = "genpois")

fld_pred_null <- glmmTMB(pred ~ 1,
                         data = fld_all,
                         offset = log(total),
                         family = "genpois")

summary(fld_pred_mod)
plot(allEffects(fld_pred_mod))

plot(residuals(fld_pred_null))
simulationOutput <- simulateResiduals(fittedModel = fld_pred_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput)

ggplot(fld_all, aes(x = Sterilized, y = prey/total)) +
  geom_boxplot() + theme_bw() +scale_y_log10()

ggplot(lab_all, aes(x = Sterilized, y = known/total)) +
  geom_boxplot() + theme_bw() +scale_y_log10()

ggplot(fld_all, aes(x = Sterilized, y = pred/total)) +
  geom_boxplot() + theme_bw() 

ggplot(lab_all, aes(x = Sterilized, y = pred/total)) +
  geom_boxplot() + theme_bw() 
