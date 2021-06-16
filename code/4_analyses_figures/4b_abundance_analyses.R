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
library(patchwork)
###########################
# Load datasets ####
############################

#these are all by potential prey species or all predator sequences, so I will
#want to summarise them to be just a value per sample of read abundance
#of each type, removing any samples with zero values? I guess? maybe not at 
#first, i'll have to think about it.

lab_prey <- read.csv(here("data", 
                          "outputs", 
                          "2c_rarefied_taxonomic_sort", 
                          "lab_known_prey_rare.csv"))

lab_prey <- lab_prey %>%
  group_by(sample, Sterilized) %>%
  summarise(known = sum(reads))

lab_all_prey <- read.csv(here("data", 
                              "outputs", 
                              "2c_rarefied_taxonomic_sort", 
                              "lab_all_prey_rare.csv"))

lab_all_prey <- lab_all_prey %>%
  group_by(sample, Sterilized) %>%
  summarise(prey = sum(reads))

lab_pred <- read.csv(here("data", 
                          "outputs", 
                          "2c_rarefied_taxonomic_sort", 
                          "lab_predator_rare.csv"))

lab_pred <- lab_pred %>%
  group_by(sample, Sterilized) %>%
  summarise(pred = sum(reads))

lab_nondiet <- read.csv(here("data", 
                             "outputs", 
                             "2c_rarefied_taxonomic_sort", 
                             "lab_nondiet_rare.csv"))

lab_nondiet <- lab_nondiet %>%
  group_by(sample, Sterilized) %>%
  summarise(nondiet = sum(reads))

lab_all <- lab_prey %>%
  left_join(lab_all_prey, by = c("sample", "Sterilized")) %>%
  left_join(lab_pred, by = c("sample", "Sterilized")) %>%
  left_join(lab_nondiet, by = c("sample", "Sterilized")) %>%
  mutate(total = 55205) 

fld_prey <- read.csv(here("data", 
                          "outputs", 
                          "2c_rarefied_taxonomic_sort", 
                          "field_prey_rare.csv"))

fld_prey <- fld_prey %>%
  group_by(sample, Sterilized) %>%
  summarise(prey = sum(reads))

fld_pred <- read.csv(here("data", 
                          "outputs", 
                          "2c_rarefied_taxonomic_sort", 
                          "field_predator_rare.csv"))

fld_pred <- fld_pred %>%
  group_by(sample, Sterilized) %>%
  summarise(pred = sum(reads))

fld_nondiet <- read.csv(here("data", 
                             "outputs", 
                             "2c_rarefied_taxonomic_sort", 
                             "field_nondiet_rare.csv"))

fld_nondiet <- fld_nondiet %>%
  group_by(sample, Sterilized) %>%
  summarise(nondiet = sum(reads))

fld_all <- fld_prey %>%
  left_join(fld_pred, by = c("sample", "Sterilized")) %>%
  left_join(fld_nondiet, by = c("sample", "Sterilized")) %>%
  mutate(total = 16004)


# Read totals summaries ---------------------------------------------------

lab_all %>%
  dplyr::select(-known) %>%
  full_join(fld_all) %>% 
  summarise(prey_prop = prey/total,
            pred_prop = pred/total,
            non_prop = nondiet/total) %>% 
  summarise(mean_prey = mean(prey_prop),
            se_prey = sd(prey_prop)/sqrt(n()),
            mean_pred = mean(pred_prop),
            se_pred = sd(pred_prop)/sqrt(n()),
            mean_non = mean(non_prop),
            se_non = sd(non_prop)/sqrt(n()))

###########################
# Lab model ####
############################

#For all models, I'm going to subset just samples that had some prey
#items detected, as what we are asking here is whether some aspect
#of the surface sterilization treatment altered DNA abundance 
#proportions/composition when prey is detected
lab_all %>%
  filter(known == 0) %>%
  ungroup() %>%
  tally()

lab_all %>%
  filter(known != 0) %>%
  ungroup() %>%
  tally()

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
fld_all %>%
  filter(prey == 0) %>%
  ungroup() %>%
  tally()

fld_all %>%
  filter(prey > 0) %>%
  ungroup() %>%
  tally()

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
# Visualizations ####
############################
#let's visualize this then
my_y_title_2 <- expression(paste("Percent ", italic("O. japonica"), " reads per sample"))

(b_2 <- ggplot(lab_all_nz, aes(x = Sterilized, y = (known/total)*100)) +
  geom_boxplot(fill = "#7B8D65") +
  scale_y_log10() +
  theme_bw() + 
  labs(x = "Surface sterilization treatment", y = my_y_title_2) +
  scale_x_discrete(labels=c("NS" = "Not S. Sterilized", 
                              "SS" = "Surface Sterilized")) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12)))

(e_2 <- ggplot(fld_all_nz, aes(x = Sterilized, y = (prey/total)*100)) +
  geom_boxplot(fill = "#F29979") +
  scale_y_log10()+
  labs(x = "Surface sterilization treatment", y = "Percent prey reads per sample") +
  theme_bw() + 
  scale_x_discrete(labels=c("NS" = "Not S. Sterilized", 
                            "SS" = "Surface Sterilized")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15)))

###########################
# Plot grid and summaries for paper ####
############################

plot_grid(b_2, e_2, nrow = 2, align = "vh")

abund_plot <- e_2 + b_2 +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = c('a'), tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 20, vjust = 2))

lab_all_nz %>%
  summarise(prop = known/total) %>%
  summarise(mean = mean(prop), total= n(), sd = sd(prop), se = sd/sqrt(total))

max(lab_all_nz$prey)/55205
min(lab_all_nz$prey)/55205

max(fld_all_nz$prey)/16004
min(fld_all_nz$prey)/16004

fld_all_nz %>%
  summarise(prop = prey/total) %>%
  summarise(mean = mean(prop), total= n(), sd = sd(prop), se = sd/sqrt(total))

###########################
# Supplement: Lab Predator model ####
############################

#These predator models are not nice...
#Basically, what I'm trying to show here is that there is no "relic" that would
#be driving the marginal difference in the mesocosm predator known prey reduction
#but rather this represents some effect of surface contamination. So yeah...

lab_pred_mod <- glmmTMB(pred ~ Sterilized,
                   data = lab_all,
                   REML = FALSE,
                   family = "gaussian")

lab_pred_null <- glmmTMB(pred ~ 1,
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
# Supplement: Field Predator model ####
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


###########################
# Supplement: Lab all Prey model ####
############################

lab_all_prey <- lab_all %>%
  filter(prey > 0)

lab_all_mod <- glmmTMB(prey ~ Sterilized,
                   data = lab_all_prey,
                   offset = log(total),
                   REML = FALSE,
                   family = "genpois")

lab_all_null <- glmmTMB(prey ~ 1,
                    data = lab_all_prey,
                    offset = log(total),
                    REML = FALSE,
                    family = "genpois")

AICc(lab_mod, lab_null)

lab_all_mod <- glmmTMB(prey ~ Sterilized,
                   data = lab_all_prey,
                   offset = log(total),
                   family = "genpois")

lab_all_null <- glmmTMB(prey ~ 1,
                    data = lab_all_prey,
                    offset = log(total),
                    family = "genpois")

summary(lab_all_mod)

plot(allEffects(lab_all_mod))

plot(residuals(lab_all_null))
simulationOutput <- simulateResiduals(fittedModel = lab_all_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput)

###########################
# Supplement: Lab nondiet models ####
############################
lab_nd_mod <- glmmTMB(nondiet ~ Sterilized,
                       data = lab_all,
                       offset = log(total),
                       REML = FALSE,
                       family = "genpois")

lab_nd_null <- glmmTMB(nondiet ~ 1,
                        data = lab_all,
                        offset = log(total),
                        REML = FALSE,
                        family = "genpois")

AICc(lab_nd_mod, lab_nd_null)

lab_nd_mod <- glmmTMB(nondiet ~ Sterilized,
                      data = lab_all,
                      offset = log(total),
                      family = "genpois")

lab_nd_null <- glmmTMB(nondiet ~ 1,
                       data = lab_all,
                       offset = log(total),
                       family = "genpois")

summary(lab_nd_mod)

plot(allEffects(lab_nd_mod))

plot(residuals(lab_nd_null))
simulationOutput <- simulateResiduals(fittedModel = lab_nd_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput)

###########################
# Supplement: Field nondiet models ####
############################
fld_nd_mod <- glmmTMB(nondiet ~ Sterilized,
                      data = fld_all,
                      offset = log(total),
                      REML = FALSE,
                      family = "genpois")

fld_nd_null <- glmmTMB(nondiet ~ 1,
                       data = fld_all,
                       offset = log(total),
                       REML = FALSE,
                       family = "genpois")

AICc(fld_nd_mod, fld_nd_null)

fld_nd_mod <- glmmTMB(nondiet ~ Sterilized,
                      data = fld_all,
                      offset = log(total),
                      family = "genpois")

fld_nd_null <- glmmTMB(nondiet ~ 1,
                       data = fld_all,
                       offset = log(total),
                       family = "genpois")

summary(fld_nd_mod)

plot(allEffects(fld_nd_mod))

plot(residuals(fld_nd_null))
simulationOutput <- simulateResiduals(fittedModel = fld_nd_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput)

###########################
# Supplement: Predator and other prey graphs ####
############################
(lab_pred_graph <- ggplot(lab_all, aes(x = Sterilized, y = pred)) +
   geom_boxplot(fill = "#7B8D65") + theme_bw() + 
    labs(x = "Surface sterilization treatment", y = "Consumer DNA reads") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) + 
    theme(legend.position = "none")
  )

(fld_pred_graph <- ggplot(fld_all, aes(x = Sterilized, y = pred)) +
    geom_boxplot(fill = "#F29979") + theme_bw() + 
    labs(x = "Surface sterilization treatment", y = "Consumer DNA reads") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) + 
    theme(legend.position = "none")
)

(lab_prey_graph <- ggplot(lab_all_prey, aes(x = Sterilized, y = prey)) +
    geom_boxplot(fill = "#7B8D65") + theme_bw() + 
    scale_y_log10() +
    labs(x = "Surface sterilization treatment", y = "Potential diet DNA reads") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) + 
    theme(legend.position = "none")
)

(fld_prey_graph <- ggplot(fld_all_nz, aes(x = Sterilized, y = prey)) +
    geom_boxplot(fill = "#F29979") +
    scale_y_log10()+
    labs(x = "Surface sterilization treatment", y = "Potential diet DNA reads") +
    theme_bw() + theme(legend.position = "none") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")))

(lab_nd_graph <- ggplot(lab_all, aes(x = Sterilized, y = nondiet)) +
  geom_boxplot(fill = "#7B8D65") + scale_y_log10() +
    theme_bw() + 
    labs(x = "Surface sterilization treatment", y = "Nondiet DNA reads") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) + 
    theme(legend.position = "none"))

(fld_nd_graph <- ggplot(fld_all, aes(x = Sterilized, y = nondiet)) +
  geom_boxplot(fill = "#F29979") +scale_y_log10() +
    labs(x = "Surface sterilization treatment", y = "Nondiet DNA reads") +
    theme_bw() + theme(legend.position = "none") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")))
    
plot_grid(lab_pred_graph, fld_pred_graph, lab_prey_graph, fld_prey_graph, 
          lab_nd_graph, fld_nd_graph, ncol = 2, align = "vh")
