###########################
# FIeld Richness Analyses of ASVs ####
# Ana Miller-ter Kuile
# May 27, 2020 
############################

# This script is looks at ASV richness for field predators, and then
#looks at species composition of this richness using  CCA analysis

# Ending is additional analyses based on abundance for field, and including
#presence-absence and abundance-based PERMANOVAs for the mesocosm spiders as well

###########################
# Required packages ####
############################

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "MuMIn", 
                  "DHARMa", "performance",
                  "vegan", "remotes",
                  "eulerr", "dummies",
                  "esc", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###########################
# Load data ####
############################

field <- read.csv(here("data",
                       "outputs", 
                       "rarefied_taxonomic_sort", 
                       "field_prey_rare.csv"))


ASV <- read.csv(here("data", 
                     "outputs", 
                     "rarefied_taxonomic_sort", 
                     "field_prey_ASVs_rare.csv"))


###########################
# Field Richness analysis ####
############################

richness <- field %>%
  group_by(sample, Sterilized) %>%
  summarise(SR = sum(reads > 0, na.rm=TRUE)) #number of species with greater than 0 reads

#create a model with and without surface sterilization effect, with poisson
#distribution since these are species count data
#set REML to false so that we can do AICc comparision
rich_mod <- glmmTMB(SR ~ Sterilized,
                    data = richness, 
                    family = poisson,
                    REML = FALSE)

rich_null <- glmmTMB(SR ~ 1,
                     data = richness, 
                     family = poisson, 
                     REML = FALSE)

#compare them with AIC, showing that the null is a better fit
AICc(rich_mod, rich_null)

#re-fit with REML
rich_null <- glmmTMB(SR ~ 1,
                     data = richness, 
                     family = poisson)

rich_mod <- glmmTMB(SR ~ Sterilized,
                    data = richness, 
                    family = poisson)

#based on this summary, surface sterilizatoin treatment is 
#non-significant. 
summary(rich_mod)

summary(rich_null)
#these diagnostics look fine!
plot(residuals(rich_null))
simulationOutput <- simulateResiduals(fittedModel = rich_null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput) 

###########################
# Plot for paper of richness ####
############################
(rich_graph <- ggplot(richness, aes(x = Sterilized, y = SR)) +
   geom_boxplot(fill = "#F29979") + theme_bw() +
   ylim(0,7) +
   labs(x = "Surface sterilization treatment", y = "Prey family richness") +
   scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
   theme(legend.position = "none",
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 15)))

richness %>%
  ungroup() %>%
  summarise(mean= mean(SR), total = n(), sd = sd(SR), se = sd/sqrt(total))
max(richness$SR)
min(richness$SR)

###########################
# Field ASV Richness analysis ####
############################

richness_ASV <- ASV %>%
  group_by(sample, Sterilized) %>%
  summarise(SR = sum(reads > 0, na.rm=TRUE)) #number of ASVs with greater than 0 reads

richness_ASV %>%
  ungroup() %>%
  summarise(rich = mean(SR),
            max = max(SR),
            se = sd(SR)/sqrt(n()))
#create a model with and without surface sterilization effect, with poisson
#distribution since these are species count data
#set REML to false so that we can do AICc comparision
rich_mod <- glmmTMB(SR ~ Sterilized,
                    data = richness_ASV, 
                    family = poisson,
                    REML = FALSE)

rich_null <- glmmTMB(SR ~ 1,
                     data = richness_ASV, 
                     family = poisson, 
                     REML = FALSE)

#compare them with AIC, showing that the null is a better fit
AICc(rich_mod, rich_null)

#re-fit with REML
rich_null <- glmmTMB(SR ~ 1,
                     data = richness, 
                     family = poisson)

rich_mod <- glmmTMB(SR ~ Sterilized,
                    data = richness, 
                    family = poisson)

#based on this summary, surface sterilizatoin treatment is 
#non-significant. 
summary(rich_mod)

summary(rich_null)
#these diagnostics look fine!
plot(residuals(rich_null))
simulationOutput <- simulateResiduals(fittedModel = rich_null, plot= T) 
testZeroInflation(simulationOutput) 
testDispersion(simulationOutput) 

###########################
# Plot for paper of richness ####
############################
(ASV_rich_graph <- ggplot(richness_ASV, aes(x = Sterilized, y = SR)) +
   geom_boxplot(fill = "#F29979") + theme_bw() +
   labs(x = "Surface sterilization treatment", y = "Prey ASV richness") +
   scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
   theme(legend.position = "none",
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 15)))

richness_ASV %>%
  ungroup() %>%
  summarise(mean= mean(SR), total = n(), sd = sd(SR), se = sd/sqrt(total))
max(richness_ASV$SR)
min(richness_ASV$SR)

# Plot Grid ---------------------------------------------------------------

rich_plot <- rich_graph + ASV_rich_graph +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = c('a'), tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 20, vjust = 2))

rich_plot
