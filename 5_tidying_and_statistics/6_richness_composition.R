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
library(vegan) #rrarefy and adonis
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
# Field presence-absence composition analysis ####
############################

#This PERMANOVA is run like a GLMM, but we need to first verify that all samples
#have greater than zero species presence in them, or the model will not work:

comp <- field %>%
  group_by(sample) %>%
  filter(sum(reads) >0) %>% #remove samples with zero prey
  ungroup() %>%
  mutate(unique_ID = as.factor(unique_ID), presence = ifelse(reads > 0, 1, 0))
#make unique_ID a factor and create a presence column

#This GLMM is run by saying, how does the fixed effect of sterilization impact
#presence, with a random effects structure with both a random interept term
#for unique_ID (let each species have a different intercept) and a random slopes 
#term for surface sterilization treatment (let each species' relationship with
#with surface sterilization differ, ie let some species increase with surface 
#sterilization, and others decrease)

comp_mod <- glmmTMB(presence ~ Sterilized + (1+Sterilized|unique_ID), 
                       data = comp1,
                       family = "binomial")

comp_null <- glmmTMB(presence ~ 1 + (1|unique_ID), 
                           data = comp1,
                           family = "binomial")

AICc(comp_mod, comp_null)

summary(comp_null)
summary(comp_mod) #sterilization term not significant

#model diagnostics look good
binned_residuals(comp_null)
simulationOutput <- simulateResiduals(fittedModel = comp_null)
fit <- plot(simulationOutput, asFactor=TRUE)

###########################
# adonis instead of GLMM for field presence ####
############################

#An alternative to the above approach is to use ADONIS from the vegan package instead

#Adonis requires a matrix-type object with species as column names and 
#samples/sites as rows, which the following pipe does for the field df
#sites with zero across break adonis, so we need to remove them first
comp1 <- field %>%
  group_by(sample) %>%
  filter(sum(reads) >0) %>%
  ungroup() %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>% #make a presence-absence value
  dplyr::select(presence, sample, unique_ID) %>% #select only the variables for matrix
  pivot_wider(names_from = unique_ID, values_from = presence) %>% 
  #make column names based on unique ID, values from presence in cells
  column_to_rownames(var = "sample") 
#set the column names to the sample, so that only numeric values are in the matrix

#adonis requires a matrix as input, so we will convert df to that
comp1 <- data.matrix(comp, rownames.force = TRUE)

#Now that those are removed, we need to remove them from the eventual metadataframe
#we are creating next, so I will create a dataframe of those sample names to subset
#that dataframe
samples <- as.data.frame(rownames(comp1))

#Adonis also requires a meta-data DF for the analysis based on the factors of interest
#for the analysis. In this analysis, we are interested in comparing surface sterilized
#to unsterilized, so we can subset these from our DF as well

meta_field <- field %>%
  dplyr::select(sample, Sterilized) %>% #select metadata on sterilization and sample
  distinct(sample, Sterilized) %>% #reduce long DF to be just one row per sample
  semi_join(samples, by = c("sample" = "rownames(comp1)")) #subset just those samples
#with non-zero values

#Now we can run adonis, which asks whether there are differences in community
#composition of diet species between the two treatment levels of sterilization
#dist = jaccard with binary = TRUE because this is presence-absence data
adonis(comp1 ~ Sterilized, data = meta_field, dist = "jaccard", binary = TRUE)
#there is no significnat difference between treatment groups 

