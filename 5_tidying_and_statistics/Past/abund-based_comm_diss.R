#Community matrix analyses
#Ana Miller-ter Kuile
#February 16, 2020

#This code looks at the following questions for each community matrix,
#broken up into sterilized, not sterilized, field, and lab:
#1. is community different (abundance-based) between sterilized and 
#unsterilized individuals? (PERMANOVA)
#1a. follow-up ecological questions here? (i.e. who are community differences,
#are they "prey", "pred", "non-diet", or "no hit"? are there different numbers
#of prey ASVs when the animal is sterilized first?)
#2. how abundant known prey ASVs in lab individuals by pipeline
#and by sterilization?

#categories to sort to:
#Groups field: predator, prey, unknown, no-hit, non-prey
#Groups lab: predator, known-prey, other-prey, unknown, no-hit, non-prey


#Load Packages ####
library(tidyverse)
library(here)
library(ggplot2)
library(ggfortify)
library(GUniFrac)
library(vegan)
library(iNEXT)
library(cowplot)
library(ecodist)
library(lme4)
library(broom)
library(MASS)
library(ggeffects)
library(DHARMa)
library(glmmTMB)

#(Now that I've moved some code to another script, need to consider here if any of the
#dataframes here need to be re-imported again)

#needed for all samples ####
d2_taxa <- read_csv(here("data", "d2_uc_tax_ass.csv"))
u3_taxa <- read_csv(here("data", "u3_uc_tax_ass.csv"))

#Rarefy for abundance-based dissimilarity####
#i'm rarefying by field/lab separately, I *think* this is okay...

#UNOISE FIELD rarefy and mixed model####
f_u3_nmds <- read.csv(here("data", "fld_u3uc_nmds.csv"))
#set rownames to sample variable
rownames(f_u3_nmds) <- f_u3_nmds$sample
#and then select just the numeric values in this dataframe
fld_u3_nmds <- f_u3_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#delete those zero columns since they represent "species" that are 
#not present in the field "communities"
fld_u3_nmds <- fld_u3_nmds[which(colSums(fld_u3_nmds)> 0)] #176 down to 145

#import the metadata for later
meta_u3_field <- read.csv(here("data", "fld_u3_nmds_meta.csv"))

raremin_u3_fld <- min(rowSums(fld_u3_nmds))
raremin_u3_fld #16004, which is substantially lower than some others...
rowSums(fld_u3_nmds)

# Rarefy, which means to standardize all samples by the lowest value (raremin). 
# change sample = to the rarefaction depth (# of reads per sample)
rarefied_u3_fld <- rrarefy(fld_u3_nmds, sample = raremin_u3_fld)
base::range(rowSums(rarefied_u3_fld)) # checks that it worked
rarefied_u3_fld <- as.data.frame(rarefied_u3_fld) # save it as a dataframe

colSums(rarefied_u3_fld) #145 total ZOTUS here
rarefied_u3_fld[which(colSums(rarefied_u3_fld) == 0)] #11 equal zero across all

#transpose dataframe so that we can make it long by sample and 
#attach metadata AND taxonomic assignments of ASVs to it
rarefied_u3_fld_t <- as.data.frame(t(rarefied_u3_fld))
rarefied_u3_fld_t$ASV <- row.names(rarefied_u3_fld_t) #make ASV column for join

#make long by sample, then join to taxonomic assignments.
rarefied_u3_fld_long <- rarefied_u3_fld_t %>%
  gather(sample, reads, HEV87:HEV83) %>%
  left_join(u3_taxa, by = "ASV")

#Groups field: predator, prey, unknown, no-hit, non-prey
rarefied_u3_fld_long$groups <- rarefied_u3_fld_long$taxonomy

rarefied_u3_fld_grp <- rarefied_u3_fld_long %>%
  group_by(sample, groups) %>%
  summarize(reads = sum(reads)) %>%
  left_join(meta_u3_field, by = "sample")

#set species ID to a factor
rarefied_u3_fld_grp$groups <- as.factor(rarefied_u3_fld_grp$groups)
rarefied_u3_fld_grp$Ster <- ifelse(rarefied_u3_fld_grp$Sterilized == "NS", 0, 1)
#model that says that the slope of the response to sterilization can vary
#randomly by Species_ID
#u3_fld_abund_mod <- glmmTMB(reads ~ Ster*groups + (1|sample),
#                            data = rarefied_u3_fld_grp,
#                            family = "nbinom2")
#summary(u3_fld_abund_mod)
hist(rarefied_u3_fld_grp$reads)

rarefied_u3_fld_grp_np <- rarefied_u3_fld_grp[which(rarefied_u3_fld_grp$groups != "predator"),]
rarefied_u3_fld_grp_np$groups <- droplevels(rarefied_u3_fld_grp_np$groups)

#this is our model, which is grouped by sample, saying:
#we are wondering what the effect of sterilization AND taxonomic group
#on reads, and let's say that sterilization may interact with groups
#such that the effect of sterilization may vary depending on group

u3_fld_abund_np_mod <- glmmTMB(reads ~ Ster*groups + (1|sample),
                            data = rarefied_u3_fld_grp_np,
                            family = "genpois",
                            REML = TRUE)

simulationOutput <- simulateResiduals(fittedModel = u3_fld_abund_np_mod) 
plot(simulationOutput, asFactor=TRUE) #these look um
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed
#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = u3_fld_abund_mod) 
#plot(simulationOutput, asFactor=TRUE) #these look um
#testZeroInflation(simulationOutput) #not zero inflated
#testDispersion(simulationOutput) #not overdispersed

#looks at the summary of this model
summary(u3_fld_abund_np_mod) #based on this summary it doesn't 
#seem like sterilization acts on any of these groups yAY
plot(residuals(u3_fld_abund_np_mod))

#this is the null model with just a random intercept by species ID
u3_fld_abund_np_mod2 <- glmmTMB(reads ~ Ster + groups + (1|sample),
                               data = rarefied_u3_fld_grp_np,
                               family = "genpois",
                               REML = TRUE)

u3_fld_abund_np_mod3 <- glmmTMB(reads ~ groups + (1|sample),
                               data = rarefied_u3_fld_grp_np,
                               family = "genpois", 
                               REML = TRUE)

u3_fld_abund_np_mod4 <- glmmTMB(reads ~ Ster + (1|sample),
                                data = rarefied_u3_fld_grp_np,
                                family = "genpois",
                                REML = TRUE)

u3_fld_abund_np_mod_null <- glmmTMB(reads ~ 1 + (1|sample), 
                            data = rarefied_u3_fld_grp_np,
                            family = "genpois", 
                            REML = TRUE)

#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(u3_fld_abund_np_mod, u3_fld_abund_np_mod2, u3_fld_abund_np_mod3, 
      u3_fld_abund_np_mod4, u3_fld_abund_np_mod_null)

u3_fld_abund_np_mod_null_ml <- glmmTMB(reads ~ 1 + (1|sample), 
                                    data = rarefied_u3_fld_grp_np,
                                    family = "genpois")

summary(u3_fld_abund_np_mod_null_ml)
simulationOutput <- simulateResiduals(fittedModel = u3_fld_abund_np_mod_null_ml) 
plot(simulationOutput, asFactor=TRUE) #these look um
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed

#UNOISE LAB RAREFY AND PERMANOVA ####
l_u3_nmds <- read.csv(here("data", "lab_u3uc_nmds.csv"))
#set rownames to sample variable
rownames(l_u3_nmds) <- l_u3_nmds$sample
#and then select just the numeric values in this dataframe
lab_u3_nmds <- l_u3_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#delete those zero columns since they represent "species" that are 
#not present in the field "communities"
lab_u3_nmds <- lab_u3_nmds[which(colSums(lab_u3_nmds)> 0)] #176 down to 75

#import the metadata for later
meta_u3_lab <- read.csv(here("data", "lab_u3uc_nmds_meta.csv"))

raremin_u3_lab <- min(rowSums(lab_u3_nmds))
rowSums(lab_u3_nmds)
raremin_u3_lab #55206, which is about half the max number

# Rarefy, which means to standardize all samples by the lowest value (raremin). 
# change sample = to the rarefaction depth (# of reads per sample)
rarefied_u3_lab <- rrarefy(lab_u3_nmds, sample = raremin_u3_lab)
base::range(rowSums(rarefied_u3_lab)) # checks that it worked
rarefied_u3_lab <- as.data.frame(rarefied_u3_lab) # save it as a dataframe

colSums(rarefied_u3_lab) #145 total ZOTUS here
rarefied_u3_lab[which(colSums(rarefied_u3_lab) == 0)] #0 equal zero across all

#transpose dataframe so that we can make it long by sample and 
#attach metadata AND taxonomic assignments of ASVs to it
rarefied_u3_lab_t <- as.data.frame(t(rarefied_u3_lab))
rarefied_u3_lab_t$ASV <- row.names(rarefied_u3_lab_t) #make ASV column for join

#make long by sample, then join to taxonomic assignments.
rarefied_u3_lab_long <- rarefied_u3_lab_t %>%
  gather(sample, reads, HEV15:HEV29) %>%
  left_join(u3_taxa, by = "ASV")

#Groups field: predator, prey, unknown, no-hit, non-prey
u3_lab_groups <- ifelse(rarefied_u3_lab_long$taxonomy == "prey" & rarefied_u3_lab_long$ID_ncbi == "Acrididae", "known prey",
                        ifelse(rarefied_u3_lab_long$taxonomy == "prey" & rarefied_u3_lab_long$ID_ncbi != "Acrididae", "other prey",
                               rarefied_u3_lab_long$taxonomy))
rarefied_u3_lab_long$groups <- u3_lab_groups

rarefied_u3_lab_grp <- rarefied_u3_lab_long %>%
  group_by(sample, groups) %>%
  summarize(reads = sum(reads)) %>%
  left_join(meta_u3_lab, by = "sample")
str(rarefied_u3_lab_grp)
#set species ID to a factor
rarefied_u3_lab_grp$groups <- as.factor(rarefied_u3_lab_grp$groups)
rarefied_u3_lab_grp$Ster <- ifelse(rarefied_u3_lab_grp$Sterilized == "NS", 0, 1)

rarefied_u3_lab_grp_np <- rarefied_u3_lab_grp[which(rarefied_u3_lab_grp$groups != "predator"),]
rarefied_u3_lab_grp_np$groups <- droplevels(rarefied_u3_lab_grp_np$groups)

#this is our model, which is grouped by sample, saying:
#we are wondering what the effect of sterilization AND taxonomic group
#on reads, and let's say that sterilization may interact with groups
#such that the effect of sterilization may vary depending on group

u3_lab_abund_np_mod <- glmmTMB(reads ~ Ster*groups + (1|sample),
                               data = rarefied_u3_lab_grp_np,
                               family = "nbinom2",
                               REML = TRUE)

simulationOutput <- simulateResiduals(fittedModel = u3_lab_abund_np_mod) 
plot(simulationOutput, asFactor=TRUE) #these look um
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed
#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = u3_fld_abund_mod) 
#plot(simulationOutput, asFactor=TRUE) #these look um
#testZeroInflation(simulationOutput) #not zero inflated
#testDispersion(simulationOutput) #not overdispersed

#looks at the summary of this model
summary(u3_lab_abund_np_mod) #based on this summary it doesn't 
#seem like sterilization acts on any of these groups yAY
plot(residuals(u3_lab_abund_np_mod))

#this is the null model with just a random intercept by species ID
u3_lab_abund_np_mod2 <- glmmTMB(reads ~ Ster + groups + (1|sample),
                                data = rarefied_u3_lab_grp_np,
                                family = "nbinom2",
                                REML = TRUE)

u3_lab_abund_np_mod3 <- glmmTMB(reads ~ groups + (1|sample),
                                data = rarefied_u3_lab_grp_np,
                                family = "nbinom2", 
                                REML = TRUE)

u3_lab_abund_np_mod4 <- glmmTMB(reads ~ Ster + (1|sample),
                                data = rarefied_u3_lab_grp_np,
                                family = "nbinom2",
                                REML = TRUE)

u3_lab_abund_np_mod_null <- glmmTMB(reads ~ 1 + (1|sample), 
                                    data = rarefied_u3_lab_grp_np,
                                    family = "nbinom2", 
                                    REML = TRUE)
#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(u3_fld_abund_np_mod, u3_fld_abund_np_mod2, u3_fld_abund_np_mod3, 
      u3_fld_abund_np_mod4, u3_fld_abund_np_mod_null)

u3_lab_abund_np_mod_null_ml <- glmmTMB(reads ~ 1 + (1|sample), 
                                       data = rarefied_u3_lab_grp_np,
                                       family = "nbinom2")
summary(u3_lab_abund_np_mod_null_ml)
simulationOutput <- simulateResiduals(fittedModel = u3_lab_abund_np_mod_null_ml) 
plot(simulationOutput, asFactor=TRUE) #these look um
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed

#DADA2 FIELD rarefy and PERMANOVA####
f_d2_nmds <- read.csv(here("data", "fld_d2uc_nmds.csv"))

#set rownames to sample variable
rownames(f_d2_nmds) <- f_d2_nmds$sample
#and then select just the numeric values in this dataframe
fld_d2_nmds <- f_d2_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#delete those zero columns since they represent "species" that are 
#not present in the field "communities"
fld_d2_nmds <- fld_d2_nmds[which(colSums(fld_d2_nmds)> 0)] #176 down to 145

#import the metadata for later
meta_d2_field <- read.csv(here("data", "fld_d2_nmds_meta.csv"))

raremin_d2_fld <- min(rowSums(fld_d2_nmds))
raremin_d2_fld #11229, which is substantially lower than some others...
rowSums(fld_d2_nmds)

# Rarefy, which means to standardize all samples by the lowest value (raremin). 
# change sample = to the rarefaction depth (# of reads per sample)
rarefied_d2_fld <- rrarefy(fld_d2_nmds, sample = raremin_d2_fld)
base::range(rowSums(rarefied_d2_fld)) # checks that it worked
rarefied_d2_fld <- as.data.frame(rarefied_d2_fld) # save it as a dataframe

colSums(rarefied_d2_fld) #145 total ZOTUS here
rarefied_d2_fld[which(colSums(rarefied_d2_fld) == 0)] #35 equal zero across all

#transpose dataframe so that we can make it long by sample and 
#attach metadata AND taxonomic assignments of ASVs to it
rarefied_d2_fld_t <- as.data.frame(t(rarefied_d2_fld))
rarefied_d2_fld_t$ASV <- row.names(rarefied_d2_fld_t) #make ASV column for join

#make long by sample, then join to taxonomic assignments.
rarefied_d2_fld_long <- rarefied_d2_fld_t %>%
  gather(sample, reads, HEV100:HEV89) %>%
  left_join(d2_taxa, by = "ASV")

#Groups field: predator, prey, unknown, no-hit, non-prey
rarefied_d2_fld_long$groups <- rarefied_d2_fld_long$taxonomy

rarefied_d2_fld_grp <- rarefied_d2_fld_long %>%
  group_by(sample, groups) %>%
  summarize(reads = sum(reads)) %>%
  left_join(meta_d2_field, by = "sample")

#set species ID to a factor
rarefied_d2_fld_grp$groups <- as.factor(rarefied_d2_fld_grp$groups)
rarefied_d2_fld_grp$Ster <- ifelse(rarefied_d2_fld_grp$Sterilized == "NS", 0, 1)

#remove predator
rarefied_d2_fld_grp_np <- rarefied_d2_fld_grp[which(rarefied_d2_fld_grp$groups != "predator"),]
rarefied_d2_fld_grp_np$groups <- droplevels(rarefied_d2_fld_grp_np$groups)

#this is our model, which is grouped by sample, saying:
#we are wondering what the effect of sterilization AND taxonomic group
#on reads, and let's say that sterilization may interact with groups
#such that the effect of sterilization may vary depending on group

d2_fld_abund_np_mod <- glmmTMB(reads ~ Ster*groups + (1|sample),
                               data = rarefied_d2_fld_grp_np,
                               family = "genpois",
                               REML = TRUE)

simulationOutput <- simulateResiduals(fittedModel = d2_fld_abund_np_mod) 
plot(simulationOutput, asFactor=TRUE) #these look um
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed
#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = u3_fld_abund_mod) 
#plot(simulationOutput, asFactor=TRUE) #these look um
#testZeroInflation(simulationOutput) #not zero inflated
#testDispersion(simulationOutput) #not overdispersed

#looks at the summary of this model
summary(d2_fld_abund_np_mod) #based on this summary it doesn't 
#seem like sterilization acts on any of these groups yAY
plot(residuals(d2_fld_abund_np_mod))

#this is the null model with just a random intercept by species ID
d2_fld_abund_np_mod2 <- glmmTMB(reads ~ Ster + groups + (1|sample),
                                data = rarefied_d2_fld_grp_np,
                                family = "genpois",
                                REML = TRUE)

d2_fld_abund_np_mod3 <- glmmTMB(reads ~ groups + (1|sample),
                                data = rarefied_d2_fld_grp_np,
                                family = "genpois", 
                                REML = TRUE)

d2_fld_abund_np_mod4 <- glmmTMB(reads ~ Ster + (1|sample),
                                data = rarefied_d2_fld_grp_np,
                                family = "nbinom1",
                                REML = TRUE)

d2_fld_abund_np_mod_null <- glmmTMB(reads ~ 1 + (1|sample), 
                                    data = rarefied_d2_fld_grp_np,
                                    family = "nbinom1", 
                                    REML = TRUE)

#ANOVA comparsion showing that the fulls is a better fit based on AIC and BIC
anova(d2_fld_abund_np_mod, d2_fld_abund_np_mod2, d2_fld_abund_np_mod3, 
      d2_fld_abund_np_mod4, d2_fld_abund_np_mod_null)


d2_fld_abund_np_mod_ml <- glmmTMB(reads ~ Ster*groups + (1|sample),
                               data = rarefied_d2_fld_grp_np,
                               family = "genpois")
#checking other model assumptions here:
simulationOutput <- simulateResiduals(fittedModel = d2_fld_abund_np_mod_ml) 
plot(simulationOutput, asFactor=TRUE) #these look um...
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed
#looks at the summary of this model
summary(d2_fld_abund_np_mod_ml)
plot(residuals(d2_fld_abund_np_mod_ml))

#DADA2 LAB RAREFY AND PERMANOVA ####
l_d2_nmds <- read.csv(here("data", "lab_d2uc_nmds.csv"))
#set rownames to sample variable
rownames(l_d2_nmds) <- l_d2_nmds$sample
#and then select just the numeric values in this dataframe
lab_d2_nmds <- l_d2_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#delete those zero columns since they represent "species" that are 
#not present in the field "communities"
lab_d2_nmds <- lab_d2_nmds[which(colSums(lab_d2_nmds)> 0)] #176 down to 75

#import the metadata for later
meta_d2_lab <- read.csv(here("data", "lab_d2uc_nmds_meta.csv"))

raremin_d2_lab <- min(rowSums(lab_d2_nmds))
rowSums(lab_d2_nmds)
raremin_d2_lab #50345, which is about half the max number

# Rarefy, which means to standardize all samples by the lowest value (raremin). 
# change sample = to the rarefaction depth (# of reads per sample)
rarefied_d2_lab <- rrarefy(lab_d2_nmds, sample = raremin_d2_lab)
base::range(rowSums(rarefied_d2_lab)) # checks that it worked
rarefied_d2_lab <- as.data.frame(rarefied_d2_lab) # save it as a dataframe

colSums(rarefied_d2_lab) #145 total ZOTUS here
rarefied_d2_lab[which(colSums(rarefied_d2_lab) == 0)] #1 equal zero across all

#transpose dataframe so that we can make it long by sample and 
#attach metadata AND taxonomic assignments of ASVs to it
rarefied_d2_lab_t <- as.data.frame(t(rarefied_d2_lab))
rarefied_d2_lab_t$ASV <- row.names(rarefied_d2_lab_t) #make ASV column for join

#make long by sample, then join to taxonomic assignments.
rarefied_d2_lab_long <- rarefied_d2_lab_t %>%
  gather(sample, reads, HEV07:HEV29) %>%
  left_join(d2_taxa, by = "ASV")

#Groups field: predator, prey, unknown, no-hit, non-prey
d2_lab_groups <- ifelse(rarefied_d2_lab_long$taxonomy == "prey" & rarefied_d2_lab_long$ID_ncbi == "Acrididae", "known prey",
                        ifelse(rarefied_d2_lab_long$taxonomy == "prey" & rarefied_d2_lab_long$ID_ncbi != "Acrididae", "other prey",
                               rarefied_d2_lab_long$taxonomy))
rarefied_d2_lab_long$groups <- d2_lab_groups

rarefied_d2_lab_grp <- rarefied_d2_lab_long %>%
  group_by(sample, groups) %>%
  summarize(reads = sum(reads)) %>%
  left_join(meta_d2_lab, by = "sample")
str(rarefied_d2_lab_grp)
#set species ID to a factor
rarefied_d2_lab_grp$groups <- as.factor(rarefied_d2_lab_grp$groups)
rarefied_d2_lab_grp$Ster <- ifelse(rarefied_d2_lab_grp$Sterilized == "NS", 0, 1)

#remove preator
rarefied_d2_lab_grp_np <- rarefied_d2_lab_grp[which(rarefied_d2_lab_grp$groups != "predator"),]
rarefied_d2_lab_grp_np$groups <- droplevels(rarefied_d2_lab_grp_np$groups)

#this is our model, which is grouped by sample, saying:
#we are wondering what the effect of sterilization AND taxonomic group
#on reads, and let's say that sterilization may interact with groups
#such that the effect of sterilization may vary depending on group

d2_lab_abund_np_mod <- glmmTMB(reads ~ Ster*groups + (1|sample),
                               data = rarefied_d2_lab_grp_np,
                               family = "genpois",
                               REML = TRUE)

simulationOutput <- simulateResiduals(fittedModel = d2_lab_abund_np_mod) 
plot(simulationOutput, asFactor=TRUE) #these look um
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed
#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = u3_fld_abund_mod) 
#plot(simulationOutput, asFactor=TRUE) #these look um
#testZeroInflation(simulationOutput) #not zero inflated
#testDispersion(simulationOutput) #not overdispersed

#looks at the summary of this model
summary(d2_lab_abund_np_mod) #based on this summary it doesn't 
#seem like sterilization acts on any of these groups yAY
plot(residuals(d2_lab_abund_np_mod))

#this is the null model with just a random intercept by species ID
d2_lab_abund_np_mod2 <- glmmTMB(reads ~ Ster + groups + (1|sample),
                                data = rarefied_d2_lab_grp_np,
                                family = "genpois",
                                REML = TRUE)

d2_lab_abund_np_mod3 <- glmmTMB(reads ~ groups + (1|sample),
                                data = rarefied_d2_lab_grp_np,
                                family = "genpois", 
                                REML = TRUE)

d2_lab_abund_np_mod4 <- glmmTMB(reads ~ Ster + (1|sample),
                                data = rarefied_d2_lab_grp_np,
                                family = "genpois",
                                REML = TRUE)

d2_lab_abund_np_mod_null <- glmmTMB(reads ~ 1 + (1|sample), 
                                    data = rarefied_d2_lab_grp_np,
                                    family = "genpois", 
                                    REML = TRUE)

#ANOVA comparsion showing that the fulls is a better fit based on AIC and BIC
anova(d2_lab_abund_np_mod, d2_lab_abund_np_mod2, d2_lab_abund_np_mod3, 
      d2_lab_abund_np_mod4, d2_lab_abund_np_mod_null)


d2_lab_abund_np_mod_ml <- glmmTMB(reads ~ Ster*groups + (1|sample),
                                  data = rarefied_d2_lab_grp_np,
                                  family = "genpois")
#checking other model assumptions here:
simulationOutput <- simulateResiduals(fittedModel = d2_lab_abund_np_mod_ml) 
plot(simulationOutput, asFactor=TRUE) #these look um...
testZeroInflation(simulationOutput) #not zero inflated
testDispersion(simulationOutput) #not overdispersed
#looks at the summary of this model
summary(d2_lab_abund_np_mod_ml)
plot(residuals(d2_lab_abund_np_mod_ml))

#Adonis PERMANOVA abundance ####
u3_fld_adonis <- rarefied_u3_fld_grp %>%
  dplyr::select(sample, groups, reads) %>%
  spread(groups, reads) %>%
  ungroup()

rownames(u3_fld_adonis) <- u3_fld_adonis$sample
u3_fld_adonis <- u3_fld_adonis %>%
  dplyr::select(-sample)

adonis(u3_fld_adonis ~ Sterilized, data = u3_meta_field, dist = "bray")

u3_lab_adonis <- rarefied_u3_lab_grp %>%
  dplyr::select(sample, groups, reads) %>%
  spread(groups, reads) %>%
  ungroup()

rownames(u3_lab_adonis) <- u3_lab_adonis$sample
u3_lab_adonis <- u3_lab_adonis %>%
  dplyr::select(-sample)

adonis(u3_lab_adonis ~ Sterilized, data = u3_meta_lab, dist = "bray")


d2_fld_adonis <- rarefied_d2_fld_grp %>%
  dplyr::select(sample, groups, reads) %>%
  spread(groups, reads) %>%
  ungroup()

rownames(d2_fld_adonis) <- d2_fld_adonis$sample
d2_fld_adonis <- d2_fld_adonis %>%
  dplyr::select(-sample)

adonis(d2_fld_adonis ~ Sterilized, data = d2_meta_field, dist = "bray")

d2_lab_adonis <- rarefied_d2_lab_grp %>%
  dplyr::select(sample, groups, reads) %>%
  spread(groups, reads) %>%
  ungroup()

rownames(d2_lab_adonis) <- d2_lab_adonis$sample
d2_lab_adonis <- d2_lab_adonis %>%
  dplyr::select(-sample)

adonis(d2_lab_adonis ~ Sterilized, data = d2_meta_lab, dist = "bray")

#Bray-Curtis Abundance based similarity ####
#field u3
bray_u3_field <- metaMDS(u3_fld_adonis, distance = "bray", k=2)
#evaluate the fit of this NMDS
stressplot(bray_u3_field)
bray_u3_field$stress #0.130425

#create a dataframe for plotting
bray_u3field_df <- data.frame(MDS1=bray_u3_field$points[,1], MDS2=bray_u3_field$points[,2])
#combining with metadata for plotting
bray_u3field_meta <- cbind(u3_meta_field, bray_u3field_df)

#plots NMDS with elipses
bray_u3_fld_grph <- ggplot(bray_u3field_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    theme_bw() +
    labs(title = "UNOISE3 Field Bray-Curtis dissimilarity") +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#lab u3
bray_u3_lab <- metaMDS(u3_lab_adonis, distance = "bray", k=2)
#evaluate the fit of this NMDS
stressplot(bray_u3_lab)
bray_u3_lab$stress #0.1361582

#create a dataframe for plotting
bray_u3lab_df <- data.frame(MDS1=bray_u3_lab$points[,1], MDS2=bray_u3_lab$points[,2])
#combining with metadata for plotting
bray_u3lab_meta <- cbind(u3_meta_lab, bray_u3lab_df)

#plots NMDS with elipses
bray_u3_lab_grph <- ggplot(bray_u3lab_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    theme_bw() +
    labs(title = "UNOISE3 Lab Bray-Curtis dissimilarity") +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#field d2
bray_d2_field <- metaMDS(d2_fld_adonis, distance = "bray", k=2)
#evaluate the fit of this NMDS
stressplot(bray_d2_field)
bray_d2_field$stress #0.1042269

#create a dataframe for plotting
bray_d2field_df <- data.frame(MDS1=bray_d2_field$points[,1], MDS2=bray_d2_field$points[,2])
#combining with metadata for plotting
bray_d2field_meta <- cbind(d2_meta_field, bray_d2field_df)

#plots NMDS with elipses
bray_d2_fld_grph <- ggplot(bray_d2field_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    theme_bw() +
    labs(title = "DADA2 Field Bray-Curtis dissimilarity") +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#lab d2
bray_d2_lab <- metaMDS(d2_lab_adonis, distance = "bray", k=2)
#evaluate the fit of this NMDS
stressplot(bray_d2_lab)
bray_d2_lab$stress #0.1276943

#create a dataframe for plotting
bray_d2lab_df <- data.frame(MDS1=bray_d2_lab$points[,1], MDS2=bray_d2_lab$points[,2])
#combining with metadata for plotting
bray_d2lab_meta <- cbind(d2_meta_lab, bray_d2lab_df)

#plots NMDS with elipses
bray_d2_lab_grph <- ggplot(bray_d2lab_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    theme_bw() +
    labs(title = "DADA2 Lab Bray-Curtis dissimilarity") +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

plot_grid(bray_u3_fld_grph, bray_d2_fld_grph,
          bray_u3_lab_grph, bray_d2_lab_grph, nrow = 2)
