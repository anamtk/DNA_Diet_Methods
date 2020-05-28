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
# Field presence-absence composition analysis (PERMANOVA) ####
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
                       data = comp,
                       family = "binomial")

comp_null <- glmmTMB(presence ~ 1 + (1|unique_ID), 
                           data = comp,
                           family = "binomial")

AICc(comp_mod, comp_null)

summary(comp_null)
summary(comp_mod) #sterilization term not significant

#model diagnostics look good
binned_residuals(comp_null)
simulationOutput <- simulateResiduals(fittedModel = comp_null)
fit <- plot(simulationOutput, asFactor=TRUE)


###########################
# Field Abundance-based composition analysis (PERMANOVA)####
############################

#We set up the GLMM in a really similar way, only this time we use read
#abundance as our dependent variable, and we use a poisson distribution,
#since these are now count data rather than binary data

#This GLMM is run by saying, how does the fixed effect of sterilization impact
#abundance, with a random effects structure with both a random interept term
#for unique_ID (let each species have a different intercept) and a random slopes 
#term for surface sterilization treatment (let each species' relationship with
#with surface sterilization differ, ie let some species increase with surface 
#sterilization, and others decrease)
#again, we set REML=FALSE for AIC comparison FIRST, and then refit with REML
#for model summaries and diagnostics

bray_mod <- glmmTMB(reads ~ Sterilized + (1+Sterilized|unique_ID), 
                    data = comp,
                    family = "genpois",
                    REML = FALSE)

bray_null <- glmmTMB(reads ~ 1 + (1|unique_ID), 
                     data = comp,
                     family = "genpois",
                     REML = FALSE)

AICc(bray_mod, bray_null)

#refit with REML
bray_mod <- glmmTMB(reads ~ Sterilized + (1+Sterilized|unique_ID), 
                    data = comp,
                    family = "genpois")

bray_null <- glmmTMB(reads ~ 1 + (1|unique_ID), 
                     data = comp,
                     family = "genpois")

summary(bray_null)
summary(bray_mod) #sterilization term not significant

#model diagnostics look good
plot(residuals(bray_null))
simulationOutput <- simulateResiduals(fittedModel = bray_null)
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput) 

###########################
# Heat map visualization of both presence and abundance ####
#(At order level for clarity)
############################
taxonomies <- read.csv(here("data", "outputs", "taxonomic_assignments", "taxonomic_assignments.csv"))

u3_comm_ord <- u3_comm_rare_f %>%
  left_join(u3_id_all, by = "ASV") %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(-taxonomy)

#UNOISE3
#make the ID columns characters for the ifelse statements below
u3_comm_ord$Order_ncbi <- as.character(u3_comm_ord$Order_ncbi)
u3_comm_ord$ID_ncbi <- as.character(u3_comm_ord$ID_ncbi)

#ifelse statements prioritizing NCBI since there are more.
Orders_u3 <- ifelse(u3_comm_ord$Order_ncbi == "", u3_comm_ord$ID_ncbi,
                    u3_comm_ord$Order_ncbi)
#Orders_u3

#make this a dataframe so we can attach it back to the community
#matrix for analyses at the sample level and later Jaccard dissimiliarity
Orders_u3 <- as.data.frame(Orders_u3)
Orders_u3$ASV <- u3_comm_ord$ASV
Orders_u3 <- rename(Orders_u3, "unique_order" = "Orders_u3")
unique_Orders_u3 <- as.data.frame(unique(Orders_u3[c("unique_order")]))
unique_Orders_u3$sp_number <- seq.int(nrow(unique_Orders_u3))

ords_comm_u3 <- Orders_u3 %>%
  left_join(unique_Orders_u3, by = "unique_order") %>%
  left_join(u3_comm_ord, by = "ASV") %>%
  gather(sample, reads, HEV65:HEV100) %>%
  left_join(metadata, by = "sample") %>%
  ungroup() %>%
  group_by(sample, unique_order, Sterilized) %>%
  summarize(reads = sum(reads)) %>%
  mutate(pipeline = "u3", presence = ifelse(reads > 0, 1, 0)) %>%
  group_by(unique_order) %>%
  filter(sum(presence) > 0) #from 518 to 407, removed 3 orders?


orders_comm <- ords_comm_u3 %>%
  ungroup() %>%
  mutate(unique_order = as.factor(unique_order), surf_ster=ifelse(Sterilized == "NS", 0, 1)) 

levels(orders_comm$unique_order)

heat_map <- orders_comm %>%
  group_by(unique_order, Sterilized) %>%
  summarise(reads = mean(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0))

heat_map_nz <- heat_map %>%
  filter(reads > 0)
quantile(heat_map_nz$reads)

heat_map <- heat_map %>%
  ungroup() %>%
  dplyr::select(Sterilized, unique_order, reads) %>%
  pivot_wider(names_from = Sterilized, values_from = reads) %>%
  arrange((NS + SS), (NS)) %>%
  mutate(unique_order=factor(unique_order, levels=unique_order)) %>%
  gather(Sterilized, reads, NS:SS) 


heat_map$quantiles <- ifelse(heat_map$reads == 0, 0,
                             ifelse(heat_map$reads > 0 & heat_map$reads <= 0.1052632, 1, 
                                    ifelse(heat_map$reads > 0.1052632 & heat_map$reads <= 0.7456140, 2,
                                           ifelse(heat_map$reads > 0.7456140 & heat_map$reads <= 2.8128655, 3,
                                                  ifelse(heat_map$reads > 2.8128655 & heat_map$reads <= 21.5614035, 4, 5)))))

heat_map$quantiles <- as.factor(heat_map$quantiles)  

pal <- c(
  '0' = "#FFFFFF",
  '1' = "#F27D72", 
  '2' = "#D26F67", 
  '3' = "#B2615C",
  '4' = "#925451",
  '5' = "#734646"
)

pal2 <- c(
  '0' = "#FFFFFF",
  '1' = "#F29979", 
  '2' = "#D2846C", 
  '3' = "#B26F5F",
  '4' = "#925B53",
  '5' = "#734646"
)

heat_map_graph <- ggplot(heat_map, aes(x = Sterilized, y = unique_order, fill=quantiles, height = 0.95, width = 0.95)) +
  geom_tile() + 
  coord_equal() +
  labs(x = "Surface sterilization treatment", y = "Diet group") +
  scale_x_discrete(labels=c("NS" = "Not Sterilized", "SS" = "Surface Sterilized")) +
  scale_fill_manual(name = "Mean read abundance\n(divided by quantiles)",
                    values = pal2,
                    limits = names(pal2),
                    labels = c("0", "< 0.1", "< 0.7", "< 2.8", "< 21.6", "< 219.5")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###########################
# Effect Sizes Presence graph ####
############################

#This is a graph of the effect sizes of presence-absence between surface sterilized 
#and non-surface sterilized predator individuals

#Effect sizes are built off of means and either SD or SE
effect <- comp %>%
  group_by(unique_ID, Sterilized) %>%
  summarize(mean = mean(presence), sd = sd(presence)) %>%
  ungroup() %>%
  group_by(unique_ID) %>%
  pivot_wider(names_from = c(Sterilized),
              values_from = c(mean, sd)) %>% #this pivots so there is an average for
  #SS and NS groups
  mutate(se_NS = sd_NS/sqrt(19), se_SS = sd_SS/sqrt(18)) #this computs the SE for groups

#this is the overall presence of each species, which we will use to sort the 
#graph visualization
pres_sort <- comp %>% 
  group_by(unique_ID) %>%
  summarise(overall = mean(presence)) #gets the overall presence of that diet item

#this computes the effect sies based on means, standared errors, and sample sizes
es_ID <- esc_mean_se(grp1m = effect$mean_NS, grp1se = effect$se_NS, grp1n = 19,
                        grp2m = effect$mean_SS, grp2se = effect$se_SS, grp2n = 18, es.type = "g")

#extract data of interest from teh es_ID object
IDs <- as.character(effect$unique_ID)
Hedges_g <- es_ID$es
Lower_CI <- es_ID$ci.lo
Upper_CI <- es_ID$ci.hi

#make into a DF
effects <- as.data.frame(cbind(IDs, Hedges_g, Lower_CI, Upper_CI))

#manipulate data types
effects$Hedges_g <- as.numeric(as.character(effects$Hedges_g))
effects$Lower_CI <- as.numeric(as.character(effects$Lower_CI))
effects$Upper_CI <- as.numeric(as.character(effects$Upper_CI))
effects$IDs <- as.factor(effects$IDs)

#join this with the overall presence DF so we can order by overall
#presence for graph visualization
effects1 <- effects %>%
  left_join(pres_sort, by = c("IDs" = "unique_ID")) %>%
  arrange(overall) %>%
  mutate(IDs=factor(IDs, levels=IDs)) 

#this is the visualization, with the IDs at the top being most present on average
#and decreasing abundance as you go to the bottom of the graph
ggplot(effects1, aes(x = IDs, y = Hedges_g)) +
  geom_point() +theme_bw() +geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "UNOISE3 effect size of average presence") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + coord_flip()

###########################
# Effect Sizes Abundance graph ####
############################

#Same as above, just with abundance data instead
#This is a graph of the effect sizes of abundance between surface sterilized 
#and non-surface sterilized predator individuals

#Effect sizes are built off of means and either SD or SE
effect_abund <- comp %>%
  group_by(unique_ID, Sterilized) %>%
  summarize(mean = mean(reads), sd = sd(reads)) %>%
  ungroup() %>%
  group_by(unique_ID) %>%
  pivot_wider(names_from = c(Sterilized),
              values_from = c(mean, sd)) %>% #this pivots so there is an average for
  #SS and NS groups
  mutate(se_NS = sd_NS/sqrt(19), se_SS = sd_SS/sqrt(18)) #this computs the SE for groups

#this is the overall presence of each species, which we will use to sort the 
#graph visualization
abund_sort <- comp %>% 
  group_by(unique_ID) %>%
  summarise(overall = mean(reads)) #gets the overall average reads of that diet item

#this computes the effect sies based on means, standared errors, and sample sizes
es_ID_abund <- esc_mean_se(grp1m = effect_abund$mean_NS, grp1se = effect_abund$se_NS, grp1n = 19,
                     grp2m = effect_abund$mean_SS, grp2se = effect_abund$se_SS, grp2n = 18, es.type = "g")

#extract data of interest from teh es_ID object
IDs <- as.character(effect_abund$unique_ID)
Hedges_g <- es_ID_abund$es
Lower_CI <- es_ID_abund$ci.lo
Upper_CI <- es_ID_abund$ci.hi

#make into a DF
effects_abund <- as.data.frame(cbind(IDs, Hedges_g, Lower_CI, Upper_CI))

#manipulate data types
effects_abund$Hedges_g <- as.numeric(as.character(effects_abund$Hedges_g))
effects_abund$Lower_CI <- as.numeric(as.character(effects_abund$Lower_CI))
effects_abund$Upper_CI <- as.numeric(as.character(effects_abund$Upper_CI))
effects_abund$IDs <- as.factor(effects_abund$IDs)

#join this with the overall presence DF so we can order by overall
#presence for graph visualization
effects_abund1 <- effects_abund %>%
  left_join(abund_sort, by = c("IDs" = "unique_ID")) %>%
  arrange(overall) %>%
  mutate(IDs=factor(IDs, levels=IDs)) 

#this is the visualization, with the IDs at the top being most abundant
#and decreasing abundance as you go to the bottom of the graph
ggplot(effects_abund1, aes(x = IDs, y = Hedges_g)) +
  geom_point() +theme_bw() +geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "UNOISE3 effect size of average abundance") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + coord_flip()

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

