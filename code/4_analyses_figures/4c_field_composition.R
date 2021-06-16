###########################
# FIeld Richness Analyses ####
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

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "MuMIn", 
                  "DHARMa", "performance",
                  "vegan", "remotes",
                  "eulerr", "dummies",
                  "esc")

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
                       "2c_rarefied_taxonomic_sort", 
                       "field_prey_rare.csv"))

ASV <- read.csv(here("data", 
                     "outputs", 
                     "2c_rarefied_taxonomic_sort", 
                     "field_prey_ASVs_rare.csv"))

###########################
# Field presence-absence composition analysis (PERMANOVA) ####
############################

#This PERMANOVA is run like a GLMM, but we need to first verify that all samples
#have greater than zero species presence in them, or the model will not work:

comp <- field %>%
  group_by(sample) %>%
  filter(sum(reads) >0) %>% #remove samples with zero prey
  ungroup() %>%
  mutate(Family_ncbi = as.factor(Family_ncbi), presence = ifelse(reads > 0, 1, 0))
#make Family_ncbi a factor and create a presence column

#This GLMM is run by saying, how does the fixed effect of sterilization impact
#presence, with a random effects structure with both a random interept term
#for Family_ncbi (let each family have a different intercept) and a random slopes 
#term for surface sterilization treatment (let each family's relationship with
#with surface sterilization differ, ie let some families increase with surface 
#sterilization, and others decrease)

comp_mod <- glmmTMB(presence ~ Sterilized + (1+Sterilized|Family_ncbi), 
                       data = comp,
                       family = "binomial")

comp_null <- glmmTMB(presence ~ 1 + (1|Family_ncbi), 
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
# Heat map visualization of both presence and abundance ####
#(At family level)
############################
#make order and ID values characters for ifelse sorting to order level below
heat <- field

heat$Family_ncbi <- as.character(heat$Family_ncbi)

#ifelse statements creating an overall "Order" column for visualization
heat$Family <- heat$Family_ncbi

#now we can summarise heat by order and sterilization treatment for
#the heatmap graph
heat_family <- heat %>%
  group_by(Family, Sterilized) %>% #group by order and sterilization treatment
  summarise(reads = mean(reads)) %>% #then find the mean abundance for each order
  ungroup() %>% #ungroup so we can make order a factor
  mutate(Family = as.factor(Family)) %>% #make order a factor
  pivot_wider(names_from = Sterilized, values_from = reads) %>% #make wide for two columns
  #for sterilization treatment so we can arrange in descending order of abundance
  arrange((NS + SS), (NS)) %>% #arrange in descending order of abundance
  mutate(Family = factor(Family, levels = Family)) %>% #reset the levels of order on this
  #new abundance-based ordering of the Order factor
  gather(Sterilized, reads, NS:SS) #gather DF back up for visualization

#Because read abundance has such a broad range, we want to set a quantile-based
#variable for these abundances for the color ramp in the graph, but also want
#to discount the effects of zero reads in this, so we will create a DF of non-zero
#values for reads, and then find quantiles of this to inform our quantiles for the figure
heat_family_nz <- heat_family %>%
  filter(reads > 0)
quantile(heat_family_nz$reads)
#0.05263158   0.37353801   1.32602339   6.77631579 220.05263158 

#now we can set a quantile variable in our DF
heat_family$quantiles <- ifelse(heat_family$reads == 0, 0,
                               ifelse(heat_family$reads > 0 & heat_family$reads <= 0.05263158, 1, 
                                      ifelse(heat_family$reads > 0.05263158 & heat_family$reads <= 0.37353801, 2,
                                             ifelse(heat_family$reads > 0.37353801 & heat_family$reads <= 1.32602339, 3,
                                                    ifelse(heat_family$reads > 1.32602339 & heat_family$reads <= 6.77631579, 4, 5)))))

#making it a factor for visualization
heat_family$quantiles <- as.factor(heat_family$quantiles)  

#these are two color palettes that could be used in this figure
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

(heat_map_family <- ggplot(heat_family, aes(x = Sterilized, y = Family, fill=quantiles, height = 0.95, width = 0.95)) +
    geom_tile() + 
    coord_equal() +
    labs(x = "Surface sterilization treatment", y = "Diet family") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
    scale_fill_manual(name = "Mean read abundance\n(divided by quantiles)",
                      values = pal2,
                      limits = names(pal2),
                      labels = c("0", "< 0.05", "< 0.37", "< 1.33", "< 6.78", "> 6.78")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

##0.05263158   0.37353801   1.32602339   6.77631579 220.05263158 
factor_levels <- heat_family %>%
  pull(Family)

###########################
# Field presence-absence composition analysis CCA ####
############################

# Matrix for CCA -----------------------------------------------
zeros <- ASV %>%
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>%
  filter(reads == 0)

ASV <- ASV %>%
  anti_join(zeros, by = "sample")

#matrix with diet as columns, samples as rows
mat <- ASV %>%
  dplyr::select(sample, reads, ASV) %>%
  pivot_wider(names_from = ASV,
              values_from = reads,
              values_fill = 0) %>% #fill missing values with 0's
  column_to_rownames(var = "sample")

mat_pres <- ASV %>%
  dplyr::select(sample, reads, ASV) %>%
  mutate(presence = ifelse(reads>0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  pivot_wider(names_from = ASV,
              values_from = presence,
              values_fill = 0) %>% #fill missing values with 0's
  column_to_rownames(var = "sample")

#metadata
meta <- ASV %>%
  dplyr::select(sample, Sterilized) %>%
  distinct(sample, Sterilized) %>%
  column_to_rownames(var = "sample")

#check that rownames are the same:
all.equal(rownames(mat), rownames(meta))

#check that rownames are the same:
all.equal(rownames(mat_pres), rownames(meta))


# CCA -----------------------------------------------------------
#working off this tutorial
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
# alternatively we could use a shorthand approach
cca <-  cca(mat_pres ~ . , data=meta)

#view the CCA loadings
cca

summary(cca)
#how much variation explained in this RDA
RsquareAdj(cca)
#ANOVA of whole model
anova(cca, permutations=10000)
#ANOVA of model terms
anova(cca, by='margin', permutations=10000)


###########################
# (OPTIONAL: adonis instead of GLMM for field presence) ####
############################
#Now we can run adonis, which asks whether there are differences in community
#composition of diet species between the two treatment levels of sterilization
#dist = jaccard with binary = TRUE because this is presence-absence data
adonis(mat_pres ~ Sterilized, data = meta, dist = "jaccard", binary = TRUE)
#there is no significnat difference between treatment groups 


# Heat Map of ASVs ----------------------------------------------------------------
heat <- ASV

taxa <- read.csv(here("data",
                      "outputs",
                      "taxonomic_assignments",
                      "taxonomic_assignments.csv"))

taxa <- taxa %>%
  dplyr::select(ASV, Family_ncbi) %>%
  filter(Family_ncbi != "")

heat <- heat %>% 
  left_join(taxa, by = "ASV") %>%
  unite("Family_ID", c(ASV, Family_ncbi), sep = " (", remove = F) %>%
  mutate(Family_ID = paste0(Family_ID, ")"))

#now we can summarise heat by order and sterilization treatment for
#the heatmap graph
heat_ASV <- heat %>%
  group_by(ASV, Family_ncbi, Family_ID, Sterilized) %>% #group by order and sterilization treatment
  summarise(reads = mean(reads)) %>% #then find the mean abundance for each order
  ungroup() %>% #ungroup so we can make order a factor
  mutate(ASV = as.factor(ASV)) %>% #make order a factor
  pivot_wider(names_from = Sterilized, values_from = reads) %>% #make wide for two columns
  #for sterilization treatment so we can arrange in descending order of abundance
  arrange((NS + SS), (NS)) %>% #arrange in descending order of abundance
  mutate(ASV = factor(ASV, levels = ASV)) %>% #reset the levels of order on this
  #new abundance-based ordering of the Order factor
  gather(Sterilized, reads, NS:SS) #gather DF back up for visualization

#Because read abundance has such a broad range, we want to set a quantile-based
#variable for these abundances for the color ramp in the graph, but also want
#to discount the effects of zero reads in this, so we will create a DF of non-zero
#values for reads, and then find quantiles of this to inform our quantiles for the figure
heat_ASV_nz <- heat_ASV %>%
  filter(reads > 0)
quantile(heat_ASV_nz$reads)
#0.05882353   0.18750000   0.64705882   6.02941176 245.94117647

#now we can set a quantile variable in our DF
heat_ASV$quantiles <- ifelse(heat_ASV$reads == 0, 0,
                                ifelse(heat_ASV$reads > 0 & heat_ASV$reads <= 0.05882353, 1, 
                                       ifelse(heat_ASV$reads > 0.05882353 & heat_ASV$reads <= 0.18750000, 2,
                                              ifelse(heat_ASV$reads > 0.18750000 & heat_ASV$reads <= 0.64705882, 3,
                                                     ifelse(heat_ASV$reads > 0.64705882 & heat_ASV$reads <= 6.02941176, 4, 5)))))

#making it a factor for visualization
heat_ASV$quantiles <- as.factor(heat_ASV$quantiles)  

#these are two color palettes that could be used in this figure
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

(heat_map_ASV <- ggplot(heat_ASV, aes(x = Sterilized, y = ASV, fill=quantiles, height = 0.95, width = 0.95)) +
    geom_tile() + 
    coord_equal() +
    labs(x = "Surface sterilization treatment", y = "Diet ASV") +
    scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
    scale_fill_manual(name = "Mean read abundance\n(divided by quantiles)",
                      values = pal2,
                      limits = names(pal2),
                      labels = c("0", "< 0.05", "< 0.19", "< 0.65", "< 6.03", "> 6.03")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

#0.05882353   0.18750000   0.64705882   6.02941176 245.94117647

###########################
# Per-sample heat map as an alternative to effects map? ####
############################

heat_nz <- heat %>%
  filter(reads > 0)
quantile(heat_nz$reads)
#1    1    6   25 4079  

#now we can set a quantile variable in our DF
field %>%
  dplyr::select(sample, Sterilized) %>%
  unique()

heat$sample <- factor(heat$sample, levels = c("HEV100", "HEV106", "HEV107", "HEV108",
                                              "HEV109", "HEV110", "HEV111", "HEV79",
                                              "HEV81", "HEV82", "HEV83", "HEV87",
                                              "HEV88", "HEV89", "HEV95", "HEV96", 
                                              "HEV97", "HEV98", "HEV99", "HEV101",
                                              "HEV102", "HEV103", "HEV104", "HEV105",
                                              "HEV65", "HEV66", "HEV67", "HEV68", "HEV70",
                                              "HEV71", "HEV74", "HEV76", "HEV90", 
                                              "HEV91", "HEV92", "HEV93", "HEV94"))
heat$quantiles <- ifelse(heat$reads == 0, 0,
                                ifelse(heat$reads > 0 & heat$reads <= 1, 1, 
                                              ifelse(heat$reads > 1 & heat$reads <= 6, 3,
                                                     ifelse(heat$reads > 6 & heat$reads <= 25, 4, 5))))

#making it a factor for visualization
heat$quantiles <- as.factor(heat$quantiles)  

#these are two color palettes that could be used in this figure
pal3 <- c(
  '0' = "#FFFFFF",
  '1' = "#F27D72", 
  '3' = "#D26F67", 
  '4' = "#B2615C",
  '5' = "#925451"
)

heat1 <- heat %>%
  filter(reads > 0)
(heat_map <- ggplot(heat, aes(x = sample, y = Family, fill=quantiles)) +
    geom_tile(color = "black") + 
    coord_equal() +
    labs(x = "Sample", y = "Diet family") +
    scale_fill_manual(name = "Mean read abundance\n(divided by quantiles)",
                      values = pal3,
                      limits = names(pal3),
                      labels = c("0", "< 1", "< 6", "< 25", "> 25")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

##0.05263158   0.37353801   1.32602339   6.77631579 220.05263158 

###########################
# Diet species names for stats ####
############################
levels(comp$Family_ncbi)
unique(comp$Order_ncbi)
###########################
# (SUPPLEMENT: Field Abundance-based composition analysis (PERMANOVA))####
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

bray_mod <- glmmTMB(reads ~ Sterilized + (1+Sterilized|Family_ncbi), 
                    data = comp,
                    family = "genpois",
                    REML = FALSE)

bray_null <- glmmTMB(reads ~ 1 + (1|Family_ncbi), 
                     data = comp,
                     family = "genpois",
                     REML = FALSE)

AICc(bray_mod, bray_null)

#refit with REML
bray_mod <- glmmTMB(reads ~ Sterilized + (1+Sterilized|Family_ncbi), 
                    data = comp,
                    family = "genpois")

bray_null <- glmmTMB(reads ~ 1 + (1|Family_ncbi), 
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
# (SUPPLEMENT: Effect Sizes Presence graph) ####
############################

#This is a graph of the effect sizes of presence-absence between surface sterilized 
#and non-surface sterilized predator individuals

#Effect sizes are built off of means and either SD or SE
effect <- comp %>%
  group_by(Family_ncbi, Sterilized) %>%
  summarize(mean = mean(presence), sd = sd(presence)) %>%
  ungroup() %>%
  group_by(Family_ncbi) %>%
  pivot_wider(names_from = c(Sterilized),
              values_from = c(mean, sd)) %>% #this pivots so there is an average for
  #SS and NS groups
  mutate(se_NS = sd_NS/sqrt(19), se_SS = sd_SS/sqrt(18)) #this computs the SE for groups

#this is the overall presence of each species, which we will use to sort the 
#graph visualization
pres_sort <- comp %>% 
  group_by(Family_ncbi) %>%
  summarise(overall = mean(presence)) #gets the overall presence of that diet item

#this computes the effect sies based on means, standared errors, and sample sizes
es_ID <- esc_mean_se(grp1m = effect$mean_NS, grp1se = effect$se_NS, grp1n = 19,
                        grp2m = effect$mean_SS, grp2se = effect$se_SS, grp2n = 18, es.type = "g")

#extract data of interest from teh es_ID object
IDs <- as.character(effect$Family_ncbi)
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
  left_join(pres_sort, by = c("IDs" = "Family_ncbi")) %>%
  arrange(overall) %>%
  mutate(IDs=factor(IDs, levels=IDs)) 

#this is the visualization, with the IDs at the top being most present on average
#and decreasing abundance as you go to the bottom of the graph
pres_effect_f <- ggplot(effects1, aes(x = IDs, y = Hedges_g)) +
  geom_point() +theme_bw() +geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "UNOISE3 effect size of average presence") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + coord_flip()

###########################
# (SUPPLEMENT) Effect Sizes Abundance graph ####
############################

#Same as above, just with abundance data instead
#This is a graph of the effect sizes of abundance between surface sterilized 
#and non-surface sterilized predator individuals

#Effect sizes are built off of means and either SD or SE
effect_abund <- comp %>%
  group_by(Family_ncbi, Sterilized) %>%
  summarize(mean = mean(reads), sd = sd(reads)) %>%
  ungroup() %>%
  group_by(Family_ncbi) %>%
  pivot_wider(names_from = c(Sterilized),
              values_from = c(mean, sd)) %>% #this pivots so there is an average for
  #SS and NS groups
  mutate(se_NS = sd_NS/sqrt(19), se_SS = sd_SS/sqrt(18)) #this computs the SE for groups

#this is the overall presence of each species, which we will use to sort the 
#graph visualization
abund_sort <- comp %>% 
  group_by(Family_ncbi) %>%
  summarise(overall = mean(reads)) #gets the overall average reads of that diet item

#this computes the effect sies based on means, standared errors, and sample sizes
es_ID_abund <- esc_mean_se(grp1m = effect_abund$mean_NS, grp1se = effect_abund$se_NS, grp1n = 19,
                     grp2m = effect_abund$mean_SS, grp2se = effect_abund$se_SS, grp2n = 18, es.type = "g")

#extract data of interest from teh es_ID object
IDs <- as.character(effect_abund$Family_ncbi)
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
  left_join(abund_sort, by = c("IDs" = "Family_ncbi")) %>%
  arrange(overall) %>%
  mutate(IDs=factor(IDs, levels=IDs)) 

#this is the visualization, with the IDs at the top being most abundant
#and decreasing abundance as you go to the bottom of the graph
abund_effect_f <- ggplot(effects_abund1, aes(x = IDs, y = Hedges_g)) +
  geom_point() +theme_bw() +geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "UNOISE3 effect size of average abundance") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + coord_flip()


###########################
# SUPPLEMENT: Plot of effects for export####
############################

plot_grid(pres_effect_f, abund_effect_f)

###########################
# (OPTIONAL: adonis instead of GLMM for field presence) ####
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
  dplyr::select(presence, sample, Family_ncbi) %>% #select only the variables for matrix
  pivot_wider(names_from = Family_ncbi, values_from = presence) %>% 
  #make column names based on unique ID, values from presence in cells
  column_to_rownames(var = "sample") 
#set the column names to the sample, so that only numeric values are in the matrix

#adonis requires a matrix as input, so we will convert df to that
comp1 <- data.matrix(comp1, rownames.force = TRUE)

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

