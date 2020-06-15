###########################
# (SUPPLEMENT: Lab All Prey Richness Analyses) ####
# Ana Miller-ter Kuile
# May 28, 2020 
############################

#This script looks at composition (both presence-absence and abundance) of 
#ALL potential prey in the mesocosm-fed spiders. This is supplementary
#to other analyses, but just to demonstrate that the only real change in 
#this set of consumers is in presence of the offered prey item (suggesting
#offered prey in the mesocosm presented some level of contamination)
#Since these are just supplemental, and because there are fewer data, which seems
#to pose a challenge for the glmmTMB funciton, rather than running as GLMMs, I'm
#going to use the adonis fucntion from vegan instead. Makes it shorter too

###########################
# Required packages ####
############################

library(here) #tidy data
library(tidyverse) #tidy data
library(ggplot2) #visualize data
library(vegan) #adonis
library(esc) #effect sizes

###########################
# Load data ####
############################

lab <- read.csv(here("data", "outputs", "rarefied_taxonomic_sort", "lab_all_prey_rare.csv"))

###########################
# Adonis presence-absence ####
############################

#Adonis requires a matrix-type object with species as column names and 
#samples/sites as rows, which the following pipe does for the lab df.
#sites with zero across break adonis, so we need to remove them first
comp_lab <- lab %>%
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
comp_lab <- data.matrix(comp_lab, rownames.force = TRUE)

#Now that those are removed, we need to remove them from the eventual metadataframe
#we are creating next, so I will create a dataframe of those sample names to subset
#that dataframe
samples <- as.data.frame(rownames(comp_lab))

#Adonis also requires a meta-data DF for the analysis based on the factors of interest
#for the analysis. In this analysis, we are interested in comparing surface sterilized
#to unsterilized, so we can subset these from our DF as well

meta_lab <- lab %>%
  dplyr::select(sample, Sterilized) %>% #select metadata on sterilization and sample
  distinct(sample, Sterilized) %>% #reduce long DF to be just one row per sample
  semi_join(samples, by = c("sample" = "rownames(comp_lab)")) #subset just those samples
#with non-zero values

#Now we can run adonis, which asks whether there are differences in community
#composition of diet species between the two treatment levels of sterilization
#dist = jaccard with binary = TRUE because this is presence-absence data
adonis(comp_lab ~ Sterilized, data = meta_lab, dist = "jaccard", binary = TRUE)
#there is no significnat difference between treatment groups 

###########################
# Adonis abundance ####
############################
#Now we can repeat this with abundance data next.

#Adonis requires a matrix-type object with species as column names and 
#samples/sites as rows, which the following pipe does for the lab df.
#sites with zero across break adonis, so we need to remove them first
abund_lab <- lab %>%
  group_by(sample) %>%
  filter(sum(reads) >0) %>%
  ungroup() %>%
  dplyr::select(reads, sample, Family_ncbi) %>% #select only the variables for matrix
  pivot_wider(names_from = Family_ncbi, values_from = reads) %>% 
  #make column names based on unique ID, values from presence in cells
  column_to_rownames(var = "sample") 
#set the column names to the sample, so that only numeric values are in the matrix

#adonis requires a matrix as input, so we will convert df to that
abund_lab <- data.matrix(abund_lab, rownames.force = TRUE)

#Now that those are removed, we need to remove them from the eventual metadataframe
#we are creating next, so I will create a dataframe of those sample names to subset
#that dataframe
samples <- as.data.frame(rownames(abund_lab))

#Adonis also requires a meta-data DF for the analysis based on the factors of interest
#for the analysis. In this analysis, we are interested in comparing surface sterilized
#to unsterilized, so we can subset these from our DF as well

meta_lab <- lab %>%
  dplyr::select(sample, Sterilized) %>% #select metadata on sterilization and sample
  distinct(sample, Sterilized) %>% #reduce long DF to be just one row per sample
  semi_join(samples, by = c("sample" = "rownames(abund_lab)")) #subset just those samples
#with non-zero values

#Now we can run adonis, which asks whether there are differences in community
#composition of diet species between the two treatment levels of sterilization
#dist = bray because this is abundance data
adonis(abund_lab ~ Sterilized, data = meta_lab, dist = "bray")
#there is a marginally significant difference in abundance-based
#community... HMM...

###########################
# Heat map visualization of both presence and abundance ####
############################
heat_lab <- lab

#make order and ID values characters for ifelse sorting to order level below
heat_lab$Family_ncbi <- as.character(heat_lab$Family_ncbi)


#now we can summarise heat by order and sterilization treatment for
#the heatmap graph
heat_lab <- heat_lab %>%
  group_by(Family_ncbi, Sterilized) %>% #group by order and sterilization treatment
  summarise(reads = mean(reads)) %>% #then find the mean abundance for each order
  ungroup() %>% #ungroup so we can make order a factor
  mutate(Order = as.factor(Family_ncbi)) %>% #make order a factor
  pivot_wider(names_from = Sterilized, values_from = reads) %>% #make wide for two columns
  #for sterilization treatment so we can arrange in descending order of abundance
  arrange((NS + SS), (NS)) %>% #arrange in descending order of abundance
  mutate(Family_ncbi = factor(Family_ncbi, levels = Family_ncbi)) %>% #reset the levels of order on this
  #new abundance-based ordering of the Order factor
  gather(Sterilized, reads, NS:SS) #gather DF back up for visualization

#Because read abundance has such a broad range, we want to set a quantile-based
#variable for these abundances for the color ramp in the graph, but also want
#to discount the effects of zero reads in this, so we will create a DF of non-zero
#values for reads, and then find quantiles of this to inform our quantiles for the figure
heat_lab_nz <- heat_lab %>%
  filter(reads > 0)
quantile(heat_lab_nz$reads)
# 0.09090909   4.81818182  22.43750000  74.06818182 698.75000000  

#now we can set a quantile variable in our DF
heat_lab$quantiles <- ifelse(heat_lab$reads == 0, 0,
                         ifelse(heat_lab$reads > 0 & heat_lab$reads <= 0.09090909 , 1, 
                                ifelse(heat_lab$reads > 0.09090909 & heat_lab$reads <= 4.81818182, 2,
                                       ifelse(heat_lab$reads > 4.81818182 & heat_lab$reads <= 22.43750000, 3,
                                              ifelse(heat_lab$reads > 22.43750000 & heat_lab$reads <= 74.06818182, 4, 5)))))

#making it a factor for visualization
heat_lab$quantiles <- as.factor(heat_lab$quantiles)  

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

heat_lab_graph <- ggplot(heat_lab, aes(x = Sterilized, y = Family_ncbi, fill=quantiles, height = 0.95, width = 0.95)) +
  geom_tile() + 
  coord_equal() +
  labs(x = "Surface sterilization treatment", y = "Diet family") +
  scale_x_discrete(labels=c("NS" = "Not Sterilized", "SS" = "Surface Sterilized")) +
  scale_fill_manual(name = "Mean read abundance\n(divided by quantiles)",
                    values = pal2,
                    limits = names(pal2),
                    labels = c("0", "< 0.1", "< 4.8", "< 22.4", "< 74.1", "> 74.1")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 0.09090909   4.81818182  22.43750000  74.06818182 698.75000000   

###########################
# Effect Sizes Presence graph ####
############################

#Effect sizes are built off of means and either SD or SE
effect <- lab %>%
  group_by(Family_ncbi, Sterilized) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  summarize(mean = mean(presence), sd = sd(presence)) %>%
  ungroup() %>%
  group_by(Family_ncbi) %>%
  pivot_wider(names_from = c(Sterilized),
              values_from = c(mean, sd)) %>% #this pivots so there is an average for
  #SS and NS groups
  mutate(se_NS = sd_NS/sqrt(19), se_SS = sd_SS/sqrt(18)) #this computs the SE for groups

#this is the overall presence of each species, which we will use to sort the 
#graph visualization
pres_sort <- lab %>% 
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  group_by(Family_ncbi) %>%
  summarise(overall = mean(presence, na.rm=TRUE)) #gets the overall presence of that diet item

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
pres_effect_l <- ggplot(effects1, aes(x = IDs, y = Hedges_g)) +
  geom_point() +theme_bw() +geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "UNOISE3 effect size of average presence") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + coord_flip()

###########################
#Effect Sizes Abundance graph ####
############################

#Same as above, just with abundance data instead
#This is a graph of the effect sizes of abundance between surface sterilized 
#and non-surface sterilized predator individuals

#Effect sizes are built off of means and either SD or SE
effect_abund <- lab %>%
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
abund_sort <- lab %>% 
  group_by(Family_ncbi) %>%
  summarise(overall = mean(reads, na.rm=TRUE)) #gets the overall average reads of that diet item

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
abund_effect_l <- ggplot(effects_abund1, aes(x = IDs, y = Hedges_g)) +
  geom_point() +theme_bw() +geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "UNOISE3 effect size of average abundance") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + coord_flip()

###########################
#Plot grid for publication ####
############################
plot_grid(pres_effect_l, abund_effect_l)
