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
                  "esc")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###########################
# Load data ####
############################

ASV <- read.csv(here("data", 
                     "outputs", 
                     "rarefied_taxonomic_sort", 
                     "field_prey_ASVs_rare.csv"))

###########################
# Field Richness analysis ####
############################

richness <- ASV %>%
  group_by(sample, Sterilized) %>%
  summarise(SR = sum(reads > 0, na.rm=TRUE)) #number of ASVs with greater than 0 reads

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
simulationOutput <- simulateResiduals(fittedModel = rich_null, plot= T) 
testZeroInflation(simulationOutput) 
testDispersion(simulationOutput) 

###########################
# Plot for paper of richness ####
############################
(rich_graph <- ggplot(richness, aes(x = Sterilized, y = SR)) +
   geom_boxplot(fill = "#F29979") + theme_bw() +
   labs(x = "Surface sterilization treatment", y = "ASV richness") +
   scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
   theme(legend.position = "none"))

richness %>%
  ungroup() %>%
  summarise(mean= mean(SR), total = n(), sd = sd(SR), se = sd/sqrt(total))
max(richness$SR)
min(richness$SR)
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

#pretty ggplot plot 
#site metadata
sites <- meta %>%
  rownames_to_column("site")
#get the site (sample) scores out and attach to site metadata
CCAscores <- scores(cca, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") %>%
  left_join(sites, by = "site")

#get the vectors out representing the loadings by species
CCAvect <- scores(cca, display = c("cn")) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") 

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = CCAscores, aes(x = CCA1, y = CA1, color = Sterilized), size = 3) +
  geom_segment(data = CCAvect, 
               aes(x = 0, y = 0, xend = CCA1, yend = CA1), 
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  geom_text(data = CCAvect, aes(x = CCA1, y = CA1, label = ID), 
            nudge_y = -0.1, nudge_x = -0.15, size = 5) +
  theme_bw() +
  labs(x = "CCA1 (1.8%)",
       y = "CA1 (9.0%)") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) 


###########################
# (OPTIONAL: adonis instead of GLMM for field presence) ####
############################
#Now we can run adonis, which asks whether there are differences in community
#composition of diet species between the two treatment levels of sterilization
#dist = jaccard with binary = TRUE because this is presence-absence data
adonis(mat_pres ~ Sterilized, data = meta, dist = "jaccard", binary = TRUE)
#there is no significnat difference between treatment groups 


# Heat Map ----------------------------------------------------------------
heat <- ASV

heat_nz <- heat %>%
  filter(reads > 0)
quantile(heat_nz$reads)
#1    1    5   20.5 4079  

#now we can summarise heat by order and sterilization treatment for
#the heatmap graph
ASV_ord <- heat %>%
  group_by(ASV) %>%
  summarise(reads = sum(reads)) %>%
  arrange(reads) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  dplyr::select(ASV) %>%
  as_vector

#now we can set a quantile variable in our DF
ASV %>%
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
                                ifelse(heat$reads > 1 & heat$reads <= 5, 3,
                                       ifelse(heat$reads > 6 & heat$reads <= 20.5, 4, 5))))

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

heat <- heat %>%
  mutate(ASV = as.factor(ASV)) %>%
  mutate(ASV = factor(ASV, levels = ASV_ord))

(heat_map <- ggplot(heat, aes(x = sample, y = ASV, fill=quantiles)) +
    geom_tile(color = "black") + 
    coord_equal() +
    labs(x = "Sample", y = "ASV") +
    scale_fill_manual(name = "Mean read abundance\n(divided by quantiles)",
                      values = pal3,
                      limits = names(pal3),
                      labels = c("0", "< 1", "< 5", "< 20", "> 20")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept = 17.5))


