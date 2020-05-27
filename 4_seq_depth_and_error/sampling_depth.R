#Sampling Depth
#Ana Miller-ter Kuile
#February 16, 2020

#this is code for looking at sampling depth across samples to ensure that I've
#sufficiently sampled each sample in this dataset. 
#this code looks at sampling depth in community matrices created by both dada2 and unoise3

#sampling depth using iNEXT, first WITH predator reads, then without (first is to
#determine actual sequencing depth, second is to determine sampling depth of prey; first
#is bioinformatics, second is ecological question)
#Packages ####
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
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(effects)

#UNOISE sequencing depth####
#needed here: ASV matrix by samples minus the controls and asv column
#import data matrix
u3_comm <- read.csv(here("data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = "\t")

#rename columns for simplicity
colnames(u3_comm) <- sapply(str_split(colnames(u3_comm), "S"), function(x){return(x[[1]])})
u3_comm <- rename(u3_comm, "ASV" = "X.OTU.ID")

#set row names to ASV labels
rownames(u3_comm) <- u3_comm$ASV
#select samples minus controls and the ASV column
u3_comm_depth <- u3_comm %>%
  dplyr::select(-ASV, -CL1, -CL4, -NEG, -QC1)

#remove any ASVs that have zeros across all samples (probably from removing control)
u3_comm_depth <- u3_comm_depth[rowSums(u3_comm_depth) != 0,] #3 ASVs disappeared

#this determines sequencing depth across all samples
u3_seq_depth <- iNEXT(u3_comm_depth, q=0, datatype="abundance") #this determines sequencing depth for each sample

#u3_seq_depth$AsyEst
u3_seq_depth$DataInfo$SC
u3_seq_depth$DataInfo
#seq_depth$iNextEst
#seq_depth$AsyEst
#seq_depth$iNextEst
#graph the interpolated and extrapolated sampling depth per sample
u3_seq_depth_graph <- ggiNEXT(u3_seq_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "UNOISE3 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))
u3_seq_depth_graph
#DADA2 Sequencing depth ####
#needed here: ASV matrix by samples minus the controls and asv column
#import data matrix
d2_comm <- read.csv(here("data", "ASV_tables", "dada2_uc_asv_tab.tsv"), sep = "\t")

#rename columns for simplicity
colnames(d2_comm) <- sapply(str_split(colnames(d2_comm), "_"), function(x){return(x[[1]])})
colnames(d2_comm) <- str_remove(colnames(d2_comm), "\\.")

#set row names to ASV labels
rownames(d2_comm) <- d2_comm$X
#select samples minus controls and the ASV column
d2_comm_depth <- d2_comm %>%
  dplyr::select(-X, -CL12, -CL42, -QC1)

#remove any ASVs that have zeros across all samples (probably from removing control)
d2_comm_depth <- d2_comm_depth[rowSums(d2_comm_depth) != 0,] #3 ASVs disappeared

#this determines sequencing depth across all samples
d2_seq_depth <- iNEXT(d2_comm_depth, q=0, datatype="abundance") #this determines sequencing depth for each sample

d2_seq_depth$DataInfo$SC
d2_seq_depth$DataInfo
#seq_depth$iNextEst
#seq_depth$AsyEst
#seq_depth$iNextEst
#graph the interpolated and extrapolated sampling depth per sample
d2_seq_depth_graph <- ggiNEXT(d2_seq_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))
d2_seq_depth_graph

#Ster vs. unster. seq depth ####
#for either DADA2 OR UNOISE3, are sterilized spiders sequencing depth different from unsterilized?
u3_dep <- u3_seq_depth$DataInfo
d2_dep <- d2_seq_depth$DataInfo

metadata <- read.csv(here("data", "Sample_Metadata.csv"))

u3_dep <- u3_dep %>%
  left_join(metadata, by =c("site" = "sample")) %>%
  mutate(pipeline = "u3")

all_dep <- d2_dep %>%
  left_join(metadata, by =c("site" = "sample")) %>%
  mutate(pipeline = "d2") %>%
  bind_rows(u3_dep)

all_lab_dep <- all_dep %>%
  filter(Source == "LAB")

all_fld_dep <- all_dep %>%
  filter(Source == "FIELD")

mod_l <- glmmTMB(n ~ Sterilized + (1|pipeline),
               data = all_lab_dep,
               family = "genpois",
               REML = FALSE)
hist(residuals(mod_l))
plot(residuals(mod_l))
mod_l_null <- glmmTMB(n ~ 1 + (1|pipeline),
                    data = all_lab_dep,
                    family = "genpois",
                    REML = FALSE)

AICc(mod_l, mod_l_null)

mod_l <- glmmTMB(n ~ Sterilized + (1|pipeline),
               data = all_lab_dep,
               family = "genpois")

hist(residuals(mod_l))
plot(residuals(mod_l))
summary(mod_l)

simulationOutput_depth <- simulateResiduals(fittedModel = mod_l) 
depth_fit <- plot(simulationOutput_depth) #look okay
depth_zi <- testZeroInflation(simulationOutput_depth) #not zero inflated
depth_od <- testDispersion(simulationOutput_depth)

mod_f <- glmmTMB(n ~ Sterilized + (1|pipeline),
                 data = all_fld_dep,
                 family = "nbinom2",
                 REML = FALSE)
hist(residuals(mod_f))
plot(residuals(mod_f))
mod_f_null <- glmmTMB(n ~ 1 + (1|pipeline),
                      data = all_fld_dep,
                      family = "nbinom2",
                      REML = FALSE)

AICc(mod_f, mod_f_null)

mod_f <- glmmTMB(n ~ Sterilized + (1|pipeline),
                 data = all_fld_dep,
                 family = "nbinom1")

hist(residuals(mod_f))
plot(residuals(mod_f))
summary(mod_f)

simulationOutput_depth <- simulateResiduals(fittedModel = mod_f) 
depth_fit <- plot(simulationOutput_depth) #look okay
depth_zi <- testZeroInflation(simulationOutput_depth) #not zero inflated
depth_od <- testDispersion(simulationOutput_depth)

ggplot(all_dep, aes(x = Sterilized, y = n)) +
  geom_boxplot() + theme_bw() +
  facet_wrap (~Source) +
  labs(x = "Surface Sterilized", y = "Sequencing depth",
       title = "Sequencing depth by surface sterilization after sequences assigned to ASVs")

#prey only: first attach taxonomies####
#now do sampling depth analyses of JUST prey ASVs, which means I will have to subset
#only prey ASVs from this based on the prey taxonomies for both bold and ncbi
#d2_comm_depth and u3_comm_depth are the community matrices to bind to
#DADA2#
#import both NCBI and BOLD taxonomies
d2_ncbi <- read.csv(here("6_taxonomic_assignment", "MEGAN", "dada2_unclean",
                         "dada2_uc_ncbi_taxa.csv"))
#this is everything that got a taxonomic assignment through MEGAN
d2_all <- read.csv(here("6_taxonomic_assignment", "MEGAN", "dada2_unclean",
                        "dada2_uc_all.csv"))
d2_bold <- read.csv(here("6_taxonomic_assignment", "BOLD", "dada2_uc", 
                         "dada2_uc_taxa_bold.csv"))
d2_bold <- rename(d2_bold, "ASV" = "Query.ID")

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it
#sort into categories for bold:
d2_bold$taxonomy <- ifelse(
  d2_bold$ID == "Heteropoda venatoria", "predator",
  ifelse(d2_bold$Kingdom == "Fungi", "non-diet", 
         ifelse(d2_bold$Phylum == "Arthropoda" & d2_bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
d2_ncbi$taxonomy <- ifelse(
  d2_ncbi$ID == "Sparassidae" | d2_ncbi$ID == "Heteropoda venatoria", "predator",
  ifelse(d2_ncbi$Phylum == "Arthropoda" & d2_ncbi$ID != "Sparassidae" & d2_ncbi$ID != "Heteropoda venatoria", 
         "prey", "non-diet"
  ))

#this is an anti-join to subset from ALL MEGAN hits for those that aren't within
#diet categories, but which still got assignments in MEGAN
d2_nond <- d2_all %>%
  anti_join(d2_ncbi, by = "ASV")
#this gives these all the category of "non-diet" to distinguish from those that didn't
#get any taxonomic assignment at all. 
d2_nond$taxonomy <- d2_nond$Category

#thus, the taxonomy variable has 4 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table)
#subset only the variables of interest from these two dataframes before merging
d2_bold_id <- d2_bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

d2_ncbi_id <- d2_ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
d2_id_both <- d2_bold_id %>%  
  full_join(d2_ncbi_id, by = c("ASV", "taxonomy"))

#join these both again with the non-diet hits
d2_id_all <- d2_id_both %>%
  full_join(d2_nond, by = c("ASV", "taxonomy"))

d2_id_all <- d2_id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y,
         ID_nondiet = ID,
         Order_nondiet = Order)

d2_taxonomy <- colnames(d2_id_all)

#Prey only seq depth####
#this binds them all to the community matrix for analyses
#and then filters only prey, and then removes all of the rows from the 
#taxonomy dataframe
d2_comm_depth_p <- d2_comm_depth %>%
  rownames_to_column(var = "ASV") %>%
  left_join(d2_id_all, by = "ASV") %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(-all_of(d2_taxonomy))
  
#UNOISE3##
#import both NCBI and BOLD taxonomies
u3_ncbi <- read.csv(here("6_taxonomic_assignment", "MEGAN", "unoise_unclean",
                         "unoise_uc_ncbi_taxa.csv"))
u3_bold <- read.csv(here("6_taxonomic_assignment", "BOLD", "usearch_uc", 
                         "unoise_uc_bold_taxa.csv"))
u3_bold <- rename(u3_bold, "ASV" = "Query.ID")

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it
#sort into categories for bold:
u3_bold$taxonomy <- ifelse(
  u3_bold$ID == "Heteropoda venatoria", "predator",
  ifelse(u3_bold$Kingdom == "Fungi" | u3_bold$Class == "Mammalia", "non-diet", 
         ifelse(u3_bold$Phylum == "Arthropoda" | u3_bold$Class == "Reptilia" & u3_bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
u3_ncbi$taxonomy <- u3_ncbi$Category

#thus, the taxonomy variable has 4 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table)
#subset only the variables of interest from these two dataframes before merging
u3_bold_id <- u3_bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

u3_ncbi_id <- u3_ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
u3_id_all <- u3_bold_id %>%  
  full_join(u3_ncbi_id, by = c("ASV", "taxonomy"))

u3_id_all <- u3_id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y)

#make a vector for following pipe to draw from
u3_taxonomy <- colnames(u3_id_all)

#this binds them all to the community matrix for analyses
#and then filters only prey, and then removes all of the rows from the 
#taxonomy dataframe
u3_comm_depth_p <- u3_comm_depth %>%
  rownames_to_column(var = "ASV") %>%
  left_join(u3_id_all, by = "ASV") %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(-all_of(u3_taxonomy))

#Now repeat sequencing depth with just prey
#this determines sequencing depth across all samples
u3_comm_depth_p1 <- u3_comm_depth_p[,colSums(u3_comm_depth_p) > 0] #removed 1 with zero counts

u3p_seq_depth <- iNEXT(u3_comm_depth_p1, q=0, datatype="abundance") #this determines sequencing depth for each sample

u3p_seq_depth$DataInfo$SC
u3_seq_depth$DataInfo
u3p_seq_depth$iNextEst
u3p_seq_depth$AsyEst
#seq_depth$iNextEst
#graph the interpolated and extrapolated sampling depth per sample
u3p_seq_depth_graph <- ggiNEXT(u3p_seq_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))
u3p_seq_depth_graph

###
d2_comm_depth_p1 <- d2_comm_depth_p[,colSums(d2_comm_depth_p) > 0] #removed 11 with zero counts

d2p_seq_depth <- iNEXT(d2_comm_depth_p1, q=0, datatype="abundance") #this determines sequencing depth for each sample

d2p_seq_depth$DataInfo$SC
d2p_seq_depth$DataInfo
#seq_depth$iNextEst
#seq_depth$AsyEst
#seq_depth$iNextEst
#graph the interpolated and extrapolated sampling depth per sample
d2p_seq_depth_graph <- ggiNEXT(d2p_seq_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))
d2p_seq_depth_graph

#RAW READ ABUNDANCE ster vs. unster for prey only ####
#for either DADA2 OR UNOISE3, are sterilized spiders sequencing depth different from unsterilized?
u3p_dep <- u3p_seq_depth$DataInfo
d2p_dep <- d2p_seq_depth$DataInfo

metadata <- read.csv(here("data", "Sample_Metadata.csv"))

u3p_dep <- u3p_dep %>%
  full_join(metadata, by =c("site" = "sample")) %>%
  mutate(pipeline = "u3") %>%
  replace_na(list(n = 0)) 

d2p_dep <- d2p_dep %>%
  full_join(metadata, by =c("site" = "sample")) %>%
  mutate(pipeline = "d2")  %>%
  replace_na(list(n = 0))

allp_dep <- d2p_dep %>%
  bind_rows(u3p_dep)

u3p_lab_dep <- u3p_dep %>%
  filter(Source == "LAB")

u3p_fld_dep <- u3p_dep %>%
  filter(Source == "FIELD")

d2p_lab_dep <- d2p_dep %>%
  filter(Source == "LAB")

d2p_fld_dep <- d2p_dep %>%
  filter(Source == "FIELD")
###
mod_p <- glm.nb(n ~ Sterilized*Source + pipeline,
             data = allp_dep)

plot(allEffects(mod_p))
simulationOutput_depth <- simulateResiduals(fittedModel = mod_p) 
depth_fit <- plot(simulationOutput_depth)
depth_zi <- testZeroInflation(simulationOutput_depth)
depth_od <- testDispersion(simulationOutput_depth)
summary(mod_p)
library(emmeans)
model.emm_depth <- emmeans(mod_p, pairwise ~ Sterilized | Source)
pairs(model.emm_depth)

#for all predators (mesocosm and field), sterilization is an important predictor
#in our model. However, sterilization does not significantly change prey 
#sequencing depth for field spiders, BUT it does significantly increase prey
#read abundance for mesocosm spiders. 
sum_pdep <- allp_dep %>%
  group_by(Sterilized, Source) %>%
  summarize(mean = mean(n), n(), se = sd(n)/sqrt(`n()`)) 

site.labs <- c("Field-collected", "Mesocosm")
names(site.labs) <- c("FIELD", "LAB")

ggplot(sum_pdep, aes(x = Sterilized, y = mean)) +
  geom_bar(stat = "identity") +theme_bw() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  facet_wrap (~Source,labeller = labeller(Source = site.labs)) +
  scale_y_log10()+
  labs(x = "Surface Sterilization Treatment", y = "Sequencing depth",
       title = "Sequencing depth by surface sterilization after sequences assigned to ASVs") +
  scale_x_discrete(labels = c("NS" = "Not S. sterilized", "SS" = "Surface sterilized"))
?scale_x_discrete

