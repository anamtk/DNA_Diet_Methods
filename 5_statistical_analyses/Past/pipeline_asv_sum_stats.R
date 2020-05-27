#ASV table Summaries
#Ana Miller-ter Kuile
#February 14, 2020

#This script will look at the ASV/ZOTU tables created by each of my pipelines for prey
#items (dada2 and unoise with clean and unclean) and summarize each based on items of
#interest, including:
#1. positive and negative control mapping

#2. total number of ASVs for each pipeline/individual (uc, prey, predator, unmapped)
#both TOTAL and per sample

#the following require binding the combined MEGAN/BOLD taxonomy ID dataframe:

#3. total number of prey ASVs for each pipeline
#both TOTAL and per sample

#4. total prey ASV read abundance (raw)
#both TOTAL and per sample

#5. proportion of prey to other ASVs (number of ASVs)

#6. proportion of prey reads to other ASVs (raw read abundance)
###a. consider here sequencing depth? Ask Austen

#7. Predator/prey ASV ratios

#8. predator/prey read abundance (raw)

#9. amount of known diet item for lab-fed (raw reads)

#10. amount of knwon diet item for lab-fed (num of ASVs)

#11. Phylogenetic diversity of prey ASVs

#Load required packages####
library(here)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(rcompanion)
library(effects)
library(lattice)
library(DHARMa)
library(MuMIn)
library(glmmTMB)
library(picante)
library(ape)
library(phytools)
library(phyloseq)


#Import ASV tables####
d2_uc <- read.csv(here("data", "ASV_tables", "dada2_uc_asv_tab.tsv"), sep = "\t")
d2_c <- read.csv(here("data", "ASV_tables", "dada2_c_prey_asv_tab.tsv"), sep = '\t')
d2_um <- read.csv(here("data", "ASV_tables", "dada2_c_um_asv_tab.tsv"), sep = '\t')
d2_pred <- read.csv(here("data", "ASV_tables", "dada2_c_pred_asv_tab.tsv"), sep = '\t')
u3_uc <- read.csv(here("data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = '\t')
u3_c <- read.csv(here("data", "ASV_tables", "unoise_c_prey_zotu_tab.txt"), sep = '\t')
u3_um <- read.csv(here("data", "ASV_tables", "unoise_c_um_zotu_tab.txt"), sep = '\t')
u3_pred <- read.csv(here("data", "ASV_tables", "unoise_c_pred_zotu_tab.txt"), sep = '\t')

#rename X to ASV in all these tables.
d2_uc <- rename(d2_uc, "ASV" = "X")
d2_c <- rename(d2_c, "ASV" = "X")
d2_um <- rename(d2_um, "ASV" = "X")
d2_pred <- rename(d2_pred, "ASV" = "X")
u3_uc <- rename(u3_uc, "ASV" = "X.OTU.ID")
u3_c <- rename(u3_c, "ASV" = "X.OTU.ID")
u3_um <- rename(u3_um, "ASV" = "X.OTU.ID")
u3_pred <- rename(u3_pred, "ASV" = "X.OTU.ID")

#remove "prey_ref_" and "pred_ref_" from the sample names in dada2 datasets:
for ( col in 1:ncol(d2_c)){
  colnames(d2_c)[col] <-  sub("prey_ref_*", "", colnames(d2_c)[col])
}

for ( col in 1:ncol(d2_pred)){
  colnames(d2_pred)[col] <-  sub("predator_ref_*", "", colnames(d2_pred)[col])
}

#rename all the sample names across dataframes for consistency
colnames(d2_uc) <- sapply(str_split(colnames(d2_uc), "_"), function(x){return(x[[1]])})
colnames(d2_uc) <- str_remove(colnames(d2_uc), "\\.")

d2_uc <- rename(d2_uc, "CL1" = "CL12")
d2_uc <- rename(d2_uc, "CL4" = "CL42")
colnames(d2_pred) <- sapply(str_split(colnames(d2_pred), "_"), function(x){return(x[[1]])})
colnames(d2_pred) <- str_remove(colnames(d2_pred), "\\.")

colnames(d2_c) <- sapply(str_split(colnames(d2_c), "_"), function(x){return(x[[1]])})
colnames(d2_c) <- str_remove(colnames(d2_c), "\\.")

colnames(d2_um) <- sapply(str_split(colnames(d2_um), "_"), function(x){return(x[[1]])})
colnames(d2_um) <- str_remove(colnames(d2_um), "\\.")

d2_um <- rename(d2_um, "CL1" = "CL12")
d2_um <- rename(d2_um, "CL4" = "CL42")
colnames(u3_uc) <- sapply(str_split(colnames(u3_uc), "S"), function(x){return(x[[1]])})
u3_uc <- rename(u3_uc, "ASV" = "A")

#these are all fine
#colnames(u3_c)
#colnames(u3_um)
#colnames(u3_pred)

#1. positive and negative control mapping ####
#dada2 uncleaned
cl1_d2uc <- d2_uc %>%
  summarize(count = sum(CL1 != 0)) %>% #1
  mutate(pipeline = "d2uc", control = "CL1") 
cl4_d2uc <- d2_uc %>%
  summarize(count = sum(CL4 != 0)) %>% 
  mutate(pipeline = "d2uc", control = "CL4") #1
qc1_d2uc <- d2_uc %>%
  summarize(count = sum(QC1 != 0)) %>%
  mutate(pipeline = "d2uc", control = "CL4") #1
#NEG is ZERO
neg_d2uc <- as.data.frame(0)
neg_d2uc <- neg_d2uc %>%
  mutate(pipeline = "d2uc", control = "NEG") %>%
  rename("count" = "0")


#usearch uncleaned
cl1u3uc <- u3_uc %>%
  summarize(count = sum(CL1 != 0)) %>% #3
  mutate(pipeline = "u3uc", control = "CL1")
cl4u3uc <- u3_uc %>%
  summarize(count = sum(CL4 != 0)) %>% #3
  mutate(pipeline = "u3uc", control = "CL4")
qc1u3uc <- u3_uc %>%
  summarize(count = sum(QC1 != 0)) %>% #3
  mutate(pipeline = "u3uc", control = "QC1")
negu3uc <- u3_uc %>%
  summarize(count = sum(NEG != 0)) %>% #1
  mutate(pipeline = "u3uc", control = "NEG")

#both cleaned map to zero
d2ccont <- as.data.frame(0)
d2ccont <- d2ccont %>%
  mutate(pipeline = "d2c", control = "ALL") %>%
  rename("count" = "0")

u3ccont <- as.data.frame(0)
u3ccont <- u3ccont %>%
  mutate(pipeline = "u3c", control = "ALL") %>%
  rename("count" = "0")
#bind all these together into a dataframe
control_counts <- bind_rows(cl1_d2uc, cl4_d2uc, qc1_d2uc,
                            neg_d2uc, cl1u3uc, cl4u3uc, qc1u3uc,
                            negu3uc, d2ccont, u3ccont)

control_count_plot <- ggplot(control_counts, aes(x = pipeline, y = count)) +
    geom_boxplot() +
    labs(x = "Pipeline", y = "Number of ASVs per control", title = "Number of ASVs mapped to each positive and negative control") +
    theme_bw()

#controls drop out of the cleaned datasets
#and they drop out of these as well, and negative goes to predator with one ASV
colnames(d2_pred)
colnames(u3_pred)
u3_pred %>%
  summarize(sum(NEG != 0)) #1

#no new ASVs get mapped in the cleaned out unmapped sequences 
u3_um %>%
  summarize(sum(CL1 != 0))#3
u3_um %>%
  summarize(sum(CL4 != 0))#3
u3_um %>%
  summarize(sum(QC1 != 0))#3

d2_um %>%
  summarize(sum(CL1 != 0))#1
d2_um %>%
  summarize(sum(CL4 != 0))#1
d2_um %>%
  summarize(sum(QC1 != 0))#1
#Summary: cleaning does not add any ASVs to the controls, and dada2 maps better to controls

#2. total number of ASVs for each pipeline/individual (uc, prey, predator, unmapped) ####
lab_spiders <- metadata %>%
  select(sample, Source) %>%
  filter(Source == "LAB") %>%
  select(sample)
lab_spiders <- as_vector(as.character(lab_spiders))
  
fld_spiders <- metadata %>%
  select(sample, Source) %>%
  filter(Source == "FIELD") %>%
  select(sample)
fld_spiders <- as_vector(fld_spiders)
fld_spiders <- as.character(fld_spiders)
#dada2 uc
nrow(d2_uc) #214 total ASVs
d2_uc_l <- d2_uc %>%
  select(lab_spiders)
d2_uc_f <- d2_uc %>%
  select(fld_spiders)
#by sample:
d2_ASVct <- d2_uc %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_d2uc, ASV:QC1, factor_key = TRUE) 
d2_lASVct <- d2_uc_l %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_d2uc, HEV07:HEV29, factor_key = TRUE) 
d2_fASVct <- d2_uc_f %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_d2uc, HEV65:HEV100, factor_key = TRUE) 
#unoise3 uc
nrow(u3_uc) #176 total ASVs
#by sample:
u3_ASVct <- u3_uc %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_u3uc, ASV:QC1, factor_key = TRUE)

#dada2 c
nrow(d2_c) #138 total ASVs
#by sample:
d2c_ASVct <- d2_c %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_d2c, ASV:HEV89, factor_key = TRUE)

#unoise3 c
nrow(u3_c) #99 total ASVs
#by sample:
u3c_ASVct <- u3_c %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_u3c, ASV:HEV99, factor_key = TRUE)

#join these two together
by_sample <- d2_ASVct %>%
  full_join(d2c_ASVct, by = "sample") %>%
  full_join(u3_ASVct, by = "sample") %>%
  full_join(u3c_ASVct, by = "sample")

metadata <- read.csv(here("data", "Sample_Metadata.csv"))

by_sample_ster <- by_sample %>%
  left_join(metadata, by = "sample")

ASV_tot_ster <- by_sample_ster %>%
  dplyr::select(ASVs_d2uc, ASVs_d2c, ASVs_u3uc, ASVs_u3c, Sterilized, Source) %>%
  gather(pipeline, value, ASVs_d2uc:ASVs_u3c) %>%
  group_by(pipeline, Sterilized)
#total ASVs
total_ASVs <- by_sample %>%
  filter(sample == "ASV") %>%
  gather(pipeline, value, ASVs_d2uc:ASVs_u3c)

total_asv_graph <- ggplot(ASV_tot_ster, aes(x = pipeline, y = value, fill = Sterilized)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Pipeline", y = "Total ASVs", 
       title = "Total ASVs produced for each pipeline") + theme_bw() +
  scale_fill_manual(values= c("#7B8D65","#F29979"))


#and remove the control samples for glmm below
by_sample <- by_sample %>%
  filter(!sample %in% c("ASV", "CL1", "CL4", "QC1", "NEG"))

by_sample_long <- by_sample %>%
  gather(pipeline, value, ASVs_d2uc:ASVs_u3c, factor_key=TRUE)

#join with sample metadata so we can split by field and lab 
#and then use sterilization in the models 
by_sample_long <- by_sample_long %>%
  left_join(metadata, by = "sample")

#create binary sterilized factor
by_sample_long$Ster <- ifelse(by_sample_long$Sterilized == "NS", 0, 1)
by_sample_long$Ster <- as.factor(as.character(by_sample_long$Ster))

#split field and lab
f_by_sample_long <- by_sample_long %>%
  filter(Source == "FIELD")

l_by_sample_long <- by_sample_long %>%
  filter(Source == "LAB")

#FIELD models
tot_ASV_mod_f <- glmmTMB(value ~ pipeline*Ster + (1|sample),
                       data = f_by_sample_long,
                       family = "genpois",
                       REML = FALSE)

tot_ASV_mod2_f <- glmmTMB(value ~ pipeline + Ster + (1|sample),
                       data = f_by_sample_long,
                       family = "genpois",
                       REML = FALSE)

tot_ASV_mod3_f <- glmmTMB(value ~ Ster + (1|sample),
                       data = f_by_sample_long,
                       family = "genpois", 
                       REML = FALSE)

tot_ASV_mod4_f <- glmmTMB(value ~ pipeline + (1|sample),
                       data = f_by_sample_long,
                       family = "genpois",
                       REML = FALSE)

tot_ASV_mod_null_f <- glmmTMB(value ~ 1 + (1|sample),
                       data = f_by_sample_long,
                       family = "genpois", 
                       REML = FALSE)

#assess the corrected AIC for these models:
AICc(tot_ASV_mod_f, tot_ASV_mod2_f, tot_ASV_mod3_f, tot_ASV_mod4_f, tot_ASV_mod_null_f)
anova(tot_ASV_mod2_f, tot_ASV_mod4_f)
#it appears that a few are within 2 AIC values of each other, including model 1,2, and 4
#the most parsimonious of these is model 4, which only includes pipeline, though
#when we look at the next most parsimonious, mod2, sterilization is marginally significant
summary(tot_ASV_mod4_f)
summary(tot_ASV_mod2_f) #p-value of 0.0609 for Sterilization

tot_ASV_mod2_f <- glmmTMB(value ~ pipeline + Ster + (1|sample),
                          data = f_by_sample_long,
                          family = "genpois")

model.emm_totASV_f <- emmeans(tot_ASV_mod2_f, "pipeline")
model.emm_totASV_s_f <- emmeans(tot_ASV_mod2_f, "Ster")
pairs(model.emm_totASV_f)
pairs(model.emm_totASV_s_f) #marginal significance here.

x_totASV_f <- residuals(tot_ASV_mod2_f)
plotNormalHistogram(x_totASV_f)
#residuals look relatively normal
plot(residuals(tot_ASV_mod2_f)) #residuals look fairly homoskedastic
plot(allEffects(tot_ASV_mod2_f))
#Unoise unclean produced significantly more ASVs per sample than any other
#pipeline. cleaning significantly decreased the number of ASVs per sample for 
#both pipelines. Sterilization marginally reduced the number of ASVs per sample
#because the interaction was non-significant, the effect of sterilization
#on total ASVs did not differe by pipeline.

#checking other model assumptions here:
simulationOutput_totASV_f <- simulateResiduals(fittedModel = tot_ASV_mod2_f) 
totASVfit_f <- plot(simulationOutput_totASV_f, asFactor=TRUE) #these look good
totASVzi_f <- testZeroInflation(simulationOutput_totASV_f) #not zero inflated
totASVod_f <- testDispersion(simulationOutput_totASV_f) #not overdispersed

#and that boxplot here
ASV_sample_plot_f <- ggplot(f_by_sample_long, aes(x = pipeline, y = value, fill = Sterilized)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Pipeline", y = "ASVs per sample", 
       title = "Total number of ASVs per sample by pipeline and sterilization") +
  scale_fill_manual(values = c("#7B8D65","#F29979"))

#cleaned out
nrow(d2_um) #52
nrow(d2_pred) #38
nrow(u3_um) #72
nrow(u3_pred) #25

#all cleaned added
nrow(d2_c) + nrow(d2_um) + nrow(d2_pred) #228
nrow(u3_c) + nrow(u3_um) + nrow(u3_pred) #196

#Intermediate: Bind BOLD and NCBI taxonomies ####
#For some MEGAN assignments, I subset two groups - one with everything that is likely diet
#(which I assigned taxonomies to) and everything else into an ALL_taxa category, which
#included mostly fungi and some other non-diet items which will be under a non-diet cateogry
#for some other MEGAN assignments, the list was small enough that I just created one
#ALL_taxa file, which will get the same types of assignments. 
#The end goal of this section is to have a dataframe which matches ASVs in the community
#matrix for each pipeline (clean and unclean) to their taxonomic categories
#DADA2 UNCLEAN#
#import both NCBI and BOLD taxonomies
d2_uc_ncbi <- read.csv(here("6_taxonomic_assignment", "MEGAN", "dada2_unclean",
                            "dada2_uc_ncbi_taxa.csv"))
#this is everything that got a taxonomic assignment through MEGAN
d2_uc_all <- read.csv(here("6_taxonomic_assignment", "MEGAN", "dada2_unclean",
                           "dada2_uc_all.csv"))
d2_uc_bold <- read.csv(here("6_taxonomic_assignment", "BOLD", "dada2_uc", 
                            "dada2_uc_taxa_bold.csv"))
d2_uc_bold <- rename(d2_uc_bold, "ASV" = "Query.ID")

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it
#sort into categories for bold:
d2_uc_bold$taxonomy <- ifelse(
  d2_uc_bold$ID == "Heteropoda venatoria", "predator",
  ifelse(d2_uc_bold$Kingdom == "Fungi", "non-diet", 
         ifelse(d2_uc_bold$Phylum == "Arthropoda" & d2_uc_bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
d2_uc_ncbi$taxonomy <- ifelse(
  d2_uc_ncbi$ID == "Sparassidae" | d2_uc_ncbi$ID == "Heteropoda venatoria", "predator",
  ifelse(d2_uc_ncbi$Phylum == "Arthropoda" & d2_uc_ncbi$ID != "Sparassidae" & d2_uc_ncbi$ID != "Heteropoda venatoria", 
         "prey", "non-diet"
  ))

#this is an anti-join to subset from ALL MEGAN hits for those that aren't within
#diet categories, but which still got assignments in MEGAN
d2_uc_nond <- d2_uc_all %>%
  anti_join(d2_uc_ncbi, by = "ASV")
#this gives these all the category of "non-diet" to distinguish from those that didn't
#get any taxonomic assignment at all. 
d2_uc_nond$taxonomy <- d2_uc_nond$Category

#thus, the taxonomy variable has 5 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table), and unknown
#subset only the variables of interest from these two dataframes before merging
d2_uc_bold_id <- d2_uc_bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

d2_uc_ncbi_id <- d2_uc_ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
d2_uc_id_both <- d2_uc_bold_id %>%  
  full_join(d2_uc_ncbi_id, by = c("ASV", "taxonomy"))

#join these both again with the non-diet hits
d2_uc_id_all <- d2_uc_id_both %>%
  full_join(d2_uc_nond, by = c("ASV", "taxonomy"))

d2_uc_id_all <- d2_uc_id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y,
         ID_nondiet = ID,
         Order_nondiet = Order)

#this binds them all to the community matrix for analyses
d2_uc_id_tab <- d2_uc %>%
  left_join(d2_uc_id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
d2_uc_id_tab$taxonomy <- replace_na(d2_uc_id_tab$taxonomy, "no hit")

#write the taxonomies to a file for use in the community dissimilarity scripts script
d2_uc_taxonomies <- d2_uc_id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, ID_nondiet, Order_nondiet, taxonomy)
write.csv(d2_uc_taxonomies, here("data", "d2_uc_tax_ass.csv"))

#DADA2 CLEAN##
#import both NCBI and BOLD taxonomies
d2_c_ncbi <- read.csv(here("6_taxonomic_assignment", "MEGAN", "dada2_clean",
                           "prey", "dada2_c_prey_ncbi_taxa.csv"))
#this is everything that got a taxonomic assignment through MEGAN
d2_c_all <- read.csv(here("6_taxonomic_assignment", "MEGAN", "dada2_clean",
                           "prey", "dada2_c_prey_ALL.csv"))
d2_c_bold <- read.csv(here("6_taxonomic_assignment", "BOLD", "dada2_c", 
                            "prey", "dada2_c_prey_taxa_bold.csv"))

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it
#sort into categories for bold:
d2_c_bold$taxonomy <- ifelse(
  d2_c_bold$ID == "Heteropoda venatoria", "predator",
  ifelse(d2_c_bold$Kingdom == "Fungi", "non-diet", 
         ifelse(d2_c_bold$Phylum == "Arthropoda" & d2_c_bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
d2_c_ncbi$taxonomy <- ifelse(
  d2_c_ncbi$ID == "Sparassidae" | d2_c_ncbi$ID == "Heteropoda venatoria", "predator",
  ifelse(d2_c_ncbi$Phylum == "Arthropoda" & d2_c_ncbi$ID != "Heteropoda venatoria" & d2_c_ncbi$ID != "Sparassidae", 
         "prey", "non-diet"
  ))

#this is an anti-join to subset from ALL MEGAN hits for those that aren't within
#diet categories, but which still got assignments in MEGAN
d2_c_nond <- d2_c_all %>%
  anti_join(d2_c_ncbi, by = "ASV")
#this gives these all the category of "non-diet" to distinguish from those that didn't
#get any taxonomic assignment at all. 
d2_c_nond$taxonomy <- d2_c_nond$Category

#thus, the taxonomy variable has 5 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table), and unknown
#subset only the variables of interest from these two dataframes before merging
d2_c_bold_id <- d2_c_bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

d2_c_ncbi_id <- d2_c_ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
d2_c_id_both <- d2_c_bold_id %>%  
  full_join(d2_c_ncbi_id, by = c("ASV", "taxonomy"))

#join these both again with the non-diet hits
d2_c_id_all <- d2_c_id_both %>%
  full_join(d2_c_nond, by = c("ASV", "taxonomy"))

d2_c_id_all <- d2_c_id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y,
         ID_nondiet = ID)

#this binds them all to the community matrix for analyses
d2_c_id_tab <- d2_c %>%
  left_join(d2_c_id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
d2_c_id_tab$taxonomy <- replace_na(d2_c_id_tab$taxonomy, "no hit")

#write the taxonomies to a file for use in the community dissimilarity scripts script
d2_c_taxonomies <- d2_c_id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, ID_nondiet, taxonomy)
write.csv(d2_c_taxonomies, here("data", "d2_c_tax_ass.csv"))

#UNOISE UNCLEAN##
#import both NCBI and BOLD taxonomies
u3_uc_ncbi <- read.csv(here("6_taxonomic_assignment", "MEGAN", "unoise_unclean",
                           "unoise_uc_ncbi_taxa.csv"))
u3_uc_bold <- read.csv(here("6_taxonomic_assignment", "BOLD", "usearch_uc", 
                           "unoise_uc_bold_taxa.csv"))
u3_uc_bold <- rename(u3_uc_bold, "ASV" = "Query.ID")

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it
#sort into categories for bold:
u3_uc_bold$taxonomy <- ifelse(
  u3_uc_bold$ID == "Heteropoda venatoria", "predator",
  ifelse(u3_uc_bold$Kingdom == "Fungi" | u3_uc_bold$Class == "Mammalia", "non-diet", 
         ifelse(u3_uc_bold$Phylum == "Arthropoda" | u3_uc_bold$Class == "Reptilia" & u3_uc_bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
u3_uc_ncbi$taxonomy <- u3_uc_ncbi$Category
  
#  ifelse(
#  d2_c_ncbi$ID == "Sparassidae", "predator",
#  ifelse(d2_c_ncbi$Phylum == "Arthropoda" & d2_c_ncbi$ID != "Heteropoda venatoria", 
#         "prey", "non-diet"
#  ))

#thus, the taxonomy variable has 5 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table), and unknown
#subset only the variables of interest from these two dataframes before merging
u3_uc_bold_id <- u3_uc_bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

u3_uc_ncbi_id <- u3_uc_ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
u3_uc_id_all <- u3_uc_bold_id %>%  
  full_join(u3_uc_ncbi_id, by = c("ASV", "taxonomy"))

u3_uc_id_all <- u3_uc_id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y)

#this binds them all to the community matrix for analyses
u3_uc_id_tab <- u3_uc %>%
  left_join(u3_uc_id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
u3_uc_id_tab$taxonomy <- replace_na(u3_uc_id_tab$taxonomy, "no hit")

#write the taxonomies to a file for use in the community dissimilarity scripts script
u3_uc_taxonomies <- u3_uc_id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, taxonomy)
write.csv(u3_uc_taxonomies, here("data", "u3_uc_tax_ass.csv"))

#UNOISE CLEAN##
#import both NCBI and BOLD taxonomies
u3_c_ncbi <- read.csv(here("6_taxonomic_assignment", "MEGAN", "unoise_clean",
                            "prey", "unoise_c_prey_ALL_ncbi_taxa.csv"))
u3_c_bold <- read.csv(here("6_taxonomic_assignment", "BOLD", "usearch_c", 
                            "prey", "unoise3_c_prey_bold_taxa.csv"))
u3_c_bold <- rename(u3_c_bold, "ASV" = "Query.ID")

#Sort by "predator", "prey", and "non-diet" before merging into one 
#dataframe for assigning to the abundance table!
#this following ifelse set will depend on the taxonomy of the set of samples
#and for this one, these are the categories that seemed to be appropriate based
#on the taxonomy list after looking at it
#sort into categories for bold:
u3_c_bold$taxonomy <- ifelse(
  u3_c_bold$ID == "Heteropoda venatoria", "predator",
  ifelse(u3_c_bold$Kingdom == "Fungi" | u3_c_bold$Class == "Mammalia", "non-diet", 
         ifelse(u3_c_bold$Phylum == "Arthropoda" & u3_c_bold$ID != "Heteropoda venatoria", 
                "prey", "non-diet"
         )))

#sort into categories for ncbi:
u3_c_ncbi$taxonomy <- u3_c_ncbi$Category

#  ifelse(
#  d2_c_ncbi$ID == "Sparassidae", "predator",
#  ifelse(d2_c_ncbi$Phylum == "Arthropoda" & d2_c_ncbi$ID != "Heteropoda venatoria", 
#         "prey", "non-diet"
#  ))

#thus, the taxonomy variable has 5 levels: predator, prey, non-diet, non-hit (which
#we will assign below after merging with full ASV table), and unknown
#subset only the variables of interest from these two dataframes before merging
u3_c_bold_id <- u3_c_bold %>%
  dplyr::select(ASV, ID, Order, taxonomy)

u3_c_ncbi_id <- u3_c_ncbi %>%
  dplyr::select(ASV, ID, Order, taxonomy)

#join them together by ASV and taxonomy
u3_c_id_all <- u3_c_bold_id %>%  
  full_join(u3_c_ncbi_id, by = c("ASV", "taxonomy"))

u3_c_id_all <- u3_c_id_all %>%
  rename(ID_bold = ID.x,
         Order_bold = Order.x,
         ID_ncbi = ID.y,
         Order_ncbi = Order.y)

#this binds them all to the community matrix for analyses
u3_c_id_tab <- u3_c %>%
  left_join(u3_c_id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
u3_c_id_tab$taxonomy <- replace_na(u3_c_id_tab$taxonomy, "no hit")

#write the taxonomies to a file for use in the community dissimilarity scripts script
u3_c_taxonomies <- u3_c_id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, taxonomy)
write.csv(u3_c_taxonomies, here("data", "u3_c_tax_ass.csv"))
#3. total number of prey ASVs for each pipeline ####
#summarize each dataframe by sample (which also gives ASVs in the ASV category)
#grouped by taxonomy, giving both a total ASV count in each category AND
#a per sample value for each category (no hit, non-diet, predator, prey, unknown)
d2_uc_ASVs <- d2_uc_id_tab %>%
  group_by(taxonomy) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, ASVs_d2_uc, ASV:QC1, factor_key=TRUE) %>%
  dplyr::select(taxonomy, sample, ASVs_d2_uc)
  
d2_c_ASVs <- d2_c_id_tab %>%
  group_by(taxonomy) %>%
  summarize_all(~sum(. !=0)) %>%
  gather(sample, ASVs_d2_c, ASV:HEV89, factor_key = TRUE) %>%
  dplyr::select(taxonomy, sample, ASVs_d2_c)

u3_uc_ASVs <- u3_uc_id_tab %>%
  group_by(taxonomy) %>%
  summarize_all(~sum(. !=0)) %>%
  gather(sample, ASVs_u3_uc, ASV:QC1, factor_key = TRUE) %>%
  dplyr::select(taxonomy, sample, ASVs_u3_uc)

u3_c_ASVs <- u3_c_id_tab %>%
  group_by(taxonomy) %>%
  summarize_all(~sum(. !=0)) %>%
  gather(sample, ASVs_u3_c, ASV:HEV99, factor_key = TRUE) %>%
  dplyr::select(taxonomy, sample, ASVs_u3_c)

#join these four together
taxa_ASVs <- d2_uc_ASVs %>%
  full_join(d2_c_ASVs, by = c("taxonomy", "sample")) %>%
  full_join(u3_uc_ASVs, by = c("taxonomy", "sample")) %>%
  full_join(u3_c_ASVs, by = c("taxonomy", "sample"))

taxa_ASVs[is.na(taxa_ASVs)] <- 0

#filter out controls for analyses of by sample reads
taxa_ASVs_samples <- taxa_ASVs %>%
  filter(!sample %in% c("ASV", "CL1", "NEG", "CL4", "QC1"))

#make this dataframe longer for analysis with lme4
taxa_ASVs_samples_long <- taxa_ASVs_samples %>%
  gather(pipeline, ASVs, ASVs_d2_uc:ASVs_u3_c)

#grouped lme by sample to see if there are differences across pipelines

#select only prey ASVs
taxa_ASVs_prey <- taxa_ASVs_samples_long %>%
  filter(taxonomy == "prey")

#spllit by field and lab by joining with metadata first
taxa_ASVs_prey <- taxa_ASVs_prey %>%
  left_join(metadata, by = "sample")

#create binary sterilized factor
taxa_ASVs_prey$Ster <- ifelse(taxa_ASVs_prey$Sterilized == "NS", 0, 1)
taxa_ASVs_prey$Ster <- as.factor(as.character(taxa_ASVs_prey$Ster))

#split field and lab
f_taxa_ASVs_prey <- taxa_ASVs_prey %>%
  filter(Source == "FIELD")

l_taxa_ASVs_prey <- taxa_ASVs_prey %>%
  filter(Source == "LAB")

#then create a model fitting the number of ASVs to pipeline grouped by sample
lme_prey_ASVs_f <- glmmTMB(ASVs ~ pipeline*Ster + (1|sample),
                 data = f_taxa_ASVs_prey,
                 family = "genpois",
                 REML = FALSE)
lme_prey_ASVs_f2 <- glmmTMB(ASVs ~ pipeline + Ster + (1|sample),
                           data = f_taxa_ASVs_prey,
                           family = "genpois",
                           REML = FALSE)
lme_prey_ASVs_f3 <- glmmTMB(ASVs ~ pipeline + (1|sample),
                            data = f_taxa_ASVs_prey,
                            family = "genpois", 
                            REML = FALSE)
lme_prey_ASVs_f4 <- glmmTMB(ASVs ~ Ster + (1|sample),
                            data = f_taxa_ASVs_prey,
                            family = "genpois",
                            REML = FALSE)
lme_prey_ASVs_null_f <- glmmTMB(ASVs ~ 1 + (1|sample),
                      data = f_taxa_ASVs_prey,
                      family = "genpois",
                      REML = FALSE)

#compare by AIC values
AICc(lme_prey_ASVs_f, lme_prey_ASVs_f2, lme_prey_ASVs_f3, lme_prey_ASVs_f4,
     lme_prey_ASVs_null_f)
#this shows that 2 and 3 are very similar to each other, with 3 being lower
#by one AIC value AND having a more parsimonious structure, we'll check 
#sterilization significance anyway below:
model.emm_preyASV_f <- emmeans(lme_prey_ASVs_f2, "pipeline")
model.emm_preyASV_f_s <- emmeans(lme_prey_ASVs_f2, "Ster")
pairs(model.emm_preyASV_f)
pairs(model.emm_preyASV_f_s)
#Sterilization is not significant, so we will stick with more parsimonious
#model for the rest:
lme_prey_ASVs_f3 <- glmmTMB(ASVs ~ pipeline + (1|sample),
                            data = f_taxa_ASVs_prey,
                            family = "genpois")
model.emm_preyASV_f3 <- emmeans(lme_prey_ASVs_f3, "pipeline")
pairs(model.emm_preyASV_f3)
#UNOISE3 had a greater number of prey ASVs for both clean and unclean
#cleaning did not significantly increase the number of prey ASVs for 
#either pipeline. Sterilization did not significantly change the number
#of prey ASVs in any pipeline

x_preyASV_f <- residuals(lme_prey_ASVs_f3)
plotNormalHistogram(x_preyASV_f)
#residuals look relatively normal
summary(lme_prey_ASVs_f3)
plot(residuals(lme_prey_ASVs_f3)) #residuals look fairly homoskedastic
plot(allEffects(lme_prey_ASVs_f3)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyASV_f <- simulateResiduals(fittedModel = lme_prey_ASVs_f3) 
preyASVfit_f <- plot(simulationOutput_preyASV_f, asFactor=TRUE) #these look good
preyASVzi_f <- testZeroInflation(simulationOutput_preyASV_f) #not zero inflated
preyASVod_f <- testDispersion(simulationOutput_preyASV_f) #not overdispersed

#boxplot visualization
prey_ASV_graph_f <- ggplot(f_taxa_ASVs_prey, aes(x = pipeline, y = ASVs, fill = Sterilized)) +
         geom_boxplot() +
         labs(x = "Pipeline", y = "Prey ASVs per sample", 
              title = "Field: Total prey ASVs assigned to each sample") +
         theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#LAB prey ASVs
#then create a model fitting the number of ASVs to pipeline grouped by sample
lme_prey_ASVs_l <- glmmTMB(ASVs ~ pipeline*Ster + (1|sample),
                           data = l_taxa_ASVs_prey,
                           family = "genpois",
                           REML = FALSE)
lme_prey_ASVs_l2 <- glmmTMB(ASVs ~ pipeline + Ster + (1|sample),
                            data = l_taxa_ASVs_prey,
                            family = "genpois",
                            REML =FALSE)
lme_prey_ASVs_l3 <- glmmTMB(ASVs ~ pipeline + (1|sample),
                            data = l_taxa_ASVs_prey,
                            family = "genpois",
                            REML = FALSE)
lme_prey_ASVs_l4 <- glmmTMB(ASVs ~ Ster + (1|sample),
                            data = l_taxa_ASVs_prey,
                            family = "genpois",
                            REML = FALSE)
lme_prey_ASVs_null_l <- glmmTMB(ASVs ~ 1 + (1|sample),
                                data = l_taxa_ASVs_prey,
                                family = "genpois",
                                REML = FALSE)

#compare by AIC values
AICc(lme_prey_ASVs_l, lme_prey_ASVs_l2, lme_prey_ASVs_l3, lme_prey_ASVs_l4,
     lme_prey_ASVs_null_l)
#this shows that 4 and null are most closely related to each other, with
#the null being the most parsimonious

#furthermore, when we do a pairwise comparision of sterilization from
#model 4, we see that the pairwise comparision is non-significant
lme_prey_ASVs_l4 <- glmmTMB(ASVs ~ Ster + (1|sample),
                            data = l_taxa_ASVs_prey,
                            family = "genpois")
model.emm_preyASV_l <- emmeans(lme_prey_ASVs_f4, "Ster")
pairs(model.emm_preyASV_l)

#so best model is null:
lme_prey_ASVs_null_l <- glmmTMB(ASVs ~ 1 + (1|sample),
                                data = l_taxa_ASVs_prey,
                                family = "genpois")
x_preyASV_l <- residuals(lme_prey_ASVs_null_l)
plotNormalHistogram(x_preyASV_l)
#residuals look relatively normal
summary(lme_prey_ASVs_null_l)
plot(residuals(lme_prey_ASVs_null_l)) #residuals look fairly homoskedastic

#checking other model assumptions here:
simulationOutput_preyASV_l <- simulateResiduals(fittedModel = lme_prey_ASVs_null_l) 
preyASVfit_l <- plot(simulationOutput_preyASV_l, asFactor=TRUE) #these look good
preyASVzi_l <- testZeroInflation(simulationOutput_preyASV_l) #not zero inflated
preyASVod_l <- testDispersion(simulationOutput_preyASV_l) #not overdispersed

#boxplot visualization
prey_ASV_graph_l <- ggplot(l_taxa_ASVs_prey, aes(x = pipeline, y = ASVs, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey ASVs per sample", 
       title = "Lab: Total prey ASVs assigned to each sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#4. total prey ASV read abundance ####
#add up the total number of reads by group in each sample for each pipeline
d2_uc_reads <- d2_uc_id_tab %>%
  dplyr::select(-ID_bold, -Order_bold, -ID_ncbi, -Order_ncbi, -ID_nondiet, -Order_nondiet, -Category) %>%
  group_by(ASV) %>%
  gather(sample, reads, CL1:QC1, factor_key = TRUE) %>%
  group_by(sample, taxonomy) %>%
  summarize(reads_d2_uc = sum(reads))

d2_c_reads <- d2_c_id_tab %>%
  dplyr::select(-ID_bold, -Order_bold, -ID_ncbi, -Order_ncbi, -ID_nondiet, -Category) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV100:HEV89, factor_key=TRUE) %>%
  group_by(sample, taxonomy) %>%
  summarize(reads_d2_c = sum(reads))

u3_uc_reads <- u3_uc_id_tab %>%
  dplyr::select(-ID_bold, -Order_bold, -ID_ncbi, -Order_ncbi) %>%
  group_by(ASV) %>%
  gather(sample, reads, CL1:QC1, factor_key = TRUE) %>%
  group_by(sample, taxonomy) %>%
  summarize(reads_u3_uc = sum(reads))

u3_c_reads <- u3_c_id_tab %>%
  dplyr::select(-ID_bold, -Order_bold, -ID_ncbi, -Order_ncbi) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV99, factor_key=TRUE) %>%
  group_by(sample, taxonomy) %>%
  summarize(reads_u3_c = sum(reads))
  
#join these all together now
taxa_reads <- d2_uc_reads %>%
  full_join(d2_c_reads, by = c("taxonomy", "sample")) %>%
  full_join(u3_uc_reads, by = c("taxonomy", "sample")) %>%
  full_join(u3_c_reads, by = c("taxonomy", "sample"))

taxa_reads[is.na(taxa_reads)] <- 0

#filter out controls for analyses of by sample reads
taxa_reads_samples <- taxa_reads %>%
  filter(!sample %in% c("CL1", "NEG", "CL4", "QC1"))

#make this dataframe longer for analysis with lme4
taxa_reads_samples_long <- taxa_reads_samples %>%
  gather(pipeline, reads, reads_d2_uc:reads_u3_c)

#select only prey reads
taxa_reads_prey <- taxa_reads_samples_long %>%
  filter(taxonomy == "prey")

taxa_reads_prey <- taxa_reads_prey %>%
  left_join(metadata, by = "sample")

taxa_reads_prey$Ster <- as.factor(ifelse(taxa_reads_prey$Sterilized == "NS", 0, 1))

taxa_reads_prey_f <- taxa_reads_prey %>%
  filter(Source == "FIELD")

taxa_reads_prey_l <- taxa_reads_prey %>%
  filter(Source == "LAB")

#field prey reads
#then create a model fitting the number of reads to pipeline grouped by sample
read_mod_f <- glmmTMB(reads ~ pipeline*Ster + (1|sample),
                 data = taxa_reads_prey_f,
                 family = "genpois",
                 REML = FALSE)
read_mod_f2 <- glmmTMB(reads ~ pipeline + Ster + (1|sample),
                      data = taxa_reads_prey_f,
                      family = "genpois",
                      REML = FALSE)
read_mod_f3 <- glmmTMB(reads ~ pipeline + (1|sample),
                      data = taxa_reads_prey_f,
                      family = "genpois",
                      REML = FALSE)
read_mod_f4 <- glmmTMB(reads ~ Ster + (1|sample),
                      data = taxa_reads_prey_f,
                      family = "genpois",
                      REML = FALSE)
read_mod_null_f <- glmmTMB(reads ~ 1 + (1|sample),
                      data = taxa_reads_prey_f,
                      family = "genpois",
                      REML = FALSE)

#compare models with AIC
AICc(read_mod_f, read_mod_f2, read_mod_f3, read_mod_f4, read_mod_null_f)
#models 2 and 3 are lowest, with 3 being the most parsimonious, which includes 
#only pipeline. 
read_mod_f2 <- glmmTMB(reads ~ pipeline + Ster + (1|sample),
                       data = taxa_reads_prey_f,
                       family = "genpois")
model.emm_preyreads_fs <- emmeans(read_mod_f2, "Ster")
pairs(model.emm_preyreads_fs)
#since pairs from sterilization non-significant, goign to use the more parsimonious
#model

read_mod_f3 <- glmmTMB(reads ~ pipeline + (1|sample),
                       data = taxa_reads_prey_f,
                       family = "genpois")
model.emm_preyreads_f <- emmeans(read_mod_f3, "pipeline")
pairs(model.emm_preyreads_f)
#cleaning did not change the number of prey reads for either DADA2 or UNOISE3.
#Unoise had significantly greater reads than DADA2

x_read_f <- residuals(read_mod_f3)
plotNormalHistogram(x_read_f)
#residuals look relatively normal
summary(read_mod_f3)
plot(residuals(read_mod_f3)) #residuals look fairly homoskedastic
plot(allEffects(read_mod_f3)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyreads_f <- simulateResiduals(fittedModel = read_mod_f3) 
preyreadsfit_f <- plot(simulationOutput_preyreads_f, asFactor=TRUE) #these look good
preyreadszi_f <- testZeroInflation(simulationOutput_preyreads_f) #looks good
preyreadsod_f <- testDispersion(simulationOutput_preyreads_f) #not overdispersed

#boxplot visualization
prey_reads_graph_f <- ggplot(taxa_reads_prey_f, aes(x = pipeline, y = reads, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey reads per sample", 
       title = "Field: Total number of prey reads per ASV per sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#LAB prey reads
#then create a model fitting the number of reads to pipeline grouped by sample
read_mod_l <- glmmTMB(reads ~ pipeline*Ster + (1|sample),
                      data = taxa_reads_prey_l,
                      family = "genpois",
                      REML = FALSE)
read_mod_l2 <- glmmTMB(reads ~ pipeline + Ster + (1|sample),
                       data = taxa_reads_prey_l,
                       family = "genpois",
                       REML = FALSE)
read_mod_l3 <- glmmTMB(reads ~ pipeline + (1|sample),
                       data = taxa_reads_prey_l,
                       family = "genpois",
                       REML = FALSE)
read_mod_l4 <- glmmTMB(reads ~ Ster + (1|sample),
                       data = taxa_reads_prey_l,
                       family = "genpois",
                       REML = FALSE)
read_mod_null_l <- glmmTMB(reads ~ 1 + (1|sample),
                           data = taxa_reads_prey_l,
                           family = "genpois",
                           REML = FALSE)

#compare models with AIC
AICc(read_mod_l, read_mod_l2, read_mod_l3, read_mod_l4, read_mod_null_l)
#model 3 is more than two AIC values better fit than any other model, so
#we will continue with this model, which only includes pipeline.
read_mod_l3 <- glmmTMB(reads ~ pipeline + (1|sample),
                       data = taxa_reads_prey_l,
                       family = "genpois")
model.emm_preyreads_l <- emmeans(read_mod_l3, "pipeline")
pairs(model.emm_preyreads_l) 
#UNOISE3 has a greater number of prey reads than DADA2 for both
#cleaned and uncleaned datasets. Cleaning did not significantly
#change the number of prey reads per sample. 


x_read_l <- residuals(read_mod_l3)
plotNormalHistogram(x_read_l)
#residuals look relatively normal
summary(read_mod_l3)
plot(residuals(read_mod_l3)) #residuals look fairly homoskedastic
plot(allEffects(read_mod_l3)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyreads_l <- simulateResiduals(fittedModel = read_mod_l3) 
preyreadsfit_l <- plot(simulationOutput_preyreads_l, asFactor=TRUE) #these look good
preyreadszi_l <- testZeroInflation(simulationOutput_preyreads_l) #looks good
preyreadsod_l <- testDispersion(simulationOutput_preyreads_l) #not overdispersed

#boxplot visualization
prey_reads_graph_l <- ggplot(taxa_reads_prey_l, aes(x = pipeline, y = reads, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey reads per sample", title = "Lab: Total number of prey reads per ASV per sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))
#5. proportion of prey ASVs to total ASVs (number of ASVs) ####
#I'm using a general poisson with an offset to model proportions.
#the tot_oth column created here along with this dataframe can be 
#used to look at prey/other ratio in Step 7
#will be using d2_uc_ASVs, d2_c_ASVs, u3_uc_ASVs, and u3_c_ASVs here
d2_uc_totals <- d2_uc_ASVs %>%
  spread(taxonomy, ASVs_d2_uc) %>%
  mutate(tot_oth = (predator + `no hit` + `non-diet`),
         tot = (predator + `no hit` + `non-diet` + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
d2_uc_totals$pipeline <- "d2_uc"

d2_c_totals <- d2_c_ASVs %>%
  spread(taxonomy, ASVs_d2_c) %>%
  mutate(tot_oth = (predator + `non-diet`),
         tot = (predator + `non-diet` + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
d2_c_totals$pipeline <- "d2_c"

u3_uc_totals <- u3_uc_ASVs %>%
  spread(taxonomy, ASVs_u3_uc) %>%
  mutate(tot_oth = (predator + `no hit` + `non-diet`),
         tot = (predator + `no hit` + `non-diet` + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
u3_uc_totals$pipeline <- "u3_uc"

u3_c_totals <- u3_c_ASVs %>%
  spread(taxonomy, ASVs_u3_c) %>%
  mutate(tot_oth = (predator + `non-diet`),
         tot = (predator + `non-diet` + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
u3_c_totals$pipeline <- "u3_c"

#join these all together now
ASV_totals <- bind_rows(d2_uc_totals, d2_c_totals,
                        u3_uc_totals, u3_c_totals)

#filter out controls for analyses of by sample ASV
ASV_totals_samples <- ASV_totals %>%
  filter(!sample %in% c("ASV", "CL1", "NEG", "CL4", "QC1"))

#total ASVS
ASV_totals %>%
  filter(sample == "ASV") 

#ASV_totals_samples$tot_scaled <- scale(ASV_totals_samples$tot)[, 1]
ASV_totals_samples <- ASV_totals_samples[which(ASV_totals_samples$tot > 0),]

ASV_totals_samples <- ASV_totals_samples %>%
  left_join(metadata, by = "sample")

ASV_totals_samples$Ster <- as.factor(ifelse(ASV_totals_samples$Sterilized == "NS", 0, 1))

ASV_totals_samples_f <- ASV_totals_samples %>%
  filter(Source == "FIELD")

ASV_totals_samples_l <- ASV_totals_samples %>%
  filter(Source == "LAB")
#FIELD
#then create a model fitting the number of ASVs to pipeline grouped by sample
prey_prop_mod_f <- glmmTMB(prey ~ pipeline*Ster + (1|sample),
                          data = ASV_totals_samples_f,
                          family = "genpois",
                         offset = log(tot),
                         REML = FALSE)
prey_prop_mod_f2 <- glmmTMB(prey ~ pipeline + Ster + (1|sample),
                           data = ASV_totals_samples_f,
                           family = "genpois",
                           offset = log(tot),
                           REML = FALSE)
prey_prop_mod_f3 <- glmmTMB(prey ~ pipeline + (1|sample),
                           data = ASV_totals_samples_f,
                           family = "genpois",
                           offset = log(tot),
                           REML = FALSE)
prey_prop_mod_f4 <- glmmTMB(prey ~ Ster + (1|sample),
                           data = ASV_totals_samples_f,
                           family = "genpois",
                           offset = log(tot),
                           REML = FALSE)
#the null model where pipeline does not matter
prey_prop_mod_null_f <- glmmTMB(prey ~ 1 + (1|sample),
                              data = ASV_totals_samples_f,
                              family = "genpois",
                              offset = log(tot),
                              REML= FALSE)
#compare these with AIC
AICc(prey_prop_mod_f, prey_prop_mod_f2, prey_prop_mod_f3, prey_prop_mod_f4,
     prey_prop_mod_null_f)
#model 3 is *just* over two AIC values lower than model 2, has the both the lowest
#AIC and the most parsimonious structure
prey_prop_mod_f3 <- glmmTMB(prey ~ pipeline + (1|sample),
                            data = ASV_totals_samples_f,
                            family = "genpois",
                            offset = log(tot))

model.emm_preyprop_f <- emmeans(prey_prop_mod_f3, "pipeline")
pairs(model.emm_preyprop_f)
#Cleaning increased proportion of prey ASVs for both DADA2 and UNOISE3
#The uncleaned datasets do not have significantly different numbers of
#prey ASVs, though DADA2 clean has more than UNOISE clean (probs because
#of using DADA2 for BBSPLIT, TBH). 

x_prop_f <- residuals(prey_prop_mod_f3)
plotNormalHistogram(x_prop_f)
#residuals look relatively normal
plot(residuals(prey_prop_mod_f3)) #residuals look fairly homoskedastic
plot(allEffects(prey_prop_mod_f3)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyprop_f <- simulateResiduals(fittedModel = prey_prop_mod_f3) 
preypropfit <- plot(simulationOutput_preyprop_f, asFactor=TRUE) #these look weird now... with offset
preypropzi <- testZeroInflation(simulationOutput_preyprop_f) #what is the opposite of zero inflation??
preypropod <- testDispersion(simulationOutput_preyprop_f) #not overdispersed

#boxplot visualization
prey_prop_graph_f <- ggplot(ASV_totals_samples_f, aes(x = pipeline, y = prey/tot, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey ASV percentage by sample", title = "Field: Percent of prey ASVs per sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#LAB
#then create a model fitting the number of ASVs to pipeline grouped by sample
prey_prop_mod_l <- glmmTMB(prey ~ pipeline*Ster + (1|sample),
                           data = ASV_totals_samples_l,
                           family = "genpois",
                           offset = log(tot),
                           REML = FALSE)
prey_prop_mod_l2 <- glmmTMB(prey ~ pipeline + Ster + (1|sample),
                            data = ASV_totals_samples_l,
                            family = "genpois",
                            offset = log(tot),
                            REML = FALSE)
prey_prop_mod_l3 <- glmmTMB(prey ~ pipeline + (1|sample),
                            data = ASV_totals_samples_l,
                            family = "genpois",
                            offset = log(tot),
                            REML = FALSE)
prey_prop_mod_l4 <- glmmTMB(prey ~ Ster + (1|sample),
                            data = ASV_totals_samples_l,
                            family = "genpois",
                            offset = log(tot),
                            REML = FALSE)
#the null model where pipeline does not matter
prey_prop_mod_null_l <- glmmTMB(prey ~ 1 + (1|sample),
                                data = ASV_totals_samples_l,
                                family = "genpois",
                                offset = log(tot),
                                REML= FALSE)
#compare these with AIC
AICc(prey_prop_mod_l, prey_prop_mod_l2, prey_prop_mod_l3, prey_prop_mod_l4,
     prey_prop_mod_null_l)
#model 1 has the lowest AIC, which means that pipeline is important and the 
#effect varies with sterilization
prey_prop_mod_l <- glmmTMB(prey ~ pipeline*Ster + (1|sample),
                            data = ASV_totals_samples_f,
                            family = "genpois",
                            offset = log(tot))

model.emm_preyprop_l <- emmeans(prey_prop_mod_l, pairwise ~ pipeline | Ster)
model.emm_preyprop_ls <- emmeans(prey_prop_mod_l, pairwise ~ Ster | pipeline)
pairs(model.emm_preyprop_l)
pairs(model.emm_preyprop_ls)
#Cleaning increased proportion of prey ASVs for both DADA2 and UNOISE3
#The uncleaned datasets do not have significantly different numbers of
#prey ASVs
#for unsterilized individuals, DADA2 produced a higher proportion of prey
#ASVs than UNOISE, though for sterilized individuals, these were not
#significantly different. 

x_prop_l<- residuals(prey_prop_mod_l)
plotNormalHistogram(x_prop_l)
#residuals look relatively normal
plot(residuals(prey_prop_mod_l)) #residuals look fairly homoskedastic
plot(allEffects(prey_prop_mod_l)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyprop_l <- simulateResiduals(fittedModel = prey_prop_mod_l) 
preypropfit_l <- plot(simulationOutput_preyprop_l, asFactor=TRUE) #these look weird now... with offset
preypropzi_l <- testZeroInflation(simulationOutput_preyprop_l) #what is the opposite of zero inflation??
preypropod_l <- testDispersion(simulationOutput_preyprop_l) #not overdispersed

#boxplot visualization
prey_prop_graph_l <- ggplot(ASV_totals_samples_l, aes(x = pipeline, y = prey/tot, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey ASV percentage by sample", title = "Lab: Percent of prey ASVs per sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#6. proportion of prey reads to total ASV reads (read abundance) ####
#will be using d2_uc_reads, d2_c_reads, u3_uc_reads, and u3_c_reads here
#as in the ASV steps above, i'm creating both a total reads and a total
#non-prey reads cateogry here that we can use in step 8 below
d2_uc_read_tot <- d2_uc_reads %>%
  spread(taxonomy, reads_d2_uc) %>%
  mutate(tot_oth = (predator + `no hit` + `non-diet`),
         tot = (predator + `no hit` + `non-diet` + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
d2_uc_read_tot$pipeline <- "d2_uc"

d2_c_read_tot <- d2_c_reads %>%
  spread(taxonomy, reads_d2_c) %>%
  mutate(tot_oth = (predator + `non-diet` ),
         tot = (predator + `non-diet`  + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
d2_c_read_tot$pipeline <- "d2_c"

u3_uc_read_tot <- u3_uc_reads %>%
  spread(taxonomy, reads_u3_uc) %>%
  mutate(tot_oth = (predator + `no hit` + `non-diet` ),
         tot = (predator + `no hit` + `non-diet`  + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
u3_uc_read_tot$pipeline <- "u3_uc"

u3_c_read_tot <- u3_c_reads %>%
  spread(taxonomy, reads_u3_c) %>%
  mutate(tot_oth = (predator + `non-diet` ),
         tot = (predator + `non-diet`  + prey)) %>%
  dplyr::select(sample, prey, tot_oth, tot)
u3_c_read_tot$pipeline <- "u3_c"


#join these all together now
read_totals <- bind_rows(d2_uc_read_tot, d2_c_read_tot,
                        u3_uc_read_tot, u3_c_read_tot)

#filter out controls for analyses of by sample ASV
read_totals_samples <- read_totals %>%
  filter(!sample %in% c("CL1", "CL4", "QC1", "NEG"))

read_totals_samples <- read_totals_samples %>%
  left_join(metadata, by = "sample")

read_totals_samples$Ster <- as.factor(ifelse(read_totals_samples$Sterilized == "NS", 0, 1))

read_totals_samples_f <- read_totals_samples %>%
  filter(Source == "FIELD")

read_totals_samples_l <- read_totals_samples %>%
  filter(Source == "LAB")

#there are three samples with ZERO reads in the total column, and when 
#i remove them the model below runs. 
read_totals_samples_f <- read_totals_samples_f[which(read_totals_samples_f$tot >0),]
#then create a model fitting the number of ASVs to pipeline grouped by sample
prey_rprop_mod_f <- glmmTMB(prey ~ pipeline*Ster + (1|sample),
                         data = read_totals_samples_f,
                         family = "genpois",
                         offset = log(tot),
                         REML =FALSE)
prey_rprop_mod_f2 <- glmmTMB(prey ~ pipeline + Ster + (1|sample),
                            data = read_totals_samples_f,
                            family = "genpois",
                            offset = log(tot),
                            REML =FALSE)
prey_rprop_mod_f3 <- glmmTMB(prey ~ pipeline + (1|sample),
                            data = read_totals_samples_f,
                            family = "genpois",
                            offset = log(tot),
                            REML =FALSE)
prey_rprop_mod_f3zi <- glmmTMB(prey ~ pipeline + (1|sample),
                             data = read_totals_samples_f,
                             family = "genpois",
                             offset = log(tot),
                             ziformula = ~1,
                             REML =FALSE)
prey_rprop_mod_f4 <- glmmTMB(prey ~ Ster + (1|sample),
                            data = read_totals_samples_f,
                            family = "genpois",
                            offset = log(tot),
                            REML =FALSE)
prey_rprop_mod_null_f <- glmmTMB(prey ~ 1 + (1|sample),
                              data = read_totals_samples_f,
                              family = "genpois",
                              offset = log(tot),
                              REML = FALSE)

AICc(prey_rprop_mod_f, prey_rprop_mod_f2, prey_rprop_mod_f3,prey_rprop_mod_f3zi, prey_rprop_mod_f4,
     prey_rprop_mod_null_f)

#based on these AIC values, model 3 is the best, though 2 is within 2 AIC values. 
prey_rprop_mod_f3zi <- glmmTMB(prey ~ pipeline + (1|sample),
                             data = read_totals_samples_f,
                             family = "genpois",
                             ziformula = ~1,
                             offset = log(tot))

model.emm_preyrprop_f <- emmeans(prey_rprop_mod_f3zi, "pipeline")
pairs(model.emm_preyrprop_f)
#cleaning significantly increased the proportion of prey reads of total
#UNOISE unclean has significantly higher proportions of prey reads than
#does DADA2

x_rprop_f <- residuals(prey_rprop_mod_f3zi)
plotNormalHistogram(x_rprop_f)
#residuals look relatively normal
plot(residuals(prey_rprop_mod_f3zi)) #residuals look fairly homoskedastic
plot(allEffects(prey_rprop_mod_f3zi)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyrprop_f <- simulateResiduals(fittedModel = prey_rprop_mod_f3zi) 
preyrpropfit_f <- plot(simulationOutput_preyrprop_f, asFactor=TRUE) #haha what is with the residuals?
preyrpropzi_f <- testZeroInflation(simulationOutput_preyrprop_f) #look fine
preyrpropod_f <- testDispersion(simulationOutput_preyrprop_f) #not overdispersed

#boxplot visualization
prey_rprop_graph_f <- ggplot(read_totals_samples_f, aes(x = pipeline, y = prey/tot, 
                                                        fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey read proportion by sample", 
       title = "Field: Perecent of prey reads in each sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#LAB
#then create a model fitting the number of ASVs to pipeline grouped by sample
prey_rprop_mod_l <- glmmTMB(prey ~ pipeline*Ster + (1|sample),
                            data = read_totals_samples_l,
                            family = "genpois",
                            offset = log(tot),
                            REML =FALSE)
prey_rprop_mod_l2 <- glmmTMB(prey ~ pipeline + Ster + (1|sample),
                             data = read_totals_samples_l,
                             family = "genpois",
                             offset = log(tot),
                             REML =FALSE)
prey_rprop_mod_l3 <- glmmTMB(prey ~ pipeline + (1|sample),
                             data = read_totals_samples_l,
                             family = "genpois",
                             offset = log(tot),
                             REML =FALSE)
prey_rprop_mod_l4 <- glmmTMB(prey ~ Ster + (1|sample),
                             data = read_totals_samples_l,
                             family = "genpois",
                             offset = log(tot),
                             REML =FALSE)
prey_rprop_mod_null_l <- glmmTMB(prey ~ 1 + (1|sample),
                                 data = read_totals_samples_l,
                                 family = "genpois",
                                 offset = log(tot),
                                 REML = FALSE)

AICc(prey_rprop_mod_l, prey_rprop_mod_l2, prey_rprop_mod_l3, prey_rprop_mod_l4,
     prey_rprop_mod_null_l)

#based on these AIC values, model 3 is the best 
prey_rprop_mod_l3 <- glmmTMB(prey ~ pipeline + (1|sample),
                               data = read_totals_samples_l,
                               family = "genpois",
                               offset = log(tot))

model.emm_preyrprop_l <- emmeans(prey_rprop_mod_l3, "pipeline")
pairs(model.emm_preyrprop_l)
#cleaning significantly increased the proportion of prey reads of total
#UNOISE and DADA2 do not have significantly different number of prey ASVs
#for either clean or unclean datasets

x_rprop_l <- residuals(prey_rprop_mod_l3)
plotNormalHistogram(x_rprop_l)
#residuals look relatively normal
plot(residuals(prey_rprop_mod_l3)) #residuals look fairly homoskedastic
plot(allEffects(prey_rprop_mod_l3)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_preyrprop_l <- simulateResiduals(fittedModel = prey_rprop_mod_l3) 
preyrpropfit_l <- plot(simulationOutput_preyrprop_l, asFactor=TRUE) #haha what is with the residuals?
preyrpropzi_l <- testZeroInflation(simulationOutput_preyrprop_l) #look fine
preyrpropod_l <- testDispersion(simulationOutput_preyrprop_l) #not overdispersed

#boxplot visualization
prey_rprop_graph_l <- ggplot(read_totals_samples_l, aes(x = pipeline, y = prey/tot, 
                                                        fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Prey read proportion by sample", 
       title = "Lab: Perecent of prey reads in each sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))
###a. consider here sequencing depth? Ask Austen
#7. Prey/OTHER ASV ratios ####
#the tot_oth column created in ASV_totals_samples 
#can be used in this analysis
#need to delete those for which total other is equal to zero
#ASV_totals_samples_ratio <- ASV_totals_samples[which(ASV_totals_samples$tot_oth > 0),]
#then create a model fitting the number of ASVs to pipeline grouped by sample
#prey_ratio_mod <- glmmTMB(prey ~ pipeline + (1|sample),
#                         data = ASV_totals_samples_ratio,
#                         family = "genpois",
#                         ziformula = ~1,
#                         offset = log(tot_oth))
#the null model where pipeline does not matter
#prey_ratio_mod_null <- glmmTMB(prey ~ 1 + (1|sample),
#                              data = ASV_totals_samples_ratio,
#                              family = "genpois",
#                              ziformula = ~1,
#                              offset = log(tot_oth))

#model.emm <- emmeans(prey_ratio_mod, "pipeline")
#pairs(model.emm)
#cleaning increased the proportion of prey ASVs for
#both dada2 and unoise
#dada2 uncleaned had higher prey ASV ratio than unoise uncleaned
#the cleaned datasets do not vary in their prey/other ASV ratio

#x_rat <- residuals(prey_ratio_mod)
#plotNormalHistogram(x_rat)
#residuals look relatively normal
#plot(residuals(prey_ratio_mod)) #residuals look fairly homoskedastic
#plot(allEffects(prey_ratio_mod)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = prey_ratio_mod) 
#WARNING: THESE BREAK R RIGHT NOW
#plot(simulationOutput, asFactor=TRUE) #these look good
#testZeroInflation(simulationOutput) #these are zero inflated
#testDispersion(simulationOutput) #not overdispersed
#hist(ASV_totals_samples$prey)
#boxplot visualization
#ggplot(ASV_totals_samples_ratio, aes(x = pipeline, y = prey/tot_oth)) +
#  geom_boxplot() +
#  labs(x = "Pipeline", y = "Prey/Other ASV ratio per sample") +
#  theme_bw()
#8. predator/prey read abundance ####
#read_totals_samples using the tot_oth category
#first, delete columns with zero in the tot_oth category
#read_totals_samples_ratio <- read_totals_samples[which(read_totals_samples$tot_oth > 0),]

#then create a model fitting the number of ASVs to pipeline grouped by sample
#prey_rratio_mod <- glmmTMB(prey ~ pipeline + (1|sample),
#                          data = read_totals_samples_ratio,
#                          family = "genpois",
#                          ziformula = ~1,
#                          offset = log(tot_oth))
#the null model where pipeline does not matter
#prey_rratio_mod_null <- glmmTMB(prey ~ 1 + (1|sample),
#                               data = read_totals_samples_ratio,
#                               family = "genpois",
#                               ziformula = ~1,
#                               offset = log(tot_oth))

#model.emm <- emmeans(prey_rratio_mod, "pipeline")
#pairs(model.emm)
#cleaning significantly increased the ratio of prey to other reads
#for both pipelines
#there is no significant difference in prey read ratios between either
#cleaned or uncleaned datasets

#x_rrat <- residuals(prey_rratio_mod)
#plotNormalHistogram(x_rrat)
#residuals look relatively normal
#plot(residuals(prey_rratio_mod)) #residuals look fairly homoskedastic
#plot(allEffects(prey_rratio_mod)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = prey_rratio_mod) 
#plot(simulationOutput, asFactor=TRUE) #these look good
#testZeroInflation(simulationOutput) #these are zero inflated...
#testDispersion(simulationOutput) #not overdispersed

#boxplot visualization
#ggplot(read_totals_samples_ratio, aes(x = pipeline, y = prey/tot_oth)) +
#  geom_boxplot() +
#  labs(x = "Pipeline", y = "Prey/Other read ratio per sample") +
#  theme_bw()
#9. amount of known diet item for lab-fed (ASVs) ####
#subset samples of interest from dataframe (lab fed) using these dataframes:
#d2_uc_id_tab
#d2_c_id_tab
#u3_uc_id_tab
#d2_ASVct <- d2_uc %>%
#  summarize_all(~sum(. != 0)) %>%
#  gather(sample, ASVs_d2uc, ASV:QC1, factor_key = TRUE) 
#u3_c_id_tab

metadata <- read.csv(here("data", "Sample_Metadata.csv"))
lab <- as.character(metadata$sample[which(metadata$Source == "LAB")])
 
d2_uc_lab <- d2_uc_id_tab %>%
  dplyr::select(all_of(lab), ID_bold, ID_ncbi, ID_nondiet) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi, -ID_nondiet) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, prey_d2uc, HEV07:HEV29, factor_key=TRUE)

d2_c_lab <- d2_c_id_tab %>%
  dplyr::select(all_of(lab), ID_bold, ID_ncbi, ID_nondiet) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi, -ID_nondiet) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, prey_d2c, HEV07:HEV29, factor_key=TRUE)

u3_uc_lab <- u3_uc_id_tab %>%
  dplyr::select(all_of(lab), ID_bold, ID_ncbi) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, prey_u3uc, HEV07:HEV29, factor_key=TRUE)

u3_c_lab <- u3_c_id_tab %>%
  dplyr::select(all_of(lab), ID_bold, ID_ncbi) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, prey_u3c, HEV07:HEV29, factor_key=TRUE)

#join these two together
known_prey_ASVs <- d2_uc_lab %>%
  full_join(d2_c_lab, by = "sample") %>%
  full_join(u3_uc_lab, by = "sample") %>%
  full_join(u3_c_lab, by = "sample")


#delete those that are zero across all pipelines (as my model is being dumb)
known_prey_ASVs <- known_prey_ASVs %>%
  filter(prey_u3uc > 0)

known_ASVs_long <- known_prey_ASVs %>%
  gather(pipeline, value, prey_d2uc:prey_u3c)
#this model is not working when i group by sample... :'(
#this is a BAD model right now...
#known_ASVs_mod <- glmmTMB(value ~ pipeline,
#                       data = known_ASVs_long,
#                       ziformula = ~1,
#                       family = "genpois")

#known_ASVs_mod_null <- glmmTMB(value ~ 1,
#                            data = known_ASVs_long,
#                            ziformula = ~1,
#                            family = "genpois")

#model.emm <- emmeans(known_ASVs_mod, "pipeline")
#pairs(model.emm)
#ACCORDING TO BAD MODEL:
#the number of ASVs assigned to known prey did not vary
#between cleaned-uncleaned pairs of the same pipeline
#however, the number varies across pipelines, both for clean
#and unclean pairs, with known prey ASVs being higher in
#both dadd2 clean and unclean pipelines

#x_ASVs <- residuals(known_ASVs_mod)
#plotNormalHistogram(x_ASVs)
#residuals look relatively normal
#summary(known_ASVs_mod)
#plot(residuals(known_ASVs_mod)) #residuals look fairly homoskedastic
#plot(allEffects(known_ASVs_mod)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
#simulationOutput <- simulateResiduals(fittedModel = known_ASVs_mod) 
#plot(simulationOutput, asFactor=FALSE) #these look good
#testZeroInflation(simulationOutput) #not zero inflated
#testDispersion(simulationOutput) #not overdispersed
#AICc(known_ASVs_mod, known_ASVs_mod_null) #model with pipeline is better than model without

#and that boxplot here
#ggplot(known_ASVs_long, aes(x = pipeline, y = value)) +
#  geom_boxplot() +
#  labs(x = "Pipeline", y = "Known ASVs per sample") +
#  theme_bw()


#10. amount of knwon diet item for lab-fed (num of reads) ####
d2_uc_rlab <- d2_uc_id_tab %>%
  dplyr::select(all_of(lab), ASV, ID_bold, ID_ncbi, ID_nondiet) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi, -ID_nondiet) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(reads_d2_uc = sum(reads))

d2_c_rlab <- d2_c_id_tab %>%
  dplyr::select(all_of(lab), ASV, ID_bold, ID_ncbi, ID_nondiet) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi, -ID_nondiet) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(reads_d2_c = sum(reads))

u3_uc_rlab <- u3_uc_id_tab %>%
  dplyr::select(all_of(lab), ASV, ID_bold, ID_ncbi) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(reads_u3_uc = sum(reads))

u3_c_rlab <- u3_c_id_tab %>%
  dplyr::select(all_of(lab), ASV, ID_bold, ID_ncbi) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(reads_u3_c = sum(reads))

#join these two together
known_prey_reads <- d2_uc_rlab %>%
  full_join(d2_c_rlab, by = "sample") %>%
  full_join(u3_uc_rlab, by = "sample") %>%
  full_join(u3_c_rlab, by = "sample")

known_reads_long <- known_prey_reads %>%
  gather(pipeline, value, reads_d2_uc:reads_u3_c)

known_reads_long <- known_reads_long %>%
  left_join(metadata, by = "sample")

known_reads_long$Ster <- as.factor(ifelse(known_reads_long$Sterilized == "NS", 0, 1))

#now run the model
known_reads_mod <- glmmTMB(value ~ pipeline*Ster + (1|sample),
                          data = known_reads_long,
                          family = "genpois",
                          REML = FALSE)
known_reads_mod2 <- glmmTMB(value ~ pipeline + Ster + (1|sample),
                           data = known_reads_long,
                           family = "genpois",
                           REML = FALSE)
known_reads_mod3 <- glmmTMB(value ~ pipeline + (1|sample),
                           data = known_reads_long,
                           family = "genpois",
                           REML = FALSE)
known_reads_mod4 <- glmmTMB(value ~ Ster + (1|sample),
                           data = known_reads_long,
                           family = "genpois",
                           REML = FALSE)

known_reads_mod_null <- glmmTMB(value ~ 1 + (1|sample),
                               data = known_reads_long,
                               family = "genpois")

AICc(known_reads_mod, known_reads_mod2, known_reads_mod3, known_reads_mod4,
     known_reads_mod_null)
#model 1 is the best fit based on AIC values.
known_reads_mod <- glmmTMB(value ~ pipeline*Ster + (1|sample),
                           data = known_reads_long,
                           family = "genpois")

model.emm_knownprey <- emmeans(known_reads_mod, pairwise ~ pipeline | Ster)
model.emm_knownpreys <- emmeans(known_reads_mod, pairwise ~ Ster | pipeline)
pairs(model.emm_knownprey)
pairs(model.emm_knownpreys)
emmeans(known_reads_mod, pairwise ~ pipeline | Ster)
#Sterilization significantly decreased the number of known prey reads for all
#pipelines
#There are significantly more reads in UNOISE datasets for both sterilized and 
#unsterilized individuals for both cleaned and uncleaned pipeilnes.

x_knownprey <- residuals(known_reads_mod)
plotNormalHistogram(x_knownprey)
#residuals look relatively normal
summary(known_reads_mod)
plot(residuals(known_reads_mod)) #eek
plot(allEffects(known_reads_mod)) #some version of this I think would be a good paper figure
#but ideally prettier

#checking other model assumptions here:
simulationOutput_knownprey <- simulateResiduals(fittedModel = known_reads_mod) 
knownpreyfit <- plot(simulationOutput_knownprey, asFactor=FALSE) #these look good
knownpreyzi <- testZeroInflation(simulationOutput_knownprey) #not zero inflated
knownpreyod <- testDispersion(simulationOutput_knownprey) #maybe overdispersed???

#and that boxplot here
known_prey_plot <- ggplot(known_reads_long, aes(x = pipeline, y = value, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Known prey reads per sample", title = "Number of known prey reads per sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979"))

#Intermediate: Phylogenetic trees of prey ASVs ####
#currently I have not included phylogenies from BOLD in this, 
#which... maybe doesn't need to happen? I dont' know sigh
#according to the taxonomy databases, there is ONE ASV of prey
#origin that BOLD matches that MEGAN doesn't ACROSS ALL pipelines
#so I think i'm okay with no BOLD

#SOME PD resources to start:
#ses.pd in the picante package:
#https://rdrr.io/rforge/picante/man/pd.html
#general phylogeny tricks here:
#http://www.phytools.org/eqg/Exercise_3.2/

#I created a meta-tree of both hits from the unclean dada2 and unoise
#pipelines in MEGAN (since cleaning didn't add any diversity, this represents
#the total possible diversity of any pipeline)
#I created a tree at the Order level and at the Family level
#then, MEGAN didn't compute branch lengths, so I did
ord_tree <- read.tree(here("6_taxonomic_assignment", "bilateria_order.tre"))
ord_tree1 <- compute.brlen(ord_tree, method = "Grafen")
ord_tree2 <- di2multi(ord_tree1) #this collapses multichotomies so that branch
#lengths aren't domianted by zeros

#same for family tree
fam_tree <- read.tree(here("6_taxonomic_assignment", "bilateria_family.tre"))
fam_tree1 <- compute.brlen(fam_tree, method = "Grafen")
fam_tree2 <- di2multi(fam_tree1)
fam_tree2$edge.length
#For community data:
#will be using a manipulated and transposed version of each
#pipeline's _ncbi dataframe attached to the community matrix "pipeline"_"c/uc"

#the general direction here is to subset all reads that include
#prey and not anything else within my community matrices
#going to be using the NCBI matches and the basic community matrix
#filtering out everything but prey

#11. Phylogenetic Diversity of Prey ASVs at ORDER LEVEL####
#DADA2 UNCLEAN at ORDER LEVEL##
#need to combine and add all ASV reads across all samples to
#get a meta-ASV count for the entire pipeline
#this following pipe does that
d2_uc_sum <- d2_uc %>%
  dplyr::select(-CL1, -CL4, -QC1) %>% #remove controls
  gather(sample, reads, HEV100:HEV89) %>% #gather by ASV and sample
  group_by(ASV) %>% #group by ASV
  summarize(reads= sum(reads)) %>% #summarize the total reads in each ASV
  filter(reads > 0) #filter only those that are greater than 0

d2_uc_ncbi$ASV <- as.character(d2_uc_ncbi$ASV)
d2_uc_sum$ASV <- as.character(d2_uc_sum$ASV)
d2_uc_phy_tab <- d2_uc_ncbi %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(ASV, ID, Order, Family) %>%
  left_join(d2_uc_sum, by = "ASV") %>%
  dplyr::select(-ASV, -ID) #55 TOTAL

#create an order dataframe first by selecting only order column and
#then deleting everything in the order column that is blank 
d2_uc_phy_tab_ord <- d2_uc_phy_tab %>%
  dplyr::select(-Family) %>%
  na_if("") %>%
  drop_na(Order) #46 observations, so 9 had no order ID

#gather by sample and group by sample and Order to summarize how many are in each order
#(we are removing duplicate orders by combining them here). Then spread back out to a 
#community matrix with samples as the columns and orders in those samples as rows
d2_uc_phy_tab_ord <- d2_uc_phy_tab_ord %>%
  group_by(Order) %>%
  summarize(reads = sum(reads))  #14 ORDERS

#set rownames to ID column, and then remove that column
rownames(d2_uc_phy_tab_ord) <- d2_uc_phy_tab_ord$Order
d2_uc_phy_tab_ord <- d2_uc_phy_tab_ord %>%
  dplyr::select(-Order)

#transpose for phylogenetic analyses
d2_uc_phy_tab_ord_t <- as.data.frame(t(d2_uc_phy_tab_ord))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_d2uc_ord <- pd(d2_uc_phy_tab_ord_t, ord_tree2, include.root=TRUE)
pd_d2uc_ord$pipeline <- "d2_uc"
#UNOISE 3 UC
u3_uc_sum <- u3_uc %>%
  dplyr::select(-CL1, -CL4, -QC1, -NEG) %>% #remove controls
  gather(sample, reads, HEV15:HEV83) %>% #gather by ASV and sample
  group_by(ASV) %>% #group by ASV
  summarize(reads= sum(reads)) %>% #summarize the total reads in each ASV
  filter(reads > 0) #filter only those that are greater than 0

u3_uc_ncbi$ASV <- as.character(u3_uc_ncbi$ASV)
u3_uc_sum$ASV <- as.character(u3_uc_sum$ASV)
u3_uc_phy_tab <- u3_uc_ncbi %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(ASV, ID, Order, Family) %>%
  left_join(u3_uc_sum, by = "ASV") %>%
  dplyr::select(-ASV, -ID) #52 TOTAL

#create an order dataframe first by selecting only order column and
#then deleting everything in the order column that is blank 
u3_uc_phy_tab_ord <- u3_uc_phy_tab %>%
  dplyr::select(-Family) %>%
  na_if("") %>%
  drop_na(Order) #43 observations

#gather by sample and group by sample and Order to summarize how many are in each order
#(we are removing duplicate orders by combining them here). Then spread back out to a 
#community matrix with samples as the columns and orders in those samples as rows
u3_uc_phy_tab_ord <- u3_uc_phy_tab_ord %>%
  group_by(Order) %>%
  summarize(reads = sum(reads, na.rm = TRUE))  #13 ORDERS

#set rownames to ID column, and then remove that column
rownames(u3_uc_phy_tab_ord) <- u3_uc_phy_tab_ord$Order
u3_uc_phy_tab_ord <- u3_uc_phy_tab_ord %>%
  dplyr::select(-Order)

#transpose for phylogenetic analyses
u3_uc_phy_tab_ord_t <- as.data.frame(t(u3_uc_phy_tab_ord))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_u3uc_ord <- pd(u3_uc_phy_tab_ord_t, ord_tree2, include.root=TRUE)
pd_u3uc_ord$pipeline <- "u3_uc"
#DADA2 CLEAN ##
#need to combine and add all ASV reads across all samples to
#get a meta-ASV count for the entire pipeline
#this following pipe does that
d2_c_sum <- d2_c %>%
  gather(sample, reads, HEV100:HEV89) %>% #gather by ASV and sample
  group_by(ASV) %>% #group by ASV
  summarize(reads= sum(reads)) %>% #summarize the total reads in each ASV
  filter(reads > 0) #filter only those that are greater than 0

d2_c_ncbi$ASV <- as.character(d2_c_ncbi$ASV)
d2_c_sum$ASV <- as.character(d2_c_sum$ASV)
d2_c_phy_tab <- d2_c_ncbi %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(ASV, ID, Order, Family) %>%
  left_join(d2_c_sum, by = "ASV") %>%
  dplyr::select(-ASV, -ID) #55 TOTAL

#create an order dataframe first by selecting only order column and
#then deleting everything in the order column that is blank 
d2_c_phy_tab_ord <- d2_c_phy_tab %>%
  dplyr::select(-Family) %>%
  na_if("") %>%
  drop_na(Order) #49 observations

#gather by sample and group by sample and Order to summarize how many are in each order
#(we are removing duplicate orders by combining them here). Then spread back out to a 
#community matrix with samples as the columns and orders in those samples as rows
d2_c_phy_tab_ord <- d2_c_phy_tab_ord %>%
  group_by(Order) %>%
  summarize(reads = sum(reads))  #14 ORDERS

#set rownames to ID column, and then remove that column
rownames(d2_c_phy_tab_ord) <- d2_c_phy_tab_ord$Order
d2_c_phy_tab_ord <- d2_c_phy_tab_ord %>%
  dplyr::select(-Order)

#transpose for phylogenetic analyses
d2_c_phy_tab_ord_t <- as.data.frame(t(d2_c_phy_tab_ord))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_d2c_ord <- pd(d2_c_phy_tab_ord_t, ord_tree2, include.root=TRUE)
pd_d2c_ord$pipeline <- "d2_c"
#UNOISE 3 UC
u3_c_sum <- u3_c %>%
  gather(sample, reads, HEV07:HEV99) %>% #gather by ASV and sample
  group_by(ASV) %>% #group by ASV
  summarize(reads= sum(reads)) %>% #summarize the total reads in each ASV
  filter(reads > 0) #filter only those that are greater than 0

u3_c_ncbi$ASV <- as.character(u3_c_ncbi$ASV)
u3_c_sum$ASV <- as.character(u3_c_sum$ASV)
u3_c_phy_tab <- u3_c_ncbi %>%
  filter(taxonomy == "prey") %>%
  dplyr::select(ASV, ID, Order, Family) %>%
  left_join(u3_c_sum, by = "ASV") %>%
  dplyr::select(-ASV, -ID) #48 TOTAL

#create an order dataframe first by selecting only order column and
#then deleting everything in the order column that is blank 
u3_c_phy_tab_ord <- u3_c_phy_tab %>%
  dplyr::select(-Family) %>%
  na_if("") %>%
  drop_na(Order) #43 observations

#gather by sample and group by sample and Order to summarize how many are in each order
#(we are removing duplicate orders by combining them here). Then spread back out to a 
#community matrix with samples as the columns and orders in those samples as rows
u3_c_phy_tab_ord <- u3_c_phy_tab_ord %>%
  group_by(Order) %>%
  summarize(reads = sum(reads, na.rm = TRUE))  #11 ORDERS

#set rownames to ID column, and then remove that column
rownames(u3_c_phy_tab_ord) <- u3_c_phy_tab_ord$Order
u3_c_phy_tab_ord <- u3_c_phy_tab_ord %>%
  dplyr::select(-Order)

#transpose for phylogenetic analyses
u3_c_phy_tab_ord_t <- as.data.frame(t(u3_c_phy_tab_ord))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_u3c_ord <- pd(u3_c_phy_tab_ord_t, ord_tree2, include.root=TRUE)
pd_u3c_ord$pipeline <- "u3_c"

#bind together
order_pd <- bind_rows(pd_d2uc_ord, pd_u3uc_ord, pd_d2c_ord, pd_u3c_ord)
order_pd <- order_pd %>%
  gather(measure, value, PD:SR)
#visualize, need to make prettier, and somehow standardize by species richness?
order_pd_graph <- ggplot(order_pd, aes(x = pipeline, y = value, color = measure)) +
  geom_point(size = 4) +
  labs(x = "Pipeline", y = "Order phylogenetic diversity and order richness", title = "Order PD and richness by pipeline") +
  theme_bw() 

#12.Phylogenetic Diversity of prey ASVs at FAMILY LEVEL #####
#using fam_tree2
#DADA2 Unclean
#create a family dataframe first by selecting only family column and
#then deleting everything in the family column that is blank 
d2_uc_phy_tab_fam <- d2_uc_phy_tab %>%
  dplyr::select(-Order) %>%
  na_if("") %>%
  drop_na(Family) #28 observations

#gather by sample and group by sample and Family to summarize how many are in each family
#(we are removing duplicate families by combining them here). Then spread back out to a 
#community matrix with samples as the columns and families in those samples as rows
d2_uc_phy_tab_fam <- d2_uc_phy_tab_fam %>%
  group_by(Family) %>%
  summarize(reads = sum(reads)) #17 families

#set rownames to ID column, and then remove that column
rownames(d2_uc_phy_tab_fam) <- d2_uc_phy_tab_fam$Family
d2_uc_phy_tab_fam <- d2_uc_phy_tab_fam %>%
  dplyr::select(-Family)

#transpose for phylogenetic analyses
d2_uc_phy_tab_fam_t <- as.data.frame(t(d2_uc_phy_tab_fam))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_d2uc_fam <- pd(d2_uc_phy_tab_fam_t, fam_tree2, include.root=TRUE)
pd_d2uc_fam$pipeline <- "d2_uc"
#UNOISE Unclean##
#create a family dataframe first by selecting only family column and
#then deleting everything in the family column that is blank 
u3_uc_phy_tab_fam <- u3_uc_phy_tab %>%
  dplyr::select(-Order) %>%
  na_if("") %>%
  drop_na(Family) #25 observations

#gather by sample and group by sample and Family to summarize how many are in each family
#(we are removing duplicate families by combining them here). Then spread back out to a 
#community matrix with samples as the columns and families in those samples as rows
u3_uc_phy_tab_fam <- u3_uc_phy_tab_fam %>%
  group_by(Family) %>%
  summarize(reads = sum(reads, na.rm=TRUE)) #16 families

#set rownames to ID column, and then remove that column
rownames(u3_uc_phy_tab_fam) <- u3_uc_phy_tab_fam$Family
u3_uc_phy_tab_fam <- u3_uc_phy_tab_fam %>%
  dplyr::select(-Family)

#transpose for phylogenetic analyses
u3_uc_phy_tab_fam_t <- as.data.frame(t(u3_uc_phy_tab_fam))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_u3uc_fam <- pd(u3_uc_phy_tab_fam_t, fam_tree2, include.root=TRUE)
pd_u3uc_fam$pipeline <- "u3_uc"
#DADA2 CLEAN
#create a family dataframe first by selecting only family column and
#then deleting everything in the family column that is blank 
d2_c_phy_tab_fam <- d2_c_phy_tab %>%
  dplyr::select(-Order) %>%
  na_if("") %>%
  drop_na(Family) #31 observations

#gather by sample and group by sample and Family to summarize how many are in each family
#(we are removing duplicate families by combining them here). Then spread back out to a 
#community matrix with samples as the columns and families in those samples as rows
d2_c_phy_tab_fam <- d2_c_phy_tab_fam %>%
  group_by(Family) %>%
  summarize(reads = sum(reads, na.rm=TRUE)) #17 families

#set rownames to ID column, and then remove that column
rownames(d2_c_phy_tab_fam) <- d2_c_phy_tab_fam$Family
d2_c_phy_tab_fam <- d2_c_phy_tab_fam %>%
  dplyr::select(-Family)

#transpose for phylogenetic analyses
d2_c_phy_tab_fam_t <- as.data.frame(t(d2_c_phy_tab_fam))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_d2c_fam <- pd(d2_c_phy_tab_fam_t, fam_tree2, include.root=TRUE)
pd_d2c_fam$pipeline <- "d2_c"
#UNOISE CLEAn##
#create a family dataframe first by selecting only family column and
#then deleting everything in the family column that is blank 
u3_c_phy_tab_fam <- u3_c_phy_tab %>%
  dplyr::select(-Order) %>%
  na_if("") %>%
  drop_na(Family) # 24 observations

#gather by sample and group by sample and Family to summarize how many are in each family
#(we are removing duplicate families by combining them here). Then spread back out to a 
#community matrix with samples as the columns and families in those samples as rows
u3_c_phy_tab_fam <- u3_c_phy_tab_fam %>%
  group_by(Family) %>%
  summarize(reads = sum(reads, na.rm=TRUE)) #15 families

#set rownames to ID column, and then remove that column
rownames(u3_c_phy_tab_fam) <- u3_c_phy_tab_fam$Family
u3_c_phy_tab_fam <- u3_c_phy_tab_fam %>%
  dplyr::select(-Family)

#transpose for phylogenetic analyses
u3_c_phy_tab_fam_t <- as.data.frame(t(u3_c_phy_tab_fam))

#finally, use this community matrix for each pipeline plus the 
#meta phylogeny from both Unoise and DADA2 to calculate Faith's PD
#for each sample
pd_u3c_fam <- pd(u3_c_phy_tab_fam_t, fam_tree2, include.root=TRUE)
pd_u3c_fam$pipeline <- "u3_c"

family_pd <- bind_rows(pd_d2uc_fam, pd_u3uc_fam, pd_d2c_fam, pd_u3c_fam)
family_pd <- family_pd %>%
  gather(measure, value, PD:SR)

family_pd_graph <- ggplot(family_pd, aes(x = pipeline, y = value, color = measure)) +
  geom_point(size = 4) +
  labs(x = "Pipeline", y = "Family phylogenetic diversity and family richness", title = "Family PD and richness by pipeline") +
  theme_bw()

#plot trees####
d2uc_tree <- d2_uc_phy_tab_ord %>%
  rename("d2uc" = "reads") %>%
  rownames_to_column("Order")

u3uc_tree <- u3_uc_phy_tab_ord %>%
  rename("u3uc" = "reads") %>%
  rownames_to_column("Order")

d2c_tree <- d2_c_phy_tab_ord %>%
  rename("d2c" = "reads") %>%
  rownames_to_column("Order")

u3c_tree <- u3_c_phy_tab_ord %>%
  rename("u3c" = "reads") %>%
  rownames_to_column("Order")

order_tree_pipe <- u3uc_tree %>%
  full_join(d2uc_tree, by = "Order") #%>%
  #left_join(u3c_tree, by = "Order") %>%
  #left_join(d2c_tree, by = "Order")

#set NA to zeros
order_tree_pipe[is.na(order_tree_pipe)] <- 0

#set rownames to order and then delete that column
rownames(order_tree_pipe) <- order_tree_pipe$Order
order_tree_pipe <- order_tree_pipe %>%
  dplyr::select(-Order)

#set to a matrix
order_tree_pipe <- as.matrix(order_tree_pipe)


#import taxonomies, set rownames, and make matrix
order_taxonomies <- read.csv(here("data", "order_phylos.csv"))
rownames(order_taxonomies) <- order_taxonomies$Order
order_taxonomies <- as.matrix(order_taxonomies)

#from phyloseq package
order_ASV <- otu_table(order_tree_pipe, taxa_are_rows = TRUE)
order_taxa <- tax_table(order_taxonomies)
physeq_order <-  phyloseq(order_ASV, order_taxa, ord_tree2)

order_phy_graph <- plot_tree(physeq_order, nodelabf=nodeplotblank, plot.margin=0.6, label.tips = "taxa_names", shape = "samples")

#Family Tree now
d2uc_treef <- d2_uc_phy_tab_fam %>%
  rename("d2uc" = "reads") %>%
  rownames_to_column("Family")

u3uc_treef <- u3_uc_phy_tab_fam %>%
  rename("u3uc" = "reads") %>%
  rownames_to_column("Family")

d2c_treef <- d2_c_phy_tab_fam %>%
  rename("d2c" = "reads") %>%
  rownames_to_column("Family")

u3c_treef <- u3_c_phy_tab_fam %>%
  rename("u3c" = "reads") %>%
  rownames_to_column("Family")

family_tree_pipe <- u3uc_treef %>%
  full_join(d2uc_treef, by = "Family") #%>%
  #left_join(u3c_treef, by = "Family") %>%
  #left_join(d2c_treef, by = "Family")

#set NA to zeros
family_tree_pipe[is.na(family_tree_pipe)] <- 0

#set rownames to order and then delete that column
rownames(family_tree_pipe) <- family_tree_pipe$Family
family_tree_pipe <- family_tree_pipe %>%
  dplyr::select(-Family)

#set to a matrix
family_tree_pipe <- as.matrix(family_tree_pipe)

#import taxonomies, set rownames, and make matrix
family_taxonomies <- read.csv(here("data", "family_phylos.csv"))
rownames(family_taxonomies) <- family_taxonomies$Family
family_taxonomies <- as.matrix(family_taxonomies)

#from phyloseq package
family_ASV <- otu_table(family_tree_pipe, taxa_are_rows = TRUE)
family_taxa <- tax_table(family_taxonomies)
physeq_family <-  phyloseq(family_ASV, family_taxa, fam_tree2)

family_phy_graph <- plot_tree(physeq_family, nodelabf = nodeplotblank, plot.margin=0.6, label.tips = "taxa_names", shape = "samples")
family_phy_graph
