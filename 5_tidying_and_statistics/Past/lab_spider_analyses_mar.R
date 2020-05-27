#Intro####
#This code looks at some patterns in mesocosm predators, specifically
#looking for the detection and abundance of prey DNA in predators 
#which we have fed a known diet item to
#Does a specific sequencing pipeline detect more prey ASVs or greater
#prey read abundances? Does sterilization alter prey ASV detection or 
#read abundance by removing surface contamination?

#we will be doing two analyses here:
#  1.	Known prey detection (ASV presence-absence)
#  2.	Known prey abundance (read abundance)


#Load required packages####
library(here) #tidy data
library(tidyverse) #tidy data
library(ggplot2) #visualize data
library(emmeans) #marginal means of models
library(rcompanion) #pairwise comparisions
library(effects) #dotplots and stuff
library(lattice) #again dotplots
library(DHARMa) #model diagnostics
library(MuMIn) #model diagnostics
library(glmmTMB) #mixed models
library(lme4) #mixed models
library(performance) #for binomail model fit
library(see) #for binomial model fit
library(ggeffects) #marginal means
library(coin)
library(MASS)

#Import ASV tables####
d2_comm <- read.csv(here("data", "ASV_tables", "dada2_uc_asv_tab.tsv"), sep = "\t")
u3_comm <- read.csv(here("data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = '\t')

#rename X to ASV in all these tables.
d2_comm <- rename(d2_comm, "ASV" = "X")
u3_comm <- rename(u3_comm, "ASV" = "X.OTU.ID")

#rename all the sample names across dataframes for consistency
colnames(d2_comm) <- sapply(str_split(colnames(d2_comm), "_"), function(x){return(x[[1]])})
colnames(d2_comm) <- str_remove(colnames(d2_comm), "\\.")

d2_comm <- rename(d2_comm, "CL1" = "CL12")
d2_comm <- rename(d2_comm, "CL4" = "CL42")

colnames(u3_comm) <- sapply(str_split(colnames(u3_comm), "S"), function(x){return(x[[1]])})
u3_comm <- rename(u3_comm, "ASV" = "A")

#sample metadata
metadata <- read.csv(here("data", "Sample_Metadata.csv"))

#ASSIGN taxonomies ####
#Now we need to assign taxonomies so we can select the taxonomies of the known prey item
#(Acrididae, Oxya)

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

#this binds them all to the community matrix for analyses
d2_id_tab <- d2_comm %>%
  left_join(d2_id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
d2_id_tab$taxonomy <- replace_na(d2_id_tab$taxonomy, "no hit")

#write this to a file for importing into the field spider analyses
d2_taxonomies <- d2_id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, ID_nondiet, Order_nondiet, taxonomy)
write.csv(d2_taxonomies, here("data", "d2_uc_tax_ass.csv"))

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

#this binds them all to the community matrix for analyses
u3_id_tab <- u3_comm %>%
  left_join(u3_id_all, by = "ASV")

#set NA taxonomies from join with "no hit" since they didn't match to any
#taxonomic assignments within our prey, predator, or non-prey categories
u3_id_tab$taxonomy <- replace_na(u3_id_tab$taxonomy, "no hit")

#write to a file for import into field analyses
u3_taxonomies <- u3_id_tab %>%
  dplyr::select(ASV, ID_bold, Order_bold, ID_ncbi, Order_ncbi, taxonomy)
write.csv(u3_taxonomies, here("data", "u3_uc_tax_ass.csv"))

## ASVs: Subset lab spiders####
lab <- as.character(metadata$sample[which(metadata$Source == "LAB")])

d2_lab <- d2_id_tab %>%
  dplyr::select(all_of(lab), ID_bold, ID_ncbi, ID_nondiet) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi, -ID_nondiet) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, prey_d2, HEV07:HEV29, factor_key=TRUE)

u3_lab <- u3_id_tab %>%
  dplyr::select(all_of(lab), ID_bold, ID_ncbi) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi) %>%
  summarize_all(~sum(. != 0)) %>%
  gather(sample, prey_u3, HEV07:HEV29, factor_key=TRUE)

#Known prey ASVs: Analyses ####
#join these two together
known_prey_ASVs <- d2_lab %>%
  full_join(u3_lab, by = "sample")

#Known ASVs ~ Pipeline ####
#First compare pipelines using a McNemar paired binary test
#from https://rcompanion.org/handbook/H_05.html
#then create a presence/absence binary variable for this
pres_D2 <- ifelse(known_prey_ASVs$prey_d2 > 0, 1, 0)
pres_U3 <- ifelse(known_prey_ASVs$prey_u3 > 0, 1, 0)

#and add to a new dataframe for the McNemar test
ASVs_mcnemar <- known_prey_ASVs %>%
  mutate(D2 = pres_D2, U3 = pres_U3) %>%
  dplyr::select(-prey_d2, -prey_u3)

#this just determines which samples have shared presence, absence, and
#differ in presence/absenece
both_pres <- ASVs_mcnemar[which(ASVs_mcnemar$D2 == ASVs_mcnemar$U3 & ASVs_mcnemar$D2 != 0),]
both_abs <- ASVs_mcnemar[which(ASVs_mcnemar$D2 == ASVs_mcnemar$U3 & ASVs_mcnemar$D2 == 0),]
d2_pres <- ASVs_mcnemar[which(ASVs_mcnemar$D2 > 0 & ASVs_mcnemar$U3 == 0),]
u3_pres <- ASVs_mcnemar[which(ASVs_mcnemar$U3 > 0 & ASVs_mcnemar$D2 == 0),]

#the values for these to populate the table below
nrow(both_pres) #15
nrow(both_abs) #3
nrow(d2_pres) #0
nrow(u3_pres) #1

#create a table
Input <- ("DADA2       UNOISE3.pres   UNOISE3.abs
        DADA2.pres     15          0
        DADA2.abs     1         3
        ")

#turn that into a matrix
Matrix.1 <-  as.matrix(read.table(textConnection(Input),
                                header=TRUE,
                                row.names=1))

#this sum should be the number of samples
sum(Matrix.1)

Matrix.1
#this is mcnemar test
nominalSymmetryTest(Matrix.1,
                    digits = 3)
#this is mcnemar test again, with p-value of 1. This means that
#the matrix is symmetrical, and the two pipelines do not differ


mcnemar.test(Matrix.1)
#there is no difference in detection between the two bioinformatics pipelines,
#let's visualize that
known_ASVs_long <- known_prey_ASVs %>%
  gather(pipeline, value, prey_d2:prey_u3)

present <- known_ASVs_long %>%
  group_by(pipeline) %>%
  filter(value > 0) %>%
  summarise(present = n()) 

total <- known_ASVs_long %>%
  group_by(pipeline) %>%
  summarise(total = n()) %>%
  dplyr::select(pipeline, total)

ASVs_known <- present %>%
  left_join(total, by = "pipeline")

ASVs_known <- ASVs_known %>%
  mutate(presence = present/total)

ggplot(ASVs_known, aes(x = pipeline, y = presence)) +
  geom_bar(stat="identity", position= "dodge") +theme_bw() +
  labs(x = "Pipeline", y = "Percent Known Prey ASV Detection")

#Known ASVs ~ Surface sterilized ####
#NOW, we can look at surface sterilization effects, which we 
#will do with a binomial GLMM
#combine with sample metadata
known_ASVs_long <- known_ASVs_long %>%
  left_join(metadata, by = "sample")

known_ASVs_long$surf_ster <- as.factor(ifelse(known_ASVs_long$Sterilized == "NS", 0, 1))
known_ASVs_long$Presence <- ifelse(known_ASVs_long$value == 0, 0, 1)

ASV_ster_mod <- glmmTMB(Presence ~ surf_ster + (1|pipeline),
                    data = known_ASVs_long, 
                    family = binomial)


simulationOutput_knownASV <- simulateResiduals(fittedModel = ASV_ster_mod) 
knownASVfit <- plot(simulationOutput_knownASV, asFactor=FALSE)
binned_residuals(ASV_ster_mod)

ASV_ster_null <- glmmTMB(Presence ~ 1 + (1|pipeline),
                    data = known_ASVs_long, 
                    family = "binomial")

AICc(ASV_ster_mod, ASV_ster_null) #model with sterilization is better fit
anova(ASV_ster_mod, ASV_ster_null)
plot(allEffects(ASV_ster_mod)) #sterilization REDUCES detection of known prey
model.emm_knownASV <- emmeans(ASV_ster_mod, "surf_ster")
pairs(model.emm_knownASV)

summary(ASV_ster_mod)

#let's visualize this then
present <- known_ASVs_long %>%
  group_by(pipeline, Sterilized) %>%
  filter(value > 0) %>%
  summarise(present = n()) 

total <- known_ASVs_long %>%
  group_by(pipeline, Sterilized) %>%
  summarise(total = n()) %>%
  dplyr::select(pipeline, Sterilized, total)

ASVs_known <- present %>%
  left_join(total, by = c("pipeline" = "pipeline", "Sterilized" = "Sterilized"))

ASVs_known <- ASVs_known %>%
  mutate(presence = present/total)

ggplot(ASVs_known, aes(x = pipeline, y = presence, fill = Sterilized)) +
  geom_bar(stat="identity", position= "dodge") +theme_bw() +
  labs(x = "Pipeline", y = "Percent Known Prey ASV Detection",
       title = "Percent of samples with known prey detected") +
  scale_fill_manual(values = c("#7B8D65","#F29979"))

#Reads: subset lab spiders####
d2_rlab <- d2_id_tab %>%
  dplyr::select(all_of(lab), ASV, ID_bold, ID_ncbi, ID_nondiet) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi, -ID_nondiet) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(reads_d2 = sum(reads))

u3_rlab <- u3_id_tab %>%
  dplyr::select(all_of(lab), ASV, ID_bold, ID_ncbi) %>%
  filter(ID_ncbi == "Acrididae") %>%
  dplyr::select(-ID_bold, -ID_ncbi) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(reads_u3 = sum(reads))

d2_tot_lab <- d2_id_tab %>%
  dplyr::select(all_of(lab), ASV) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(total = sum(reads)) %>%
  mutate(pipeline = "reads_d2")

u3_tot_lab <- u3_id_tab %>%
  dplyr::select(all_of(lab), ASV) %>%
  group_by(ASV) %>%
  gather(sample, reads, HEV07:HEV29, factor_key = TRUE) %>%
  group_by(sample) %>%
  summarize(total = sum(reads)) %>%
  mutate(pipeline = "reads_u3")

total_reads <- d2_tot_lab %>%
  bind_rows(u3_tot_lab)
# Known prey reads: analyses ####

#Known reads ~ Pipeline ####
#join these two together
known_prey_reads <- d2_rlab %>%
  full_join(u3_rlab, by = "sample")

wilcoxsign_test(known_prey_reads$reads_d2 ~ known_prey_reads$reads_u3)
wilcoxsign_test(known_prey_reads$reads_d2 ~ known_prey_reads$reads_u3,
                zero.method = "Wilcoxon")

#According to this, there *is* a difference between pipelines in the number
#of known prey reads detected. UNOISE3 has a significantly higher value
#for prey reads across samples. 

#let's visualize this now:
known_reads_long <- known_prey_reads %>%
  gather(pipeline, value, reads_d2:reads_u3) %>%
  left_join(metadata, by = "sample")

ggplot(known_reads_long, aes(x = pipeline, y = value)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Pipeline", y = "Read abundance per sample", 
       title = "Known prey read abundance by bioinformatics pipeline")

#or:
sum_reads <- known_reads_long %>%
  group_by(pipeline) %>%
  summarize(mean = mean(value), se = sd(value)/(sqrt(n_distinct(sample))))

ggplot(sum_reads, aes(x = pipeline, y = mean)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.1, size = 0.75) +
  theme_bw() +
  labs(x = "Pipeline", y = "Mean reads per sample",
       title = "Mean number of reads per sample by pipeline") +
  scale_x_discrete(labels=c("reads_d2" = "DADA2", "reads_u3" = "UNOISE3"))
  
#there are significantly different numbers of reads by pipeline, with
#UNOISE3 producing more reads per sample of known prey than DADA2

# Known reads ~ Sterilization ####
#sterilization, split by pipeline for this.
known_reads_long$Ster <- ifelse(known_reads_long$Sterilized == "NS", 0, 1)

reads_ster_mod <- glmmTMB(value ~ Ster + (1|pipeline),
                    data = known_reads_long,
                    family = "genpois",
                    REML = FALSE)

reads_ster_mod_null <- glmmTMB(value ~ 1 + (1|pipeline),
                      data = known_reads_long,
                      family = "genpois",
                      REML = FALSE)

AICc(reads_ster_mod, reads_ster_mod_null) #with sterilizatin is a better model

reads_ster_mod <- glmmTMB(value ~ Ster + (1|pipeline),
                    data = known_reads_long,
                    family = "genpois")

simulationOutput_knownprey <- simulateResiduals(fittedModel = reads_ster_mod) 
knownpreyfit <- plot(simulationOutput_knownprey, asFactor=FALSE) #these look good
knownpreyzi <- testZeroInflation(simulationOutput_knownprey) #not zero inflated
knownpreyod <- testDispersion(simulationOutput_knownprey) 

x_knownprey <- residuals(reads_ster_mod)
plotNormalHistogram(x_knownprey) #eek look yeah...
plot(residuals(reads_ster_mod)) #two weirdos but otherwise okay
#nb did not fix this problem, and not zero inflated so i dont' know

model.emm_knownprey <- emmeans(reads_ster_mod, "Ster")
pairs(model.emm_knownprey)
plot(allEffects(reads_ster_mod))

#sterilization seems to matter, and reduces known prey read abundance. HUH.

#and that boxplot here
ggplot(known_reads_long, aes(x = pipeline, y = value, fill = Sterilized)) +
  geom_boxplot() +
  labs(x = "Pipeline", y = "Known prey reads per sample", title = "Number of known prey reads per sample") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979")) +
  scale_y_log10()


reads_sum <- known_reads_long %>%
  group_by(pipeline, Sterilized) %>%
  summarize(mean = mean(value), se = sd(value)/n_distinct(sample))

dodge <- position_dodge(width=0.9)

ggplot(reads_sum, aes(x = pipeline, y = mean, fill = Sterilized)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + scale_fill_manual(values = c("#7B8D65","#F29979")) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), position = dodge, width= 0.25) +
  labs(x = "Pipeline", y = "Mean read abundance", 
       title = "Mean read abundance by pipeline and sterilization") +
  scale_x_discrete(labels=c("reads_d2" = "DADA2", "reads_u3" = "UNOISE3"))
  



