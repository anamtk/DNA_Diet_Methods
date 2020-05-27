#Intro ####
#This code looks at patterns in field-collected predators, specifically
#looking at the detection, diversity, and abundance of prey ASVs and 
#reads, looking at how well different bioinformatic pipelines pick up 
#prey DNA and also how sterilization might alter the prey detection,
#diversity, abundance, and composition of samples. I'll be taking a courser
#look at just the broad diversity and abundance of the samples, and then
#I'll be zooming in to look at which species may change between communities
#and whether our detection/abundance of core prey (ie core microbiome) is altered by
#either bioinformatics pipeline or sterilization.

#Field-collected Predators, 
#or: for predators with unknown diet, how do we ensure that prey 
#detected by metabarcoding methods is, indeed, prey, and not surface
#contamination?

#1.	Big picture
#i.	Intermediate: rarefy by lowest sample read first
#ii.	Intermediate: Concatenate by unique taxonomic assignments first
#b.	Total species richness of all prey taxonomies (total and by sample)
#i.	Total species (value only, no stats)
#ii.	By sample (glm(SR ~ Sterilization *algorithm, family = poisson))
#c.	Phylogenetic diversity of prey ASVs (just a value)
#i.	Not sure I will include this after all?
#  d.	Prey ASV read abundance (by sample)
#2.	Zooming in
#a.	Jaccard presence-absence community analyses of prey ASVs
#i.	Presence-absence PERMANOVA on raw data (glmer/glmmTMB(presence ~ sterilization+algorithm + (1+sterilization | species), family = binomial)
#                                           b.	Core prey (ala. Core microbiome) presence-absence analyses
#                                           i.	Presence-absence PERMANOVA on rarefied data (glmer/glmmTMB(presence ~ sterilization+algorithm + (1+sterilization | species), family = binomial)
                                                                                           

#Load required packages ####
library(here) #tidy data
library(tidyverse) #tidy data
library(ggplot2) #visualization
library(emmeans) #model marginal means
library(rcompanion) #model pairwise comparisons
library(effects) #model effect visualizations
library(lattice) #model effects visulaizations
library(DHARMa) #model diagnostics
library(MuMIn) #model diagnostics
library(glmmTMB) #mixed models
library(picante) #phylogenetic analyses
library(ape) #phylogenetic analyses
library(phytools) #phylogenetic analyses
library(phyloseq) #phylogenetic analyses
library(hciR) #for as_matrix function

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

#view sample sizes
metadata %>%
  filter(Source == "FIELD") %>%
  group_by(Sterilized) %>%
  summarize(sum = n())

#Assign taxonomies ####
#now we need to assign taxonomies so we can subset prey ASVs, and then we will
#be combining them by unique taxonomy for later analyses.

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

#Subset field spiders ####
#create a character vector of the names of the field spiders
fld <- as.character(metadata$sample[which(metadata$Source == "FIELD")])

#now subset the tables with taxonomic assignments for the columns
#that match these samples
d2_fld <- d2_id_tab %>% #214 when not separated into prey only
  dplyr::select(ASV, all_of(fld), ID_bold, ID_ncbi, ID_nondiet, taxonomy) %>%
  filter(taxonomy == "prey") #55

u3_fld <- u3_id_tab %>% #176 when not separated into prey only
  dplyr::select(ASV, all_of(fld), ID_bold, ID_ncbi, taxonomy) %>%
  filter(taxonomy == "prey") #42 prey


#Total prey ASV Count (with duplicates) ####
d2_tot <- d2_fld %>%
  gather(sample, reads, HEV65:HEV100) %>%
  left_join(metadata, by = "sample") %>%
  group_by(Sterilized) %>%
  summarise(counts = sum(reads > 0, na.rm = TRUE)) %>%
  mutate(pipeline = "d2")

u3_tot <- u3_fld %>%
  gather(sample, reads, HEV65:HEV100) %>%
  left_join(metadata, by = "sample") %>%
  group_by(Sterilized) %>%
  summarise(counts = sum(reads > 0, na.rm = TRUE)) %>%
  mutate(pipeline = "u3")

tot_ASVs <- u3_tot %>%
  bind_rows(d2_tot)

tot_ASVs <- tot_ASVs %>%
  group_by(Sterilized) %>%
  summarize(mean = mean(counts), se = sd(counts)/sqrt(2))

ggplot(tot_ASVs, aes(x = Sterilized, y = mean, fill = Sterilized)) +
  geom_bar(stat="identity") +theme_bw() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  labs(x = "Surface sterilization", y = "Number of ASVs", 
       title = "Number of prey ASVs detected comparing sterilized to unsterilized 
       (duplicate taxonomic matches not combined)") +
  scale_fill_manual(values = c("#7B8D65","#F29979")) +
  scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
  theme(legend.position = "none")

#based on raw ASV counts, non surface-sterilized spiders have greater ASV counts for
#both pipelines

#Total species richness of prey ASVs####
#(first, combine uniques) ####
#use d2_fld and u3_fld with uniques instructions below
#UNOISE3
#make the ID columns characters for the ifelse statements below
u3_fld$ID_bold <- as.character(u3_fld$ID_bold)
u3_fld$ID_ncbi <- as.character(u3_fld$ID_ncbi)

#ifelse statement that assigns each ASV to a specific taxonomy
#I want the BOLD ids to be the last option because these are to
#species in many cases where NCBI/Genbank blasted to order, family, or genus
#therefore, I'm going to say that if there is no BOLD id, give it the NCBI
#ID, otherwise, give it the BOLD ID.
IDs_u3 <- ifelse(is.na(u3_fld$ID_bold), u3_fld$ID_ncbi, 
                                               u3_fld$ID_bold)
IDs_u3
#make this a dataframe so we can attach it back to the community
#matrix for analyses at the sample level and later Jaccard dissimiliarity
IDs_u3 <- as.data.frame(IDs_u3)
IDs_u3$ASV <- u3_fld$ASV
IDs_u3 <- rename(IDs_u3, "unique_ID" = "IDs_u3")
unique_ID_u3 <- as.data.frame(unique(IDs_u3[c("unique_ID")]))
unique_ID_u3$sp_number <- seq.int(nrow(unique_ID_u3))

richness_u3 <- IDs_u3 %>%
  left_join(unique_ID_u3, by = "unique_ID") %>%
  left_join(u3_fld, by = "ASV") %>%
  gather(sample, reads, HEV65:HEV100) %>%
  left_join(metadata, by = "sample")

#richness_u3$sp_number <- as.factor(richness_u3$sp_number)
tot_rich_u3 <- richness_u3 %>%
  filter(reads > 0) %>%
  group_by(Sterilized) %>%
  summarize(rich = n_distinct(sp_number)) %>%
  mutate(pipeline = "u3")

#DADA2
#make the ID columns characters for the ifelse statements below
d2_fld$ID_bold <- as.character(d2_fld$ID_bold)
d2_fld$ID_ncbi <- as.character(d2_fld$ID_ncbi)
d2_fld$ID_nondiet <- as.character(d2_fld$ID_nondiet)
#ifelse statement that assigns each ASV to a specific taxonomy
#I want the BOLD ids to be the last option because these are to
#species in many cases where NCBI/Genbank blasted to order, family, or genus
#therefore, I'm going to say that if there is no BOLD id, give it the NCBI
#ID, otherwise, give it the BOLD ID.
IDs_d2 <- ifelse(is.na(d2_fld$ID_bold) & is.na(d2_fld$ID_ncbi), d2_fld$ID_nondiet, 
                 ifelse(is.na(d2_fld$ID_bold), d2_fld$ID_ncbi, 
                        d2_fld$ID_bold))
IDs_d2
#make this a dataframe so we can attach it back to the community
#matrix for analyses at the sample level and later Jaccard dissimiliarity
IDs_d2 <- as.data.frame(IDs_d2)
IDs_d2$ASV <- d2_fld$ASV
IDs_d2 <- rename(IDs_d2, "unique_ID" = "IDs_d2")
unique_ID_d2 <- unique(IDs_d2[c("unique_ID")]) 
unique_ID_d2$sp_number <- seq.int(nrow(unique_ID_d2))

richness_d2 <- IDs_d2 %>%
  left_join(unique_ID_d2, by = "unique_ID") %>%
  left_join(d2_fld, by = "ASV") %>%
  gather(sample, reads, HEV65:HEV100) %>%
  left_join(metadata, by = "sample")

#richness_u3$sp_number <- as.factor(richness_u3$sp_number)
tot_rich_d2 <- richness_d2 %>%
  filter(reads > 0) %>%
  group_by(Sterilized) %>%
  summarize(rich = n_distinct(sp_number)) %>%
  mutate(pipeline = "d2")

tot_rich <- tot_rich_u3 %>%
  bind_rows(tot_rich_d2)

tot_rich_sum <- tot_rich %>%
  group_by(Sterilized) %>%
  summarize(mean = mean(rich), se = sd(rich)/sqrt(2))

ggplot(tot_rich_sum, aes(x = Sterilized, y = mean, fill = Sterilized)) +
  geom_bar(stat = "identity") + theme_bw() +
  scale_fill_manual(values = c("#7B8D65","#F29979")) +
  geom_errorbar(aes(ymin= mean-se, ymax = mean+se), width = 0.2) +
  labs(x = "Surface sterilization", y = "Species richness", 
       title = "Prey richness detected comparing sterilized to unsterilized") +
  scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
  theme(legend.position = "none")

#again, surface sterilization seems to slightly reduce species richness. However,
#this is a bit misleading as there is one more sample in the NS group

#By sample richness ####
#we will be using richness_u3 and richness_d2 for this analysis
samples <- metadata %>%
  filter(Source == "FIELD") %>%
  dplyr::select(sample)

samp_rich_u3 <- richness_u3 %>%
  filter(reads > 0) %>%
  group_by(sample) %>%
  summarize(rich = n_distinct(sp_number)) %>%
  mutate(pipeline = "u3") %>%
  right_join(samples, by = "sample")

samp_rich_d2 <- richness_d2 %>%
  filter(reads > 0) %>%
  group_by(sample) %>%
  summarize(rich = n_distinct(sp_number)) %>%
  right_join(samples, by = "sample") %>%
  replace_na(list(rich = 0)) %>%
  mutate(pipeline = "d2")

#some samples have zero counts of prey species richness, so when we bind with metadata
#need to be sure to get those samples back, so we can set the NAs to zeros

samp_rich <- samp_rich_u3 %>%
  bind_rows(samp_rich_d2)

samp_rich <- samp_rich %>%
  left_join(metadata, by = "sample") 
  
#model: does by-sample richness change with sterilization?
#going to do a glmm with Poisson distribution of 
#richness ~ Sterilization + (1|Pipeline)

samp_rich$Ster <- ifelse(samp_rich$Sterilized == "NS", 0, 1)

rich_mod <- glmmTMB(rich ~ Ster + (1|pipeline),
                    data = samp_rich,
                    family = "genpois",
                    REML = FALSE)

rich_null <- glmmTMB(rich ~ 1 + (1|pipeline),
                     data = samp_rich,
                     family = "genpois",
                     REML = FALSE)

anova(rich_mod, rich_null) #no sig difference between these models, so the 
#null is probably best
AICc(rich_mod, rich_null)

rich_null <- glmmTMB(rich ~ 1 + (1|pipeline),
                     data = samp_rich,
                     family = "genpois")

simulationOutput_rich <- simulateResiduals(fittedModel = rich_null) 
rich_fit <- plot(simulationOutput_rich) #look okay
rich_zi <- testZeroInflation(simulationOutput_rich) #not zero inflated
rich_od <- testDispersion(simulationOutput_rich) 

plot(residuals(rich_null)) #look fine!

#richness by sample does not change with sterilization!! WOO

ggplot(samp_rich, aes(x = Sterilized, y = rich, fill = Sterilized)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Surface sterilization", y = "Per sample prey richness", 
       title = "Per sample prey richness by sterilization") +
  scale_fill_manual(values = c("#7B8D65","#F29979")) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) 

#Phylogenetic diversity of prey ASVs
#I'M NOT SURE that i need phylogenetic diversity as a main analysis - why would
#it matter for surface sterilization if phylo changed? if we get the same info 
#from doing a presence-absence multidimensional (PERMANOVA) analysis.
#Intermediate: Phylogenetic trees of prey ASVs ####
#this is a phylogenetic tree of the combined taxonomies from both NCBI and BOLD

#to create trees
#following this strategy:
#https://rdrr.io/cran/ape/man/as.phylo.formula.html
#data(carnivora)
#frm <- ~SuperFamily/Family/Genus/Species
#tr <- as.phylo(frm, data = carnivora)
#tr1 <- as.phylo(frm, data=carnivora, collapse = FALSE)
#plot(tr)
#Nnode(tr)
#Nnode(tr1)
#plot(tr1)
#plot(tr)
## compare with:
#Nnode(as.phylo(frm, data = carnivora, collapse = FALSE))

#SOME PD resources to start:
#ses.pd in the picante package:
#https://rdrr.io/rforge/picante/man/pd.html
#general phylogeny tricks here:
#http://www.phytools.org/eqg/Exercise_3.2/
#UNOISE3
u3_bold_pr <- u3_bold %>%
  filter(taxonomy == "prey")

u3_ncbi_pr <- u3_ncbi %>%
  filter(taxonomy == "prey")

u3_tax_pr <- u3_bold_pr %>%
  full_join(u3_ncbi_pr, by = "ASV") %>%
  dplyr::select(-BLAST_ID, -BLAST_Category, -BLAST_Level, -BLAST_Kingdom,
                -BLAST_Phylum, -BLAST_Class, -BLAST_Order, -BLAST_Family,
                -BLAST_Genus, -BLAST_Species, -taxonomy.x, -taxonomy.y, 
                -Level.x, -Level.y, -Category)

u3_tax_pr <- u3_tax_pr %>%
  mutate_if(is.factor, as.character) %>%
  glimpse()

u3_tax_pr$Kingdom <- ifelse(is.na(u3_tax_pr$Kingdom.x), u3_tax_pr$Kingdom.y, u3_tax_pr$Kingdom.x)
u3_tax_pr$Phylum <- ifelse(is.na(u3_tax_pr$Phylum.x), u3_tax_pr$Phylum.y, u3_tax_pr$Phylum.x)
u3_tax_pr$Class <- ifelse(is.na(u3_tax_pr$Class.x), u3_tax_pr$Class.y, u3_tax_pr$Class.x)
u3_tax_pr$Order <- ifelse(is.na(u3_tax_pr$Order.x), u3_tax_pr$Order.y, u3_tax_pr$Order.x)
u3_tax_pr$Family <- ifelse(is.na(u3_tax_pr$Family.x), u3_tax_pr$Family.y, u3_tax_pr$Family.x)
u3_tax_pr$Genus <- ifelse(is.na(u3_tax_pr$Genus.x), u3_tax_pr$Genus.y, u3_tax_pr$Genus.x)
u3_tax_pr$Species <- ifelse(is.na(u3_tax_pr$Species.x), u3_tax_pr$Species.y, u3_tax_pr$Species.x)

#DADA2
d2_bold_pr <- d2_bold %>%
  filter(taxonomy == "prey")

d2_ncbi_pr <- d2_ncbi %>%
  filter(taxonomy == "prey")

d2_tax_pr <- d2_bold_pr %>%
  full_join(d2_ncbi_pr, by = "ASV") %>%
  dplyr::select(-BLAST.ID, -blast_lev, -blast_k, -blast_p, -blast_c, 
                -blast_o, -taxonomy.x, -taxonomy.y, -Level.x, -Level.y)
d2_tax_pr <- d2_tax_pr %>%
  mutate_if(is.factor, as.character) %>%
  glimpse()

d2_tax_pr$Kingdom <- ifelse(is.na(d2_tax_pr$Kingdom.x), d2_tax_pr$Kingdom.y, d2_tax_pr$Kingdom.x)
d2_tax_pr$Phylum <- ifelse(is.na(d2_tax_pr$Phylum.x), d2_tax_pr$Phylum.y, d2_tax_pr$Phylum.x)
d2_tax_pr$Class <- ifelse(is.na(d2_tax_pr$Class.x), d2_tax_pr$Class.y, d2_tax_pr$Class.x)
d2_tax_pr$Order <- ifelse(is.na(d2_tax_pr$Order.x), d2_tax_pr$Order.y, d2_tax_pr$Order.x)
d2_tax_pr$Family <- ifelse(is.na(d2_tax_pr$Family.x), d2_tax_pr$Family.y, d2_tax_pr$Family.x)
d2_tax_pr$Genus <- ifelse(is.na(d2_tax_pr$Genus.x), d2_tax_pr$Genus.y, d2_tax_pr$Genus.x)
d2_tax_pr$Species <- ifelse(is.na(d2_tax_pr$Species.x), d2_tax_pr$Species.y, d2_tax_pr$Species.x)

#now we are going to just select variables we're interested in and then
#get only the distinct values based on these.
u3_tax_pr <- u3_tax_pr %>%
  mutate_if(is.character, as.factor) %>%
  dplyr::select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=T)

#now we are going to just select variables we're interested in and then
#get only the distinct values based on these.
d2_tax_pr <- d2_tax_pr %>%
  mutate_if(is.character, as.factor) %>%
  dplyr::select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=T)

all_ASV_dat <- u3_tax_pr %>%
  full_join(d2_tax_pr, by = c("Kingdom" = "Kingdom", "Phylum" = "Phylum", 
                              "Class" = "Class", "Order" = "Order",
                              "Family" = "Family", "Genus" = "Genus",
                              "Species" = "Species"
                              ))

#a couple were duplicated... so going to remove them
all_ASV_dat <- all_ASV_dat %>%
  filter(!ASV.y %in% c("ASV_100", "ASV_170")) %>%
  mutate_if(is.character, as.factor)


frm <- ~Kingdom/Phylum/Class/Order/Family/Genus/Species
tr <- as.phylo(frm, data = all_ASV_dat, collapse = FALSE)

#compute branch lengths for each tree
tr1 <- compute.brlen(tr, method = "Grafen")
tr1 <- di2multi(tr1) #this collapses multichotomies so that branch
#lengths aren't domianted by zeros
plot(tr1)

#for community data:
#need the dataframe with all ASVs with unique taxonomies assigned from both GenBank and BOLD
#so something like richness_u3 and richness_d2
#need to get a DF with pipeline, Sterilized, and species presence
PD_u3 <- richness_u3 %>%
  filter(reads > 0) %>%
  dplyr::select(unique_ID, ASV, sp_number, reads, Sterilized, sample) %>%
  group_by(Sterilized) %>%
  distinct(unique_ID) 

PD_ss_u3 <- PD_u3 %>%
  filter(Sterilized == "SS") %>%
  mutate(presence = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

#set rownames for calling PD function
row.names(PD_ss_u3) <- PD_ss_u3$unique_ID
PD_ss_u3 <- PD_ss_u3 %>%
  dplyr::select(presence) #select only presence from DF once rownames set
PD_ss_u3_t <- as.data.frame(t(PD_ss_u3)) #transpose this dataframe
pd_ss_u3 <- pd(PD_ss_u3_t, tr1, include.root=FALSE)
pd_ss_u3 <- pd_ss_u3 %>%
  mutate(pipeline = "u3", Sterilized = "SS")

PD_ns_u3 <- PD_u3 %>%
  filter(Sterilized == "NS") %>%
  mutate(presence = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

row.names(PD_ns_u3) <- PD_ns_u3$unique_ID
PD_ns_u3 <- PD_ns_u3 %>%
  dplyr::select(presence) #select only presence from DF once rownames set
PD_ns_u3_t <- as.data.frame(t(PD_ns_u3)) #transpose this dataframe
pd_ns_u3 <- pd(PD_ns_u3_t, tr1, include.root=FALSE)
pd_ns_u3 <- pd_ns_u3 %>%
  mutate(pipeline = "u3", Sterilized = "NS")

PD_d2 <- richness_d2 %>%
  filter(reads > 0) %>%
  dplyr::select(unique_ID, ASV, sp_number, reads, Sterilized, sample) %>%
  group_by(Sterilized) %>%
  distinct(unique_ID)

PD_ss_d2 <- PD_d2 %>%
  filter(Sterilized == "SS") %>%
  mutate(presence = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

row.names(PD_ss_d2) <- PD_ss_d2$unique_ID
PD_ss_d2 <- PD_ss_d2 %>%
  dplyr::select(presence) #select only presence from DF once rownames set
PD_ss_d2_t <- as.data.frame(t(PD_ss_d2)) #transpose this dataframe
pd_ss_d2 <- pd(PD_ss_d2_t, tr1, include.root=FALSE)
pd_ss_d2 <- pd_ss_d2 %>%
  mutate(pipeline = "d2", Sterilized = "SS")

PD_ns_d2 <- PD_d2 %>%
  filter(Sterilized == "NS") %>%
  mutate(presence = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

row.names(PD_ns_d2) <- PD_ns_d2$unique_ID
PD_ns_d2 <- PD_ns_d2 %>%
  dplyr::select(presence) #select only presence from DF once rownames set
PD_ns_d2_t <- as.data.frame(t(PD_ns_d2)) #transpose this dataframe
pd_ns_d2 <- pd(PD_ns_d2_t, tr1, include.root=FALSE)
pd_ns_d2 <- pd_ns_d2 %>%
  mutate(pipeline = "d2", Sterilized = "NS")

pd_all <- pd_ss_u3 %>%
  bind_rows(pd_ns_u3) %>%
  bind_rows(pd_ss_d2) %>%
  bind_rows(pd_ns_d2) 

pd_all_sum <- pd_all %>%
  group_by(Sterilized) %>%
  summarize(mean = mean(PD), se = sd(PD)/sqrt(2))

ggplot(pd_all_sum, aes(x = Sterilized, y = mean, fill = Sterilized)) +
  geom_bar(stat="identity") +theme_bw() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  labs(x = "Surface sterilization", y = "Prey phylogenetic diversity (Faith's PD)", 
       title = "Phylogenetic diversity of prey") +
  scale_fill_manual(values = c("#7B8D65","#F29979")) +
  scale_x_discrete(labels=c("NS" = "Not S. Sterilized", "SS" = "Surface Sterilized")) +
  theme(legend.position = "none")

#phylogenetic diversity looks higher for non-sterilized. Again, keep in mind that
#this is one extra sample as well. 

#not working Phylogeny graph Section ####
#for phyloseq tree graphing, I need: 
#1. an OTU table which should just be a dataframe with 
#PD_d2 and PD_u3 with each row as a pipeline-sterilization pair
PD_d2_ns <- PD_d2 %>%
  filter(Sterilized == "NS") %>%
  mutate(NS_d2 = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

PD_d2_ss <- PD_d2 %>%
  filter(Sterilized == "SS") %>%
  mutate(SS_d2 = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

PD_u3_ns <- PD_u3 %>%
  filter(Sterilized == "NS") %>%
  mutate(NS_u3 = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

PD_u3_ss <- PD_u3 %>%
  filter(Sterilized == "SS") %>%
  mutate(SS_u3 = 1) %>%
  ungroup() %>%
  dplyr::select(-Sterilized)

otu_tab <- PD_u3_ns %>%
  full_join(PD_u3_ss, by = "unique_ID") %>%
  full_join(PD_d2_ns, by = "unique_ID") %>%
  full_join(PD_d2_ss, by = "unique_ID") %>%
  replace(., is.na(.), 0) 
rownames(otu_tab) <- otu_tab$unique_ID
otu_tab <- otu_tab %>%
  dplyr::select(-unique_ID)

#2. taxonomies, which can be d2_tax_pr bound to u3_tax_pr
taxonomies_ph <- d2_tax_pr %>%
  full_join(u3_tax_pr, by = c("Kingdom", "Phylum", "Class", "Order", "Family",
                              "Genus", "Species")) %>%
  filter(!ASV.x %in% c("ASV_100", "ASV_170")) %>%
  dplyr::select(-ASV.x, -ASV.y) 

taxonomies_ph <- taxonomies_ph%>%
  mutate_all(na_if,"")

taxonomies_ph$Species <- ifelse(is.na(taxonomies_ph$Order), taxonomies_ph$Class,
                                ifelse(is.na(taxonomies_ph$Family), taxonomies_ph$Order,
                                             ifelse(is.na(taxonomies_ph$Genus), taxonomies_ph$Family,
                                                          ifelse(is.na(taxonomies_ph$Species), taxonomies_ph$Genus,
                                                                 taxonomies_ph$Species))))

frm <- ~Kingdom/Phylum/Class/Order/Family/Genus/Species
taxonomies_ph <- taxonomies_ph %>%
  mutate_if(is.character, as.factor) %>%
  glimpse()
  
tr2 <- as.phylo(frm, data = taxonomies_ph, collapse = FALSE)
plot(tr2)
taxonomies_ph <- as.matrix(taxonomies_ph, rownames.value = c(nrow(taxonomies_ph)))


#3. taxonomic tree, which in this case is tr1
#trouble with the phyloseq object, and need to figure out
#where this error is coming from:
#Error in `taxa_names<-`(`*tmp*`, value = gsub("\"", "", taxa_names(x),  : 
#taxa_names<-: You are attempting to assign duplicated taxa_names
#https://github.com/joey711/phyloseq/issues/785
#from phyloseq package
ASVs <- otu_table(otu_tab, taxa_are_rows = TRUE)
taxa <- tax_table(taxonomies_ph)



#compute branch lengths for each tree
#tr1 <- compute.brlen(tr, method = "Grafen")
#tr1 <- di2multi(tr1) #this collapses multichotomies so that branch
#lengths aren't domianted by zeros
#plot(tr1)

physeq_graph <-  phyloseq(ASVs, taxa, tr1)

family_phy_graph <- plot_tree(physeq_family, nodelabf = nodeplotblank, plot.margin=0.6, label.tips = "taxa_names", shape = "samples")
family_phy_graph


#end of not working phyl. section####

#### end here Apr1 ##### 
#OUtline to flesh out####

  #e.	Prey ASV read abundance (rarefied first)
  #i. per sample (stats)
  