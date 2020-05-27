#Presence-based community dissimilarity
#Ana Miller-ter Kuile
#February 2020

#load packages ####
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

#This code looks at the following questions for each community matrix,
#broken up into sterilized, not sterilized, field, and lab:
#1. is community different (presence/absence) between sterilized and 
#unsterilized individuals? (PERMANOVA)
#1a. follow-up ecological questions here? (i.e. who are community differences,
#are they "prey", "pred", "non-diet", or "no hit"? are there different numbers
#of prey ASVs when the animal is sterilized first?)

#Thoughts March 3###
#after some issues with treating every ASV as a unique ID, I have decided I am
#going to combine things by unique taxonomies - so if it matched to the same
#thing in NCBI/BOLD, combine those reads for subsequent analyses. Furthermore, while
#I'm stickign with that taxonomic level for things that *could* be prey, I'm going
#to go with a higher grouping variable for things that are not prey, or unknown, so
#anything with an "unknown" or "non-diet" or "no hit" taxonomy will be lumped into
#one thing. For the purposes of these analyses, the important question is whether these
#things are present more often for unsterilized vs. surface sterilized individuals, and
#their actual specific identities are less important. For the paper, I will likely
#include soem supplementary data indicating what taxonomies fit in these categories.

#note on NMDS stress:
#if stress is high, reposition the points in 2 dimensions in the direction of 
#decreasing stress, and repeat until stress is below some threshold.
#**A good rule of thumb: stress < 0.05 provides an excellent 
#representation in reduced dimensions, < 0.1 is great, < 0.2 is good/ok, 
#and stress < 0.3 provides a poor representation.** 
#To reiterate: high stress is bad, low stress is good!

#Data for all analyses ####
samples <- read.csv(here("data", "Sample_Metadata.csv"))

#unoise: ASV tables and metadata####

#divide into "lab" and "field" spider datasets using samples dataframe
#first transpose dataframe so samples are on the rows
u3_comm <- read.csv(here("data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = "\t")

#rename columns for simplicity
colnames(u3_comm) <- sapply(str_split(colnames(u3_comm), "S"), function(x){return(x[[1]])})
u3_comm <- rename(u3_comm, "ASV" = "X.OTU.ID")

#set row names to ASV labels
rownames(u3_comm) <- u3_comm$ASV
#select samples minus controls and the ASV column
u3_comm <- u3_comm %>%
  dplyr::select(-CL1, -CL4, -QC1, -ASV, -NEG)

#transpose the matrix (I think) so that ASVs on columns, samples on rows
u3_comm_t <- as.data.frame(t(u3_comm))
#set all na to 0, if there are NA values
#d2_uc_comm_t[is.na(d2_uc_comm_t)] <- 0
#d2_uc_comm[is.na(d2_uc_comm)] <- 0
#remove any ASVs that have zeros across all samples (probably from removing control)
u3_comm <- u3_comm[rowSums(u3_comm) != 0,] #3 ASVs disappeared

u3_comm_t$sample <- rownames(u3_comm_t)

#combine by sample name with sample dataframe
u3_nmds_all <- u3_comm_t %>%
  left_join(samples, by = "sample")

#unoise: Field tables and metadata####
f_u3_nmds <- u3_nmds_all %>%
  dplyr::filter(Source == "FIELD")

#write this to a file so we can import it later quickly into the abundance based 
#community dissimilarity code
write_csv(f_u3_nmds, here("data", "fld_u3uc_nmds.csv"))


#and then select just the numeric values in this dataframe
rownames(f_u3_nmds) <- f_u3_nmds$sample
fld_u3_nmds <- f_u3_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#any zeros will mess up the analysis, so we need to check this first
rowSums(fld_u3_nmds)
colSums(fld_u3_nmds) #there are quite a few ASVs not present in this subset of data

#delete those zero columns for analyses
fld_u3_nmds <- fld_u3_nmds[which(colSums(fld_u3_nmds)> 0)] #176 down to 145

#metadata for field NMDS
#need to look at sample metadata DF to set the meta df for both field and lab
u3_meta_field <- f_u3_nmds %>%
  dplyr::select(sample, Method, Island, Year, ID, Date.Collected,
                dada2_sample, unoise_sample, Size, Extract..Date, Sterilized, Source)

#again, going to write this so it is easier to import later
write_csv(u3_meta_field, here("data", "fld_u3_nmds_meta.csv"))

#name the taxonomic assignment dataframe from the pipeline code set
u3_taxa <- read_csv(here("data", "u3_uc_tax_ass.csv"))

#unoise: uniques for community analyses####
#I am going to run a GLMER PERMANOVA with unique species combined, rather than by ASV
#to see if it fixes my singular model problems and simplifies that data analyses.

#need to bond taxonomy with DF with ASVs on rows and separated into field and lab
#transpose the community matrix:
u3_f_glmm <- as.data.frame(t(fld_u3_nmds))
#set an ASV column for combining with the taxonomic assignment dataframe
u3_f_glmm$ASV <- rownames(u3_f_glmm)

#combine our matrix of samples and ASVs with the corresponding taxonomies
u3_f_glmm <- u3_f_glmm %>%
  left_join(u3_taxa, by = "ASV")

#creates a new vector of the unique species IDs for these data
#i'm going to categorize all "non-diet in one in the hopes of simplifying
#the model below, and then treating all diet as unique potential "diet" items
#difference being - I *know* the fungi are not diet, but I don't know anything
#that could be diet is or isn't actually diet or contamination
uniques_f_u3 <- ifelse(u3_f_glmm$taxonomy == "non-diet", u3_f_glmm$taxonomy,
                         ifelse(u3_f_glmm$taxonomy != "non-diet" & !is.na(u3_f_glmm$ID_bold), 
                                u3_f_glmm$ID_bold,
                                ifelse(u3_f_glmm$taxonomy != "non-diet" & !is.na(u3_f_glmm$NCBI_ID),
                                       u3_f_glmm$NCBI_ID,
                                ifelse(u3_f_glmm$taxonomy != "non-diet" & !is.na(u3_f_glmm$ID_ncbi),
                                       u3_f_glmm$ID_ncbi, "no_hit"))))


# 146 uniques
#makes all predators map to Sparassidae, resetting some of the weirdness b/w BOLD and NCBI
uniques_f_u3[uniques_f_u3 == "Heteropoda venatoria"] <- "Sparassidae"

#makes these unique taxonomies a dataframe                 
uniques_f_u3 <- as.data.frame(uniques_f_u3)  
#add X variable for joining with the sample dataframe
uniques_f_u3$X <- u3_f_glmm$X1
uniques_f_u3 <- rename(uniques_f_u3, "uniques" = "uniques_f_u3")
#join with sample dataframe
u3_f_glmm <- u3_f_glmm %>%
  left_join(uniques_f_u3, by = c("X1" = "X"))

#create a unique ID for each species for the Permanova GLMER
unique_field_u3 <- unique(uniques_f_u3[c("uniques")]) #26 unique taxonomies
#give these numbers here:
unique_field_u3$species_ID <- seq.int(nrow(unique_field_u3))

#join this with our sample dataframe
u3_f_glmm <- u3_f_glmm %>%
  left_join(unique_field_u3, by = "uniques")

#make this a long dataframe for GLMer
u3_field_glmm <- u3_f_glmm %>%
  gather(sample, reads, HEV87:HEV83)

#combines reads from each sample from the same unique "species" and adds up the reads
u3_field_glmm <- u3_field_glmm %>%
  group_by(sample, uniques) %>% 
  summarise(reads=sum(reads)) %>%
  ungroup()

u_rare_u3 <- u3_field_glmm %>%
  group_by(uniques) %>%
  tally(reads > 0)
#there are about 10 uniques that are only in one sample 

u_ct_u3 <- u3_field_glmm %>%
  group_by(sample) %>%
  tally(reads > 0)
#every sample has 3 or more unique IDs

#unoise: field GLMER PERMANOVA ####
#add sample meta-data from sample dataframe
u3_field_glmm <- u3_field_glmm %>%
  inner_join(samples, by = "sample")

#set a field for presence/absence for analysis with PERMANOVA
u3_field_glmm$Presence <- ifelse(u3_field_glmm$reads >0,1,0)

#adds the unique sample ID value to this dataframe so we can call it as a 
#random effect below
u3_field_glmm <- u3_field_glmm %>%
  left_join(unique_field_u3, by = "uniques")

#set species ID to a factor
u3_field_glmm$species_ID <- as.factor(u3_field_glmm$species_ID)
u3_field_glmm$Ster <- ifelse(u3_field_glmm$Sterilized == "NS", 0, 1)
#model that says that the slope of the response to sterilization can vary
#randomly by Species_ID
u3_field_model <- glmer(Presence ~ Ster + (1+Ster|species_ID), 
                     data = u3_field_glmm,
                     family = "binomial")
#check if it is singular based on warning from running model

#looks at the summary of this model
summary(u3_field_model)
plot(residuals(u3_field_model))

#dotplot(ranef(field_model,
#              condVar = TRUE))

#this is the null model with just a random intercept by species ID
u3_field_no_ster <- glmer(Presence ~ 1 + (1|species_ID), 
                       data = u3_field_glmm, 
                       family = "binomial")

#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(u3_field_model, u3_field_no_ster)

#unoise: Field NMDS ####
#i'm going to do an NMDS of this to visualize the community
#using d2uc_field_glmm_uniques
#select columns of interest and then spread back to a community
#matrix for NMDS with the samples as row names
perm_field <- u3_field_glmm %>%
  dplyr::select(sample, uniques, reads) %>%
  spread(uniques, reads) %>%
  column_to_rownames("sample")

#run NMDS based on jaccard dissimilarity
u3_nmds_field <- metaMDS(perm_field, distance = "jaccard",binary =TRUE, trymax=1000, k=2) 
#evaluate the fit of this NMDS
stressplot(u3_nmds_field) 
u3_nmds_field$stress #0.1814212, not *great*, but not bad
#create a dataframe for plotting
u3_nmds_field_df <- data.frame(MDS1=u3_nmds_field$points[,1], MDS2=u3_nmds_field$points[,2])
#combining with metadata for plotting
u3_nmds_field_meta <- cbind(u3_meta_field, u3_nmds_field_df)

#plots NMDS with elipses
u3_f_nmds_graph <- ggplot(u3_nmds_field_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    labs(title = "UNOISE3 Field Jaccard Dissimilarity") +
    theme_bw() +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#I combined all the "non-diet" ASVs into categories of "non-diet" and "unknown" 
#because they were messing up the model, which is still singular... but I'm still
#okay with the way this model is run I think. SO I think this is the level at which
#I will be running this. YAY

#unoise: Lab tables and metadata ####
#subset lab collected samples from these  
l_u3_nmds <- u3_nmds_all %>%
  dplyr::filter(Source == "LAB")

#again write for later import
write_csv(l_u3_nmds, here("data", "lab_u3uc_nmds.csv"))

#set rownames to sample variable
rownames(l_u3_nmds) <- l_u3_nmds$sample

#and then select just the numeric values in this dataframe
lab_u3_nmds <- l_u3_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#any zeros will mess up the analysis, so we need to check this first
#rowSums(lab_u3_nmds)
#colSums(lab_u3_nmds) #there are quite a few ASVs not present in this subset of data

#delete those zero columns for analyses
lab_u3_nmds <- lab_u3_nmds[which(colSums(lab_u3_nmds)> 0)] #176 to 75

#metadata for lab NMDS
#need to look at sample metadata DF to set the meta df for both field and lab
u3_meta_lab <- l_u3_nmds %>%
  dplyr::select(sample, Method, Island, Year, ID, Date.Collected,
                dada2_sample, unoise_sample, Size, Extract..Date, Sterilized, Source)
#and write for later import
write_csv(u3_meta_lab, here("data", "lab_u3uc_nmds_meta.csv"))

#unoise: uniques for community analyses####
#I am going to run a GLMER PERMANOVA with unique species combined, rather than by ASV
#to see if it fixes my singular model problems and simplifies that data analyses.

#need to bond taxonomy with DF with ASVs on rows and separated into field and lab
#transpose the community matrix:
u3_l_glmm <- as.data.frame(t(lab_u3_nmds))
#set an ASV column for combining with the taxonomic assignment dataframe
u3_l_glmm$ASV <- rownames(u3_l_glmm)

#combine our matrix of samples and ASVs with the corresponding taxonomies
u3_l_glmm <- u3_l_glmm %>%
  left_join(u3_taxa, by = "ASV")

#creates a new vector of the unique species IDs for these data
u3_uniques_lab <- ifelse(u3_l_glmm$taxonomy == "non-diet", u3_l_glmm$taxonomy,
                             ifelse(u3_l_glmm$taxonomy != "non-diet" & !is.na(u3_l_glmm$ID_bold), u3_l_glmm$ID_bold,
                                    ifelse(u3_l_glmm$taxonomy != "non-diet" & !is.na(u3_l_glmm$NCBI_ID), u3_l_glmm$NCBI_ID,
                                    ifelse(u3_l_glmm$taxonomy != "non-diet" & !is.na(u3_l_glmm$ID_ncbi), u3_l_glmm$ID_ncbi,
                                           "no_hit"))))

#u3_uniques_lab

#makes all predators map to Sparassidae, resetting some of the weirdness b/w BOLD and NCBI
u3_uniques_lab[u3_uniques_lab == "Heteropoda venatoria"] <- "Sparassidae"

#makes these unique taxonomies a dataframe  
#77 unique taxonomies               
u3_uniques_lab <- as.data.frame(u3_uniques_lab)  
#add X variable for joining with the sample dataframe
u3_uniques_lab$X <- u3_l_glmm$X1

u3_uniques_lab <- rename(u3_uniques_lab, "uniques" = "u3_uniques_lab")
#join with sample dataframe
u3_l_glmm <- u3_l_glmm %>%
  left_join(u3_uniques_lab, by =c("X1" = "X"))

#create a unique ID for each species for the Permanova GLMER
u3_unique_lab <- unique(u3_uniques_lab[c("uniques")]) #15 unique taxonomies
#give these numbers here:
u3_unique_lab$species_ID <- seq.int(nrow(u3_unique_lab))

#join this with our sample dataframe
u3_l_glmm <- u3_l_glmm %>%
  left_join(u3_unique_lab, by = "uniques")

#make this a long dataframe for GLMer
u3_lab_glmm <- u3_l_glmm %>%
  gather(sample, reads, HEV15:HEV29)

#combines reads from each sample from the same unique "species" and adds up the reads
u3_lab_glmm <- u3_lab_glmm %>%
  group_by(sample, uniques) %>% 
  summarise(reads=sum(reads)) %>%
  ungroup()

#unoise: lab GLMER PERMANOVA ####

#add sample meta-data from sample dataframe
u3_lab_glmm <- u3_lab_glmm %>%
  inner_join(samples, by = "sample")

#set a field for presence/absence for analysis with PERMANOVA
u3_lab_glmm$Presence <- ifelse(u3_lab_glmm$reads >0,1,0)

#adds the unique sample ID value to this dataframe so we can call it as a 
#random effect below
u3_lab_glmm <- u3_lab_glmm %>%
  left_join(u3_unique_lab, by = "uniques")

#set species ID to a factor
u3_lab_glmm$species_ID <- as.factor(u3_lab_glmm$species_ID)
u3_lab_glmm$Ster <- ifelse(u3_lab_glmm$Sterilized == "NS", 0, 1)
#model that says that the slope of the response to sterilization can vary
#randomly by Species_ID
u3_lab_model <- glmer(Presence ~ Ster + (1+Ster|species_ID), 
                   data = u3_lab_glmm, 
                   glmerControl(optimizer = "bobyqa"),
                   family = binomial)

#looks at the summary of this model
summary(u3_lab_model)

#dotplots with species which differ between communities being shown in the Ster
#dotplot without error bars crossing zero
#dotplot(ranef(u3_lab_model,
#              condVar =T))


#this is the null model with just a random intercept by species ID
u3_lab_no_ster <- glmer(Presence ~ 1 + (1|species_ID), 
                     data = u3_lab_glmm, 
                     glmerControl(optimizer = "bobyqa"),
                     family = binomial)

#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(u3_lab_model, u3_lab_no_ster)

#unoise: Lab NMDS####
#i'm going to do an NMDS of this to visualize the community
#using u3_lab_glmm
#select columns of interest and then spread back to a community
#matrix for NMDS with the samples as row names
u3_perm_lab <- u3_lab_glmm %>%
  dplyr::select(sample, uniques, reads) %>%
  spread(uniques, reads) %>%
  column_to_rownames("sample")

#run NMDS based on jaccard dissimilarity
u3_nmds_lab <- metaMDS(u3_perm_lab, distance = "jaccard", trymax = 1000,binary=TRUE, k=2) 
#evaluate the fit of this NMDS
stressplot(u3_nmds_lab) 
u3_nmds_lab$stress #0.1009513 looks good!
#create a dataframe for plotting
u3_nmds_lab_df <- data.frame(MDS1=u3_nmds_lab$points[,1], MDS2=u3_nmds_lab$points[,2])
#combining with metadata for plotting
u3_nmds_lab_meta <- cbind(u3_meta_lab, u3_nmds_lab_df)

#plots NMDS with elipses
u3_l_nmds_graph <- ggplot(u3_nmds_lab_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    labs(title = "UNOISE3 Lab Jaccard Dissimilarity") +
    theme_bw() +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#adonis Permanova presence/absence ####
#for permanova with adonis
#with predator ASVs
adonis(fld_u3_nmds ~ Sterilized, data = u3_meta_field, dist = "jaccard", binary=TRUE)

adonis(lab_u3_nmds ~ Sterilized, data = u3_meta_lab, dist = "jaccard", binary = TRUE)

#dada2: ASV tables and metadata####

#divide into "lab" and "field" spider datasets using samples dataframe
#first transpose dataframe so samples are on the rows
d2_comm <- read.csv(here("data", "ASV_tables", "dada2_uc_asv_tab.tsv"), sep = "\t")

#rename columns for simplicity
colnames(d2_comm) <- sapply(str_split(colnames(d2_comm), "_"), function(x){return(x[[1]])})
colnames(d2_comm) <- str_remove(colnames(d2_comm), "\\.")
d2_comm <- rename(d2_comm, "ASV" = "X")

#set row names to ASV labels
rownames(d2_comm) <- d2_comm$ASV
#select samples minus controls and the ASV column
d2_comm <- d2_comm %>%
  dplyr::select(-CL12, -CL42, -QC1, -ASV)

#transpose the matrix (I think) so that ASVs on columns, samples on rows
d2_comm_t <- as.data.frame(t(d2_comm))
#set all na to 0, if there are NA values
#d2_uc_comm_t[is.na(d2_uc_comm_t)] <- 0
#d2_uc_comm[is.na(d2_uc_comm)] <- 0
#remove any ASVs that have zeros across all samples (probably from removing control)
d2_comm <- d2_comm[rowSums(d2_comm) != 0,] #3 ASVs disappeared

d2_comm_t$sample <- rownames(d2_comm_t)

#combine by sample name with sample dataframe
d2_nmds_all <- d2_comm_t %>%
  left_join(samples, by = "sample")

#dada2: Field tables and metadata####
f_d2_nmds <- d2_nmds_all %>%
  dplyr::filter(Source == "FIELD")

#write this to a file so we can import it later quickly into the abundance based 
#community dissimilarity code
write_csv(f_d2_nmds, here("data", "fld_d2uc_nmds.csv"))

#set rownames to sample variable
rownames(f_d2_nmds) <- f_d2_nmds$sample

#and then select just the numeric values in this dataframe
fld_d2_nmds <- f_d2_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#any zeros will mess up the analysis, so we need to check this first
#rowSums(fld_d2_nmds)
#colSums(fld_d2_nmds) #there are quite a few ASVs not present in this subset of data

#delete those zero columns for analyses
fld_d2_nmds <- fld_d2_nmds[which(colSums(fld_d2_nmds)> 0)] #214 to 145

#metadata for field NMDS
#need to look at sample metadata DF to set the meta df for both field and lab
d2_meta_field <- f_d2_nmds %>%
  dplyr::select(sample, Method, Island, Year, ID, Date.Collected,
                dada2_sample, unoise_sample, Size, Extract..Date, Sterilized, Source)

#again, going to write this so it is easier to import later
write_csv(d2_meta_field, here("data", "fld_d2_nmds_meta.csv"))

#name the taxonomic assignment dataframe from the pipeline code set
d2_taxa <- read_csv(here("data", "d2_uc_tax_ass.csv"))

#dada2: uniques for community analyses####
#I am going to run a GLMER PERMANOVA with unique species combined, rather than by ASV
#to see if it fixes my singular model problems and simplifies that data analyses.

#need to bond taxonomy with DF with ASVs on rows and separated into field and lab
#transpose the community matrix:
d2_f_glmm <- as.data.frame(t(fld_d2_nmds))
#set an ASV column for combining with the taxonomic assignment dataframe
d2_f_glmm$ASV <- rownames(d2_f_glmm)

#combine our matrix of samples and ASVs with the corresponding taxonomies
d2_f_glmm <- d2_f_glmm %>%
  left_join(d2_taxa, by = "ASV")

#creates a new vector of the unique species IDs for these data
#i'm going to categorize all "non-diet in one in the hopes of simplifying
#the model below, and then treating all diet as unique potential "diet" items
#difference being - I *know* the fungi are not diet, but I don't know anything
#that could be diet is or isn't actually diet or contamination
uniques_f_d2 <- ifelse(d2_f_glmm$taxonomy == "non-diet", d2_f_glmm$taxonomy,
                              ifelse(d2_f_glmm$taxonomy != "non-diet" & !is.na(d2_f_glmm$ID_bold), 
                                     d2_f_glmm$ID_bold,
                                     ifelse(d2_f_glmm$taxonomy != "non-diet" & !is.na(d2_f_glmm$NCBI_ID), d2_f_glmm$NCBI_ID,
                                     ifelse(d2_f_glmm$taxonomy != "non-diet" &  !is.na(d2_f_glmm$ID_ncbi),
                                            d2_f_glmm$ID_ncbi, "no_hit"))))


# 154 uniques
#makes all predators map to Sparassidae, resetting some of the weirdness b/w BOLD and NCBI
uniques_f_d2[uniques_f_d2 == "Heteropoda venatoria"] <- "Sparassidae"

#makes these unique taxonomies a dataframe                 
uniques_f_d2 <- as.data.frame(uniques_f_d2)  
#add X variable for joining with the sample dataframe
uniques_f_d2$X <- d2_f_glmm$X1
uniques_f_d2 <- rename(uniques_f_d2, "uniques" = "uniques_f_d2")
#join with sample dataframe
d2_f_glmm <- d2_f_glmm %>%
  left_join(uniques_f_d2, by = c("X1" = "X"))

#create a unique ID for each species for the Permanova GLMER
unique_field_d2 <- unique(uniques_f_d2[c("uniques")]) #30 unique taxonomies
#give these numbers here:
unique_field_d2$species_ID <- seq.int(nrow(unique_field_d2))

#join this with our sample dataframe
d2_f_glmm <- d2_f_glmm %>%
  left_join(unique_field_d2, by = "uniques")

#make this a long dataframe for GLMer
d2_field_glmm <- d2_f_glmm %>%
  gather(sample, reads, HEV100:HEV89)

#combines reads from each sample from the same unique "species" and adds up the reads
d2_field_glmm <- d2_field_glmm %>%
  group_by(sample, uniques) %>% 
  summarise(reads=sum(reads)) %>%
  ungroup()

u_rare_d2 <- d2_field_glmm %>%
  group_by(uniques) %>%
  tally(reads > 0)
#there are about 17 uniques that are only in one sample 

u_ct_d2 <- d2_field_glmm %>%
  group_by(sample) %>%
  tally(reads > 0)
#a couple samples only have 1 unique ID, and several others have only 2

#dada2: field GLMER PERMANOVA ####
#add sample meta-data from sample dataframe
d2_field_glmm <- d2_field_glmm %>%
  inner_join(samples, by = "sample")

#set a field for presence/absence for analysis with PERMANOVA
d2_field_glmm$Presence <- ifelse(d2_field_glmm$reads >0,1,0)

#adds the unique sample ID value to this dataframe so we can call it as a 
#random effect below
d2_field_glmm <- d2_field_glmm %>%
  left_join(unique_field_d2, by = "uniques")

#set species ID to a factor
d2_field_glmm$species_ID <- as.factor(d2_field_glmm$species_ID)
d2_field_glmm$Ster <- ifelse(d2_field_glmm$Sterilized == "NS", 0, 1)
#model that says that the slope of the response to sterilization can vary
#randomly by Species_ID
d2_field_model <- glmer(Presence ~ Ster + (1+Ster|species_ID), 
                     data = d2_field_glmm,
                     family = "binomial")
#check if it is singular based on warning from running model
#isSingular(d2_field_model) #isn't singular... so dunno
#looks at the summary of this model
summary(d2_field_model)
plot(residuals(d2_field_model))

#this is the null model with just a random intercept by species ID
d2_field_no_ster <- glmer(Presence ~ 1 + (1|species_ID), 
                       data = d2_field_glmm, 
                       family = "binomial")

#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(d2_field_model, d2_field_no_ster)

#this is a dot graph to visualize sample variation
#dotplot(ranef(d2_field_model, condVar=T))

#dada2: Field NMDS ####
#i'm going to do an NMDS of this to visualize the community
#using d2uc_field_glmm_uniques
#select columns of interest and then spread back to a community
#matrix for NMDS with the samples as row names
perm_field <- d2_field_glmm %>%
  dplyr::select(sample, uniques, reads) %>%
  spread(uniques, reads) %>%
  column_to_rownames("sample")

#run NMDS based on jaccard dissimilarity
d2_nmds_field <- metaMDS(perm_field, distance = "jaccard",binary =TRUE, trymax=1000, k=2) 
#evaluate the fit of this NMDS
stressplot(d2_nmds_field) 
d2_nmds_field$stress #0.1620101
#create a dataframe for plotting
d2_nmds_field_df <- data.frame(MDS1=d2_nmds_field$points[,1], MDS2=d2_nmds_field$points[,2])
#combining with metadata for plotting
d2_nmds_field_meta <- cbind(d2_meta_field, d2_nmds_field_df)

#plots NMDS with elipses
d2_f_nmds_graph <- ggplot(d2_nmds_field_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    labs(title = "DADA2 Field Jaccard Dissimilarity") +
    theme_bw() +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#I combined all the "non-diet" ASVs into categories of "non-diet" and "unknown" 
#because they were messing up the model, which is still singular... but I'm still
#okay with the way this model is run I think. SO I think this is the level at which
#I will be running this. YAY

#dada2: Lab tables and metadata ####
#subset lab collected samples from these  
l_d2_nmds <- d2_nmds_all %>%
  dplyr::filter(Source == "LAB")

#again write for later import
write_csv(l_d2_nmds, here("data", "lab_d2uc_nmds.csv"))
#set rownames to sample variable
rownames(l_d2_nmds) <- l_d2_nmds$sample

#and then select just the numeric values in this dataframe
lab_d2_nmds <- l_d2_nmds %>%
  dplyr::select(-sample, -Method, -Island, -Year, -ID, -Date.Collected,
                -dada2_sample, -unoise_sample, -Size, -Extract..Date, -Sterilized, -Source)

#any zeros will mess up the analysis, so we need to check this first
#rowSums(lab_d2_nmds)
#colSums(lab_d2_nmds) #there are quite a few ASVs not present in this subset of data

#delete those zero columns for analyses
lab_d2_nmds <- lab_d2_nmds[which(colSums(lab_d2_nmds)> 0)] #214 to 82

#metadata for lab NMDS
#need to look at sample metadata DF to set the meta df for both field and lab
d2_meta_lab <- l_d2_nmds %>%
  dplyr::select(sample, Method, Island, Year, ID, Date.Collected,
                dada2_sample, unoise_sample, Size, Extract..Date, Sterilized, Source)
#and write for later import
write_csv(d2_meta_lab, here("data", "lab_d2uc_nmds_meta.csv"))

#dada2: uniques for community analyses####
#I am going to run a GLMER PERMANOVA with unique species combined, rather than by ASV
#to see if it fixes my singular model problems and simplifies that data analyses.

#need to bond taxonomy with DF with ASVs on rows and separated into field and lab
#transpose the community matrix:
d2_l_glmm <- as.data.frame(t(lab_d2_nmds))
#set an ASV column for combining with the taxonomic assignment dataframe
d2_l_glmm$ASV <- rownames(d2_l_glmm)

#combine our matrix of samples and ASVs with the corresponding taxonomies
d2_l_glmm <- d2_l_glmm %>%
  left_join(d2_taxa, by = "ASV")

#creates a new vector of the unique species IDs for these data
d2_uniques_lab <- ifelse(d2_l_glmm$taxonomy == "non-diet", d2_l_glmm$taxonomy,
                                ifelse(d2_l_glmm$taxonomy != "non-diet" & !is.na(d2_l_glmm$ID_bold), d2_l_glmm$ID_bold, 
                                       ifelse(d2_l_glmm$taxonomy != "non-diet" & !is.na(d2_l_glmm$NCBI_ID), d2_l_glmm$NCBI_ID,
                                       ifelse(d2_l_glmm$taxonomy != "non-diet" & !is.na(d2_l_glmm$ID_ncbi), d2_l_glmm$ID_ncbi,
                                              "no_hit"))))

#d2_uniques_lab

#makes all predators map to Sparassidae, resetting some of the weirdness b/w BOLD and NCBI
d2_uniques_lab[d2_uniques_lab == "Heteropoda venatoria"] <- "Sparassidae"

#makes these unique taxonomies a dataframe                 
d2_uniques_lab <- as.data.frame(d2_uniques_lab)  
#add X variable for joining with the sample dataframe
d2_uniques_lab$X <- d2_l_glmm$X1

d2_uniques_lab <- rename(d2_uniques_lab, "uniques" = "d2_uniques_lab")
#join with sample dataframe
d2_l_glmm <- d2_l_glmm %>%
  left_join(d2_uniques_lab, by =c("X1" = "X"))

#create a unique ID for each species for the Permanova GLMER
d2_unique_lab <- unique(d2_uniques_lab[c("uniques")]) #12 unique taxonomies
#give these numbers here:
d2_unique_lab$species_ID <- seq.int(nrow(d2_unique_lab))

#join this with our sample dataframe
d2_l_glmm <- d2_l_glmm %>%
  left_join(d2_unique_lab, by = "uniques")

#make this a long dataframe for GLMer
d2_lab_glmm <- d2_l_glmm %>%
  gather(sample, reads, HEV07:HEV29)

#combines reads from each sample from the same unique "species" and adds up the reads
d2_lab_glmm <- d2_lab_glmm %>%
  group_by(sample, uniques) %>% 
  summarise(reads=sum(reads)) %>%
  ungroup()

#dada2: lab GLMER PERMANOVA ####

#add sample meta-data from sample dataframe
d2_lab_glmm <- d2_lab_glmm %>%
  inner_join(samples, by = "sample")

#set a field for presence/absence for analysis with PERMANOVA
d2_lab_glmm$Presence <- ifelse(d2_lab_glmm$reads >0,1,0)

#adds the unique sample ID value to this dataframe so we can call it as a 
#random effect below
d2_lab_glmm <- d2_lab_glmm %>%
  left_join(d2_unique_lab, by = "uniques")

#set species ID to a factor
d2_lab_glmm$species_ID <- as.factor(d2_lab_glmm$species_ID)
d2_lab_glmm$Ster <- ifelse(d2_lab_glmm$Sterilized == "NS", 0, 1)
#model that says that the slope of the response to sterilization can vary
#randomly by Species_ID
d2_lab_model <- glmer(Presence ~ Ster + (1+Ster|species_ID), 
                      data = d2_lab_glmm, 
                      glmerControl(optimizer = "bobyqa"),
                      family = binomial)

#looks at the summary of this model
summary(d2_lab_model)
plot(residuals(d2_lab_model))
#isSingular(d2_lab_model) #FUCK, it's singular...
#dotplots with species which differ between communities being shown in the Ster
#dotplot without error bars crossing zero
#dotplot(ranef(d2_lab_model,
#              condVar =T))


#this is the null model with just a random intercept by species ID
d2_lab_no_ster <- glmer(Presence ~ 1 + (1|species_ID), 
                        data = d2_lab_glmm, 
                        glmerControl(optimizer = "bobyqa"),
                        family = binomial)

#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(d2_lab_model, d2_lab_no_ster)

#dada2: Lab NMDS####
#i'm going to do an NMDS of this to visualize the community
#using d2_lab_glmm
#select columns of interest and then spread back to a community
#matrix for NMDS with the samples as row names
d2_perm_lab <- d2_lab_glmm %>%
  dplyr::select(sample, uniques, reads) %>%
  spread(uniques, reads) %>%
  column_to_rownames("sample")

#run NMDS based on jaccard dissimilarity
d2_nmds_lab <- metaMDS(d2_perm_lab, distance = "jaccard", trymax = 1000,binary=TRUE, k=2) 
#evaluate the fit of this NMDS
stressplot(d2_nmds_lab) 
d2_nmds_lab$stress #0.0873047 looks good!
#create a dataframe for plotting
d2_nmds_lab_df <- data.frame(MDS1=d2_nmds_lab$points[,1], MDS2=d2_nmds_lab$points[,2])
#combining with metadata for plotting
d2_nmds_lab_meta <- cbind(d2_meta_lab, d2_nmds_lab_df)

#plots NMDS with elipses
d2_l_nmds_graph <- ggplot(d2_nmds_lab_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    labs(title = "DADA2 Lab Jaccard Dissimilarity") +
    theme_bw() +
    scale_color_manual(values = c("#7B8D65","#F29979")) 

#adonis Permanova presence/absence ####
#for permanova with adonis
#with predator ASVs
adonis(fld_u3_nmds ~ Sterilized, data = u3_meta_field, dist = "jaccard", binary=TRUE)

adonis(lab_u3_nmds ~ Sterilized, data = u3_meta_lab, dist = "jaccard", binary = TRUE)

adonis(fld_d2_nmds ~ Sterilized, data = d2_meta_field, dist = "jaccard", binary=TRUE)

adonis(lab_d2_nmds ~ Sterilized, data = d2_meta_lab, dist = "jaccard", binary = TRUE)
#these adonis calls are ALSO not significant!

#I was thinking that a full model with pipeline*sterilized would be good since 
#that would compare whether the pipliens produce different results. However, then
#I rememeberd that the list of ASV "species" for dada2 is different from unoise,
#so I feel okay with analyzing them separately given that. 

#dotplot visulizations in ggplot: https://stackoverflow.com/questions/13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-how-to-mak
#https://bbolker.github.io/morelia_2018/notes/mixedlab.html

#Plot Visulaizations together####
#NMDS Plots####
plot_grid(u3_f_nmds_graph, d2_f_nmds_graph, align = "h")
plot_grid(u3_l_nmds_graph, d2_l_nmds_graph, align = "h")

plot_grid(u3_f_nmds_graph, d2_f_nmds_graph,
          u3_l_nmds_graph, d2_l_nmds_graph, nrow = 2)

#Dotplots of random effects ####
#dotplot of random effects
re_u3f <- ranef(u3_field_model) #make a random effect object
dd_u3f <- as.data.frame(re_u3f) #make that into a dataframe
dd_u3f <- dd_u3f %>% #select only the random intercept grouper, removing "intercept"
  filter(term == "Ster")
unique_field_u3$species_ID <- as.factor(as.character(unique_field_u3$species_ID))
dd_u3f <- dd_u3f %>%
  left_join(unique_field_u3, by=c("grp" = "species_ID"))

#plot this, which shows by-species variation within the community
(u3f_re_plot <- ggplot(dd_u3f, aes(y=uniques, x = condval)) +
  geom_point() + 
  geom_errorbarh(aes(xmin = condval - 2*condsd,
                     xmax = condval + 2*condsd), height =0) +
  labs(title = "UNOISE Field: Effect of sterilization by species") +
  geom_vline(xintercept = 0, color = "gray") + theme_bw())

#dotplot of random effects
re_d2f <- ranef(d2_field_model) #make a random effect object
dd_d2f <- as.data.frame(re_d2f) #make that into a dataframe
dd_d2f <- dd_d2f %>% #select only the random intercept grouper, removing "intercept"
  filter(term == "Ster")
unique_field_d2$species_ID <- as.factor(as.character(unique_field_d2$species_ID))
dd_d2f <- dd_d2f %>%
  left_join(unique_field_d2, by=c("grp" = "species_ID"))
#plot this, which shows by-species variation within the community
(d2f_re_plot <- ggplot(dd_d2f, aes(y=uniques, x = condval)) +
    geom_point() + 
    geom_errorbarh(aes(xmin = condval - 2*condsd,
                       xmax = condval + 2*condsd), height =0) +
    labs(title = "DADA2 Field: Effect of sterilization by species") +
    geom_vline(xintercept = 0, color = "gray") + theme_bw())

plot_grid(u3f_re_plot, d2f_re_plot, align = "h")

#dotplot of random effects
re_u3l <- ranef(u3_lab_model) #make a random effect object
dd_u3l <- as.data.frame(re_u3l) #make that into a dataframe
dd_u3l <- dd_u3l %>% #select only the random intercept grouper, removing "intercept"
  filter(term == "Ster")
u3_unique_lab$species_ID <- as.factor(as.character(u3_unique_lab$species_ID))
dd_u3l <- dd_u3l %>%
  left_join(u3_unique_lab, by=c("grp" = "species_ID"))

#plot this, which shows by-species variation within the community
(u3l_re_plot <- ggplot(dd_u3l, aes(y=uniques, x = condval)) +
    geom_point() + 
    geom_errorbarh(aes(xmin = condval - 2*condsd,
                       xmax = condval + 2*condsd), height =0) +
    labs(title = "UNOISE Lab: Effect of sterilization by species") +
    geom_vline(xintercept = 0, color = "gray") + theme_bw())

#dotplot of random effects
re_d2l <- ranef(d2_lab_model) #make a random effect object
dd_d2l <- as.data.frame(re_d2l) #make that into a dataframe
dd_d2l <- dd_d2l %>% #select only the random intercept grouper, removing "intercept"
  filter(term == "Ster")
d2_unique_lab$species_ID <- as.factor(as.character(d2_unique_lab$species_ID))
dd_d2l <- dd_d2l %>%
  left_join(d2_unique_lab, by=c("grp" = "species_ID"))
#plot this, which shows by-species variation within the community
(d2l_re_plot <- ggplot(dd_d2l, aes(y=uniques, x = condval)) +
    geom_point() + 
    geom_errorbarh(aes(xmin = condval - 2*condsd,
                       xmax = condval + 2*condsd), height =0) +
    labs(title = "DADA2 Lab: Effect of sterilization by species") +
    geom_vline(xintercept = 0, color = "gray") + theme_bw())

plot_grid(u3l_re_plot, d2l_re_plot, align = "h")

#line plots of random effects slopes####
#this is a line graph to visualize the per sample slopes
#unoise field: sig groups:
#1, 6, 3, 7, 11, 10, 8 ,14
dotplot(ranef(u3_field_model, condVar=TRUE))
me_f_u3 <- ggpredict(u3_field_model, terms = c("Ster", "species_ID"), type = "re")
me_f_u3 <- me_f_u3 %>%
  left_join(unique_field_u3, by = c("group" = "species_ID"))

me_f_u3$sig <- if_else(me_f_u3$group %in% c(1, 6, 3, 7, 11, 10, 8 ,14), "sig", "non-sig")

me_f_u3sig <- me_f_u3[which(me_f_u3$sig == "sig"),]
pres <- me_f_u3sig %>%
  rename("Sterilized" = "x") %>%
 group_by(uniques, Sterilized) %>%
  summarize(Presence = mean(predicted)) 

me_f_u3sig %>%
  group_by(x) %>%
  summarize(se = sd(predicted)/sqrt(nrow(me_f_u3sig))) %>%
  summarize(mean= mean(se))

pres <- pres %>%
  spread(Sterilized, Presence) %>%
  mutate("Difference" = (`0` - `1`)) %>%
  rename("Unsterilized Presence" = `0`) %>%
  rename("Sterilized Presence" = `1`)

knitr::kable(pres, row.names=FALSE)

#plot these 
(me_fu3_line <- ggplot(me_f_u3, aes(x, predicted, color = uniques)) +
  geom_smooth(aes(linetype = sig), method = "lm") +
  scale_linetype_manual(values=c("dashed", "solid")) +
  labs(x = "Sterilization", y = "Predicted Presence", 
       title = "UNOISE Field marginal means") +
  scale_x_continuous(breaks = c(0,1)) +
  theme_bw())

#unoise lab: 
#dotplot(ranef(u3_lab_model, condVar=TRUE))
#3,1,2,6,4
#this is a line graph to visualize the per sample slopes
#me_l_u3 <- ggpredict(u3_lab_model, terms = c("Ster", "species_ID"), type = "re")
#me_l_u3$group <- as.numeric(as.character(me_l_u3$group))
#me_l_u3 <- me_l_u3 %>%
#  left_join(u3_unique_lab, by = c("group" = "species_ID"))

#me_l_u3$sig <- "non-sig"

#plot these 
#(me_lu3_line <- ggplot(me_l_u3, aes(x, predicted, color = uniques)) +
#    geom_smooth(aes(linetype = sig), method = "lm", linetype = "dashed") +
#    labs(x = "Sterilization", y = "Predicted Presence", 
#         title = "UNOISE Lab marginal means") +
#    scale_x_continuous(breaks = c(0,1)) +
#    theme_bw())

#dada2 FIELD
#sig groups: 1, 5, 15, 14, 12, 11, 10, 3
#me_f_d2 <- ggpredict(d2_field_model, terms = c("Ster", "species_ID"), type = "re")

#me_f_d2$group <- as.numeric(as.character(me_f_d2$group))

#me_f_d2 <- me_f_d2 %>%
#  left_join(unique_field_d2, by = c("group" = "species_ID"))

#me_f_d2$sig <- if_else(me_f_d2$group %in% c(1, 5, 15, 14, 12, 11, 10, 3), "sig", "non-sig")

#plot these 
#(me_fd2_line <- ggplot(me_f_d2, aes(x, predicted, color = uniques)) +
#    geom_smooth(aes(linetype = sig), method = "lm") +
#    scale_linetype_manual(values=c("dashed", "solid")) +
#    labs(x = "Sterilization", y = "Predicted Presence", 
#         title = "DADA2 Field marginal means") +
#    scale_x_continuous(breaks = c(0,1)) +
#    theme_bw())

#me_f_d2sig <- me_f_d2[which(me_f_d2$sig == "sig"),]
#me_f_d2sig %>%
#  group_by(x) %>%
#  summarize(mean = mean(predicted)) 
#0.436-0.372
#me_f_d2sig %>%
#  group_by(x) %>%
#  summarize(se = sd(predicted)/sqrt(nrow(me_f_d2sig))) %>%
 # summarize(mean = mean(se))

#dada2 Lab
#sig groups 1, 5, 4
#me_l_d2 <- ggpredict(d2_lab_model, terms = c("Ster", "species_ID"), type = "re")

#me_l_d2$group <- as.numeric(as.character(me_l_d2$group))

#me_l_d2 <- me_l_d2 %>%
#  left_join(d2_unique_lab, by = c("group" = "species_ID"))

#me_l_d2$sig <- if_else(me_l_d2$group %in% c(1, 5, 4), "sig", "non-sig")

#plot these 
#(me_ld2_line <- ggplot(me_l_d2, aes(x, predicted, color = uniques)) +
#    geom_smooth(aes(linetype = sig), method = "lm") +
#    scale_linetype_manual(values=c("dashed", "solid")) +
#    labs(x = "Sterilization", y = "Predicted Presence", 
#         title = "DADA2 Lab marginal means") +
#    scale_x_continuous(breaks = c(0,1)) +
#    theme_bw())

#me_l_d2sig <- me_l_d2[which(me_l_d2$sig == "sig"),]
#me_l_d2sig %>%
#  group_by(x) %>%
#  summarize(mean = mean(predicted)) 
#0.856-0.822
#me_l_d2sig %>%
#  group_by(x) %>%
#  summarize(se = sd(predicted)/sqrt(nrow(me_l_d2sig))) %>%
#  summarize(mean = mean(se))

#graph: known Diet presence####
u3_lab_glmm_known <- u3_lab_glmm[which(u3_lab_glmm$uniques == "Oxya"),]
d2_lab_glmm_known <- d2_lab_glmm[which(d2_lab_glmm$uniques == "Acrididae"),]

absent_u3 <- u3_lab_glmm_known %>%
  group_by(Sterilized) %>%
  tally(Presence <= 0) %>%
  mutate(pipeline = "u3", total=nrow(u3_lab_glmm_known))

absent_d2 <- d2_lab_glmm_known %>%
  group_by(Sterilized) %>%
  tally(Presence <= 0) %>%
  mutate(pipeline = "d2", total=nrow(d2_lab_glmm_known))

absent <- bind_rows(absent_u3, absent_d2)

absent_plot <- ggplot(absent, aes(x = Sterilized, y = (total-n)/total*100, fill = pipeline)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Percent of samples", title = "Percent of samples with known diet present") +
  scale_fill_manual(values = c("#734646","#F27D72")) +
  theme_bw()
absent_plot
