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
#2. how many known prey ASVs/ASV abundance of known prey in lab individuals by pipeline
#and by sterilization?

#Thoughts February 20####
#here, I need to think a bit more about what taxonomic level things need to be 
#broken up into. Species ID seems to be misleading after going through, so I'm
#thinking about ORDER for everything that is here, and then for the remaining
#ASVs that got matched in MEGAN, giving categories such as "higher-taxonomy potential
#prey", "higher-taxonomy non-prey", and "no hits"

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

#(Now that I've moved some code to another script, need to consider here if any of the
#dataframes here need to be re-imported again)

#Rarefy for abundance-based dissimilarity####
#field, will be working with fld_d2uc_nmds
fld_d2uc_nmds <- read.csv(here("data", "fld_d2uc_nmds.csv"))
meta_field <- read.csv(here("data", "fld_d2_nmds_meta.csv"))
raremin_fld <- min(rowSums(fld_d2uc_nmds))
raremin_fld
# Rarefy, which means to standardize all samples by the lowest value (raremin). 
# change sample = to the rarefaction depth (# of reads per sample)
rarefied_fld <- rrarefy(fld_d2uc_nmds, sample = raremin_fld)
base::range(rowSums(rarefied_fld)) # checks that it worked
rarefied_fld <- as.data.frame(rarefied_fld) # save it as a dataframe
rowSums(rarefied_fld)
colSums(rarefied_fld)
# Remove OTUs with a read abundance of 0
rarefied_fld <- rarefied_fld[, colSums(rarefied_fld) > 0] #from 154 to 116

#lab will be working with lab_d2uc_nmds
lab_d2uc_nmds <- read.csv(here("data", "lab_d2uc_nmds.csv"))
meta_lab <- read.csv(here("data", "lab_d2uc_nmds_meta.csv"))
raremin_lab <- min(rowSums(lab_d2uc_nmds))
raremin_lab

rarefied_lab <- rrarefy(lab_d2uc_nmds, sample = raremin_lab)
base::range(rowSums(rarefied_lab))
rarefied_lab <- as.data.frame(rarefied_lab)
rowSums(rarefied_lab)

#are any zero across all samples?
colSums(rarefied_lab)

#remove zero columns
rarefied_lab <- rarefied_lab[, colSums(rarefied_lab) > 0] #from 82 to 82, none removed

#Bray-Curtis Abundance based similarity ####
#field
bray_nmds_field <- metaMDS(rarefied_fld, distance = "bray", k=2)
#evaluate the fit of this NMDS
stressplot(bray_nmds_field)
bray_nmds_field$stress

#create a dataframe for plotting
bray_field_df <- data.frame(MDS1=bray_nmds_field$points[,1], MDS2=bray_nmds_field$points[,2])
#combining with metadata for plotting
bray_field_meta <- cbind(meta_field, bray_field_df)

#plots NMDS with elipses
(f <- ggplot(bray_field_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    theme_bw() +
    labs(title = "Field abundance-based dissimilarity (bray-curtis)") +
    scale_color_manual(values = c("#7B8D65","#F29979")) +
    theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20)))

plot_grid(b, f, align="h")

#lab
bray_nmds_lab <- metaMDS(rarefied_lab, distance = "bray", k=2)
#evaluate the fit of this NMDS
stressplot(bray_nmds_lab)
bray_nmds_lab$stress

#create a dataframe for plotting
bray_lab_df <- data.frame(MDS1=bray_nmds_lab$points[,1], MDS2=bray_nmds_lab$points[,2])
#combining with metadata for plotting
bray_lab_meta <- cbind(meta_lab, bray_lab_df)

#plots NMDS with elipses
(g <- ggplot(bray_lab_meta, aes(x=MDS1,y=MDS2, color=Sterilized))+
    geom_point(size = 6)+
    stat_ellipse(level = 0.60, size = 1)+
    coord_fixed()+
    theme_bw() +
    labs(title = "Lab abundance-based dissimilarity (bray-curtis)") +
    scale_color_manual(values = c("#7B8D65","#F29979")) +
    theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20)))

plot_grid(a, g, align="h")

#Adonis PERMANOVA abundance ####
adonis(rarefied_fld ~ Sterilized, data = meta_field, dist = "bray")
adonis(rarefied_lab ~ Sterilized, data = meta_lab, dist = "bray")
#no abundance-based community differences based on ADONIS, will do GLMMs next

#GLMM PERMANOVA abundance ####
#need to bond taxonomy with DF with ASVs on rows
#transpose the community matrix:
rare_f_glmm <- as.data.frame(t(rarefied_fld))
#set an ASV column for combining with the taxonomic assignment dataframe
rare_f_glmm$ASV <- rownames(rare_f_glmm)

#import the taxonomic assignment datafram and then add a row for ASV
taxa <- read_csv(here("data", "d2_uc_tax_ass.csv"))
taxa$ASV <- rownames(taxa)

#combine our matrix of samples and ASVs with the corresponding taxonomies
rare_f_glmm <- rare_f_glmm %>%
  left_join(taxa, by = "ASV")

#rename these variables for ifelse below, and change to characters rather than factors
rare_f_glmm <- rare_f_glmm %>%
  rename(ID_bold = ID.x,
         ID_ncbi = ID.y,
         ID_nondiet = ID)
rare_f_glmm$ID_bold <- as.character(rare_f_glmm$ID_bold)
rare_f_glmm$ID_ncbi <- as.character(rare_f_glmm$ID_ncbi)
rare_f_glmm$ID_nondiet <- as.character(rare_f_glmm$ID_nondiet)

#creates a new vector of the unique species IDs for these data
uniques <- ifelse(!is.na(rare_f_glmm$ID_bold), rare_f_glmm$ID_bold, 
                  ifelse(!is.na(rare_f_glmm$ID_ncbi), rare_f_glmm$ID_ncbi, 
                         ifelse(!is.na(rare_f_glmm$ID_nondiet), rare_f_glmm$ID_nondiet, "no_hit")))

#makes all predators map to Sparassidae, resetting some of the weirdness b/w BOLD and NCBI
uniques[uniques == "Heteropoda venatoria"] <- "Sparassidae"

#makes these unique taxonomies a dataframe                 
uniques <- as.data.frame(uniques)  
#add X variable for joining with the sample dataframe
uniques$X <- rare_f_glmm$X

#join with sample dataframe
rare_f_glmm <- rare_f_glmm %>%
  left_join(uniques, by = "X")

#create a unique ID for each species for the Permanova GLMER
unique_f_rare <- unique(uniques[c("uniques")]) #59 unique taxonomies
#give these numbers here:
unique_f_rare$species_ID <- seq.int(nrow(unique_f_rare))

#join this with our sample dataframe
rare_f_glmm <- rare_f_glmm %>%
  left_join(unique_f_rare, by = "uniques")

#make this a long dataframe for GLMer
rare_field_glmm <- rare_f_glmm %>%
  gather(sample, reads, HEV.100_S22:HEV89_S32)

#combines reads from each sample from the same unique "species" and adds up the reads
rare_field_glmm_uniques <- rare_field_glmm %>%
  group_by(sample, uniques) %>% 
  summarise(reads=sum(reads)) %>%
  ungroup()

#add sample meta-data from sample dataframe
rare_field_glmm_uniques <- rare_field_glmm_uniques %>%
  inner_join(samples, by = c("sample" = "dada2_sample"))

#adds the unique sample ID value to this dataframe so we can call it as a 
#random effect below
rare_field_glmm_uniques <- rare_field_glmm_uniques %>%
  left_join(unique_field, by = "uniques")

#set species ID to a factor
rare_field_glmm_uniques$species_ID <- as.factor(rare_field_glmm_uniques$species_ID)
rare_field_glmm_uniques$Ster <- ifelse(rare_field_glmm_uniques$Sterilized == "NS", 0, 1)
#model that says that the slope of the response to sterilization can vary
#randomly by Species_ID
field_model_ab <- glmer(reads ~ Ster + (1+Ster|species_ID), 
                     data = rare_field_glmm_uniques, glmerControl(optimizer = "bobyqa"),
                     family = poisson)

#looks at the summary of this model
summary(field_model_ab)

#dotplots with species which differ between communities being shown in the SpeciesSS
#dotplot without error bars crossing zero
dotplot(ranef(field_model_ab,
              condVar =T))
#looks like there are a few

#this is the null model with just a random intercept by species ID
field_no_ster_ab <- glmer(reads ~ Ster + (1|species_ID), 
                       data = rare_field_glmm_uniques, glmerControl(optimizer = "bobyqa"),
                       family = poisson)

#ANOVA comparsion showing that the null is a better fit based on AIC and BIC
anova(field_model_ab, field_no_ster_ab)

#this is a line graph to visualize the per sample slopes
me_f <- ggpredict(field_model_ab, terms = c("Ster", "species_ID"), type = "re")

#these are the species IDs which have a signifcant difference bw NS and SS
#Pos: 1, 13,3,17,11,14, 24, 8,23, 15, 53, 21, 40, 41, 37, 28, 
#Neg: 5, 6, 10, 19, 9, 13, 20, 16, 
#minus 1 because it throws off graph
me_pos <- me_f %>%
  filter(group %in% c(13,3,17,11,14, 24, 8,23, 15, 53, 21, 40, 41, 37, 28))

#plot these 
ggplot(me_pos, aes(x, predicted, color = group)) +
  geom_smooth(method = "lm") 

#subset only negative
me_neg <- me_f %>%
  filter(group %in% c(5, 6, 10, 19, 9, 13, 20, 16))

ggplot(me_neg, aes(x, predicted, color = group)) +
  geom_smooth(method = "lm") 

me_sig_f <- me_f %>%
  filter(group %in% c(1,3,5,6,11,12,13,17,19))

#To keep working on ####
#dotplot visulizations in ggplot: https://stackoverflow.com/questions/13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-how-to-mak
#https://bbolker.github.io/morelia_2018/notes/mixedlab.html
#Kyle Edwards lecture 22 search dot
#consider grouping ASVs by "prey" "predator" "non-diet" and "no-hit" AND THEN doing both
#abundance based dissimilarity/PERMANOVA. I think an issue is that
#because so many species are rare based on a species-level approach, that the
#PERMANOVA has a hard time picking up detection variation.
ranef(field_model)
