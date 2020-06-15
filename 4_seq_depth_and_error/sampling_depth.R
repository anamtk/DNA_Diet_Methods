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
#library(ggfortify)
#library(GUniFrac)
#library(vegan)
library(iNEXT)
#library(cowplot)
##library(ecodist)
#library(lme4)
#library(broom)
#library(MASS)
#library(ggeffects)
#library(glmmTMB)
#library(DHARMa)
#library(MuMIn)
#library(effects)

#UNOISE sequencing depth####
#needed here: ASV matrix by samples minus the controls and asv column
#import data matrix
u3_comm <- read.csv(here("data", "denoised_data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = "\t")

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
  labs(x = "Sequencing Depth", y = "ASV Richness") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))
u3_seq_depth_graph

