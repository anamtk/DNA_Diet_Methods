###########################
#Sampling Depth
#Ana Miller-ter Kuile
#February 16, 2020
##########################

#this is code for looking at sampling depth across samples to ensure that I've
#sufficiently sampled each sample in this dataset. 
#this code looks at sampling depth in community matrices created by both dada2 and unoise3

#######################
#Load packages ####
######################

library(tidyverse)
library(here)
library(ggplot2)
library(iNEXT)
library(cowplot)

#######################
#Load and tidy data ####
######################
#needed here: ASV matrix by samples minus the controls and asv column
#import data matrix
u3_comm <- read.csv(here("data", 
                         "data_raw",
                         "1_denoised_data", 
                         "ASV_tables", 
                         "unoise_uc_zotu_tab.txt"), sep = "\t")

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

#######################
#iNEXT sequencing depth####
######################

#this determines sequencing depth across all samples
u3_seq_depth <- iNEXT(u3_comm_depth, q=0, datatype="abundance") #this determines sequencing depth for each sample

#u3_seq_depth$AsyEst
u3_seq_depth$DataInfo$SC #look at sampling completeness
u3_seq_depth$DataInfo #gets per sample richness, etc
#other parts of the iNEXT object that could be useful:
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

#######################
#Split field and lab####
######################

meta <- read.csv(here("data", 
                      "data_raw",
                      "2_sample_data",
                      "Sample_Metadata.csv"))

lab <- meta %>%
  filter(Source == "LAB") %>%
  dplyr::select(sample) %>%
  as_vector()

field <- meta %>%
  filter(Source == "FIELD") %>%
  dplyr::select(sample) %>%
  as_vector()

comm_lab <- u3_comm_depth %>%
  dplyr::select(all_of(lab))

comm_field <- u3_comm_depth %>%
  dplyr::select(all_of(field))

#######################
#Field iNEXT####
######################

#this determines sequencing depth across all samples
field_depth <- iNEXT(comm_field, q=0, datatype="abundance") #this determines sequencing depth for each sample


field_depth$DataInfo$SC #look at sampling completeness
field_depth$DataInfo #gets per sample richness, etc
#other parts of the iNEXT object that could be useful:

#graph the interpolated and extrapolated sampling depth per sample
df <- fortify(field_depth, type=1)
head(df)

df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method,
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))

fld_seq_depth <- ggplot(df, aes(x=x, y=y, colour=site)) +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
  theme_bw() +
  labs(x="Sequencing Depth", y="ASV diversity") +
  ylim(0, 35) +
  geom_vline(xintercept = 16004, linetype = "dashed") +
  theme(legend.position = "none",
        text=element_text(size=18), axis.title.y = element_blank()) 

#######################
#Lab iNEXT####
######################

#this determines sequencing depth across all samples
lab_depth <- iNEXT(comm_lab, q=0, datatype="abundance") #this determines sequencing depth for each sample


lab_depth$DataInfo$SC #look at sampling completeness
lab_depth$DataInfo #gets per sample richness, etc
#other parts of the iNEXT object that could be useful:

df <- fortify(lab_depth, type=1)
head(df)

df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method,
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))

lab_seq_depth <- ggplot(df, aes(x=x, y=y, colour=site)) +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
  theme_bw() +
  ylim(0, 35) +
  labs(x="Sequencing Depth", y="ASV diversity") +
  geom_vline(xintercept = 55205, linetype = "dashed") +
  theme(legend.position = "none",
        text=element_text(size=18)) 

#######################
#paper output####
######################
plot_grid(lab_seq_depth, fld_seq_depth,  align = "vh")
