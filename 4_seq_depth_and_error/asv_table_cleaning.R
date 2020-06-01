#Correcting Error using negatives (from Jerde)####
#Ana Miller-ter Kuile
#February 5, 2020
library(here)
library(tidyverse)
library(MASS) #fitdistr function
library(fitdistrplus)
library(DHARMa)

#Since UNOISE3 was best for ecological assessments in the pipeline comparisons,
#I'm going to overlook the error in this pipeline for now and try to correct it
#using the method from Chris Jerde in the eDNA workshop

#load the dataframe
u3_error_check <- read.csv(here("data", "denoised_data", "ASV_tables", "unoise_uc_zotu_tab.txt"), sep = "\t")
#rename columns for simplicity
colnames(u3_error_check) <- sapply(str_split(colnames(u3_error_check), "S"), function(x){return(x[[1]])})
u3_error_check <- rename(u3_error_check, "ASV" = "X.OTU.ID")

#Error Analyses####
#Method 1: this is from Jerde in eDNA workshop on error rates and figuring out errors that may
#throw off species richness estimates
#look at a histogram of values
error_check <- u3_error_check %>%
  dplyr::select(ASV, NEG)

hist(error_check$NEG)

#ggplot(error_check, aes(x = NEG)) +
#         geom_density(fill = "white", alpha = 0.5, position = "identity") +
#  theme_bw()

model1 <- glm(NEG ~ 1, family = poisson, data = error_check)
#plot(residuals(model1))
simulate_poi <- simulateResiduals(fittedModel = model1)
mod1fit <- plot(simulate_poi, asFactor=TRUE) #these look good
dev.off()
mod1zi <- testZeroInflation(simulate_poi) #not zero inflated
#mod1od <- testDispersion(simulate_poisson) #maybe overdispersed???

model2 <- glm.nb(NEG ~ 1, data = error_check)
#plot(residuals(model2))
BIC(model1, model2) #this seems to indicate poisson is better, but I have tons 
#of zero inflation, not really sure what to do with this - look into it

fitdistr(error_check$NEG, "Poisson")
#lambda 0.005681818
#fitdistr(error_check$NEG, "negative binomial")
#size = 100, mu = .005681818

ppois(1,lambda=0.005681818,lower=FALSE) # answer is <<0.0001
#ppois(2,lambda=0.005681818,lower=FALSE) # answer is <<0.0001
#ppois(0, lambda =0.005681818, lower=FALSE ) #answer is 0.006
#?pnbinom
#pnbinom(1, size = 100, mu = 0.00568181818, lower = FALSE) #<<0.0001
#pnbinom(0, size = 100, mu = 0.00568181818, lower = FALSE) #0.006

#nothing comes from the error distribution? seems like it...
#going to go ahead and say that a 1 is a real value in this dataset
