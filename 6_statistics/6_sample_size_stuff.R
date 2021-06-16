library(simr)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(pwr)
library(pscl)
library(here)
library(lmtest)
library(pbkrtest)

###########################
# Load dataset ####
############################

#the lab dataframe is already in the format it needs to be for analysis, with
#the addition of a presence-absence column
lab_detect <- read.csv(here("data", 
                            "outputs", 
                            "rarefied_taxonomic_sort", 
                            "lab_known_prey_rare.csv"))

###########################
# Create presence-absence column in DFs ####
############################
lab_detect <- lab_detect %>%
  mutate(presence = ifelse(reads > 0, 1, 0))

###########################
# Mesocosm Known prey detection ####
############################           

lab_detect_mod <- glm(presence ~ Sterilized,
                          data = lab_detect, 
                          family = binomial)

lab_detect_null <- glm(presence ~ 1,
                      data = lab_detect, 
                      family = binomial)

lrtest(lab_detect_mod, lab_detect_null)

p <- coef(summary(lab_detect_mod))[,"Pr(>|z|)"]
p <- p[2]
###########################
# Summary stats to get out variable for power test ####
############################ 

summary(lab_detect_mod)

u <- 1 #number of df from model
v <- 19-u-1 #sample size minus u minus 1
r2 <- pR2(lab_detect_mod) #want McFadden
R2 <- r2[4] #select McFadden pseudo R^2 from list
f2 <- R2/(1-R2)

###########################
# Power test of observed ####
############################ 
#https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html
pwr.f2.test(u = u,
            v=v,
            sig.level = 0.05,
            f2 = f2)

pwr.f2.test(u = u,
            sig.level = 0.05,
            f2 = f2,
            power = 0.8)

v2 <- 34.01722
n <- v2 + u + 1 #sample sizes estimate from v
#36


#############################################
#P-value function for a bunch of random binomial dists ####
#############################################
#simulate the binomial distribution
#https://stats.stackexchange.com/questions/258175/simulate-data-based-on-output-from-a-generalized-linear-model-in-r

sim_p <- function(mod = lab_detect_mod){
  y <- rbinom(19,1,fitted(mod,method="response")) #simulate binomial from model
  x <- c(rep("SS",8), rep("NS", 11)) #create treatment groups
  #x <- c(rep("SS",15), rep("NS", 21))
  df <- as.data.frame(cbind(x,y)) #make this a df
  
  df <- df %>%
    mutate(y = as.numeric(y)) #make y numeric
  
  mod <- glm(y ~ x, #full model
             data = df,
             family = binomial)
  
  p <- coef(summary(mod))[,"Pr(>|z|)"] #pull out p-value
  p <- p[2] #get p-value for predictor
  
  return(p) #return the p-value for the LR (H0: models are the same)
}


# Simulate over 1000 iterations -------------------------------------------

## Number of simulations
nsims <- 1000
## Pre-allocate emptyvector to memory
pVals <- rep(NA, nsims)
## Run simulations
for(i in 1:nsims){
  pVals[i] <- sim_p(lab_detect_mod)
}

pVals <- as.data.frame(pVals)
str(pVals)


# Distribution of p-values ------------------------------------------------


(small <- ggplot(pVals, aes(x = pVals)) +
  geom_histogram() +
  labs(x= "Sterilization effect p-value", y = "Count") +
  geom_vline(xintercept = 0.05, linetype = "dashed", size = 0.7) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20)) +
  annotate(geom = "text", x = 0.23, y = 200, label = "p-value = 0.05", size = 5))

pVals %>%
  mutate(reject = ifelse(pVals <= 0.05, "reject", "no_reject")) %>%
  group_by(reject) %>% 
  tally()
140/(140+860)
pVals <- pVals %>% 
  mutate(sim = "small")

# n = 36 function -------------------------------------------------------------

sim_p2 <- function(mod = lab_detect_mod){
  y <- rbinom(36,1,fitted(mod,method="response")) #simulate binomial from model
  x <- c(rep("SS",15), rep("NS", 21)) #create treatment groups
  #x <- c(rep("SS",15), rep("NS", 21))
  df <- as.data.frame(cbind(x,y)) #make this a df
  
  df <- df %>%
    mutate(y = as.numeric(y)) #make y numeric
  
  mod <- glm(y ~ x, #full model
             data = df,
             family = binomial)
  
  p <- coef(summary(mod))[,"Pr(>|z|)"] #pull out p-value
  p <- p[2] #get p-value for predictor
  
  return(p) #return the p-value for the LR (H0: models are the same)
  
}


# Simulate over 1000 iterations -------------------------------------------

## Number of simulations
nsims <- 1000
## Pre-allocate emptyvector to memory
pVals2 <- rep(NA, nsims)
## Run simulations
for(i in 1:nsims){
  pVals2[i] <- sim_p2(lab_detect_mod)
}

pVals2 <- as.data.frame(pVals2)
pVals2 <- pVals2 %>%
  mutate(sim = "big") %>%
  rename("pVals" = "pVals2")

str(pVals2)

# Distribution of p-values ------------------------------------------------


big <- ggplot(pVals2, aes(x = pVals2)) +
  geom_histogram() +
  geom_vline(xintercept = 0.07, linetype = "dashed")+
  theme_bw()

pVals2 %>%
  mutate(reject = ifelse(pVals2 <= 0.07, "reject", "no_reject")) %>%
  group_by(reject) %>% 
  tally()

22/(22+978)

p_v <- pVals %>%
  bind_rows(pVals2)
  
ggplot(p_v, aes(x = pVals, fill = sim)) +
  geom_histogram(position = "dodge") +
  geom_vline(xintercept = 0.07, linetype = "dashed") +
  theme_bw()

small
big
