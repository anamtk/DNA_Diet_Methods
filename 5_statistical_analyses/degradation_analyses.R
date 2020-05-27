dept <- u3_commr %>%
  gather(sample, reads, HEV07:HEV29) %>%
  group_by(sample) %>%
  summarize(reads = sum(reads)) %>% 
  left_join(metadata, by = "sample") 

ggplot(dept, aes(x = Sterilized, y = reads)) +
  geom_boxplot()

mod <- glmmTMB(reads ~ Sterilized,
               data=dept, 
               family = "genpois")

mod_null <- glmmTMB(reads ~ 1,
               data=dept, 
               family = "genpois")

AICc(mod, mod_null)

plot(allEffects(mod))

simulationOutput <- simulateResiduals(fittedModel = mod) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput) 

dept_f <- u3_commf %>%
  gather(sample, reads, HEV65:HEV100) %>%
  group_by(sample) %>%
  summarize(reads = sum(reads)) %>%
  left_join(metadata, by = "sample")

ggplot(dept_f, aes(x = Sterilized, y = reads)) +
  geom_boxplot()

mod1 <- glmmTMB(reads ~ Sterilized,
               data=dept_f, 
               family = "genpois")

mod_1null <- glmmTMB(reads ~ 1,
                    data=dept_f, 
                    family = "genpois")

AICc(mod1, mod_1null)

plot(allEffects(mod1))
plot(residuals(mod1))
simulationOutput <- simulateResiduals(fittedModel = mod_1null) 
fit <- plot(simulationOutput, asFactor=TRUE)
zi <- testZeroInflation(simulationOutput) 
od <- testDispersion(simulationOutput) 

summary(mod1)
