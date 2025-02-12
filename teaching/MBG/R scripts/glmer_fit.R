rm(list = ls())

library(RiskMap)
library(lme4)
data("liberia")

# Create the ID of the location
liberia$ID_loc <- 1:nrow(liberia)

# Binomial mixed model with log-elevation
fit_glmer_lib <- glmer(cbind(npos, ntest) ~ log(elevation) + (1|ID_loc), family = binomial,
                       data = liberia, nAGQ = 25)

summary(fit_glmer_lib)


