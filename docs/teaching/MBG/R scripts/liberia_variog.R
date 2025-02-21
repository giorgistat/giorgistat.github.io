rm(list = ls())

library(RiskMap)
library(lme4)
library(sf)
data("liberia")

# Create the ID of the location
liberia$ID_loc <- 1:nrow(liberia)

# Binomial mixed model with log-elevation
fit_glmer_lib <- glmer(cbind(npos, ntest) ~ log(elevation) + (1|ID_loc), family = binomial,
                       data = liberia, nAGQ = 25)

summary(fit_glmer_lib)

liberia$Z_hat <- ranef(fit_glmer_lib)$ID_loc[,1]

liberia <- st_as_sf(liberia, coords = c("long", "lat"),
                    crs=4326)
crs_lb <- propose_utm(liberia)
liberia <- st_transform(liberia, crs = crs_lb)

liberia_variog <- s_variogram(data = liberia,
                              variable = "Z_hat",
                              bins = c(15, 30, 40, 80, 120,
                                       160, 200, 250, 300, 350),
                              scale_to_km = TRUE,
                              n_permutation = 10000)

plot_s_variogram(liberia_variog)


plot_s_variogram(liberia_variog,
                 plot_envelope = TRUE)
