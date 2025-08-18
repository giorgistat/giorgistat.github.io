rm(list = ls())  # Clear workspace
library(RiskMap)  # Load RiskMap package for spatial modeling

# Load Required Packages
library(ggplot2)  # For data visualization
library(sf)  # For spatial data handling
library(RiskMap)  # Spatial risk mapping functions

# Load data-sets
# These datasets contain spatial and epidemiological data
data(liberia)
data(malkenya)
data(anopheles)

# Exploring associations with risk factors using count data

# Calculate prevalence as proportion of positive cases
tibble::glimpse(liberia)  # View dataset structure
liberia$prev <- liberia$npos/liberia$ntest

# Plot prevalence against elevation
ggplot(liberia, aes(x = elevation, y = prev)) + geom_point() +
  labs(x="Elevation (meters)",y="Prevalence")

# Compute the empirical logit transformation to stabilize variance
liberia$elogit <- log((liberia$npos+0.5)/
                        (liberia$ntest-liberia$npos+0.5))

# Visualizing empirical logit transformation
ggplot(liberia, aes(x = elevation, y = elogit)) + geom_point() +

  # Adding a smoothing spline to explore trends
  labs(x="Elevation (meters)",y="Empirical logit") +
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+

  # Adding linear regression fit with log-transformed elevation
  stat_smooth(method = "lm", formula = y ~ log(x),
              col="green",lty="dashed",se=FALSE) +

  # Adding linear regression fit with change point in 150 meters
  stat_smooth(method = "lm", formula = y ~ x + pmax(x-150, 0),
              col="red",lty="dashed",se=FALSE)


####################

library(RiskMap)
library(lme4)
data("liberia")

# Create the ID of the location
liberia$ID_loc <- 1:nrow(liberia)

# Binomial mixed model with log-elevation
fit_glmer_lib <- glmer(cbind(npos, ntest) ~ log(elevation) + (1|ID_loc), family = binomial,
                       data = liberia, nAGQ = 25)

summary(fit_glmer_lib)


# Empirical variogram


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
