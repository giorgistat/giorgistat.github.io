rm(list = ls())

library(RiskMap)
library(sf)
library(ggplot2)

data("liberia")

# Convert to an sf object
liberia_sf <- st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
crs_lb <- propose_utm(liberia_sf)
liberia_sf <- st_transform(liberia_sf, crs = crs_lb)


# Fitting a Binomial geostatistical model

lb_fit <- glgpm(npos ~ log(elevation) + gp(),
                den = ntest,
                family = "binomial",
                data = liberia_sf)

summary(lb_fit)

# Using an alternative link function

lb_fit_probit <- glgpm(npos ~ log(elevation) + gp(),
                den = ntest,
                invlink = function(x) pnorm(x),
                family = "binomial",
                data = liberia_sf)

summary(lb_fit_probit)
