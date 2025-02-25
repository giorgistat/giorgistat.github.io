rm(list = ls())

library(RiskMap)
library(sf)
library(ggplot2)

data("liberia")

# Convert to an sf object
liberia_sf <- st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
crs_lb <- propose_utm(liberia_sf)
liberia_sf <- st_transform(liberia_sf)


# Fitting a Binomial geostatistical model

lb_fit <- glgpm(npos ~ log(elevation) + gp(),
                den = ntest,
                family = "binomial",
                data = liberia_sf)

check_mcmc(lb_fit)
check_mcmc(lb_fit, check_mean = FALSE, component = sample(1:nrow(liberia_sf), 1))

# Updating the MCMCL
theta0 <- coef(lb_fit)

lb_fit <- glgpm(npos ~ log(elevation) + gp(),
                den = ntest,
                par0 = theta0,
                family = "binomial",
                data = liberia_sf)

summary(lb_fit)


#####################################################################################
#####################################################################################

data("anopheles")

anopheles_sf <- st_as_sf(anopheles, coords = c("web_x", "web_y"), crs = 3857)


# Fitting a Poisson geostatistical model

an_fit <- glgpm(An.gambiae ~ elevation + gp(),
                family = "poisson",
                data = anopheles_sf)

check_mcmc(an_fit)
check_mcmc(an_fit, check_mean = FALSE, component = sample(1:nrow(anopheles_sf), 1))

# Updating the MCMCL
theta0 <- coef(an_fit)

an_fit <- glgpm(An.gambiae ~ elevation + gp(),
                par0 = theta0,
                family = "poisson",
                data = anopheles_sf)


summary(an_fit)
