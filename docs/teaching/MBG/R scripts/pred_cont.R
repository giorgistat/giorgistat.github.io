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


# Obtaining the grid

library(rgeoboundaries)

shp_lb <- geoboundaries(country = "liberia", adm_lvl = "adm0")
shp_lb <- st_transform(shp_lb, crs=crs_lb)


# Create the grid

grid_lb <- create_grid(shp_lb, spat_res = 5)


# Extract elevation

library(elevatr)
elevation <- get_elev_point(st_as_sf(grid_lb), prj = crs_lb, src = "aws")$elevation

pred_lb <- data.frame(elevation=elevation)

# Predict S(x)

pred_S_lb <-
pred_over_grid(lb_fit, grid_pred = grid_lb,
               predictors = pred_lb)

# Predict T(x)

pred_T_lb <-
  pred_target_grid(pred_S_lb,
                   f_target = list(prev = function(x) exp(x)/(1+exp(x))),
                   pd_summary = list(mean = mean,
                                     lower_lim = function(x) quantile(x, 0.025),
                                     upper_lim = function(x) quantile(x, 0.975),
                                     exceed20 = function(x) mean(x > 0.2),
                                     exceed30 = function(x) mean(x > 0.3)))


plot(pred_T_lb, which_target = "prev", which_summary = "mean",
     main = "Predictive mean")

plot(pred_T_lb, which_target = "prev", which_summary = "lower_lim",
     main = "Quantile 0.025")

plot(pred_T_lb, which_target = "prev", which_summary = "upper_lim",
     main = "Quantile 0.975")

plot(pred_T_lb, which_target = "prev", which_summary = "exceed20",
     main = "Exceedance probability (L = 0.2)")

plot(pred_T_lb, which_target = "prev", which_summary = "exceed30",
     main = "Exceedance probability (L = 0.3)")
