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
shp_lb_0 <- geoboundaries(country = "liberia", adm_lvl = "adm0")
shp_lb_0 <- st_transform(shp_lb_0, crs=crs_lb)

shp_lb_1 <- geoboundaries(country = "liberia", adm_lvl = "adm1")
shp_lb_1 <- st_transform(shp_lb_1, crs=crs_lb)


# Create the grid

grid_lb <- create_grid(shp_lb_0, spat_res = 5)


# Extract elevation

library(elevatr)
elevation <- get_elev_point(st_as_sf(grid_lb), prj = crs_lb, src = "aws")$elevation

pred_lb <- data.frame(elvation=elevation)

# Predict S(x)

pred_S_lb <-
  pred_over_grid(lb_fit, grid_pred = grid_lb,
                 predictors = pred_lb, type="joint")


## Unweighted regional prevalence

pred_area <- pred_target_shp(pred_S_lb, shp = shp_lb_1,
                            shp_target = function(Tx) mean(Tx),
                            f_target = list(prev =
                                              function(lp) exp(lp)/(1+exp(lp))),
                            pd_summary = list(mean = mean,
                                              exceed20 = function(x) mean(x > 0.2)),
                            col_names = "shapeName")

# Plot point predictions of average prevalence
plot(pred_area, which_target = "prev", which_summary = "mean",
     palette = "RdYlGn",
     limits = c(0.1, 0.30),
     breaks =  seq(0.1, 0.30, by = 0.05)) +
  guides(fill=guide_legend(title="Prevalence")) +
  ggtitle("Average prevalence \n (no weights)") +
  theme(plot.title = element_text(size = 15))

# Plot of the exceedance probabilities
plot(pred_area, which_target = "prev", which_summary = "exceed20",
     palette = "RdYlGn",
     limits = c(0, 1),
     breaks = seq(0,1, by = 0.1)) +
  guides(fill=guide_legend(title="Probability")) +
  ggtitle("Exceedance probability (L = 0.2) \n (no weights)") +
  theme(plot.title = element_text(size = 15))


# Obtaining population density
library(wpgpDownloadR)
lbr_url <- wpgpGetCountryDataset(ISO3 = "LBR", covariate = "ppp_2014")
library(terra)
lbr_pop <- rast(lbr_url)
lbr_pop <- project(lbr_pop, "EPSG:32629")

# Extra pop density weights at the prediction grid
weights_pred <- extract(lbr_pop, st_as_sf(grid_lb))$lbr_ppp_2014


# Prediction of the population weighted regional average prevalence

pred_area_w <- pred_target_shp(pred_S_lb, shp = shp_lb_1,
                              shp_target = function(Tx) sum(Tx),
                              f_target = list(prev =
                                                function(lp) exp(lp)/(1+exp(lp))),
                              pd_summary = list(mean = mean,
                                                exceed20 = function(x) mean(x > 0.2)),
                              weights = weights_pred,
                              standardize_weights = TRUE,
                              col_names = "shapeName")

# Plot point predictions of average prevalence (weighted)
plot(pred_area_w, which_target = "prev", which_summary = "mean",
     palette = "RdYlGn",
     limits = c(0.1, 0.30),
     breaks =  seq(0.1, 0.30, by = 0.05)) +
  guides(fill=guide_legend(title="Prevalence")) +
  ggtitle("Average prevalence \n (population weighted)") +
  theme(plot.title = element_text(size = 15))

# Plot of the exceedance probabilities (weighted)
plot(pred_area_w, which_target = "prev", which_summary = "exceed20",
     palette = "RdYlGn",
     limits = c(0, 1),
     breaks = seq(0,1, by = 0.1)) +
  guides(fill=guide_legend(title="Probability")) +
  ggtitle("Exceedance probability (L = 0.2) \n (population weighted)") +
  theme(plot.title = element_text(size = 15))

