# 1) LOAD R PACKAGES -----------------------------------------------------------

# Load R packages
library(wpgpDownloadR) # Download population data from worldpop
library(elevatr)       # Download elevation data
library(terra)         # Raster handling

# 2) CREATE PREDICTION GRID ----------------------------------------------------

# Convert the admin boundaries to the same CRS of our data
liberia_admin0 <- st_transform(liberia_admin0, crs = st_crs(liberia_utm))

# Create the grid (5km  resolution)
grid_lb <- create_grid(liberia_admin0, spat_res = 5)

# Visualise the grid
plot(liberia_admin0$geometry)
plot(grid_lb, cex = .5, add = T, col = "grey")

# 3) EXTRACT COVARIATE VALUES --------------------------------------------------

# Extract elevation
elevation <- get_elev_point(st_as_sf(grid_lb),
                            prj = st_crs(grid_lb),
                            src = "aws")$elevation

pred_lb <- as.data.frame(elevation)

# 4) PIXEL LEVEL PREDICTIONS ---------------------------------------------------

# Predict S(x) (marginal predictions)
pred_S_lb <- pred_over_grid(lb_fit, grid_pred = grid_lb, predictors = pred_lb)

# Predict targes of interest T(x)
pred_T_lb <- pred_target_grid(
  pred_S_lb,
  f_target = list(prev = function(x) exp(x) / (1 + exp(x))),
  pd_summary = list(mean = mean,
                    lower_lim = function(x) quantile(x, 0.025),
                    upper_lim = function(x) quantile(x, 0.975),
                    exceed20 = function(x) mean(x > 0.2),
                    exceed30 = function(x) mean(x > 0.3))
  )

# Visualise the predictions
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

# 4) AREA LEVEL PREDICTIONS (UNWEIGHTED) ---------------------------------------

# Retrieve admin 1 boundaries for liberia
liberia_admin1 <- gb_adm1("Liberia")

# Change CRS
liberia_admin1 <- st_transform(liberia_admin1, crs = st_crs(liberia_admin0))

plot(liberia_admin1$geometry)

# Predict S(x) (joint predictions)
pred_S_lb <- pred_over_grid(lb_fit, grid_pred = grid_lb, predictors = pred_lb,
                            type = "joint")


# Unweighted regional prevalence
pred_area <- pred_target_shp(
  pred_S_lb,
  shp = liberia_admin1,
  shp_target = function(x) mean(x),
  f_target = list(prev = function(x) exp(x) / (1 + exp(x))),
  pd_summary = list(mean = mean,
                    exceed20 = function(x) mean(x > 0.2)),
  col_names = "shapeName"
  )

# Plot point predictions of average prevalence
ggplot(pred_area$shp) +
  geom_sf(aes(fill = prev_mean * 100), col = "black") +
  scale_fill_distiller("Average Prevalence (%)", palette = "YlOrRd", direction = 1) +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(legend.position = "top", legend.key.width = unit(1.5, ( "cm")))

# Plot of the exceedance probabilities
ggplot(pred_area$shp) +
  geom_sf(aes(fill = prev_exceed20), col = "black") +
  scale_fill_gradient2("Exceedance probability (L = 20%)",
                       midpoint = 0.5, low = "blue4", high = "red4") +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(legend.position = "top", legend.key.width = unit(1.5, ( "cm")))

# 5) AREA LEVEL PREDICTIONS (WEIGHTED) -----------------------------------------

# Obtaining population

# Search for datasets available for Liberia
# usign the ISO3 country code
lbr_datasets <- wpgpListCountryDatasets(ISO3 = "LBR")

# Download population and load as a raster
lbr_url <- wpgpGetCountryDataset(ISO3 = "LBR", covariate = "ppp_2014")
lbr_pop <- rast(lbr_url)

# Aggregate raster at the same resolution of our grid
lbr_pop_5km <- aggregate(lbr_pop, fact = 50, fun = "sum", na.rm = T)

# Extra pop weights at the prediction grid
grid_lb_wgs84 <- st_transform(st_as_sf(grid_lb), crs = 4326)
weights_pred <- extract(lbr_pop, grid_lb_wgs84)$lbr_ppp_2014

# Prediction of the population weighted regional average prevalence
pred_area_w <- pred_target_shp(
  pred_S_lb,
  shp = liberia_admin1,
  shp_target = function(x) sum(x),
  f_target = list(prev = function(x) exp(x) / (1 + exp(x))),
  pd_summary = list(mean = mean,
                    exceed20 = function(x) mean(x > 0.2)),
  weights = weights_pred,
  standardize_weights = TRUE,
  col_names = "shapeName")

# Plot point predictions of average prevalence (weighted)
ggplot(pred_area_w$shp) +
  geom_sf(aes(fill = prev_mean * 100), col = "black") +
  scale_fill_distiller("Average Prevalence (%)", palette = "YlOrRd", direction = 1) +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(legend.position = "top", legend.key.width = unit(1.5, ( "cm")))

# Plot of the exceedance probabilities (weighted)
ggplot(pred_area_w$shp) +
  geom_sf(aes(fill = prev_exceed20), col = "black") +
  scale_fill_gradient2("Exceedance probability (L = 20%)",
                       midpoint = 0.5, low = "blue4", high = "red4") +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(legend.position = "top", legend.key.width = unit(1.5, ( "cm")))
