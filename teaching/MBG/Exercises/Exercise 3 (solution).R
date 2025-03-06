rm(list = ls())

library(RiskMap)
library(lme4)

data("tz_malaria")
tz_malaria <- st_as_sf(tz_malaria, coords = c("utm_x", "utm_y"), crs = 32736)

# Fitting a model with covariates
glgm_fit <- glgpm(Pf ~ Temperature + pmax(Temperature - 33, 0) +
                    EVI + gp(nugget = NULL),
                  den = Ex,
                  data = tz_malaria,
                  family = "binomial")

# Shapefile of Tanzania
library(rgeoboundaries)
shp_tz_0 <- geoboundaries(country = "tanzania", adm_lvl = "adm0")
shp_tz_0 <- st_transform(shp_tz_0, crs = 32736)

shp_tz_1 <- geoboundaries(country = "tanzania", adm_lvl = "adm1")
shp_tz_1 <- st_transform(shp_tz_1, crs = 32736)

# Create the grid
grid_tza <- create_grid(shp_tz_0, spat_res = 10)


# EVI for Tanzania
tza_evi <- rast("teaching/MBG/Data/Tanzania_Annual_EVI_2015.tif")
tza_evi <- project(tza_evi, "EPSG:32736")

# Temperature for Tanzania
tza_temp <- rast("teaching/MBG/Data/Tanzania_Annual_LST_2015.tif")
tza_temp <- project(tza_temp, "EPSG:32736")

# Population density for Tanzania
library(terra)
tza_pop <- rast("teaching/MBG/Data/Tanzania_population_density_2015.tif")
tza_pop <- project(tza_pop, "EPSG:32736")

# Extract covariates over the grid
pred_tza <- data.frame(EVI = extract(tza_evi, st_as_sf(grid_tza))[,2],
                       Temperature = extract(tza_temp, st_as_sf(grid_tza))[,2])


# Prediction of S(x)

pred_S_tza_cov <-
  pred_over_grid(glgm_fit, predictors = pred_tza,
                 grid_pred = grid_tza, type="joint")

# Prediction of prevalence over the grid

pred_T_grid <-
  pred_target_grid(pred_S_tza_cov,
                  f_target = list(prev = function(x) exp(x)/(1+exp(x))),
                  pd_summary = list(mean = mean,
                                    below10 = function(x) mean(x < 0.05),
                                    btw10_30 = function(x) mean(x > 0.1 &
                                                                x < 0.3),
                                    class = function(x) {
                                      v1 <- mean(x < 0.05)
                                      v2 <- mean(x > 0.1 &
                                                  x < 0.3)
                                      v3 <- mean(x > 0.3)
                                      (1:3)[which.max(c(v1,v2,v3))]
                                    }))

plot(pred_T_grid, which_target = "prev", which_summary = "mean",
     main = "Predictive mean")

plot(pred_T_grid, which_target = "prev", which_summary = "below10")

plot(pred_T_grid, which_target = "prev", which_summary = "btw10_30")

plot(pred_T_grid, which_target = "prev", which_summary = "class")


pred_T_adm1 <-
  pred_target_shp(pred_S_tza_cov,
                   f_target = list(prev = function(x) exp(x)/(1+exp(x))),
                   shp = shp_tz_1,
                   shp_target = mean,
                   pd_summary = list(mean = mean,
                                     below10 = function(x) mean(x < 0.05),
                                     btw10_30 = function(x) mean(x > 0.1 &
                                                                   x < 0.2),
                                     class = function(x) {
                                       v1 <- mean(x < 0.05)
                                       v2 <- mean(x > 0.1 &
                                                    x < 0.2)
                                       v3 <- mean(x > 0.2)
                                       (1:3)[which.max(c(v1,v2,v3))]
                                     }),
                  col_names = "shapeName")

plot(pred_T_adm1, which_target = "prev", which_summary = "mean",
     palette = "RdYlGn",
     limits = c(0, 0.30),
     breaks =  seq(0, 0.30, by = 0.05)) +
  guides(fill=guide_legend(title="Prevalence")) +
  ggtitle("Average prevalence") +
  theme(plot.title = element_text(size = 15))

plot(pred_T_adm1, which_target = "prev", which_summary = "class",
     palette = "RdYlGn",
     breaks = 1:3) +
  guides(fill=guide_legend(title="Class")) +
  ggtitle("Predicted district class") +
  theme(plot.title = element_text(size = 15))


## Estimating the total number of malaria cases

tza_weights <- 10*extract(tza_pop, pred_S_tza_cov$grid_pred)[,2]
pred_T_adm1 <-
  pred_target_shp(pred_S_tza_cov,
                  f_target = list(prev = function(x) exp(x)/(1+exp(x))),
                  shp = shp_tz_0,
                  weights = tza_weights,
                  standardize_weights = FALSE,
                  shp_target = sum,
                  pd_summary = list(mean = mean,
                                    lower = function(x) quantile(x, 0.025),
                                    upper = function(x) quantile(x, 0.975)),
                  col_names = "shapeName")

pred_T_adm1$target$Tanzania$prev



################
################
data("anopheles")
anopheles_sf <- st_as_sf(anopheles, coords = c("web_x", "web_y"), crs = 3857)

an_fit <- glgpm(An.gambiae ~ elevation + gp(),
                family = "poisson",
                data = anopheles_sf)

# Create grid from convex hull
shp_ch <- convex_hull_sf(an_fit$data_sf)
an_grid <- create_grid(shp_ch, spat_res = 2)


an_elev <- get_elev_point(st_as_sf(an_grid), prj = 3857, src = "aws")$elevation


predictors <- data.frame(elevation= an_elev)

pred_an_S <- pred_over_grid(an_fit, grid = an_grid,
                            predictors = predictors,
                            type = "joint")

pred_n_mosq_grid <-
  pred_target_grid(pred_an_S,
                   f_target = list(n_mosq = function(lp) exp(lp)),
                   pd_summary = list(mean = function(Tx) mean(Tx)))

an_weights <- 1*(predictors$elevation > 400 & predictors$elevation < 830)


pred_n_mosq_shp <-
  pred_target_shp(pred_an_S, shp = shp_ch,
                  weights = an_weights,
                  shp_target = sum,
                  f_target = list(n_mosq = function(lp) exp(lp)),
                  pd_summary = list(mean = function(Tx) mean(Tx),
                                    q025 = function(Tx) quantile(Tx, 0.025),
                                    q075 = function(Tx) quantile(Tx, 0.975)))

pred_n_mosq_shp$target$reg1$n_mosq
