# Point 1

library(ggplot2)

# Scatter plot for elogit vs Temperature
p1 <- ggplot(tz_malaria, aes(x = Temperature, y = elogit)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "blue") +
  labs(title = "Empirical logit vs Temperature",
       x = "Temperature (°C)",
       y = "Empirical logit") +
  theme_minimal()

# Scatter plot for elogit vs Precipitation
p2 <- ggplot(tz_malaria, aes(x = Precipitation, y = elogit)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") +
  labs(title = "Empirical logit vs Precipitation",
       x = "Precipitation (mm)",
       y = "Empirical logit") +
  theme_minimal()

# Scatter plot for elogit vs EVI
p3 <- ggplot(tz_malaria, aes(x = EVI, y = elogit)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "green") +
  labs(title = "Empirical logit vs EVI",
       x = "EVI",
       y = "Empirical logit") +
  theme_minimal()

# Scatter plot for elogit vs NTL
p4 <- ggplot(tz_malaria, aes(x = log(NTL), y = elogit)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "purple") +
  labs(title = "Empirical logit vs NTL",
       x = "log-NTL",
       y = "Empirical logit") +
  theme_minimal()


# Point 2
library(boot)

# Fit the binomial GLM
fit_bin <- glm(cbind(Pf, Ex - Pf) ~ Temperature + I(Temperature^2),
               data = tz_malaria,
               family = binomial)

# Function to calculate peak temperature
peak_temp <- function(data) {
  fit <- glm(cbind(Pf, Ex - Pf) ~ Temperature + I(Temperature^2),
             data = data,
             family = binomial)
  coef_fit <- coef(fit)
  a <- coef_fit["I(Temperature^2)"]
  b <- coef_fit["Temperature"]
  peak_temp <- -b / (2 * a)
  return(peak_temp)
}

# Parametric bootstrapping
parametric_bootstrap <- function(fit, estimator, nsim = 1000) {
  original_data <- fit$data
  out <- numeric(nsim)
  for (i in 1:nsim) {
    simulated_response <- simulate(fit)$sim_1
    new_data <- original_data
    new_data$Pf <- simulated_response[,1]
    out[i] <- estimator(new_data)
  }
  return(out)
}

# Run bootstrapping
set.seed(123)
peak_temps <- parametric_bootstrap(fit_bin, estimator = peak_temp, nsim = 1000)

# Calculate peak temperature and CI
peak_temp_estimate <- peak_temp(tz_malaria)
ci_lower <- quantile(peak_temps, 0.025)
ci_upper <- quantile(peak_temps, 0.975)

cat("Temperature at which prevalence is highest:", peak_temp_estimate, "\n")
cat("95% Confidence Interval: [", ci_lower, ",", ci_upper, "]\n")

# Point 3

# Estimate change point
estimate_change_point <- function(data, change_points = seq(25, 35, by = 0.1)) {
  log_likelihoods <- numeric(length(change_points))
  for (i in seq_along(change_points)) {
    change_point <- change_points[i]
    fit <- glm(cbind(Pf, Ex - Pf) ~ Temperature + pmax(Temperature - change_point, 0),
               data = data,
               family = binomial)
    log_likelihoods[i] <- logLik(fit)
  }
  optimal_change_point <- change_points[which.max(log_likelihoods)]
  return(optimal_change_point)
}

opt_cp <- estimate_change_point(tz_malaria)

# Fit linear spline model
fit_bin_sp <- glm(cbind(Pf, Ex - Pf) ~ Temperature + pmax(Temperature - opt_cp, 0),
                  data = tz_malaria,
                  family = binomial)

# Bootstrapping for change point
cp_boot <- parametric_bootstrap(fit = fit_bin_sp,
                                estimator = estimate_change_point,
                                nsim = 100)

ci_lower_cp <- quantile(cp_boot, 0.025)
ci_upper_cp <- quantile(cp_boot, 0.975)
cat("Estimate of the change point in temperature:", opt_cp, "\n")
cat("95% Confidence Interval: [", ci_lower_cp, ",", ci_upper_cp, "]\n")

# Plot empirical logit with GAM and linear spline
tz_malaria$predicted_sp <- predict(fit_bin_sp)
p_spl <- ggplot(tz_malaria, aes(x = Temperature, y = elogit)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "blue") +
  geom_line(aes(y = predicted_sp), color = "red") +
  labs(title = "Empirical logit vs Temperature",
       x = "Temperature (°C)",
       y = "Empirical logit") +
  theme_minimal()

# Point 4

library(lme4)

# Fit GLMM
fit2_glmer <- glmer(cbind(Pf, Ex - Pf) ~ Temperature + pmax(Temperature - opt_cp, 0) +
                      EVI + (1 | cluster.number),
                    data = tz_malaria, nAGQ = 100,
                    family = binomial)

summary(fit2_glmer)

# Point 5

library(sf)
library(rgeoboundaries)
library(terra)

# Convert to sf object
tz_sf <- st_as_sf(tz_malaria, coords = c("utm_x", "utm_y"), crs = 32736)

# Download Tanzania boundary
shp_tz <- geoboundaries(country = "tanzania", adm_lvl = "adm0")
shp_tz <- st_transform(shp_tz, crs = 32736)

# Create a 5 km grid
grid_utm <- create_grid(shp_tz, spat_res = 5)

# Read and project raster files
r_temp <- rast("path/to/local/Tanzania_Annual_LST_2015.tif")
r_temp <- terra::project(r_temp, "EPSG:32736")
r_evi <- rast("path/to/local/Tanzania_Annual_EVI_2015.tif")
r_evi <- terra::project(r_evi, "EPSG:32736")

# Extract values and predict prevalence
pred_data_frame <- data.frame(
  Temperature = extract(r_temp, st_coordinates(grid_utm))[,1],
  EVI = extract(r_evi, st_coordinates(grid_utm))[,1]
)

pred_glmer <- predict(fit2_glmer, newdata = pred_data_frame, type = "response", re.form = NA)

# Create raster plot
raster_pred <- data.frame(
  x = st_coordinates(grid_utm)[, 1],
  y = st_coordinates(grid_utm)[, 2],
  prev = pred_glmer
)

ggplot(data = raster_pred) +
  geom_raster(aes(x = x, y = y, fill = prev)) +
  scale_fill_viridis_c() +
  coord_cartesian() +
  theme_minimal() +
  labs(title = "Predictions", fill = "Prevalence")
