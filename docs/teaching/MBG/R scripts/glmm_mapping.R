# Clear the workspace
rm(list = ls())

# Load necessary libraries
library(RiskMap)  # For RiskMap functionality
library(lme4)    # For mixed-effects models
data("liberia")   # Load the Liberia dataset

# Create a unique ID for each location in the Liberia dataset
liberia$ID_loc <- 1:nrow(liberia)

# Fit a binomial mixed-effects model with log-elevation as a fixed effect and location as a random effect
fit_glmer_lib <- glmer(cbind(npos, ntest) ~ log(elevation) + (1|ID_loc),
                       family = binomial, data = liberia, nAGQ = 25)

# Summarize the model results
summary(fit_glmer_lib)

### Propose UTM projection for Liberia
library(sf)
# Convert Liberia data to an sf object with WGS84 (EPSG:4326) CRS
liberia_sf <- st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
# Propose a UTM CRS for Liberia
crs_lb <- propose_utm(liberia_sf)

### Get the shapefile for Liberia
library(rgeoboundaries)
# Download Liberia's administrative boundary (ADM0 level)
shp_lb <- geoboundaries(country = "liberia", adm_lvl = "adm0")
# Transform the shapefile to the proposed UTM CRS
shp_lb <- st_transform(shp_lb, crs = crs_lb)
# Create a grid over Liberia with a spatial resolution of 5 km
grid_utm <- create_grid(shp_lb, spat_res = 5)
# Transform the grid back to WGS84 (EPSG:4326) for compatibility with other data
grid_ll <- st_transform(grid_utm, crs = 4326)

### Get elevation raster for Liberia
library(elevatr)
# Download elevation data for Liberia at zoom level 5 and clip it to the country boundary
elev_lb <- get_elev_raster(locations = shp_lb, z = 5, clip = "locations")

# Convert the elevation raster to a data frame for use with ggplot2
library(terra)
raster_df <- as.data.frame(elev_lb, xy = TRUE)
colnames(raster_df)[3] <- "elev"  # Rename the elevation column

# Plot the elevation raster using ggplot2
library(ggplot2)
ggplot(data = raster_df) +
  geom_raster(aes(x = x, y = y, fill = elev)) +  # Plot elevation as a raster
  scale_fill_viridis_c() +  # Use viridis color scale
  coord_cartesian() +  # Use Cartesian coordinates
  theme_minimal() +  # Use a minimal theme
  labs(title = "Elevation", fill = "Density") +  # Add title and legend label
  geom_sf(data = shp_lb, col = 2, lwd = 2, fill = NA)  # Overlay Liberia boundary

### Download and process temperature data
library(geodata)
# Download monthly minimum temperature data for Liberia from WorldClim
clim_lb <- worldclim_country(country = "Liberia", var = "tmin", path = tempdir())
# Mask the temperature data to Liberia's boundary
a <- mask(clim_lb, st_transform(shp_lb, crs = 4326))
# Plot the mean minimum temperature
terra::plot(mean(a), plg = list(title = "Min. temperature (C)"))

### Predict prevalence using elevation data
# Extract elevation values at grid locations
elev_data_frame <- data.frame(elevation = extract(elev_lb, st_coordinates(grid_utm)))

# Fit a binomial GLM model with log-elevation as a predictor
glm_fit <- glm(cbind(npos, ntest - npos) ~ log(elevation),
               data = liberia, family = "binomial")

# Predict prevalence using the fitted model
pred_glm <- predict(glm_fit, newdata = elev_data_frame, type = "response")

# Create a data frame with predicted prevalence and grid coordinates
raster_pred <- data.frame(x = st_coordinates(grid_utm)[, 1],
                          y = st_coordinates(grid_utm)[, 2],
                          prev = pred_glm)

# Plot the predicted prevalence using ggplot2
ggplot(data = raster_pred) +
  geom_raster(aes(x = x, y = y, fill = prev)) +  # Plot prevalence as a raster
  scale_fill_viridis_c() +  # Use viridis color scale
  coord_cartesian() +  # Use Cartesian coordinates
  theme_minimal() +  # Use a minimal theme
  labs(title = "Predictions", fill = "Prevalence")  # Add title and legend label
