# 1) LOAD R PACKAGES AND DATA --------------------------------------------------

# Load R packages
library(ggplot2)        # For data visualization
library(sf)             # For spatial data handling
library(RiskMap)        # Spatial risk mapping functions
library(rgeoboundaries) # Retrieve administrative boundaries
library(mapview)        # Interactive viz of spatial data
library(ggspatial)      # Tools for spatial viz
library(lme4)           # Fit random effects model

# Load Liberia Onchocerciasis dataset
data(liberia)

# 2) EXPLORATORY DATA ANALYSIS -------------------------------------------------

## 2.1) DATA VISUALISATION -----------------------------------------------------

# View dataset structure
tibble::glimpse(liberia)

# Compute prevalence
liberia$prev <- liberia$npos/liberia$ntest

# Convert to sf object before visualisation
liberia_sf <- st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)

# Quick viz
plot(liberia_sf$geometry)
plot(liberia_sf["prev"])

# Quick interactive viz with mapview
mapview(liberia_sf)
mapview(liberia_sf, zcol = "prev")

# Map with ggplot
map <- ggplot(liberia_sf) +
  geom_sf(aes(fill = prev * 100), shape = 21, size = 2.5, col = "black") +
  scale_fill_distiller("Prevalence (%)", type = "div", palette = 5) +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(legend.position = "top", legend.key.width = unit(1.5, ( "cm")))

map

# Download administrative boundaries for Liberia (level 0: country)
liberia_admin0 <- gb_adm0("Liberia")

# Add country border
map_border <- map + geom_sf(data = liberia_admin0, col = "black", fill = NA)

map_border

# Add scale bar and north arrow
map_border +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr")

## 2.2) EXPLORE ASSOCIATIONS ---------------------------------------------------

# Exploring associations with risk factors using count data

# Plot prevalence against elevation
ggplot(liberia, aes(x = elevation, y = prev)) +
  geom_point() +
  labs(x = "Elevation (meters)", y = "Prevalence")

# Compute the empirical logit transformation to stabilize variance
liberia$elogit <- log((liberia$npos + 0.5) / (liberia$ntest - liberia$npos + 0.5))

# Visualizing empirical logit transformation
ggplot(liberia, aes(x = elevation, y = elogit)) +
  geom_point() +
  # Adding a smoothing spline to explore trends
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE)+
  # Adding linear regression fit with log-transformed elevation
  stat_smooth(method = "lm", formula = y ~ log(x),
              col = "green", lty = "dashed", se = FALSE) +
  # Adding linear regression fit with change point in 150 meters
  stat_smooth(method = "lm", formula = y ~ x + pmax(x-150, 0),
              col = "red", lty = "dashed", se = FALSE) +
  labs(x="Elevation (meters)", y = "Empirical logit")

# 3) FIT NON-SPATIAL MODEL -----------------------------------------------------

# Create the ID of the location
liberia$ID_loc <- 1:nrow(liberia)

# Binomial mixed model with log-elevation
fit_glmer_lib <- glmer(cbind(npos, ntest) ~ log(elevation) + (1|ID_loc),
                       family = binomial,
                       data = liberia,
                       nAGQ = 25)

summary(fit_glmer_lib)

# 4) CHECK RESIDUAL SPATIAL VARIATION ------------------------------------------

# Extract estimated random effects Z_i
liberia_sf$Z_hat <- ranef(fit_glmer_lib)$ID_loc[,1]

# When computing distances is better to convert to metric system
crs_utm_lb <- propose_utm(liberia_sf)
liberia_utm <- st_transform(liberia_sf, crs = crs_utm_lb)

# Compute variogram
liberia_variog <- s_variogram(data = liberia_utm,
                              variable = "Z_hat",
                              bins = seq(15, 250, l = 10),
                              scale_to_km = TRUE,
                              n_permutation = 10000)

# Visualise variogram
plot_s_variogram(liberia_variog)

# Add envelope
plot_s_variogram(liberia_variog, plot_envelope = TRUE)
