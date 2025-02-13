rm(list = ls())

library(RiskMap)
library(lme4)
data("liberia")

# Create the ID of the location
liberia$ID_loc <- 1:nrow(liberia)

# Binomial mixed model with log-elevation
fit_glmer_lib <- glmer(cbind(npos, ntest) ~ log(elevation) + (1|ID_loc), family = binomial,
                       data = liberia, nAGQ = 25)

summary(fit_glmer_lib)

### Poropose UTM
library(sf)
liberia_sf <- st_as_sf(liberia, coords = c("long", "lat"), crs = 4326)
crs_lb <- propose_utm(liberia_sf)

### Getting the shapefile for Liberia
library(rgeoboundaries)
shp_lb <- geoboundaries(country = "liberia", adm_lvl = "adm0")



### Getting elevation raster for Liberia
