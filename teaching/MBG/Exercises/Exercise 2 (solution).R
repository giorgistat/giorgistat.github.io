rm(list = ls())

library(RiskMap)
library(lme4)
library(sf)

data("tz_malaria")
tz_malaria <- st_as_sf(tz_malaria, coords = c("utm_x", "utm_y"), crs = 32736)


# Point 1

glmer_fit_null <- glmer(cbind(Pf, Ex - Pf) ~ (1 | cluster.number),
                   data = tz_malaria, nAGQ = 100,
                   family = binomial)

tz_malaria$Z_hat_null <- ranef(glmer_fit_null)$cluster.number[,1]


tz_vari_null <-
  s_variogram(data = tz_malaria, variable = "Z_hat_null",
              bins = seq(1,500, length = 15),
              scale_to_km = TRUE,
              n_permutation = 1000)

plot_s_variogram(tz_vari_null, plot_envelope = TRUE)


# Point 2

glmer_fit <- glmer(cbind(Pf, Ex - Pf) ~ Temperature + pmax(Temperature - 33, 0) +
                      EVI + (1 | cluster.number),
                    data = tz_malaria, nAGQ = 100,
                    family = binomial)

tz_malaria$Z_hat <- ranef(glmer_fit)$cluster.number[,1]

tz_vari <-
s_variogram(data = tz_malaria, variable = "Z_hat",
            bins = seq(1,500, length = 15),
            scale_to_km = TRUE,
            n_permutation = 1000)

plot_s_variogram(tz_vari, plot_envelope = TRUE)




### Point 3

# Fitting an intercept only model
glgm_fit_null <- glgpm(Pf ~ gp(nugget = NULL),
                       den = Ex,
                      data = tz_malaria,
                      family = "binomial")

summary(glgm_fit_null)

glgm_fit_null <- glgpm(Pf ~ gp(nugget = NULL),
                       den = Ex,
                       par0 = coef(glgm_fit_null),
                       data = tz_malaria,
                       family = "binomial")


# Fitting a model with covariates
glgm_fit <- glgpm(Pf ~ Temperature + pmax(Temperature - 33, 0) +
                    EVI + gp(nugget = NULL),
                       den = Ex,
                       data = tz_malaria,
                       family = "binomial")

summary(glgm_fit)


### Point 4

library(ggplot2)
library(dplyr)

# Extract parameters
theta_hat_null <- coef(glgm_fit_null)
theta_hat <- coef(glgm_fit)

# Define a sequence of distances
distances <- seq(0, 1000, length.out = 200)

# Compute variogram values for both models
variogram_df <- data.frame(
  Distance = rep(distances, 2),
  Variogram = c(
    theta_hat_null$tau2 + theta_hat_null$sigma2 * (1 - exp(-distances / theta_hat_null$phi)),
    theta_hat$tau2 + theta_hat$sigma2 * (1 - exp(-distances / theta_hat$phi))
  ),
  Model = factor(rep(c("No covariates", "With covariates"), each = length(distances)))
)

# Plot using ggplot2
ggplot(variogram_df, aes(x = Distance, y = Variogram, color = Model)) +
  geom_line(size = 1.2) +
  labs(
    title = "Theoretical Variogram Comparison",
    x = "Distance",
    y = "Variogram"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  theme(
    legend.title = element_blank(), # Removes legend title
    legend.position = "bottom" # Places legend at the bottom for clarity
  )

