rm(list = ls())

library(RiskMap)
library(sf)
library(ggplot2)

data("galicia")

# Convert to an sf object
galicia_sf <- st_as_sf(galicia, coords = c("x", "y"), crs = 32629)

ggplot(data = galicia_sf) +
  geom_sf(aes(color = lead, size = lead)) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "",
       color = "Lead conc.")

lm_fit <- lm(log(lead) ~ 1, data = galicia_sf)
galicia_sf$residuals <- lm_fit$residuals

# Compute the variogram, using the residuals from the linear model fit,
# and the 95% confidence level envelope for spatial independence
galicia_variog <- s_variogram(galicia_sf, variable = "residuals",
                                scale_to_km = TRUE,
                                bins = seq(10, 140, length = 15),
                                n_permutation = 10000)

# Plotting the results
plot_s_variogram(galicia_variog, plot_envelope = TRUE)

# Fit a linear geostatistical model
fit_galicia <-
  glgpm(log(lead) ~ gp(x, y, kappa = 1.5), data=galicia, family = "gaussian",
        crs = 32629, scale_to_km = TRUE, messages = FALSE)

summary(fit_galicia)

# Maximum likelihood estimates
par_hat <- coef(fit_galicia)

# Given parameters
sigma2 <- par_hat$sigma2
phi <- par_hat$phi

# Compute the practical range where correlation = 0.05
practical_range <- -phi * log(0.05)

# Define distance values for plotting
h <- seq(0, practical_range * 1.2, length.out = 100)
variogram <- sigma2 * (1 - exp(-h / phi))

# Plot the variogram
plot(h, variogram, type = "l", col = "blue", lwd = 2,
     xlab = "Distance (h)", ylab = "Variogram γ(h)",
     main = "Theoretical Variogram (Exponential Model)")

# Add vertical line for practical range
abline(v = practical_range, col = "red", lty = 2)

# Add legend
legend("topleft", legend = c("Theoretical Variogram", "Practical Range"),
       col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1))
