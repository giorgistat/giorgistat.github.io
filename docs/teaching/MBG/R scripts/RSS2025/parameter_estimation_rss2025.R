# 1) FIT A BINOMIAL GEOSTATISTICAL MODEL ---------------------------------------

lb_fit <- glgpm(npos ~ log(elevation) + gp(kappa = 0.5, nugget = 0),
                den = ntest,
                family = "binomial",
                data = liberia_utm)

summary(lb_fit)

# Update MCML until convergence

while(lb_fit$log.lik > 1) {

  # Extract estimated parameters
  theta0 <- coef(lb_fit)

  # Re-fit the model
  lb_fit <- glgpm(npos ~ log(elevation) + gp(),
                  den = ntest,
                  par0 = theta0,
                  family = "binomial",
                  data = liberia_utm,
                  messages = F)

  print(lb_fit$log.lik)
}

summary(lb_fit)

# Visualise estimated correlation function
udist <- seq(0, 250, l = 1000)
corr <- exp(-udist / coef(lb_fit)$phi)

df <- data.frame(udist, corr)

ggplot(df) +
  geom_line(aes(x = udist, y = corr)) +
  labs(x = "Distance (km)", y = "Estimated correlation") +
  geom_vline(xintercept = coef(lb_fit)$phi * 3, linetype = 2, col = "red")

# 2) FIT A PROBIT MODEL --------------------------------------------------------

# Using an alternative link function
lb_fit_probit <- glgpm(npos ~ log(elevation) + gp(),
                       den = ntest,
                       invlink = function(x) pnorm(x),
                       family = "binomial",
                       data = liberia_sf)

summary(lb_fit_probit)
