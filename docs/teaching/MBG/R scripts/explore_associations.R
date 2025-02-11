rm(list = ls())  # Clear workspace
library(RiskMap)  # Load RiskMap package for spatial modeling

# Load Required Packages
library(ggplot2)  # For data visualization
library(sf)  # For spatial data handling
library(RiskMap)  # Spatial risk mapping functions

# Load data-sets
# These datasets contain spatial and epidemiological data
data(liberia)
data(malkenya)
data(anopheles)

# Exploring associations with risk factors using count data

# Calculate prevalence as proportion of positive cases
tibble::glimpse(liberia)  # View dataset structure
liberia$prev <- liberia$npos/liberia$ntest

# Plot prevalence against elevation
ggplot(liberia, aes(x = elevation, y = prev)) + geom_point() +
  labs(x="Elevation (meters)",y="Prevalence")

# Compute the empirical logit transformation to stabilize variance
liberia$elogit <- log((liberia$npos+0.5)/
                        (liberia$ntest-liberia$npos+0.5))

# Visualizing empirical logit transformation
ggplot(liberia, aes(x = elevation, y = elogit)) + geom_point() +

  # Adding a smoothing spline to explore trends
  labs(x="Elevation (meters)",y="Empirical logit") +
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+

  # Adding linear regression fit with log-transformed elevation
  stat_smooth(method = "lm", formula = y ~ log(x),
              col="green",lty="dashed",se=FALSE) +

  # Adding linear regression fit with change point in 150 meters
  stat_smooth(method = "lm", formula = y ~ x + pmax(x-150, 0),
              col="red",lty="dashed",se=FALSE)

# Anopheles mosquito count transformation and visualization
anopheles$log_counts <- log(anopheles$An.gambiae)
ggplot(anopheles, aes(x = elevation, y = log_counts)) + geom_point() +

  # Adding a smoothing spline to explore trends
  labs(x="Elevation (meters)",y="Log number of An. gambiae mosquitoes") +
  stat_smooth(method = "lm", formula = y ~ x, se=FALSE)

# Community survey data
# Extract subset for community survey
malkenya_comm <- malkenya[malkenya$Survey=="community", ]

# Age-group Analysis
# Group ages into predefined classes
tibble::glimpse(malkenya_comm)  # View dataset structure
malkenya_comm$Age_class <- cut(malkenya_comm$Age,
                               breaks = c(0, 5, 10, 15, 30, 40, 50, 100),
                               include.lowest = TRUE)

# Compute empirical logit within each age group
tibble::glimpse(malkenya_comm)
age_class_data <- aggregate(RDT ~ Age_class + Gender,
                            data = malkenya_comm,
                            FUN = function(y)
                              log((sum(y)+0.5)/(length(y)-sum(y)+0.5)))

# Compute the average age within each group
age_class_data$age_mean_point <- aggregate(Age ~ Age_class + Gender,
                                           data = malkenya_comm,
                                           FUN = mean)$Age

# Compute number of individuals per age group by gender
age_class_data$n_obs <-  aggregate(Age ~ Age_class + Gender,
                                   data = malkenya_comm,
                                   FUN = length)$Age

# Visualizing empirical logit by age groups
ggplot(age_class_data, aes(x = age_mean_point, y = RDT,
                           size = n_obs,
                           colour = Gender)) +
  geom_point() +
  labs(x="Age (years)",y="Empirical logit")

# Logistic regression models to analyze age and gender interactions
glm_age_gender_interaction <- glm(RDT ~ Age + Gender:Age +
                                    pmax(Age-15, 0) + Gender:pmax(Age-15, 0) +
                                    pmax(Age-40, 0) + Gender:pmax(Age-40, 0),
                                  data = malkenya_comm, family = binomial)

summary(glm_age_gender_interaction)

# Alternative model without interaction terms
glm_age_gender_no_interaction <- glm(RDT ~ Age +  pmax(Age-15, 0) + pmax(Age-40, 0),
                                     data = malkenya_comm, family = binomial)

# Compare models using Chi-square test
anova(glm_age_gender_no_interaction, glm_age_gender_interaction, test = "Chisq")

# Elevation Analysis
# Group elevation into quantiles to analyze relationship with RDT
tibble::glimpse(malkenya_comm)
malkenya_comm$elevation_class <- cut(malkenya_comm$elevation,
                                     breaks = quantile(malkenya_comm$elevation, seq(0, 1, by = 0.1)),
                                     include.lowest = TRUE)

# Compute empirical logit by elevation class
elev_class_data <- aggregate(RDT ~ elevation_class,
                             data = malkenya_comm,
                             FUN = function(y)
                               log((sum(y)+0.5)/(length(y)-sum(y)+0.5)))

# Compute average elevation within each class
elev_class_data$elevation_mean <- aggregate(elevation ~ elevation_class,
                                            data = malkenya_comm,
                                            FUN = mean)$elevation

# Visualizing empirical logit by elevation
ggplot(elev_class_data, aes(x = elevation_mean, y = RDT),
       size = n_obs) +
  geom_point() +
  labs(x="Elevation (meters)",y="Empirical logit")
