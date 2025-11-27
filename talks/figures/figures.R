# Load necessary libraries
library(ggplot2)

# Function to generate data
generate_mixture_data <- function(N, pi1, pi2, mu, sigma) {
  set.seed(123)
  c1 <- rbinom(N, 1, pi1) + 1
  c2 <- rbinom(N, 1, pi2) + 1
  y1 <- rnorm(N, mean = mu[c1, 1], sd = sigma[c1, 1])
  y2 <- rnorm(N, mean = mu[c2, 2], sd = sigma[c2, 2])
  cluster <- paste0("(", c1, ",", c2, ")")
  data.frame(y1 = y1, y2 = y2, cluster = factor(cluster))
}

# Parameters for Example 1 (balanced mixing)
mu <- matrix(c(0, 1, 5, 6), nrow = 2, byrow = TRUE)
sigma <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)
data1 <- generate_mixture_data(500, 0.5, 0.5, mu, sigma)

# Plot Example 1
ggplot(data1, aes(y1, y2, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  ggtitle("Example 1: Balanced Mixing") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("example1.png", width = 6, height = 4)

# Parameters for Example 2 (skewed mixing)
data2 <- generate_mixture_data(500, 0.3, 0.7, mu, sigma)

# Plot Example 2
ggplot(data2, aes(y1, y2, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  ggtitle("Example 2: Skewed Mixing") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("example2.png", width = 6, height = 4)