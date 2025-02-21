library(shiny)
library(ggplot2)
library(gstat)
library(sp)
library(gridExtra)
library(sf)

# Matern correlation function
matern_corr <- function(h, kappa, phi) {
  if (h == 0) return(1)
  factor <- (2^(1 - kappa)) / gamma(kappa)
  scaled_h <- h / phi
  factor * (scaled_h^kappa) * besselK(scaled_h, kappa)
}

# UI
ui <- fluidPage(
  titlePanel("Matérn GP Simulator"),
  sidebarLayout(
    sidebarPanel(
      numericInput("phi", "Phi (scale parameter):", 0.1, min = 0.01, max = 1, step = 0.01),
      numericInput("kappa", "Kappa (smoothness):", 2.5, min = 0.5, max = 5, step = 0.1),
      numericInput("sigma2", "Sigma^2 (variance):", 1, min = 0.1, max = 5, step = 0.1),
      actionButton("simulate", "Generate Plots")
    ),
    mainPanel(
      plotOutput("maternPlot"),
      plotOutput("simulationPlot")
    )
  )
)

# Server function
server <- function(input, output) {
  observeEvent(input$simulate, {
    h_vals <- seq(0, sqrt(2), length.out = 100)
    corr_vals <- sapply(h_vals, matern_corr, kappa = input$kappa, phi = input$phi)

    output$maternPlot <- renderPlot({
      ggplot(data.frame(h = h_vals, corr = corr_vals), aes(x = h, y = corr)) +
        geom_line() +
        labs(title = "Matérn Correlation Function", x = "Distance h", y = "Correlation")
    })

    # Simulating a Spatial Field
    grid_size <- 50
    x <- seq(0, 1, length.out = grid_size)
    y <- seq(0, 1, length.out = grid_size)
    grid <- expand.grid(x = x, y = y)
    coordinates(grid) <- ~x + y
    grid <- st_as_sf(grid)
    st_crs(grid) <- 3857

    # Define the variogram model (Matérn covariance)
    sim_res <-
      glgpm_sim(n_sim = 1,
                formula = ~ gp(kappa = input$kappa),
                data = grid,
                family = "gaussian",
                crs = 3857,
                sim_pars = list(beta = 0, sigma2 = input$sigma2,
                                phi = input$phi,
                                tau2 = 0, sigma2_me = 0),
                scale_to_km = FALSE)

    # Convert the result to a data frame for ggplot
    sim_df <- data.frame(x = st_coordinates(grid)[,1],
                         y = st_coordinates(grid)[,2],
                         value = sim_res$data_sim$`_sim1`)

    output$simulationPlot <- renderPlot({
      ggplot(sim_df, aes(x = x, y = y, fill = value)) +
        geom_tile() +
        scale_fill_viridis_c() +
        theme_minimal() +
        labs(title = "Simulated Matérn GP Surface", fill = "Value")
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
