library(ggplot2) 
library(dplyr)
library(bayesplot)
library(tibble)
library(purrr)
library(tidybayes)
library(loo)
library(tidyr)
library(gridExtra)


##### DATA CLEANING ######

# Load and prepare data
data <- read.csv('data_men.csv')
data$event_type <- factor(data$event_type)
data$competitor_id <- as.integer(factor(data$competitor_id))

# Create nationality groups
data$nationality_group <- as.numeric(factor(data$nationality))



# Separate marathon and half-marathon data
marathon_data <- filter(data, event_type == "Men's Marathon")
half_marathon_data <- filter(data, event_type == "Men's Half Marathon")

# Function to prepare data for each event
prepare_event_data <- function(event_data) {
  event_data %>%
    group_by(event_id) %>%
    summarize(
      time_period = first(time_index),
      n_competitors = n(),
      competitor_ids = list(competitor_id),
      places = list(place)
    ) %>%
    ungroup()
}

marathon_prepared <- prepare_event_data(marathon_data)
half_marathon_prepared <- prepare_event_data(half_marathon_data)

# Find the maximum number of competitors in any event
max_competitors <- max(c(marathon_prepared$n_competitors, half_marathon_prepared$n_competitors))

# Pad the competitor_ids and places lists with 0 (instead of NA)
pad_list <- function(x, n) {
  c(x, rep(0, n - length(x)))
}

marathon_prepared <- marathon_prepared %>%
  mutate(
    competitor_ids = map(competitor_ids, ~pad_list(.x, max_competitors)),
    places = map(places, ~pad_list(.x, max_competitors))
  )

half_marathon_prepared <- half_marathon_prepared %>%
  mutate(
    competitor_ids = map(competitor_ids, ~pad_list(.x, max_competitors)),
    places = map(places, ~pad_list(.x, max_competitors))
  )

# Get first appearance time for each competitor
first_appearance <- data %>%
  group_by(competitor_id) %>%
  summarize(first_time = min(time_index))


##### PREPARING STAN DATA #####

# Prepare data for Stan
stan_data <- list(
  T = max(data$time_index),
  M = max(data$competitor_id),
  N_marathon = nrow(marathon_prepared),
  N_half_marathon = nrow(half_marathon_prepared),
  max_competitors = max_competitors,
  time_period_marathon = marathon_prepared$time_period,
  time_period_half_marathon = half_marathon_prepared$time_period,
  n_competitors_marathon = marathon_prepared$n_competitors,
  n_competitors_half_marathon = half_marathon_prepared$n_competitors,
  competitor_id_marathon = do.call(rbind, marathon_prepared$competitor_ids),
  competitor_id_half_marathon = do.call(rbind, half_marathon_prepared$competitor_ids),
  place_marathon = do.call(rbind, marathon_prepared$places),
  place_half_marathon = do.call(rbind, half_marathon_prepared$places),
  first_appearance = first_appearance$first_time,
  B = 1  # Assuming time index starts at 1
)

# Ensure all elements are integers
stan_data$competitor_id_marathon <- matrix(as.integer(stan_data$competitor_id_marathon), 
                                           nrow=nrow(stan_data$competitor_id_marathon))
stan_data$competitor_id_half_marathon <- matrix(as.integer(stan_data$competitor_id_half_marathon), 
                                                nrow=nrow(stan_data$competitor_id_half_marathon))
stan_data$place_marathon <- matrix(as.integer(stan_data$place_marathon), 
                                   nrow=nrow(stan_data$place_marathon))
stan_data$place_half_marathon <- matrix(as.integer(stan_data$place_half_marathon), 
                                        nrow=nrow(stan_data$place_half_marathon))


#### Load fitted stan data ####

fit_ab<- readRDS("SD_w0.rds")


#### Divergences check ####

rstan::check_hmc_diagnostics(fit_ab)
sampler_params <- rstan::get_sampler_params(fit_ab, inc_warmup = FALSE)
# Check for divergences
divergences <- sapply(sampler_params, function(x) sum(x[,"divergent__"]))
print(divergences)


  
##### Plotting thetas and ethas #####

# Function to get competitor name by competitor_id
get_competitor_name <- function(competitor_id) {
  competitor_name <- data %>%
    filter(competitor_id == !!competitor_id) %>%
    pull(name) %>%
    unique()
  
  # Return the name or a placeholder if the competitor is not found
  if (length(competitor_name) > 0) {
    return(competitor_name)
  } else {
    return(paste("Competitor", competitor_id))
  }
}



plot_theta_eta_evolution <- function(fit, competitor_id) {
  
  # Extract theta and eta samples
  theta_samples <- rstan::extract(fit)$theta
  eta_samples <- rstan::extract(fit)$eta
  
  # Extract samples for the specified competitor across all time periods
  theta_comp <- theta_samples[, , competitor_id]
  eta_comp <- eta_samples[, , competitor_id]
  
  # Get the first appearance time for the competitor
  first_time <- stan_data$first_appearance[competitor_id]
  
  # Create a sequence of all time points
  time_points <- 1:ncol(theta_comp)
  
  # Create mean and credible intervals (5th and 95th percentiles) for theta and eta
  theta_mean <- apply(theta_comp, 2, mean)
  theta_ci_5th <- apply(theta_comp, 2, quantile, probs = 0.05)
  theta_ci_95th <- apply(theta_comp, 2, quantile, probs = 0.95)
  
  eta_mean <- apply(eta_comp, 2, mean)
  eta_ci_5th <- apply(eta_comp, 2, quantile, probs = 0.05)
  eta_ci_95th <- apply(eta_comp, 2, quantile, probs = 0.95)
  
  # Set values before first appearance to NA (they won't be plotted)
  theta_mean[time_points < first_time] <- NA
  theta_ci_5th[time_points < first_time] <- NA
  theta_ci_95th[time_points < first_time] <- NA
  
  eta_mean[time_points < first_time] <- NA
  eta_ci_5th[time_points < first_time] <- NA
  eta_ci_95th[time_points < first_time] <- NA
  
  # Create data frames for plotting
  theta_stats <- data.frame(
    Time = time_points,
    Mean = theta_mean,
    CI_5th = theta_ci_5th,
    CI_95th = theta_ci_95th
  )
  
  eta_stats <- data.frame(
    Time = time_points,
    Mean = eta_mean,
    CI_5th = eta_ci_5th,
    CI_95th = eta_ci_95th
  )
  
  # Get competitor's name
  competitor_name <- get_competitor_name(competitor_id)
  
  # Plot evolution of theta
  p_theta <- ggplot(theta_stats, aes(x = Time)) +
    geom_line(aes(y = Mean), color = "blue") +
    geom_ribbon(aes(ymin = CI_5th, ymax = CI_95th), fill = "blue", alpha = 0.2) +
    labs(title = paste("Evolution of Theta for Competitor", competitor_id, "(", competitor_name, ")"),
         x = "Time", y = "Theta") +
    theme_minimal()
  
  # Plot evolution of eta
  p_eta <- ggplot(eta_stats, aes(x = Time)) +
    geom_line(aes(y = Mean), color = "red") +
    geom_ribbon(aes(ymin = CI_5th, ymax = CI_95th), fill = "red", alpha = 0.2) +
    labs(title = paste("Evolution of Eta for Competitor", competitor_id, "(", competitor_name, ")"),
         x = "Time", y = "Eta") +
    theme_minimal()
  
  # Arrange the two plots side by side
  grid.arrange(p_theta, p_eta, ncol = 2)
}

# Example usage:

plot_theta_eta_evolution(fit_ab, competitor_id = 7)
plot_theta_eta_evolution(fit_ab, competitor_id = 12)
plot_theta_eta_evolution(fit_ab, competitor_id = 22)
plot_theta_eta_evolution(fit_ab, competitor_id = 23)
plot_theta_eta_evolution(fit_ab, competitor_id = 55)       
plot_theta_eta_evolution(fit_ab, competitor_id = 71)       





##### MODEL CHECKING #####

# 1. Trace plots and autocorrelation
# Choose a few key parameters to check
params_to_check <-c("a", "b", "sigma", "theta[23,7]", "eta[23,7]", "theta[26,7]", "eta[26,7]","theta[20,55]", "eta[20,12]", "theta[25,55]", "eta[25,55]", "theta[22,71]")

# Trace plots
traceplot <- mcmc_trace(fit_ab, pars = params_to_check)
print(traceplot)
#ggsave("traceplot.png", traceplot, width = 12, height = 8)

# Autocorrelation plots
autocorr <- mcmc_acf(fit_ab, pars = params_to_check)
print(autocorr)
#ggsave("autocorrelation.png", autocorr, width = 12, height = 8)

# Function to find the first lag where autocorrelation is smaller than 0.1
find_lag_below_threshold <- function(fit, params, threshold = 0.1, max_lags = 100) {
  results <- list()
  
  # Loop over each parameter
  for (param in params) {
    # Get autocorrelation values for the parameter up to max_lags
    autocorr_values <- bayesplot::mcmc_acf(fit, pars = param, lags = max_lags)
    
    # Extract autocorrelation data
    autocorr_data <- as.data.frame(autocorr_values$data)
    
    # Find the first lag where autocorrelation is smaller than threshold (0.1)
    lag_below_threshold <- autocorr_data %>%
      filter(AC < threshold) %>%
      slice(1) %>%
      pull(Lag)
    
    # If no lag is found, return "No lag found" message
    if (length(lag_below_threshold) == 0) {
      lag_below_threshold <- "No lag found below threshold"
    }
    
    # Store the result for the parameter
    results[[param]] <- lag_below_threshold
  }
  
  return(results)
}


# Call the function to find lags where autocorrelation is smaller than 0.1
lags_below_threshold <- find_lag_below_threshold(fit_ab, params_to_check, threshold = 0.1)

# Print the results
print(lags_below_threshold)




# 2. Effective sample size (ESS)
fit_summary <- rstan::summary(fit_ab)$summary

# Now you should be able to extract n_eff
n_eff <- fit_summary[, "n_eff"]

# Print the effective sample sizes for the parameters of interest
print(n_eff[params_to_check])

# 3. Gelman-Rubin statistic (R-hat)
rhat <- rhat(fit_ab)
print(summary(rhat))
hist(rhat, main="Histogram of R-hat values", xlab="R-hat")



# 4. Posterior predictive checks
# For marathon
y_rep_marathon <- as.matrix(fit_ab, pars = "y_rep_marathon")
ppc_marathon <- ppc_dens_overlay(y = stan_data$place_marathon[,1], 
                                 yrep = y_rep_marathon[1:100,])
print(ppc_marathon)
#ggsave("ppc_marathon.png", ppc_marathon, width = 10, height = 6)

# For half-marathon
y_rep_half_marathon <- as.matrix(fit_ab, pars = "y_rep_half_marathon")
ppc_half_marathon <- ppc_dens_overlay(y = stan_data$place_half_marathon[,1], 
                                      yrep = y_rep_half_marathon[1:100,])
print(ppc_half_marathon)
#ggsave("ppc_half_marathon.png", ppc_half_marathon, width = 10, height = 6)

# 5. Energy plots

nuts_params <- rstan::get_sampler_params(fit_ab, inc_warmup = FALSE)
nuts_df <- do.call(rbind, nuts_params) %>% as.data.frame()
# Add chain and iteration information
n_chains <- length(nuts_params)
n_iterations <- nrow(nuts_df) / n_chains

# Create a new data frame in the required format
nuts_df <- nuts_df %>%
  mutate(Chain = rep(1:n_chains, each = n_iterations),
         Iteration = rep(1:n_iterations, times = n_chains)) %>%
  pivot_longer(cols = -c(Chain, Iteration),
               names_to = "Parameter",
               values_to = "Value")
print(nuts_df)
# Validate the new data frame
if (!all(c("Chain", "Iteration", "Parameter", "Value") %in% names(nuts_df))) {
  stop("NUTS parameters are not in the expected format.")
}
energy <- mcmc_nuts_energy(nuts_df)
print(energy)
#ggsave("energy_plot.png", energy, width = 10, height = 6)


# 6. Density plots
mcmc_dens_overlay(fit_ab, pars = params_to_check) +
  ggtitle("Density Plots for Selected Parameters")


#### Exploration of a and b ####

# Assuming `fit_abcd` is your combined Stan fit object
posterior_samples <- rstan::extract(fit_ab, pars = c("a", "b"))

# Extracting samples for 'a' and 'b'
a_samples <- posterior_samples$a
b_samples <- posterior_samples$b


# Create a data frame with the samples
posterior_df <- data.frame(a = a_samples, b = b_samples)

# Scatter plot with density contours
ggplot(posterior_df, aes(x = a, y = b)) +
  geom_point(alpha = 0.1, color = "blue") +  # Scatter plot for joint distribution
  geom_density2d(color = "red") +            # Add density contours
  labs(title = "Joint Posterior Distribution of a and b",
       x = "Parameter a",
       y = "Parameter b") +
  theme_minimal()


# Calculate posterior correlation
cor_ab <- cor(a_samples, b_samples)
print(cor_ab)

# Posterior summary of a and b

print(mean(a_samples ))
print(mean(b_samples))

print(median(a_samples ))
print(median(b_samples))

# Calculate quantiles
quantiles_a <- quantile(a_samples, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantiles_b <- quantile(b_samples, c(0.05, 0.25, 0.5, 0.75, 0.95))
print(quantiles_a)
print(quantiles_b)



#### Printing thetas and ethas for the last period ####

library(dplyr)
library(tidyr)
library(xtable)
library(rstan)

# Extract posterior samples for theta and eta
theta_samples <- rstan::extract(fit_ab)$theta
eta_samples <- rstan::extract(fit_ab)$eta

# Calculate mean and 95% CI for theta and eta for period t = 26
theta_stats <- apply(theta_samples[, 26, ], 2, function(x) {
  c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))
})
eta_stats <- apply(eta_samples[, 26, ], 2, function(x) {
  c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))
})

# Combine into a single data frame
theta_df <- as.data.frame(t(theta_stats))
eta_df <- as.data.frame(t(eta_stats))

# Add labels for easier identification
theta_df <- theta_df %>%
  rownames_to_column(var = "Competitor") %>%
  mutate(Parameter = "Theta")

eta_df <- eta_df %>%
  rownames_to_column(var = "Competitor") %>%
  mutate(Parameter = "Eta")

# Rename columns for clarity
theta_df <- theta_df %>%
  rename(
    Mean = mean,
    `CI Lower` = `2.5%`,
    `CI Upper` = `97.5%`
  ) %>%
  select(Competitor, Mean, `CI Lower`, `CI Upper`)

eta_df <- eta_df %>%
  rename(
    Mean = mean,
    `CI Lower` = `2.5%`,
    `CI Upper` = `97.5%`
  ) %>%
  select(Competitor, Mean, `CI Lower`, `CI Upper`)

# Add athlete names by searching in the data dataframe
theta_df <- theta_df %>%
  mutate(name = data$name[match(Competitor, data$competitor_id)]) %>%
  select(name, Competitor, Mean, `CI Lower`, `CI Upper`)

eta_df <- eta_df %>%
  mutate(name = data$name[match(Competitor, data$competitor_id)]) %>%
  select(name, Competitor, Mean, `CI Lower`, `CI Upper`)


# Create LaTeX table using xtable
latex_table_theta <- xtable(theta_df, caption = "Posterior Means and 95\\% Credible Intervals for Theta  (t = 26)")
latex_table_eta <- xtable(eta_df, caption = "Posterior Means and 95\\% Credible Intervals for Eta (t = 26)")

# Print LaTeX table to console
print(latex_table_theta, type = "latex", include.rownames = FALSE, caption.placement = "top", booktabs = TRUE)
print(latex_table_eta, type = "latex", include.rownames = FALSE, caption.placement = "top", booktabs = TRUE)
