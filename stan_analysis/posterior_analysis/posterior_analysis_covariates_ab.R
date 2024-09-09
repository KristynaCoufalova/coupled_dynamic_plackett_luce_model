library(rstan)
library(ggplot2)
library(dplyr)
library(bayesplot)
library(tibble)
library(purrr)
library(tidybayes)
library(loo)
library(tidyr)
library(gridExtra)
library(xtable)
library(lubridate)

# Load fitted model
fit_ab <- readRDS(file = "fit_covariates_ab.rds")

#### Data loading ####

# Read and preprocess data
data <- read.csv('data_men.csv')
data$event_type <- factor(data$event_type)
data$competitor_id <- as.integer(factor(data$competitor_id))
data$nationality_group <- as.numeric(factor(data$nationality))
data$event_date <- as.Date(data$event_date)
data$birth_date <- dmy(data$birth_date)

data <- data %>%
  mutate(
    age = interval(birth_date, event_date) / years(1),
    age_int = as.integer(age)
  )

# Separate marathon and half-marathon data
marathon_data <- filter(data, event_type == "Men's Marathon")
half_marathon_data <- filter(data, event_type == "Men's Half Marathon")

# Function to prepare event data
prepare_event_data <- function(event_data) {
  event_data %>%
    group_by(event_id) %>%
    summarize(
      time_period = first(time_index),
      n_competitors = n(),
      competitor_ids = list(competitor_id),
      places = list(place),
      ages = list(age_int)
    ) %>%
    ungroup()
}

marathon_prepared <- prepare_event_data(marathon_data)
half_marathon_prepared <- prepare_event_data(half_marathon_data)

max_competitors <- max(c(marathon_prepared$n_competitors, half_marathon_prepared$n_competitors))

# Function to pad lists
pad_list <- function(x, n) {
  c(x, rep(0, n - length(x)))
}

# Pad data
marathon_prepared <- marathon_prepared %>%
  mutate(
    competitor_ids = map(competitor_ids, ~pad_list(.x, max_competitors)),
    places = map(places, ~pad_list(.x, max_competitors)),
    ages = map(ages, ~pad_list(.x, max_competitors))
  )

half_marathon_prepared <- half_marathon_prepared %>%
  mutate(
    competitor_ids = map(competitor_ids, ~pad_list(.x, max_competitors)),
    places = map(places, ~pad_list(.x, max_competitors)),
    ages = map(ages, ~pad_list(.x, max_competitors))
  )

# Calculate first appearance for each competitor
first_appearance <- data %>%
  group_by(competitor_id) %>%
  summarize(first_time = min(time_index))


# Prepare age data
M <- max(data$competitor_id)
T <- max(data$time_index)



# Prepare Stan data
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
  B = 1,
  age_marathon = do.call(rbind, marathon_prepared$ages),
  age_half_marathon = do.call(rbind, half_marathon_prepared$ages)
)







#### Divergences check ####
rstan::check_hmc_diagnostics(fit_ab)
sampler_params <- rstan::get_sampler_params(fit_ab, inc_warmup = FALSE)
divergences <- sapply(sampler_params, function(x) sum(x[,"divergent__"]))
print(divergences)

##### Plotting thetas and ethas #####
plot_theta_eta_evolution <- function(fit, competitor_id) {
  theta_samples <- rstan::extract(fit)$theta
  eta_samples <- rstan::extract(fit)$eta
  
  theta_comp <- theta_samples[, , competitor_id]
  eta_comp <- eta_samples[, , competitor_id]
  
  theta_stats <- data.frame(
    Time = 1:ncol(theta_comp),
    Mean = apply(theta_comp, 2, mean),
    CI_5th = apply(theta_comp, 2, quantile, probs = 0.05),
    CI_95th = apply(theta_comp, 2, quantile, probs = 0.95)
  )
  
  eta_stats <- data.frame(
    Time = 1:ncol(eta_comp),
    Mean = apply(eta_comp, 2, mean),
    CI_5th = apply(eta_comp, 2, quantile, probs = 0.05),
    CI_95th = apply(eta_comp, 2, quantile, probs = 0.95)
  )
  
  p_theta <- ggplot(theta_stats, aes(x = Time)) +
    geom_line(aes(y = Mean), color = "blue") +
    geom_ribbon(aes(ymin = CI_5th, ymax = CI_95th), fill = "blue", alpha = 0.2) +
    labs(title = paste("Evolution of Theta for Competitor", competitor_id),
         x = "Time", y = "Theta") +
    theme_minimal()
  
  p_eta <- ggplot(eta_stats, aes(x = Time)) +
    geom_line(aes(y = Mean), color = "red") +
    geom_ribbon(aes(ymin = CI_5th, ymax = CI_95th), fill = "red", alpha = 0.2) +
    labs(title = paste("Evolution of Eta for Competitor", competitor_id),
         x = "Time", y = "Eta") +
    theme_minimal()
  
  grid.arrange(p_theta, p_eta, ncol = 2)
}

# Example usage
plot_theta_eta_evolution(fit_ab, competitor_id = 1)
plot_theta_eta_evolution(fit_ab, competitor_id = 20)
plot_theta_eta_evolution(fit_ab, competitor_id = 67)

##### MODEL CHECKING #####
params_to_check <- c("a", "b", "beta_age", "beta_age_sq", "gamma_age", "gamma_age_sq", "sigma", "theta[23,20]", "eta[23,20]", 
                     "theta[23,12]", "eta[23,12]")

traceplot <- mcmc_trace(fit_ab, pars = params_to_check)
print(traceplot)

autocorr <- mcmc_acf(fit_ab, pars = params_to_check)
print(autocorr)

fit_summary <- rstan::summary(fit_ab)$summary
n_eff <- fit_summary[, "n_eff"]
print(n_eff[params_to_check])

rhat <- rhat(fit_ab)
print(summary(rhat))
hist(rhat, main="Histogram of R-hat values", xlab="R-hat")

y_rep_marathon <- as.matrix(fit_ab, pars = "y_rep_marathon")
ppc_marathon <- ppc_dens_overlay(y = stan_data$place_marathon[,1], yrep = y_rep_marathon[1:100,])
print(ppc_marathon)

y_rep_half_marathon <- as.matrix(fit_ab, pars = "y_rep_half_marathon")
ppc_half_marathon <- ppc_dens_overlay(y = stan_data$place_half_marathon[,1], yrep = y_rep_half_marathon[1:100,])
print(ppc_half_marathon)

nuts_params <- rstan::get_sampler_params(fit_ab, inc_warmup = FALSE)
nuts_df <- do.call(rbind, nuts_params) %>% as.data.frame()
n_chains <- length(nuts_params)
n_iterations <- nrow(nuts_df) / n_chains

nuts_df <- nuts_df %>%
  mutate(Chain = rep(1:n_chains, each = n_iterations),
         Iteration = rep(1:n_iterations, times = n_chains)) %>%
  pivot_longer(cols = -c(Chain, Iteration),
               names_to = "Parameter",
               values_to = "Value")
energy <- mcmc_nuts_energy(nuts_df)
print(energy)

mcmc_dens_overlay(fit_ab, pars = params_to_check) +
  ggtitle("Density Plots for Selected Parameters")

#### Exploration of a and b ####
posterior_samples <- rstan::extract(fit_ab, pars = c("a", "b"))
a_samples <- posterior_samples$a
b_samples <- posterior_samples$b

posterior_df <- data.frame(a = a_samples, b = b_samples)

ggplot(posterior_df, aes(x = a, y = b)) +
  geom_point(alpha = 0.1, color = "blue") +
  geom_density2d(color = "red") +
  labs(title = "Joint Posterior Distribution of a and b",
       x = "Parameter a",
       y = "Parameter b") +
  theme_minimal()

cor_ab <- cor(a_samples, b_samples)
print(cor_ab)

print(mean(a_samples))
print(mean(b_samples))

print(median(a_samples))
print(median(b_samples))

quantiles_a <- quantile(a_samples, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantiles_b <- quantile(b_samples, c(0.05, 0.25, 0.5, 0.75, 0.95))
print(quantiles_a)
print(quantiles_b)

#### Printing thetas and ethas for the last period ####
theta_samples <- rstan::extract(fit_ab)$theta
eta_samples <- rstan::extract(fit_ab)$eta

theta_stats <- apply(theta_samples[, 26, ], 2, function(x) {
  c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))
})
eta_stats <- apply(eta_samples[, 26, ], 2, function(x) {
  c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))
})

theta_df <- as.data.frame(t(theta_stats))
eta_df <- as.data.frame(t(eta_stats))

theta_df <- theta_df %>%
  rownames_to_column(var = "Competitor") %>%
  mutate(Parameter = "Theta")

eta_df <- eta_df %>%
  rownames_to_column(var = "Competitor") %>%
  mutate(Parameter = "Eta")

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

theta_df <- theta_df %>%
  mutate(name = data$name[match(Competitor, data$competitor_id)]) %>%
  select(name, Competitor, Mean, `CI Lower`, `CI Upper`)

eta_df <- eta_df %>%
  mutate(name = data$name[match(Competitor, data$competitor_id)]) %>%
  select(name, Competitor, Mean, `CI Lower`, `CI Upper`)

latex_table_theta <- xtable(theta_df, caption = "Posterior Means and 95\\% Credible Intervals for Theta  (t = 26)")
latex_table_eta <- xtable(eta_df, caption = "Posterior Means and 95\\% Credible Intervals for Eta (t = 26)")

print(latex_table_theta, type = "latex", include.rownames = FALSE, caption.placement = "top", booktabs = TRUE)
print(latex_table_eta, type = "latex", include.rownames = FALSE, caption.placement = "top", booktabs = TRUE)


##### CHECKING THE AGE EFFECT #####

# Assuming 'age' was included as a covariate in your Stan model,
# the following code will help you explore its posterior distribution and effect.

# Extract posterior samples for the age coefficient
age_effect_marathon_samples <- rstan::extract(fit_ab)$beta_age  
age_effect_squared_marathon_samples <- rstan::extract(fit_ab)$beta_age_sq  
age_effect_half_marathon_samples <- rstan::extract(fit_ab)$gamma_age 
age_effect_squared_half_marathon_samples <- rstan::extract(fit_ab)$gamma_age_sq


# Function to summarize and plot age effects
analyze_age_effect <- function(samples, title) {
  mean_effect <- mean(samples)
  median_effect <- median(samples)
  quantiles_effect <- quantile(samples, c(0.05, 0.25, 0.5, 0.75, 0.95))
  
  # Print summary statistics
  cat(paste0(title, ":\n"))
  cat("  Mean: ", mean_effect, "\n")
  cat("  Median: ", median_effect, "\n")
  cat("  Quantiles: \n")
  print(quantiles_effect)
  
  # Plot the posterior distribution
  ggplot(data.frame(effect = samples), aes(x = effect)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = paste("Posterior Distribution of", title),
         x = paste(title, "Coefficient"),
         y = "Density") +
    theme_minimal()
}

# Analyze and plot the age effects
# Marathon age effect
analyze_age_effect(age_effect_marathon_samples, "Marathon Age Effect")

# Marathon age squared effect
analyze_age_effect(age_effect_squared_marathon_samples, "Marathon Age Squared Effect")

# Half-marathon age effect
analyze_age_effect(age_effect_half_marathon_samples, "Half-Marathon Age Effect")

# Half-marathon age squared effect
analyze_age_effect(age_effect_squared_half_marathon_samples, "Half-Marathon Age Squared Effect")

# Joint plot for marathon age and age squared effects
posterior_df_marathon <- data.frame(
  age_effect = age_effect_marathon_samples, 
  age_squared_effect = age_effect_squared_marathon_samples
)

ggplot(posterior_df_marathon, aes(x = age_effect, y = age_squared_effect)) +
  geom_point(alpha = 0.1, color = "blue") +
  geom_density2d(color = "red") +
  labs(title = "Joint Posterior Distribution of Marathon Age and Age Squared Effects",
       x = "Age Effect Coefficient",
       y = "Age Squared Effect Coefficient") +
  theme_minimal()

# Joint plot for half-marathon age and age squared effects
posterior_df_half_marathon <- data.frame(
  age_effect = age_effect_half_marathon_samples, 
  age_squared_effect = age_effect_squared_half_marathon_samples
)

ggplot(posterior_df_half_marathon, aes(x = age_effect, y = age_squared_effect)) +
  geom_point(alpha = 0.1, color = "purple") +
  geom_density2d(color = "red") +
  labs(title = "Joint Posterior Distribution of Half-Marathon Age and Age Squared Effects",
       x = "Age Effect Coefficient",
       y = "Age Squared Effect Coefficient") +
  theme_minimal()

# Correlation between age effect and age squared effect for marathon
cor_marathon <- cor(age_effect_marathon_samples, age_effect_squared_marathon_samples)
cat("Correlation between Marathon Age Effect and Age Squared Effect: ", cor_marathon, "\n")

# Correlation between age effect and age squared effect for half-marathon
cor_half_marathon <- cor(age_effect_half_marathon_samples, age_effect_squared_half_marathon_samples)
cat("Correlation between Half-Marathon Age Effect and Age Squared Effect: ", cor_half_marathon, "\n")

##### FIT LINEAR COVARIATES #####

fit_ns_ab <- readRDS(file = "fit_non_sq_covariates_ab.rds")



##### MODEL COMPARISON #####

# Extract log-likelihoods
log_lik1 <- rstan::extract(fit_ab, pars = "log_lik")$log_lik
log_lik2 <- rstan::extract(fit_ns_ab, pars = "log_lik")$log_lik




# Compute LOO-CV
loo1 <- loo(log_lik1)
loo2 <- loo(log_lik2)



# Print LOO-CV results
print(loo1)
print(loo2)


# Compare the two models using LOO-CV
loo_comparison <- loo_compare(loo1, loo2)
print(loo_comparison)

# Compute WAIC
waic1 <- waic(log_lik1)
waic2 <- waic(log_lik2)



# Print WAIC results
print(waic1)
print(waic2)


# Compare the two models using WAIC
waic_comparison <- loo_compare(waic1, waic2)
print(waic_comparison)


