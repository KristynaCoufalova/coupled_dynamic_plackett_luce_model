library(tidyverse)
library(rstan)
library(bayesplot)
# Set options for rstan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

# Parameters
N_competitors <-80
T_periods <- 26
N_marathon <- 87
N_half_marathon <- 57
avg_competitors_per_race <- 5
a <- 0.4  # Coupling parameter
b <- 0.2  # Coupling parameter

sigma <- 1  # Scale parameter

# Create Q_M matrix
Q_M <- diag(N_competitors) - matrix(1/N_competitors, N_competitors, N_competitors)

# 1. Create competitors with assigned abilities
create_initial_abilities <- function(N) {
  epsilon_B <- rnorm(N,0,sigma)
  abilities <- Q_M %*% epsilon_B  # Ensure initial abilities sum to zero
}


evolve_abilities <- function(theta_prev, eta_prev, a,b, Q_M, sigma) {
  epsilon_1 <- rnorm(N_competitors, 0, sqrt(2*a*(1-a))*sigma)
  epsilon_2 <- rnorm(N_competitors, 0, sqrt(2*b*(1-b))*sigma)
  
  theta_new <- (1-a)* theta_prev + a * eta_prev + Q_M %*% epsilon_1
  eta_new <- (1-b) * eta_prev + b * theta_prev + Q_M %*% epsilon_2
  
  list(theta = theta_new, eta = eta_new)
}

# Initialize abilities
theta <- matrix(0, nrow = N_competitors, ncol = T_periods)
eta <- matrix(0, nrow = N_competitors, ncol = T_periods)
theta[,1] <- create_initial_abilities(N_competitors)
eta[,1] <- create_initial_abilities(N_competitors)

# Evolve abilities over time
for (t in 2:T_periods) {
  evolved <- evolve_abilities(theta[,t-1], eta[,t-1], a,b, Q_M, sigma)
  theta[,t] <- evolved$theta
  eta[,t] <- evolved$eta
}

# 2. Generate competitions (unchanged)
generate_competitions <- function(N_races, T_periods, N_competitors, avg_competitors) {
  races <- tibble(
    race_id = 1:N_races,
    time_period = sample(1:T_periods, N_races, replace = TRUE),
    n_competitors = pmin(rpois(N_races, avg_competitors) + 1, N_competitors)  # Ensure n_competitors does not exceed N_competitors
  )
  
  race_competitors <- races %>%
    mutate(competitors = map(n_competitors, ~sample(1:N_competitors, .x, replace = FALSE))) %>%
    unnest(competitors)
  
  return(race_competitors)
}

marathon_races <- generate_competitions(N_marathon, T_periods, N_competitors, avg_competitors_per_race)
half_marathon_races <- generate_competitions(N_half_marathon, T_periods, N_competitors, avg_competitors_per_race)

# 3. Simulate rankings using Plackett-Luce model (unchanged)
simulate_rankings <- function(race_data, abilities) {
  race_data %>%
    group_by(race_id) %>%
    mutate(
      ability = abilities[cbind(competitors, time_period)],
      prob = exp(ability) / sum(exp(ability)),
      place = sample(1:n(), size = n(), replace = FALSE, prob = prob)
    ) %>%
    ungroup()
}




marathon_results <- simulate_rankings(marathon_races, theta)
half_marathon_results <- simulate_rankings(half_marathon_races, eta)

# 4. Prepare data for Stan model (unchanged)

library(tidyverse)

# Function to prepare data for each event
prepare_event_data <- function(event_data) {
  event_data %>%
    group_by(race_id) %>%
    summarize(
      event_id = first(race_id),
      time_period = first(time_period),
      n_competitors = n(),
      competitor_ids = list(competitors),
      places = list(place)
    ) %>%
    ungroup()
}

# Prepare marathon and half-marathon data
marathon_prepared <- prepare_event_data(marathon_results)
half_marathon_prepared <- prepare_event_data(half_marathon_results)

# Find the maximum number of competitors in any event
max_competitors <- max(c(marathon_prepared$n_competitors, half_marathon_prepared$n_competitors))

# Function to pad lists with 0
pad_list <- function(x, n) {
  c(x, rep(0, n - length(x)))
}

# Pad the competitor_ids and places lists with 0
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
first_appearance <- bind_rows(marathon_results, half_marathon_results) %>%
  group_by(competitors) %>%
  summarize(first_time = min(time_period))



# Prepare data for Stan
stan_data <- list(
  T = T_periods,
  M = N_competitors,
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
  B=1
)

# Ensure all elements are integers where necessary
stan_data$competitor_id_marathon <- matrix(as.integer(stan_data$competitor_id_marathon), 
                                           nrow=nrow(stan_data$competitor_id_marathon))
stan_data$competitor_id_half_marathon <- matrix(as.integer(stan_data$competitor_id_half_marathon), 
                                                nrow=nrow(stan_data$competitor_id_half_marathon))
stan_data$place_marathon <- matrix(as.integer(stan_data$place_marathon), 
                                   nrow=nrow(stan_data$place_marathon))
stan_data$place_half_marathon <- matrix(as.integer(stan_data$place_half_marathon), 
                                        nrow=nrow(stan_data$place_half_marathon))


# You can now use stan_data with your Stan model
# Fit the Stan model
artificial_fit_ab <- stan(
  file = 'artificial.stan',
  data = stan_data,
  iter = 30000,
  chains = 2,
  seed = 123
)
# Save the final combined result
saveRDS(artificial_fit_ab, file = "artificial_fit_ab_no_sigma.rds")
#artificial_fit <- readRDS("artificial_fit_ab.rds")

print(artificial_fit)
library(rstan)
library(bayesplot)

# Extract posterior samples
posterior <- rstan::extract(artificial_fit_ab)

print(dim(posterior$theta))
print(dim(posterior$eta))

# Calculate posterior means for theta and eta
#theta_post_mean <- apply(posterior$theta, c(2,3), mean)
#eta_post_mean <- apply(posterior$eta, c(2,3), mean)

# Calculate posterior means for theta and eta
theta_post_mean <- rowMeans(aperm(posterior$theta, c(2, 3, 1)), dims = 2)
eta_post_mean <- rowMeans(aperm(posterior$eta, c(2, 3, 1)), dims = 2)



# True values
theta_true <- t(theta)
eta_true <- t(eta)
print(dim(theta))
print(dim(eta))

# Transpose the posterior mean matrices to match the true values
#theta_post_mean <- t(theta_post_mean)
#eta_post_mean <- t(eta_post_mean)


# Calculate the probabilities
compute_probabilities <- function(theta_matrix) {
  exp_theta <- exp(theta_matrix)  # Compute exponentials of theta
  sum_exp_theta <- rowSums(exp_theta)  # Sum of exponentials for each time period
  probabilities <- exp_theta / sum_exp_theta  # Divide by the sum to get probabilities
  
  return(probabilities)
}

# Apply the function
probabilities_theta_true <- compute_probabilities(theta_true)
probabilities_theta_post <- compute_probabilities(theta_post_mean)



# Prepare the data for plotting
time_periods <- 1:nrow(probabilities_theta_true)
data <- data.frame(
  Time = rep(time_periods, each = ncol(probabilities_theta_true)),
  Competitor = rep(1:ncol(probabilities_theta_true), times = nrow(probabilities_theta_true)),
  True_Probability = as.vector(probabilities_theta_true),
  Post_Probability = as.vector(probabilities_theta_post)
)

# Gather the data into long format for ggplot2
data_long <- data %>%
  pivot_longer(cols = c(True_Probability, Post_Probability), 
               names_to = "Probability_Type", 
               values_to = "Probability")

# Create the plot
ggplot(data_long, aes(x = Time, y = Probability, color = Probability_Type, linetype = Probability_Type)) +
  geom_line() +
  facet_wrap(~ Competitor, scales = "free_y") +
  labs(title = "True vs Posterior Probabilities Over Time",
       x = "Time Period",
       y = "Probability",
       color = "Probability Type",
       linetype = "Probability Type") +
  theme_minimal()



# Filter the data for the first 20 competitors
data_long_filtered <- data_long %>%
  filter(Competitor <= 20)

# Create the plot for the first 20 competitors
ggplot(data_long_filtered, aes(x = Time, y = Probability, color = Probability_Type, linetype = Probability_Type)) +
  geom_line() +
  facet_wrap(~ Competitor, scales = "free_y") +
  labs(title = "True vs Posterior Probabilities of Winning Marathon Race Over Time (First 20 Competitors)",
       x = "Time Period",
       y = "Probability",
       color = "Probability Type",
       linetype = "Probability Type") +
  theme_minimal()
###### Eta

# Calculate the probabilities for eta
probabilities_eta_true <- compute_probabilities(eta_true)
probabilities_eta_post <- compute_probabilities(eta_post_mean)

# Prepare the data for plotting for eta
time_periods <- 1:nrow(probabilities_eta_true)
data_eta <- data.frame(
  Time = rep(time_periods, each = ncol(probabilities_eta_true)),
  Competitor = rep(1:ncol(probabilities_eta_true), times = nrow(probabilities_eta_true)),
  True_Probability = as.vector(probabilities_eta_true),
  Post_Probability = as.vector(probabilities_eta_post)
)

# Gather the data into long format for ggplot2 for eta
data_eta_long <- data_eta %>%
  pivot_longer(cols = c(True_Probability, Post_Probability), 
               names_to = "Probability_Type", 
               values_to = "Probability")

# Filter the data for the first 20 competitors for eta
data_eta_long_filtered <- data_eta_long %>%
  filter(Competitor <= 4)

# Create the plot for the first 20 competitors for eta
ggplot(data_eta_long_filtered, aes(x = Time, y = Probability, color = Probability_Type, linetype = Probability_Type)) +
  geom_line() +
  facet_wrap(~ Competitor, scales = "free_y") +
  labs(title = "True vs Posterior Probabilities of Winning Half Marathon Over Time (First 20 Competitors for eta)",
       x = "Time Period",
       y = "Probability",
       color = "Probability Type",
       linetype = "Probability Type") +
  theme_minimal()



#####
# Compute the correlation for each competitor
correlations <- sapply(1:ncol(probabilities_theta_true), function(i) {
  cor(probabilities_theta_true[, i], probabilities_theta_post[, i])
})

# Create a data frame to store the results
correlation_results <- data.frame(
  Competitor = 1:ncol(probabilities_theta_true),
  Correlation = correlations
)

# View the results
print(correlation_results)




compute_placement_probabilities <- function(theta_matrix) {
  n_time_periods <- nrow(theta_matrix)
  n_competitors <- ncol(theta_matrix)
  
  # Initialize a list to store the probabilities for each time period
  placement_probs <- vector("list", n_time_periods)
  
  for (t in 1:n_time_periods) {
    remaining_competitors <- 1:n_competitors  # Start with all competitors
    probs <- matrix(0, nrow = n_competitors, ncol = n_competitors)  # To store placement probabilities
    
    for (j in 1:n_competitors) {
      exp_theta <- exp(theta_matrix[t, remaining_competitors])
      sum_exp_theta <- sum(exp_theta)
      probs[j, remaining_competitors] <- exp_theta / sum_exp_theta
      
      # Select the competitor with the highest probability for the j-th place
      max_prob_index <- which.max(probs[j, remaining_competitors])
      
      # Remove the selected competitor from the remaining competitors
      remaining_competitors <- remaining_competitors[-max_prob_index]
    }
    
    placement_probs[[t]] <- probs
  }
  
  return(placement_probs)
}

# Apply the corrected function to compute placement probabilities
placement_probs_true <- compute_placement_probabilities(theta_true)

# Example of how to access the probabilities for a specific time period
placement_probs_true[[4]]  # Access probabilities for time period 1


#Placement probabilities for time period 1:
#1st place: Competitor 1: 0.5, Competitor 2: 0.3, Competitor 3: 0.2
#2nd place: Competitor 1: 0.4, Competitor 2: 0.4, Competitor 3: 0.2
#3rd place: Competitor 1: 0.6, Competitor 2: 0.3, Competitor 3: 0.1





a_simulated <- posterior$a
b_simulated <- posterior$b


print(mean(a_simulated))
print(mean(b_simulated))




sigma_simulated <- posterior$sigma


print(mean(sigma_simulated))
print(median(sigma_simulated))




params_to_check <-c("a", "b", "theta[23,7]", "eta[23,7]", "theta[26,7]", "eta[26,7]","theta[20,55]", "eta[20,12]", "theta[25,55]", "eta[25,55]", "theta[22,71]")
# 2. Effective sample size (ESS)
fit_summary <- rstan::summary(artificial_fit_ab)$summary

# Now you should be able to extract n_eff
n_eff <- fit_summary[, "n_eff"]



# Print the effective sample sizes for the parameters of interest
print(n_eff[params_to_check])

# Trace plots
traceplot <- mcmc_trace(artificial_fit_ab, pars = params_to_check)
print(traceplot)


###### WITH SIGMA ####

artificial_fit_ab_sigma <- stan(
  file = 'stan_ab.stan',
  data = stan_data,
  iter = 30000,
  chains = 2,
  seed = 123
)

#saveRDS(artificial_fit_ab_sigma, file = "artificial_fit_ab_sigma.rds")


# Extract posterior samples
posterior <- rstan::extract(artificial_fit_ab_sigma)

print(dim(posterior$theta))
print(dim(posterior$eta))

# Calculate posterior means for theta and eta
#theta_post_mean <- apply(posterior$theta, c(2,3), mean)
#eta_post_mean <- apply(posterior$eta, c(2,3), mean)

# Calculate posterior means for theta and eta
theta_post_mean <- rowMeans(aperm(posterior$theta, c(2, 3, 1)), dims = 2)
eta_post_mean <- rowMeans(aperm(posterior$eta, c(2, 3, 1)), dims = 2)



# True values
theta_true <- t(theta)
eta_true <- t(eta)
print(dim(theta))
print(dim(eta))

# Transpose the posterior mean matrices to match the true values
#theta_post_mean <- t(theta_post_mean)
#eta_post_mean <- t(eta_post_mean)


# Calculate the probabilities
compute_probabilities <- function(theta_matrix) {
  exp_theta <- exp(theta_matrix)  # Compute exponentials of theta
  sum_exp_theta <- rowSums(exp_theta)  # Sum of exponentials for each time period
  probabilities <- exp_theta / sum_exp_theta  # Divide by the sum to get probabilities
  
  return(probabilities)
}

# Apply the function
probabilities_theta_true <- compute_probabilities(theta_true)
probabilities_theta_post <- compute_probabilities(theta_post_mean)



# Prepare the data for plotting
time_periods <- 1:nrow(probabilities_theta_true)
data <- data.frame(
  Time = rep(time_periods, each = ncol(probabilities_theta_true)),
  Competitor = rep(1:ncol(probabilities_theta_true), times = nrow(probabilities_theta_true)),
  True_Probability = as.vector(probabilities_theta_true),
  Post_Probability = as.vector(probabilities_theta_post)
)

# Gather the data into long format for ggplot2
data_long <- data %>%
  pivot_longer(cols = c(True_Probability, Post_Probability), 
               names_to = "Probability_Type", 
               values_to = "Probability")

# Create the plot
ggplot(data_long, aes(x = Time, y = Probability, color = Probability_Type, linetype = Probability_Type)) +
  geom_line() +
  facet_wrap(~ Competitor, scales = "free_y") +
  labs(title = "True vs Posterior Probabilities Over Time",
       x = "Time Period",
       y = "Probability",
       color = "Probability Type",
       linetype = "Probability Type") +
  theme_minimal()



# Filter the data for the first 20 competitors
data_long_filtered <- data_long %>%
  filter(Competitor <= 4)

# Create the plot for the first 20 competitors
ggplot(data_long_filtered, aes(x = Time, y = Probability, color = Probability_Type, linetype = Probability_Type)) +
  geom_line() +
  facet_wrap(~ Competitor, scales = "free_y") +
  labs(title = "True vs Posterior Probabilities of Winning Marathon Race Over Time (First 4 Competitors)",
       x = "Time Period",
       y = "Probability",
       color = "Probability Type",
       linetype = "Probability Type") +
  theme_minimal()



params_to_check <-c("a", "b","sigma", "theta[23,7]", "eta[23,7]", "theta[26,7]", "eta[26,7]","theta[20,55]", "eta[20,12]", "theta[25,55]", "eta[25,55]", "theta[22,71]")
# 2. Effective sample size (ESS)
fit_summary <- rstan::summary(artificial_fit_ab_sigma)$summary

# Now you should be able to extract n_eff
n_eff <- fit_summary[, "n_eff"]


# Print the effective sample sizes for the parameters of interest
print(n_eff[params_to_check])

# Trace plots
traceplot <- mcmc_trace(artificial_fit_ab_sigma, pars = params_to_check)
print(traceplot)
