library(rstan)
library(ggplot2)
library(dplyr)
library(bayesplot)
library(tibble)
library(purrr)
library(tidybayes)
library(loo)
library(tidyr)


# Set options for rstan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)


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


##### FITTIG MODEL ab####
# Fit the Stan model
fit_ab <- stan(
  file = 'stan_ab.stan',
  data = stan_data,
  iter = 30000,
  chains = 2,
  seed = 123
)

saveRDS(fit_ab, file = "fit_ab.rds")


##### FITTIG MODELS abcd ####
# Fit the Stan model
fit_abcd <- stan(
  file = 'stan_abcd.stan',
  data = stan_data,
  iter = 30000,
  chains = 2,
  seed = 123
)

saveRDS(fit_abcd, file = "fit_abcd.rds")
