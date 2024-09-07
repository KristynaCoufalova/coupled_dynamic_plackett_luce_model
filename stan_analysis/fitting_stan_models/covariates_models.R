library(rstan)
library(ggplot2)
library(dplyr)
library(bayesplot)
library(tibble)
library(purrr)
library(lubridate)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

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


###### FITTING MODELS - SQUARED AGE EFFECT#####

fit_covariates_ab <- stan(
  file = 'covariates_ab_corrected.stan',
  data = stan_data,
  iter = 30000,
  chains = 2,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

saveRDS(fit_covariates_ab, file = "Fit_covariates_ab.rds")





###### FITTING MODELS - LINEAR AGE EFFECT #####

fit_non_sq_covariates_ab <- stan(
  file = 'non_sq_covariates_ab.stan',
  data = stan_data,
  iter = 30000,
  chains = 2,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

saveRDS(fit_non_sq_covariates_ab, file = "Fit_non_sq_covariates_ab.rds")

