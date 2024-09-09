
//

functions {
  real plackett_luce_lpmf(array[] int place, array[] int competitor_id, vector ability) {
    int n = size(place);
    real log_prob = 0;
    array[n] int valid_indices;
    int valid_count = 0;
    
    for (i in 1:n) {
      if (place[i] != 0 && competitor_id[i] != 0) {  // Assuming 0 represents NA
        valid_count += 1;
        valid_indices[valid_count] = i;
      }
    }
    
    array[valid_count] int order = sort_indices_asc(place[valid_indices[1:valid_count]]);
    
    for (i in 1:(valid_count-1)) {
      int idx = valid_indices[order[i]];
      vector[valid_count - i + 1] subability;
      for (j in i:valid_count) {
        subability[j-i+1] = ability[competitor_id[valid_indices[order[j]]]];
      }
      log_prob += ability[competitor_id[idx]] - log_sum_exp(subability);
    }
    return log_prob;
  }
  
  array[] int plackett_luce_rng(vector ability) {
    int n = num_elements(ability);
    array[n] int rank;
    vector[n] softmax_probs = softmax(ability);
    
    for (i in 1:n) {
      int chosen = categorical_rng(softmax_probs);
      rank[i] = chosen;
      softmax_probs[chosen] = 0;
      softmax_probs = softmax_probs / sum(softmax_probs);
    }
    return rank;
  }
}

data {
  int<lower=0> T;  // Number of time periods
  int<lower=1> M;  // Total number of competitors
  int<lower=0> N_marathon;  // Number of marathon events
  int<lower=0> N_half_marathon;  // Number of half-marathon events
  int<lower=0> max_competitors;  // Maximum number of competitors in any event
  array[N_marathon] int<lower=1, upper=T> time_period_marathon;
  array[N_half_marathon] int<lower=1, upper=T> time_period_half_marathon;
  array[N_marathon] int<lower=1, upper=max_competitors> n_competitors_marathon;
  array[N_half_marathon] int<lower=1, upper=max_competitors> n_competitors_half_marathon;
  array[N_marathon, max_competitors] int competitor_id_marathon;
  array[N_half_marathon, max_competitors] int competitor_id_half_marathon;
  array[N_marathon, max_competitors] int place_marathon;
  array[N_half_marathon, max_competitors] int place_half_marathon;
  array[M] int<lower=1, upper=T> first_appearance;  // First appearance time for each competitor
  int<lower=1> B;  // Initial time period
}

parameters {
  matrix[T, M] theta_raw;  // Raw marathon abilities
  matrix[T, M] eta_raw;    // Raw half-marathon abilities
  real<lower=0, upper=1> a;  // Coupling parameter
  real<lower=0, upper=1> b;  // Coupling parameter
}

transformed parameters {
  matrix[T, M] theta;  // Marathon abilities
  matrix[T, M] eta;    // Half-marathon abilities
  matrix[M, M] Q_M = diag_matrix(rep_vector(1, M)) - rep_matrix(1.0 / M, M, M);
  
  
  // Initialize at first appearance time for each competitor
  for (m in 1:M) {
    int t = first_appearance[m];
    theta[t, m] = dot_product(Q_M[m], theta_raw[t]);
    eta[t, m] = dot_product(Q_M[m], eta_raw[t]);
    
    // Set abilities to 0 before first appearance
    // Note: These zeros don't influence model results as they're not used in likelihood
    // Note: If we modeled the whole period for all athletes, the model would be overparameterized
    // Note: We could leave this blank instead of zeros, but that would return NAs in the posterior,
    //       which is not good for model checking and posterior analysis
    for (prev_t in B:(t-1)) {
      theta[prev_t, m] = 0;
      eta[prev_t, m] = 0;
    }
    
    // Evolution for t > first_appearance[m]
    // Note: This models how abilities change over time after first appearance
    for (next_t in (t+1):T) {
      theta[next_t, m] = (1-a) * theta[next_t-1, m] + a * eta[next_t-1, m] + 
                         dot_product(Q_M[m], theta_raw[next_t]) * sqrt(2*a*(1-a));
      eta[next_t, m] = (1-b) * eta[next_t-1, m] + b * theta[next_t-1, m] + 
                       dot_product(Q_M[m], eta_raw[next_t]) * sqrt(2*b*(1-b));
    }
  }
}

model {
  // Priors
  a ~ uniform(0,1);  
  b ~ uniform(0,1); 
  to_vector(theta_raw) ~ normal(0,1);
  to_vector(eta_raw) ~ normal(0,1);
  
  // Likelihood for marathon
  for (n in 1:N_marathon) {
    int t = time_period_marathon[n];
    target += plackett_luce_lpmf(place_marathon[n] | 
                                 competitor_id_marathon[n], 
                                 theta[t]');
  }
  
  // Likelihood for half-marathon
  for (n in 1:N_half_marathon) {
    int t = time_period_half_marathon[n];
    target += plackett_luce_lpmf(place_half_marathon[n] | 
                                 competitor_id_half_marathon[n], 
                                 eta[t]');
  }
}

generated quantities {
  vector[N_marathon + N_half_marathon] log_lik;
  
  // Log-likelihood for marathon
  for (n in 1:N_marathon) {
    int t = time_period_marathon[n];
    log_lik[n] = plackett_luce_lpmf(place_marathon[n] | 
                                    competitor_id_marathon[n], 
                                    theta[t]');
  }
  
  // Log-likelihood for half-marathon
  for (n in 1:N_half_marathon) {
    int t = time_period_half_marathon[n];
    log_lik[N_marathon + n] = plackett_luce_lpmf(place_half_marathon[n] | 
                                                 competitor_id_half_marathon[n], 
                                                 eta[t]');
  }
  array[N_marathon] int y_rep_marathon;
  array[N_half_marathon] int y_rep_half_marathon;
  
  // Generate predicted ranks for marathon
  for (n in 1:N_marathon) {
    int t = time_period_marathon[n];
    int n_comp = n_competitors_marathon[n];
    array[n_comp] int simulated_rank = plackett_luce_rng(theta[t, competitor_id_marathon[n, 1:n_comp]]');
    y_rep_marathon[n] = simulated_rank[1];
  }
  
  // Generate predicted ranks for half-marathon
  for (n in 1:N_half_marathon) {
    int t = time_period_half_marathon[n];
    int n_comp = n_competitors_half_marathon[n];
    array[n_comp] int simulated_rank = plackett_luce_rng(eta[t, competitor_id_half_marathon[n, 1:n_comp]]');
    y_rep_half_marathon[n] = simulated_rank[1];
  }
}


