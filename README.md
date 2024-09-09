**Coupled Dynamic Plackett-Luce Model for Marathon and Half Marathon Performance**

This repository contains the code and data for the implementation of a Coupled Dynamic Plackett-Luce Model, as described in the dissertation "A Coupled Dynamic Plackett-Luce Model for Predicting Olympic Marathon Winners" by Kristyna Coufalova.

**Overview**
This project presents a novel approach to modeling and predicting performance in marathon and half marathon events. It extends the traditional Plackett-Luce model to incorporate time-varying abilities of athletes, coupling between marathon and half marathon performances, and age effects on performance.

**Repository Structure**
*explanatory_data_analysis/*

analysis.R: This script performs initial data exploration and visualization. It includes:
            Loading and cleaning of marathon and half marathon data
            Creation of network graphs to visualize competition structures
            Calculation of centrality measures for athletes
            Generation of distribution plots for race participation



*stan_analysis/*
artificial_data/

artificial_analysis_ab.R: This script generates and analyzes simulated data for the model with coupling parameters a and b. It includes:
                          Simulation of athlete abilities over time
                          Fitting of the model to simulated data
                          Comparison of estimated parameters with true values
                          Visualization of results from simulated data



*fitting_stan_models/*

base_models.R: This script implements and fits the base coupled dynamic Plackett-Luce model. It includes:
              Data preparation for Stan
              Fitting of the model with coupling parameters a and b
              Fitting of the model with coupling parameters a, b, c, and d


covariates_models.R: This script extends the base model to include covariates, particularly age effects. It includes:
                    Implementation of models with linear and quadratic age effects
                    Fitting of these extended models to the data



*posterior_analysis/*

posterior_analysis_base_model_ab.R: This script analyzes the posterior distributions from the base model with parameters a and b. It includes:
                                    Trace plots and convergence diagnostics
                                    Analysis of coupling parameters a and b
                                    Visualization of athlete ability evolution over time
                                    Model validation checks

posterior_analysis_base_model_abcd.R: Similar to the above, but for the model with parameters a, b, c, and d.

posterior_analysis_covariates_ab.R: This script analyzes the posterior distributions from the model including age covariates. It includes:
                                    Analysis of age effects on performance
                                    Comparison of models with and without quadratic age terms
                                    Model selection using LOO-CV and WAIC



stan_files/

covariates_ab.stan: Stan model specification for the coupled dynamic Plackett-Luce model with age covariates.
artificial.stan: Stan model specification for simulation checks (with sigma fixed to 1)
non_sq_covariates_ab.stan: Stan model specification for the model with linear age effects only.
stan_ab.stan: Stan model specification for the base model with coupling parameters a and b.
stan_abcd.stan: Stan model specification for the extended model with coupling parameters a, b, c, and d.

Usage

Clone the repository:
Copygit clone https://github.com/KristynaCoufalova/coupled_dynamic_plackett_luce_model.git

Install required R packages:
RCopyinstall.packages(c("rstan", "tidyverse", "bayesplot", "loo", "igraph", "ggraph"))

Run the scripts in the following order:

explanatory_data_analysis/analysis.R
stan_analysis/fitting_stan_models/base_models.R
stan_analysis/fitting_stan_models/covariates_models.R
stan_analysis/posterior_analysis/posterior_analysis_base_model_ab.R
stan_analysis/posterior_analysis/posterior_analysis_base_model_abcd.R
stan_analysis/posterior_analysis/posterior_analysis_covariates_ab.R
stan_analysis/artificial_data/artificial_analysis_ab.R

