#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: check-models
## Desc: perform posterior predictive checks and get outputs for SI
## Date created: 2023-07-18

# Packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("brms")
library("patchwork")
library("stringr")
library("bayesplot")
library("bayestestR")

mod <-
  readRDS(here::here("output", "models",
                     "full_conn_binom_nseeds_abund.rds"))


# Calculate contrast ------------------------------------------------------

# Choose the context at which to evaluate the contrast
# Because predictors are scaled, 0 = mean value
newdat <- tibble(
  conn_RC_sc = c(1, 0),
  conn_RH_sc = c(0, 1),
  conn_NRC_sc = 0,
  total_seeds = 1,
  log_total_seeds_sc = 0,
  log_median_abundance_sc = 0,
  quadrat = NA,
  trap = NA,
  year = NA,
  sp4 = NA
)

# Posterior draws of the linear predictor
lp <- posterior_linpred(
  mod,
  newdata = newdat,
  re_formula = NA,
  transform = FALSE
)

# Contrast: RC minus RH
contrast_draws <- lp[, 1] - lp[, 2]

quantile(contrast_draws, probs = c(0.025, 0.5, 0.975))
median(contrast_draws)


# Posterior predictive checks ---------------------------------------------

# Check zero inflation
pp_check(mod, type = "stat", stat = function(y) mean(y == 0), ndraws = 500)

# observed zeros
obs_zero <- as.integer(mod$data$abscised_seeds == 0)

# compute predicted overall zero proportion per draw
pp_mod <- posterior_predict(mod, ndraws = 1000)
pred_zero_rate <- apply(pp_mod, 1, function(draw) mean(draw == 0))

# summary
obs_zero_rate <- mean(obs_zero)
median_pred <- median(pred_zero_rate)
ci_pred <- quantile(pred_zero_rate, c(0.025, 0.975))

tibble(
  observed = obs_zero_rate,
  pred_median = median_pred,
  pred_lo = ci_pred[1],
  pred_hi = ci_pred[2]
)


# MCMC diagnostics --------------------------------------------------------

plot_mcmc_check <- function(model) {
  mcmc_trace(model, regex_pars = "b_",
             iter1 = 1000,
             facet_args = list(ncol = 2)) +
    scale_x_continuous(breaks = seq(1000, 4000, by = 1000)) +
    theme_classic(base_size = 15)
}

plot_mcmc_check(mod)


# Get posterior param estimates -------------------------------------------

names <-
  c("Intercept",
    "Reproductive conspecific density",
    "Reproductive heterospecific density",
    "Non-reproductive conspecific density",
    "log Total conspecific seeds",
    "log Species abundance",
    "Reproductive conspecific density:log Total conspecific seeds",
    "Reproductive heterospecific density:log Total conspecific seeds",
    "Non-reproductive conspecific density:log Total conspecific seeds",
    "Reproductive conspecific density:log Species abundance",
    "Reproductive heterospecific density:log Species abundance",
    "Non-reproductive conspecific density:log Species abundance")

values <-
  c("b_Intercept",
    "b_conn_RC_sc",
    "b_conn_RH_sc",
    "b_conn_NRC_sc",
    "b_log_total_seeds_sc",
    "b_log_median_abundance_sc",
    "b_conn_RC_sc:log_total_seeds_sc",
    "b_conn_RH_sc:log_total_seeds_sc",
    "b_conn_NRC_sc:log_total_seeds_sc",
    "b_conn_RC_sc:log_median_abundance_sc",
    "b_conn_RH_sc:log_median_abundance_sc",
    "b_conn_NRC_sc:log_median_abundance_sc")

fixed_eff_out <-
  mod %>%
  bayestestR::describe_posterior(effects = "fixed",
                                 component = "all",
                                 ci = 0.95,
                                 ci_method = "HDI",
                                 centrality = "median",
                                 test = FALSE,
                                 diagnostic = c("ESS", "ESS_bulk", "Rhat")) %>%
    mutate(Parameter = str_replace(Parameter, values, names))

fixed_eff_out %>%
  mutate(across(!c(Rhat, Parameter), ~ formatC(.x, format = "f", digits = 2))) %>%
  mutate(Rhat = formatC(Rhat, format = "f", digits = 3)) %>%
  write_csv(here::here("output", "results", "describe_posterior_fixed.csv"))

lookup <- setNames(names, values)
lookup_r <- setNames(names, gsub("^b_", "", values))

rand_eff_out <-
  mod %>%
  bayestestR::describe_posterior(effects = "random",
                                 component = "all",
                                 ci = 0.95,
                                 ci_method = "HDI",
                                 centrality = "median",
                                 test = FALSE,
                                 diagnostic = c("ESS", "ESS_bulk", "Rhat")) %>%
  mutate(Parameter = str_replace_all(Parameter, lookup_r))

rand_eff_out %>%
  mutate(across(!c(Rhat, Parameter), ~ formatC(.x, format = "f", digits = 2))) %>%
  mutate(Rhat = formatC(Rhat, format = "f", digits = 3)) %>%
  write_csv(here::here("output", "results", "describe_posterior_random.csv"))
