#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit-model
## Desc: Fit models
## takes ~11 hours to refit each model on HPC cluster
## with 4 CPUs and 50G memory per CPU

options(mc.cores = 4)

# the fit is always loaded if it exists and fitting is skipped
# change to "always" to refit models
options(brms.file_refit = "on_change")

# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")
library("extraDistr", lib.loc = "~/local/rlibs")


# Get data ----------------------------------------------------------------

test_data <-
  readRDS("data/clean/trap_connect_buffer.rds")

test_data_preds <-
  readRDS("data/clean/trap_connect_buffer_preds.rds")


# Define models -----------------------------------------------------------

# main model formula
binom_mod_abund <-
  bf(
    abscised_seeds | trials(total_seeds) ~
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_total_seeds_sc +
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_median_abundance_sc +
      (1 | quadrat/trap) +
      (1 | year) +
      (1 + conn_RC_sc + conn_RH_sc + conn_NRC_sc | sp4),
    family = beta_binomial(link = "logit", link_phi = "log")
  )

# formula including seed predator presence term
binom_mod_abund_preds <-
  bf(
    abscised_seeds | trials(total_seeds) ~
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_total_seeds_sc +
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_median_abundance_sc +
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * SeedPred_pres +
      (1 | quadrat/trap) +
      (1 | year) +
      (1 + conn_RC_sc + conn_RH_sc + conn_NRC_sc | sp4),
    family = beta_binomial(link = "logit", link_phi = "log")
  )


# Set priors --------------------------------------------------------------

priors <- c(
  set_prior("normal(0, 1)", class = "b")
)


# Fit models --------------------------------------------------------------

# Fit main model
fit_abund <-
  brm(
    formula = binom_mod_abund,
    data = test_data,
    prior = priors,
    sample_prior = "yes",
    chains = 4,
    iter = 5000,
    control = list(adapt_delta = 0.95),
    cores = 4,
    seed = 123,
    init_r = 0.1,
    file = "output/models/full_conn_binom_nseeds_abund"
  )

brms::add_criterion(x = fit_abund, criterion = "loo")

print(fit_abund$criteria$loo)

# Fit model excluding dominant species, Alseis blackiana
fit_abund_no_alsb <-
  brm(
    formula = binom_mod_abund,
    data = filter(test_data, sp4 != "alsb"),
    prior = priors,
    sample_prior = "yes",
    chains = 4,
    iter = 5000,
    control = list(adapt_delta = 0.95),
    cores = 4,
    seed = 123,
    init_r = 0.1,
    file = "output/models/full_conn_binom_nseeds_abund_no_alsb"
  )

brms::add_criterion(x = fit_abund_no_alsb, criterion = "loo")

print(fit_abund_no_alsb$criteria$loo)

# Fit model with seed predator presence term
fit_abund_preds <-
  brm(
    formula = binom_mod_abund_preds,
    data = test_data_preds,
    prior = priors,
    sample_prior = "yes",
    chains = 4,
    iter = 5000,
    control = list(adapt_delta = 0.95),
    cores = 4,
    seed = 123,
    init_r = 0.1,
    file = "output/models/full_conn_binom_nseeds_abund_preds"
  )

brms::add_criterion(x = fit_abund_preds, criterion = "loo")

print(fit_abund_preds$criteria$loo)

# Fit model to compare with seed predator presence term
# only a subset of species have seed predator data
# so comparison is restricted to the shared set of observed data
fit_abund_predscomp <-
  brm(
    formula = binom_mod_abund,
    data = test_data_preds,
    prior = priors,
    sample_prior = "yes",
    chains = 4,
    iter = 5000,
    control = list(adapt_delta = 0.95),
    cores = 4,
    seed = 123,
    init_r = 0.1,
    file = "output/models/full_conn_binom_nseeds_abund_predscomp"
  )

brms::add_criterion(x = fit_abund_predscomp, criterion = "loo")

print(fit_abund$criteria$loo)

loo_compare(fit_abund, fit_abund_predscomp)
