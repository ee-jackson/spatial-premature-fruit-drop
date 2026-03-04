#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit-zoib
## Desc: Fit final model

options(mc.cores = 4)

# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")


# Get data ----------------------------------------------------------------

test_data <- readRDS("data/clean/trap_connect.rds")


# Define model ------------------------------------------------------------

binom_mod <-
  bf(
    abscised_seeds | trials(total_seeds) ~
      conn_RC_sc + conn_RH_sc + conn_NRC_sc + log_total_seeds_sc +
      (1|quadrat/trap) +
      (1|year) +
      (1 + conn_RC_sc + conn_RH_sc + conn_NRC_sc | sp4),
    family = beta_binomial(link = "logit", link_phi = "log")
  )

binom_mod_int <-
  bf(
    abscised_seeds | trials(total_seeds) ~
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_total_seeds_sc +
      (1|quadrat/trap) +
      (1|year) +
      (1 + conn_RC_sc + conn_RH_sc + conn_NRC_sc | sp4),
    family = beta_binomial(link = "logit", link_phi = "log")
  )

# Set priors --------------------------------------------------------------

priors <- c(
  set_prior("normal(0, 1)", class = "b")
)


# Fit model ---------------------------------------------------------------

fit <-
  brm(
    formula = binom_mod,
    data = test_data,
    prior = priors,
    sample_prior = "yes",
    chains = 4,
    iter = 5000,
    control = list(adapt_delta = 0.95),
    cores = 4,
    seed = 123,
    file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds",
    save_pars = save_pars(all = TRUE),
    file_refit = "on_change"
  )

# brms::add_criterion(x = fit, criterion = "loo",
#                     newdata = test_data)
#
# print(fit$criteria$loo)


fit_int <-
  brm(
    formula = binom_mod_int,
    data = test_data,
    prior = priors,
    sample_prior = "yes",
    chains = 4,
    iter = 5000,
    control = list(adapt_delta = 0.95),
    cores = 4,
    seed = 123,
    file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds_interact",
    save_pars = save_pars(all = TRUE),
    file_refit = "always"
  )

# brms::add_criterion(x = fit_int, criterion = "loo",
#                     newdata = test_data)
#
# print(fit_int$criteria$loo)

