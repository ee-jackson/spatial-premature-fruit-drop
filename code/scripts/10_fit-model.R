#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit-model
## Desc: Fit final model

options(mc.cores = 4)
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


# Fit model ---------------------------------------------------------------

# fit <-
#   brm(
#     formula = binom_mod,
#     data = test_data,
#     prior = priors,
#     sample_prior = "yes",
#     chains = 4,
#     iter = 5000,
#     control = list(adapt_delta = 0.95),
#     cores = 4,
#     seed = 123,
#     init_r = 0.1,
#     file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds"
#   )
#
# brms::add_criterion(x = fit, criterion = "loo",
#                     overwrite = FALSE)
#
# print(fit$criteria$loo)
#
# fit_int <-
#   brm(
#     formula = binom_mod_int,
#     data = test_data,
#     prior = priors,
#     sample_prior = "yes",
#     chains = 4,
#     iter = 5000,
#     control = list(adapt_delta = 0.95),
#     cores = 4,
#     seed = 123,
#     init_r = 0.1,
#     file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds_interact"
#   )
#
# brms::add_criterion(x = fit_int, criterion = "loo",
#                     overwrite = FALSE)
#
# print(fit_int$criteria$loo)

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
    file_refit = "always",
    file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds_abund"
  )

brms::add_criterion(x = fit_abund, criterion = "loo")

print(fit_abund$criteria$loo)

# fit_abund <-
#   brm(
#     formula = binom_mod_abund,
#     data = test_data_preds,
#     prior = priors,
#     sample_prior = "yes",
#     chains = 4,
#     iter = 5000,
#     control = list(adapt_delta = 0.95),
#     cores = 4,
#     seed = 123,
#     init_r = 0.1,
#     file_refit = "never",
#     file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds_abund_predscomp"
#   )
#
# brms::add_criterion(x = fit_abund, criterion = "loo")
#
# print(fit_abund$criteria$loo)
#
# fit_abund_preds <-
#   brm(
#     formula = binom_mod_abund_preds,
#     data = test_data_preds,
#     prior = priors,
#     sample_prior = "yes",
#     chains = 4,
#     iter = 5000,
#     control = list(adapt_delta = 0.95),
#     cores = 4,
#     seed = 123,
#     init_r = 0.1,
#     file_refit = "never",
#     file = "output/models/pheno-repro-adjust/full_conn_binom_nseeds_abund_preds"
#   )
#
# brms::add_criterion(x = fit_abund_preds, criterion = "loo")
#
# print(fit_abund_preds$criteria$loo)
#
# loo_compare(fit_abund, fit_abund_preds)
