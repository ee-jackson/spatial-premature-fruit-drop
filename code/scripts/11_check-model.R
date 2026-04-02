#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: check-models
## Desc: perform posterior predictive checks
## Date created: 2023-07-18

# packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("brms")
library("patchwork")
library("stringr")
library("bayesplot")
library("bayestestR")
library("gt")

mod <-
  readRDS(here::here("output", "models", "pheno-repro-adjust",
                     "full_conn_binom_nseeds_abund.rds"))

# Posterior predictive checks ---------------------------------------------

# For binomial data, plots of y and yrep show the proportion of successes
# rather than the raw count

plot_pp_check <- function(model) {
  pp_check(model, ndraws = 500) +
    labs(x = "Proportion of immature seeds", y = "Density") +
    theme_classic(base_size = 20)
}


# MCMC diagnostics --------------------------------------------------------

plot_mcmc_check <- function(model) {
  mcmc_trace(model, regex_pars = "b_",
             iter1 = 1000,
             facet_args = list(ncol = 2)) +
    scale_x_continuous(breaks = seq(1000, 4000, by = 1000)) +
    theme_classic(base_size = 15)
}


# Get posterior param estimates -------------------------------------------

names <-
  c("Intercept",
    "Reproductive conspecific density",
    "Reproductive heterospecific density",
    "Non-reproductive conspecific density",
    "log total seeds",
    "log abundance",
    "Reproductive conspecific density:log total seeds",
    "Reproductive heterospecific density:log total seeds",
    "Non-reproductive conspecific density:log total seeds",
    "Reproductive conspecific density:log species abundance",
    "Reproductive heterospecific density:log species abundance",
    "Non-reproductive conspecific density:log species abundance")

values <-
  c("(Intercept)",
    "b_conn_RH_sc",
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


get_table <- function(model) {
  bayestestR::describe_posterior(model,
                                 effects = "fixed",
                                 component = "all",
                                 ci = 0.95,
                                 ci_method = "HDI",
                                 centrality = "median",
                                 test = FALSE) %>%
    mutate(across(!Rhat & !Parameter, round, 2)) %>%
    mutate(Parameter = str_replace(Parameter, values, names)) %>%
    gt()
}


# Assemble figure ---------------------------------------------------------

rc_t <- get_table(mod)
gtsave(rc_t, here::here("output", "results", "describe_posterior.png"))
rc_t_png <- png::readPNG(here::here("output", "results", "describe_posterior.png"),
                         native = TRUE)

rc_pp <- plot_pp_check(mod)

(rc_pp / rc_t_png) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "pp_check_si.png"),
  width = 500,
  height = 500,
  units = "px"
)

png(
  here::here("output", "figures", "total_con_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)
