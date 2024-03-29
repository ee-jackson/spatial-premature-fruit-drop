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

model_rc <- readRDS(
  here::here("output", "models", "202308","zoib_20m.rds"))
model_tc <- readRDS(
  here::here("output", "models", "202308","zoib_total_c_20m.rds"))
model_h <- readRDS(
  here::here("output", "models", "202308","zoib_hetero_20m.rds"))

# Posterior predictive checks ---------------------------------------------

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

get_table <- function(model) {
  bayestestR::describe_posterior(model,
                                 ci = 0.95,
                                 ci_method = "HDI",
                                 centrality = "median",
                                 test = FALSE) %>%
    mutate(across(!Rhat & !Parameter, round, 2)) %>%
    gt()
}

# plot for Reproductive conspecifics model --------------------------------

rc_t <- get_table(model_rc)
gtsave(rc_t, here::here("output", "results", "repro_con_20m.png"))
rc_t_png <- png::readPNG(here::here("output", "results", "repro_con_20m.png"),
                         native = TRUE)

rc_pp <- plot_pp_check(model_rc)
rc_mcmc <- plot_mcmc_check(model_rc)

((rc_pp / rc_t_png) | rc_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "repro_con_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)

# plot for Reproductive conspecifics model --------------------------------

tc_t <- get_table(model_tc)
gtsave(tc_t, here::here("output", "results", "total_con_20m.png"))
tc_t_png <- png::readPNG(here::here("output", "results", "total_con_20m.png"),
                         native = TRUE)

tc_pp <- plot_pp_check(model_tc)
tc_mcmc <- plot_mcmc_check(model_tc)

((tc_pp / tc_t_png) | tc_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "total_con_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)

# plot for Reproductive heterospecifics model --------------------------------

h_t <- get_table(model_h)
gtsave(h_t, here::here("output", "results", "het_20m.png"))
h_t_png <- png::readPNG(here::here("output", "results", "het_20m.png"),
                        native = TRUE)

h_pp <- plot_pp_check(model_h)
h_mcmc <- plot_mcmc_check(model_h)

((h_pp / h_t_png) | h_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "repro_het_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)
