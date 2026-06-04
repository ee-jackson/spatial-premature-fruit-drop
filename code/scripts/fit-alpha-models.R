#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit-zoib
## Desc: Fit final model

options(mc.cores = 4)
options(brms.file_refit = "always")


# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")
library("extraDistr", lib.loc = "~/local/rlibs")


# Get data ----------------------------------------------------------------

test_data <-
  readRDS("data/clean/connect_all_alpha_buffer.rds") %>%
  mutate(log_total_seeds = log(total_seeds))

alphas <- unique(test_data$alpha)

common_traps <-
  test_data %>%
  filter(alpha == alphas[[1]]) %>%
  filter(location == "interior") %>%
  select(trap) %>%
  distinct() %>%
  pull(trap) %>%
  droplevels()

test_data <-
  test_data %>%
  filter(trap %in% common_traps)


# Define model ------------------------------------------------------------

binom_mod <-
  bf(
    abscised_seeds | trials(total_seeds) ~
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_total_seeds_sc +
      (conn_RC_sc + conn_RH_sc + conn_NRC_sc) * log_median_abundance_sc +
      (1 | quadrat/trap) +
      (1 | year) +
      (1 + conn_RC_sc + conn_RH_sc + conn_NRC_sc | sp4),
    family = beta_binomial(link = "logit", link_phi = "log")
  )

# Set priors --------------------------------------------------------------

priors <- c(
  set_prior("normal(0, 1)", class = "b")
)


# alpha = 0.07 ------------------------------------------------------------

test_data07 <-
  filter(test_data, alpha == alphas[[1]]) %>%
  filter(location == "interior") %>%
  mutate(
    conn_RC_sc = as.numeric(scale(conn_RC)),
    conn_RH_sc = as.numeric(scale(conn_RH)),
    conn_NRC_sc = as.numeric(scale(conn_NRC)),
    log_total_seeds_sc = as.numeric(scale(log_total_seeds))
  )

fit07 <-
  brm(
    formula = binom_mod,
    data = test_data07,
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123,
    init_r = 0.1,
    control = list(adapt_delta = 0.95),
    file = "output/models/alpha-tests-subset-abund/bbinom_0.07"
  )

# interior20 <-
#   test_data20 %>%
#   filter(trap %in% common_traps)

brms::add_criterion(x = fit07, criterion = "loo",
                    # newdata = interior20,
                    overwrite = TRUE)

rm(fit07, test_data07)


# alpha = 0.099 -----------------------------------------------------------

test_data09 <-
  filter(test_data, alpha == alphas[[2]]) %>%
  filter(location == "interior") %>%
  mutate(
    conn_RC_sc = as.numeric(scale(conn_RC)),
    conn_RH_sc = as.numeric(scale(conn_RH)),
    conn_NRC_sc = as.numeric(scale(conn_NRC)),
    log_total_seeds_sc = as.numeric(scale(log_total_seeds))
  )

fit09 <-
  brm(
    formula = binom_mod,
    data = test_data09,
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123,
    init_r = 0.1,
    control = list(adapt_delta = 0.95),
    file = "output/models/alpha-tests-subset-abund/bbinom_0.099"
  )

# interior15 <-
#   test_data15 %>%
#   filter(trap %in% common_traps)

brms::add_criterion(x = fit09, criterion = "loo",
                    # newdata = interior15,
                    overwrite = TRUE)

rm(fit09, test_data09)


# alpha = 0.14 ------------------------------------------------------------

test_data14 <-
  filter(test_data, alpha == alphas[[3]]) %>%
  filter(location == "interior") %>%
  mutate(
    conn_RC_sc = as.numeric(scale(conn_RC)),
    conn_RH_sc = as.numeric(scale(conn_RH)),
    conn_NRC_sc = as.numeric(scale(conn_NRC)),
    log_total_seeds_sc = as.numeric(scale(log_total_seeds))
  )

fit14 <-
  brm(
    formula = binom_mod,
    data = test_data14,
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123,
    init_r = 0.1,
    control = list(adapt_delta = 0.95),
    file = "output/models/alpha-tests-subset-abund/bbinom_0.14"
  )

# interior10 <-
#   test_data10 %>%
#   filter(trap %in% common_traps)

brms::add_criterion(x = fit14, criterion = "loo",
                    # newdata = interior10,
                    overwrite = TRUE)

rm(fit14, test_data14)


# alpha = 0.199 -----------------------------------------------------------

test_data199 <-
  filter(test_data, alpha == alphas[[4]]) %>%
  filter(location == "interior") %>%
  mutate(
    conn_RC_sc = as.numeric(scale(conn_RC)),
    conn_RH_sc = as.numeric(scale(conn_RH)),
    conn_NRC_sc = as.numeric(scale(conn_NRC)),
    log_total_seeds_sc = as.numeric(scale(log_total_seeds))
  )

fit199 <-
  brm(
    formula = binom_mod,
    data = test_data199,
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123,
    init_r = 0.1,
    control = list(adapt_delta = 0.95),
    file = "output/models/alpha-tests-subset-abund/bbinom_0.199"
  )

# interior05 <-
#   test_data05 %>%
#   filter(trap %in% common_traps)

brms::add_criterion(x = fit199, criterion = "loo",
                    # newdata = interior05,
                    overwrite = TRUE)

rm(fit199, test_data199)


# alpha = 0.28 ------------------------------------------------------------

test_data28 <-
  filter(test_data, alpha == alphas[[5]]) %>%
  filter(location == "interior") %>%
  mutate(
    conn_RC_sc = as.numeric(scale(conn_RC)),
    conn_RH_sc = as.numeric(scale(conn_RH)),
    conn_NRC_sc = as.numeric(scale(conn_NRC)),
    log_total_seeds_sc = as.numeric(scale(log_total_seeds))
  )

fit28 <-
  brm(
    formula = binom_mod,
    data = test_data28,
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123,
    init_r = 0.1,
    control = list(adapt_delta = 0.95),
    file = "output/models/alpha-tests-subset-abund/bbinom_0.28"
  )

# interior05 <-
#   test_data05 %>%
#   filter(trap %in% common_traps)

brms::add_criterion(x = fit28, criterion = "loo",
                    # newdata = interior05,
                    overwrite = TRUE)

rm(fit28, test_data28)


# alpha = 0.40 ------------------------------------------------------------

test_data40 <-
  filter(test_data, alpha == alphas[[6]]) %>%
  filter(location == "interior") %>%
  mutate(
    conn_RC_sc = as.numeric(scale(conn_RC)),
    conn_RH_sc = as.numeric(scale(conn_RH)),
    conn_NRC_sc = as.numeric(scale(conn_NRC)),
    log_total_seeds_sc = as.numeric(scale(log_total_seeds))
  )

fit40 <-
  brm(
    formula = binom_mod,
    data = test_data40,
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123,
    init_r = 0.1,
    control = list(adapt_delta = 0.95),
    file = "output/models/alpha-tests-subset-abund/bbinom_0.40"
  )

# interior05 <-
#   test_data05 %>%
#   filter(trap %in% common_traps)

brms::add_criterion(x = fit40, criterion = "loo",
                    # newdata = interior05,
                    overwrite = TRUE)

rm(fit40, test_data40)

