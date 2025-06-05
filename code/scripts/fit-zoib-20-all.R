#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: 03_fit-growth-model.R
## Desc: Fit final growth model


# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")


# Get data ----------------------------------------------------------------

consp_nonrepro <-
  readRDS("data/clean/trap_connect_nonrepro_consp_20m.rds")

consp_repro <-
  readRDS("data/clean/trap_connect_repro_consp_20m.rds")

hetero_repro <-
  readRDS("data/clean/trap_connect_repro_hetero_20m.rds")

# combine
# there are some NAs for connectivity - this is when there are no trees
# of that category in the plot in that year, so replace with zero
# only happens for sp4 attb and tri6
trap_connect <-
  consp_nonrepro %>%
  rename(d_conspecific_nonreproductive = connectivity) %>%
  full_join(consp_repro) %>%
  rename(d_conspecific_reproductive = connectivity) %>%
  full_join(hetero_repro) %>%
  rename(d_heterospecific_reproductive = connectivity) %>%
  mutate(across(contains("d_"), ~ replace_na(., 0)))

# don't include traps < 20m from the edge of the plot
# centre and scale connectivity
test_data <-
  trap_connect %>%
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) %>%
  select(- x, - y, - capsules) %>%
  transform(d_conspecific_nonreproductive_sc =
              scale(d_conspecific_nonreproductive)) %>%
  transform(d_conspecific_reproductive_sc =
              scale(d_conspecific_reproductive)) %>%
  transform(d_heterospecific_reproductive_sc =
              scale(d_heterospecific_reproductive)) %>%
  filter(sum_parts >= 3) %>% # do we need to do this?
  mutate_at(c("sp4", "year", "trap", "quadrat"), ~as.factor(.))


# Define model ------------------------------------------------------------

zoib_mod <-
  bf(
    proportion_abscised ~
      d_conspecific_nonreproductive_sc +
      d_conspecific_reproductive_sc +
      d_heterospecific_reproductive_sc +
      (1|quadrat/trap) +
      (1|year) +
      (1 +
         d_conspecific_nonreproductive_sc +
         d_conspecific_reproductive_sc +
         d_heterospecific_reproductive_sc | sp4),
    phi ~
      d_conspecific_nonreproductive_sc +
      d_conspecific_reproductive_sc +
      d_heterospecific_reproductive_sc +
      (1|quadrat/trap) +
      (1|year) +
      (1 +
         d_conspecific_nonreproductive_sc +
         d_conspecific_reproductive_sc +
         d_heterospecific_reproductive_sc | sp4),
    zoi ~
      d_conspecific_nonreproductive_sc +
      d_conspecific_reproductive_sc +
      d_heterospecific_reproductive_sc +
      (1|quadrat/trap) +
      (1|year) +
      (1 +
         d_conspecific_nonreproductive_sc +
         d_conspecific_reproductive_sc +
         d_heterospecific_reproductive_sc | sp4),
    coi ~
      d_conspecific_nonreproductive_sc +
      d_conspecific_reproductive_sc +
      d_heterospecific_reproductive_sc +
      (1|quadrat/trap) +
      (1|year) +
      (1 +
         d_conspecific_nonreproductive_sc +
         d_conspecific_reproductive_sc +
         d_heterospecific_reproductive_sc | sp4),
    family = zero_one_inflated_beta()
  )

zero_one_inflated_beta( link = "logit",
                        link_phi = "log",
                        link_zoi = "logit",
                        link_coi = "logit"
)


# Set priors --------------------------------------------------------------

priors <-
  set_prior(
    prior = "normal(0, 10)",
    class = "b")

# Fit model ---------------------------------------------------------------

fit <-
  brm(
    formula = zoib_mod,
    data = test_data,
    prior = priors,
    chains = 4,
    #control = list(max_treedepth = 12, adapt_delta = 0.99),
    cores = 4,
    seed = 123,
    file = "output/models/all_dens_20m"
  )

add_criterion(
  fit,
  criterion = c("loo", "waic")
)

print(fit$criteria$loo)
print(fit$criteria$waic)
