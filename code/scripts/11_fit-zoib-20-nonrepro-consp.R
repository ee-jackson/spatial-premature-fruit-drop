#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: 03_fit-growth-model.R
## Desc: Fit final growth model


# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")


# Get data ----------------------------------------------------------------

trap_connect <-
  readRDS("data/clean/trap_connect_nonrepro_consp_20m_dioecious.rds")

not_wind_disp_species <-
  read.csv("data/clean/species_list.csv") #%>%
  #filter(dsp_wind != TRUE)

# don't include traps < 20m from the edge of the plot
# centre and scale connectivity
test_data <-
  trap_connect %>%
  filter(sp4 %in% not_wind_disp_species$sp4) %>% # drops 30 species
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) %>%
  select(- x, - y, - capsules) %>%
  transform(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) %>%
  mutate_at(c("sp4", "year", "trap", "quadrat"), ~as.factor(.))


# Define model ------------------------------------------------------------

zoib_mod <-
  bf(
    proportion_abscised ~ connectivity_sc + (1|quadrat/trap) + (1|year) + (1 + connectivity_sc|sp4),
    phi ~ connectivity_sc + (1|quadrat/trap) + (1|year) + (1 + connectivity_sc|sp4),
    zoi ~ connectivity_sc + (1|quadrat/trap) + (1|year) + (1 + connectivity_sc|sp4),
    coi ~ connectivity_sc + (1|quadrat/trap) + (1|year) + (1 + connectivity_sc|sp4),
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
    iter = 10000,
    #control = list(max_treedepth = 12, adapt_delta = 0.99),
    cores = 4,
    seed = 123,
    file = "output/models/nonrepro_consp_20m__yesdioecious"
  )

add_criterion(
  fit,
  criterion = c("loo", "waic")
)

print(fit$criteria$loo)
print(fit$criteria$waic)
