#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit_binomial_models.R
## Desc: fit 8 models, ZOIB and ZIB for 4 species
## Date: March 2021

rm(list = ls())

library("tidyverse")
library("rstan")
library("brms")
library("tidybayes")

load("/home/users/ft840275/SpatialPatterns/data/trapConnect.RData")

trapConnect$abscised_seeds <- as.integer(trapConnect$abscised_seeds)
trapConnect$total_seeds <- as.integer(trapConnect$total_seeds)
trapConnect$year <- as.factor(trapConnect$year)
trapConnect$trap <- as.factor(trapConnect$trap)

trapConnect <- transform(trapConnect, 
    CI_pred.sc = scale(CI_pred)
    )

testdat <- subset(trapConnect, sum_parts > 10) # it seems like this is key, wouldn't run when i subset to > 5

testdat$OLRE <- seq_len(nrow(testdat))

# priors
prior1 <- c(
    prior(normal(0, 10), class = Intercept, coef = ""), # the intercept
    prior(cauchy(0, 10), class = sd), # the variance
    prior(normal(0, 10), class = b) # the mean ?
    )

## Zero-inflated binomial
zib_mod <- bf(
  formula = abscised_seeds | trials(total_seeds) ~ 
  CI_pred.sc + (1|trap) + (1|year) + (1|OLRE),
  zi ~ CI_pred.sc + (1|trap) + (1|year) + (1|OLRE),
  family = zero_inflated_binomial(logit))


## Zero-One inflated beta
zoib_mod <- bf(
  proportion_abscised ~ CI_pred.sc + (1|trap) + (1|year), 
  phi ~ CI_pred.sc + (1|trap) + (1|year), 
  zoi ~ CI_pred.sc + (1|trap) + (1|year), 
  coi ~ CI_pred.sc + (1|trap) + (1|year), 
  family = zero_one_inflated_beta()
)

zero_one_inflated_beta( link = "logit", 
  link_phi = "log", 
  link_zoi = "logit", 
  link_coi = "logit"
)

FARO_d <- subset(testdat, SP4 == "FARO")
JACC_d <- subset(testdat, SP4 == "JACC")
HYBP_d <- subset(testdat, SP4 == "HYBP")
BEIP_d <- subset(testdat, SP4 == "BEIP")

# I have 4 datasets, 2 models to fit to each dataset and 4 chains per model
# I will the models (within a dataset) sequentially, and 4 cores per model in parallel.
# tried to parallelise the run_my_models function with futures but couldn't get it to work

run_my_models <- function(data, name) {
  ZIBfit <- brm(formula = zib_mod,  
              data = data,
              chains = 4,
              warmup = 2000, iter = 5000,
              prior = prior1,
              cores = 4)

  ZOIBfit <- brm(formula = zoib_mod,  
              data = data,
              chains = 4,
              warmup = 2000, iter = 5000,
              prior = prior1,
              cores = 4)

  tibble(
    SP4 = name, 
    model = c("ZIB", "ZOIB"),
    fit = list(ZIBfit, ZOIBfit)
  )
}

fits1 <- run_my_models(data=FARO_d, name="FARO") 
fits2 <- run_my_models(data=JACC_d, name="JACC") 
fits3 <- run_my_models(data=HYBP_d, name="HYBP")
fits4 <- run_my_models(data=BEIP_d, name="BEIP") 

all_models <- bind_rows(fits1, fits2, fits3, fits4)

save(all_models, file = "/home/users/ft840275/SpatialPatterns/output/mods_4sp.RData")