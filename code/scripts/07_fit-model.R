#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit_binomial_models.R
## Desc: fit ZOIB and ZIB on full dataset with species as a random effect
## Date: April 2021

rm(list = ls())

library("tidyverse")
library("rstan")
library("brms")
library("tidybayes")
library("purrr")
library("future")

load("/home/users/ft840275/SpatialPatterns/data/cleanData.RData")

trapConnect$abscised_seeds <- as.integer(trapConnect$abscised_seeds)
trapConnect$total_seeds <- as.integer(trapConnect$total_seeds)
trapConnect$year <- as.factor(trapConnect$year)
trapConnect$trap <- as.factor(trapConnect$trap)

trapConnect <- transform(trapConnect, 
    CI_pred.sc = scale(CI_pred)
    )

testdat <- subset(trapConnect, sum_parts > 10) # it seems like this is key, wouldn't run when I subset to > 5

testdat$OLRE <- seq_len(nrow(testdat))

# use default priors
# prior1 <- c(
#     prior(normal(0, 10), class = Intercept, coef = ""), # the intercept
#     prior(cauchy(0, 10), class = sd), # the variance
#     prior(normal(0, 10), class = b) # the mean ?
#     )

## Zero-inflated binomial (ZIB)
zib_mod <- bf(
  formula = abscised_seeds | trials(total_seeds) ~ 
  CI_pred.sc + (1|trap) + (1|year) + (1 + CI_pred.sc|SP4) + (1|OLRE) + (1|quadrat),
  zi ~ CI_pred.sc + (1|trap) + (1|year)+ (1 + CI_pred.sc|SP4) + (1|OLRE) + (1|quadrat),
  family = zero_inflated_binomial(logit))


## Zero-One inflated beta (ZOIB)
zoib_mod <- bf(
  proportion_abscised ~ CI_pred.sc + (1|trap) + (1|year) + (1 + CI_pred.sc|SP4) + (1|quadrat), 
  phi ~ CI_pred.sc + (1|trap) + (1|year) + (1 + CI_pred.sc|SP4) + (1|quadrat), 
  zoi ~ CI_pred.sc + (1|trap) + (1|year) + (1 + CI_pred.sc|SP4) + (1|quadrat), 
  coi ~ CI_pred.sc + (1|trap) + (1|year) + (1 + CI_pred.sc|SP4) + (1|quadrat), 
  family = zero_one_inflated_beta()
)

zero_one_inflated_beta( link = "logit", 
  link_phi = "log", 
  link_zoi = "logit", 
  link_coi = "logit"
)


# the models will run sequentially, 4 cores per model in parallel

run_my_models <- function(data, name) {
  ZIBfit <- brm(formula = zib_mod,  
              data = data,
              chains = 4,
              warmup = 2000, iter = 5000,
              control = list(max_treedepth = 15),
              cores = 4, seed = 123, 
              file = paste0(name, sep = "_", "ZIB") )

  ZOIBfit <- brm(formula = zoib_mod,  
              data = data,
              chains = 4,
              warmup = 2000, iter = 5000,
              control = list(max_treedepth = 15),
              cores = 4, seed = 123, 
              file = paste0(name, sep = "_", "ZOIB") )
  ZIBfit <- brms::add_criterion(ZIBfit, "loo")
  ZOIBfit <- brms::add_criterion(ZOIBfit, "loo")
  tibble(
    SP4 = name, 
    model = c("ZIB", "ZOIB"),
    fit = list(ZIBfit, ZOIBfit)
  )
}

all_models_quad_rslope <- run_my_models(data = testdat, name = "full_sp_quad_rslope")

#save on cluster
save(all_models_quad_rslope, file = "/home/users/ft840275/SpatialPatterns/output/mods_full_quad_rslope.RData")