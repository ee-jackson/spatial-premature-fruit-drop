#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: ZOIB.R
## Desc: first attempts at using a Zero-One Inflated Beta model
## Date: October 2020

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size=20))
library("brms")

load("trapConnect.RData")


# species that show heavy zero inflation:
# jac1co, tri2tu, alsebl, luehse, hybapr

# subset to one of these
dat.jac <- subset(trapConnect, sp == "jac1co")

# for some reason brm needs to be explicity told that these are integers when doing a binomial model
dat.jac$abscised_seeds <- as.integer(dat.jac$abscised_seeds)
dat.jac$total_seeds <- as.integer(dat.jac$total_seeds)

# simple intercept test to check rstan is working
intercept.test <- brm(formula = abscised_seeds | trials(total_seeds) ~ 1 + (1|trap) + (1|year),  
              data = dat.jac,
              warmup = 1000, iter = 3000, 
              cores = 4, chains = 2, 
              seed = 123, family = binomial(logit))  

# "abscised_seeds | trials(total_seeds)" is the equivelent of cbind(abscised_seeds, viable_seeds) in lme4. 

# this would be my 'ideal' model, though perhaps 3 random effects is too complex?
zoib_model <- bf(proportion_abscised ~ CI + (1|trap) + (1|year) + (1|quadratID),
  phi ~ CI + (1|trap) + (1|year) + (1|quadratID),
  zoi ~ CI + (1|trap) + (1|year) + (1|quadratID),
  coi ~ CI + (1|trap) + (1|year) + (1|quadratID), 
  family = zero_one_inflated_beta()
)

# this won't run
model1 <- brm(formula = zoib_model,  
              data = dat.jac,
              warmup = 1000, iter = 3000, 
              cores = 4, chains = 4,
              seed = 123, save_model = "stan_model") 
