
library("dplyr")
library("rstan")
library("brms")

readRDS("/home/users/ft840275/spatial_patterns/data/clean/total_trap_connect_20m.rds") -> trap_connect

edge_traps <- subset(trap_connect, x > 980 | x < 20 | y > 480 | y < 20)

edge_traps %>%
  select(trap) %>%
  distinct() -> edge_traps_list

trap_connect %>%
  filter(!trap %in% edge_traps_list$trap) %>%
  select(- x, - y, - capsules) %>%
  transform(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat_20m

testdat_20m$year <- as.factor(testdat_20m$year)
testdat_20m$trap <- as.factor(testdat_20m$trap)
testdat_20m$quadrat <- as.factor(testdat_20m$quadrat)

zoib_mod <- bf(
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

brm(
  formula = zoib_mod,
  data = testdat_20m,
  chains = 4,
  warmup = 2000,
  iter = 5000,
  control = list(max_treedepth = 12, adapt_delta = 0.99),
  cores = 4,
  seed = 123,
  file = "zoib_total_c_20m"
) -> fit

add_criterion(
  fit,
  criterion = c("loo", "waic")
)

print(fit$criteria$loo)
print(fit$criteria$waic)
