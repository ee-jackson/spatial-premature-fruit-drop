#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: make-figures
## Desc:
## Date created: 2021-07-27

# packages ----------------------------------------------------------------

library("tidyverse")
library("brms")
library("ggdist")
library("tidybayes")
library("patchwork")

# load data ---------------------------------------------------------------

data <- readRDS(here::here("data", "clean", "data_for_model_20m.rds"))
model <- readRDS(here::here("output", "models", "202307","zoib_capsules_20m.rds"))
sp_data <- read.csv(here::here("data", "clean", "species_list.csv"))


# format / transform data -------------------------------------------------

data %>%
  mutate(abscised_seeds = as.integer(abscised_seeds),
         total_seeds = as.integer(total_seeds),
         year = as.factor(year),
         trap = as.factor(trap)) %>%
  mutate(connectivity_sc = scale(connectivity_pred)) %>%
  filter(sum_parts >= 3) -> testdat

# plot posterior prediction -----------------------------------------------

# Compute posterior samples of the expected value/mean of the posterior predictive distribution
# giving posterior draws from the expectation of the posterior predictive; i.e. posterior distributions of conditional means

conditional_effects(model, re_formula = NA)-> cond_ef

testdat %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 101)
  ) %>%
  add_epred_draws(model, ndraws = 1000, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat$connectivity_sc, 'scaled:scale') +
      attr(testdat$connectivity_sc, 'scaled:center')
  ) %>%
  ggplot(aes(
    x = connectivity_us,
    y = .epred,
    fill_ramp = after_stat(.width)
  )) +
  stat_lineribbon(.width = ppoints(50), fill = "forestgreen") +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "forestgreen")) +
  geom_point(
    data = testdat,
    aes(
      x = connectivity_pred,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.7, size = 2.5,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 30) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003)) +
  xlab("Connectivity") +
  ylab("Rate of premature seed abscission") +
  theme(legend.position = "none") -> fig

png(
  here::here("output", "figures", "overall-effect-20m.png"),
  width = 945,
  height = 945,
  units = "px",
  type = "cairo"
)
fig
dev.off()

# plot parameter estimates ------------------------------------------------

# Extract posterior samples of specified parameters

model %>%
  brms::posterior_samples() %>%
  select(b_connectivity_pred_sc, b_coi_connectivity_pred_sc, b_zoi_connectivity_pred_sc, b_phi_connectivity_pred_sc) %>%
  pivot_longer(cols = contains("connectivity_pred_sc"), names_to = "parameters") %>%
  mutate(parameters = recode(parameters, b_connectivity_pred_sc = "mu", b_coi_connectivity_pred_sc = "coi",
                             b_zoi_connectivity_pred_sc = "zoi", b_phi_connectivity_pred_sc = "phi")) %>%
  mutate(parameters = as.character(parameters)) -> posterior_slope_params

pdf("post_params_epred.pdf", width = 3.15, height = 3.15)
ggplot(posterior_slope_params, aes(y = reorder(parameters, abs(value)), x = as.numeric(value))) +
  ggdist::stat_halfeye(normalize = "xy", .width = 0.95, fill = "forestgreen",
                       interval_size =0.005, stroke = 0.5,
                       shape = 21, slab_alpha = 0.6, point_size = 1,
                       point_fill = "white") +
  labs(x = "Effect of connectivity on posterior means", y = "") +
  geom_vline(xintercept = 0, linetype = 1, size = 0.25, colour = "blue") +
  scale_y_discrete(labels = c("mu" = substitute(Beta),
                              "phi"  = substitute(phi),
                              "coi" = substitute(gamma),
                              "zoi" = substitute(alpha))) +
  theme_classic(base_size = 7 )+
  theme(axis.text.y = element_text(face = "italic", size = 12, colour = "black"))
dev.off()


