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

# load data ---------------------------------------------------------------

data <- readRDS(here::here("data", "data_for_model.rds"))
model <- readRDS(here::here("output", "models", "zoib_capsules_quad_rslope_nest.rds"))
sp_data <- read.csv(here::here("data", "species_list.csv"))

model <- readRDS("output/models/full_sp_quad_rslope_ZOIB.rds")
model <- readRDS("output/models/zoib_capsules_quad_rslope_nest.rds")

# format / transform data -------------------------------------------------

data %>%
  mutate(abscised_seeds = as.integer(abscised_seeds),
         total_seeds = as.integer(total_seeds),
         year = as.factor(year),
         trap = as.factor(trap)) %>%
  mutate(connectivity_pred_sc = scale(connectivity_pred)) %>%
  filter(sum_parts >= 3) -> testdat

# plot posterior prediction -----------------------------------------------

# Compute posterior samples of the expected value/mean of the posterior predictive distribution
# giving posterior draws from the expectation of the posterior predictive; i.e. posterior distributions of conditional means

testdat %>%
  tidybayes::add_epred_draws(model, ndraws = 100, re_formula = NA) -> test_fit_na


test_fit_na %>%
  ggplot(aes(x = connectivity_pred_sc)) +
  geom_point(data = testdat, aes(y = proportion_abscised, size = total_seeds),
             colour = "black", alpha = 0.4, pch = 16, size = 0.5) +
  stat_lineribbon(aes(y = .epred), alpha = 0.7, size = 0.25, colour = "darkgreen", show.legend = FALSE) +
  theme_classic(base_size = 7) +
  xlab("Connectivity") +
  ylab("Rate of premature seed abscission") +
  scale_x_continuous(expand = c(0 , 0)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_brewer(palette = "Greens") -> fitted

pdf("epred_draws.pdf", width = 3.15, height = 3.15)
fitted
dev.off()

test_fit_na %>%
  ggplot(aes(x = connectivity_pred_sc)) +
  geom_point(data = testdat, aes(y = proportion_abscised, size = total_seeds),
             colour = "black", alpha = 0.4, pch = 16, size = 0.5) +
  geom_line(aes(y = .epred, group = .draw), alpha = 0.7, colour = "darkgreen", show.legend = FALSE) +
  theme_classic(base_size = 7) +
  xlab("Connectivity") +
  ylab("Rate of premature seed abscission") +
  scale_x_continuous(expand = c(0 , 0)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_brewer(palette = "Greens") -> fitted2

pdf("epred_draws_spag.pdf", width = 3.15, height = 3.15)
fitted2
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


# plot species estimates --------------------------------------------------

# Extract posterior samples of species-level mu estimate, slope not intercept
# this will show the effect of connectivity on the mean of non- zero or one values (mu), per species

model %>%
  brms::posterior_samples() %>%
  select(contains("r_sp4[")) %>%
  select(contains(",connectivity_pred_sc]")) %>%
  pivot_longer(cols = contains("r_sp4["), names_to = "sp4") %>%
  mutate(sp4 = as.character(sp4)) %>%
  mutate(sp4 = gsub("r_sp4[", "", sp4, fixed = TRUE)) %>%
  mutate(sp4 = gsub(",connectivity_pred_sc]", "", sp4, fixed = TRUE)) -> sp_ests

# add species  names
sp_data %>%
  mutate(taxa = paste(genus, species, sep = " ")) -> sp_data

sp_data %>%
  select(sp4, taxa) %>%
  unique() %>%
  right_join(sp_ests) -> sp_ests_n

# plot
pdf("sp_ests.pdf", width = 3.15, height = 3.15)
ggplot(sp_ests_n, aes(y = reorder(taxa, value), x = as.numeric(value))) +
  ggdist::stat_gradientinterval(.width = 0.95, fill = "forestgreen",
                                interval_size = 0.005, stroke = 0.5,
                                shape = 21,
                                point_fill = "white") +
  labs(x = expression(paste("Effect of connectivity on\nmean rate of premature seed abscission, ", Beta)), y = "") +
  geom_vline(xintercept = 0, linetype = 1, size = 0.25, colour = "blue") +
  coord_cartesian(xlim = c(-1, 1), expand = FALSE) +
  theme_classic(base_size = 7 ) +
  theme(axis.text.y = element_text(face = "italic", size = 5, colour = "black"),
        axis.title.x = element_text(hjust = 0, vjust = 1, margin = margin(t = 10)),
        plot.margin = margin(2, 3, 2, 2, "mm"))
dev.off()

# Estimated mu parameters from the ZOIB fit, as points and error bars (95% CIs)
