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

load("data/cleanData.RData")
model <- readRDS("output/models/full_sp_quad_rslope_ZOIB.rds")

# format / transform data -------------------------------------------------

trapConnect$abscised_seeds <- as.integer(trapConnect$abscised_seeds)
trapConnect$total_seeds <- as.integer(trapConnect$total_seeds)
trapConnect$year <- as.factor(trapConnect$year)
trapConnect$trap <- as.factor(trapConnect$trap)

trapConnect <- transform(trapConnect,
                         CI_pred.sc = scale(CI_pred)
)

testdat <- subset(trapConnect, sum_parts > 10)


# plot posterior prediction -----------------------------------------------

# Compute posterior samples of the expected value/mean of the posterior predictive distribution
# giving posterior draws from the expectation of the posterior predictive; i.e. posterior distributions of conditional means

testdat %>%
  add_epred_draws(model, ndraws = 1000, re_formula = NA) -> test_fit_na


test_fit_na %>%
  ggplot(aes(x = CI_pred.sc)) +
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

# plot parameter estimates ------------------------------------------------

# Extract posterior samples of specified parameters

model %>%
  brms::posterior_samples() %>%
  select(b_CI_pred.sc, b_coi_CI_pred.sc, b_zoi_CI_pred.sc, b_phi_CI_pred.sc) %>%
  pivot_longer(cols = contains("CI_pred.sc"), names_to = "parameters") %>%
  mutate(parameters = recode(parameters, b_CI_pred.sc = "mu", b_coi_CI_pred.sc = "coi",
                             b_zoi_CI_pred.sc = "zoi", b_phi_CI_pred.sc = "phi")) %>%
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
  select(contains("r_SP4[")) %>%
  select(contains(",CI_pred.sc]")) %>%
  pivot_longer(cols = contains("r_SP4["), names_to = "SP4") %>%
  mutate(SP4 = as.character(SP4)) %>%
  mutate(SP4 = gsub("r_SP4[", "", SP4, fixed = TRUE)) %>%
  mutate(SP4 = gsub(",CI_pred.sc]", "", SP4, fixed = TRUE)) -> sp_ests

# add species  names
testdat %>%
  mutate(taxa = paste(GENUS, SPECIES, sep = " ")) -> testdat

testdat %>%
  select(SP4, taxa) %>%
  unique() %>%
  right_join(sp_ests) -> sp_ests_n

# plot
pdf("sp_ests.pdf", width = 3.15, height = 3.15)
ggplot(sp_ests_n, aes(y = reorder(taxa, value), x = as.numeric(value))) +
  ggdist::stat_gradientinterval(.width = 0.95, fill = "forestgreen",
                                interval_size =0.005, stroke = 0.5,
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
