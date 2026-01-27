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

data <- readRDS(here::here("data", "clean", "trap_connect_repro_consp_20m_dioecious.rds"))
model <- readRDS(here::here("output", "models", "repro_consp_20m_yesdioecious.rds"))

data_t <- readRDS(here::here("data", "clean", "trap_connect_nonrepro_consp_20m_dioecious.rds"))
model_t <- readRDS(here::here("output", "models", "nonrepro_consp_20m_yesdioecious.rds"))

data_h <- readRDS(here::here("data", "clean", "trap_connect_repro_hetero_20m_dioecious.rds"))
model_h <- readRDS(here::here("output", "models", "repro_hetero_20m_yesdioecious.rds"))

sp_data <- read.csv(here::here("data", "clean", "species_list.csv"))


# format / transform data -------------------------------------------------

data %>%
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) %>%
  select(- x, - y, - capsules) %>%
  filter(sum_parts >= 3) %>%
  mutate(connectivity_sc = scale(connectivity)) %>%
  mutate(
         year = as.factor(year),
         trap = as.factor(trap)) -> testdat

##

data_t %>%
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) %>%
  select(- x, - y, - capsules) %>%
  filter(sum_parts >= 3) %>%
  mutate(connectivity_sc = scale(connectivity)) %>%
  mutate(
    year = as.factor(year),
    trap = as.factor(trap)) -> testdat_t

##

data_h %>%
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) %>%
  select(- x, - y, - capsules) %>%
  filter(sum_parts >= 3) %>%
  mutate(connectivity_sc = scale(connectivity)) %>%
  mutate(
    year = as.factor(year),
    trap = as.factor(trap)) -> testdat_h


# -------------------------------------------------------------------------
# get predictions

testdat %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 50)
  ) %>%
  add_epred_draws(model, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat$connectivity_sc, 'scaled:scale') +
      attr(testdat$connectivity_sc, 'scaled:center')
  ) -> repro_mod

testdat_t %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 50)
  ) %>%
  add_epred_draws(model_t, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat_t$connectivity_sc, 'scaled:scale') +
      attr(testdat_t$connectivity_sc, 'scaled:center')
  ) -> total_mod

testdat_h %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 50)
  ) %>%
  add_epred_draws(model_h, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat_h$connectivity_sc, 'scaled:scale') +
      attr(testdat_h$connectivity_sc, 'scaled:center')
  ) -> hetro_mod

# -------------------------------------------------------------------------
# big 4 panel plot

repro_mod %>%
  ggplot(aes(
    x = connectivity_us,
    y = .epred,
    fill_ramp = after_stat(.width)
  )) +
  ggdist::stat_lineribbon(.width = ppoints(50), fill = "#E69F00", linewidth = 0.25) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "#E69F00")) +
  geom_point(
    data = testdat,
    aes(
      x = connectivity,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.5, size = 0.01,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0,1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") -> p1


total_mod %>%
  ggplot(aes(
    x = connectivity_us,
    y = .epred,
    fill_ramp = after_stat(.width)
  )) +
  ggdist::stat_lineribbon(.width = ppoints(50), fill = "#56B4E9", linewidth = 0.25) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "#56B4E9")) +
  geom_point(
    data = testdat_t,
    aes(
      x = connectivity,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.5, size = 0.01,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0,1)) +
  xlab("Neighbourhood density") +
  ylab("Proportion of immature seeds") +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") -> p2


hetro_mod %>%
  ggplot(aes(
    x = connectivity_us,
    y = .epred,
    fill_ramp = after_stat(.width)
  )) +
  ggdist::stat_lineribbon(.width = ppoints(50), fill = "#009E73", linewidth = 0.25) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "#009E73")) +
  geom_point(
    data = testdat_h,
    aes(
      x = connectivity,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.5, size = 0.01,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0,1)) +
  xlab("Neighbourhood density") +
  ylab("Proportion of immature seeds") +
  theme(legend.position = "none",
        axis.title.y = element_text(hjust = 0),
        axis.title.x = element_text(hjust = 0)) -> p3



# pd ----------------------------------------------------------------------
testdat %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 70, expand =2)
  ) %>%
  add_epred_draws(model, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat$connectivity_sc, 'scaled:scale') +
      attr(testdat$connectivity_sc, 'scaled:center')
  ) -> repro_mod2

testdat_t %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 70, expand =2)
  ) %>%
  add_epred_draws(model_t, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat_t$connectivity_sc, 'scaled:scale') +
      attr(testdat_t$connectivity_sc, 'scaled:center')
  ) -> total_mod2

testdat_h %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 70, expand = 2)
  ) %>%
  add_epred_draws(model_h, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat_h$connectivity_sc, 'scaled:scale') +
      attr(testdat_h$connectivity_sc, 'scaled:center')
  ) -> hetro_mod2

bind_rows(repro_mod2, total_mod2, hetro_mod2,  .id = "dataset") %>%
  ggplot(aes(x = connectivity_us, y = .epred, colour = dataset, fill = dataset)) +
  stat_lineribbon(.width = .95, alpha = 0.5, linewidth = 0.5) +
  stat_lineribbon(.width = 0, alpha = 1, linewidth = 0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("Reproductive sized conspecifics",
                               "Non-reproductive sized conspecifics",
                               "Reproductive sized heterospecifics")) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_classic(base_size = 10) +
  coord_cartesian(xlim = c(0,2), ylim = c(0, 1), expand = FALSE) +
  xlab("Neighbourhood density") +
  ylab("Proportion of immature seeds") +
  guides(colour = "none",
         fill = guide_legend(override.aes = list(alpha = 1, linewidth = 0))) +
  theme(legend.title = element_blank(), legend.position = c(0.35, 0.9),
        axis.title.y = element_text(hjust = 0),
        axis.title.x = element_text(hjust = 0)) -> p4

((p1 / p2 / p3) | p4 ) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag.position  = c(.935, .96))

# save as png
png(
  here::here("output", "figures", "figure2.png"),
  width = 6.81102,
  height = 4.6145,
  units = "in",
  type = "cairo",
  res = 600
)

((p1 / p2 / p3) | p4 ) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag.position  = c(.935, .96))

dev.off()
