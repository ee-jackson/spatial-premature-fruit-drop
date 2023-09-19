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

data <- readRDS(here::here("data", "clean", "trap_connect_20m.rds"))
model <- readRDS(here::here("output", "models", "202308", "zoib_20m.rds"))

data_t <- readRDS(here::here("data", "clean", "total_trap_connect_20m.rds"))
model_t <- readRDS(here::here("output", "models", "202308", "zoib_total_c_20m.rds"))

data_h <- readRDS(here::here("data", "clean", "hetero_trap_connect_20m.rds"))
model_h <- readRDS(here::here("output", "models", "202308", "zoib_hetero_20m.rds"))

sp_data <- read.csv(here::here("data", "clean", "species_list.csv"))


# format / transform data -------------------------------------------------

edge_traps <- subset(data, x > 980 | x < 20 | y > 480 | y < 20)

edge_traps %>%
  select(trap) %>%
  distinct() -> edge_traps_list

data %>%
  filter(!trap %in% edge_traps_list$trap) %>%
  select(- x, - y, - capsules) %>%
  transform(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat_20m

testdat_20m %>%
  mutate(
         year = as.factor(year),
         trap = as.factor(trap)) %>%
  mutate(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat

##

data_t %>%
  filter(!trap %in% edge_traps_list$trap) %>%
  select(- x, - y, - capsules) %>%
  transform(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat_20m

testdat_20m %>%
  mutate(
    year = as.factor(year),
    trap = as.factor(trap)) %>%
  mutate(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat_t

##

data_h %>%
  filter(!trap %in% edge_traps_list$trap) %>%
  select(- x, - y, - capsules) %>%
  transform(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat_20m

testdat_20m %>%
  mutate(
    year = as.factor(year),
    trap = as.factor(trap)) %>%
  mutate(connectivity_sc = scale(connectivity)) %>%
  filter(sum_parts >= 3) -> testdat_h


# -------------------------------------------------------------------------
# get predictions

testdat %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 101)
  ) %>%
  add_epred_draws(model, ndraws = 1000, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat$connectivity_sc, 'scaled:scale') +
      attr(testdat$connectivity_sc, 'scaled:center')
  ) -> repro_mod

testdat %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 101)
  ) %>%
  add_epred_draws(model_t, ndraws = 1000, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat$connectivity_sc, 'scaled:scale') +
      attr(testdat$connectivity_sc, 'scaled:center')
  ) -> total_mod

testdat %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 101)
  ) %>%
  add_epred_draws(model_h, ndraws = 1000, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat$connectivity_sc, 'scaled:scale') +
      attr(testdat$connectivity_sc, 'scaled:center')
  ) -> hetro_mod

# -------------------------------------------------------------------------
# big 4 panel plot

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
  ggdist::stat_lineribbon(.width = ppoints(50), fill = "#E69F00", linewidth = 0.5) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "#E69F00")) +
  geom_point(
    data = testdat,
    aes(
      x = connectivity,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.6, size = 1,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 25) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0,1)) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") -> p1

testdat_t %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 101)
  ) %>%
  add_epred_draws(model_t, ndraws = 1000, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat_t$connectivity_sc, 'scaled:scale') +
      attr(testdat_t$connectivity_sc, 'scaled:center')
  ) %>%
  ggplot(aes(
    x = connectivity_us,
    y = .epred,
    fill_ramp = after_stat(.width)
  )) +
  ggdist::stat_lineribbon(.width = ppoints(50), fill = "#56B4E9", linewidth = 0.5) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "#56B4E9")) +
  geom_point(
    data = testdat_t,
    aes(
      x = connectivity,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.6, size = 1,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 25) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0,1)) +
  xlab("Connectivity") +
  ylab("Proportion of immature seeds") +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") -> p2

testdat_h %>%
  modelr::data_grid(
    connectivity_sc = modelr::seq_range(connectivity_sc, n = 101)
  ) %>%
  add_epred_draws(model_h, ndraws = 1000, re_formula = NA) %>%
  mutate(
    connectivity_us = connectivity_sc *
      attr(testdat_h$connectivity_sc, 'scaled:scale') +
      attr(testdat_h$connectivity_sc, 'scaled:center')
  ) %>%
  ggplot(aes(
    x = connectivity_us,
    y = .epred,
    fill_ramp = after_stat(.width)
  )) +
  ggdist::stat_lineribbon(.width = ppoints(50), fill = "#009E73", linewidth = 0.5) +
  scale_fill_ramp_continuous(range = c(1, 0), guide = guide_rampbar(to = "#009E73")) +
  geom_point(
    data = testdat_h,
    aes(
      x = connectivity,
      y = proportion_abscised
    ),
    inherit.aes = FALSE,
    alpha = 0.6, size = 1,
    shape = 16, colour = "grey25"
  ) +
  theme_classic(base_size = 25) +
  scale_x_continuous(expand = c(0.003, 0.003)) +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0,1)) +
  xlab("Connectivity") +
  ylab("Proportion of immature seeds") +
  theme(legend.position = "none",
        axis.title.y = element_text(hjust = 0),
        axis.title.x = element_text(hjust = 0)) -> p3



bind_rows(repro_mod, total_mod, hetro_mod,  .id = "dataset") %>%
  ggplot(aes(x = connectivity_us, y = .epred, color = dataset, fill = dataset)) +
  stat_lineribbon(.width = .95, alpha = 0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("Reproductive conspecifics",
                               "All conspecifics",
                               "Reproductive heterospecifics")) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  stat_lineribbon(.width = 0, alpha = 1) +
  theme_classic(base_size = 30) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits= c(0,1), expand = c(0, 0)) +
  xlab("Connectivity") +
  ylab("Proportion of immature seeds") +
  guides(colour = FALSE,
         fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.title = element_blank(), legend.position = c(0.25, 0.9),
        axis.title.y = element_text(hjust = 0),
        axis.title.x = element_text(hjust = 0)) -> p4

((p1 / p2 / p3) | p4 ) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag.position  = c(.935, .96))

png(
  here::here("output", "figures", "all-models-predict-big-labs.png"),
  width = 1476,
  height = 1000,
  units = "px",
  type = "cairo"
)

((p1 / p2 / p3) | p4 ) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag.position  = c(.935, .96))

dev.off()
