#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: make-fig-04
## Desc: make figure 4

# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")
library("patchwork")
library("ggdist")
library("ggtext")

# set the ggplot theme
theme_set(
  theme_bw(base_size = 8) +
    theme(
      panel.spacing = unit(0, "lines"))
)

# Get model ---------------------------------------------------------------

mod <-
  readRDS(here::here("output", "models", "pheno-repro-adjust",
                     "full_conn_binom_nseeds_abund.rds"))

# function to return variables to measurement scale
mod_data <-
  readRDS("data/clean/trap_connect_buffer.rds")

unscale <- function(x, var_name) {
  x * attr(mod_data[[var_name]], "scaled:scale") +
    attr(mod_data[[var_name]], "scaled:center")
}


# Make panel a ------------------------------------------------------------

cond_eff <-
  conditional_effects(mod, effects = "log_median_abundance_sc")

p1 <-
  cond_eff[["log_median_abundance_sc"]] %>%
  mutate(log_median_abundance =
           unscale(log_median_abundance_sc, "log_median_abundance_sc_mat")) %>%
  ggplot(aes(x = log_median_abundance, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = "#56B4E9", alpha = 0.5) +
  geom_line(linewidth = 1, colour = "#2F6380") +
  labs(x = "Species abundance
       <span style='font-size:5pt'>(log basal area m<sup>2</sup>)</span>",
       y = "Proportion of immature seeds") +
  coord_cartesian(ylim = c(0,1),
                  expand = 0) +
  theme(axis.title.x = element_markdown())


# Make pabel b ------------------------------------------------------------

draws <-
  as_draws_df(mod)

draws %>%
  mutate(
    abund_effect_on_RC_vs_RH =
      `b_conn_RC_sc:log_median_abundance_sc` -
      `b_conn_RH_sc:log_median_abundance_sc`,
    .keep = "none"
  ) %>%
  summarise(
    mean = mean(abund_effect_on_RC_vs_RH),
    sd = sd(abund_effect_on_RC_vs_RH),
    q025 = quantile(abund_effect_on_RC_vs_RH, 0.025),
    q50  = quantile(abund_effect_on_RC_vs_RH, 0.5),
    q975 = quantile(abund_effect_on_RC_vs_RH, 0.975),
    p_gt_0 = mean(abund_effect_on_RC_vs_RH > 0),
    p_lt_0 = mean(abund_effect_on_RC_vs_RH < 0)
  )

# there is one abundance value per species
species_abund <-
  mod$data %>%
  distinct(sp4, log_median_abundance_sc)

# extract species random slopes for RC and RH
sp_re <-
  draws %>%
  select(.draw, matches("^r_sp4\\[")) %>%
  pivot_longer(
    cols = -.draw,
    names_to = "param",
    values_to = "value"
  ) %>%
  filter(str_detect(param, ",conn_RC_sc\\]$|,conn_RH_sc\\]$")) %>%
  mutate(
    sp4 = str_match(param, "^r_sp4\\[(.*?),")[, 2],
    term = str_match(param, ",(conn_RC_sc|conn_RH_sc)\\]$")[, 2]
  ) %>%
  select(.draw, sp4, term, value) %>%
  pivot_wider(names_from = term, values_from = value)

# extract the fixed effects for RC and RH
fixef_rc_rh <-
  draws %>%
  mutate(
    .draw,
    b_RC = .data[["b_conn_RC_sc"]],
    b_RH = .data[["b_conn_RH_sc"]],
    b_RC_seed = .data[["b_conn_RC_sc:log_total_seeds_sc"]],
    b_RH_seed = .data[["b_conn_RH_sc:log_total_seeds_sc"]],
    b_RC_abund = .data[["b_conn_RC_sc:log_median_abundance_sc"]],
    b_RH_abund = .data[["b_conn_RH_sc:log_median_abundance_sc"]],
    .keep = "none"
  )

# fix log total seeds at zero (the mean)
L0 <- 0

# compute the RC-RH contrast per draw
species_rc_rh_draws <-
  fixef_rc_rh %>%
  left_join(sp_re, by = ".draw") %>%
  left_join(species_abund, by = "sp4") %>%
  mutate(
    rc_rh_contrast =
      (b_RC - b_RH) +
      (b_RC_seed - b_RH_seed) * L0 +
      (b_RC_abund - b_RH_abund) * log_median_abundance_sc +
      (conn_RC_sc - conn_RH_sc)
  ) %>%
  mutate(log_median_abundance =
           unscale(log_median_abundance_sc, "log_median_abundance_sc_mat")) %>%
  select(.draw, sp4, log_median_abundance, rc_rh_contrast)


p2 <-
  species_rc_rh_draws %>%
  ggplot(aes(x = log_median_abundance, group = sp4,
             y = rc_rh_contrast)) +
  stat_interval(aes(interval_alpha = after_stat(level)),
                size = 1.5, slab_alpha = 0.8,
                interval_colour = "#56B4E9", .width = c(0.95, 0.8, 0.5)) +
  stat_pointinterval(size = 0.15, colour = "#2F6380",
                      .width = c(0.95), interval_alpha = 0.8,
                     shape = 21, stroke = 0.5,
                     point_fill = "white") +
  labs(y = "Stabilising CNDD effect",
       x = "Species abundance
       <span style='font-size:5pt'>(log basal area m<sup>2</sup>)</span>") +
  geom_hline(yintercept = 0, linetype = 1, colour = "#D55E00", linewidth = 0.5) +
  coord_cartesian(ylim = c(-2.5, 5)) +
  theme(legend.position = "none",
        axis.title.x = element_markdown())

p2


# Combine panels and save -------------------------------------------------

png(
  here::here("output", "figures", "figure_04.png"),
  width = 110,
  height = 90,
  units = "mm",
  type = "cairo",
  res = 600
)

p1 + inset_element(p2, left = 0.05, bottom = 0.4, right = 0.7, top = 0.95,
                   align_to = "panel") +
  plot_annotation(tag_levels = "a")

dev.off()
