#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: make-fig-03
## Desc: make figure 3

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

sample_size <-
  mod$data %>%
  group_by(sp4) %>%
  summarise(n = n())

sp_names <-
  read_csv(here::here("data", "clean", "species_list.csv")) %>%
  left_join(sample_size) %>%
  mutate(genus_species = paste("<i>", genus, species, "</i>", n, sep = " ")) %>%
  select(sp4, genus_species)

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
  select(.draw, sp4, log_median_abundance_sc, rc_rh_contrast)

species_rc_rh_draws_medians <-
  species_rc_rh_draws %>%
  group_by(sp4) %>%
  point_interval(rc_rh_contrast) %>%
  select(sp4, rc_rh_contrast) %>%
  rename(median = rc_rh_contrast)

p <-
  species_rc_rh_draws %>%
  left_join(sp_names) %>%
  left_join(species_rc_rh_draws_medians) %>%
  group_by(genus_species) %>%
  ggplot(aes(y = reorder(genus_species, median),
             x = rc_rh_contrast,
             group = genus_species)) +
  ggdist::stat_gradientinterval(.width = 0.95, fill = "#56B4E9",
                                stroke = 0.5, linewidth = 0.5,
                                shape = 21, fatten_point = 0.7,
                                point_fill = "white",
                                fill_type = "gradient",
                                point_interval = median_qi,
                                slab_linewidth = 0, normalize = "groups"
                                ) +
  scale_slab_alpha_continuous(range = c(0, 1), limits = c(0, 0.95)) +
  labs(x = "Stabilising CNDD effect", y = "") +
  coord_cartesian(xlim = c(-2.5, 5)) +
  geom_vline(xintercept = 0, linetype = 1, colour = "#D55E00", linewidth = 0.25) +
  theme(axis.text.y = element_markdown())

png(
  here::here("output", "figures", "figure_03.png"),
  width = 82,
  height = 200,
  units = "mm",
  type = "cairo",
  res = 600
)

p

dev.off()
