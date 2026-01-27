#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script:
## Desc:
## Date created: 2023-07-18


# Packages ----------------------------------------------------------------

library("tidyverse")
library("patchwork")
library("here")
library("marginaleffects")
library("tidybayes")
library("ggtext")
library("ggpubr")


# Get abundance data ------------------------------------------------------

tree_data <-
  readRDS(here::here("data", "clean", "tree_data.rds")) %>%
  filter(year %in% c("1990", "1995", "2000", "2005", "2010",
                     "2010", "2015", "2022"))

trap_data <- readRDS(here::here("data", "clean", "trap_data.rds"))


monoecious_species <-
  read.csv(here::here("data", "clean", "species_list.csv"))

# only keep species which appear in both datasets
shared_sp <-
  trap_data %>%
  select(sp4) %>%
  distinct() %>%
  inner_join(
    y = tree_data %>%
      select(sp4) %>%
      distinct()
  ) #%>%
  #filter(sp4 %in% monoecious_species$sp4) # drops 21 dioecious species

tree_data <-
  tree_data %>%
  filter(sp4 %in% shared_sp$sp4)

abun <-
  tree_data %>%
  group_by(year, sp4, genus, species) %>%
  summarise(abundance = n_distinct(tree, na.rm = TRUE)) %>%
  group_by(sp4, genus, species) %>%
  summarise(median_abundance = median(abundance, na.rm = TRUE),
            mean_abundance = mean(abundance, na.rm = TRUE),
            max_abundance = max(abundance, na.rm = TRUE)) %>%
  mutate(log_median_abundance = log(median_abundance))


# abundance distribution panel --------------------------------------------

q25 <-
  quantile(abun$log_median_abundance) %>%
  pluck(2)

q50 <-
  quantile(abun$log_median_abundance) %>%
  pluck(3)

q75 <-
  quantile(abun$log_median_abundance) %>%
  pluck(4)

abdun_dist_fig <-
  abun %>%
  mutate(grouped_abun = case_when(
    log_median_abundance < q25 ~ "1Rare",
    log_median_abundance > q25  & log_median_abundance < q50 ~ "2Uncommon",
    log_median_abundance > q50  & log_median_abundance < q75 ~ "3Common",
    log_median_abundance > q75 ~ "4Abundant"
  )) %>%
  ggplot(aes(x = log_median_abundance)) +
  ggdist::stat_halfeye(
    aes(fill = after_stat(case_when(
      x < q25 ~ "q1",
      x > q25  & x < q50 ~ "q2",
      x > q50  & x < q75 ~ "q3",
      x > q75 ~ "q4"
    ))),
    slab_colour = "black",
    point_interval = median_qi,
    adjust = .5,
    width = .6,
    .width =  0,
    justification = -.3,
    point_colour = NA,
    slab_alpha = 0.6,slab_linewidth = 0.7) +
  geom_boxplot(
    width = .25,
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_vline(aes(
    xintercept = q25
    ),
    linetype = 2) +
  geom_vline(aes(
    xintercept = q50
    ),
    linetype = 2) +
  geom_vline(aes(
    xintercept = q75
    ),
    linetype = 2) +
  labs(x = "Species abundance<br>log count") +
  theme_classic(base_size = 7) +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer("Blues",
                    direction = 1)


# predicted fruit drop panel ----------------------------------------------

# get connectivity data
consp_repro <-
  readRDS(here::here("data", "clean", "trap_connect_repro_consp_20m_dioecious.rds"))

# don't include traps < 20m from the edge of the plot
# centre and scale connectivity
test_data <-
  consp_repro %>%
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) %>%
  select(- x, - y, - capsules) %>%
  transform(connectivity_sc =
              scale(connectivity)) %>%
  filter(sum_parts >= 3) %>%
  mutate_at(c("sp4", "year", "trap", "quadrat"), ~as.factor(.))

model <-
  readRDS(here::here("output", "models",
                     "repro_consp_20m_yesdioecious.rds"))

model_h <-
  readRDS(here::here("output", "models",
                     "repro_hetero_20m_yesdioecious.rds"))

# get predictions
preds_sp_con <-
  test_data %>%
  modelr::data_grid(
    connectivity_sc =
      c(0, 0.5, 1, 1.5, 2),
    sp4 = unique(test_data$sp4),
    .model = model
  ) %>%
  tidybayes::add_epred_draws(model,
                             re_formula =
                               ~ (1 +
                                    connectivity_sc | sp4))

preds_sp_het <-
  test_data %>%
  modelr::data_grid(
    connectivity_sc =
      c(0, 0.5, 1, 1.5, 2),
    sp4 = unique(test_data$sp4),
    .model = model_h
  ) %>%
  tidybayes::add_epred_draws(model_h,
                             re_formula =
                               ~ (1 +
                                    connectivity_sc | sp4))


# panel 2 -----------------------------------------------------------------

pred_fig <-
  preds_sp_con %>%
  point_interval() %>%
  bind_rows(point_interval(preds_sp_het),
            .id = "type") %>%
  mutate(type = ifelse(type == 1, "Reproductive conspecifics",
                       "Reproductive heterospecifics")) %>%
  left_join(abun) %>%
  mutate(grouped_abun = case_when(
    log_median_abundance < q25 ~ "1Rare",
    log_median_abundance > q25  & log_median_abundance < q50 ~ "2Uncommon",
    log_median_abundance > q50  & log_median_abundance < q75 ~ "3Common",
    log_median_abundance > q75 ~ "4Abundant"
  )) %>%
  ggplot(aes(x = connectivity_sc,
             y = .epred,
             colour = type
  )) +
  facet_grid(type~grouped_abun) +
  geom_rect(data = tibble(grouped_abun =
                            c("1Rare", "2Uncommon", "3Common", "4Abundant")),
            aes(fill = grouped_abun),
            inherit.aes = FALSE, show.legend = FALSE, alpha = 0.5,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(alpha = 0.6,
             size = 0.5,
             shape = 16) +
  geom_smooth(method = "lm",
              show.legend = FALSE,linewidth = 0.5) +
  stat_regline_equation(colour = "grey20",
                        size = 2) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Predicted immature seed drop",
       x = "Neighbourhood density") +
  scale_fill_brewer("Blues",
                    direction = 1) +
  scale_colour_manual(values = c("#E69F00", "#009E73"))


# panel 3 -----------------------------------------------------------------

preds_difference <-
    preds_sp_con %>%
    left_join(preds_sp_het,
              by = c("connectivity_sc", "sp4",
                     ".row", ".chain",
                     ".iteration", ".draw"),
              suffix = c("_conspecific", "_heterospecific")) %>%
    mutate(difference =  .epred_conspecific - .epred_heterospecific) %>%
    point_interval()

pred_fig_diff <-
  preds_difference %>%
  left_join(abun) %>%
  mutate(grouped_abun = case_when(
    log_median_abundance < q25 ~ "1Rare",
    log_median_abundance > q25  & log_median_abundance < q50 ~ "2Uncommon",
    log_median_abundance > q50  & log_median_abundance < q75 ~ "3Common",
    log_median_abundance > q75 ~ "4Abundant"
  )) %>%
  ggplot(aes(x = connectivity_sc,
             y = difference
  )) +
  facet_wrap(~grouped_abun, nrow = 1) +
  geom_rect(data = tibble(grouped_abun = c("1Rare", "2Uncommon", "3Common", "4Abundant")),
            aes(fill = grouped_abun),
            inherit.aes = FALSE, show.legend = FALSE, alpha = 0.5,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(alpha = 0.5,
             shape = 16,
             size = 0.5) +
  geom_smooth(method = "lm",linewidth = 0.5) +
  stat_regline_equation(colour = "grey20",
                        size = 2) +
  geom_hline(yintercept = 0,
             colour = "red",
             linetype = 2) +
  theme_classic(base_size = 7) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Difference in\npredicted immature seed drop",
       x = "Neighbourhood density") +
  scale_fill_brewer("Blues",
                    direction = 1)


abdun_dist_fig /
  pred_fig /
  pred_fig_diff +
  plot_annotation(tag_levels = "a") +
  patchwork::plot_layout(heights = c(1,2,1))

pdf(
  here::here("output", "figures", "figure3.pdf"),
  width = 4.33071,
  height = 5.90551
)
abdun_dist_fig /
  pred_fig /
  pred_fig_diff +
  plot_annotation(tag_levels = "a") +
  patchwork::plot_layout(heights = c(1,2,1))
dev.off()
