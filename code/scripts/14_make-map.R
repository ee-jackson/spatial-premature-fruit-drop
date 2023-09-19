#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: make-map
## Desc:
## Date created: 2023-07-18


# packages ----------------------------------------------------------------

library("here")
library("tidyverse")
library("ggmap")
library("ggsn")
library("sf")
library("patchwork")

traps <- read.csv(here::here("data", "clean", "trap_locations.csv"))


# 50ha plot map -----------------------------------------------------------

traps %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(
    size = 1,
    alpha = 0.7,
    shape = 16,
    colour = "forestgreen") +
  theme_classic(base_size = 10) +
  labs(x = "1000m", y = "500m") +
  theme(
    line = element_blank(),
    axis.title.y = element_text(angle = 90),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_rect(colour = "blue", fill=NA, size = 1)
    ) +
  theme(plot.background = element_rect(colour = NA,
                                    fill = NA,
                                    linewidth = 1)) -> traps_map


# BCI map -----------------------------------------------------------------

read_sf(here::here("data", "raw", "BCI_Plot_50ha", "BCI_Plot_50ha.shp")) %>%
  st_transform(plot_50ha, crs = st_crs(4326)) -> plot_50ha

bbox <- make_bbox(c(-79.875, -79.75),
                  c(9.10, 9.185))

bci_basemap <- ggmap::get_map(bbox, source = "stamen",
                              force = TRUE, maptype = "toner-lite")

bci_basemap %>%
  ggmap() +
  geom_sf(data = plot_50ha, inherit.aes = FALSE, fill = NA,
          colour = "blue", linewidth = 0.5) +
  theme_void() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black",
                                fill = NA,
                                linewidth = 0.5)) -> bci_plot


# combine -----------------------------------------------------------------

bci_plot + inset_element(traps_map, left = 0.35, bottom = 0.005,
                           right = 0.995, top = 0.6,
                           align_to = "plot", clip = TRUE) -> final_map

ggsave(here::here("output", "figures", "bci_map.png"),
       plot = final_map, device = "png", dpi = 600,
       width = 110, height = 150, units = "mm")

trap_png <- png::readPNG(here::here("output", "figures", "trap.png"),
                                     native = TRUE)

map_png <- png::readPNG(here::here("output", "figures", "bci_map.png"),
                         native = TRUE)

(ggplot() + theme_void()) / map_png / trap_png +
  plot_layout(heights = c(0, 2, 2)) +
  plot_annotation(tag_levels = list(c('', 'a', 'b')))

ggsave(here::here("output", "figures", "bci_map_trap.png"),
       device = "png", dpi = 600, width = 110, height = 150, units = "mm")
