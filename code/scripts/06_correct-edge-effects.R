#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: correct-edge-effects
## Desc: adjust values of connectivity for traps that are at the edge
## Date created: 2020-11-05

# Load packages ---------------------------

library("tidyverse") # v 1.3.1
library("here") # v 1.0.1
library("sf")

# Load data ---------------------------

readRDS(here::here("data", "clean", "trap_connect.rds")) %>%
  rename(connectivity_obs = connectivity) -> trap_connect

edge_traps <- subset(trap_connect, x > 880 | x < 120 | y > 380 | y < 120)

correct_for_edge <- function(species, yr){
	# draw the fdp
	fdplot <- rbind(c(0,0), c(0, 500), c(1000, 500), c(1000, 0), c(0,0))
	# make it a polygon
	fdplot.sf <- sf::st_polygon(list(fdplot))
	# subset edge_traps to sp i and yr i
	dplyr::filter(edge_traps,
		(.data[["sp4"]] == species) & (.data[["year"]] == yr)) %>%
	# make traps sf points
	sf::st_as_sf(coords = c("x", "y")) %>%
	# create a radius around the traps
	mutate(buf = st_buffer(geometry, dist = 120)) %>%
	# calculate the area of intersection between each radius & the fdp
	mutate(area = sf::st_area(
		sf::st_intersection(.data[["buf"]], fdplot.sf) )) %>%
	# calculate our prediction
	mutate(connectivity_pred = (connectivity_obs/area) * (pi * 120^2))
}

# pull out every combo of yr and sp in edge_traps
edge_traps %>%
	count(year, sp4) %>%
	select(sp4, year) %>%
	rename(species = sp4, yr = year) -> keys

# map function over rows in keys
edge_traps_adjusted <- pmap(keys, correct_for_edge)

# bind output into one df
edge_traps_adjusted_b <- bind_rows(edge_traps_adjusted)

# Merge and save dataset ---------------------------

as.data.frame(edge_traps_adjusted_b) %>%
	select(year, sp4, trap, connectivity_pred) -> edge_traps_adjusted_s

# connectivity_pred for non-edge traps are NA -> fill with connectivity_obs
trap_connect %>%
  left_join(edge_traps_adjusted_s, by = c("year", "sp4", "trap")) %>%
	mutate(connectivity_pred = coalesce(connectivity_pred, connectivity_obs)) -> trap_connect_corrected

trap_connect_corrected %>%
  select(- x, - y, - capsules) %>%
  saveRDS(here::here("data", "clean", "data_for_model.rds"))
