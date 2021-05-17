#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: correct-edge-effects
## Desc: adjust values of CI for traps that are at the edge
## Date created: 2020-11-05

library("tidyverse"); theme_set(theme_bw(base_size=8))
library("sf")
load(here:here("data", "clean", "trapConnect.RData"))

trapConnect %>% rename(CI_obs = CI) -> trapConnect

edgeTraps <- subset(trapConnect, X > 900 | X < 100 | Y > 400 | Y < 100)

correct_for_edge <- function(species, yr){
	# draw the fdp
	fdplot <- rbind(c(0,0), c(0, 500), c(1000, 500), c(1000, 0), c(0,0))
	# make it a polygon
	fdplot.sf <- sf::st_polygon(list(fdplot))
	# subset edgeTraps to sp i and yr i
	dplyr::filter(edgeTraps, 
		(.data[["SP4"]]==species) & (.data[["year"]]==yr)) %>%
	# make traps sf points
	sf::st_as_sf(coords = c("X", "Y")) %>% 
	# create a radius around the traps
	mutate(buf = st_buffer(geometry, dist = 100)) %>% 
	# calculate the area of intersection between each radius & the fdp
	mutate(area = sf::st_area( 
		sf::st_intersection(.data[["buf"]], fdplot.sf) )) %>% 
	# calculate our prediction
	mutate(CI_pred = (CI_obs/area) * (pi * 100^2))
}

# pull out every combo of yr and sp in edgeTraps
edgeTraps %>% 
	count(year, SP4) %>% 
	select(SP4, year) %>% 
	rename(species = SP4, yr = year) -> keys

# map function over rows in keys
edgeTraps_adjusted <- pmap(keys, correct_for_edge)

# bind output into one df
edgeTraps_adjusted_b <- bind_rows(edgeTraps_adjusted)

###############################################################################
## Merge and save dataset
###############################################################################

as.data.frame(edgeTraps_adjusted_b) %>%
	select(year, SP4, trap, CI_pred) -> edgeTraps_adjusted_s

# CI_pred for non-edge traps are NA -> fill with CI_obs
left_join(trapConnect, edgeTraps_adjusted_s, by = c("year", "SP4", "trap")) %>%
	mutate(CI_pred = coalesce(CI_pred, CI_obs)) -> trapConnect

save(trapConnect, sp.list, here:here("data", "clean", "cleanData.RData"))
