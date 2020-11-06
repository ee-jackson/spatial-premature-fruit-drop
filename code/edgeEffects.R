#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: edgeEffects
## Desc: adjust values of CI for traps that are at the edge
## Date: November 2020

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size=10))
library("sf")
load("../data/clean/trapConnect2.RData")


###############################################################################
## test
###############################################################################

# trial on one yr and sp, subset
dat.jac <- subset(trapConnect, sp == "jac1co" & year == "2015")

# subset to traps within 50m of the edge
dat.jac.edgetraps <- subset(dat.jac, X > 950 | X < 50 | Y > 450 | Y < 50)

# make the traps sf points
trap.sf <- st_as_sf(dat.jac.edgetraps, coords = c("X", "Y"))

# make polygons in a 50m radius around the traps
trap.sf$buf <- st_buffer(trap.sf$geometry, dist = 50)

# draw the fdp and make it a polygon
fdplot <- rbind(c(0,0), c(0, 500), c(1000, 500), c(1000, 0),c(0,0))
fdplot.sf <- st_polygon(list(fdplot))

# what's the intersection between each radius and the fdp?
st_intersection(trap.sf$buf, fdplot.sf) %>%
	ggplot() + geom_sf()

###############################################################################
## make iterable
###############################################################################

edgeTraps <- subset(trapConnect, X > 950 | X < 50 | Y > 450 | Y < 50)

correct_for_edge <- function(species, yr){
	# draw the fdp
	fdplot <- rbind(c(0,0), c(0, 500), c(1000, 500), c(1000, 0),c(0,0))
	fdplot.sf <- st_polygon(list(fdplot)) # make it a polygon
	dplyr::filter(edgeTraps, (.data[["sp"]]==species) & (.data[["year"]]==yr)) %>% # subset edgeTraps to sp i and yr i
	sf::st_as_sf(coords = c("X", "Y")) %>% # make traps sf points
	# create a 50m radius around the traps
	mutate(buf = st_buffer(geometry, dist = 50)) %>% 
	# calculate the area of intersection between each trap & the fdp
	mutate(area = sf::st_area( sf::st_intersection(.data[["buf"]], fdplot.sf) )) %>%
	# calculate
	mutate(CI.edge.adjusted = (CI/area) * (pi * 50^2) )
}

# pull out every combo of yr and sp in edgeTraps
edgeTraps %>% 
	count(year, sp) %>% 
	select(sp, year) %>% 
	rename(species = sp, yr = year) -> keys

# map function over rows in keys
edgeTraps.adjusted <- pmap(keys, correct_for_edge)

# bind output into one df
edgeTraps.adjusted.b <- bind_rows(edgeTraps.adjusted)

###############################################################################
## Merge and save dataset
###############################################################################

edgeTraps.adjusted.b %>%
	select(year, sp, trap, CI.edge.adjusted) -> edgeTraps.adjusted.s

left_join(trapConnect, edgeTraps.adjusted.s, by = c("year", "sp", "trap")) %>%
	mutate(CI.edge.adjusted = coalesce(CI.edge.adjusted, CI)) -> trapConnect.a2

save(trapConnect.a2, sp.list, file = "../data/clean/trapConnect2.edge.RData")
