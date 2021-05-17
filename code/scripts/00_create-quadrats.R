#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: create-quadrats
## Desc: assign each trap to a 50m2 quadrat
## Date created: 2020-10-07

library("tidyverse"); theme_set(theme_bw(base_size = 8))
library("sf")

# co-ordinates of the original 200 traps
trap_loc <- read.csv(here::here("data", "raw", "trapLocations.csv")) 

# paired traps - angles at 2m from the original traps
pair_traps <- read.csv(here::here("data", "raw", "pairedTraps.csv"))

# format trap ids
trap_loc$trap <- formatC(trap_loc$trap, width = 3, format = "d", flag = "0")
trap_loc$trap <- paste("trap", trap_loc$trap, sep="_")

pair_traps$trap <- formatC(pair_traps$trap, 
	width = 3, format = "d", flag = "0")
pair_traps$trap <- paste("trap", pair_traps$trap, sep="_")

pair_traps$paired_trap <- formatC(pair_traps$paired_trap, 
	width = 3, format = "d", flag = "0")
pair_traps$paired_trap <- paste("trap", pair_traps$paired_trap, sep="_")

# add coordinates from trap_loc to paired traps
pair_traps %>%
	left_join(trap_loc, by = "trap") %>% 
	select(-trap) %>%
	rename(trap = paired_trap) -> pair_traps_loc

# first convert from degrees to radians
# then calculate coordinates of paired traps
# and then bind
pair_traps_loc %>% 
	mutate(angle = angle * (pi/180)) %>%
	mutate(X = X + cos(angle) * 2, Y = Y + sin(angle) * 2) %>% 
	select(-angle) %>%
	bind_rows(trap_loc) -> all_trap_loc

# convert to sf object
trap_sf <- st_as_sf(all_trap_loc, coords = c("X", "Y"))

# change bounding box
new_bb <- c(0, 0, 1000, 500)
names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
attr(new_bb, "class") = "bbox"
attr(st_geometry(trap_sf), "bbox") = new_bb

# create quadrats of 50m2
quads <- st_make_grid(
  trap_sf,
  cellsize = 50,
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)

#convert polygons to sf object and add id column
quads %>% 
	st_sfc() %>% 
	st_sf() %>% 
	mutate(quadrat = paste(seq(1:200))) -> quads_sf 

# calculate what traps are in what quadrats
joined <- st_intersection(trap_sf, quads_sf)

# plot to check
ggplot() + 
	geom_sf(data = quads_sf) +
	geom_sf(data = joined, aes(col=quadrat)) 

# convert to dataframe
joined %>%
	st_drop_geometry() %>%
	left_join(all_trap_loc, by = "trap") -> trap_loc_quad

# format quadrat id
trap_loc_quad$quadrat <- formatC(as.numeric(trap_loc_quad$quadrat),width = 3, format = "d", flag = "0")

trap_loc_quad$quadrat <- paste("quadrat", trap_loc_quad$quadrat, sep="_")

# if traps are on the boundary between 2 quadrats they are duplicated
# remove duplicates
trap_loc_quad[duplicated(trap_loc_quad$trap),]
trap_loc_quad <- trap_loc_quad[!duplicated(trap_loc_quad$trap),]

save(trap_loc_quad, file = here:here("data", "clean", "trapLocations.RData"))
