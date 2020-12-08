#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: createQuadrats
## Desc: assign each trap to a 50m2 quadrat
## Date: 20201007

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size=10))
library("sf")

trapLoc <- read.csv("../data/clean/trapLocations.csv")

trapLoc$trap <- formatC(trapLoc$trap, width = 3, format = "d", flag = "0")
trapLoc$trap <- paste("trap", trapLoc$trap, sep="_")

# convert to sf object
trap.sf <- st_as_sf(trapLoc, coords = c("X", "Y"))

# change bounding box
new_bb <- c(0, 0, 1000, 500)
names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
attr(new_bb, "class") = "bbox"
attr(st_geometry(trap.sf), "bbox") = new_bb

# create quadrats of 50m2
st_make_grid(
  trap.sf,
  cellsize = 50,
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE
) -> quads

#convert polygons to sf object and add id column
quads.sf <- quads %>% st_sfc() %>% st_sf(geom=.) %>% 
	mutate(quadratID = paste("quadrat", seq(1:200), sep="_")) 

# calculate what traps are in what quadrats
joined <- st_intersection(trap.sf, quads.sf)

# plot to check
ggplot() + geom_sf(data=quads.sf) +
geom_sf(data=joined ,aes(col=quadratID)) 

# convert to dataframe
joined %>%
	as.data.frame() %>%
	select(trap, quadratID) -> dat 

dat[duplicated(dat$trap),]
dat <- dat[!duplicated(dat$trap),]

#write.csv(dat, "../data/clean/trapQuadrats.csv", row.names=FALSE)

#join to trap data
trapDat <- read.csv("../data/clean/proportionAbscisedPerTrap.csv") 

trapDat$trap <- formatC(trapDat$trap, width = 3, format = "d", flag = "0")
trapDat$trap <- paste("trap", trapDat$trap, sep="_")

trapData <- left_join(trapDat, dat, by= "trap")

write.csv(trapData, "../data/clean/trapData.csv", row.names=FALSE)