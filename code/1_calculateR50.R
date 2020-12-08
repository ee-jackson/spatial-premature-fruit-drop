#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculateR50.R
## Desc: Calculate the reproductive size threshold (50%) for each sp. R 50 =  1/2 dmax 
## Date: August 2020

rm(list = ls())

require("tidyverse")

# cd "/Users/eleanorjackson/OneDrive - University of Reading/SpatialPatterns/data/bci.tree"
file_names = as.list(dir(pattern = "bci.tree*"))
lapply(file_names, load, environment())

#cd "/Users/eleanorjackson/OneDrive - University of Reading/SpatialPatterns/code"

bci.tree <- bind_rows(bci.tree3,bci.tree4,bci.tree5,bci.tree6,bci.tree7,bci.tree8, .id = "df")

# Dmax as mean of 6 the individuals with the highest dbh
bci.tree %>%
	group_by(treeID) %>%
	slice_max(dbh, n = 1) %>% # use the largest dbh for each treeID
	ungroup() %>%
	group_by(sp) %>%
	slice_max(dbh, n = 6) %>% 
	summarise(dbh=mean(dbh, na.rm = TRUE)) %>%
	mutate(R50 = dbh/2) %>%
	rename(dmax="dbh") -> tree.max

write.csv(tree.max, "../data/clean/R50.csv", row.names=FALSE)