#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculateR50.R
## Desc: Calculate the reproductive size threshold (50%) for each sp. R 50 =  1/2 dmax 
## Date: August 2020

rm(list = ls())

require(tidyverse)

# cd "/Users/eleanorjackson/OneDrive - University of Reading/SpatialPatterns/data/bci.tree"
file_names = as.list(dir(pattern = "bci.tree*"))
lapply(file_names, load, environment())

#cd "/Users/eleanorjackson/OneDrive - University of Reading/SpatialPatterns/code"

bci.tree <- bind_rows(bci.tree1,bci.tree2,bci.tree3,bci.tree4,bci.tree5,bci.tree6,bci.tree7,bci.tree8, .id = "id")

# Dmax as mean of 6 the individuals with the highest dbh
bci.tree %>%
	group_by(sp) %>%
	arrange(desc(dbh)) %>% 
	slice(6) %>% 
	summarise(dbh=mean(dbh)) %>%
	mutate(R50 = dbh/2) %>%
	ungroup() %>%
	rename(dmax="dbh") -> tree.max

# Dmax = the dbh of the individual with the largest value
bci.tree %>%
	group_by(sp) %>%
	slice(which.max(dbh)) %>%
	mutate(R50 = dbh/2) %>%
	ungroup() %>%
	select("sp", "dbh", "R50") %>%
	rename(dmax="dbh") -> tree.max1

write.csv(tree.max, "../output/tables/R50.csv")
