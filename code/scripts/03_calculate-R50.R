#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-R50.R
## Desc: Calculate the reproductive size threshold (50%) for each sp. R 50 =  1/2 dmax
## Date created: 2020-08-10

library("tidyverse") # v 1.3.1
library("here") # v 1.0.1

# load all bci files
file_names <- as.list(dir(path = here::here("data", "raw", "bci.tree"),
                         pattern = "bci.tree*", full.names = TRUE))

lapply(file_names, load, environment())

bci.tree <- bind_rows(bci.tree3, bci.tree4, bci.tree5, bci.tree6,
                      bci.tree7, bci.tree8, .id = "df")

# Dmax as mean of 6 the individuals with the highest dbh
bci.tree %>%
	group_by(treeID) %>%
	slice_max(dbh, n = 1) %>% # use the largest dbh for each treeID
	ungroup() %>%
	group_by(sp) %>%
	slice_max(dbh, n = 6) %>%
	summarise(dbh = mean(dbh, na.rm = TRUE), .groups = "drop") %>%
	mutate(r50 = dbh/2) %>%
	rename(dmax = "dbh") -> tree_max

write.csv(tree_max, here::here("data", "clean", "R50.csv"),
          row.names = FALSE)
