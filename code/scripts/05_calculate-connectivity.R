#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-connectivity.R
## Desc: Script used to calculate connectivity indicies for each trap in every year
## Date created: 2020-09-30

library("tidyverse") # v 1.3.1
library("here") # v 1.0.1
library("rdist") # v 0.0.5

# Load data ---------------------------
trap_data <- readRDS(here::here("data", "clean", "trap_data.rds"))

tree_data <- readRDS(here::here("data", "clean", "tree_data.rds"))

# Select species ---------------------------
# trap_data %>%
# 	group_by(sp4, year) %>%
# 	mutate(n_traps = length(unique(trap))) %>%
# 	arrange(desc(n_traps)) %>%
# 	ungroup() %>%
# 	group_by(sp4) %>%
# 	mutate(median_traps = median(n_traps, na.rm = TRUE)) %>%
# 	select(sp4, median_traps) %>%
# 	distinct() %>%
# 	as.data.frame() -> trap_summary
#
#
# tree_data %>%
# 	group_by(sp4, survey_year) %>%
# 	mutate(n_trees = length(unique(tree))) %>%
# 	arrange(desc(n_trees)) %>%
# 	ungroup() %>%
# 	group_by(sp4) %>%
# 	mutate(median_trees = median(n_trees, na.rm = TRUE)) %>%
# 	select(median_trees, sp4) %>%
# 	distinct() %>%
# 	as.data.frame() -> tree_summary
#
# # inner join  = sp only included if occur in both lists
# summary <- inner_join(trap_summary, tree_summary, by = "sp4")
#
# summary %>%
#   filter(median_trees >15 & median_traps >15) %>%
# 	pull(sp4) -> sp_list

# Calculate euclidean distances ---------------------------

calculate_dist <- function (species, yr) {
	# subset tree data
	bd <- dplyr::filter(tree_data,
	                    (tree_data[,"sp4"] == species) & (tree_data[,"year"] == yr) )

	# subset trap data
	td <- dplyr::filter(trap_data,
	                    (trap_data[,"sp4"] == species) & (trap_data[,"year"] == yr) )

	# create distance matrix using x and y co-ordinates
	dists <- cdist(bd[,c("x", "y")], td[,c("x", "y")], metric = "euclidean")
	dists_df <- as.data.frame(dists)

	# label col names with trap IDs
	colnames(dists_df) <- unlist(td$trap)

	# bind with the tree data
	cbind(bd, dists_df)
}

# create lists of sp and year to pass to function
trap_data %>%
  distinct(year) %>%
  inner_join(distinct(tree_data, year)) %>%
  pull(year) -> year_list

trap_data %>%
  distinct(sp4) %>%
  inner_join(distinct(tree_data, sp4)) %>%
  pull(sp4) -> sp_list

# apply function to all pairwise combinations of year and species
# will return list of dfs
all_dists <- outer(sp_list, year_list, FUN = Vectorize(calculate_dist))

bci_dists <- dplyr::bind_rows(all_dists)

# Hanski's Connectivity index ---------------------------

# function to calculate CI
calculate_CI <- function (trap) {
	drop_na(bci_dists, trap) %>%
	group_by(tree, year) %>%
	mutate(a = dbh * exp( (-0.02) * eval(parse(text = trap)) ) ) %>%
	ungroup() %>%
	dplyr::select(year, sp4, trap, a) %>%
	group_by(year, sp4) %>%
	summarise(trap = paste(trap), connectivity = sum(a), .groups = "drop")
}

trap_data %>%
  distinct(trap) %>%
  pull(trap) -> trap_list

CI_data <- parallel::mclapply(trap_list, calculate_CI, mc.cores = 3)

CI_data_b <- bind_rows(CI_data)

# Merge and save dataset ---------------------------

CI_data_b  %>%
  left_join(trap_data, by = c("trap", "year", "sp4")) -> trap_connect

saveRDS(trap_connect,
	file = here::here("data", "clean", "trap_connect.rds"))
