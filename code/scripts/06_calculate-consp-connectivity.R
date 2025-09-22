#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-connectivity.R
## Desc: Calculate connectivity to reproductive conspecifics
## Date created: 2020-09-30

library("tidyverse")
library("here")
library("rdist")


# Load data ---------------------------

trap_data <- readRDS(here::here("data", "clean", "trap_data.rds"))

tree_data <-
  readRDS(here::here("data", "clean", "tree_data.rds")) %>%
  filter(dbh_mm >= r50)  %>% # only reproductive-sized
  select(sp4, year, tree, x, y, basal_area_m2)

monoecious_species <-
  read.csv(here::here("data", "clean", "species_list.csv")) %>%
  filter(dioecious != TRUE)


# Filter species ----------------------------------------------------------

# only keep species which appear in both datasets
shared_sp <-
  trap_data %>%
  select(sp4) %>%
  distinct() %>%
  inner_join(
    y = tree_data %>%
      select(sp4) %>%
      distinct()
  ) #%>%
  #filter(sp4 %in% monoecious_species$sp4) # drops 21 dioecious species

trap_data <-
  trap_data %>%
  filter(sp4 %in% shared_sp$sp4)

tree_data <-
  tree_data %>%
  filter(sp4 %in% shared_sp$sp4)


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
keys <-
  expand(trap_data, sp4, year)

# apply function to all pairwise combinations of year and species
# will return list of dfs
all_dists <-
  purrr::map2(
    .x = keys$sp4,
    .y = keys$year,
    .f = calculate_dist,
    .progress = TRUE
  )

bci_dists <- dplyr::bind_rows(all_dists)

# Hanski's Connectivity index ---------------------------

# function to calculate CI
calculate_CI <- function (trap) {
	drop_na(bci_dists, trap) %>%
	group_by(tree, year) %>%
	mutate(a = exp( (-1/20) * eval(parse(text = trap)) ) * basal_area_m2 ) %>%
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
	file = here::here("data", "clean", "trap_connect_repro_consp_20m_dioecious.rds"))
