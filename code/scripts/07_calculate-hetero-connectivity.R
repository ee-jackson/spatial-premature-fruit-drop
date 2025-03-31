#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-connectivity.R
## Desc: Calculate connectivity indicies for each trap in every year
## Date created: 2023-08-02

library("dplyr") # v 1.3.1
library("rdist") # v 0.0.5
library("tidyr")

# Load data ---------------------------
trap_data <- readRDS("/home/users/ft840275/spatial_patterns/data/clean/trap_data.rds")

tree_data <- readRDS("/home/users/ft840275/spatial_patterns/data/clean/tree_data.rds")

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
calculate_hetero_CI <- function (trap_id, sp_id) {
  drop_na(bci_dists, trap_id) %>%
    group_by(tree, year) %>%
    mutate(a = exp( (-1/20) * eval(parse(text = trap_id)) ) * dbh^(0.5) ) %>%
    ungroup() %>%
    dplyr::filter(sp4 != sp_id) %>%
    dplyr::select(year, all_of(trap_id), a) %>%
    group_by(year) %>%
    summarise(trap = paste(trap_id), sp4 = paste(sp_id),
              connectivity = sum(a, na.rm = TRUE), .groups = "drop")
}

trap_data %>%
  filter(sp4 %in% sp_list) %>%
  select(trap, sp4) %>%
  distinct() %>%
  rename(trap_id = trap, sp_id = sp4) -> trap_sp_df

CI_data_h <-
  purrr::map2(.f = calculate_hetero_CI,
              .x= trap_sp_df$trap_id,
              .y = trap_sp_df$sp_id)

CI_data_h_b <- bind_rows(CI_data_h)

# Merge and save dataset ---------------------------

CI_data_h_b  %>%
  inner_join(trap_data, by = c("trap", "year", "sp4")) -> trap_connect

saveRDS(trap_connect,
        file = "/home/users/ft840275/spatial_patterns/data/clean/hetero_trap_connect_20m.rds")
