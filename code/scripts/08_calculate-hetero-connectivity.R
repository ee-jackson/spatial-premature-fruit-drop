#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-connectivity.R
## Desc: Calculate connectivity to reproductive, co-fruiting heterospecifics
## Date created: 2023-08-02

library("dplyr")
library("rdist", lib.loc = "~/local/rlibs")
library("tidyr")


# Load data ---------------------------

trap_data <-
  readRDS("data/clean/trap_data.rds")

tree_data <-
  readRDS("data/clean/tree_data.rds") %>%
  filter(dbh_mm >= r50)  %>% # only reproductive-sized
  select(sp4, year, tree, x, y, basal_area_m2)

fruiting_data <-
  readRDS("data/clean/cofruit_data.rds")

monoecious_species <-
  read.csv("data/clean/species_list.csv") %>%
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
  ) %>%
  filter(sp4 %in% monoecious_species$sp4) # drops 21 dioecious species

trap_data <-
  trap_data %>%
  filter(sp4 %in% shared_sp$sp4)

tree_data <-
  tree_data %>%
  filter(sp4 %in% shared_sp$sp4)


# Calculate euclidean distances ---------------------------

# distance between each tree and every trap, in each year

calculate_dist <- function (species, yr, fruiting, tree_df, trap_df) {

  # get list of species which fruit at the same time as species i
  fd <- fruiting %>%
    dplyr::filter(sp4 == species) %>%
    tidyr::unnest(co_fruit_sp) %>%
    dplyr::pull(co_fruit_sp)

  # only keep heterospecific trees i.e. not species i
  # and only heterospecifics which fruit at the same time as species i
  bd <- tree_df %>%
    dplyr::filter(sp4 != species  &
                    sp4 %in% fd &
                    year == yr)

  # only keep traps with seeds of species i
  td <- trap_df %>%
    dplyr::filter(sp4 == species) %>%
    dplyr::filter(year == yr)

  # create distance matrix using x and y co-ordinates
  dists <- cdist(bd[,c("x", "y")], td[,c("x", "y")], metric = "euclidean")
  dists_df <- as.data.frame(dists)

  # label col names with trap IDs
  colnames(dists_df) <- unlist(td$trap)

  # bind with the tree data
  cbind(bd, dists_df) %>%
    mutate(sp_i = species) -> out

  return(out)
}

# get every possible comb of trap, yr and sp to pass to function
keys <-
  expand(trap_data, sp4, year)

# apply function
# will return list of dfs
all_dists <-
  purrr::map2(
    .x = keys$sp4,
    .y = keys$year,
    .f = calculate_dist,
    fruiting = fruiting_data,
    tree_df = tree_data,
    trap_df = trap_data,
    .progress = TRUE
  )

bci_dists <-
  dplyr::bind_rows(all_dists)


# Hanski's Connectivity index ---------------------------

# function to calculate CI
calculate_hetero_CI <- function (trap_id, sp_id) {
  bci_dists %>%
    drop_na(trap_id) %>%
    group_by(tree, year) %>%
    mutate(a = exp( (-1/20) * eval(parse(text = trap_id)) ) * basal_area_m2 ) %>%
    ungroup() %>%
    dplyr::filter(sp_i == sp_id) %>%
    dplyr::select(year, all_of(trap_id), a) %>%
    group_by(year) %>%
    summarise(trap = paste(trap_id), sp4 = paste(sp_id),
              connectivity = sum(a, na.rm = TRUE), .groups = "drop")
}

keys_traps <-
  expand(trap_data, sp4, trap) %>%
  rename(trap_id = trap, sp_id = sp4)

CI_data_h <-
  purrr::map2(.f = calculate_hetero_CI,
              .x = keys_traps$trap_id,
              .y = keys_traps$sp_id)

CI_data_h_b <- bind_rows(CI_data_h)


# Merge and save dataset ---------------------------

CI_data_h_b  %>%
  inner_join(trap_data, by = c("trap", "year", "sp4")) -> trap_connect

saveRDS(trap_connect,
        file = "data/clean/trap_connect_repro_hetero_20m.rds")
