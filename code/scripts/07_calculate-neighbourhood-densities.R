#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-neighbourhood-densities.R
## Desc: Calculate neighbourhood densities using trap and tree datasets
## Can take > 20 mins to run locally - would advise running on HPC cluster
## Date created: 2026-03-02

library("tidyverse")
library("rdist", lib.loc = "~/local/rlibs")


# Load data ---------------------------------------------------------------

trap_data <-
  readRDS("data/clean/trap_data.rds")

tree_data <-
  readRDS("data/clean/tree_data.rds")

fruiting_data <-
  readRDS("data/clean/cofruit_data.rds")


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
  )

trap_data <-
  trap_data %>%
  filter(sp4 %in% shared_sp$sp4) %>%
  arrange(sp4, year, trap)

tree_data <-
  tree_data %>%
  filter(sp4 %in% shared_sp$sp4) %>%
  arrange(sp4, year)


# Partition trees into reproductive vs non repro --------------------------

tree_repro <- tree_data %>%
  filter(dbh_mm >= repro_dbh)

tree_nonrepro <- tree_data %>%
  filter(dbh_mm < repro_dbh)


# Helper function to compute neighbourhood density ------------------------

# Given:
#  bd: tibble of trees (rows: trees) with columns x,y,basal_area_m2
#  td: tibble of traps (rows: traps) with columns x,y,trap (trap order matters)
# returns numeric vector length nrow(td) of connectivity contributions

calculate_connectivity <- function(bd, td, alpha = 1/5) {
  # no traps -> empty numeric named vector
  if (nrow(td) == 0) return(set_names(numeric(0), character(0)))

  # no trees -> zeros for every trap
  if (nrow(bd) == 0) {
    z <- numeric(nrow(td))
    if (!is.null(td$trap)) names(z) <- td$trap
    return(z)
  }

  # compute pairwise distances: rows = trees, cols = traps
  D <- rdist::cdist(
    X = as.matrix(bd[, c("x", "y")]),
    Y = as.matrix(td[, c("x", "y")]),
    metric = "euclidean")

  contribs <- sweep(exp(- alpha * D), 1, bd$basal_area_m2, `*`)

  # sum by trap
  conn <- colSums(contribs, na.rm = TRUE)

  # remove big objects
  rm(D, contribs)

  # name the vector by trap to be explicit
  if (!is.null(td$trap)) names(conn) <- td$trap

  conn
}


# Function to compute all 3 densities per trap x sp x year ----------------

compute_all_three <-
  function(species, yr, fruiting_df, tree_repro_df, tree_nonrepro_df, trap_df, alpha = 1/5) {

  # traps for focal species-year
  td <- trap_df %>% filter(sp4 == !!species, year == !!yr)

  # if no traps for this species-year, return empty tibble
  if (nrow(td) == 0) {
    return(tibble(year = character(0), trap = character(0), sp4 = character(0),
                  conn_RH = double(0),
                  conn_RC = double(0),
                  conn_NRC = double(0)))
  }

  # get list of co-fruiting species
  co_sp <- fruiting_df %>%
    filter(sp4 == !!species) %>%
    tidyr::unnest(co_fruit_sp) %>%
    pull(co_fruit_sp) %>%
    as.character()

  # heterospecific reproductive trees that co-fruit
  # (same year, co-fruiting, not focal species)
  bd_RH <- tree_repro_df %>%
    filter(sp4 != !!species, sp4 %in% co_sp, year == !!yr)

  # conspecific reproductive trees
  # (same species, same year)
  bd_RC <- tree_repro_df %>%
    filter(sp4 == !!species, year == !!yr)

  # conspecific non-reproductive trees
  # (same species, same year)
  bd_NRC <- tree_nonrepro_df %>%
    filter(sp4 == !!species, year == !!yr)

  # calculate connectivity vectors
  conn_RH <- calculate_connectivity(bd_RH, td, alpha = alpha)
  conn_RC <- calculate_connectivity(bd_RC, td, alpha = alpha)
  conn_NRC <- calculate_connectivity(bd_NRC, td, alpha = alpha)

  # assemble tibble; preserve trap order from td
  tibble(
    year = !!yr,
    trap = td$trap,
    sp4 = !!species,
    alpha = !!alpha,
    conn_RH = as.numeric(conn_RH[td$trap]),
    conn_RC = as.numeric(conn_RC[td$trap]),
    conn_NRC = as.numeric(conn_NRC[td$trap])
  )
}


# Iterate over species x year combos --------------------------------------

keys <- trap_data %>%
  expand(sp4, year) %>%
  arrange(sp4, year)

alphas <- exp(seq(log(0.07), log(0.40), length.out = 6))

all_connectivities <- purrr::map2(
  .x = keys$sp4,
  .y = keys$year,
  .f = ~ compute_all_three(.x, .y,
                           fruiting_df = fruiting_data,
                           tree_repro_df = tree_repro,
                           tree_nonrepro_df = tree_nonrepro,
                           trap_df = trap_data,
                           alpha = alphas[[3]]),
  .progress = TRUE) %>%
  list_rbind()


# Alternative if testing different alpha values:

# keys <- trap_data %>%
#   expand(sp4, year, alpha = exp(seq(log(0.07), log(0.40), length.out = 6))) %>%
#   arrange(sp4, year)
#
# trap_connect_all3 <- purrr::pmap(
#   .l = list(species = keys$sp4,
#             yr = keys$year,
#             alpha = keys$alpha),
#   .f = compute_all_three,
#   fruiting_df = fruiting_data,
#   tree_repro_df = tree_repro,
#   tree_nonrepro_df = tree_nonrepro,
#   trap_df = trap_data) %>%
#   list_rbind() %>%
#   inner_join(trap_data, by = c("trap", "year", "sp4")) %>%
#   mutate_at(c("sp4", "trap", "quadrat", "year"), ~as.factor(.))


# Add seed rain and trap metadata -----------------------------------------

trap_connect_all3 <-
  all_connectivities %>%
  inner_join(trap_data, by = c("trap", "year", "sp4")) %>%
  mutate_at(c("sp4", "trap", "quadrat", "year"), ~as.factor(.))

saveRDS(trap_connect_all3,
        file = "data/clean/trap_connect.rds")

# saveRDS(trap_connect_all3,
#         file = "data/clean/connect_all_alpha.rds")
