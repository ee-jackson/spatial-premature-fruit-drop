#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-abundance
## Desc: create tree_data
## Date created: 2026-02-20


# Packages ----------------------------------------------------------------

library("tidyverse")
library("here")


# Get data ----------------------------------------------------------------

tree_data <-
  readRDS(here::here("data", "clean", "tree_data.rds")) %>%
  filter(year %in% c("1990", "1995", "2000", "2005", "2010",
                     "2010", "2015", "2022")) # not replicate years

trap_data <-
  readRDS(here::here("data", "clean", "trap_data.rds"))

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

tree_data <-
  tree_data %>%
  filter(sp4 %in% shared_sp$sp4)


# Calculate abundance -----------------------------------------------------

abun <-
  tree_data %>%
  group_by(year, sp4, genus, species) %>%
  summarise(abundance = n_distinct(tree, na.rm = TRUE)) %>%
  group_by(sp4, genus, species) %>%
  summarise(median_abundance = median(abundance, na.rm = TRUE),
            mean_abundance = mean(abundance, na.rm = TRUE),
            max_abundance = max(abundance, na.rm = TRUE))

write_csv(abun,
          here::here("data", "clean", "species_abundance.csv"))
