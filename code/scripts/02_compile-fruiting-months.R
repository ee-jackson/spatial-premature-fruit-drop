#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-traps.R
## Desc: create dataset of co-fruiting species
## Date created: 2025-04-29

library("tidyverse")
library("here")
library("foreign")

# seed rain data
foreign::read.dbf(
  here::here("data", "raw", "BCI_TRAP200_20250224.DBF")) %>%
  mutate(MASS = na_if(MASS, -9)) %>%
  rename_with(tolower) %>%
  select(-c(mass, seedeqs)) -> seed_rain

# species list
read.csv(here::here("data", "clean", "species_list.csv")) -> species_list

# clean seed rain data ----------------------------------------------------

# tidy and combine with species data to drop species that we can't use
seed_rain %>%
  mutate(date = as.Date(date, "%Y-%m-%d"),
         year = format(as.Date(date), "%Y"),
         month = format(as.Date(date), "%m")) %>%
  mutate(species = tolower(species),
         trap = formatC(
           trap,
           width = 3,
           format = "d",
           flag = "0"
         )) %>%
  mutate(trap = paste("trap", trap, sep = "_")) %>%
  rename(sp4 = species) %>%
  inner_join(species_list, by = "sp4") -> seed_dat


abs_dat_comb <-
  seed_dat %>%
  filter(part == 1 | part == 2 | part == 3 | part == 4 | part == 5 | part == 7) %>%
  rowwise() %>%
  mutate(seed_equiv = case_when(

    # mature fruit
    part == 1 ~ quantity * seeds_per_fruit,

    # single diaspores
    part == 2 ~ quantity,

    # capsules
    part == 3 ~ quantity * seeds_per_fruit,

    # fragments
    part == 4 ~ quantity * seeds_per_fruit,

    # immature fruits
    part == 5 ~ quantity * seeds_per_fruit,

    # insect damaged fruits - only recorded for select species
    part == 7 ~ quantity * seeds_per_fruit))


# calculate viable seeds --------------------------------------------------

# estimate number of fruits based on fruits + capsules,
abs_dat_comb %>%
  filter(part == 1 | part == 3| part == 5 | part == 7) %>%
  group_by(sp4, year, month, trap, capsules) %>%
  summarise(capsules_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> sum_caps_seeds

# estimate number of fruits based on fruits + fragments,
abs_dat_comb %>%
  filter(part == 1 | part == 4 | part == 5 | part == 7) %>%
  group_by(sp4, year, month, trap, capsules) %>%
  summarise(frags_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> sum_frag_seeds

# and based on fruits + seeds
abs_dat_comb %>%
  filter(part == 1 | part == 2| part == 5| part == 7) %>%
  group_by(sp4, year, month, trap, capsules) %>%
  summarise(fruits_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> sum_fruits_seeds

full_join(sum_caps_seeds, sum_fruits_seeds)  %>%
  full_join(sum_frag_seeds) %>%
  # if we can't estimate fruits from capsules for a species, make NA
  mutate(capsules_seeds = case_when(
    capsules == FALSE | is.na(capsules) ~ NA,
    .default = capsules_seeds)) %>%
  filter(!if_all(c(capsules_seeds,
                   fruits_seeds,
                   frags_seeds), is.na)) %>%
  # use the largest estimate
  rowwise() %>%
  mutate(n_seeds = max(
    capsules_seeds,
    fruits_seeds,
    frags_seeds,
    na.rm = TRUE)
    ) %>%
  ungroup() %>%
  select(-capsules_seeds, -fruits_seeds) -> abs_dat

# quick look
abs_dat %>%
  group_by(sp4, month) %>%
  summarise(median = mean(n_seeds, na.rm = TRUE)) %>%
  ggplot(aes(x = month, y = median)) +
  geom_col() +
  facet_wrap(~sp4, scales = "free_y")


# Get co-fruiting species -------------------------------------------------

# list of unique year-months for each species
fruit_months <-
  abs_dat %>%
  drop_na(n_seeds) %>%
  select(sp4, month, year) %>%
  distinct() %>%
  mutate(yr_month = paste(year, month, sep = "_"))

# list of species
sp4_list <-
  abs_dat %>%
  select(sp4) %>%
  distinct() %>%
  pull(sp4)

# function to get species with overlapping year-months for species i
get_cofruit_species <-
  function(species_id, fruit_months_data) {
  sp_i_month_yrs <-
    fruit_months %>%
    filter(sp4 == species_id) %>%
    select(yr_month)

  sp_i_co_sp <-
    fruit_months %>%
    filter(sp4 != species_id) %>%
    filter(yr_month %in% sp_i_month_yrs$yr_month) %>%
    select(sp4) %>%
    distinct() %>%
    rename(co_fruit_sp = sp4) %>%
    mutate(sp4 = species_id)

  return(sp_i_co_sp)
  }

# run function for every species
cofruit_data <-
  lapply(X = sp4_list,
       FUN = get_cofruit_species,
       fruit_months_data = fruit_months) %>%
  bind_rows() %>%
  group_by(sp4) %>%
  nest(co_fruit_sp = co_fruit_sp)

saveRDS(cofruit_data,
        here::here("data", "clean", "cofruit_data.rds"))
