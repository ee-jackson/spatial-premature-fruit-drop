#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-traps.R
## Desc: create dataset fruiting months for each species
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

    # damaged fruits
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

# check with a plot
abs_dat %>%
  group_by(sp4, month) %>%
  summarise(median = mean(n_seeds, na.rm = TRUE)) %>%
  ggplot(aes(x = month, y = median)) +
  geom_col() +
  facet_wrap(~sp4, scales = "free_y")

fruit_months <-
  abs_dat %>%
  select(sp4, month) %>%
  distinct()

write_csv(fruit_months,
        here::here("data", "clean", "fruit_months.csv"))
