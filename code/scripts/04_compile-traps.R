#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-traps.R
## Desc: create dataset of proportion of seeds abscised per trap per sp per year
## Date created: 2020-08-05

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
species_list <-
  read_csv(here::here("data", "clean", "species_list.csv"))

# trap locations
trap_locs <-
  read_csv(here::here("data", "clean", "trap_locations.csv"))

# phenological year starts
pheno_years <-
  read_csv(here::here("data", "clean", "phenological_years.csv"))


# clean seed rain data ----------------------------------------------------

# tidy and combine with species data to drop species that we can't use
seed_rain %>%
  mutate(date = as.Date(date, "%Y-%m-%d"),
         calendar_year = format(as.Date(date), "%Y"),
         day_of_year = yday(date) - 1) %>%
  # before Nov 1989 immature fruits also contained aborted fruits
  filter(!calendar_year %in% c("1987", "1988", "1989")) %>%
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

# add pheno year data
seed_pheno <-
  seed_dat %>%
  left_join(pheno_years, by = "sp4") %>%
  mutate(calendar_year = as.numeric(calendar_year)) %>%
  # if the census date is after the pheno start then use calendar year
  # if census date is before the pheno start, use the previous calendar year
  mutate(
    pheno_year = if_else(
      day_of_year >= pheno_year_start,
      calendar_year,
      calendar_year - 1
    )
  ) %>%
  mutate(
    pheno_year = if_else(
      is.na(pheno_year),
      calendar_year,
      pheno_year
    )
  ) %>%
  mutate(pheno_year = as.character(pheno_year)) %>%
  # pheno year 2024 is incomplete
  # will only be from pheno start date to 2024-12-30
  filter(pheno_year != "2024")

# sum weekly counts of parts into yearly counts per trap per species
seed_sums <-
  seed_pheno %>%
  filter(part == 1 | part == 2 | part == 3 | part == 4 | part == 5) %>%
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
    part == 5 ~ quantity * seeds_per_fruit))

# count of parts
sum_dat <-
  seed_pheno %>%
  filter(part == 1 | part == 2 | part == 3 | part == 4 | part == 5) %>%
  group_by(sp4, pheno_year, trap) %>%
  summarise(sum_parts = sum(quantity, na.rm = TRUE),
            .groups = "drop")


# calculate viable seeds --------------------------------------------------

# estimate number of mature fruits based on fruits + capsules,
seed_sums %>%
  filter(part == 1 | part == 3) %>%
  group_by(sp4, pheno_year, trap, capsules) %>%
  summarise(capsules_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> sum_fruits_caps

# estimate number of mature fruits based on fruits + fragments,
seed_sums %>%
  filter(part == 1 | part == 4) %>%
  group_by(sp4, pheno_year, trap, capsules) %>%
  summarise(frags_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> sum_fruits_frags

# and based on fruits + seeds
seed_sums %>%
  filter(part == 1 | part == 2) %>%
  group_by(sp4, pheno_year, trap, capsules) %>%
  summarise(fruits_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> sum_fruits_seeds

#combine
full_join(sum_fruits_caps, sum_fruits_seeds)  %>%
  full_join(sum_fruits_frags) %>%
  # if we can't estimate fruits from capsules for a species, make NA
  mutate(capsules_seeds = case_when(
    capsules == FALSE | is.na(capsules) ~ NA,
    .default = capsules_seeds)) %>%
  filter(!if_all(c(capsules_seeds,
                   fruits_seeds,
                   frags_seeds), is.na)) %>%
  # use the largest estimate
  rowwise() %>%
  mutate(viable_seeds = max(
    capsules_seeds,
    fruits_seeds,
    frags_seeds,
    na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-capsules_seeds, -fruits_seeds, -frags_seeds) -> abs_dat_viable


# calculate abscised seeds ------------------------------------------------

# abscised seeds,  part 5 = immature fruit
seed_sums %>%
  filter(part == 5) %>%
  group_by(sp4, pheno_year, trap, capsules) %>%
  summarise(abscised_seeds = sum(seed_equiv, na.rm = TRUE),
            .groups = "drop") -> abs_dat_abscised


# calculate proportion abscised -------------------------------------------

# join abscised and viable seeds data, where there is no match will = NA, change to zero
full_join(abs_dat_abscised, abs_dat_viable) %>%
  mutate(viable_seeds = replace_na(viable_seeds, 0),
         abscised_seeds = replace_na(abscised_seeds, 0)) -> abs_dat_abscised_viable

# calculate proportion abscised
# use non-integers for calculating proportion_abscised
# but round values for estimates of seeds
abs_dat_abscised_viable %>%
  rowwise() %>%
  mutate(proportion_abscised =
           abscised_seeds / sum(abscised_seeds, viable_seeds, na.rm = TRUE)) %>%
  mutate_at(c("abscised_seeds", "viable_seeds"), ~ round(.)) %>%
  mutate(total_seeds = abscised_seeds + viable_seeds) %>%
  ungroup() %>%
  filter(total_seeds != 0) -> prop_dat

prop_dat %>%
  left_join(sum_dat, by = c("sp4", "pheno_year", "trap")) %>%
  left_join(trap_locs, by = "trap") %>%
  rename(year = pheno_year) %>%
  # don't include traps < 20m from the edge of the plot
  filter(x < 980 & x > 20) %>%
  filter(y < 480 & y > 20) -> trap_dat

saveRDS(trap_dat,
          here::here("data", "clean", "trap_data.rds"))
