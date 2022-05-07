#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-traps.R
## Desc: create dataset of proportion of seeds abscised per trap per sp per year
## Date created: 2020-08-05

library("tidyverse") # v 1.3.1
library("here") # v 1.0.1

# seed rain data
read.table(
  here::here("data", "raw", "BCI_TRAP200_20190215_spcorrected.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
) -> seed_rain

# seed traits from Joe
read.csv(
  here::here("data", "raw", "20120227_seedsMassForTraits.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  select(LIFEFORM, N_SEEDFULL, SP4) -> seed_trait

read.csv(here::here("data", "clean", "species_list.csv")) -> species_list

read.csv(here::here("data", "clean", "trap_locations.csv")) -> trap_locs

# clean seed rain data ----------------------------------------------------

seed_rain %>%
  mutate(fecha = as.character(fecha)) %>%
  mutate(fecha = as.Date(fecha, "%Y-%m-%d"),
         year = format(as.Date(fecha), "%Y")) %>%
  left_join(seed_trait, by = c("sp" = "SP4")) %>%

  # don't have the full data for these years
  filter(year != "1987" & year != "2019") %>%

  # keep only woody stems
  filter(
    LIFEFORM == "SHRUB" |
      LIFEFORM == "UNDERSTORY" |
      LIFEFORM == "MIDSTORY" | LIFEFORM == "TREE"
  ) %>%
  drop_na("N_SEEDFULL") %>%
  mutate(sp = tolower(sp),
         trap = formatC(
           trap,
           width = 3,
           format = "d",
           flag = "0"
         )) %>%
  mutate(trap = paste("trap", trap, sep = "_")) %>%
  rename(sp4 = sp, n_seedfull = N_SEEDFULL) -> seed_dat

# sum weekly counts of parts into yearly counts per trap per species
seed_dat %>%
  filter(part == 1 | part == 2 | part == 5 | part == 3) %>%
  group_by(part, sp4, year, trap, n_seedfull) %>%
  summarise(quantity_sum = sum(quantity, na.rm = TRUE),
            .groups = "drop") -> abs_dat

# combine with species data to drop species that we can't use
species_list %>%
  filter(dioecious == FALSE &
           ((animal_disp == FALSE | is.na(animal_disp)) | capsules == TRUE)) %>%
  select("sp4", "capsules", "capsules_per_fruit") %>%
  inner_join(abs_dat, by = "sp4") -> abs_dat_comb

# calculate viable seeds --------------------------------------------------

abs_dat_comb %>%
  subset(part == 1 | part == 2 | part == 3) %>%

  rowwise() %>%
  mutate(parts_avgseeds = case_when(

    # mature fruit
    part == 1 ~ quantity_sum * n_seedfull,

    # single diaspores
    part == 2 ~ quantity_sum,

    # capsules
    part == 3 ~ (quantity_sum / capsules_per_fruit) * n_seedfull)) %>%

  # estimate number of mature fruits based on fruits + capsules, and based on fruits + seeds
  group_by(sp4, year, trap, capsules) %>%
  summarise(capsules_seeds = ifelse(
    part == 1 | part == 3, sum(parts_avgseeds, na.rm = TRUE), 0),

    fruits_seeds = ifelse(
      part == 1 | part == 2, sum(parts_avgseeds, na.rm = TRUE), 0),
    .groups = "drop") %>%

  # if fruits + seeds estimate is larger, use that rather than fruits + capsules
  mutate(capsules_seeds =
           ifelse(fruits_seeds > capsules_seeds, fruits_seeds, capsules_seeds)) %>%

  # if we can't estimate fruits from capsules for a species, use fruits + seeds
  mutate(viable_seeds = case_when(
    capsules == FALSE | is.na(capsules) ~ fruits_seeds,
    capsules == TRUE ~ capsules_seeds)) %>%

  select(-capsules_seeds, -fruits_seeds) -> abs_dat_viable

# calculate abscised seeds ------------------------------------------------

# abscised seeds,  part 5 = immature fruit
abs_dat_comb %>%
  subset(part == 5) %>%
  rowwise() %>%
  mutate(abscised_seeds = quantity_sum*n_seedfull) %>%
  group_by(sp4, year, trap) %>%
  summarise(abscised_seeds = sum(abscised_seeds, na.rm = TRUE),
            .groups = "drop") -> abs_dat_abscised

# calculate proportion abscised -------------------------------------------

# join abscised and viable seeds data, where there is no match will = NA, change to zero
full_join(abs_dat_abscised, abs_dat_viable, by= c("sp4", "year", "trap")) %>%
  mutate_at(vars(viable_seeds, abscised_seeds), ~replace(., is.na(.), 0)) -> abs_dat_abscised_viable

# calculate proportion abscised
abs_dat_abscised_viable %>%
  rowwise() %>%
  mutate(total_seeds = sum(abscised_seeds, viable_seeds, na.rm = TRUE),
         proportion_abscised = abscised_seeds / total_seeds) %>%
  ungroup() -> prop_dat

abs_dat %>%
  group_by(sp4, year, trap) %>%
  summarise(sum_parts = sum(quantity_sum, na.rm = TRUE), .groups = "drop") -> sum_dat

prop_dat %>%
  left_join(sum_dat, by = c("sp4", "year", "trap")) %>%
  left_join(trap_locs, by = "trap") -> trap_dat

saveRDS(trap_dat,
          here::here("data", "clean", "trap_data.rds"))
