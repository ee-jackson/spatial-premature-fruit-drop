#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-traps.R
## Desc: Compile list of species in the seed rain dataset and add data on reproduction/dispersal modes
## Date created: 2021-04-20

library("tidyverse") # v 1.3.1
library("WorldFlora") # v 1.1
library("readr") # v 2.0.2
library("here") # v 1.0.1

seedRain <- read.table(here::here("data", "raw",
                                  "BCI_TRAP200_20190215_spcorrected.txt"),
                       header = TRUE, stringsAsFactors = FALSE)

read.csv(here::here("data", "raw", "20120305_nomenclature.csv")) %>%
  rename_with(tolower) -> nomenclature

read.csv(here::here("data", "raw", "Metadata_Capsulas_Part3_OsvaldoCalderon.csv")) -> capsulas

read.csv(here::here("data", "raw", "species_20070228_DispersalModes.csv")) -> dispersal


# Get list of unique species in seed rain data ----------------------------

# some cleaning
seedRain %>%
  mutate(sp = str_squish(sp)) %>%
  filter(!(str_detect(sp, "\\d"))) %>%
  filter(sp != "dust",sp != "branch",sp != "petiole",sp != "pedicel-fruit",sp != "pedicel-flower",sp != "") %>%
  rename(sp4 = sp) %>%
  distinct(sp4) -> species

# join with nomenclature data
species %>%
  left_join(nomenclature, by = "sp4") %>%
  select(sp4, sp6, family, genus, species) -> species


# Update species names ----------------------------------------------------

# update species names according to The World Flora Online (http://www.worldfloraonline.org)
# data v.2021.01 valid at 2021-12-09

species %>%
  select(sp4, genus, species) %>%
  mutate(genus_species = paste0(genus, sep = " ", species)) %>%
  select(-c(genus, species)) -> orig_sp

wfo <- WFO.match(spec.data = orig_sp,
                 WFO.file = here::here("data", "raw", "classification.txt"),
                 no.dates = TRUE, spec.name = "genus_species")

# finds one unique matching name for each submitted name
WFO.one(wfo, priority = "Accepted") %>%
  select(sp4, family, genus,
         specificEpithet)  -> wfo_sp

left_join(species, wfo_sp, by = "sp4") %>%
  select(-c(family.x, genus.x, species)) %>%
  rename(
    species = specificEpithet,
    genus = genus.y,
    family = family.y
  ) %>%
  mutate(across(c("family", "genus", "species", ),
            na_if, "")) -> species_names

# GUA1 is listed as a synonym of Guarea guidonia (GUA2) in the World Flora Online database
# Here we adopt the name Guarea megantha for GUA1
species_names %>%
  mutate(species = ifelse(sp4 == "GUA1", "megantha", species)) -> species_names


# Clean capsule data ------------------------------------------------------

# Using dissection-based estimates of number of fruits per capsule for species
# for which they are available, and Osvaldo’s estimates for the rest

capsulas %>%
  mutate(
    capsules = ifelse(
      Capsulas_por_fruta_conocimiento_expert == "No" |
        Capsulas_por_fruta_conocimiento_expert == "Na - species not known to Osvaldo" |
        Capsulas_por_fruta_conocimiento_expert == "Comida por animal. Categoria 10." |
        Capsulas_por_fruta_conocimiento_expert == "----",
      FALSE,
      TRUE
    )
  ) %>%
  mutate(
    Capsulas_por_fruta_conocimiento_expert = as.character(Capsulas_por_fruta_conocimiento_expert),
    Capsulas_por_fruta_disecciones = as.character(Capsulas_por_fruta_disecciones)
  ) %>%
  mutate(
    capsules_per_fruit =
      case_when(
        is.na(Capsulas_por_fruta_disecciones) ~ Capsulas_por_fruta_conocimiento_expert,
        Capsulas_por_fruta_disecciones == NaN ~ Capsulas_por_fruta_conocimiento_expert,
        Capsulas_por_fruta_disecciones == 0 ~ Capsulas_por_fruta_conocimiento_expert,
        TRUE ~ Capsulas_por_fruta_disecciones
      )
  ) %>%
  mutate(capsules_per_fruit = ifelse(capsules == FALSE, NA, capsules_per_fruit)) %>%
  mutate(capsules_per_fruit = gsub(" a ", ",", capsules_per_fruit)) %>%
  separate(
    capsules_per_fruit,
    into = c("min", "max"),
    sep = ",",
    remove = FALSE
  ) %>%
  mutate_at(c("capsules_per_fruit", "min", "max"), readr::parse_number) %>%
  rowwise() %>%
  mutate(median = median(c(min, max))) %>%
  ungroup() %>%
  mutate(capsules_per_fruit = ifelse(is.na(median), capsules_per_fruit, median)) %>%
  mutate(capsules_per_fruit = round(capsules_per_fruit, digits = 2)) %>%
  select(sp4, capsules, capsules_per_fruit) -> capsulas_clean

species_names %>%
  left_join(capsulas_clean, by = "sp4") -> species_capsulas

# Clean dispersal mode data -----------------------------------------------

key <- species_names %>% select(sp4, sp6)

dispersal %>%
  left_join(key, by = c("code6" = "sp6")) %>%
  drop_na(sp4) %>%
  mutate(
    animal_disp = case_when(
      dsp_ant == TRUE |
        dsp_bat == TRUE |
        dsp_bird == TRUE |
        dsp_bbird == TRUE |
        dsp_mam == TRUE ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  select(sp4, animal_disp) -> dispersal_clean

species_capsulas %>%
  left_join(dispersal_clean, by = "sp4") -> species_capsulas_dispersal

# Add dioecious variable --------------------------------------------------

dioecious_sp <- c("ade1tr", "alchco", "alchla", "alibed", "amaico", "ast2gr",
  "brosal", "cecrin", "cecrob", "cha1te", "coccma", "coccco", "dio2ar",
  "drypst", "gar2in", "guapst", "guargu", "guarsp", "guargr", "hampap",
  "hyeral", "maclti", "maquco", "margno", "neeaam", "ocotwh", "ocotce",
  "ocotob", "perexa", "picrla", "poutre", "poutst", "poutfo", "protpa",
  "protte", "protco", "randar", "randfo", "simaam", "sipapa", "soroaf",
  "stylst", "tet2pa", "tratas", "tri2tu", "tri2pa", "tripcu", "tropra",
  "olmeas", "urerba", "virose", "virosu", "virosp", "xyl2ol", "zantbe",
  "zantp1", "zantpr", "zantse")

species_capsulas_dispersal %>%
  mutate_at(c("sp6", "sp4"), tolower) %>%
  mutate(dioecious = case_when(sp6 %in% dioecious_sp ~ TRUE,
                               TRUE ~ FALSE)
         ) -> species_capsulas_dispersal_dioecious


# Write clean data --------------------------------------------------------

write.csv(species_capsulas_dispersal_dioecious,
          here::here("data", "clean", "species_list.csv"),
          row.names = FALSE)
