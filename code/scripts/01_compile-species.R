#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-traps.R
## Desc: Compile list of species in the seed rain dataset and add data on reproduction/dispersal modes
## Date created: 2021-04-20

library("tidyverse")
library("WorldFlora")
library("readr")
library("here")
library("foreign")

# trap data
seedRain <-
  foreign::read.dbf(
  here::here("data", "raw", "BCI_TRAP200_20250224.DBF")) %>%
  mutate(MASS = na_if(MASS, -9)) %>%
  rename_with(tolower)

# species ID codes + nomenclature
read.csv(here::here("data", "raw", "20210412_nomenclature.csv")) %>%
  rename_with(tolower) %>%
  select(sp4, sp6, genus, species) -> nomenclature1

read.csv(here::here("data", "raw", "20120305_nomenclature.csv")) %>%
  rename_with(tolower) %>%
  select(sp4, sp6, genus, species) %>%
  drop_na(sp4) %>%
  filter(!sp4 %in% nomenclature1$sp4) -> nomenclature2

nomenclature <-
  bind_rows(nomenclature1, nomenclature2)

# capsules per fruit as defined by expert
read.csv(here::here("data",
                    "raw",
                    "Metadata_Capsulas_Part3_OsvaldoCalderon.csv")) -> capsulas

# species dispersal modes
read.csv(here::here("data",
                    "raw",
                    "species_20070228_DispersalModes.csv")) %>%
  select(sp, code6, starts_with("dsp")) %>%
  mutate(sp = na_if(sp, "")) %>%
  mutate(code6 = na_if(code6, "")) -> dispersal

# seed traits
seed_trait <-
  read_tsv(
    here::here(
      "data",
      "raw",
      "seed-masses",
      "Fruit_Seed_masses_20210111_BanisteriopsisCorrectedToBactrismajor.txt"
    )
  ) %>%
  rename_with(tolower) %>%
  mutate(id = paste(sp4, ind.unique, fruit.unique, sep = "_")) %>%
  select(sp4, id, n_seedfull, n_capsules, lifeform) %>%
  distinct() %>%
  group_by(sp4, lifeform) %>%
  summarise(seeds_per_fruit = mean(n_seedfull, na.rm = TRUE),
            capsules_per_fruit = mean(n_capsules, na.rm = TRUE))


# Get list of unique species in seed rain data ----------------------------

# some cleaning
seedRain %>%
  mutate(species = str_squish(species)) %>%
  filter(species != "dust", species != "branch",
         species != "petiole", species != "pedicel-fruit",
         species != "pedicel-flower", species != "",
         species != "NADA", species != "AUXX",
         species != "UNK?") %>%
  rename(sp4 = species) %>%
  distinct(sp4) -> seedRain_sp

# join with nomenclature and seed trait data, only woody stems
seedRain_sp %>%
  left_join(nomenclature, by = "sp4") %>%
  left_join(seed_trait, by = "sp4") %>%
  filter(
    lifeform == "S" |
      lifeform == "U" |
      lifeform == "M" | lifeform == "T"
  ) -> species


# Update species names ----------------------------------------------------

# update species names according to The World Flora Online (http://www.worldfloraonline.org)
# data v.2024.12 valid at 2025-04-23 https://doi.org/10.5281/zenodo.14538251

species %>%
  select(sp4, genus, species) %>%
  filter(!is.na(species)) %>%
  filter(!is.na(genus)) %>%
  mutate(genus_species = paste0(genus, sep = " ", species)) %>%
  select(-c(genus, species)) -> orig_sp

wfo <- WFO.match(spec.data = orig_sp,
                 WFO.file = here::here("data", "raw", "classification.csv"),
                 no.dates = TRUE, spec.name = "genus_species")

# finds one unique matching name for each submitted name
WFO.one(wfo, priority = "Accepted") %>%
  select(sp4, family, genus,
         specificEpithet, scientificNameAuthorship)  -> wfo_sp

species %>%
  select(-c(genus, species)) %>%
  left_join(wfo_sp, by = "sp4") %>%
  rename(species = specificEpithet,
         name_authorship = scientificNameAuthorship) %>%
  mutate(across(c(family, genus, species),
            na_if, "")) -> species_names


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
    capsules_per_fruit_oc =
      case_when(
        is.na(Capsulas_por_fruta_disecciones) ~ Capsulas_por_fruta_conocimiento_expert,
        Capsulas_por_fruta_disecciones == NaN ~ Capsulas_por_fruta_conocimiento_expert,
        Capsulas_por_fruta_disecciones == 0 ~ Capsulas_por_fruta_conocimiento_expert,
        TRUE ~ Capsulas_por_fruta_disecciones
      )
  ) %>%
  mutate(capsules_per_fruit_oc = ifelse(capsules == FALSE, NA, capsules_per_fruit_oc)) %>%
  mutate(capsules_per_fruit_oc = gsub(" a ", ",", capsules_per_fruit_oc)) %>%
  separate(
    capsules_per_fruit_oc,
    into = c("min", "max"),
    sep = ",",
    remove = FALSE
  ) %>%
  mutate_at(c("capsules_per_fruit_oc", "min", "max"), readr::parse_number) %>%
  rowwise() %>%
  mutate(median = median(c(min, max))) %>%
  ungroup() %>%
  mutate(capsules_per_fruit_oc = ifelse(is.na(median), capsules_per_fruit_oc, median)) %>%
  mutate(capsules_per_fruit_oc = round(capsules_per_fruit_oc, digits = 2)) %>%
  select(sp4, capsules, capsules_per_fruit_oc) -> capsulas_clean

species_names %>%
  left_join(capsulas_clean, by = "sp4") %>%
  mutate(capsules_per_fruit = case_when(
    is.na(capsules_per_fruit) ~ capsules_per_fruit_oc,
    .default = capsules_per_fruit
  )) %>%
  select(- capsules_per_fruit_oc) %>%
  mutate(capsules = case_when(
    is.na(capsules) & capsules_per_fruit >  0 ~ TRUE,
    .default = capsules)
    ) -> species_capsulas


# Clean dispersal mode data -----------------------------------------------

# fill any missing 4 digit IDs in dispersal data
key <-
  species_names %>%
  select(sp4, sp6) %>%
  drop_na()

missing_sp4 <-
  dispersal %>%
  filter(is.na(sp)) %>%
  left_join(key, by = c("code6" = "sp6")) %>%
  mutate(sp = coalesce(sp4, sp)) %>%
  select(-sp4) %>%
  drop_na(sp)

dispersal <-
  bind_rows(dispersal, missing_sp4) %>%
  drop_na(sp) %>%
  distinct() %>%
  filter(code6 != "CHR1EC") %>% # duplicated sp4
  select(-code6)

# create column indicating animal dispersal
species_capsulas %>%
  left_join(y = dispersal, by = c("sp4" = "sp")) %>%
  mutate(
    animal_disp = case_when(
        dsp_bat == TRUE |
        dsp_bird == TRUE |
        dsp_bbird == TRUE |
        dsp_mam == TRUE ~ TRUE,
      TRUE ~ FALSE
    )
  ) -> dispersal_clean


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

dispersal_clean %>%
  mutate_at(c("sp6", "sp4"), tolower) %>%
  mutate(dioecious = case_when(sp6 %in% dioecious_sp ~ TRUE,
                               TRUE ~ FALSE)
         ) -> species_capsulas_dispersal_dioecious


# Write clean data --------------------------------------------------------

# only keep animal dispersed species if we can estimate fruits from capsules

species_capsulas_dispersal_dioecious %>%
  filter(
    case_when(
      animal_disp == TRUE & (capsules == FALSE | is.na(capsules)) ~ F,
      .default = T
    ) ) %>%
  filter(!is.na(seeds_per_fruit)) %>%
  write.csv(here::here("data", "clean", "species_list.csv"),
          row.names = FALSE)
