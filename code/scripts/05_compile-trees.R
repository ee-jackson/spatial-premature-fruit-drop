#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-trees.R
## Desc: create tree_data
## Date created: 2020-12-08

library("tidyverse")
library("here")

# most recent census is in an excel sheet
bci10 <-
  readxl::read_excel(here::here("data", "raw","bci-50ha", "bci2022_2023.xlsx"),
                     na = c("", "NULL")) %>%
  rename_with(tolower) %>%
  mutate(census = 10,
         date = as.Date(exactdate)) %>%
  select(-c(plotname, family, genus, speciesname, subspecies, speciesid,
            quadratname, qy, qx, subspecies,
            listoftsm, highhom, exactdate, dbhid)) %>%
  rename(sp6 = mnemonic) %>%
  mutate(across(c(sp6, status, stemid, stemtag, tag, treeid), as.factor)) %>%
  select(sort(tidyselect::peek_vars()))

# older censuses are txt files - one for each census
file_names <- as.list(dir(path = here::here("data", "raw", "bci-50ha"),
                         pattern = "*bci_all_stems.txt", full.names = TRUE))

# census 9 was abandoned due to covid
file_names[9] <- NULL

bci_data_list <-
  lapply(file_names, read_tsv)

names(bci_data_list) <-
  lapply(file_names, basename)

species_list <-
  read.csv(here::here("data", "clean", "species_list.csv"))

R50 <-
  read.csv(here::here("data", "clean","R50.csv"))

tidy_bci <- function(bci_data, sp_data) {
  bci_data %>%
    rename_with(tolower) %>%
    select(-c(latin, subspecies, quadrat, stem,
              codes, no.)) %>%
    rename(sp6 = mnemonic) %>%
    filter(sp6 %in% sp_data$sp6) %>%
    mutate(across(c(sp6, status, stemid, stemtag, tag, treeid), as.factor)) %>%
    select(sort(tidyselect::peek_vars()))
}

bci_data_list <-
  lapply(bci_data_list, tidy_bci, sp_data = species_list)


# select the largest stem for multi stemmed individuals
# status AD = when a tree was noted as dead in one census but was found alive in a later census
# status AR = when a stem was broken in one census but alive in a later census
get_trees <- function(bci_data, sp_data) {
  bci_data %>%
    filter(sp6 %in% sp_data$sp6) %>%
    filter(status == "A" | status == "AD"| status == "AR" |
             status == "alive") %>%
    group_by(treeid) %>%
    slice_max(order_by = dbh, with_ties = TRUE) %>%
    slice_max(order_by = hom, with_ties = FALSE) %>%
    ungroup()
}

bci10_onestem <-
  get_trees(bci_data = bci10, sp_data = species_list)

bci_data_onestem_list <-
  lapply(bci_data_list, get_trees, sp_data = species_list)

# expand the bci survey data so we can match it to nearest years of the trap data
# we are repeating each row of the bci datasets 5 times and assigning 2 years either side

bci3 <-
  bci_data_onestem_list$`03_bci_all_stems.txt`[rep(seq_len(nrow(
    bci_data_onestem_list$`03_bci_all_stems.txt`
  )), each = 3), ]
bci3$year <- rep(c("1990", "1991", "1992"),
                 nrow(bci3) / 3)

bci4 <-
  bci_data_onestem_list$`04_bci_all_stems.txt`[rep(seq_len(nrow(
    bci_data_onestem_list$`04_bci_all_stems.txt`
  )), each = 5), ]
bci4$year <- rep(c("1993", "1994", "1995", "1996", "1997"),
                 nrow(bci4) / 5)

bci5 <-
  bci_data_onestem_list$`05_bci_all_stems.txt`[rep(seq_len(nrow(
    bci_data_onestem_list$`05_bci_all_stems.txt`
  )), each = 5), ]
bci5$year <- rep(c("2002", "2001", "2000", "1999", "1998"),
                 nrow(bci5) / 5)

bci6 <-
  bci_data_onestem_list$`06_bci_all_stems.txt`[rep(seq_len(nrow(
    bci_data_onestem_list$`06_bci_all_stems.txt`
  )), each = 5), ]
bci6$year <- rep(c("2003", "2004", "2005", "2006", "2007"),
                 nrow(bci6) / 5)

bci7 <-
  bci_data_onestem_list$`07_bci_all_stems.txt`[rep(seq_len(nrow(
    bci_data_onestem_list$`07_bci_all_stems.txt`
  )), each = 5), ]
bci7$year <- rep(c("2008", "2009", "2010", "2011", "2012"),
                 nrow(bci7) / 5)

bci8 <-
  bci_data_onestem_list$`08_bci_all_stems.txt`[rep(seq_len(nrow(
    bci_data_onestem_list$`08_bci_all_stems.txt`
  )), each = 6), ]
bci8$year <- rep(c("2013", "2014", "2015", "2016", "2017", "2018"),
                 nrow(bci8) / 6)

bci10 <-
  bci10_onestem[rep(seq_len(nrow(bci10_onestem)), each = 6), ]
bci10$year <- rep(c("2019", "2020", "2021", "2022", "2023", "2024"),
                          nrow(bci10) / 6)

# bind all the bci datasets together
rm(bci_data_list)
rm(bci_data_onestem_list)
rm(bci10_onestem)
bind_rows(mget(ls(pattern="^bci*")), .id='df') -> bci_bind

# join up data and format
bci_bind %>%
  left_join(species_list, by = "sp6") %>%
  left_join(R50, by = "sp6") %>%
  mutate(treeid =  formatC(
    treeid,
    width = 6,
    format = "d",
    flag = "0"
  )) %>%
  mutate(treeid = paste("tree", treeid, sep = "_")) %>%
  rename(tree = treeid,
         x = px, y = py,
         dbh_mm = dbh) %>%
  mutate(dbh_m = dbh_mm / 1000) %>% # calculate basal area
  mutate(radius_m = dbh_m / 2) %>%
  mutate(basal_area_m2 = pi * radius_m^2) %>%
  select(sp4, sp6, genus, species, family,
         capsules, dioecious, animal_disp, dsp_wind, r50,
         tree, x, y, dbh_mm, basal_area_m2, census, year) -> tree_data

saveRDS(tree_data,
        here::here("data", "clean", "tree_data.rds"))
