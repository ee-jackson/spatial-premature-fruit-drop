#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-trees.R
## Desc: create tree_data
## Date created: 2020-12-08

library("tidyverse") # v 1.3.1
library("here") # v 1.0.1

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

file_names <- as.list(dir(path = here::here("data", "raw", "bci-50ha"),
                         pattern = "*bci_all_stems.txt", full.names = TRUE))

bci_data_list <-
  lapply(file_names, read_tsv)

names(bci_data_list) <-
  lapply(file_names, basename)

species_list <- read.csv(here::here("data", "clean", "species_list.csv"))

R50 <- read.csv(here::here("data", "clean","R50.csv"))

tidy_bci <- function(bci_data) {
  bci_data %>%
    rename_with(tolower) %>%
    select(-c(latin, subspecies, quadrat, stem,
              codes, no.)) %>%
    rename(sp6 = mnemonic) %>%
    mutate(across(c(sp6, status, stemid, stemtag, tag, treeid), as.factor)) %>%
    select(sort(tidyselect::peek_vars()))
}

bci_data_list <-
  lapply(bci_data_list, tidy_bci)

# expand the bci survey data so we can match it to nearest years of the trap data
# we are repeating each row of the bci datasets 5 times and assigning 2 years either side

bci3 <-
  bci_data_list$`03_bci_all_stems.txt`[rep(seq_len(nrow(bci_data_list$`03_bci_all_stems.txt`)), each = 5), ]
bci3$year <- rep(c("1990","1991","1992"), nrow(bci3)/3)

bci4 <-
  bci_data_list$`04_bci_all_stems.txt`[rep(seq_len(nrow(bci_data_list$`04_bci_all_stems.txt`)), each = 5), ]
bci4$year <- rep(c("1993","1994","1995","1996","1997"), nrow(bci4)/5)

bci5 <-
  bci_data_list$`05_bci_all_stems.txt`[rep(seq_len(nrow(bci_data_list$`05_bci_all_stems.txt`)), each = 5), ]
bci5$year <- rep(c("2002","2001","2000","1999","1998"), nrow(bci5)/5)

bci6 <-
  bci_data_list$`06_bci_all_stems.txt`[rep(seq_len(nrow(bci_data_list$`06_bci_all_stems.txt`)), each = 5), ]
bci6$year <- rep(c("2003","2004","2005","2006","2007"), nrow(bci6)/5)

bci7 <-
  bci_data_list$`07_bci_all_stems.txt`[rep(seq_len(nrow(bci_data_list$`07_bci_all_stems.txt`)), each = 5), ]
bci7$year <- rep(c("2008","2009","2010","2011","2012"), nrow(bci7)/5)

bci8 <-
  bci_data_list$`08_bci_all_stems.txt`[rep(seq_len(nrow(bci_data_list$`08_bci_all_stems.txt`)), each = 6), ]
bci8$year <- rep(c("2013","2014","2015","2016","2017","2018"), nrow(bci8)/6)

bci10 <-
  bci10[rep(seq_len(nrow(bci10)), each = 6), ]
bci10$year <- rep(c("2019","2020","2021","2022","2023","2024"), nrow(bci10)/6)

# bind all the bci datasets together
rm(bci_data_list)
bind_rows(mget(ls(pattern="^bci*")), .id='df') -> bci_bind

## TO DO multiple stems!!

# subset to the species in the trap data and clean up
species_list %>%
  select(sp4, sp6) %>%
  left_join(bci_bind, by = c("sp6" = "sp")) %>%
  left_join(R50, by = c("sp6" = "sp")) %>%
  # status AD = when a tree was noted as dead in one census but was found alive in a later census
  filter((status == "A" | status == "AD") & dbh >= r50) %>%
  mutate(treeID =  formatC(
    treeID,
    width = 6,
    format = "d",
    flag = "0"
  )) %>%
  mutate(treeID = paste("tree", treeID, sep = "_")) %>%
  rename(tree = treeID, x = gx, y = gy) %>%
  select(sp4, tree, x, y, dbh, year, survey_year) -> tree_data

saveRDS(tree_data,
        here::here("data", "clean", "tree_data.rds"))
