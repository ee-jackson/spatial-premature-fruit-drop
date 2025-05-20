#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-R50.R
## Desc: Calculate the reproductive size threshold (50%) for each sp. R 50 =  1/2 dmax
## Date created: 2020-08-10

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

bci_data_list <-
  lapply(file_names, read_tsv)

names(bci_data_list) <-
  lapply(file_names, basename)

# make format consistent with bci10 file
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

# select the largest stem for multi stemmed individuals
get_trees <- function(bci_data) {
  bci_data %>%
    filter(status == "A" | status == "AD"| status == "AR" |
             status == "alive") %>%
    group_by(treeid) %>%
    slice_max(order_by = dbh, with_ties = TRUE) %>%
    slice_max(order_by = hom, with_ties = FALSE) %>%
    ungroup()
}

bci_data_onestem_list <-
  lapply(bci_data_list, get_trees)

bci_data_onestem_list$bci10 <-
  get_trees(bci_data = bci10)

# bind files together
bci_bind <-
  bind_rows(mget(ls(pattern="^bci*")), .id='df')

# Dmax as mean of 6 the individuals with the highest dbh
bci_bind %>%
	group_by(treeid) %>%
	slice_max(dbh, with_ties = FALSE) %>% # use the largest dbh for each treeID
	ungroup() %>%
	group_by(sp6) %>%
	slice_max(dbh, n = 6, with_ties = FALSE) %>% # 6 largest trees per sp
	summarise(dbh = mean(dbh, na.rm = TRUE), .groups = "drop") %>%
	mutate(r50 = dbh/2) %>%
	rename(dmax = "dbh") -> tree_max

write.csv(tree_max, here::here("data", "clean", "R50.csv"),
          row.names = FALSE)
