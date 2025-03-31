#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-trees.R
## Desc: create tree_data
## Date created: 2020-12-08

library("tidyverse") # v 1.3.1
library("here") # v 1.0.1

file_names <- as.list(dir(path = here::here("data", "raw", "bci.tree"),
                         pattern = "bci.tree*", full.names = TRUE))

lapply(file_names, load, environment())

species_list <- read.csv(here::here("data", "clean", "species_list.csv"))

R50 <- read.csv(here::here("data", "clean","R50.csv"))


# expand the bci survey data so we can match it to nearest years of the trap data
# we are repeating each row of the bci datasets 5 times and assigning 2 years either side

e.bci.tree3 <- bci.tree3[rep(seq_len(nrow(bci.tree3)), each = 5), ]
e.bci.tree3$year <- rep(c("1988","1989","1990","1991","1992"))

e.bci.tree4 <- bci.tree4[rep(seq_len(nrow(bci.tree4)), each = 5), ]
e.bci.tree4$year <- rep(c("1993","1994","1995","1996","1997"))

e.bci.tree5 <- bci.tree5[rep(seq_len(nrow(bci.tree5)), each = 5), ]
e.bci.tree5$year <- rep(c("2002","2001","2000","1999","1998"))

e.bci.tree6 <- bci.tree6[rep(seq_len(nrow(bci.tree6)), each = 5), ]
e.bci.tree6$year <- rep(c("2003","2004","2005","2006","2007"))

e.bci.tree7 <- bci.tree7[rep(seq_len(nrow(bci.tree7)), each = 5), ]
e.bci.tree7$year <- rep(c("2008","2009","2010","2011","2012"))

e.bci.tree8 <- bci.tree8[rep(seq_len(nrow(bci.tree8)), each = 5), ]
e.bci.tree8$year <- rep(c("2013","2014","2015","2016","2017"))

# bind all the bci datasets together and add column for survey year
bind_rows(mget(ls(pattern="e.bci.tree*")), .id='df') %>%
    mutate(survey_year = case_when(df == "e.bci.tree8" ~ "2015",
                        df == "e.bci.tree7" ~ "2010",
                        df == "e.bci.tree6" ~ "2005",
                        df == "e.bci.tree5" ~ "2000",
                        df == "e.bci.tree4" ~ "1995",
                        df == "e.bci.tree3" ~ "1990")) -> bci_bind

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
