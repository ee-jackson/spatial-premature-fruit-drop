#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: compile-trees.R
## Desc: create treeData
## Date created: 2020-12-08

library("tidyverse")

# find and load all bci files in one go

file_names <- as.list(dir(path = here::here("data", "raw", "bci.tree"),
	pattern="bci.tree*"))

file_paths <- paste("data", "raw", "bci.tree", file_names, sep="/")

lapply(file_paths, load, environment())

# read in other data
load(here::here("data", "raw", "bci.spptable.rdata")) # species data
RST <- read.csv(here::here("data", "clean","R50.csv")) # reproductive size threshold
nomenclature <- read.csv(here::here("data", "raw","20120305_nomenclature.csv")) #nomenclature

# expand the bci survey data so we can match it to nearest years of the trap data
# basically we are repeating each row of the bci datasets 5 times and assigning 2 years either side
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
                        df == "e.bci.tree3" ~ "1990")) -> bci.bind

rm(list = ls(pattern='*bci.tree*')) # remove some of those big files

# join bci with sp data and reproductive size threshold data
bci <- left_join(bci.bind, bci.spptable, by = "sp")
bci <- left_join(bci, RST, by = "sp")

# subset to only trees that are alive and are above reproductive size threshold
bci <- bci[!is.na(bci$R50), ]
bci <- dplyr::filter(bci, (bci[,"status"]=="A") |
		(bci[,"status"]=="AD") & # when a tree was noted as dead in one census but was found alive in a later census
		(bci[,"dbh"]>=bci[,"R50"]))

# add id tags to id values so they don't just look like numbers
bci$treeID <- formatC(bci$treeID, width = 3, format = "d", flag = "0")
bci$treeID <- paste("tree", bci$treeID, sep="_")

# remove unidentified sp
bci <- subset(bci, bci$sp != "uniden")

# add SP4 id codes
nomenclature %>%
  mutate(SP4 = tolower(SP4), SP6 = tolower(SP6)) %>%
  select(SP4, SP6) -> codes

left_join(bci, codes, by = c("sp" = "SP6")) -> bci

# save
save(bci, file = here::here("data", "clean", "treeData.RData"))