#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculateConnectivity.R
## Desc: Script used to calculate connectivity indicies for each trap in every year
## Date: October 2020

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size=10))
library("rdist")
library("lme4")
library("ggeffects")

###############################################################################
# read in and compile data
###############################################################################

######## TREE DATA ########
# find and load all bci files in one go

file_names <- as.list(dir(
	path="/home/users/ft840275/SpatialPatterns/data/bci.tree", 
	pattern="bci.tree*"))
file_paths <- paste("/home/users/ft840275/SpatialPatterns/data/bci.tree/", file_names, sep="")

lapply(file_paths, load, environment())

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

# read in other data
load("../data/bci.spptable.rdata") # species data
RST <- read.csv("../data/R50.csv") # reproductive size threshold

# join bci with sp data and reproductive size threshold data
bci <- left_join(bci.bind, bci.spptable, by = "sp")
RST$sp <- tolower(RST$sp)
bci <- left_join(bci, RST, by = "sp")

# subset to only trees that are alive and are above reproductive size threshold
bci <- bci[!is.na(bci$R50), ]
bci <- dplyr::filter(bci, (bci[,"status"]=="A") & 
		(bci[,"dbh"]>bci[,"R50"]))

# add id tags to id values so they don't just look like numbers
bci$treeID <- paste("tree", bci$treeID, sep="_")

# will need year to be a character when we merge later
bci$year <- as.character(bci$year)

######## TRAP DATA ########
trapDat <- read.csv("../data/trapData.csv")
trapDat$year <- as.character(trapDat$year)

# save these so we don't have to do it every time
#save(trapDat, bci, file = "cleanTrapandTree.RData")
#load("cleanTrapandTree.RData")

####### CHOOSE SPECIES #########
trapDat %>%
	group_by(SP6, year) %>%
	mutate(n_traps = length(unique(trap))) %>%
	arrange(desc(n_traps)) %>%
	ungroup() %>%
	group_by(SP6) %>%
	mutate(median_traps = median(n_traps,na.rm = TRUE)) %>%
	select(SP6, median_traps) %>%
	distinct() %>%
	as.data.frame() -> trapSummary


bci %>%
	group_by(sp, survey_year) %>%
	mutate(n_trees = length(unique(treeID))) %>%
	arrange(desc(n_trees)) %>%
	ungroup() %>%
	group_by(sp) %>%
	mutate(median_trees = median(n_trees,na.rm = TRUE)) %>%
	select(median_trees, sp) %>%
	distinct() %>%
	as.data.frame() -> treeSummary

summary <- full_join(trapSummary, treeSummary, by= c("SP6" = "sp"))

subset(summary, median_trees >15 & median_traps >15) %>%
	pull(SP6) -> sp.list # this gives 41 species

# Subset to traps where there is more than 10 parts per year per sp
trapDat <- subset(trapDat, sum_parts>10)

###############################################################################
## Calculate euclidean distances
###############################################################################

calculate_dist <- function (species, yr) {
	# subset tree data
	bd <- dplyr::filter(bci, (bci[,"status"]=="A") & (bci[,"sp"]==species) & 
		(bci[,"dbh"]>bci[,"R50"]) & (bci[,"year"]==yr))

	# subset trap data
	td <- dplyr::filter(trapDat, (trapDat[,"year"]==yr) & 
		(trapDat[,"SP6"]==species))

	# create distance matrix using x and y co-ordinates
	dists <- cdist(bd[,c("gx", "gy")], td[,c("X", "Y")], metric="euclidean")
	dists.df <- as.data.frame(dists)

	# label col and row names with trap and tree IDs
	rownames(dists.df) <- unlist(bd$treeID)
	colnames(dists.df) <- unlist(td$trap)

	# bind with the bci data
	cbind(bd, dists.df)
}

# create lists of sp and year to pass to function
#sp.list <- unique(bci[["sp"]])
#sp.list <- c("cordbi", "mourmy")
year.list <- unique(bci[["year"]]) 

# apply function to all pairwise combinations of year and species
all.dists <- outer(sp.list, year.list, FUN = Vectorize(calculate_dist))
bci.dists <- bind_rows(all.dists)

###############################################################################
## Conspecific neighbourhood fecundity index
###############################################################################

# function to calculate NFI
calculate_NFI <- function (trap) {
	dplyr::filter(bci.dists, bci.dists[,trap] <= 100) %>% # subset to only trees within a 100m radius of the trap
	group_by(treeID, year, sp) %>% 
	mutate(a = dbh / eval(parse(text=trap))) %>% # divide dbh by dist
	ungroup() %>%
	dplyr::select(year, sp, trap, a) %>%
	group_by(year, sp) %>%
	summarise(trap = paste(trap), NFI = sum(a), .groups = "drop") # sum over all trees to = 1 value per trap
}

# create list of trap IDs to pass through the function
trap.list <- colnames(select(bci.dists, matches("trap_")))

# apply function to each trap and parallelize
numCores <- detectCores()
NFIdat <- mclapply(trap.list, calculate_NFI, mc.cores = numCores)

# bind output into one dataframe
NFIdat.b <- bind_rows(NFIdat)

# merge with trap data
trapConnect <- left_join(NFIdat.b, trapDat,  
	by = c("trap", "year", "sp" = "SP6"))

head(NFIdat.b) #take a look at NFI values

###############################################################################
## Hanski's Connectivity index
###############################################################################

# function to calculate CI
calculate_CI <- function (trap) {
	dplyr::filter(bci.dists, bci.dists[,trap] <= 100) %>% 
	group_by(treeID, year, sp) %>% 
	mutate(a = dbh * exp( (-1/50) * eval(parse(text=trap)) ) ) %>% 
	ungroup() %>%
	dplyr::select(year, sp, trap, a) %>%
	group_by(year, sp) %>%
	summarise(trap = paste(trap), CI = sum(a), .groups = "drop")
}

head(CIdat.b) #take a look at CI values

###############################################################################
## Merge and save dataset
###############################################################################

# merge together into one big dataset
CIdat <- mclapply(trap.list, calculate_CI, mc.cores = numCores)
CIdat.b <- bind_rows(CIdat)
trapConnect <- left_join(trapConnect, CIdat.b, by = c("trap", "year", "sp"))

# save it
save(trapConnect, file = "trapConnect.RData")
