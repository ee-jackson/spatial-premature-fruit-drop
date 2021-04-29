#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculateConnectivity.R
## Desc: Script used to calculate connectivity indicies for each trap in every year
## Date: October 2020

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size = 8))
library("rdist")
library("parallel")

######## LOAD DATA ########
load("../data/clean/trapData.RData")
load("../data/clean/treeData.RData")

trapDat$year <- as.character(trapDat$year)


####### CHOOSE SPECIES #########
trapDat %>%
	group_by(SP4, year) %>%
	mutate(n_traps = length(unique(trap))) %>%
	arrange(desc(n_traps)) %>%
	ungroup() %>%
	group_by(SP4) %>%
	mutate(median_traps = median(n_traps, na.rm = TRUE)) %>%
	select(SP4, median_traps) %>%
	distinct() %>%
	as.data.frame() -> trapSummary


bci %>%
	group_by(SP4, survey_year) %>%
	mutate(n_trees = length(unique(treeID))) %>%
	arrange(desc(n_trees)) %>%
	ungroup() %>%
	group_by(SP4) %>%
	mutate(median_trees = median(n_trees, na.rm = TRUE)) %>%
	select(median_trees, SP4) %>%
	distinct() %>%
	as.data.frame() -> treeSummary

# I think the NAs in the summary happen if a sp is in bci but not in trapDat and vice versa... but why would this happen?
summary <- full_join(trapSummary, treeSummary, by= "SP4")

subset(summary, median_trees >15 & median_traps >15) %>%
	pull(SP4) -> sp.list # this gives 40 species

###############################################################################
## Calculate euclidean distances
###############################################################################

calculate_dist <- function (species, yr) {
	# subset tree data
	bd <- dplyr::filter(bci, (bci[,"SP4"]==species) & (bci[,"year"]==yr))

	# subset trap data
	td <- dplyr::filter(trapDat, (trapDat[,"SP4"]==species) & 
		(trapDat[,"year"]==yr) )

	# create distance matrix using x and y co-ordinates
	dists <- cdist(bd[,c("gx", "gy")], td[,c("X", "Y")], metric="euclidean")
	dists.df <- as.data.frame(dists)

	# label col names with trap IDs
	colnames(dists.df) <- unlist(td$trap)

	# bind with the bci data
	cbind(bd, dists.df)
}

# create lists of sp and year to pass to function
year.list <- unique(bci[["year"]]) 

# apply function to all pairwise combinations of year and species
# will return list of dfs
all.dists <- outer(sp.list, year.list, FUN = Vectorize(calculate_dist))

bci.dists <- dplyr::bind_rows(all.dists)

###############################################################################
## Conspecific neighbourhood fecundity index
###############################################################################

# function to calculate NFI
calculate_NFI <- function (trap) {
	dplyr::filter(bci.dists, bci.dists[,trap] <= 50) %>% # subset to only trees within a 50m radius of the trap
	group_by(treeID, year) %>% 
	mutate(a = dbh / eval(parse(text=trap))) %>% # divide dbh by dist
	ungroup() %>%
	dplyr::select(year, SP4, trap, a) %>%
	group_by(year, SP4) %>%
	summarise(trap = paste(trap), NFI = sum(a), .groups = "drop") # sum over all trees to = 1 value per trap
}

# create list of trap IDs to pass through the function
trap.list <- unique(trapDat$trap)

# apply function to each trap and parallelize
NFIdat <- mclapply(trap.list, calculate_NFI, mc.cores = 4)

# bind output into one dataframe
NFIdat.b <- bind_rows(NFIdat)

head(NFIdat.b) #take a look at NFI values

# merge with trap data
trapConnect <- left_join(NFIdat.b, trapDat,  
	by = c("trap", "year", "SP4"))

###############################################################################
## Hanski's Connectivity index
###############################################################################

# function to calculate CI
calculate_CI <- function (trap) {
	drop_na(bci.dists, trap) %>%
	group_by(treeID, year) %>% 
	mutate(a = dbh * exp( (-0.02) * eval(parse(text=trap)) ) ) %>% 
	ungroup() %>%
	dplyr::select(year, SP4, trap, a) %>%
	group_by(year, SP4) %>%
	summarise(trap = paste(trap), CI = sum(a), .groups = "drop")
}

CIdat <- parallel::mclapply(trap.list, calculate_CI, mc.cores = 4)

CIdat.b <- bind_rows(CIdat)

head(CIdat.b) #take a look at CI values

###############################################################################
## Merge and save dataset
###############################################################################

# merge together into one big dataset
trapConnect <- left_join(trapConnect, CIdat.b, by = c("trap", "year", "SP4"))

# save it
save(trapConnect, sp.list, file = "../data/clean/trapConnect.RData")
