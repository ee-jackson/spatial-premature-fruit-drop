#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: connectivity2
## Desc: Calculate connectivity indicies for multiple years and sp
## Date: August 2020

rm(list = ls())

library("tidyverse")
library("rdist")
library("parallel")
library("lme4")

# load all bci files (they must be in the wd)
#file_names <- as.list(dir(pattern="bci.tree"))
#lapply(file_names, load, environment())

load("../data/bci.tree/bci.tree8.rdata")
load("../data/bci.tree/bci.tree7.rdata")

load("../data/bci.spptable.rdata")
trapDat <- read.csv("../output/tables/proportionAbscisedPerTrap.csv") 
RDA <- read.csv("../output/tables/R50.csv") 

#bci <- rbind(file_names)
bci <- rbind(bci.tree8, bci.tree7)

bci.tree <- bind_rows(file_names, .id = "df")
#add survey year
bci.tree %>%
    mutate(year = case_when(df = bci.tree8 ~ "2015",
                            df = bci.tree7 ~ "2010",
                            df = bci.tree6 ~ "2005",
                        	df = bci.tree5 ~ "2000",
                        	df = bci.tree4 ~ "1995",
                        	df = bci.tree3 ~ "1990",
                        	df = bci.tree2 ~ "1985",
                        	df = bci.tree1 ~ "1982"))


bci <- left_join(bci, bci.spptable, by = "sp")
bci <- left_join(bci, RDA, by = "sp")

# add tags to id values so they don't just look like numbers
trapDat$trap <- paste("trap", trapDat$trap, sep="_")
bci$treeID <- paste("tree", bci$treeID, sep="_")

###############################################################################
## Calculate euclidean distances
###############################################################################

#sp.list <- unique(bci[["sp"]])
#year.list <- unique(bci[["year"]])
sp.list <- list("cordbi", "oenoma")
year.list <- list("2015", "2010")

calculate_dist <- function (species, yr) {
	# subset tree data
	bd <- dplyr::filter(bci, (bci[,"status"]=="A") & (bci[,"sp"]==species) & (bci[,"dbh"]>bci[,"R50"]) & (bci[,"year"]==yr))

	# subset trap data
	td <- dplyr::filter(trapDat, (trapDat[,"year"]==yr) & (trapDat[,"SP6"]==species))

	dists <- cdist(bd[,c("gx", "gy")], td[,c("X", "Y")], metric="euclidean")
	dists.df <- as.data.frame(dists)
	rownames(dists.df) <- unlist(bd$treeID)
	colnames(dists.df) <- unlist(td$trap)
	cbind(bd, dists.df)
}

# apply function to all pairwise combinations of year and species
all.dists <- outer(sp.list, year.list, FUN = Vectorize(calculate_dist))
bci.dists <- bind_rows(all.dists, .id = "df")

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
	summarise(trap = paste(trap), NFI = sum(a)) %>% # sum over all trees to = 1 value per trap
	ungroup()
}

# create list of trap IDs to pass through the function
trap.list <- colnames(select(bci.dists, matches("trap_")))

# apply function to each trap
numCores <- detectCores()
NFIdat <- mclapply(trap.list, calculate_NFI, mc.cores = numCores)

# bind output into one dataframe
NFIdat.b <- do.call(rbind, NFIdat)

# merge with trap data
trapConnect <- left_join(NFIdat.b, trapDat,  by = c("trap", "year", "sp" = "SP6"))
head(trapConnect) #take a look

###############################################################################
## Hanski's Connectivity index
###############################################################################

# function to calculate CI
calculate_CI <- function (trap) {
	dplyr::filter(bci.dists, bci.dists[,trap] <= 100) %>%
	group_by(treeID, year, sp) %>% 
	mutate(a = dbh * exp( (-1/100) * eval(parse(text=trap)) ) ) %>% 
	ungroup() %>%
	dplyr::select(year, sp, trap, a) %>%
	group_by(year, sp) %>%
	summarise(trap = paste(trap), CI = sum(a))
}

CIdat <- mclapply(trap.list, calculate_CI)
CIdat.b <- do.call(rbind, CIdat) #bind_rows()
trapConnect <- left_join(trapConnect, CIdat.b, by = c("trap", "year", "sp"))
head(trapConnect) #take a look

###############################################################################
## Build and run some models
###############################################################################

m1 <- lme4::glmer(cbind(abscised_seeds, viable_seeds) ~ 
	NFI + (1|year) + (1|sp) + (1|year:sp), family = binomial(logit), 
	data = trapConnect)

m2 <- lme4::glmer(cbind(abscised_seeds, viable_seeds) ~ 
	CI + (1|year) + (1|sp) + (1|year:sp), family = binomial(logit), 
	data = trapConnect)

# compare
models <- list(m1, m2)
purrr::map_dfr(models, broom::tidy, 
							conf.int = TRUE, .id = "vars")
