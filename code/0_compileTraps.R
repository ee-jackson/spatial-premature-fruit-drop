#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: 0_compileTraps.R
## Desc: create dataset of proportion of seeds abscised per trap per sp per year
## Date: August 2020 

rm(list = ls())

require(tidyverse)

seedRain<-read.table("../../PrematureFruitAbscission/data/BCI_TRAP200_20190215_spcorrected.txt", header=TRUE, stringsAsFactors = FALSE)

trapLoc <- read.csv("../data/clean/trapLocations.csv")

seedTrait<-read.csv("../../PrematureFruitAbscission/data/20120227_seedsMassForTraits.csv", header=TRUE, stringsAsFactors = FALSE)

seedTrait <- seedTrait %>%
   select(LIFEFORM, N_SEEDFULL, SP4, SP6, GENUS, SPECIES)


seedRain$fecha <- as.character(seedRain$fecha)
seedRain$fecha <- as.Date(seedRain$fecha, "%Y-%m-%d")
seedRain$year <- format(as.Date(seedRain$fecha), "%Y")

# only use these 250 traps
seedRain %>%
	dplyr::filter(trap <= 200 | trap >=300 & trap <= 349) -> seedRain

seedDat <- left_join(seedRain, seedTrait, by = c("sp" = "SP4"))

seedDat <- subset(seedDat, part==1|part==2|part==5)
seedDat <- subset(seedDat,LIFEFORM== "LIANA"|LIFEFORM== "MIDSTORY"|LIFEFORM== "SHRUB"|LIFEFORM== "TREE"|LIFEFORM=="UNDERSTORY")

seedDat <- seedDat %>% 
	drop_na("N_SEEDFULL")

# don't have the full data for these years
seedDat <- subset(seedDat, seedDat$year != "1987" & seedDat$year != "2019")

absdat <- seedDat %>% 
	group_by(part, sp, SP6, year, trap, N_SEEDFULL, GENUS, SPECIES) %>%
	summarise(quantity_sum= sum(quantity, na.rm = TRUE)) %>%
	ungroup()

# Viable seeds

absdat.v <- subset(absdat, part== 1 | part== 2) %>% 
	replace(., is.na(.), 0) %>%
	rowwise() %>% 
	mutate(viable_seeds = ifelse(part==1, quantity_sum*N_SEEDFULL, quantity_sum)) %>%
	ungroup() %>%
	group_by(sp, SP6, year, trap, GENUS, SPECIES) %>%
	summarise(viable_seeds= sum(viable_seeds, na.rm = TRUE)) %>%
	ungroup()

# Abscised seeds

absdat.a <- subset(absdat, part== 5) %>% 
	replace(., is.na(.), 0) %>%
	rowwise() %>% 
	mutate(abscised_seeds = quantity_sum*N_SEEDFULL) %>%
	ungroup() %>%
	group_by(sp, SP6, year, trap,GENUS, SPECIES) %>%
	summarise(abscised_seeds= sum(abscised_seeds, na.rm = TRUE)) %>%
	ungroup()

# join and calculate proportion abcised
propDat <- full_join(absdat.a, absdat.v, by= c("sp", "SP6", "year", "trap", "GENUS", "SPECIES"))

propDat <- propDat %>%
	replace(., is.na(.), 0) %>%
	rowwise() %>% 
	mutate(total_seeds = sum(abscised_seeds, viable_seeds, na.rm = TRUE), proportion_abscised = abscised_seeds / total_seeds) %>%
	ungroup()

sumdat <- absdat %>% 
	group_by(sp, year, trap) %>%
	summarise(sum_parts= sum(quantity_sum, na.rm = TRUE)) %>%
	ungroup()

trapDat <- left_join(propDat, trapLoc, by = c("trap" = "trap"))
trapDat2 <- left_join(trapDat, sumdat, by = c("sp", "year", "trap"))
trapDat2$SP6 <- tolower(trapDat2$SP6)

write.csv(trapDat2,"../data/clean/proportionAbscisedPerTrap.csv")
