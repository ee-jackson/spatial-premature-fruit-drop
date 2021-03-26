#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: 0_compileTraps.R
## Desc: create dataset of proportion of seeds abscised per trap per sp per year
## Date: August 2020 

rm(list = ls())

require("tidyverse")

# seed rain data
seedRain <- read.table("../data/raw/BCI_TRAP200_20190215_spcorrected.txt",
 header=TRUE, stringsAsFactors = FALSE)

# seed traits from joe
seedTrait <- read.csv("../data/raw/20120227_seedsMassForTraits.csv", 
	header=TRUE, stringsAsFactors = FALSE)

# trap locations
load("../data/clean/trapLocations.RData")

## sort out names
seedRain <- rename(seedRain, SP4 = sp)

seedTrait <- seedTrait %>%
   select(LIFEFORM, N_SEEDFULL, SP4, GENUS, SPECIES)

seedRain$fecha <- as.character(seedRain$fecha)
seedRain$fecha <- as.Date(seedRain$fecha, "%Y-%m-%d")
seedRain$year <- format(as.Date(seedRain$fecha), "%Y")

seedDat <- left_join(seedRain, seedTrait, by = c("SP4"))

seedDat <- subset(seedDat, part == 1|part == 2|part == 5)

seedDat <- subset(seedDat, LIFEFORM != "EPIPHYTE"|LIFEFORM != "HEMIEPIPHYTE"|LIFEFORM != "HERB"|LIFEFORM != "VINE")

seedDat <- seedDat %>% 
	drop_na("N_SEEDFULL")

# don't have the full data for these years
seedDat <- subset(seedDat, seedDat$year != "1987" & seedDat$year != "2019")

# format trap id
seedDat$trap <- formatC(seedDat$trap, width = 3, format = "d", flag = "0")
seedDat$trap <- paste("trap", seedDat$trap, sep="_")

absdat <- seedDat %>% 
	group_by(part, SP4, year, trap, N_SEEDFULL, GENUS, SPECIES) %>%
	summarise(quantity_sum= sum(quantity, na.rm = TRUE)) %>%
	ungroup()

# Viable seeds

absdat_v <- subset(absdat, part == 1 | part == 2) %>% 
	replace(., is.na(.), 0) %>%
	rowwise() %>% 
	mutate(viable_seeds = ifelse(part==1, quantity_sum*N_SEEDFULL, quantity_sum)) %>%
	ungroup() %>%
	group_by(SP4, year, trap, GENUS, SPECIES) %>%
	summarise(viable_seeds= sum(viable_seeds, na.rm = TRUE)) %>%
	ungroup()

# Abscised seeds

absdat_a <- subset(absdat, part == 5) %>% 
	replace(., is.na(.), 0) %>%
	rowwise() %>% 
	mutate(abscised_seeds = quantity_sum*N_SEEDFULL) %>%
	ungroup() %>%
	group_by(SP4, year, trap,GENUS, SPECIES) %>%
	summarise(abscised_seeds= sum(abscised_seeds, na.rm = TRUE)) %>%
	ungroup()

# join and calculate proportion abcised
propdat <- full_join(absdat_a, absdat_v, by= c("SP4", "year", "trap", "GENUS", "SPECIES"))

propdat <- propdat %>%
	replace(., is.na(.), 0) %>%
	rowwise() %>% 
	mutate(total_seeds = sum(abscised_seeds, viable_seeds, na.rm = TRUE), proportion_abscised = abscised_seeds / total_seeds) %>%
	ungroup()

sumdat <- absdat %>% 
	group_by(SP4, year, trap) %>%
	summarise(sum_parts= sum(quantity_sum, na.rm = TRUE)) %>%
	ungroup()

propdat_loc <- left_join(trap_loc_quad, propdat, by = c("trap" = "trap"))

propdat_loc %>%
	left_join(sumdat, by = c("SP4", "year", "trap")) %>%
	filter(!sum_parts == 0) -> trapDat

save(trapDat, file = "../data/clean/trapData.RData")
