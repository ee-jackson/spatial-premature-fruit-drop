#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: 
## Desc: 
## Date: 

rm(list = ls())

require(tidyverse)

load("../data/bci.tree/bci.tree8.rdata")
load("../data/bci.spptable.rdata")
trapDat <- read.csv("../output/tables/proportionAbscisedPerTrap.csv") 
RDA <- read.csv("../output/tables/R50.csv") 

bci8 <- left_join(bci.tree8, bci.spptable, by = "sp")
bci8 <- left_join(bci8, RDA, by = "sp")

bci8$gs <- paste(bci8$Genus, bci8$Species)
trapDat$gs <- paste(trapDat$GENUS, trapDat$SPECIES)

###############################################################################
# Trichilia tuberculata
bci8.TT <- subset(bci8, Latin=="Trichilia tuberculata")
bci8.TT <- subset(bci8.TT, status=="A")
bci8.TT <- subset(bci8.TT, dbh>R50) #322.9795918
trapDat.TT <- subset(trapDat, year=="2010" & sp=="TRI3")

ggplot() +
	geom_point(data = trapDat.TT, aes(x= X, y= Y, colour=proportion_abscised), 
		size=2) +
	geom_point(data = bci8.TT, aes(x= gx, y= gy), 
		size=1, shape=2) +
	ggtitle("Trichilia tuberculata") +
	theme_bw()

ggsave("../output/plots/treesAndTraps_TT.png")
# R50 is 322.9795918

###############################################################################
# write it into a function

# select species to loop over
species = c("Cordia alliodora","Trichilia tuberculata","Hybanthus prunifolius",
	"Faramea occidentalis", "Quararibea asterolepis", "Alseis blackiana", 
	"Oenocarpus mapora", "Cordia bicolor", "Beilschmiedia pendula", "Tabebuia rosea", "Anacardium excelsum")

species = c("cordal","tri2tu","hybapr", "faraoc", "quaras", "alsebl", 
	"oenoma", "cordbi", "beilpe", "tab1ro", "anacex")

species2 = c("hyeral","hirttr","mourmy","protte","casear","guapst","laetpr","apeime","dendar","des2pa","diptpa","guatdu")

species3 = c("heisco","pri2co","simaam","termob","paligu","brosal","cupasy","lacmpa","micoar","sponmo","tab2ar","virosu")

plot_sp_map <- function (species, speciesname) {
	bci <- subset(bci8, sp==species & status=="A" & dbh>=R50)
	trap <- subset(trapDat, year=="2010" & SP6==species & sum_parts>5)
	ggplot() +
		geom_point(data = trap, aes(x= X, y= Y, colour=proportion_abscised), 
			size=3) +
		geom_point(data = bci, aes(x= gx, y= gy), 
			size=2, shape=2) +
		ggtitle(bci$Latin) +
		theme_bw()
}

myplots <- imap(species3, plot_sp_map)
ggpubr::ggarrange(plotlist = myplots, common.legend=TRUE)
ggsave("../output/plots/treesAndTraps_2010_3.png")
