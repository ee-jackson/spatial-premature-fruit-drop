#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: 
## Desc: 
## Date: 

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size=10))

load("../data/trapConnect.RData")

trapConnect %>% dplyr::filter(sum_parts > 10) -> dat

ggplot(subset(trapDat, SP6 == "cordbi"| SP6=="oenoma" | SP6=="faraoc"| SP6=="mourmy"& sum_parts > 20)) +
	geom_line(aes(x=year, y=sum_parts, colour = trap, group = trap)) +
	theme(legend.position="none") +
	facet_wrap(~SP6)

trapDat %>% group_by(trap) %>% dplyr::filter(median(sum_parts) > 7 & SP6 == "cordbi") %>% drop_na() %>%
ggplot() +
	geom_boxplot(aes(x=trap, y=sum_parts)) +
	ggtitle("cordbi") -> p1

trapDat %>% group_by(trap) %>% dplyr::filter(median(sum_parts) > 6 & SP6 == "oenoma") %>% drop_na() %>%
ggplot() +
	geom_boxplot(aes(x=trap, y=sum_parts)) +
	ggtitle("oenoma") -> p2

trapDat %>% group_by(trap) %>% dplyr::filter(median(sum_parts) > 8 & SP6 == "faraoc") %>% drop_na() %>%
ggplot() +
	geom_boxplot(aes(x=trap, y=sum_parts)) +
	ggtitle("faraoc") -> p3

trapDat %>% group_by(trap) %>% dplyr::filter(median(sum_parts) > 6 & SP6 == "mourmy") %>% drop_na() %>%
ggplot() +
	geom_boxplot(aes(x=trap, y=sum_parts)) +
	ggtitle("mourmy") ->p4

library("patchwork")

(p1 +p2) / (p3 + p4)

ggplot(subset(trapDat, SP6 == "cordbi"| SP6=="oenoma" | SP6=="faraoc"| SP6=="mourmy"& sum_parts > 200)) +
	geom_boxplot(aes(x=trap, y=sum_parts)) +
	facet_wrap(~SP6)

trapDat <- read.csv("../output/tables/proportionAbscisedPerTrap.csv") 
trapDat$trap <- paste("trap", trapDat$trap, sep="_")
trapDat$year <- as.character(trapDat$year)
# Subset to traps where there is more than 1 part per year per sp
trapDat <- subset(trapDat, sum_parts>1)

#########################################

#is there correlation between traps that are close to each other?

# need: a measure of distance between traps, connectivity. 
# calculate euclidean distances between different traps
# distance on x asix and semivariance on y? Variograms

load("../output/tables/trapConnect.RData")
library("gstat")

subset(trapConnect, sp=="mourmy") %>%
	group_by(trap) %>%
	mutate(median_CI = median(CI)) -> dat

sp::coordinates(dat) = ~X+Y

alltraps <- trapConnect[!duplicated(trapConnect$trap),]


variogram(CI~1, data = dat) -> vgram
plot(vgram)

fit.variogram(vgram, vgm(c("Exp", "Mat", "Sph")), fit.kappa = TRUE)

library("spatstat")
MC<- moran.mc(s1$Income, lw, nsim=599)
moran.test(s1$Income,lw)

trap_pp <- ppp(dat$X, dat$Y, owin(c(0,1000),c(0,500)))
unitname(trap_pp) <- c("meter","meter")
plot(trap_pp)
QC<-quadratcount(trap_pp, xbreaks = seq(0,1000,50), ybreaks = seq(0,500,50))
plot(QC,add=T,col="red")

