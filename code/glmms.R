#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: glmms.R
## Desc: binomial GLMMs of proportion of seeds abscised ~ connectivity
## Date: October 2020

rm(list = ls())

library("tidyverse"); theme_set(theme_bw(base_size=10))
library("lme4")
library("ggeffects")
library("gridExtra")

load("../data/clean/trapConnect2.RData")

###############################################################################
## Build and run some models
###############################################################################

trapConnect <- transform(trapConnect, 
    CI.sc = scale(CI),
    NFI.sc = scale(NFI)
    )

# a loop to run univariate binomial models in lme4
CImodels <- mclapply(setNames(sp.list, sp.list), function(var) {
	spDat <- as.data.frame(dplyr::filter(trapConnect, sp == var))
	formula <- paste("cbind(abscised_seeds, viable_seeds) ~ CI.sc + (1|trap) + (1|year)")
	lme4::glmer(formula, family = binomial(logit), data=spDat)
	},
	mc.cores = numCores
)

CImodels_q <- mclapply(setNames(sp.list, sp.list), function(var) {
	spDat <- as.data.frame(dplyr::filter(trapConnect, sp == var))
	formula <- paste("cbind(abscised_seeds, viable_seeds) ~ CI.sc + (1|trap) + (1|year) + (1|quadratID)")
	lme4::glmer(formula, family = binomial(logit), data=spDat)
	},
	mc.cores = numCores
)


NFImodels <- mclapply(setNames(sp.list, sp.list), function(var) {
	spDat <- as.data.frame(dplyr::filter(trapConnect, sp == var))
	formula <- paste("cbind(abscised_seeds, viable_seeds) ~ NFI.sc + (1|trap) + (1|year)")
	lme4::glmer(formula, family = binomial(logit), data=spDat)
	},
	mc.cores = numCores
)

NFImodels_q <- mclapply(setNames(sp.list, sp.list), function(var) {
	spDat <- as.data.frame(dplyr::filter(trapConnect, sp == var))
	formula <- paste("cbind(abscised_seeds, viable_seeds) ~ NFI.sc + (1|trap) + (1|year) + (1|quadratID)")
	lme4::glmer(formula, family = binomial(logit), data=spDat)
	},
	mc.cores = numCores
)

# compare
purrr::map_dfr(NFImodels_q, broom::tidy, 
							conf.int = TRUE, .id = "model")-> NFImodels_q_summary

purrr::map_dfr(CImodels_q, broom::tidy, 
							conf.int = TRUE, .id = "model") -> CImodels_q_summary

write.csv(NFImodels_q_summary, "NFImodels_q_summary.csv")

###############################################################################
## Plotting
###############################################################################

# unscale for plotting
# trapConnect$NFI.sc <- trapConnect$NFI.sc * 
# 	attr(trapConnect$NFI.sc, 'scaled:scale') + 
# 	attr(trapConnect$NFI.sc, 'scaled:center')

# function to plot all model fits CI
plot_model_fit.CI <- function (model, modelname) {
    output <- as.data.frame(ggpredict(model, terms = "CI.sc"))
    spDat <- as.data.frame(dplyr::filter(trapConnect, sp == modelname))
    ggplot() +
        geom_point(data= spDat, aes_string(x = "CI.sc", 
        					y = "proportion_abscised", size="sum_parts"), alpha= 0.5) +
        geom_line(data= output, aes_string(x = "x", y = "predicted"), 
        					colour= "blue") +
        geom_ribbon(data= output, aes_string(x= "x", ymin="conf.low", 
        					ymax="conf.high"), alpha= 0.1, fill="blue") +
        xlab("CI") +
        ylab("Proportion of seeds abscised") +
        theme_bw(base_size=10) +
        ggtitle(modelname)
}

# function to plot all model fits NFI
plot_model_fit.NFI <- function (model, modelname) {
    output <- as.data.frame(ggpredict(model, terms = "NFI.sc"))
    spDat <- as.data.frame(dplyr::filter(trapConnect, sp == modelname))
    ggplot() +
        geom_point(data= spDat, aes_string(x = "NFI.sc", 
        					y = "proportion_abscised", size="sum_parts"), alpha= 0.5) +
        geom_line(data= output, aes_string(x = "x", y = "predicted"), 
        					colour= "blue") +
        geom_ribbon(data= output, aes_string(x= "x", ymin="conf.low", 
        					ymax="conf.high"), alpha= 0.1, fill="blue") +
        xlab("NFI") +
        ylab("Proportion of seeds abscised") +
        theme_bw(base_size=10) +
        ggtitle(modelname)
}

myplots_NFI <- imap(NFImodels, plot_model_fit.NFI)
myplots_NFI_q <- imap(NFImodels_q, plot_model_fit.NFI)

myplots_CI <- imap(CImodels, plot_model_fit.CI)
myplots_CI_q <- imap(CImodels_q, plot_model_fit.CI)

# n <- length(myplots)
# nCol <- floor(sqrt(n))

# png("NFIplots.png", width = 2000, height = 1500)
# do.call("grid.arrange", c(myplots, ncol=nCol))
# dev.off()

png("NFIplots1_CI_q.png", width = 2000, height = 1500)
do.call(what = "grid.arrange", args = c(myplots_CI_q[c(1:16)], ncol = 4))
dev.off()

png("NFIplots2.png", width = 2000, height = 1500)
do.call(what = "grid.arrange", args = c(myplots[c(17:19)], ncol = 4))
dev.off()

png("NFIplots3.png", width = 2000, height = 1500)
do.call(what = "grid.arrange", args = c(myplots[c(20:28)], ncol = 4))
dev.off()

pdf("NFIplots_q.pdf", width = 10, height = 7)
marrangeGrob(grobs=myplots_NFI_q, nrow=3, ncol=3)
dev.off()

pdf("NFIplots.pdf", width = 10, height = 7)
marrangeGrob(grobs=myplots_NFI, nrow=3, ncol=3)
dev.off()