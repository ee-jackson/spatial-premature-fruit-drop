#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: make-fig-02
## Desc: make figure 2

# Packages ----------------------------------------------------------------

library("tidyverse")
library("brms")
library("patchwork")
library("thematic")
library("scales")

# set the ggplot theme
theme_set(
  theme_classic(base_size = 7) +
    theme(
      panel.spacing = unit(0, "lines"))
  )


# Get model ---------------------------------------------------------------

mod <-
  readRDS(here::here("output", "models", "pheno-repro-adjust",
                     "full_conn_binom_nseeds_abund.rds"))

# function to return variables to measurement scale
mod_data <-
  readRDS("data/clean/trap_connect_buffer.rds")

unscale <- function(x, var_name) {
  x * attr(mod_data[[var_name]], "scaled:scale") +
    attr(mod_data[[var_name]], "scaled:center")
}


# Make panel a ------------------------------------------------------------

# get predictions using posterior_epred() for each density index
# while keeping other predictors at mean
overlay_predictor_curves <- function(
    fit,
    predictors,
    x_range = c(-2, 2),
    resolution = 100,
    probs = c(0.025, 0.975),
    re_formula = NA,
    ndraws = NULL
) {
  dat <- fit$data

  x_seq <- seq(x_range[1], x_range[2], length.out = resolution)

  out <- map_dfr(predictors, function(pred) {

    # start with one row per x value
    newdata <- tibble(.x = x_seq)

    # fill all model variables with fixed values
    for (v in names(dat)) {
      if (v %in% predictors) {
        newdata[[v]] <- 0
      } else if (v == "total_seeds") {
        newdata[[v]] <- 1
      } else if (is.numeric(dat[[v]])) {
        newdata[[v]] <- mean(dat[[v]], na.rm = TRUE)
      } else if (is.factor(dat[[v]])) {
        newdata[[v]] <- factor(levels(dat[[v]])[1], levels = levels(dat[[v]]))
      } else if (is.character(dat[[v]])) {
        newdata[[v]] <- dat[[v]][which(!is.na(dat[[v]]))[1]]
      } else if (is.logical(dat[[v]])) {
        newdata[[v]] <- FALSE
      }
    }

    # vary only the focal predictor over the common grid
    newdata[[pred]] <- x_seq
    newdata$.x <- NULL

    epred <- posterior_epred(
      fit,
      newdata = newdata,
      re_formula = re_formula,
      ndraws = ndraws
    )

    tibble(
      x = x_seq,
      predictor = pred,
      estimate__ = apply(epred, 2, median),
      lower__    = apply(epred, 2, quantile, probs = probs[1]),
      upper__    = apply(epred, 2, quantile, probs = probs[2])
    )
  })

  out
}

# find common observed range of neighbourhood density indices
lower <- max(
  min(mod$data$conn_RC_sc, na.rm = TRUE),
  min(mod$data$conn_NRC_sc, na.rm = TRUE),
  min(mod$data$conn_RH_sc, na.rm = TRUE)
)

upper <- min(
  max(mod$data$conn_RC_sc, na.rm = TRUE),
  max(mod$data$conn_NRC_sc, na.rm = TRUE),
  max(mod$data$conn_RH_sc, na.rm = TRUE)
)

# run function over common observed range
ce3 <-
  overlay_predictor_curves(
    fit = mod,
    predictors = c("conn_RC_sc", "conn_NRC_sc", "conn_RH_sc"),
    x_range = round(c(lower, upper))
  ) %>%
  mutate(predictor = factor(predictor,
                            levels = c("conn_RC_sc",
                                       "conn_NRC_sc",
                                       "conn_RH_sc"))) %>%
  mutate(predictor = recode_factor(predictor,
                                   "conn_RC_sc" = "Reproductive conspecifics",
                                   "conn_NRC_sc" = "Non-reproductive conspecifics",
                                   "conn_RH_sc" = "Reproductive heterospecifics"))

# make plot
p1 <-
  ggplot(
    ce3,
    aes(x = x, y = estimate__, colour = predictor, fill = predictor)
  ) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, linewidth = 0) +
  geom_line(linewidth = 1) +
  labs(x = "Neighbourhood density", y = "Immature seed mortality",
       colour = "", fill = "") +
  scale_colour_manual(
    aesthetics = c("colour", "fill"),
    values = okabe_ito(3)
  ) +
  coord_cartesian(ylim = c(0,1), expand = 0) +
  theme(legend.position  = "inside",
        legend.position.inside = c(0.25, 0.9),
        legend.text = element_text(size = 7)) +
  scale_y_continuous(labels = scales::label_percent())


# Colour helper function --------------------------------------------------

# mix two colours in RGB space
mix_col <- function(col1, col2, p = 0.5) {
  rgb1 <- grDevices::col2rgb(col1)
  rgb2 <- grDevices::col2rgb(col2)

  mixed <- round((1 - p) * rgb1 + p * rgb2)

  grDevices::rgb(
    mixed[1], mixed[2], mixed[3],
    maxColorValue = 255
  )
}

# make 3 shades with the supplied colour fixed as the middle shade
shade_3 <- function(
    mid,
    light_mix = 0.55,  # more = lighter
    dark_mix  = 0.45   # more = darker
) {
  c(
    "mean - sd"    = mix_col(mid, "#FFFFFF", light_mix),
    "mean" = mid,
    "mean + sd"   = mix_col(mid, "#000000", dark_mix)
  )
}


# Plotting function for conditional effects -------------------------------

# will make panels b, c, d
make_plot <- function(data, effect, x, x_name, mid_colour) {
  data <-
    data[[effect]] %>%
    mutate(
           effect2__ = case_when(
             effect2__ == -1 ~ "mean - sd",
             effect2__ == 0 ~ "mean",
             effect2__ == 1 ~ "mean + sd"
           ),
           effect2__ = factor(effect2__,
                              levels = c("mean + sd", "mean", "mean - sd"))) %>%
    mutate(x = unscale(x = effect1__, var_name = x))

  vals <- shade_3(mid = mid_colour)

  ggplot(data,
         aes(x = x, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__,
                    fill = as.factor(effect2__)), alpha = 0.3) +
    geom_line(aes(colour = as.factor(effect2__)), linewidth = 0.5) +
    labs(x = x_name, y = "Immature seed mortality",
         colour = "Total conspecific seeds",
         fill = "Total conspecific seeds") +
    scale_colour_manual(
      aesthetics = c("colour", "fill"),
      values = vals
    ) +
    coord_cartesian(ylim = c(0,1),
                    expand = 0) +
    scale_y_continuous(labels = scales::label_percent())
}

cond_eff <- conditional_effects(mod)

p2 <- make_plot(data = cond_eff, effect = "conn_RC_sc:log_total_seeds_sc",
                x = "conn_RC_sc_mat",
                x_name = "Density of reproductive conspecifics",
                mid_colour = okabe_ito(3)[1])

p3 <- make_plot(data = cond_eff, effect = "conn_NRC_sc:log_total_seeds_sc",
                x = "conn_NRC_sc_mat",
                x_name = "Density of non-reproductive conspecifics",
                mid_colour = okabe_ito(3)[2])

p4 <- make_plot(data = cond_eff, effect = "conn_RH_sc:log_total_seeds_sc",
                x = "conn_RH_sc_mat",
                x_name = "Density of reproductive heterospecifics",
                mid_colour = okabe_ito(3)[3])


# Combine panels ----------------------------------------------------------

png(
  here::here("output", "figures", "figure_02.png"),
  width = 173,
  height = 150,
  units = "mm",
  type = "cairo",
  res = 600
)

(p1 | (p2 / p3 / p4) ) +
  plot_layout(widths = c(2, 0.75)) +
  plot_annotation(tag_levels = "a")

dev.off()
