#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: remove-edge-traps.R
## Desc: Filter out traps from trap_connect.rds which are too close to the edge
## of the forest dynamics plot
## Date created: 2026-03-19


# Packages ----------------------------------------------------------------

library("tidyverse")
library("sf")


# Data --------------------------------------------------------------------

data <-
  readRDS(here::here("data", "clean", "trap_connect.rds"))


# Define edge plots -------------------------------------------------------

# define effective radius
eps <- 0.05 #captures 95% of influence

# define alpha
alphas <- exp(seq(log(0.07), log(0.40), length.out = 6))
alpha <- alphas[[3]]

# draw the 50-ha plot
fdplot <- rbind(c(0,0), c(0, 500), c(1000, 500), c(1000, 0), c(0,0))

# make it a polygon
plot_polygon <- sf::st_polygon(list(fdplot))

data_sf <- sf::st_as_sf(data, coords = c("x", "y"),
                       remove = FALSE)

# calculate distance from each trap to plot boundary
dist_to_edge_m <-
  sf::st_distance(data_sf,
                  sf::st_boundary(plot_polygon),
                  which = "Euclidean")

data_sf$dist_to_edge_m <- c(dist_to_edge_m)

data_edge <-
  data_sf %>%
  mutate(r_eff =  (- log(eps)) / alpha) %>%
  mutate(location = ifelse(dist_to_edge_m >= r_eff, "interior",
                           "edge")) %>%
  sf::st_drop_geometry()


# Check -------------------------------------------------------------------

data_edge %>%
  sf::st_drop_geometry() %>%
  select(trap, x, y , location) %>%
  distinct() %>%
  ggplot(aes(x = x, y = y, colour = location)) +
  geom_point()

data_edge %>%
  select(trap, location) %>%
  distinct() %>%
  group_by(location) %>%
  summarise(n_distinct(trap))

# 32 of 450 traps removed


# Download seed predator data  --------------------------------------------

# download from Dryad
trait_download <- rdryad::dryad_download(doi = "10.5061/dryad.230j5ch")

# function to put downloaded files in the right place
move_files <- function(old_path, new_path){
  file.rename( from = file.path(old_path) ,
               to = file.path(new_path, basename(old_path)) )
}

# apply the function to all files
lapply(trait_download, move_files, new_path = here::here("data", "raw"))


# Save --------------------------------------------------------------------

# format data for modelling and save

# add abundance data
data_abun <-
  read_csv(here::here("data", "clean", "species_abundance.csv")) %>%
  select(sp4, median_abundance)

# add seed pred data
pred_data <-
  read_csv(here::here("data", "raw", "TidyTrait.csv")) %>%
  mutate(sp4 = str_to_lower(Codigo)) %>%
  select(sp4, SeedPred_pres)

# with seed pred data
data_edge %>%
  filter(location == "interior") %>%
  select(- c(dist_to_edge_m, r_eff, location)) %>%
  left_join(data_abun) %>%
  left_join(pred_data) %>%
  drop_na(SeedPred_pres) %>%
  # transform variables for modelling
  mutate(log_total_seeds = log(total_seeds),
         log_median_abundance = log(median_abundance)) %>%
  mutate(
    conn_RC_sc_mat = scale(conn_RC),
    conn_RH_sc_mat = scale(conn_RH),
    conn_NRC_sc_mat = scale(conn_NRC),
    log_total_seeds_sc_mat = scale(log_total_seeds),
    log_median_abundance_sc_mat = scale(log_median_abundance)
  ) %>%
  mutate(
    conn_RC_sc = as.numeric(conn_RC_sc_mat),
    conn_RH_sc = as.numeric(conn_RH_sc_mat),
    conn_NRC_sc = as.numeric(conn_NRC_sc_mat),
    log_total_seeds_sc = as.numeric(log_total_seeds_sc_mat),
    log_median_abundance_sc = as.numeric(log_median_abundance_sc_mat),
    SeedPred_pres = as.factor(SeedPred_pres)
  ) %>%
  saveRDS(here::here("data", "clean", "trap_connect_buffer_preds.rds"))

# without seed pred data
data_edge %>%
  filter(location == "interior") %>%
  select(- c(dist_to_edge_m, r_eff, location)) %>%
  left_join(data_abun) %>%
  # transform variables for modelling
  mutate(log_total_seeds = log(total_seeds),
         log_median_abundance = log(median_abundance)) %>%
  mutate(
    conn_RC_sc_mat = scale(conn_RC),
    conn_RH_sc_mat = scale(conn_RH),
    conn_NRC_sc_mat = scale(conn_NRC),
    log_total_seeds_sc_mat = scale(log_total_seeds),
    log_median_abundance_sc_mat = scale(log_median_abundance)
  ) %>%
  mutate(
    conn_RC_sc = as.numeric(conn_RC_sc_mat),
    conn_RH_sc = as.numeric(conn_RH_sc_mat),
    conn_NRC_sc = as.numeric(conn_NRC_sc_mat),
    log_total_seeds_sc = as.numeric(log_total_seeds_sc_mat),
    log_median_abundance_sc = as.numeric(log_median_abundance_sc_mat)
  ) %>%
  saveRDS(here::here("data", "clean", "trap_connect_buffer.rds"))
