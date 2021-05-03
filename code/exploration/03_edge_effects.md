Testing edge effect mitigation method
================
Eleanor Jackson
03 May, 2021

When a trap is close to the edge of the forest dynamics plot (FDP), we
do not get a complete picture of the connectivity of that trap. One
method to correct for this effect would be to scale the *observed*
connectivity index to a larger area - *predicting* connectivity. This
assumes that connectivity is homogeneous across a certain radius.

A potential problem with this method is the choice of radius to scale
to. At what point is our predicted value closest to the true value?

In this document I’m going to try and answer the following questions:

-   Is correcting for edge effects using this method better than just
    using the observed value?
-   What value of radius gives us the best prediction of the truth?

``` r
library("tidyverse"); theme_set(theme_bw(base_size=10))
library("sf")
library("rdist")
library("knitr")
library("patchwork")
library("here")

load(here::here("data", "clean", "treeData.RData"))
load(here::here("data", "clean", "trapConnect.RData"))
```

## Sample the data

So that I can compare predicted values to true values, I’m going to take
traps from the middle of the FDP, for which we already have “true”
measures of connectivity, and artificially place them at the edge. To do
this I will cut the FDP in half and discard all trees in that half. This
will mean that some of the traps which were in the middle of the plot
are now on the edge.

``` r
## TRAP DATA

# subset to traps which are on the LHS of the plot 
# and not near the 'real' edges
centreTraps <- subset(trapConnect, X < 500 & X > 100 & Y < 400 & Y > 100)

# subset to species which fulfil inclusion criteria
testTraps <- subset(centreTraps, SP4 %in% sp.list) 

# randomly select 1000 of these traps
testTraps <- sample_n(centreTraps, 1000)

testTraps %>% rename(CI_truth = CI) -> testTraps

## TREE DATA

# subset to species which fulfill inclusion criteria
subset(bci, SP4 %in% sp.list) -> bci50

# subset tree data to only include trees on the LHS of the plot
subset(bci50, gx < 500) -> bci50
```

## Observed connectivity

Re-calculate connectivity for each trap in the sample dataset, imagining
the plot is half as big. This will give us the *observed* connectivity.

``` r
## First claculate eulcidian distances

calculate_dist <- function (species, yr) {
    # subset tree data
    bd <- dplyr::filter(bci50, (bci50[,"SP4"]==species) & (bci50[,"year"]==yr))

    # subset trap data
    td <- dplyr::filter(testTraps, (testTraps[,"SP4"]==species) & 
        (testTraps[,"year"]==yr) )

    # create distance matrix using x and y co-ordinates
    dists <- rdist::cdist(bd[,c("gx", "gy")], td[,c("X", "Y")], 
        metric="euclidean")
    dists.df <- as.data.frame(dists)

    # label col names with trap IDs
    colnames(dists.df) <- unlist(td$trap)

    # bind with the bci data
    cbind(bd, dists.df)
}

# pull out every combo of yr and sp in edgeTraps to pass through function
testTraps %>% 
    count(year, SP4) %>% 
    select(SP4, year) %>% 
    rename(species = SP4, yr = year) -> dist_keys

# run function to calculate euclidean distances
dists <- pmap(dist_keys, calculate_dist)

# bind output dfs together
bci_dists <- dplyr::bind_rows(dists)

## Now calculate Hanski's connectivity index

calculate_CI <- function (trap) {
    drop_na(bci_dists, trap) %>%
    group_by(treeID, year) %>% 
    mutate(a = dbh * exp( (-0.02) * eval(parse(text=trap)) ) ) %>% 
    ungroup() %>%
    dplyr::select(year, SP4, trap, a) %>%
    group_by(year, SP4) %>%
    summarise(trap = paste(trap), CI_obs = sum(a), .groups = "drop")
}

# list of every trap to pass to function
trap_list <- unique(testTraps$trap)

# run function
CIdat <- lapply(trap_list, calculate_CI)

# bind output dfs together
CIdat_b <- bind_rows(CIdat)

# join our values of observed CI back to our values of true CI
left_join(CIdat_b, testTraps, by= c("year", "SP4", "trap")) -> sample
```

## Predict connectivity

Now we try and predict connectivity based on our observed values. I’ve
written a function to calculate the area of the intersection between the
radius around a trap and the FDP. It then calculates *predicted* CI by
dividing the *observed* CI by the area of the intersection and then
multiplying this by the area of the complete radius around the trap.

I run the function for 10 different values for ‘radius’ ranging from 20
to 200.

``` r
correct_for_edge <- function(species, yr, radius){
    # draw the fdp
    fdplot <- rbind(c(0,0), c(0, 500), c(500, 500), c(500, 0), c(0,0))
    # make it a polygon
    fdplot.sf <- sf::st_polygon(list(fdplot))
    # subset edgeTraps to sp i and yr i
    dplyr::filter(sample, 
        (.data[["SP4"]]==species) & (.data[["year"]]==yr)) %>%
    # make traps sf points
    sf::st_as_sf(coords = c("X", "Y")) %>% 
    # create a radius around the traps
    mutate(buf = st_buffer(geometry, dist = radius)) %>% 
    # calculate the area of intersection between each radius & the fdp
    mutate(area = sf::st_area( 
        sf::st_intersection(.data[["buf"]], fdplot.sf) )) %>% 
    # subset to traps that have an intersection / are on the edge
    filter(area > 0 & area < ( (pi * radius^2) - 1) ) %>%
    # calculate our prediction
    mutate(CI_pred = (CI_obs/area) * (pi * radius^2), radius = radius )
}

# create df of every combination of year sp and radius
sample %>% 
    count(year, SP4) %>% 
    select(SP4, year) %>% 
    rename(species = SP4, yr = year)    -> sample_keys

radius_keys <- tibble(radius = c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200))

sample_radius_keys <- merge(sample_keys, radius_keys)

# pass to function
testTraps_adjusted <- pmap(sample_radius_keys, correct_for_edge)
```

## Compare radii

Let’s take a look at how good our predictions are. I’m going to
calculate error as the difference between the *truth* and our
*prediction*. I’m also going to use the difference between the *truth*
and the *observed* connectivity as a control i.e. what would the error
be if we did nothing to correct for edge effects?

``` r
# bind output dfs and calculate error
bind_rows(testTraps_adjusted)%>%
    mutate(CI_error = CI_pred - CI_truth, radius = as.factor(radius)) -> output

# add a 'control' 
sample %>%
    mutate(CI_error = CI_obs - CI_truth, radius = as.factor("control")) %>%
    select(year, SP4, trap, CI_error, radius)  -> control

bind_rows(output, control) -> plotting_data

# plot error as jitter
ggplot(plotting_data) +
  geom_jitter(aes(x = radius, y = CI_error)) +
    labs(x = "radius", y = "error") +
  geom_hline(yintercept  = 0, colour = "red", linetype = 2) -> p1

# as violins
ggplot(plotting_data) +
  geom_violin(aes(x = radius, y = CI_error)) +
  labs(x = "radius", y = "error") +
  geom_hline(yintercept  = 0, colour = "red", linetype = 2) -> p2

# with x-axis ordered by decreasing median error
ggplot(plotting_data) +
  geom_jitter(aes(x = reorder(radius, abs(CI_error), median), y = CI_error)) +
  labs(x = "radius", y = "error") +
  geom_hline(yintercept  = 0, colour = "red", linetype = 2) -> p3

ggplot(plotting_data) +
  geom_violin(aes(x = reorder(radius, abs(CI_error), median), y = CI_error)) +
  labs(x = "radius", y = "error") +
  geom_hline(yintercept  = 0, colour = "red", linetype = 2) -> p4

(p1 + p2) / (p3 + p4)
```

![](figures/03_edge_effects/unnamed-chunk-5-1.png)<!-- -->

Correcting for edge effects using this method gives a lower median error
than using the observed value. The radii giving the lowest median errors
are between 100 - 140 m (the difference in error between these is
negligible).
