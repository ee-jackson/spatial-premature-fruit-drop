Abundance ~ species size?
================
Eleanor Jackson
17 July, 2025

``` r
library("tidyverse"); theme_set(theme_bw(base_size = 10))
library("patchwork")
library("here")
```

Sofia: \> BTW, do you think abundance will covary with tree size/growth
form? I suspect it will (large canopy trees would have much higher basal
area per individual)? I have been wondering if we might be more likely
to see stronger CNDD in smaller trees (understory trees) than canopy
trees, but if tree growth form and summed basal area covary then the
current results suggest that would not be the case.

``` r
R50 <- read.csv(here::here("data", "clean", "R50.csv"))
```

## Get abundance data

``` r
tree_data <- 
  readRDS(here::here("data", "clean", "tree_data.rds")) %>% 
  filter(year %in% c("1990", "1995", "2000", "2005", "2010",
                     "2010", "2015", "2022"))
```

``` r
# brms::parnames(model)
abun <- 
  tree_data %>% 
  group_by(year, sp4, sp6, genus, species) %>% 
  summarise(abundance = sum(basal_area_m2, na.rm = TRUE)) %>% 
  group_by(sp4, sp6, genus, species) %>% 
  summarise(median_abundance = median(abundance, na.rm = TRUE),
            mean_abundance = mean(abundance, na.rm = TRUE),
            max_abundance = max(abundance, na.rm = TRUE))
```

    ## `summarise()` has grouped output by 'year', 'sp4', 'sp6', 'genus'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'sp4', 'sp6', 'genus'. You can override
    ## using the `.groups` argument.

``` r
  abun %>% 
  inner_join(R50) %>% 
  left_join(read_csv(here::here("data", "clean", "species_list.csv"))) %>% 
  ggplot(aes(x = log(median_abundance), y = dmax,
             colour = lifeform)) +
  geom_point()  
```

    ## Joining with `by = join_by(sp6)`
    ## Rows: 103 Columns: 19
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (6): sp4, sp6, species, lifeform, family, genus dbl (2): seeds_per_fruit,
    ## capsules_per_fruit lgl (11): capsules, dsp_ant, dsp_bat, dsp_bird, dsp_bbird,
    ## dsp_explo, dsp_ma...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Joining with `by = join_by(sp4, sp6, genus, species)`

![](figures/22_abundance-size/unnamed-chunk-5-1.png)<!-- -->
