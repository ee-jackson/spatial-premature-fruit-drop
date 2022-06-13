Testing different values of alpha
================
Eleanor Jackson
10 June, 2022

In calculating Hanski’s connectivity index we had to include a measure
of the hypothesised dispersal ability of host-specific enemies
associated with seeds of the focal tree species. alpha = (1/average
dispersal distance in m). Dispersal distances of tropical seed predators
are unknown. We chose **alpha = 0.02**, which corresponds to an average
migration distance of 50 m. This choice was based on empirical data on
the average dispersal distance of a temperate leaf mining moth
associated with oak trees (Gripenberg and Roslin 2005).

T.O.: *If it is not too computationally intensive, I would be tempted to
pick a few values and re-run the models, then pick the one that has the
best fit and use that in main paper. This sensitivity analysis to find
the best alpha value can be described in Supp Info and shouldn’t be too
arduous? Will avoid criticism from reviews that this choice is
arbitrary.*

Lets do 20m, 30m, 40m, 50m, 60m, 70m. So alpha = 0.05, 0.033, 0.025,
0.016, 0.014.

## Compare the predictive accuracy of the models using LOO-CV

``` r
model_20m <- readRDS(here::here("output", "models", "zoib_capsules_20m.rds"))
model_30m <- readRDS(here::here("output", "models", "zoib_capsules_30m.rds"))
model_40m <- readRDS(here::here("output", "models", "zoib_capsules_40m.rds"))
model_50m <- readRDS(here::here("output", "models", "zoib_capsules_quad_rslope_nest.rds"))
model_60m <- readRDS(here::here("output", "models", "zoib_capsules_60m.rds"))
model_70m <- readRDS(here::here("output", "models", "zoib_capsules_70m.rds"))

library("tidyverse"); theme_set(theme_bw(base_size = 10))
library("broom.mixed")
library("brms")
library("ggdist")
library("loo") # v 2.4.1
library("patchwork")
options(mc.cores = 3)
```

## Compare the predictive accuracy of the models using Leave-One-Out Cross Validation

``` r
comp <- loo_compare(model_20m, model_30m, model_40m, model_50m, model_60m, model_70m)

print(comp, digits = 3)
```

    ##           elpd_diff se_diff
    ## model_20m   0.000     0.000
    ## model_30m -18.362     4.619
    ## model_40m -30.714     6.960
    ## model_50m -39.711     8.613
    ## model_60m -46.138     9.674
    ## model_70m -50.620    10.494

``` r
comp %>% 
  data.frame() %>% 
  rownames_to_column(var = "model_name") %>% 
  ggplot(aes(x    = model_name, 
             y    = elpd_loo, 
             ymin = elpd_loo - se_elpd_loo, 
             ymax = elpd_loo + se_elpd_loo)) +
  geom_pointrange(shape = 21, fill = "white") +
  coord_flip() +
  labs(x = NULL, title = "expected log predictive density (ELPD)") 
```

![](figures/12_testing-alpha-values/loo-compare-1.png)<!-- -->

## Compare estimates

``` r
tibble(model = c("model_20m", "model_30m", "model_40m", "model_50m", "model_60m", "model_70m")) %>% 
  mutate(fit   = purrr::map(model, get)) %>% 
  mutate(tidy  = purrr::map(fit, tidy)) %>% 
  unnest(tidy) %>% 
  filter(effect == "fixed" & !grepl("(Intercept)", term)) -> my_coef_tab

my_coef_tab
```

    ## # A tibble: 24 × 10
    ##    model     fit       effect component group term   estimate std.error conf.low
    ##    <chr>     <list>    <chr>  <chr>     <chr> <chr>     <dbl>     <dbl>    <dbl>
    ##  1 model_20m <brmsfit> fixed  cond      <NA>  conne…  0.0297     0.0549  -0.0775
    ##  2 model_20m <brmsfit> fixed  cond      <NA>  phi_c…  0.0108     0.0841  -0.169 
    ##  3 model_20m <brmsfit> fixed  cond      <NA>  zoi_c… -0.560      0.117   -0.816 
    ##  4 model_20m <brmsfit> fixed  cond      <NA>  coi_c…  0.678      0.240    0.217 
    ##  5 model_30m <brmsfit> fixed  cond      <NA>  conne…  0.0461     0.0620  -0.0779
    ##  6 model_30m <brmsfit> fixed  cond      <NA>  phi_c… -0.00162    0.108   -0.237 
    ##  7 model_30m <brmsfit> fixed  cond      <NA>  zoi_c… -0.609      0.154   -0.937 
    ##  8 model_30m <brmsfit> fixed  cond      <NA>  coi_c…  0.761      0.312    0.211 
    ##  9 model_40m <brmsfit> fixed  cond      <NA>  conne…  0.0610     0.0721  -0.0859
    ## 10 model_40m <brmsfit> fixed  cond      <NA>  phi_c… -0.00811    0.131   -0.288 
    ## # … with 14 more rows, and 1 more variable: conf.high <dbl>

``` r
my_coef_tab %>% 
  mutate(term = recode(term, connectivity_pred_sc = "mu", coi_connectivity_pred_sc = "coi",
                             zoi_connectivity_pred_sc = "zoi", phi_connectivity_pred_sc = "phi")) %>% 
  ggplot(aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(shape = 21, fill = "white") +
  labs(x = NULL,
       y = NULL) +
  geom_hline(yintercept = 0,  color = "blue") +
  coord_flip() +
  theme_classic() +
  facet_wrap(~term, ncol = 1)
```

![](figures/12_testing-alpha-values/estimate-compare-1.png)<!-- -->
