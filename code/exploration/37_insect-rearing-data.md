Insect rearing data
================
Eleanor Jackson
25 September, 2026

``` r
library("tidyverse")
library("patchwork")
library("glmmTMB")
```

Reviewers are concerned that the causality between immature seed
mortality and natural enemies is not established.

Sofia: We could perhaps tap into the data collected in the context of
[our food web
study](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13359). We
have insect rearing data and fruit/seed dissection data which is
assembled separately for mature and immature fruits. The data sets are
not perfect, but should still allow us to demonstrate that for most
species, immature fruits are *way* more likely to yield insects (or show
signs of insect attack) than mature fruits, which could perhaps help us
in building the case.

1.  Seeds_saved_for_rearing(…).xlsx: Contains a column called
    Nr_units_predated which shows data from dissections (rather
    conservative – sometimes the material was in quite a bad state after
    3 months of rearing, so hard to see evidence of predation). Metadata
    worksheet provides column descriptions.

2.  Insects_final_cleaned(…).xlsx: Each row shows the number of insects
    of a given morphospecies collected from a given pot on a given day
    (and therefore placed in the same vial for storage) – if lots of
    insects there might be multiple vials (and rows in the file) for a
    morphospecies x rearing pot x day combination. The column Quantity
    gives the number of individuals. SeedPot column is the same as
    Pot_ID in the previous file. (You will need to link the two to get
    information on whether seeds/fruits were mature or immature on
    collection.

3.  Seed predator summary(…).xlsx: I can’t remember why I did this
    summary and hope that it does not contain any errors, but it looks
    like some info here might be potentially useful. Worksheet
    IndsPerMorphospPerMaturityStage seems potentially useful, and
    IndsPerMorphosp_seedpreds only has filtered out morphospecies that
    we believe are feeding on seeds (not the fleshy parts of the fruit
    like e.g. most flies and some moth families do). Could possibly
    remove morphospecies that don’t appear in this worksheet from data
    summaries?

With this data, we can ask if more insects are reared from immature than
from mature fruits/seeds. Can also restrict insects to seed predators
and exclude fruit-feeding insects. Our response is: insect emergence
rate from pots containing immature versus mature fruits/seeds
i.e. number of emerged insects divided by number of seeds.

Read in the data:

``` r
data_seed <- readxl::read_xlsx(
  here::here("data", "raw", "insect-data", 
             "Seeds_saved_for_rearing_Cleaned_June2014.xlsx"),
  na = c("", "NA")) |> 
  janitor::clean_names()
```

    ## New names:
    ## • `Month` -> `Month...16`
    ## • `Month` -> `Month...20`

``` r
data_insect <- readxl::read_xlsx(
  here::here("data", "raw", "insect-data", 
             "Insects_final_cleaned_September2015.xlsx"),
  na = c("", "NA")) |> 
  janitor::clean_names()

seed_pred_morphos <- readxl::read_xlsx(
  sheet = 4,
  here::here("data", "raw", "insect-data", 
             "Seed_predator_summary_16July2014.xlsx")) |> 
  janitor::clean_names() |> 
  filter(seed_predator == "Yes")
```

It would be good to estimate the number of seeds within collected
fruits.

``` r
# Joe's fruit dissection data
seed_trait <-
  read_tsv(
    here::here(
      "data",
      "raw",
      "seed-masses",
      "Fruit_Seed_masses_20210111_BanisteriopsisCorrectedToBactrismajor.txt"
    )
  ) %>%
  rename_with(tolower) %>%
  mutate(id = paste(sp4, ind.unique, fruit.unique, sep = "_")) %>%
  select(sp4, id, n_seedfull, n_capsules) %>%
  distinct() %>%
  group_by(sp4) %>%
  summarise(seeds_per_fruit = mean(n_seedfull, na.rm = TRUE)) |> 
  select(sp4, seeds_per_fruit)
```

    ## Rows: 19087 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (7): sp4, lugar, sp6, family, genus, species, lifeform
    ## dbl (19): ind.unique, fruit.unique, date, IND, fruit_frsh, fruit_dry, N_seed...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Calculate seeds per pot:

``` r
data_seed <- 
  data_seed |> 
  left_join(seed_trait, by = c("codigo" = "sp4")) |> 
  rowwise() |> 
  mutate(nr_equ_seeds = 
           ifelse(unit == "Fruit", nr_fruits * seeds_per_fruit, nr_seeds)) 
```

There 124 species out of 483 for which we do not have data on “seeds per
fruit”.

``` r
# join pot data with insect data
insects_by_pot <- 
  data_insect %>%
  filter(!is.na(seed_pot)) %>%
  group_by(seed_pot) %>%
  summarise(
    n_insects = sum(quantity, na.rm = TRUE),
    .groups = "drop"
  )

data_join <-
  data_seed |> 
  left_join(insects_by_pot, by = join_by("pot_id" == "seed_pot")) |> 
  mutate(n_insects = replace_na(n_insects, 0)) |> 
  filter(nr_equ_seeds > 0) 

data_join |> 
  ggplot(aes(x = n_insects, fill = maturity_type)) +
  geom_density(alpha = 0.4) +
  ylab("N pots") +
  coord_cartesian(xlim = c(0,20))
```

![](figures/37_insect-rearing-data/unnamed-chunk-5-1.png)<!-- -->

Most pots had zero insects.

``` r
data_join |> 
  group_by(maturity_type) |> 
  summarise(sum_seeds = sum(nr_equ_seeds, na.rm = TRUE),
            n_insects = sum(n_insects, na.rm = TRUE),
            sum_fruits = sum(nr_fruits, na.rm = TRUE)) |> 
  mutate(seed_emergence_rate = n_insects/sum_seeds,
         fruit_emergence_rate = n_insects/sum_fruits) |> 
  print(width = Inf)
```

    ## # A tibble: 2 × 6
    ##   maturity_type sum_seeds n_insects sum_fruits seed_emergence_rate
    ##   <chr>             <dbl>     <dbl>      <dbl>               <dbl>
    ## 1 inmaduro       1368844.      9329     100604             0.00682
    ## 2 maduro          580245.     13266      39280             0.0229 
    ##   fruit_emergence_rate
    ##                  <dbl>
    ## 1               0.0927
    ## 2               0.338

Emergence rate is higher for mature material, this would mean 0.02
insects emerged per mature seed on average vs 0.007 per immature seed
and 0.33 per mature fruit vs 0.09 per immature fruit.

Try restricting to only seed predators:

``` r
insects_by_pot_seedpreds <- 
  data_insect %>%
  filter(morphospecies %in% seed_pred_morphos$morphospecies) |> 
  filter(!is.na(seed_pot)) %>%
  group_by(seed_pot) %>%
  summarise(
    n_insects = sum(quantity, na.rm = TRUE),
    .groups = "drop"
  ) 

data_join_seedpreds <-
  data_seed |> 
  left_join(insects_by_pot_seedpreds, by = join_by("pot_id" == "seed_pot")) |> 
  mutate(n_insects = replace_na(n_insects, 0))|> 
  filter(nr_equ_seeds > 0)

data_join_seedpreds |> 
  group_by(maturity_type) |> 
  summarise(sum_seeds = sum(nr_equ_seeds, na.rm = TRUE),
            n_insects = sum(n_insects, na.rm = TRUE),
            sum_fruits = sum(nr_fruits, na.rm = TRUE)) |> 
  mutate(seed_emergence_rate = n_insects/sum_seeds,
         fruit_emergence_rate = n_insects/sum_fruits) |> 
  print(width = Inf)
```

    ## # A tibble: 2 × 6
    ##   maturity_type sum_seeds n_insects sum_fruits seed_emergence_rate
    ##   <chr>             <dbl>     <dbl>      <dbl>               <dbl>
    ## 1 inmaduro       1368844.      4143     100604             0.00303
    ## 2 maduro          580245.      6705      39280             0.0116 
    ##   fruit_emergence_rate
    ##                  <dbl>
    ## 1               0.0412
    ## 2               0.171

Try fitting a model to account for differences in sampling effort and
species.

First for fruits only:

``` r
data_fit <- data_join_seedpreds %>%
  mutate(
    maturity_type = factor(maturity_type),
    maturity_type = relevel(maturity_type, ref = "maduro")
  )

mfruit <- glmmTMB(
  n_insects ~ maturity_type +
    offset(log(nr_fruits)) +
    (1 | codigo),
  family = nbinom2,
  data = filter(data_fit, nr_fruits > 0)
)

summary(mfruit)
```

    ##  Family: nbinom2  ( log )
    ## Formula:          
    ## n_insects ~ maturity_type + offset(log(nr_fruits)) + (1 | codigo)
    ## Data: filter(data_fit, nr_fruits > 0)
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##    9026.7    9053.6   -4509.3    9018.7      6200 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  codigo (Intercept) 14.32    3.784   
    ## Number of obs: 6204, groups:  codigo, 346
    ## 
    ## Dispersion parameter for nbinom2 family (): 0.172 
    ## 
    ## Conditional model:
    ##                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)            -5.7835     0.4482 -12.903  < 2e-16 ***
    ## maturity_typeinmaduro  -0.7239     0.1251  -5.785 7.27e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
exp(confint(mfruit, parm = "maturity_typeinmaduro"))
```

    ##                          2.5 %    97.5 %  Estimate
    ## maturity_typeinmaduro 0.379395 0.6196318 0.4848559

The fitted insect emergence rate per fruit in immature material is about
52% lower than in mature fruits.

Now including seeds plus the estimated seeds within fruits
(nr_equ_seeds):

``` r
m1 <- glmmTMB(
  n_insects ~ maturity_type +
    offset(log(nr_equ_seeds)) +
    (1 | codigo),
  family = nbinom2,
  data = drop_na(data_fit, seeds_per_fruit)
)

summary(m1)
```

    ##  Family: nbinom2  ( log )
    ## Formula:          n_insects ~ maturity_type + offset(log(nr_equ_seeds)) + (1 |  
    ##     codigo)
    ## Data: drop_na(data_fit, seeds_per_fruit)
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##   10932.9   10961.2   -5462.4   10924.9      8670 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  codigo (Intercept) 14.73    3.838   
    ## Number of obs: 8674, groups:  codigo, 359
    ## 
    ## Dispersion parameter for nbinom2 family (): 0.172 
    ## 
    ## Conditional model:
    ##                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)            -7.2642     0.4012 -18.104  < 2e-16 ***
    ## maturity_typeinmaduro  -0.5461     0.1103  -4.951 7.38e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
exp(confint(m1, parm = "maturity_typeinmaduro"))
```

    ##                           2.5 %    97.5 %  Estimate
    ## maturity_typeinmaduro 0.4665846 0.7189746 0.5791912

The estimated insect emergence rate for immature material is 0.58 times
the rate for mature material. i.e., the emergence rate is about 42%
lower in immature material, and it’s a significant effect

``` r
data_seed %>%
  filter(nr_equ_seeds > 0) |> # only pots containing seeds
  group_by(codigo) %>%
  summarise(
    n_maturity = n_distinct(maturity_type),
    immature = any(maturity_type == "inmaduro"),
    mature   = any(maturity_type == "maduro"),
    .groups = "drop"
  ) %>%
  count(n_maturity)
```

    ## # A tibble: 2 × 2
    ##   n_maturity     n
    ##        <int> <int>
    ## 1          1   126
    ## 2          2   279

Not all species contain both maturity types. If species and
`maturity_type` are largely confounded, much of the difference between
immature and mature material is probably associated with which plant
species occur in each maturity category, rather than maturity alone.

We could try restricting the analysis to only the species which have
both immature and mature fruits.

``` r
species_both <- 
  data_fit %>%
  group_by(codigo) %>%
  filter(n_distinct(maturity_type) == 2) %>%
  ungroup()
```

``` r
m_both <- glmmTMB(
  n_insects ~ maturity_type +
    offset(log(nr_equ_seeds)) +
    (1 | codigo),
  family = nbinom2,
  data =  drop_na(species_both, seeds_per_fruit)
)

summary(m_both)
```

    ##  Family: nbinom2  ( log )
    ## Formula:          n_insects ~ maturity_type + offset(log(nr_equ_seeds)) + (1 |  
    ##     codigo)
    ## Data: drop_na(species_both, seeds_per_fruit)
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##   10591.5   10619.6   -5291.7   10583.5      8391 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  codigo (Intercept) 13.48    3.672   
    ## Number of obs: 8395, groups:  codigo, 274
    ## 
    ## Dispersion parameter for nbinom2 family (): 0.17 
    ## 
    ## Conditional model:
    ##                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)            -6.8935     0.3666  -18.80  < 2e-16 ***
    ## maturity_typeinmaduro  -0.5386     0.1111   -4.85 1.23e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
exp(fixef(m_both)$cond)
```

    ##           (Intercept) maturity_typeinmaduro 
    ##           0.001014396           0.583549771

``` r
exp(confint(m_both, parm = "maturity_typeinmaduro"))
```

    ##                           2.5 %    97.5 %  Estimate
    ## maturity_typeinmaduro 0.4694084 0.7254457 0.5835498

This is similar to the previous model result: the estimated insect
emergence rate for immature material is 0.58 times the rate for mature
material.

What about seeds which show signs of insect attack/ emergence holes?

The metadata say: “Nr_units_predated: Number of units (seeds or fruits)
with clear signs of predation (scored when sample discarded). Unless the
following column says otherwise - all units (seeds or fruits) were
examined.”

``` r
data_seed_attack <- 
  data_seed |> 
  drop_na(seeds_per_fruit) |> 
  # convert n examined fruit to est seeds
  mutate(nr_eq_seeds_examined_for_predation = 
           ifelse(unit == "Fruit", 
                  nr_units_examined_for_predation * seeds_per_fruit,
                  nr_units_examined_for_predation)) |> 
  # if NA, all units were examined
  mutate(nr_eq_seeds_examined_for_predation = 
           coalesce(
             nr_eq_seeds_examined_for_predation, 
             nr_equ_seeds)
         ) |>
  # convert n predated fruit to est seeds
  mutate(nr_eq_seeds_predated = 
           ifelse(
             unit == "Fruit", 
             nr_units_predated * seeds_per_fruit, 
             nr_units_predated)
         ) |> 
  filter(nr_eq_seeds_predated <= nr_eq_seeds_examined_for_predation) |> 
  mutate(
    nr_eq_seeds_examined_for_predation = round(nr_eq_seeds_examined_for_predation),
    nr_eq_seeds_predated = round(nr_eq_seeds_predated)) |> 
  filter(nr_eq_seeds_examined_for_predation > 0) |> 
  mutate(
    maturity_type = factor(maturity_type),
    maturity_type = relevel(maturity_type, ref = "maduro")
  ) |> 
  drop_na(nr_eq_seeds_predated)

data_seed_attack |> 
  group_by(maturity_type) |> 
  summarise(n_predated = sum(nr_eq_seeds_predated, na.rm = TRUE),
            sum_seeds = sum(nr_eq_seeds_examined_for_predation, na.rm = TRUE)) |> 
  mutate(proportion_predated = (n_predated/sum_seeds)*100) 
```

    ## # A tibble: 2 × 4
    ##   maturity_type n_predated sum_seeds proportion_predated
    ##   <fct>              <dbl>     <dbl>               <dbl>
    ## 1 maduro             52110    430503               12.1 
    ## 2 inmaduro           43431   1038095                4.18

``` r
m_attack <- glmmTMB(
    cbind(
        nr_eq_seeds_predated,
        nr_eq_seeds_examined_for_predation - nr_eq_seeds_predated
    ) ~ maturity_type +
        (1 | codigo),
    family = betabinomial(link = "logit"),
    data = data_seed_attack
)

m_attack
```

    ## Formula:          
    ## cbind(nr_eq_seeds_predated, nr_eq_seeds_examined_for_predation -  
    ##     nr_eq_seeds_predated) ~ maturity_type + (1 | codigo)
    ## Data: data_seed_attack
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##  18182.12  18210.15  -9087.06  18174.12      8165 
    ## Random-effects (co)variances:
    ## 
    ## Conditional model:
    ##  Groups Name        Std.Dev.
    ##  codigo (Intercept) 1.027   
    ## 
    ## Number of obs: 8169 / Conditional model: codigo, 355
    ## 
    ## Dispersion parameter for betabinomial family (): 0.885 
    ## 
    ## Fixed Effects:
    ## 
    ## Conditional model:
    ##           (Intercept)  maturity_typeinmaduro  
    ##               -2.9037                -0.1847

``` r
exp(confint(m_attack, parm = "maturity_typeinmaduro"))
```

    ##                           2.5 %    97.5 %  Estimate
    ## maturity_typeinmaduro 0.7433819 0.9297163 0.8313448

The model estimates that the odds of a seed showing signs of predation
in immature material are about 17% lower than in mature material – the
opposite of what we would expect – but a large confidence interval.
