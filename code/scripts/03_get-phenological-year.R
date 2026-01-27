#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: get-phenological-year.R
## Desc: calculate phenological year for each species
## Date created: 2026-01-21

library("tidyverse")
library("here")

species_list <-
  read_csv(here::here("data", "clean", "species_list.csv"))

# seed rain data - fruits only
trap_dat <-
  foreign::read.dbf(
  here::here("data", "raw", "BCI_TRAP200_20250224.DBF")) %>%
  mutate(MASS = na_if(MASS, -9)) %>%
  rename_with(tolower) %>%
  select(-c(mass, seedeqs)) %>%
  mutate(date = as.Date(date, "%Y-%m-%d"),
         calendar_year = format(as.Date(date), "%Y"),
         day_of_year = yday(date) - 1) %>%
  # before Nov 1989 immature fruits category also contained aborted fruits
  filter(!calendar_year %in% c("1987", "1988", "1989")) %>%
  mutate(species = tolower(species),
         trap = formatC(
           trap,
           width = 3,
           format = "d",
           flag = "0"
         )) %>%
  mutate(trap = paste("trap", trap, sep = "_")) %>%
  rename(sp4 = species) %>%
  inner_join(species_list, by = "sp4")

# sum fruits per species on each day of year
summed_fruits <-
  trap_dat %>%
  filter(part == 1 | part == 5| part == 3) %>% # only fruits (mature, immature and capsules)
  group_by(sp4, day_of_year) %>%
  summarise(quantity_sum = sum(quantity, na.rm = TRUE)) %>%
  ungroup()

# sum flowers per species on each day of year
summed_flowers <-
  trap_dat %>%
  filter(part == 6 | part == 9) %>% # only flowers (male and female)
  group_by(sp4, day_of_year) %>%
  summarise(quantity_sum = sum(quantity, na.rm = TRUE)) %>%
  ungroup()


#function
rolling_variance_by_species <- function(df,
                                        day_col = day_of_year,
                                        qty_col = quantity_sum,
                                        species_col = sp4,
                                        days = 0:364) {
  # enquo() captures the column expressions so they can be evaluated later
  # inside dplyr verbs
  day_col <- enquo(day_col)
  qty_col <- enquo(qty_col)
  species_col <- enquo(species_col)

  # define helper function
  # The phenological year could start on *any* calendar day.
  # For each candidate start day (i = 0..364), we:
  #   1) re-index all observations relative to that start day
  #   2) compute how compact fruiting is for each species
  #
  # This helper performs that calculation for ONE candidate start day (i).
  # It will later be mapped across all possible days.
  compute_for_shift <- function(i) {
    df %>%

      # Day-of-year is circular (Dec 31 → Jan 1).
      # To evaluate a candidate phenological year start (i), we:
      #   - subtract i from every observation
      #   - wrap values using modulo so they stay within 0..364
      #
      # After this transformation:
      #   x = 0 means "exactly at the phenological year start"
      #   small x means "soon after the start"
      #   large x means "late in the phenological year"
      #
      # This converts a circular timing problem into a linear one.
      mutate(
        # re-index day of year
        x = ((!!day_col - i) %% 365),
        # get weighted sums of fruit per day
        x2 = x^2,
        xsum = x * (!!qty_col),
        x2sum = x2 * (!!qty_col)
      ) %>%
      # summarise by species
      group_by(!!species_col) %>%
      summarise(
        quantity = sum(!!qty_col, na.rm = TRUE),
        xsum = sum(xsum, na.rm = TRUE),
        x2sum = sum(x2sum, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      # Compute the variance
      mutate(
        variance = (x2sum / quantity) - (xsum / quantity) ^ 2,
        beginyr = i
      )
  }

  # run compute_for_shift for each day and bind rows
  result <- map_dfr(days, compute_for_shift) %>%
    arrange(sp4, beginyr) %>%
    select(beginyr, sp4, quantity, variance)

  result
}

# Fun function for both flowers and fruits
out_fl <- rolling_variance_by_species(summed_flowers)
out_fr <- rolling_variance_by_species(summed_fruits)

# now find the day with minimum variance for each species
# if multiple days with the same variance, calculate the median day
min_vars_fl <-
  out_fl %>%
  group_by(sp4) %>%
  slice_min(variance, with_ties = TRUE) %>%
  summarise(min_var_day = floor(median(beginyr)))

min_vars_fr <-
  out_fr %>%
  group_by(sp4) %>%
  slice_min(variance, with_ties = TRUE) %>%
  summarise(min_var_day = floor(median(beginyr)))

# check with plots

make_sp_plot <- function(fruit_data,
                         fruit_variance_data,
                         flower_variance_data,
                         species_code,
                         species_name) {

  fruit_line_df <- fruit_variance_data %>%
    filter(sp4 == species_code, !is.na(min_var_day))

  flower_line_df <- flower_variance_data %>%
    filter(sp4 == species_code, !is.na(min_var_day))

  fruit_data %>%
    filter(sp4 == species_code) %>%
    ggplot(aes(x = day_of_year, y = quantity_sum)) +
    geom_col() +
    geom_vline(
      data = fruit_line_df,
      aes(xintercept = min_var_day),
      colour = "blue",
      linetype = 2
    ) +
    geom_vline(
      data = flower_line_df,
      aes(xintercept = min_var_day),
      colour = "red",
      linetype = 2
    ) +
    labs(
      title = species_name,
      y = "Summed fruits (mature, immature capsules)"
    ) +
    xlim(c(0, 365))
}

all_sp <-
  trap_dat %>%
  drop_na(sp4) %>%
  filter(sp4 != "nada") %>%
  distinct(sp4) %>%
  left_join(species_list) %>%
  mutate(name = paste(sp4, "-", genus, species, sep = " ")) %>%
  arrange(sp4)

# create many plots
plot_list <- purrr::pmap(
  .f = make_sp_plot,
  fruit_data = summed_fruits,
  fruit_variance_data = min_vars_fr,
  flower_variance_data = min_vars_fl,
  list(
    species_code = all_sp$sp4,
    species_name = all_sp$name),
  .progress = TRUE
)

# split into groups of 4 plots
plot_list_split <-
  split(plot_list, rep(
    seq_along(plot_list),
    each = 4,
    length.out = length(plot_list)
  ))

# wrap each group of 4 plots for 1 pdf page
plots_wrapped <-
  purrr::map(
    .x = plot_list_split,
    .f = patchwork::wrap_plots,
    nrow = 2,
    ncol = 2
  )

# make pdf
pdf(
  width = 12,
  height = 8,
  file = here::here("docs", "pheno_years.pdf")
)
plots_wrapped
dev.off()

# save fruiting pheno data
min_vars_fr %>%
  rename(pheno_year_start = min_var_day) %>%
  write_csv(here::here("data", "clean", "phenological_years.csv"))
