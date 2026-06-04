#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: get-phenological-year.R
## Desc: Calculate phenological year start for each species
## Date created: 2026-01-21

library("tidyverse")
library("here")

species_list <-
  read_csv(here::here("data", "clean", "species_list.csv"))

# seed rain data
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
  inner_join(species_list, by = "sp4") # only keep species in our dataset

# Sum parts per species on each day of year
# mature fruits, immature fruits, single diaspores and capsules
summed_parts <-
  trap_dat %>%
  filter(part == 1 | part == 2| part == 3| part == 5) %>%
  group_by(sp4, day_of_year) %>%
  summarise(quantity_sum = sum(quantity, na.rm = TRUE)) %>%
  ungroup()

# Create function to calculate variance for each species & day
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

# Run function
out_parts <- rolling_variance_by_species(summed_parts)

# Now find the day with minimum variance for each species
# if multiple days with the same variance, calculate the median day
min_vars <-
  out_parts %>%
  group_by(sp4) %>%
  slice_min(variance, with_ties = TRUE) %>%
  summarise(min_var_day = floor(median(beginyr)))


# Check results with a figure for each species ----------------------------

# sum by parts for plotting
summed_parts <-
  trap_dat %>%
  filter(part == 1 | part == 2| part == 3| part == 5) %>%
  group_by(sp4, day_of_year, part) %>%
  summarise(quantity_sum = sum(quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(part = case_when(
    part == 1 ~ "1 - Mature fruit",
    part == 5 ~ "5 - Immature fruit",
    part == 3 ~ "3 - Capsules",
    part == 2 ~ "2 - Seeds"
  )) %>%
  mutate(part = as.ordered(part))

# Function to make a single plot
make_sp_plot <- function(fruit_data,
                         variance_data,
                         species_code,
                         species_name) {

  line_df <- variance_data %>%
    filter(sp4 == species_code, !is.na(min_var_day))

  fruit_data %>%
    filter(sp4 == species_code) %>%
    ggplot(aes(x = day_of_year, y = quantity_sum,
               colour = part, fill = part)) +
    geom_col(position = "stack") +
    geom_vline(
      data = line_df,
      aes(xintercept = min_var_day),
      colour = "red",
      linetype = 2
    ) +
    labs(
      title = species_name,
      y = "Number of parts"
    ) +
    xlim(c(0, 365)) +
    scale_colour_viridis_d(drop = FALSE) +
    scale_fill_viridis_d(drop = FALSE)
}

# Get a list of species with latin names
all_sp <-
  trap_dat %>%
  drop_na(sp4) %>%
  filter(sp4 != "nada") %>%
  distinct(sp4) %>%
  left_join(species_list) %>%
  mutate(name = paste(sp4, "-", genus, species, sep = " ")) %>%
  arrange(sp4)

# Create a plot for each species
plot_list <- purrr::pmap(
  .f = make_sp_plot,
  fruit_data = summed_parts,
  variance_data = min_vars,
  list(
    species_code = all_sp$sp4,
    species_name = all_sp$name),
  .progress = TRUE
)

# Split into groups of 4 plots
plot_list_split <-
  split(plot_list, rep(
    seq_along(plot_list),
    each = 4,
    length.out = length(plot_list)
  ))

# Wrap each group of 4 plots for 1 pdf page
plots_wrapped <-
  purrr::map(
    .x = plot_list_split,
    .f = patchwork::wrap_plots,
    guides = "collect",
    nrow = 2,
    ncol = 2
  )

# Make pdf
pdf(
  width = 12,
  height = 8,
  file = here::here("docs", "pheno_years.pdf")
)
plots_wrapped
dev.off()


# Save output data --------------------------------------------------------

min_vars %>%
  rename(pheno_year_start = min_var_day) %>%
  write_csv(here::here("data", "clean", "phenological_years.csv"))
