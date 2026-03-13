## --------------------------------------------------------------#
## Script name: script01-03_gap_analysis_width.R
##
## Purpose:
##   - gap-analysis table.
##   - For each species:
##       * n (sample size)
##       * width range (min/max)
##       * counts above/below 5 cm (50 mm)
##       * proportional representation
##       * meaningful gap flags:
##            - TRUE if below/above 5 cm is missing or too weak
##            - uses MIN_REQUIRED and LOW_N_THRESHOLD rules
##
## Author: Marcus Rizzuto
## Date Created: 3/12/2026
##
## --------------------------------------------------------------#

##### Setup --------------------------------------------------#
WIDTH_THRESHOLD_MM  <- 50   # 5 cm
MIN_REQUIRED        <- 5    # need ≥5 fish to say “represented”
LOW_N_THRESHOLD     <- 10   # species with n<10 flagged as unreliable

if (!exists("df_all") && !exists("combined_all")) {
 stop("Need df_all or combined_all first.")
}

library(dplyr)
library(tibble)

##### Build unified df ---------------------------------------#
if (exists("combined_all")) {
 df_width <- purrr::imap_dfr(
  combined_all,
  ~ tibble(
   Species  = .y,
   Width_mm = .x$Width_mm
  )
 )
} else {
 df_width <- df_all %>% select(Species, Width_mm)
}

df_width <- df_width %>% filter(!is.na(Width_mm))

##### Build final table --------------------------------------#
df_width_gap <- df_width %>%
 group_by(Species) %>%
 summarise(
  n = n(),
  min_width_mm = min(Width_mm),
  max_width_mm = max(Width_mm),
  
  n_below_5cm  = sum(Width_mm < WIDTH_THRESHOLD_MM),
  n_above_5cm  = sum(Width_mm >= WIDTH_THRESHOLD_MM),
  
  pct_below = n_below_5cm / n,
  pct_above = n_above_5cm / n,
  
  # TRUE = there is a real gap (missing or very weak)
  gap_below_5cm = n_below_5cm < MIN_REQUIRED,
  gap_above_5cm = n_above_5cm < MIN_REQUIRED,
  
  # override: species too small overall to make claims
  low_n_flag = n < LOW_N_THRESHOLD,
  
  .groups = "drop"
 ) %>%
 arrange(Species)

# ---- standard deviation  ----
df_width_sd <- df_width %>%
 dplyr::group_by(Species) %>%
 dplyr::summarise(
  sd_width_mm = sd(Width_mm),
  .groups = "drop"
 )

# Append SD to the existing gap table
df_width_gap <- df_width_gap %>%
 dplyr::left_join(df_width_sd, by = "Species")
# ---- End SD chunk ----

print(df_width_gap)
cat("\nGap table () complete.\n")

## Optional export:
 write.csv(df_width_gap, "03_outputs/01_tables/width_gap_table.csv", row.names = FALSE)